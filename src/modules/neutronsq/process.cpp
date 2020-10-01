// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Team Dissolve and contributors

#include "classes/atomtype.h"
#include "classes/box.h"
#include "classes/configuration.h"
#include "classes/neutronweights.h"
#include "classes/species.h"
#include "genericitems/listhelper.h"
#include "io/export/data1d.h"
#include "main/dissolve.h"
#include "math/filters.h"
#include "math/ft.h"
#include "modules/import/import.h"
#include "modules/neutronsq/neutronsq.h"
#include "modules/rdf/rdf.h"
#include "modules/sq/sq.h"

// Run set-up stage
bool NeutronSQModule::setUp(Dissolve &dissolve, ProcessPool &procPool)
{
    /*
     * Load and set up reference data (if a file/format was given)
     */
    if (referenceFQ_.hasValidFileAndFormat())
    {
        // Load the data
        Data1D referenceData;
        if (!referenceFQ_.importData(referenceData, &procPool))
        {
            Messenger::error("Failed to load reference data '{}'.\n", referenceFQ_.filename());
            return false;
        }

        // Truncate data beyond QMax
        const auto qMax = keywords_.asDouble("QMax") < 0.0 ? 30.0 : keywords_.asDouble("QMax");
        if (referenceData.constXAxis().lastValue() < qMax)
            Messenger::warn("Qmax limit of {:e} Angstroms**-1 for calculated NeutronSQ ({}) is beyond limit of "
                            "reference data (Qmax = {:e} Angstroms**-1).\n",
                            qMax, uniqueName(), referenceData.constXAxis().lastValue());
        else
            while (referenceData.constXAxis().lastValue() > qMax)
                referenceData.removeLastPoint();

        // Remove first point?
        if (keywords_.asBool("ReferenceIgnoreFirst"))
        {
            referenceData.removeFirstPoint();
            Messenger::print("Removed first point from supplied reference data - new Qmin = {:e} Angstroms**-1.\n",
                             referenceData.constXAxis().firstValue());
        }

        // Remove normalisation factor from data
        auto normType = keywords_.enumeration<NeutronSQModule::NormalisationType>("ReferenceNormalisation");
        if (normType != NeutronSQModule::NoNormalisation)
        {
            // We need the summed Weights in order to do the normalisation
            NeutronWeights summedWeights;
            if (!calculateSummedWeights(summedWeights))
            {
                Messenger::error("Couldn't get summed Weights for reference data in NeutronSQ module '{}', and "
                                 "so can't perform requested normalisation.\n",
                                 uniqueName());
                return false;
            }

            // Remove normalisation of data
            if (normType == NeutronSQModule::AverageOfSquaresNormalisation)
            {
                referenceData.values() *= summedWeights.boundCoherentAverageOfSquares();
                Messenger::print("NeutronSQ: Removed <b>**2 normalisation from reference data ('{}'), factor = {}.\n",
                                 uniqueName(), summedWeights.boundCoherentAverageOfSquares());
            }
            else if (normType == NeutronSQModule::SquareOfAverageNormalisation)
            {
                referenceData.values() *= summedWeights.boundCoherentSquareOfAverage();
                Messenger::print("NeutronSQ: Removed <b**2> normalisation from reference data ('{}'), factor = {}.\n",
                                 uniqueName(), summedWeights.boundCoherentSquareOfAverage());
            }
        }

        // Get window function to use for transformation of S(Q) to g(r)
        const WindowFunction &referenceWindowFunction =
            keywords_.retrieve<WindowFunction>("ReferenceWindowFunction", WindowFunction());
        if (referenceWindowFunction.function() == WindowFunction::NoWindow)
            Messenger::print("No window function will be applied in Fourier transform of S(Q) to g(r).");
        else
            Messenger::print("Window function to be applied in Fourier transform of reference data is {} ({}).",
                             WindowFunction::functionType(referenceWindowFunction.function()),
                             referenceWindowFunction.parameterSummary());

        // Store the reference data in processing
        referenceData.setName(uniqueName());
        auto &storedData = GenericListHelper<Data1D>::realise(dissolve.processingModuleData(), "ReferenceData", uniqueName(),
                                                              GenericItem::ProtectedFlag);
        storedData.setObjectTag(fmt::format("{}//ReferenceData", uniqueName()));
        storedData = referenceData;

        // Calculate and store the FT of the reference data in processing
        referenceData.setName(uniqueName());
        auto &storedDataFT = GenericListHelper<Data1D>::realise(dissolve.processingModuleData(), "ReferenceDataFT",
                                                                uniqueName(), GenericItem::ProtectedFlag);
        storedDataFT.setObjectTag(fmt::format("{}//ReferenceDataFT", uniqueName()));
        storedDataFT = referenceData;
        double rho = nTargetConfigurations() == 0 ? 0.1 : RDFModule::summedRho(this, dissolve.processingModuleData());
        Fourier::sineFT(storedDataFT, 1.0 / (2.0 * PI * PI * rho), 0.0, 0.05, 30.0, referenceWindowFunction);
        if (nTargetConfigurations() == 0)
            Messenger::warn("No configurations associated to module, so Fourier transform of reference data will "
                            "use assumed atomic density of 0.1.\n");

        // Save data?
        if (keywords_.asBool("SaveReference"))
        {
            if (procPool.isMaster())
            {
                Data1DExportFileFormat exportFormat(fmt::format("{}-ReferenceData.q", uniqueName()));
                if (!exportFormat.exportData(storedData))
                    return procPool.decideFalse();
                Data1DExportFileFormat exportFormatFT(fmt::format("{}-ReferenceData.r", uniqueName()));
                if (!exportFormatFT.exportData(storedDataFT))
                    return procPool.decideFalse();
                procPool.decideTrue();
            }
            else if (!procPool.decision())
                return false;
        }
    }

    return true;
}

// Run main processing
bool NeutronSQModule::process(Dissolve &dissolve, ProcessPool &procPool)
{
    /*
     * Calculate neutron structure factors from existing g(r) data
     *
     * This is a serial routine, with each process constructing its own copy of the data.
     * Partial calculation routines called by this routine are parallel.
     */

    // Check for zero Configuration targets
    if (targetConfigurations_.nItems() == 0)
        return Messenger::error("No configuration targets set for module '{}'.\n", uniqueName());

    const SQModule *sqModule = keywords_.retrieve<const SQModule *>("SourceSQs", nullptr);
    if (!sqModule)
        return Messenger::error("A source SQ module must be provided.\n");
    const RDFModule *rdfModule = sqModule->keywords().retrieve<const RDFModule *>("SourceRDFs", nullptr);
    if (!rdfModule)
        return Messenger::error("A source RDF module (in the SQ module) must be provided.\n");
    auto normalisation = keywords_.enumeration<NeutronSQModule::NormalisationType>("Normalisation");
    const auto &qBroadening = keywords_.retrieve<BroadeningFunction>("QBroadening", BroadeningFunction());
    const auto qDelta = keywords_.asDouble("QDelta");
    const auto qMin = keywords_.asDouble("QMin");
    double qMax = keywords_.asDouble("QMax");
    if (qMax < 0.0)
        qMax = 30.0;
    const bool saveUnweighted = keywords_.asBool("SaveUnweighted");
    const bool saveWeighted = keywords_.asBool("SaveWeighted");
    const auto &windowFunction = keywords_.retrieve<WindowFunction>("WindowFunction", WindowFunction());

    // Print argument/parameter summary
    Messenger::print("NeutronSQ: Calculating S(Q)/F(Q) over {} < Q < {} Angstroms**-1 using step size of {} Angstroms**-1.\n",
                     qMin, qMax, qDelta);
    Messenger::print("NeutronSQ: Source unweighted S(Q) will be taken from module '{}'.\n", sqModule->uniqueName());
    if (windowFunction.function() == WindowFunction::NoWindow)
        Messenger::print("NeutronSQ: No window function will be applied in Fourier transforms of g(r) to S(Q).");
    else
        Messenger::print("NeutronSQ: Window function to be applied in Fourier transforms is {} ({}).",
                         WindowFunction::functionType(windowFunction.function()), windowFunction.parameterSummary());
    const WindowFunction &referenceWindowFunction =
        keywords_.retrieve<WindowFunction>("ReferenceWindowFunction", WindowFunction());
    if (referenceWindowFunction.function() == WindowFunction::NoWindow)
        Messenger::print("No window function will be applied when calculating representative g(r) from S(Q).");
    else
        Messenger::print("Window function to be applied when calculating representative g(r) from S(Q) is {} ({}).",
                         WindowFunction::functionType(referenceWindowFunction.function()),
                         referenceWindowFunction.parameterSummary());
    if (normalisation == NeutronSQModule::NoNormalisation)
        Messenger::print("NeutronSQ: No normalisation will be applied to total F(Q).\n");
    else if (normalisation == NeutronSQModule::AverageOfSquaresNormalisation)
        Messenger::print("NeutronSQ: Total F(Q) will be normalised to <b>**2");
    else if (normalisation == NeutronSQModule::SquareOfAverageNormalisation)
        Messenger::print("NeutronSQ: Total F(Q) will be normalised to <b**2>");
    if (qBroadening.function() == BroadeningFunction::NoFunction)
        Messenger::print("NeutronSQ: No broadening will be applied to calculated S(Q).");
    else
        Messenger::print("NeutronSQ: Broadening to be applied in calculated S(Q) is {} ({}).",
                         BroadeningFunction::functionType(qBroadening.function()), qBroadening.parameterSummary());
    if (saveUnweighted)
        Messenger::print("NeutronSQ: Unweighted partials and totals will be saved.\n");
    if (saveWeighted)
        Messenger::print("NeutronSQ: Weighted partials and totals will be saved.\n");
    Messenger::print("\n");

    /*
     * Loop over target Configurations and Fourier transform their UnweightedGR into the corresponding UnweightedSQ.
     */

    bool created;

    for (Configuration *cfg : targetConfigurations_)
    {
        // Set up process pool - must do this to ensure we are using all available processes
        procPool.assignProcessesToGroups(cfg->processPool());

        // Get unweighted g(r) for this Configuration - we don't supply a specific Module prefix, since the unweighted
        // g(r) may come from one of many RDF-type modules
        if (!cfg->moduleData().contains("UnweightedGR", rdfModule->uniqueName()))
            return Messenger::error("Couldn't locate UnweightedGR for Configuration '{}'.\n", cfg->name());
        const auto &unweightedgr = GenericListHelper<PartialSet>::value(cfg->moduleData(), "UnweightedGR");

        // Does a PartialSet for the unweighted S(Q) already exist for this Configuration?
        auto &unweightedsq = GenericListHelper<PartialSet>::realise(cfg->moduleData(), "UnweightedSQ", "",
                                                                    GenericItem::InRestartFileFlag, &created);
        if (created)
            unweightedsq.setUpPartials(unweightedgr.atomTypes(), fmt::format("{}-{}", cfg->niceName(), uniqueName()),
                                       "unweighted", "sq", "Q, 1/Angstroms");

        // Is the PartialSet already up-to-date? Do we force its calculation anyway?
        auto &forceCalculation = GenericListHelper<bool>::retrieve(cfg->moduleData(), "_ForceNeutronSQ", "", false);
        const bool sqUpToDate = DissolveSys::sameString(
            unweightedsq.fingerprint(), fmt::format("{}", cfg->moduleData().version("UnweightedGR")));
        if ((!forceCalculation) && sqUpToDate)
        {
            Messenger::print("Unweighted partial S(Q) are up-to-date for Configuration '{}'.\n", cfg->name());
            continue;
        }
        forceCalculation = false;

        // Transform g(r) into S(Q)
        if (!SQModule::calculateUnweightedSQ(procPool, cfg, unweightedgr, unweightedsq, qMin, qDelta, qMax,
                                             cfg->atomicDensity(), windowFunction, qBroadening))
            return false;

        // Set names of resources (Data1D) within the PartialSet, and tag it with the fingerprint from the source
        // unweighted g(r)
        unweightedsq.setObjectTags(fmt::format("{}//{}//{}", cfg->niceName(), "NeutronSQ", "UnweightedSQ"));
        unweightedsq.setFingerprint(fmt::format("{}", cfg->moduleData().version("UnweightedGR")));

        // Save data if requested
        if (saveUnweighted && (!MPIRunMaster(procPool, unweightedsq.save())))
            return false;

        // Construct weights matrix based on Isotopologue specifications and the populations of AtomTypes in the
        // Configuration
        NeutronWeights weights;
        if (isotopologues_.contains(cfg))
        {
            // Get the set...
            auto topeSetData = isotopologues_.getIsotopologueSet(cfg);
            if (!topeSetData)
                return Messenger::error("Could not locate IsotopologueSet");
            const IsotopologueSet &topeSet = *topeSetData;

            // Iterate over Species present in the set
            for (auto &topes : topeSet.constIsotopologues())
            {
                // Find the referenced Species in our SpeciesInfo list
                auto *spInfo = cfg->usedSpeciesInfo(topes.species());
                if (!spInfo)
                    return Messenger::error("Couldn't locate SpeciesInfo for '{}' in the Configuration '{}'.\n",
                                            topes.species()->name(), cfg->niceName());

                // Add defined isotopologues, in the relative isotopic proportions defined, to the weights.
                for (const auto &isoWeight : topes.constMix())
                    weights.addIsotopologue(spInfo->species(), spInfo->population(), isoWeight.isotopologue(),
                                            isoWeight.weight());
            }
        }

        // We will automatically use the natural isotopologue for a species if it is not explicitly covered by at least
        // one Isotopologue definition
        ListIterator<SpeciesInfo> speciesInfoIterator(cfg->usedSpecies());
        while (SpeciesInfo *spInfo = speciesInfoIterator.iterate())
            if (!weights.containsIsotopologues(spInfo->species()))
            {
                Messenger::print("Isotopologue specification for Species '{}' in Configuration '{}' is missing "
                                 "- natural isotopologue will be used.\n",
                                 spInfo->species()->name(), cfg->name());

                Species *sp = spInfo->species();
                weights.addIsotopologue(sp, spInfo->population(), sp->naturalIsotopologue(), 1.0);
            }

        // Create, print, and store weights
        Messenger::print("Isotopologue and isotope composition for Configuration '{}':\n\n", cfg->name());
        weights.createFromIsotopologues(exchangeableTypes_);
        weights.print();

        // Does a PartialSet for the unweighted S(Q) already exist for this Configuration?
        auto &weightedsq = GenericListHelper<PartialSet>::realise(cfg->moduleData(), "WeightedSQ", "",
                                                                  GenericItem::InRestartFileFlag, &created);
        if (created)
            weightedsq.setUpPartials(unweightedsq.atomTypes(), fmt::format("{}-{}", cfg->niceName(), uniqueName()), "weighted",
                                     "sq", "Q, 1/Angstroms");

        // Calculate weighted S(Q)
        calculateWeightedSQ(unweightedsq, weightedsq, weights, normalisation);

        // Set names of resources (Data1D) within the PartialSet
        weightedsq.setObjectTags(fmt::format("{}//{}//{}", cfg->niceName(), uniqueName_, "WeightedSQ"));

        // Save data if requested
        if (saveWeighted && (!MPIRunMaster(procPool, weightedsq.save())))
            return false;
    }

    /*
     * Construct the weighted sum of Configuration partials and store in the processing module data list
     */

    // Create/retrieve PartialSet for summed partial S(Q)
    auto &summedUnweightedSQ = GenericListHelper<PartialSet>::realise(dissolve.processingModuleData(), "UnweightedSQ",
                                                                      uniqueName_, GenericItem::InRestartFileFlag);

    // Sum the unweighted S(Q) from the associated Configurations
    if (!SQModule::sumUnweightedSQ(procPool, this, sqModule, dissolve.processingModuleData(), summedUnweightedSQ))
        return false;

    // Save data if requested
    if (saveUnweighted && (!MPIRunMaster(procPool, summedUnweightedSQ.save())))
        return false;

    // Create/retrieve PartialSet for summed partial S(Q)
    auto &summedWeightedSQ = GenericListHelper<PartialSet>::realise(dissolve.processingModuleData(), "WeightedSQ", uniqueName_,
                                                                    GenericItem::InRestartFileFlag);
    summedWeightedSQ = summedUnweightedSQ;
    summedWeightedSQ.setFileNames(uniqueName_, "weighted", "sq");
    summedWeightedSQ.setObjectTags(fmt::format("{}//{}", uniqueName_, "WeightedSQ"));

    // Calculate weighted S(Q)
    Messenger::print("Isotopologue and isotope composition over all Configurations used in '{}':\n\n", uniqueName_);
    NeutronWeights summedWeights;
    if (!calculateSummedWeights(summedWeights))
        return false;
    summedWeights.print();
    GenericListHelper<NeutronWeights>::realise(dissolve.processingModuleData(), "FullWeights", uniqueName_,
                                               GenericItem::InRestartFileFlag) = summedWeights;
    calculateWeightedSQ(summedUnweightedSQ, summedWeightedSQ, summedWeights, normalisation);

    // Create/retrieve PartialSet for summed unweighted g(r)
    auto &summedUnweightedGR = GenericListHelper<PartialSet>::realise(dissolve.processingModuleData(), "UnweightedGR",
                                                                      uniqueName_, GenericItem::InRestartFileFlag);

    // Sum the partials from the associated Configurations
    if (!RDFModule::sumUnweightedGR(procPool, this, rdfModule, dissolve.processingModuleData(), summedUnweightedGR))
        return false;

    // Create/retrieve PartialSet for summed weighted g(r)
    auto &summedWeightedGR = GenericListHelper<PartialSet>::realise(dissolve.processingModuleData(), "WeightedGR", uniqueName_,
                                                                    GenericItem::InRestartFileFlag, &created);
    if (created)
        summedWeightedGR.setUpPartials(summedUnweightedSQ.atomTypes(), uniqueName_, "weighted", "gr", "r, Angstroms");
    summedWeightedGR.setObjectTags(fmt::format("{}//{}", uniqueName_, "WeightedGR"));

    // Calculate weighted g(r)
    calculateWeightedGR(summedUnweightedGR, summedWeightedGR, summedWeights, normalisation);

    // Calculate representative total g(r) from FT of calculated S(Q)
    auto &repGR = GenericListHelper<Data1D>::realise(dissolve.processingModuleData(), "RepresentativeTotalGR", uniqueName_,
                                                     GenericItem::InRestartFileFlag);
    repGR = summedWeightedSQ.total();
    double rho = nTargetConfigurations() == 0 ? 0.1 : RDFModule::summedRho(this, dissolve.processingModuleData());
    Fourier::sineFT(repGR, 1.0 / (2.0 * PI * PI * rho), qMin, qDelta, qMax, referenceWindowFunction);
    repGR.setObjectTag(fmt::format("{}//RepresentativeTotalGR", uniqueName_));

    // Save data if requested
    if (saveWeighted && (!MPIRunMaster(procPool, summedWeightedSQ.save())))
        return false;

    return true;
}
