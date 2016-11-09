/*
	*** Partials Module
	*** src/modules/partials.cpp
	Copyright T. Youngs 2012-2016

	This file is part of dUQ.

	dUQ is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	dUQ is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with dUQ.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "main/duq.h"
#include "modules/partials.h"
#include "classes/box.h"
#include "classes/grain.h"
#include "classes/species.h"
#include "base/sysfunc.h"
#include "base/processpool.h"
#include "base/timer.h"

// Static Members
List<Module> Partials::instances_;
List<PartialSet> Partials::partialSets_;

/*
 * Constructor / Destructor
 */

// Constructor
Partials::Partials() : Module()
{
	// Setup variables / control parameters
	addVariable("Save", false, "Whether to save partials to disk after calculation");
}	

// Destructor
Partials::~Partials()
{
}

/*
 * Instances
 */

// Create instance of this module
List<Module>& Partials::instances()
{
	return instances_;
}

// Create instance of this module
Module* Partials::createInstance()
{
	return new Partials;
}

/*
 * Definition
 */

// Return name of module
const char* Partials::name()
{
	return "Partials";
}

// Return brief description of module
const char* Partials::brief()
{
	return "Calculate partial and total RDFs";
}

// Return instance type for module
Module::InstanceType Partials::instanceType()
{
	return Module::MultipleInstance;
}

// Whether the Module has a pre-processing stage
bool Partials::hasPreProcessing()
{
	return false;
}

// Whether the Module has a processing stage
bool Partials::hasProcessing()
{
	return true;
}

// Whether the Module has a post-processing stage
bool Partials::hasPostProcessing()
{
	return false;
}

/*
 * Targets
 */

// Return the maximum number of Configurations the Module can target (or -1 for any number)
int Partials::nTargetableConfigurations()
{
	return 1;
}

// Return the maximum number of Samples the Module can target (or -1 for any number)
int Partials::nTargetableSamples()
{
	return 0;
}

/*
 * Method
 */

// Perform setup tasks for module
bool Partials::setup(ProcessPool& procPool)
{
	return true;
}

// Execute pre-processing stage
bool Partials::preProcess(DUQ& duq, ProcessPool& procPool)
{
	return false;
}

// Execute Method
bool Partials::process(DUQ& duq, ProcessPool& procPool)
{
	/*
	 * Calculate weighted or unweighted partials and total RDFs for the Sample or Configuration
	 * 
	 * This is a serial routine, with each process constructing its own copy of the data.
	 * Partial calculation routines called by this routine are parallel.
	 */

	// If there is a Sample target, then we calculate the weighted structure factors for it (using the supplied Configurations)
	// Otherwise take the Configuration targets and calculate unweighted structure factors for them.

	if (targetSamples_.nItems() > 0)
	{
	}
	else
	{
		RefListIterator<Configuration,bool> configIterator(targetConfigurations_);
		while (Configuration* cfg = configIterator.iterate())
		{
			// Retrieve control parameters from Configuration
			const bool saveData = variableAsBool("Save");

			// Print argument/parameter summary
			Messenger::print("Partials: Save data is %s.\n", DUQSys::onOff(saveData));

			Partials::calculateUnweighted(cfg, procPool);

			// Save data?
			if (saveData)
			{
				// Only the pool master saves the data
				if (procPool.isMaster())
				{
					// Find PartialSet for this Configuration
					PartialSet* partials = partialSet(cfg);
					if (partials->save()) procPool.proceed();
					else
					{
						procPool.stop();
						return false;
					}
				}
				else if (!procPool.decision()) return false;
			}
		}
	}

	return true;
}

// Execute post-processing stage
bool Partials::postProcess(DUQ& duq, ProcessPool& procPool)
{
	return false;
}

/*
 * Static Functions
 */

// Calculate partial RDFs with simple double-loop
bool Partials::calculateSimple(PartialSet* partialSet, ProcessPool& procPool)
{
	// Variables
	int n, m, nTypes, typeI, typeJ, i, j, nPoints;

	// Construct local arrays of atom type positions
	Messenger::printVerbose("Constructing local partial working arrays.\n");
	nTypes = partialSet->nTypes();
	Configuration* targetConfiguration = partialSet->targetConfiguration();
	const Box* box = targetConfiguration->box();
	Vec3<double>* r[nTypes];
	int maxr[nTypes], nr[nTypes];
	int* binss[nTypes], *bins;

	n = 0;
	for (AtomTypeData* atd = targetConfiguration->usedAtomTypes(); atd != NULL; atd = atd->next)
	{
		maxr[n] = atd->population();
		nr[n] = 0;
		r[n] = new Vec3<double>[maxr[n]];
		binss[n] = new int[maxr[n]];
		++n;
	}

	// Loop over Atoms and construct arrays
	Atom* atoms = targetConfiguration->atoms();
	for (n=0; n<targetConfiguration->nAtoms(); ++n)
	{
		m = atoms[n].localTypeIndex();
		r[m][nr[m]++] = atoms[n].r();
	}

	Messenger::printVerbose("Ready..\n");

	// Loop over assigned Atoms
	Vec3<double> centre, *ri, *rj, mim;
	int* histogram;
	double rbin = 1.0 / targetConfiguration->rdfBinWidth();

	// Loop context is to use all processes in Pool as one group
	int start = procPool.interleavedLoopStart(ProcessPool::OverPoolProcesses);
	int stride = procPool.interleavedLoopStride(ProcessPool::OverPoolProcesses);

	Messenger::printVerbose("Self terms..\n");
	
	// Self terms
	for (typeI = 0; typeI<nTypes; ++typeI)
	{
		ri = r[typeI];
		histogram = partialSet->partial(typeI,typeI).histogram();
		bins = binss[typeI];
		nPoints = partialSet->partial(typeI,typeI).nBins();
		for (i=start; i < maxr[typeI]; i += stride)
		{
			centre = ri[i];
			for (j = i+1; j < maxr[typeI]; ++j) bins[j] = box->minimumDistance(centre, ri[j]) * rbin;
			for (j = i+1; j < maxr[typeI]; ++j) if (bins[j] < nPoints) ++histogram[bins[j]];
		}
// 		printf("For types %i-%i, count = %i\n", typeI, typeI, count);
	}

	Messenger::printVerbose("Cross terms..\n");

	// Cross terms
	for (typeI = 0; typeI<nTypes; ++typeI)
	{
		ri = r[typeI];
		for (typeJ = 0; typeJ<nTypes; ++typeJ)
		{
			// Skip if typeI == typeJ, or if the number of atoms in typeI is greater than typeJ (since it is less efficient)
			if (typeI == typeJ) continue;
			if (nr[typeI] > nr[typeJ]) continue;
			if ((nr[typeI] == nr[typeJ]) && (typeI > typeJ)) continue;

			rj = r[typeJ];
			histogram = partialSet->partial(typeI,typeJ).histogram();
			bins = binss[typeJ];
			nPoints = partialSet->partial(typeI,typeJ).nBins();
			for (i=start; i < maxr[typeI]; i += stride)
			{
				centre = ri[i];
				for (j = 0; j < maxr[typeJ]; ++j) bins[j] = box->minimumDistance(centre, rj[j]) * rbin;
				for (j = 0; j < maxr[typeJ]; ++j) if (bins[j] < nPoints) ++histogram[bins[j]];
			}
// 			printf("For types %i-%i count = %i, time = %s\n", typeI, typeJ, count, timer.timeString());
		}
	}

	// Delete temporary arrays
	for (int n=0; n<nTypes; ++n)
	{
		delete[] r[n];
		delete[] binss[n];
	}

	return true;
}

// Return PartialSet for specified Configuration (if it exists)
PartialSet* Partials::partialSet(Configuration* cfg)
{
	// Search existing list
	for (PartialSet* ps = partialSets_.first(); ps != NULL; ps = ps->next) if (ps->targetConfiguration() == cfg) return ps;

	return NULL;
}

// Calculate unweighted partials for the specified Configuration
bool Partials::calculateUnweighted(Configuration* cfg, ProcessPool& procPool, int method)
{
	// Does a PartialSet already exist for this Configuration?
	PartialSet* partialSet = Partials::partialSet(cfg);
	if (partialSet == NULL)
	{
		// No match, so create new
		partialSet = partialSets_.add();
		partialSet->setup(cfg, "unweighted", "rdf");
	}

	// Is the PartialSet already up-to-date?
	if (partialSet->upToDate())
	{
		Messenger::print("No need to calculate partials for Configuration '%s' - nothing has changed since the last calculation.\n", cfg->name());
		return true;
	}

	Messenger::print("Calculating partial RDFs for Configuration '%s'...\n", cfg->name());

	/*
	 * Reset any existing data
	 */

	partialSet->reset();

	/*
	 * Calculate full (intra+inter) partials
	 */

	// if (method == 0) ...
	Timer timer;
	calculateSimple(partialSet, procPool);
	timer.stop();
	Messenger::print("--> Finished calculation of partials (%s elapsed, %s comms).\n", timer.timeString(), procPool.accumulatedTimeString());

	/*
	 * Calculate intramolecular partials
	 */

	double distance;
	int start, stride;
	const Box* box = cfg->box();

	// Set start/skip for parallel loop (pool solo)
	start = procPool.interleavedLoopStart(ProcessPool::OverPoolProcesses);
	stride = procPool.interleavedLoopStride(ProcessPool::OverPoolProcesses);

	// Loop over molecules...
	Atom* i, *j, *k;
	for (int m=start; m<cfg->nMolecules(); m += stride)
	{
		Molecule* mol = cfg->molecule(m);

		// Bonds
		for (SpeciesBond* b = mol->species()->bonds(); b != NULL; b = b->next)
		{
			i = mol->atom(b->indexI());
			j = mol->atom(b->indexJ());
			if (cfg->useMim(i->cell(), j->cell())) distance = box->minimumDistance(i, j);
			else distance = (i->r() - j->r()).magnitude();
			partialSet->boundPartial(i->localTypeIndex(), j->localTypeIndex()).add(distance);
		}

		// Angles
		for (SpeciesAngle* a = mol->species()->angles(); a != NULL; a = a->next)
		{
			i = mol->atom(a->indexI());
			j = mol->atom(a->indexJ());
			k = mol->atom(a->indexK());
			
			// Determine whether we need to apply minimum image between 'j-i' and 'j-k'
			if (cfg->useMim(i->grain()->cell(), k->grain()->cell())) distance = box->minimumDistance(i, k);
			else distance = (i->r() - k->r()).magnitude();
			partialSet->boundPartial(i->localTypeIndex(), k->localTypeIndex()).add(distance);
		}
	}

	/*
	 * Partials are up-to-date
	 */

	partialSet->setUpToDate();

	return true;
}




/*
 * From Configuration::calcualtePartails()
 */

// 	
// 	// Calculate full partials
// 	Timer timer;
// 	if (!calculatePartialRDFs(procPool)) return false;
// 
// 	// Update partials index
// 	partialsIndex_ = coordinateIndex_;
// 
// 	// Collect all processes together
// 	if (!procPool.wait(ProcessPool::Pool)) return false;
// 
// 	timer.stop();

// 
// 	// Calculate intramolecular partials
// 	timer.start();
// 	if (!calculateIntramolecularRDFs(procPool)) return false;
// 	timer.stop();
// 	Messenger::print("--> Finished calculation of intramolecular partials (%s elapsed, %s comms).\n", timer.timeString(), procPool.accumulatedTimeString());
// 
// 	// Perform summation of partial data
// 	// Note that merging/summation of cross-term data (i.e. [n][m] with [m][n]) is not necessary since the partials matrix knows
// 	// that (i,j) == (j,i) as it is stored as a half-matrix in an Array2D object.
// 	int typeI, typeJ;
// 	procPool.resetAccumulatedTime();
// 	timer.start();
// 	double rho = atomicDensity();
// 	for (typeI=0; typeI<usedAtomTypes_.nItems(); ++typeI)
// 	{
// 		for (typeJ=typeI; typeJ<usedAtomTypes_.nItems(); ++typeJ)
// 		{
// 			// Sum histogram data from all processes
// 			if (!pairRDFMatrix_.ref(typeI,typeJ).allSum(procPool)) return false;
// 			if (!boundRDFMatrix_.ref(typeI,typeJ).allSum(procPool)) return false;
// 
// 			// Create unbound histogram from total and bound data
// 			unboundRDFMatrix_.ref(typeI, typeJ) = pairRDFMatrix_.ref(typeI,typeJ);
// 			unboundRDFMatrix_.ref(typeI, typeJ).addHistogramData(boundRDFMatrix_.ref(typeI,typeJ), -1.0);
// 
// 			// Finalise (normalise) partials
// 			pairRDFMatrix_.ref(typeI,typeJ).finalise();
// 			boundRDFMatrix_.ref(typeI,typeJ).finalise();
// 			unboundRDFMatrix_.ref(typeI,typeJ).finalise();
// 
// 			// Smooth partials if requested
// 			if (rdfSmoothing_ > 0)
// 			{
// 				pairRDFMatrix_.ref(typeI,typeJ).normalisedData().smooth(rdfSmoothing_*2+1);
// 				boundRDFMatrix_.ref(typeI,typeJ).normalisedData().smooth(rdfSmoothing_*2+1);
// 				unboundRDFMatrix_.ref(typeI,typeJ).normalisedData().smooth(rdfSmoothing_*2+1);
// 			}
// 
// 			// Copy RDF data ready for Fourier transform
// 			// -- Copy RDF data
// 			pairSQMatrix_.ref(typeI,typeJ) = pairRDFMatrix_.ref(typeI,typeJ).normalisedData();
// 			pairSQMatrix_.ref(typeI,typeJ).arrayY() -= 1.0;
// 			boundSQMatrix_.ref(typeI,typeJ) = boundRDFMatrix_.ref(typeI,typeJ).normalisedData();
// // 			boundSQMatrix_.ref(typeI,typeJ).arrayY() -= 1.0;
// 			unboundSQMatrix_.ref(typeI,typeJ) = unboundRDFMatrix_.ref(typeI,typeJ).normalisedData();
// 			unboundSQMatrix_.ref(typeI,typeJ).arrayY() -= 1.0;
// 		}
// 	}
// 	timer.stop();
// 	Messenger::print("--> Finished summation and normalisation of partial RDF data (%s elapsed, %s comms).\n", timer.timeString(), procPool.accumulatedTimeString());
// 	// Perform FT of partial g(r) into S(Q)
// 	// No instrumental broadening is applied in this case - the Configuration-based S(Q) are 'pure' in that sense
// 	// TODO Parallelise this
// 	procPool.resetAccumulatedTime();
// 	timer.start();
// 	for (typeI=0; typeI<usedAtomTypes_.nItems(); ++typeI)
// 	{
// 		for (typeJ=typeI; typeJ<usedAtomTypes_.nItems(); ++typeJ)
// 		{
// 			if (!pairSQMatrix_.ref(typeI,typeJ).transformBroadenedRDF(rho, qDelta, qMax, 0.0, 0.0, windowFunction)) return false;
// 			if (!boundSQMatrix_.ref(typeI,typeJ).transformBroadenedRDF(rho, qDelta, qMax, 0.0, 0.0, windowFunction)) return false;
// 			if (!unboundSQMatrix_.ref(typeI,typeJ).transformBroadenedRDF(rho, qDelta, qMax, 0.0, 0.0, windowFunction)) return false;
// 		}
// 	}
// 	timer.stop();
// 	Messenger::print("--> Finished Fourier transform of partial g(r) into partial S(Q) (%s elapsed, %s comms).\n", timer.timeString(), procPool.accumulatedTimeString());
// 
// 	// Calculate Bragg partials (if requested)
// 	if (braggOn)
// 	{
// 		procPool.resetAccumulatedTime();
// 		timer.start();
// 		if (!calculateBraggContributions(procPool)) return false;
// // 		if (!calculateBraggSQ(procPool)) return false;
// 		timer.stop();
// 		Messenger::print("--> Finished calculation of partial Bragg S(Q) (%s elapsed, %s comms).\n", timer.timeString(), procPool.accumulatedTimeString());
// 	}
// 
// 	// Generate final partial S(Q) combining pair correlations and Bragg partials
// 	for (typeI=0; typeI<usedAtomTypes_.nItems(); ++typeI)
// 	{
// 		for (typeJ=typeI; typeJ<usedAtomTypes_.nItems(); ++typeJ)
// 		{
// 			// Grab references:
// 			// -- Full atom pair-pair structure factors (bound + unbound)
// 			Data2D& pairSQ = pairSQMatrix_.ref(typeI,typeJ);
// 			// -- Bragg partial structure factors
// 			Data2D& braggSQ = braggSQMatrix_.ref(typeI,typeJ);
// 			// -- Full partial structure factors (destination)
// 			Data2D& partialSQ = partialSQMatrix_.ref(typeI,typeJ);
// 
// 			// Copy atomic pair S(Q) information
// 			partialSQ = pairSQ;
// 
// 			// Combine Bragg(Q) data if it was calculated
// 			if (braggOn)
// 			{
// 				double xRange = 0.1, x;
// 				double xMin = braggMaximumQ_ - xRange - qDelta*0.5;
// 				for (int n=0; n<braggSQ.nPoints(); ++n)
// 				{
// 					x = (braggSQ.x(n) <= xMin ? 0.0 : (braggSQ.x(n) - xMin) / xRange);
// 					partialSQ.setY(n, pairSQ.y(n)*x + braggSQ.y(n)*(1.0-x));
// 					// TEST - Straight replacement of original pair data
// // 					partialSQ.setY(n, braggSQ.y(n));
// 				}
// 			}
// 		}
// 	}
// 
// 	// Calculate total unweighted g(r) and F(Q)
// 	totalFQ_ = pairSQMatrix_.ref(0,0);
// 	totalFQ_.arrayY() = 0.0;
// 	double factor, braggMax;
// 	for (typeI=0; typeI<usedAtomTypes_.nItems(); ++typeI)
// 	{
// 		for (typeJ=typeI; typeJ<usedAtomTypes_.nItems(); ++typeJ)
// 		{
// 			factor = usedAtomTypes_[typeI]->fraction() * usedAtomTypes_[typeJ]->fraction() * (typeI == typeJ ? 1.0 : 2.0);
// 
// 			// Total RDF
// 			totalRDF_.addY(pairRDFMatrix_.ref(typeI,typeJ).normalisedData().arrayY(), factor);
// 			// TODO Does not include contributions from Bragg partials
// 
// 			// F(Q)
// 			totalFQ_.addY(partialSQMatrix_.ref(typeI,typeJ).arrayY(), factor);
// 		}
// 	}
// 	timer.stop();
// 	Messenger::print("--> Finished summation and FT of partials, and generation of total unweighted RDF/F(Q) (%s elapsed, %s comms).\n", timer.timeString(), procPool.accumulatedTimeString());
// 
// 	return true;
// }

