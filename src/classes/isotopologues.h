// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2021 Team Dissolve and contributors

#pragma once

#include "classes/isotopologueweight.h"
#include "genericitems/base.h"
#include <vector>

// Forward Declarations
class Species;
class Isotopologue;
class ProcessPool;
class LineParser;

// Isotopologues
class Isotopologues : public GenericItemBase
{
    public:
    Isotopologues(const Species *species = nullptr, int speciesPopulation = 0);
    ~Isotopologues();

    /*
     * Isotopologue Mixture Data
     */
    private:
    // Associated Species
    const Species *species_;
    // Integer population of associated Species
    int speciesPopulation_;
    // Weighted Isotopologue mixture
    std::vector<IsotopologueWeight> mix_;

    public:
    // Set associated Species and population
    void setSpecies(const Species *sp, int population);
    // Return associated Species
    const Species *species() const;
    // Return associated Species population
    int speciesPopulation() const;
    /// Prune defunct Isotopologue entries
    void pruneMissing();
    // Add next available Isotopologue to list
    bool addNext();
    // Add specific Isotopologue to list
    void add(const Isotopologue *iso, double relativeWeight);
    // Set Isotopologue component in list
    bool set(const Isotopologue *iso, double relativeWeight);
    // Remove references to the specified Isotopologue
    void remove(const Isotopologue *iso);
    // Remove the specified IsotopologueWeight
    void remove(IsotopologueWeight *isoWeight);
    // Return whether the mix contains the specified Isotopologue
    bool contains(const Isotopologue *iso) const;
    // Return Isotopologue/weight mix
    std::vector<IsotopologueWeight> &mix();
    const std::vector<IsotopologueWeight> &mix() const;
    // Return number of Isotopologues in mix
    int nIsotopologues() const;
    // Return total relative population
    double totalRelative() const;
    // Normalise total relative population to 1.0
    void normalise();

    /*
     * GenericItemBase Implementations
     */
    public:
    // Return class name
    static std::string_view itemClassName();
    // Read data through specified LineParser
    bool read(LineParser &parser, CoreData &coreData);
    // Write data through specified LineParser
    bool write(LineParser &parser);

    /*
     * Parallel Comms
     */
    public:
    // Broadcast data
    bool broadcast(ProcessPool &procPool, const int root, const CoreData &coreData);
    // Check item equality
    bool equality(ProcessPool &procPool);
};
