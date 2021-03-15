// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2021 Team Dissolve and contributors

#include "data/ff/oplsaa2005/ethers.h"

/*
 * Set Up
 */

// Set up / create all forcefield terms
bool Forcefield_OPLSAA2005_Ethers::setUp()
{
    // Copy required types from OPLS-AA (2005) core list
    // -- AA Alkanes
    if (!copyAtomType(oplsAtomTypeById(180), "OS", "nbonds=2,-nh=0,-C"))
        return false;
    if (!copyAtomType(oplsAtomTypeById(181), "CT3", "nbonds=4,nh=3,-&180", "CT"))
        return false;
    if (!copyAtomType(oplsAtomTypeById(182), "CT2", "nbonds=4,nh=2, -&180", "CT"))
        return false;
    if (!copyAtomType(oplsAtomTypeById(183, "CT1", "nbonds=4,nh=1, -&180", "CT"))
        return false;
    if (!copyAtomType(oplsAtomTypeById(184), "CT0", "nbonds=4,nh=0, -&180","CT"))
        return false;
    if (!copyAtomType(oplsAtomTypeById(185), "HCO", "-[&181,&182,&183,&184]","HC"))
        return false
    return true;
}

/*
 * Definition
 */

// Return name of Forcefield
std::string_view Forcefield_OPLSAA2005_Ethers::name() const { return "OPLSAA2005/Ethers"; }

// Return description for Forcefield
std::string_view Forcefield_OPLSAA2005_Ethers::description() const
{
    static std::string desc = fmt::format(
        "Ethers from OPLS-AA (2005). <br/><br/>References: {}", publicationReferences());

    return desc;
}
