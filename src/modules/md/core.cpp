// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Team Dissolve and contributors

#include "modules/md/md.h"

/*
 * Constructor / Destructor
 */

MDModule::MDModule() : Module(nRequiredTargets())
{
    // Initialise Module - set up keywords etc.
    initialise();
}

MDModule::~MDModule() {}

/*
 * Instances
 */

// Create instance of this module
Module *MDModule::createInstance() const { return new MDModule; }
