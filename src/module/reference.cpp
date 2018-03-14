/*
	*** Module Reference
	*** src/module/reference.cpp
	Copyright T. Youngs 2012-2018

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

#include "module/reference.h"

// Constructor
ModuleReference::ModuleReference()
{
	module_ = NULL;
	configuration_ = NULL;
}

// Destructor
ModuleReference::~ModuleReference()
{
}

// Conversion to Module*
ModuleReference::operator Module*()
{
	return module_;
}

/*
 * Data
 */

// Set Module and location
void ModuleReference::set(Module* module, Configuration* cfg)
{
	module_ = module;
	configuration_ = cfg;
}

// Return referenced Module
Module* ModuleReference::module()
{
	return module_;
}

// Return referenced Configuration
Configuration* ModuleReference::configuration()
{
	return configuration_;
}

// Return whether the Module is a processing Module (i.e. is not in a Configuration)
bool ModuleReference::isProcessingModule() const
{
	return (configuration_ == NULL);
}