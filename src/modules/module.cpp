/*
	*** dUQ Module Interface
	*** src/modules/module.cpp
	Copyright T. Youngs 2012-2017

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

#include "modules/module.h"
#include "classes/configuration.h"
#include "base/lineparser.h"
#include "base/sysfunc.h"

// Static Members
int Module::staticLogPoint_ = 0;
int Module::staticBroadcastPoint_ = 0;

// Constructor
Module::Module()
{
	frequency_ = 1;
	enabled_ = true;
	configurationLocal_ = true;

	// LogPoints
	logPoint_ = 0;
	broadcastPoint_ = 0;
}

// Destructor
Module::~Module()
{
}

/*
 * Instances
 */

// Find instance with unique name specified
Module* Module::findInstance(const char* uniqueName)
{
	for (Module* instance = instances().first(); instance != NULL; instance = instance->next) if (DUQSys::sameString(instance->uniqueName(), uniqueName)) return instance;
	return NULL;
}

// Delete all instances of this Module
void Module::deleteInstances()
{
	while (instances().first()) instances().removeFirst();
}

/*
 * Definition
 */

// Set unique name of Module
void Module::setUniqueName(const char* uniqueName)
{
	uniqueName_ = uniqueName;
}

// Return unique name of Module
const char* Module::uniqueName()
{
	return uniqueName_.get();
}

// Add dependent Module to this Module
void Module::addDependentModule(Module* module, bool autoAdded)
{
	dependentModules_.add(module, autoAdded);
}

// Return pointer for specified dependent Module
Module* Module::dependentModule(const char* name)
{
	RefListIterator<Module,bool> iterator(dependentModules_, true);
	while (Module* module = iterator.iterate()) if (DUQSys::sameString(name, module->name())) return module;

	return NULL;
}

// Update targets for any auto-added dependent Modules with those of this Module
void Module::updateDependentTargets()
{
	RefListIterator<Module,bool> iterator(dependentModules_);
	while (Module* module = iterator.iterate())
	{
		// If the Module was *not* added automatically, then do not update the targets since it is the user's responsibility to ensure consistency
		if (!iterator.currentData()) continue;

		// Copy target Configurations
		module->copyTargetConfigurations(this);

		// Set frequency of module to match that of the parent
		module->setFrequency(frequency_);

		// Perform any set up on the dependent module
		setupDependentModule(module);

		// Also may need to update dependent targets of this Module...
		module->updateDependentTargets();
	}
}

/*
 * Keywords
 */

// Parse keyword line, returning true (1) on success, false (0) for recognised but failed, and -1 for not recognised
int Module::parseKeyword(LineParser& parser, DUQ* duq, GenericList& targetList, const char* prefix)
{
	// The LineParser currently contains a parsed line from the input file

	// Do we recognise the first item (the 'keyword')?
	ModuleKeywordBase* keyword = keywords_.find(parser.argc(0));
	if (!keyword) return -1;

	// We recognised the keyword - what should we try to do with it?
	if (keyword->type() == ModuleKeywordBase::ComplexData)
	{
		// It's a 'complex' keyword, one that either sets up a complicated object, or does something specific within the Module
		return parseComplexKeyword(keyword, parser, duq, targetList, prefix);
	}
	else
	{
		// It's a keyword that sets a simple POD or class
		// Check the number of arguments we have against the min / max for the keyword
		if ((parser.nArgs() - 1) < keyword->minArguments())
		{
			Messenger::error("Not enough arguments given to Module keyword '%s'.\n", keyword->keyword());
			return 0;
		}
		if ((keyword->maxArguments() >= 0) && ((parser.nArgs() - 1) > keyword->maxArguments()))
		{
			Messenger::error("Too many arguments given to Module keyword '%s'.\n", keyword->keyword());
			return 0;
		}

		// All OK, so pass the current line to it for parsing
		if (!keyword->parseArguments(parser, 1))
		{
			Messenger::error("Failed to parse arguments for Module keyword '%s'.\n", keyword->keyword());
			return 0;
		}

		// Add a copy of the data to the target list
		keyword->duplicateInList(targetList, prefix);
	}

	return true;
}

/*
 * Control
 */

// Frequency with which to run Module (relative to master simulation loop counter)
void Module::setFrequency(int freq)
{
	frequency_ = freq;
}

// Frequency with which to run Module (relative to master simulation loop counter)
int Module::frequency()
{
	return frequency_;
}

// Return whether the Module should run this iteration
bool Module::runThisIteration(int iteration)
{
	if ((frequency_ < 1) || (!enabled_)) return false;
	else if ((iteration%frequency_) == 0) return true;
	else return false;
}

// Return short descriptive text relating frequency to supplied iteration number
const char* Module::frequencyDetails(int iteration)
{
	static CharString result;

	if (frequency_ < 0) return "NEGATIVE?";
	else if (frequency_ == 0) return "disabled";
	else if (frequency_ == 1) return "every time";
	else if ((iteration%frequency_) == 0) return "this iteration";
	else
	{
		// Calculate number of steps necessary to get to next multiple of the frequency_
		int nToGo = frequency_ - (iteration - frequency_*(iteration/frequency_));
		if (nToGo == 1) return "next iteration";

		result.sprintf("in %i steps time", nToGo);
		return result.get();
	}
}

// Set whether the Module is enabled
void Module::setEnabled(bool b)
{
	enabled_ = b;
}

// Return whether the Module is enabled
bool Module::enabled()
{
	return enabled_;
}

/*
 * Targets
 */

// Add Configuration target
bool Module::addConfigurationTarget(Configuration* cfg)
{
	// Check how many Configurations we accept before we do anything else
	if ((nTargetableConfigurations() == -1) || (targetConfigurations_.nItems() < nTargetableConfigurations()))
	{
		targetConfigurations_.add(cfg);
		return true;
	}
	else
	{
		if (nTargetableConfigurations() == 0) Messenger::error("Can't add Configuration '%s' as a target to Module '%s' since it doesn't accept any such targets.\n", cfg->name(), name());
		else Messenger::error("Can't add Configuration '%s' as a target to Module '%s' since the maximum number (%i) has already been reached.\n", cfg->name(), name(), nTargetableConfigurations());
	}

	return false;
}

// Return number of targeted Configurations
int Module::nConfigurationTargets()
{
	return targetConfigurations_.nItems();
}

// Return first targeted Configuration
RefList<Configuration,bool>& Module::targetConfigurations()
{
	return targetConfigurations_;
}

// Return if the specified Configuration is in the targets list
bool Module::isTargetConfiguration(Configuration* cfg)
{
	return targetConfigurations_.contains(cfg);
}

// Copy Configuration targets from specified Module
void Module::copyTargetConfigurations(Module* sourceModule)
{
	// First, check if this module actually accepts target Configurations
	if ((nTargetableConfigurations() < sourceModule->nConfigurationTargets()) && (nTargetableConfigurations() != -1))
	{
		Messenger::warn("Dependent Module '%s' does not accept Configuration targets, but the source Module '%s' lists %i.\n", name(), sourceModule->name());
		return;
	}
	RefListIterator<Configuration,bool> configIterator(sourceModule->targetConfigurations());
	while (Configuration* cfg = configIterator.iterate()) addConfigurationTarget(cfg);
}

// Set whether this module is a local Module in a Configuration
void Module::setConfigurationLocal(bool b)
{
	configurationLocal_ = b;
}

// Return whether this module is a local Module in a Configuration
bool Module::configurationLocal()
{
	return configurationLocal_;
}

// Add Data target
bool Module::addDataTarget(Data* data)
{
	targetData_.add(data);
	return true;
}

// Return number of targeted Data
int Module::nDataTargets()
{
	return targetData_.nItems();
}

// Return first targeted Data
RefList<Data,bool>& Module::targetData()
{
	return targetData_;
}

/*
 * Data Reaping
 */

// Reap data for all instances of this Module
void Module::reapAllInstances()
{
	for (Module* instance = instances().first(); instance != NULL; instance = instance->next) instance->reapData();
}
