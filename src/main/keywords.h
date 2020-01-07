/*
	*** Keyword Definitions
	*** src/main/keywords.h
	Copyright T. Youngs 2012-2019

	This file is part of Dissolve.

	Dissolve is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	Dissolve is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with Dissolve.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DISSOLVE_KEYWORDS_H
#define DISSOLVE_KEYWORDS_H

#include "base/enumoptions.h"

// Forward Declarations
class LineParser;
class CoreData;
class Dissolve;
class Sample;
class Configuration;
class Species;
class Module;
class ModuleLayer;
class GenericList;
class Data;
class SpeciesInfo;
class SpeciesSite;

/*
 * Block Keywords
 */
namespace BlockKeywords
{
	// Block Keyword Enum
	enum BlockKeyword
	{
		ConfigurationBlockKeyword,		/* 'Configuration' - Defines a single Configuration for use in the simulation */
		LayerBlockKeyword,			/* 'Layer' - Defines a sequence of Modules in a processing layer */
		MasterBlockKeyword,			/* 'Master' - Contains master intramolecular terms for use in Species */
		ModuleBlockKeyword,			/* 'Module' - Sets up a Module to run after Configuration processing */
		PairPotentialsBlockKeyword,		/* 'PairPotentials' - Contains definitions of the PairPotentials for the simulation */
		SimulationBlockKeyword,			/* 'Simulation' - Setting of simulation variables affecting the calculation */
		SiteBlockKeyword,			/* 'Site' - Defines an analysis site within a Species */
		SpeciesBlockKeyword,			/* 'Species' - Begins a definition of a Species */
		nBlockKeywords				/* Number of defined BlockKeyword keywords */
	};
	// Return enum option info for BlockKeyword
	EnumOptions<BlockKeywords::BlockKeyword> keywords();
};


/*
 * Configuration Block Keywords
 */
namespace ConfigurationBlock
{
	// Configuration Block Keyword Enum
	enum ConfigurationKeyword
	{
		CellDivisionLengthKeyword,	/* 'CellDivisionLength' - Set the requested side length for regions when partitioning the unit cell */
		EndConfigurationKeyword,	/* 'EndConfiguration' - Signals the end of the Configuration block */
		GeneratorKeyword,		/* 'Generator' - Define the generator procedure for the Configuration */
		InputCoordinatesKeyword,	/* 'InputCoordinates' - Specifies the file which contains the starting coordinates */
		ModuleKeyword,			/* 'Module' - Starts the set up of a Module for this configuration */
		SizeFactorKeyword,		/* 'SizeFactor' - Scaling factor for Box lengths, Cell size, and Molecule centres-of-geometry */
		TemperatureKeyword,		/* 'Temperature' - Defines the temperature of the simulation */
		nConfigurationKeywords		/* Number of keywords defined for this block */
	};
	// Return enum option info for ConfigurationKeyword
	EnumOptions<ConfigurationBlock::ConfigurationKeyword> keywords();
	// Parse Configuration block
	bool parse(LineParser& parser, Dissolve* dissolve, Configuration* cfg);
};


/*
 * Layer Block Keywords
 */
namespace LayerBlock
{
	// Layer Block Keyword Enum
	enum LayerKeyword
	{
		DisabledKeyword,		/* 'Disabled' - Specify that the layer is currently disabled */
		EndLayerKeyword,		/* 'EndLayer' - Signals the end of the Layer block */
		FrequencyKeyword,		/* 'Frequency' - Frequency at which the layer is executed, relative to the main iteration counter */
		ModuleKeyword,			/* 'Module' - Begin a Module definition within this layer */
		nLayerKeywords			/* Number of keywords defined for this block */
	};
	// Return enum option info for LayerKeyword
	EnumOptions<LayerBlock::LayerKeyword> keywords();
	// Parse Layer block
	bool parse(LineParser& parser, Dissolve* dissolve, ModuleLayer* layer);
};


/*
 * Master Block Keywords
 */
namespace MasterBlock
{
	// Master Block Keyword Enum
	enum MasterKeyword
	{
		AngleKeyword,			/* 'Angle' - Define master Angle parameters that can be referred to */
		BondKeyword,			/* 'Bond' - Define master Bond parameters that can be referred to */
		EndMasterKeyword,		/* 'EndMaster' - Signals the end of the Master block */
		ImproperKeyword,		/* 'Improper' - Define master Improper parameters that can be referred to */
		TorsionKeyword,			/* 'Torsion' - Define master Torsion parameters that can be referred to */
		nMasterKeywords			/* Number of keywords defined for this block */
	};
	// Return enum option info for MasterKeyword
	EnumOptions<MasterBlock::MasterKeyword> keywords();
	// Parse Master block
	bool parse(LineParser& parser, CoreData& coreData);
};


/*
 * Module Block Keywords
 */
namespace ModuleBlock
{
	// Module Block Keyword Enum
	enum ModuleKeyword
	{
		ConfigurationKeyword,		/* 'Configuration' - Associates the specified Configuration to this Module */
		DisableKeyword,			/* 'Disable' - Disables the module, preventing it from running */
		EndModuleKeyword,		/* 'EndModule' - Signals the end of the Module block */
		FrequencyKeyword,		/* 'Frequency' - Frequency at which the Module is run */
		nModuleKeywords			/* Number of keywords defined for this block */
	};
	// Return enum option info for ModuleKeyword
	EnumOptions<ModuleBlock::ModuleKeyword> keywords();
	// Parse Module block
	bool parse(LineParser& parser, Dissolve* dissolve, Module* module, GenericList& targetList, bool moduleInConfiguration);
};


/*
 * PairPotential Block Keywords
 */
namespace PairPotentialsBlock
{
	// PairPotential Block Keyword Enum
	enum PairPotentialsKeyword
	{
		CoulombTruncationKeyword,		/* 'CoulombTruncation' - Truncation scheme to apply to Coulomb potential */
		DeltaKeyword,				/* 'Delta' - Gives the spacing between points in the tabulated potentials */
		EndPairPotentialsKeyword,		/* 'EndPairPotentials' - Signals the end of the PairPotentials block */
		GenerateKeyword,			/* 'Generate' - Generates a single PairPotential with the specified contributions */
		IncludeCoulombKeyword,			/* 'IncludeCoulomb' - Include Coulomb term in tabulated pair potentials" */
		ParametersKeyword,			/* 'Parameters' - Sets or re-sets the short-range and charge parameters for a specific AtomType */
		RangeKeyword,				/* 'Range' - Specifies the total range (inc. truncation width) over which to generate potentials */
		ShortRangeTruncationKeyword,		/* 'ShortRangeTruncation' - Truncation scheme to apply to short-range potential */
		ShortRangeTruncationWidthKeyword,	/* 'ShortRangeTruncationWidth' - Width of potential tail over which to reduce short-range term to zero */
		nPairPotentialsKeywords			/* Number of keywords defined for this block */
	};
	// Return enum option info for PairPotentialsKeyword
	EnumOptions<PairPotentialsBlock::PairPotentialsKeyword> keywords();
	// Parse PairPotentials block
	bool parse(LineParser& parser, Dissolve* dissolve);
};


/*
 * Simulation Block Keywords
 */
namespace SimulationBlock
{
	// Simulation Block Keyword Enum
	enum SimulationKeyword
	{
		EndSimulationKeyword,		/* 'EndSimulation' - Signals the end of the Simulation block */
		ParallelStrategyKeyword,	/* 'ParallelStrategy' - Determines the distribution of processes across Configurations */
		ParallelGroupPopulationKeyword,	/* 'ParallelGroupPopulation' - Controls the maximum number of groups to split processes in a pool in to */
		SeedKeyword,			/* 'Seed' - Random seed to use */
		nSimulationKeywords		/* Number of keywords defined for this block */
	};
	// Return enum option info for SimulationKeyword
	EnumOptions<SimulationBlock::SimulationKeyword> keywords();
	// Parse Simulation block
	bool parse(LineParser& parser, Dissolve* dissolve);
};

#endif
