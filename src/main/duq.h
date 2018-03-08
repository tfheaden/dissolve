/*
	*** dUQ Main Structure
	*** src/main/duq.h
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

#ifndef DUQ_DUQ_H
#define DUQ_DUQ_H

#include "base/ptable.h"
#include "module/module.h"
#include "classes/configuration.h"
#include "classes/data.h"
#include "classes/pairpotential.h"
#include "classes/potentialmap.h"
#include "classes/masterintra.h"

// Forward Declarations
class Atom;
class Box;
class Cell;
class Grain;
class Isotopologue;
class Molecule;
class ChangeStore;
class ModuleProxyBase;

// dUQ Main Class
class DUQ
{
	public:
	// Constructor
	DUQ();
	// Destructor
	~DUQ();


	/*
	 * Core
	 */
	private:
	// Whether we are set up, ready for simulation
	bool setUp_;

	public:
	// Clear all data
	void clear();
	// Register GenericItems
	void registerGenericItems();
	// Set up everything needed to run the simulation
	bool setUp();
	// Return whether the simulation has been set up
	bool isSetUp() const;


	/*
	 * Atom Types
	 */
	private:
	// Master list of allowed AtomTypes for all Species (referenced from List<T>::masterInstance_)
	List<AtomType>& atomTypes_;
	
	public:
	// Add AtomType
	AtomType* addAtomType(int el);
	// Return number of AtomTypes in list
	int nAtomTypes() const;
	// Return first AtomType in list
	AtomType* atomTypes() const;
	// Return AtomTypes list
	const List<AtomType>& atomTypeList() const;
	// Return nth AtomType in list
	AtomType* atomType(int n);
	// Generate unique AtomType name with base name provided
	const char* uniqueAtomTypeName(const char* base, AtomType* exclude = 0) const;
	// Search for AtomType by name
	AtomType* findAtomType(const char* name) const;


	/*
	 * Master Terms for Species
	 */
	private:
	// List of master Bond parameters for Species
	List<MasterIntra> masterBonds_;
	// List of master Angles parameters for Species
	List<MasterIntra> masterAngles_;
	// List of master Torsions parameters for Species
	List<MasterIntra> masterTorsions_;

	public:
	// Add new master Bond parameters
	MasterIntra* addMasterBond(const char* name);
	// Return number of master Bond parameters in list
	int nMasterBonds() const;
	// Return list of master Bond parameters
	const List<MasterIntra>& masterBonds() const;
	// Return nth master Bond 
	MasterIntra* masterBond(int n);
	// Return whether named master Bond parameters exist
	MasterIntra* hasMasterBond(const char* name) const;
	// Add new master Angle parameters
	MasterIntra* addMasterAngle(const char* name);
	// Return number of master Angles parameters in list
	int nMasterAngles() const;
	// Return list of master Angle parameters
	const List<MasterIntra>& masterAngles() const;
	// Return nth master Angle parameters
	MasterIntra* masterAngle(int n);
	// Return whether named master Angle parameters exist
	MasterIntra* hasMasterAngle(const char* name) const;
	// Add new master Torsion parameters
	MasterIntra* addMasterTorsion(const char* name);
	// Return number of master Torsions parameters in list
	int nMasterTorsions() const;
	// Return list of master Torsion parameters
	const List<MasterIntra>& masterTorsions() const;
	// Return nth master Torsion parameters
	MasterIntra* masterTorsion(int n);
	// Return whether named master Torsion parameters exist
	MasterIntra* hasMasterTorsion(const char* name) const;


	/*
	 * Species Definitions
	 */
	private:
	// Master List of defined Species (referenced from List<T>::masterInstance_)
	List<Species>& species_;
	
	public:
	// Add a new Species to the list
	Species* addSpecies();
	// Remove the specified Species from the list
	void removeSpecies(Species* sp);
	// Return number of defined Species
	int nSpecies() const;
	// Return Species list
	List<Species>& species();
	// Return nth Species in the list
	Species* species(int n);
	// Search for Species by name
	Species* findSpecies(const char* name) const;
	// Update Species (or all) Isotopologues (or specified)
	void updateIsotopologues(Species* species = NULL, Isotopologue* iso = NULL);
	// Remove Isotopologue from Species
	void removeSpeciesIsotopologue(Species* species, Isotopologue* iso);


	/*
	 * Pair Potentials
	 */
	private:
	// Maximum distance for tabulated PairPotentials
	double pairPotentialRange_;
	// Maximum squared distance for tabulated PairPotentials
	double pairPotentialRangeSquared_;
	// Delta to use in tabulation
	double pairPotentialDelta_;
	// Whether Coulomb terms are included in generated PairPotentials
	bool pairPotentialsIncludeCoulomb_;
	// Simulation PairPotentials
	List<PairPotential> pairPotentials_;
	// Map for PairPotentials
	PotentialMap potentialMap_;

	public:
	// Set maximum distance for tabulated PairPotentials
	void setPairPotentialRange(double range);
	// Return maximum distance for tabulated PairPotentials
	double pairPotentialRange() const;
	// Return maximum squared distance for tabulated PairPotentials
	double pairPotentialRangeSquared() const;
	// Set delta to use in tabulations
	void setPairPotentialDelta(double delta);
	// Return delta to use in tabulations
	double pairPotentialDelta() const;
	// Set whether Coulomb term should be included in generated PairPotentials
	void setPairPotentialsIncludeCoulomb(bool b);
	// Return whether Coulomb term should be included in generated PairPotentials
	bool pairPotentialsIncludeCoulomb();
	// Return index of specified PairPotential
	int indexOf(PairPotential* pp);
	// Return number of defined PairPotentials
	int nPairPotentials() const;
	// Add new pair potential to list
	PairPotential* addPairPotential(AtomType* at1, AtomType* at2);
	// Return PairPotentials list
	const List<PairPotential>& pairPotentials() const;
	// Return nth PairPotential in list
	PairPotential* pairPotential(int n);
	// Return whether specified PairPotential is defined
	PairPotential* pairPotential(AtomType* at1, AtomType* at2) const;
	// Return whether specified PairPotential is defined
	PairPotential* pairPotential(const char* at1, const char* at2) const;
	// Save all PairPotentials
	bool savePairPotentials(const char* baseName) const;
	// Return map for PairPotentials
	const PotentialMap& potentialMap();
	// Regenerate all PairPotentials, replacing those currently defined
	void regeneratePairPotentials(PairPotential::ShortRangeType srType);
	// Update all currently-defined PairPotentials
	void updateCurrentPairPotentials();
	// Generate any missing PairPotentials using the supplied short-range form
	void generateMissingPairPotentials(PairPotential::ShortRangeType srType);


	/*
	 * Configurations
	 */
	private:
	// List of all atomic configurations
	List<Configuration> configurations_;

	public:
	// Add new Configuration
	Configuration* addConfiguration();
	// Return first Configuration in list
	const List<Configuration>& configurations() const;
	// Find configuration by name
	Configuration* findConfiguration(const char* name, bool useNiceName = false) const;
	// Write Configuration through specified LineParser
	bool writeConfiguration(Configuration* cfg, LineParser& parser);
	// Read Configuration through specified LineParser
	bool readConfiguration(Configuration* cfg, LineParser& parser);


	/*
	 * Data
	 */
	private:
	List<Data> data_;

	public:
	// Add new Data
	Data* addData();
	// Find Data by name
	Data* findData(const char* name, bool useNiceName = false) const;


	/*
	 * Simulation
	 */
	private:
	// List of Modules with pre-processing tasks to perform
	RefList<Module,bool> preProcessingTasks_;
	// List of main processing Modules to run
	ModuleList processingModules_;
	// Data associated with main processing Modules
	GenericList processingModuleData_;
	// List of Modules with post-processing tasks to perform
	RefList<Module,bool> postProcessingTasks_;
	// Current simulation step
	int iteration_;
	// Number of iterations performed
	int nIterationsPerformed_;
	// Main loop timer
	Timer mainLoopTimer_;

	private:
	// Find first occurrence of named Module in pre-processing tasks
	Module* findPreProcessingTask(const char* name);
	// Find first occurrence of named Module in post-processing tasks
	Module* findPostProcessingTask(const char* name);

	public:
	// Return list of main processing Modules to run
	ModuleList& processingModules();
	// Return data associated with main processing Modules
	GenericList& processingModuleData();
	// Iterate main simulation
	bool iterate(int nIterations = -1);
	// Return current simulation step
	int iteration() const;
	// Print timing information
	void printTiming();


	/*
	 * Simulation
	 */
	private:
	// Number of test points to use when calculating Box normalisation arrays
	int nBoxNormalisationPoints_;
	// Random seed
	int seed_;
	// Frequency at which to write restart file
	int restartFileFrequency_;

	public:
	// Set number of test points to use when calculating Box normalisation arrays
	void setNBoxNormalisationPoints(int nPoints);
	// Return number of test points to use when calculating Box normalisation arrays
	int nBoxNormalisationPoints() const;
	// Set random seed
	void setSeed(int seed);
	// Return random seed
	int seed() const;
	// Set frequency with which to write various iteration dat
	void setRestartFileFrequency(int n);
	// Return frequency with which to write restart file
	int restartFileFrequency() const;
	// Set up all simulation data, checking it as we go
	bool setUpSimulation();


	/*
	 * I/O
	 */
	private:
	// Filename of current input file
	CharString filename_;
	// Whether to automatically add dependent Modules if they have not been defined
	bool autoAddDependentModules_;
	// Accumulated timing information for saving restart file
	SampledDouble saveRestartTimes_;

	public:
	// Load datafiles
	bool loadDataFiles();
	// Load Species from specified file
	bool loadSpecies(const char* filename);
	// Load input file
	bool loadInput(const char* filename);
	// Save input file
	bool saveInput(const char* filename);
	// Load restart file
	bool loadRestart(const char* filename);
	// Save restart file
	bool saveRestart(const char* filename);
	// Return whether a filename has been set
	bool hasInputFileName() const;
	// Return filename of current input file
	const char* inputFilename() const;
	// Set whether to automatically add dependent Modules if they have not been defined
	void setAutoAddDependentModules(bool b);
	// Return whether to automatically add dependent Modules if they have not been defined
	bool autoAddDependentModules();


	/*
	 * Parallel Comms
	 */
	public:
	// Parallel Strategy
	enum ParallelStrategy { SequentialConfigStrategy, EvenStrategy, nParallelStrategies };
	// Convert string to ParallelStrategy
	static ParallelStrategy parallelStrategy(const char* s);

	private:
	// Parallel strategy for Configuration work
	ParallelStrategy parallelStrategy_;
	// Pool containing all available processes
	ProcessPool worldPool_;

	public:
	// Set parallel strategy for Configuration work
	void setParallelStrategy(ParallelStrategy ps);
	// Return parallel strategy for Configuration work
	ParallelStrategy parallelStrategy();
	// Return world process pool
	ProcessPool& worldPool();
	// Set up local MPI pools
	bool setUpMPIPools();
};

#endif
