/*
	*** AtomShake Module
	*** src/modules/atomshake.cpp
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
#include "modules/atomshake.h"
#include "classes/box.h"
#include "classes/cell.h"
#include "classes/changestore.h"
#include "classes/energykernel.h"
#include "base/comms.h"
#include "base/timer.h"
#include "math/matrix4.h"

/*
 * Constructor / Destructor
 */

// Constructor
AtomShake::AtomShake() : Module()
{
	// Setup variables / control parameters
	setVariable("cutoffDistance", -1.0);
	setVariable("nShakesPerAtom", 1);
	setVariable("targetAcceptanceRate", 0.33);
	setVariable("stepSize", 0.05);
}

// Destructor
AtomShake::~AtomShake()
{
}

// Create instance of this module
Module* AtomShake::createInstance()
{
	return new AtomShake;
}

/*
 * Definition
 */

// Return name of module
const char* AtomShake::name()
{
	return "AtomShake";
}

// Return brief description of module
const char* AtomShake::brief()
{
	return "Perform atomic Monte Carlo on all atoms";
}

// Return type of module
Module::ModuleType AtomShake::type()
{
	return Module::EvolutionModule;
}

// Number of Configurations that this module requires to run
int AtomShake::nConfigurationsRequired()
{
	return 1;
}

/*
 * Method
 */

// Execute Method
bool AtomShake::execute(DUQ& duq)
{
	/*
	* Perform an Atom shake
	* 
	* This is a parallel routine, with processes operating in groups.
	*/

	// Get target Configuration
	Configuration* cfg = targetConfigurations_.firstItem();

	// Retrieve control parameters from Configuration
	const double cutoffDistance = variableAsDouble("cutoffDistance") < 0.0 ? duq.pairPotentialRange() : variableAsDouble("cutoffDistance");
	const int nShakesPerAtom = variableAsInt("nShakesPerAtom");
	const double targetAcceptanceRate = variableAsDouble("targetAcceptanceRate");
	double stepSize = variableAsDouble("stepSize");
	const double termScale = 1.0;
	const double rRT = 1.0/(.008314472*cfg->temperature());

	// Print argument/parameter summary
	Messenger::print("AtomShake: Cutoff distance is %f\n", cutoffDistance);
	Messenger::print("AtomShake: Performing %i shake(s) per Atom\n", nShakesPerAtom);
	Messenger::print("AtomShake: Translation step is %f Angstroms, target acceptance rate is %f.\n", stepSize, targetAcceptanceRate);

	// Initialise the Cell distributor
	const bool willBeModified = true, allowRepeats = false;
	cfg->initialiseCellDistribution();

	// Create a local ChangeStore and EnergyKernel
	ChangeStore changeStore;
	EnergyKernel kernel(cfg, duq.potentialMap(), cutoffDistance);

	// Initialise the random number buffer
	Comm.initialiseRandomBuffer(DUQComm::Group);

	// Enter calculation loop until no more Cells are available
	int cellId, shake, n, nbr;
	int nTries = 0, nAccepted = 0;
	bool accept;
	double currentEnergy, intraEnergy, newEnergy, newIntraEnergy, delta, totalDelta = 0.0;
	Cell* cell;
	Grain* grainI;
	Vec3<double> centre, rDelta;
	Atom** cellAtoms;

	Timer timer;
	Comm.resetAccumulatedTime();
	while (cellId = cfg->nextAvailableCell(willBeModified, allowRepeats), cellId != Cell::AllCellsComplete)
	{
		// Check for valid cell
		if (cellId == Cell::NoCellsAvailable)
		{
			// No valid cell, but still need to enter into change distribution with other processes
			changeStore.distributeAndApply(cfg);
			cfg->finishedWithCell(willBeModified, cellId);
			continue;
		}
		cell = cfg->cell(cellId);
		Messenger::printVerbose("Cell %i now the target, containing %i Grains interacting with %i neighbours.\n", cellId, cell->nGrains(), cell->nTotalCellNeighbours());

		/*
		 * Calculation Begins
		 */

		// Set current atom targets in ChangeStore (entire cell contents)
		changeStore.add(cell);

		// Loop over atoms in this cell
		cellAtoms = cell->atoms().objects();
		for (n = 0; n < cell->atoms().nItems(); ++n)
		{
			// Check for unused atom
			Atom* i = cellAtoms[n];

			// Get the atom's grain pointer
			grainI = i->grain();

			// Calculate reference intramolecular energy for atom, including intramolecular terms through the atom's grain
			currentEnergy = kernel.energy(i, DUQComm::Group);
			intraEnergy = kernel.fullIntraEnergy(grainI, termScale);

			// Loop over number of shakes per atom
			for (shake=0; shake<nShakesPerAtom; ++shake)
			{
				// Create a random translation vector
				rDelta.set(Comm.randomPlusMinusOne()*stepSize, Comm.randomPlusMinusOne()*stepSize, Comm.randomPlusMinusOne()*stepSize);

				// Translate atom and calculate new energy
				i->translateCoordinates(rDelta);
				newEnergy = kernel.energy(i, DUQComm::Group);
				newIntraEnergy = kernel.fullIntraEnergy(grainI, termScale);
				
				// Trial the transformed atom position
				delta = (newEnergy + newIntraEnergy) - (currentEnergy + intraEnergy);
// 				printf("delta = %f\n", delta);
				accept = delta < 0 ? true : (Comm.random() < exp(-delta*rRT));

				if (accept)
				{
// 					Messenger::print("Accepts move with delta %f\n", delta);
					// Accept new (current) position of target Atom
					changeStore.updateAtom(n);
					currentEnergy = newEnergy;
					intraEnergy = newIntraEnergy;
					totalDelta += delta;
					++nAccepted;
				}
				else changeStore.revert(n);
				
				++nTries;
			}

		}

		// Store modifications to Atom positions ready for broadcast later
		changeStore.storeAndReset();

		/*
		 * Calculation End
		 */

		// Distribute coordinate changes to all processes
		changeStore.distributeAndApply(cfg);
		changeStore.reset();

		// Must unlock the Cell when we are done with it!
		cfg->finishedWithCell(willBeModified, cellId);
	}
	timer.stop();

	// Grains have moved, so refold and update locations
	cfg->updateGrains();

	// Collect statistics from process group leaders
	if (!Comm.allSum(&nAccepted, 1, DUQComm::Leaders)) return false;
	if (!Comm.allSum(&nTries, 1, DUQComm::Leaders)) return false;
	if (!Comm.allSum(&totalDelta, 1, DUQComm::Leaders)) return false;
	if (Comm.processGroupLeader())
	{
		double rate = double(nAccepted)/nTries;

		Messenger::print("AtomShake: Overall acceptance rate was %4.2f% (%i of %i attempted moves) (%s work, %s comms)\n", 100.0*rate, nAccepted, nTries, timer.timeString(), Comm.accumulatedTimeString());
		Messenger::print("AtomShake: Total energy delta was %10.4e kJ/mol.\n", totalDelta);

		// Adjust step size
		stepSize *= rate/targetAcceptanceRate;
// 		if (stepSize_ < 0.05) stepSize_ = 0.05;
// 		else if (stepSize_ > maxTranslationStep_) stepSize_ = maxTranslationStep_;
// 		if (rotationStep_ < 3.0) rotationStep_ = 3.0;
	}

	// Store updated parameter values
	if (!Comm.broadcast(&stepSize, 1, 0, DUQComm::Group)) return false;
// 	stepSizeParam->setValue(stepSize); TODO
	Messenger::print("AtomShake: Updated translation step is %f Angstroms.\n", stepSize);
	
	// Increment configuration changeCount_
	if (nAccepted > 0) cfg->incrementCoordinateIndex();

	// Update total energy
	cfg->registerEnergyChange(totalDelta);
	cfg->accumulateEnergyChange();

	return true;
}

// // Perform a world atom shake
// bool DUQ::worldAtomShake(Configuration& cfg, double cutoffDistance, int nShakes, double targetAcceptanceRate, double stepSize)
// {
// 	// Control Parameters
// 	const double termScale = 1.0;	
// 	if (cutoffDistance < 0.0) cutoffDistance = pairPotentialRange_;
// 	const double rRT = 1.0/(.008314472*cfg->temperature());
// 	
// 	// Print argument/parameter summary
// 	Messenger::print("WorldAtomShake: Cutoff distance is %f\n", cutoffDistance);
// 	Messenger::print("WorldAtomShake: Translation step is %f Angstroms, target acceptance rate is %f.\n", stepSize, targetAcceptanceRate);
// 
// 	// Initialise the Cell distributor
// 	const bool willBeModified = true, allowRepeats = false;
// 	cfg->initialiseCellDistribution();
// 
// 	// Create a local ChangeStore and EnergyKernel
// 	ChangeStore changeStore;
// 	EnergyKernel kernel(cfg, potentialMap_, cutoffDistance);
// 
// 	// Initialise the random number buffer
// 	Comm.initialiseRandomBuffer(DUQComm::World);
// 
// 	// Enter calculation loop until no more Cells are available
// 	int cellId, shake, n, i;
// 	int nTries = 0, nAccepted = 0;
// 	bool accept;
// 	double currentEnergy, newEnergy, delta, totalDelta = 0.0;
// 	Cell* cell;
// 	Vec3<double> rDelta;
// 	int nAtoms = cfg->nAtoms();
// 	Atom* atoms = cfg->atoms();
// 
// 	Timer timer;
// 	Comm.resetAccumulatedTime();
// 
// 	// Add entire cell contents to a ChangeStore
// 	for (n=0; n<nAtoms; ++n) changeStore.add(&atoms[n]);
// 
// 	// Calculate reference energy
// 	currentEnergy = interatomicEnergy(cfg) + intramolecularEnergy(cfg);
// 
// 	// Loop over nShakes
// 	for (int n=0; n<nShakes; ++n)
// 	{
// 		// Randomly displace each atom in the system
// 		for (i=0; i<nAtoms; ++i)
// 		{
// 			// Create a random translation vector
// 			rDelta.set(Comm.randomPlusMinusOne()*stepSize, Comm.randomPlusMinusOne()*stepSize, Comm.randomPlusMinusOne()*stepSize);
// 
// 			// Translate atom
// 			atoms[i].translateCoordinates(rDelta);
// 		}
// 
// 		// Update neighbour lists
// 		cfg->recreateCellAtomNeighbourLists(pairPotentialRange_);
// 
// 		// Calculate new energy
// 		newEnergy = interatomicEnergy(cfg) + intramolecularEnergy(cfg);
// 
// 		delta = newEnergy - currentEnergy;
// 		printf("delta = %f\n", delta);
// 		accept = delta < 0 ? true : (Comm.random() < exp(-delta*rRT));
// 
// 		if (accept)
// 		{
// // 			Messenger::print("Accepts move with delta %f\n", delta);
// 			// Accept new (current) positions of atoms
// 			changeStore.updateAll();
// 			currentEnergy = newEnergy;
// 			totalDelta += delta;
// 			++nAccepted;
// 		}
// 		else
// 		{
// 			changeStore.revertAll();
// 			cfg->recreateCellAtomNeighbourLists(pairPotentialRange_);
// 		}
// 	}
// 
// 	// Collect statistics from process group leaders
// 	if (!Comm.allSum(&nAccepted, 1, DUQComm::Leaders)) return false;
// 	if (!Comm.allSum(&nTries, 1, DUQComm::Leaders)) return false;
// 	if (!Comm.allSum(&totalDelta, 1, DUQComm::Leaders)) return false;
// 	if (Comm.processGroupLeader())
// 	{
// 		double rate = double(nAccepted)/nTries;
// 
// 		Messenger::print("WorldAtomShake: Overall acceptance rate was %4.2f% (%i of %i attempted moves) (%s work, %s comms)\n", 100.0*rate, nAccepted, nTries, timer.timeString(), Comm.accumulatedTimeString());
// 		Messenger::print("WorldAtomShake: Total energy delta was %10.4e kJ/mol.\n", totalDelta);
// 
// 		// Adjust step size
// 		stepSize *= rate/targetAcceptanceRate;
// // 		if (stepSize_ < 0.05) stepSize_ = 0.05;
// // 		else if (stepSize_ > maxTranslationStep_) stepSize_ = maxTranslationStep_;
// // 		if (rotationStep_ < 3.0) rotationStep_ = 3.0;
// 	}
// 
// 	// Store updated parameter values
// 	if (!Comm.broadcast(&stepSize, 1, 0, DUQComm::Group)) return false;
// // 	stepSizeParam->setValue(stepSize); TODO
// 	Messenger::print("WorldAtomShake: Updated translation step is %f Angstroms.\n", stepSize);
// 	
// 	// Increment configuration changeCount_
// 	if (nAccepted > 0) cfg->incrementCoordinateIndex();
// 
// 	// Update total energy
// 	cfg->registerEnergyChange(totalDelta);
// 	cfg->accumulateEnergyChange();
// }

