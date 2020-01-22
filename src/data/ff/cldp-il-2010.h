/*
        *** Canongia-Lopes, Deschamps & Padua Ionic Liquids (version 2010/09/16)
        *** src/data/ff/cldp-il-2010.h
	Copyright T. Youngs 2020

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

#ifndef DISSOLVE_FORCEFIELD_CLDP_IL_2010_H
#define DISSOLVE_FORCEFIELD_CLDP_IL_2010_H

#include "data/ff.h"

// Forward Declarations
/* none */

// Canongia-Lopes, Deschamps & Padua Ionic Liquids (2005) Forcefield
class Forcefield_CLDP2010 : public Forcefield
{
	public:
	// Constructor / Destructor
	Forcefield_CLDP2010();
	~Forcefield_CLDP2010();


	/*
	 * Definition
	 */
	public:
	// Return name of Forcefield
	const char* name() const;
	// Return description for Forcefield
	const char* description() const;
	// Return short-range interaction style for AtomTypes
	Forcefield::ShortRangeType shortRangeType() const;
};

#endif
