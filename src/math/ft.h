/*
	*** Fourier Transforoms
	*** src/math/ft.h
	Copyright T. Youngs 2012-2020

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

#ifndef DISSOLVE_FOURIERTRANSFORM_H
#define DISSOLVE_FOURIERTRANSFORM_H

#include "math/broadeningfunction.h"
#include "math/windowfunction.h"

// Forward Declarations
class Data1D;

// Fourier Transforms
class Fourier
{
	public:
	// Perform Fourier sine transform of supplied data, over range specified, and with specified window and broadening functions applied
	static bool sineFT(Data1D& data, double normFactor, double wMin, double wStep, double wMax, WindowFunction windowFunction = WindowFunction(), BroadeningFunction broadening = BroadeningFunction());
};

#endif
