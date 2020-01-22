/*
        *** Canongia-Lopes, Deschamps & Padua Ionic Liquids (version 2010/09/16)
        *** src/data/ff/cldp-il-2010.cpp
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

#include "data/ff/cldp-il-2010.h"
#include "data/ffangleterm.h"
#include "data/ffatomtype.h"
#include "data/ffbondterm.h"
#include "data/ffparameters.h"
#include "data/fftorsionterm.h"
#include "data/ffimproperterm.h"
#include "classes/atomtype.h"
#include "classes/speciesatom.h"
#include "base/sysfunc.h"

/*
 * Implements "Molecular force field for ionic liquids (version 2010/09/16)"
 * This forcefield is the work of K. Shimizu,	A. Podgorsek,	F. Hammami,	L. Gontrani,	J. Canongia Lopes,	and A. Padua.
 * Questions regarding the forcefield itself to: agilio.padua@univ-bpclermont.fr
 * Latest version donwloadable at: http://therm10.univ-bpclermont.fr/~apadua/
 * Reformatted for use in Dissolve
 * Notes:
 * Any inconsistencies between the forcefield as implemented here and the original work are the sole responsibility of TGAY.
 * All energy values are in kJ/mol.
 */

// Constructor / Destructor
Forcefield_CLDP2010::Forcefield_CLDP2010()
{
	static ForcefieldParameters parameters[] =
	{
		//           Name    Form(ifnot default)      Param1(Epsilon)     Param2(Sigma)
			// dialkylimidazolium  JPCB108:2038(2004)
			{this,	"CT",	0.27614,	3.50 },
			{this,	"CR",	0.29288,	3.55 },
			{this,	"CW",	0.29288,	3.55 },
			{this,	"HA",	0.12552,	2.42 },
			{this,	"HC",	0.12552,	2.50 },
			{this,	"H1",	0.12552,	2.50 },
			{this,	"NA",	0.71128,	3.25 },
			// alkenylimidazolium  OPLS-AA:unpublished
			{this,	"CM",	0.31798,	3.55 },
			{this,	"HCM",	0.12552,	2.42 },
			{this,	"HAL",	0.12552,	2.42 },
			// alcohols  OPLS-AA:JACS118:11225(1996);JPC100:18010(1996)
			{this,	"H1O",	0.12552,	2.50 },
			{this,	"OH",	0.71128,	3.21 },
			{this,	"HO",	0.00000,	0.00 },
			// esters  OPLS-AA:JCC22:1340(2001)
			{this,	"CC",	0.43932,	3.75 },
			{this,	"OC",	0.87864,	2.96 },
			{this,	"OS",	0.71176,	3.00 },
			{this,	"HCO",	0.06276,	2.42 },
			// ether  OPLS-AA
			{this,	"OE",	0.87864,	2.96 },
			// pyridinium  JCSPerkin2:2365(1999)
			{this,	"HAP",	0.12552,	2.42 },
			// ammonium(also_pyrrolidinium)  OPLS-AA:JACS121:4827(1999),AMBER
			{this,	"N3",	0.71128,	3.25 },
			{this,	"H",	0.12552,	2.50 },
			// phosphonium  OPLS-AA
			{this,	"P3",	0.8368,	3.74 },
			// tetrafluoroborate_anion  ""
			{this,	"B",	0.3975,	3.58 },
			{this,	"FB",	0.2552,	3.12 },
			// hexafluorophosphate  JCSPerkin2:2365(1999)
			{this,	"P",	0.8368,	3.74 },
			{this,	"F",	0.2552,	3.12 },
			// nitrate  JPCB108:2038(2004)
			{this,	"ON",	0.3380,	3.06 },
			{this,	"NN",	0.6100,	2.77 },
			// chloride  JPCB108:2038(2004)
			{this,	"Cl",	0.8300,	3.65 },
			// bromide  JPCB110:19586(2006)
			{this,	"Br",	0.86,	3.97 },
			// bistriflylamide  JPCB108:16893(2004)
			{this,	"SBT",	1.04600,	3.55 },
			{this,	"OBT",	0.87864,	2.96 },
			{this,	"NBT",	0.71128,	3.25 },
			// perfluoroalkanes  OPLS-AA
			{this,	"CSF",	0.27614,	3.50 },
			{this,	"F",	0.22175,	2.95 },
			// dicyanamide  JPCB110:19586(2006)
			{this,	"N3A",	0.71128,	3.25 },
			{this,	"CZA",	0.27614,	3.30 },
			{this,	"NZA",	0.71128,	3.20 },
			// alkylsulfates  JPCB112:5039(2008)
			{this,	"OS4",	0.83700,	3.15 },
			{this,	"OC4",	0.58600,	2.90 },
			{this,	"HS4",	0.12600,	2.50 },
			// alkylsulfonates  JPCB112:5039(2008)
			{this,	"OS3",	0.83700,	3.15 },
			{this,	"HS3",	0.12600,	2.50 },
			// acetate  OPLS-AA
			{this,	"O2",	0.87860,	2.96 },
			{this,	"CO2",	0.43930,	3.75 }
	};

	static ForcefieldAtomType atomTypes[] =
	{
		//      Z     El          FFID    Name      Type               Description
		//                                                q      Epsilon   Sigma
			// dialkylimidazolium  JPCB108:2038(2004)
		{ this,	ELEMENT_C,	1,	"C1",	"nbonds=4,-N,-H,-H",	"First tetrahedral carbon in R group on ring nitrogen",
								-0.17,	"CT" },
		{ this,	ELEMENT_C,	2,	"C2",	"nbonds=4,-C(ring(none),-[N,P]),-H(n=2),-C",	"Second carbon in R group (> ethyl)",
								0.01,	"CT" },
		{ this,	ELEMENT_C,	3,	"CE",	"nbonds=4,-C(ring(none),-[N,P]),-H(n=3)",	"Terminal carbon of ethyl R group",
								-0.05,	"CT" },
		{ this,	ELEMENT_C,	4,	"CS",	"nbonds=4,-C(n=2),-H(n=2)",	"General carbon in R group",
								-0.12,	"CT" },
		{ this,	ELEMENT_C,	5,	"CT",	"nbonds=4,-C(-C(-H(n=2)),-C),-H(n=3)",	"Terminal carbon in R group (> ethyl)",
								-0.18,	0.27614,	3.50 },
		{ this,	ELEMENT_C,	6,	"CR",	"nbonds=3,-H,-&13,-[&13,&14]",	"Carbon adjacent to two nitrogens with alkyl groups",
								-0.11,	0.29288,	3.55 },
		{ this,	ELEMENT_C,	7,	"CW",	"nbonds=3,-H,-C,-[&13,&14,&19]",	"Carbon adjacent to nitrogen (with alkyl group) and carbon",
								-0.13,	0.29288,	3.55 },
		{ this,	ELEMENT_H,	9,	"HCR",	"-&6",	"H on carbon adjacent to two N",
								0.21,	"HA" },
		{ this,	ELEMENT_H,	10,	"HCW",	"-&7",	"H on carbon adjacent to N and C",
								0.21,	"HA" },
		{ this,	ELEMENT_H,	11,	"HC",	"",	"Hydrogen on other carbon in R group",
								0.06,	0.12552,	2.50 },
		{ this,	ELEMENT_H,	12,	"H1",	"-C(-[N,P,S])",	"Hydrogen on first carbon in R group",
								0.13,	"HC" },
		{ this,	ELEMENT_N,	13,	"NA",	"nbonds=3,ring(size=5,-N,-C,-N,-C,-C),-C(nh=1,n=2),-C",	"Nitrogen in dialkylimidazolium with alkyl group attached",
								0.15,	0.71128,	3.25 },
		// alkylimidazolium  JPCB110:19586(2006)
		{ this,	ELEMENT_N,	14,	"NAH",	"nbonds=3,ring(size=5,-N,-C,-N,-C,-C),-C(nh=1,n=2),-H",	"Nitrogen with hydrogen attached",
								-0.21,	"NA" },
		{ this,	ELEMENT_C,	15,	"CRH",	"nbonds=3,-&13,-&21,-H",	"Carbon adjacent to two nitrogens",
								0.00,	"CR" },
		{ this,	ELEMENT_C,	16,	"CWH",	"nbonds=3,-H,-C,-&21",	"Carbon adjacent to N with H attached",
								-0.03,	"CW" },
		{ this,	ELEMENT_H,	17,	"HAN",	"-&14",	"Hydrogen on ring nitrogen",
								0.37,	0.00000,	0.00 },
		// dialkylmethylmidazolium  JPCB112:5039(2008)
		{ this,	ELEMENT_C,	18,	"CRM",	"nbonds=3,-&19(n=2)",	"Methylated carbon adjacent to two nitrogens in 5-membered ring",
								0.19,	"CR" },
		{ this,	ELEMENT_N,	19,	"NAM",	"ring(size=5,-N,-C,-N,-C,-C),-C(nh=1),-C(ring(none)),-C(-C(ring(none),nh=3))",	"Ring nitrogen in dialkylmethylimidazolium",
								0.04,	"NA" },
		{ this,	ELEMENT_C,	20,	"CCR",	"-&18",	"Methyl carbon at C2 position",
								-0.26,	"CT" },
		// alkenylimidazolium  OPLS-AA:unpublished
		{ this,	ELEMENT_C,	21,	"CM",	"",	"",
								-0.115,	0.31798,	3.55 },
		{ this,	ELEMENT_C,	22,	"CMT",	"",	"",
								-0.230,	"CM" },
		{ this,	ELEMENT_H,	23,	"HCM",	"",	"",
								0.115,	0.12552,	2.42 },
		{ this,	ELEMENT_C,	24,	"CAL",	"",	"",
								0.170,	"CM" },
		{ this,	ELEMENT_H,	25,	"HAL",	"",	"",
								-0.020,	"HC" },
		// alcohols  OPLS-AA:JACS118:11225(1996);JPC100:18010(1996)
		{ this,	ELEMENT_C,	26,	"CTO",	"-O(-H)",	"Alcohol carbon",
								0.145,	"CT" },
		{ this,	ELEMENT_H,	27,	"H1O",	"-&26",	"Hydrogen on alcoholic carbon",
								0.040,	0.12552,	2.50 },
		{ this,	ELEMENT_O,	28,	"OH",	"-&26,-H",	"Alcohol oxygen",
								-0.683,	0.71128,	3.21 },
		{ this,	ELEMENT_H,	29,	"HO",	"-&28",	"Alcohol hydrogen",
								0.418,	0.00000,	0.00 },
		// alkoxy  imidazolium
		{ this,	ELEMENT_C,	30,	"C2O",	"-C(nh=2,-N),-O(-H)",	"Carbon next to OH in hydroxyethyl group",
								0.275,	"CT" },
		// esters  OPLS-AA:JCC22:1340(2001)
		{ this,	ELEMENT_C,	31,	"CC",	"-C,=O,-O(nbonds=2)",	"Ester carbon",
								0.51,	0.43932,	3.75 },
		{ this,	ELEMENT_O,	32,	"OC",	"=&31",	"Carbonyl on ester carbon",
								-0.43,	0.87864,	2.96 },
		{ this,	ELEMENT_O,	33,	"OS",	"-&31",	"Alkoxy oxygen on ester carbon",
								-0.33,	0.71176,	3.00 },
		{ this,	ELEMENT_H,	34,	"HCO",	"-&35",	"Hydrogen on first alkoxy carbon",
								0.03,	0.06276,	2.42 },
		{ this,	ELEMENT_C,	35,	"CMO",	"nbonds=4,-&33",	"First alkoxy carbon",
								0.16,	"CT" },
		{ this,	ELEMENT_C,	36,	"CEO",	"-&35,nh=3",	"Carbon at end of ethoxy chain",
								0.19,	"CT" },
		// ether  OPLS-AA
		{ this,	ELEMENT_O,	37,	"OE",	"",	"",
								0.00,	0.87864,	2.96 },
		// fluoroalkylimidazolium  JPCB114:3592(2010)
		{ this,	ELEMENT_C,	38,	"C1H",	"",	"",
								-0.05,	"CT" },
		{ this,	ELEMENT_C,	39,	"C3F",	"",	"",
								0.12,	"CF" },
		// pyridinium  JCSPerkin2:2365(1999)
		{ this,	ELEMENT_N,	40,	"NAP",	"nbonds=3,ring(size=6,-C(n=5))",	"Nitrogen in alkylpyridinium",
								0.15,	"NA" },
		{ this,	ELEMENT_C,	41,	"CA1",	"nbonds=3,ring(size=6),-N,-C,-H",	"Carbon in ring,	adjacent to ring nitrogen (ortho)",
								0.00,	"CA" },
		{ this,	ELEMENT_C,	42,	"CA2",	"nbonds=3,ring(size=6),-C,-H,-C",	"Carbon in ring,	adjacent to two carbons (meta)",
								-0.07,	"CA" },
		{ this,	ELEMENT_C,	43,	"CA3",	"nbonds=3,ring(size=6),-C(n=2,-C(-N))",	"Carbon in ring,	adjacent to two carbons (para)",
								0.02,	"CA" },
		{ this,	ELEMENT_H,	44,	"HAP",	"nbonds=1,-C(ring(size=6,-C(n=5),-N),-C)",	"Hydrogen on ring carbon",
								0.15,	0.12552,	2.42 },
		// ammonium(also_pyrrolidinium)  OPLS-AA:JACS121:4827(1999),AMBER
		{ this,	ELEMENT_N,	45,	"N3",	"",	"",
								0.12,	0.71128,	3.25 },
		{ this,	ELEMENT_H,	46,	"H",	"",	"",
								0.13,	0.12552,	2.50 },
		// phosphonium  OPLS-AA
		{ this,	ELEMENT_P,	47,	"P3",	"nbonds=4,-C(n=4)",	"Phosphorus with four alkyl groups",
								0.68,	0.8368,	3.74 },
		{ this,	ELEMENT_C,	48,	"C1P",	"nbonds=4,-P(-C(n=3)),-H,-H,-C",	"First tetrahedral carbon in R group on phosphorus",
								-0.31,	"CT" },
		// tetrafluoroborate_anion  ""
		{ this,	ELEMENT_B,	49,	"B",	"nbonds=4,-F(n=4)",	"Boron with four fluorines attached",
								0.96,	0.3975,	3.58 },
		{ this,	ELEMENT_F,	50,	"FB",	"nbonds=1,-B(-F(n=4))",	"Fluorine attached to boron in BF4",
								-0.49,	0.2552,	3.12 },
		// hexafluorophosphate  JCSPerkin2:2365(1999)
		{ this,	ELEMENT_P,	51,	"P",	"nbonds=6,-F(n=6)",	"Phosphorus with six fluorines attahed",
								1.34,	0.8368,	3.74 },
		{ this,	ELEMENT_F,	52,	"F",	"nbonds=1,-P(-F(n=6))",	"Fluorine attached to phosphorus in PF6",
								-0.39,	0.2552,	3.12 },
		// nitrate  JPCB108:2038(2004)
		{ this,	ELEMENT_O,	53,	"ON",	"nbonds=1,-N(-O,-O)",	"Oxygen attached to N in NO3",
								0.95,	0.3380,	3.06 },
		{ this,	ELEMENT_N,	54,	"NN",	"nbonds=3,-O(n=3)",	"Nitrogen with three oxygens attached",
								-0.65,	0.6100,	2.77 },
		// chloride  JPCB108:2038(2004)
		{ this,	ELEMENT_CL,	55,	"Cl",	"unbound",	"Chloride ion",
								-1.00,	0.8300,	3.65 },
		// bromide  JPCB110:19586(2006)
		{ this,	ELEMENT_BR,	56,	"Br",	"unbound",	"Bromide ion",
								-1.00,	0.86,	3.97 },
		// bistriflylamide  JPCB108:16893(2004)
		{ this,	ELEMENT_F,	57,	"F1",	"-&58",	"Fluorine in TfO or NTf2",
								-0.16,	"F" },
		{ this,	ELEMENT_C,	58,	"CBT",	"-F(n=3),-S(-O(n>=2))",	"Fluorinated carbon in TfO or NTf2",
								0.35,	"CF" },
		{ this,	ELEMENT_S,	59,	"SBT",	"-O(n>=2),-[C,F]",	"Sulfur in TfO,	NTf2,	or bis(fluorosulfonyl)imide",
								1.02,	1.04600,	3.55 },
		{ this,	ELEMENT_O,	60,	"OBT",	"-S(-O(n=2))",	"Oxygen in bistrifylamide",
								-0.53,	0.87864,	2.96 },
		{ this,	ELEMENT_N,	61,	"NBT",	"-&59(n=2),nbonds=2",	"Nitrogen in bistrifylamide",
								-0.66,	0.71128,	3.25 },
		// triflate  anion
		{ this,	ELEMENT_O,	62,	"OTF",	"-S(-O(n=3))",	"Oxygen in triflate",
								-0.63,	"OBT" },
		// bis(fluorosulfonyl)amide  ""
		{ this,	ELEMENT_F,	63,	"FSI",	"-&59",	"Fluorine on sulfur in bis(fluorosulfonyl)amide",
								-0.13,	"F" },
		// longer_perfluoroalkanerosulfonylamides  ""
		{ this,	ELEMENT_C,	64,	"C1F",	"-&59,-C",	"First fluorinated carbon in chain on sulfur",
								0.19,	"CF" },
		// perfluoroalkanes  OPLS-AA
		{ this,	ELEMENT_C,	65,	"CTF",	"-F(n=3),-C",	"Terminal carbon in fluorinated chain",
								0.36,	"CF" },
		{ this,	ELEMENT_C,	66,	"CSF",	"-F(n=2),-C(n=2)",	"Middle carbon in fluorinated chain",
								0.24,	"CF" },
		{ this,	ELEMENT_F,	67,	"F",	"-[&64,&65,&66]",	"Fluorine in fluorinated chain",
								-0.12,	0.22175,	2.95 },
		// dicyanamide  JPCB110:19586(2006)
		{ this,	ELEMENT_N,	68,	"N3A",	"nbonds=2,-&69(n=2)",	"Central nitrogen in dicyanamide",
								-0.76,	0.71128,	3.25 },
		{ this,	ELEMENT_C,	69,	"CZA",	"nbonds=2,-N(bond=triple),-N",	"Carbon in dicyanamide",
								0.64,	0.27614,	3.30 },
		{ this,	ELEMENT_N,	70,	"NZA",	"nbonds=1,-C(bond=triple)",	"Terminal nitrogen in dicyanamide",
								-0.76,	0.71128,	3.20 },
		// alkylsulfates  JPCB112:5039(2008)
		{ this,	ELEMENT_S,	71,	"SO4",	"-O(n=3,nbonds=1),-O(-C)",	"Sulfur in alkylsulfate",
								1.18,	"SBT" },
		{ this,	ELEMENT_O,	72,	"OS4",	"-&71,nbonds=1",	"Oxygen in alkylsulfate",
								-0.65,	0.83700,	3.15 },
		{ this,	ELEMENT_O,	73,	"OC4",	"-&71,nbonds=2,-C",	"Linking oxygen in alkylsulfate",
								-0.45,	0.58600,	2.90 },
		{ this,	ELEMENT_C,	74,	"CS4",	"-&73",	"First carbon in chain in alkylsulfate",
								0.22,	"CT" },
		{ this,	ELEMENT_H,	75,	"HS4",	"-&74",	"Hydrogen on first carbon in chain in alkylsulfate",
								0.00,	"HS4" },
		// alkylsulfonates  JPCB112:5039(2008)
		{ this,	ELEMENT_S,	76,	"SO3",	"-O(n=3,nbonds=1),-C(-[H,C])",	"Sulfur in alkylsulfonate",
								1.18,	"SBT" },
		{ this,	ELEMENT_O,	77,	"OS3",	"-&76,nbonds=1",	"Oxygen in alkylsulfonate",
								-0.68,	0.83700,	3.15 },
		{ this,	ELEMENT_C,	78,	"CS3",	"-&76",	"First carbon in chain in alkylsulfonate",
								-0.14,	"CT" },
		{ this,	ELEMENT_H,	79,	"HS3",	"-&78",	"Hydrogen on first carbon in chain in alkylsulfonate",
								0.00,	"HC" },
		// acetate  OPLS-AA
		{ this,	ELEMENT_O,	80,	"O2",	"-&81",	"Acetate oxygen",
								-0.80,	0.87860,	2.96 },
		{ this,	ELEMENT_C,	81,	"CO2",	"-C(nh=3),-O(n=2,nbonds=1)",	"Carbon with two oxygens in acetate",
								0.70,	0.43930,	3.75 },
		{ this,	ELEMENT_C,	82,	"CTA",	"nh=3,-&81",	"Methyl carbon in acetate",
								-0.28,	"CT" },
		// Extra/Fixes  ""
		{ this,	ELEMENT_C,	83,	"CT4",	"-&74,nh=3",	"Terminal carbon in R=ethyl group of alkylsulfonate",
								-0.18,	"CT" }
	};

	static ForcefieldBondTerm bondTerms[] =
	{
		//      i       j       Type (Harmonic)                 k               eq
		// alkanes  OPLS-AA
		{ this,	"HC*",	"CT",	SpeciesBond::HarmonicForm,	2845,	1.090 },
		{ this,	"CT",	"CT",	SpeciesBond::HarmonicForm,	2242,	1.529 },
		// dialkylimidazolium  JPCB108:2038(2004)
		{ this,	"CR",	"HA*",	SpeciesBond::HarmonicForm,	2845,	1.080 },
		{ this,	"CW",	"HA*",	SpeciesBond::HarmonicForm,	2845,	1.080 },
		{ this,	"CR",	"NA",	SpeciesBond::HarmonicForm,	3992,	1.315 },
		{ this,	"CW",	"NA",	SpeciesBond::HarmonicForm,	3574,	1.378 },
		{ this,	"CW",	"CW",	SpeciesBond::HarmonicForm,	4352,	1.341 },
		{ this,	"NA",	"CT",	SpeciesBond::HarmonicForm,	2820,	1.466 },
		// alkylimidazolium  OPLS-AA:JPCB110:19586(2006)
		{ this,	"NA",	"HA*",	SpeciesBond::HarmonicForm,	3632,	1.010 },
		// alkenylimidazolium  OPLS-AA
		{ this,	"HCM",	"CM",	SpeciesBond::HarmonicForm,	2845.12,	1.080 },
		{ this,	"CT",	"CM",	SpeciesBond::HarmonicForm,	2652.65,	1.510 },
		{ this,	"CM",	"CM",	SpeciesBond::HarmonicForm,	4594.03,	1.340 },
		// alcohols  OPLS-AA:JACS118:11225(1996);JPC100:18010(1996)
		{ this,	"HO",	"OH",	SpeciesBond::HarmonicForm,	4628,	0.945 },
		{ this,	"CT",	"OH",	SpeciesBond::HarmonicForm,	2678,	1.410 },
		// pyridinium  OPLS-AA:Theochem424:145(1998)
		{ this,	"CA",	"HA*",	SpeciesBond::HarmonicForm,	3071,	1.080 },
		{ this,	"CA",	"CA",	SpeciesBond::HarmonicForm,	3925,	1.380 },
		{ this,	"CA",	"NA",	SpeciesBond::HarmonicForm,	4042,	1.340 },
		{ this,	"NA",	"CT",	SpeciesBond::HarmonicForm,	2820,	1.480 },
		// ammonium(also_pyrrolidinium)  OPLS-AA:JACS121:4827(1999),AMBER
		{ this,	"H",	"N3",	SpeciesBond::HarmonicForm,	3632,	1.010 },
		{ this,	"N3",	"CT",	SpeciesBond::HarmonicForm,	3071,	1.471 },
		// dialkylimethylmidazolium  ""
		{ this,	"CR",	"CT",	SpeciesBond::HarmonicForm,	2653,	1.510 },
		// esters  OPLS-AA
		{ this,	"CC",	"OC",	SpeciesBond::HarmonicForm,	4770.0,	1.200 },
		{ this,	"CT",	"CC",	SpeciesBond::HarmonicForm,	2653.0,	1.522 },
		{ this,	"CC",	"OS",	SpeciesBond::HarmonicForm,	1791.0,	1.344 },
		{ this,	"CT",	"OS",	SpeciesBond::HarmonicForm,	2678.0,	1.437 },
		// ethers  OPLS-AA
		{ this,	"CT",	"OE",	SpeciesBond::HarmonicForm,	2678.0,	1.437 },
		// phosphonium  OPLS-AA:JPCB110:19586(2006)
		{ this,	"P3",	"CT",	SpeciesBond::HarmonicForm,	3550,	1.81 },
		// tetrafluoroborate  ""
		{ this,	"B",	"FB",	SpeciesBond::HarmonicForm,	3235,	1.394 },
		// hexafluorophosphate  JCSPerkin2:2365(1999)
		{ this,	"P",	"F",	SpeciesBond::HarmonicForm,	3100,	1.606 },
		// nitrate  JPCB108:2038(2004)
		{ this,	"NN",	"ON",	SpeciesBond::HarmonicForm,	5307,	1.226 },
		// triflate_and_bistriflylamide  JPCB108:16893(2004)
//                { this,	"CBT",	"FBT",	SpeciesBond::HarmonicForm,	harm,	1.323 },
		{ this,	"CF",	"SBT",	SpeciesBond::HarmonicForm,	1950,	1.818 },
		{ this,	"SBT",	"OBT",	SpeciesBond::HarmonicForm,	5331,	1.437 },
		{ this,	"NBT",	"SBT",	SpeciesBond::HarmonicForm,	3137,	1.570 },
		// bis(fluorosulfonyl)amide  ""
		{ this,	"F",	"SBT",	SpeciesBond::HarmonicForm,	1879,	1.575 },
		// longer_perfluoroalkanerosulfonylamides  ""
		// perfluoroalkanes  OPLS-AA:JPCA105:4118(2001)
		{ this,	"F",	"CF",	SpeciesBond::HarmonicForm,	3071.06,	1.332 },
		{ this,	"CF",	"CF",	SpeciesBond::HarmonicForm,	2242.62,	1.529 },
		{ this,	"CT",	"CF",	SpeciesBond::HarmonicForm,	2242.62,	1.529 },
		// dicyanamide  JPCB110:19586(2006)
		{ this,	"N3A",	"CZA",	SpeciesBond::HarmonicForm,	4206,	1.310 },
		{ this,	"CZA",	"NZA",	SpeciesBond::HarmonicForm,	7746,	1.157 },
		// alkylsulfates  JPCB112:5039(2008)
		{ this,	"CT",	"OC4",	SpeciesBond::HarmonicForm,	745.8,	1.402 },
		{ this,	"OS4",	"SO",	SpeciesBond::HarmonicForm,	5331,	1.455 },
		{ this,	"OC4",	"SO",	SpeciesBond::HarmonicForm,	1789.6,	1.633 },
		// alkylsulfonates  JPCB112:5039(2008)
		{ this,	"CT",	"SO",	SpeciesBond::HarmonicForm,	1970,	1.792 },
		{ this,	"OS3",	"SO",	SpeciesBond::HarmonicForm,	5331,	1.455 },
		// acetate  OPLS-AA:JPCB(2004)
		{ this,	"CO2",	"O2",	SpeciesBond::HarmonicForm,	5489.0,	1.250 },
		{ this,	"CT",	"CO2",	SpeciesBond::HarmonicForm,	2653.0,	1.522 }
	};

	static ForcefieldAngleTerm angleTerms[] =
	{
		//      i       j       k       Type (Harmonic)                 k            eq
		// alkanes  OPLS-AA:JACS118:11225(1996);JPC100:18010(1996)
		{ this,	"CT",	"CT",	"CT",	SpeciesAngle::HarmonicForm,	488.3,	112.7 },
		{ this,	"CT",	"CT",	"HC*",	SpeciesAngle::HarmonicForm,	313.8,	110.7 },
		{ this,	"HC*",	"CT",	"HC*",	SpeciesAngle::HarmonicForm,	276.1,	107.8 },
		// dialkylimidazolium  JPCB108:2038(2004)
		{ this,	"CW",	"NA",	"CR",	SpeciesAngle::HarmonicForm,	585.8,	108.0 },
		{ this,	"CW",	"NA",	"CT",	SpeciesAngle::HarmonicForm,	585.8,	125.6 },
		{ this,	"CR",	"NA",	"CT",	SpeciesAngle::HarmonicForm,	585.8,	126.4 },
		{ this,	"NA",	"CR",	"HA*",	SpeciesAngle::HarmonicForm,	292.9,	125.1 },
		{ this,	"NA",	"CR",	"NA",	SpeciesAngle::HarmonicForm,	585.8,	109.8 },
		{ this,	"NA",	"CW",	"CW",	SpeciesAngle::HarmonicForm,	585.8,	107.1 },
		{ this,	"NA",	"CW",	"HA*",	SpeciesAngle::HarmonicForm,	292.9,	122.0 },
		{ this,	"CW",	"CW",	"HA*",	SpeciesAngle::HarmonicForm,	292.9,	130.9 },
		{ this,	"NA",	"CT",	"HC*",	SpeciesAngle::HarmonicForm,	313.8,	110.7 },
		{ this,	"NA",	"CT",	"CT",	SpeciesAngle::HarmonicForm,	488.3,	112.7 },
		// alkylimlidazolium  JPCB110:19586(2006)
		{ this,	"CR",	"NA",	"HA*",	SpeciesAngle::HarmonicForm,	292.9,	125.4 },
		{ this,	"CW",	"NA",	"HA*",	SpeciesAngle::HarmonicForm,	292.9,	126.6 },
		// dialkylimethylmidazolium  JPCB112:5039(2008)
		{ this,	"CT",	"CR",	"NA",	SpeciesAngle::HarmonicForm,	585.8,	125.8 },
		{ this,	"CR",	"CT",	"HC*",	SpeciesAngle::HarmonicForm,	313.8,	110.7 },
		// alkenylimidazolium  OPLS-AA
		{ this,	"CT",	"CT",	"CM",	SpeciesAngle::HarmonicForm,	527.18,	111.1 },
		{ this,	"CM",	"CM",	"CT",	SpeciesAngle::HarmonicForm,	585.76,	124.0 },
		{ this,	"HCM",	"CM",	"CM",	SpeciesAngle::HarmonicForm,	292.88,	120.0 },
		{ this,	"HCM",	"CM",	"CT",	SpeciesAngle::HarmonicForm,	292.88,	117.0 },
		{ this,	"HCM",	"CT",	"CM",	SpeciesAngle::HarmonicForm,	292.88,	109.5 },
		{ this,	"CM",	"CT",	"HC*",	SpeciesAngle::HarmonicForm,	292.88,	109.5 },
		{ this,	"HCM",	"CM",	"HCM",	SpeciesAngle::HarmonicForm,	292.88,	114.3 },
		// alcohols  OPLS-AA
		{ this,	"CT",	"CT",	"OH",	SpeciesAngle::HarmonicForm,	418.40,	109.50 },
		{ this,	"HC*",	"CT",	"OH",	SpeciesAngle::HarmonicForm,	292.88,	109.50 },
		{ this,	"CT",	"OH",	"HO",	SpeciesAngle::HarmonicForm,	460.24,	108.50 },
		// esters  OPLS-AA
		{ this,	"CT",	"CC",	"OC",	SpeciesAngle::HarmonicForm,	669.00,	125.00 },
		{ this,	"HC*",	"CT",	"CC",	SpeciesAngle::HarmonicForm,	313.80,	109.50 },
		{ this,	"OS",	"CC",	"OC",	SpeciesAngle::HarmonicForm,	694.50,	125.00 },
		{ this,	"OS",	"CC",	"CT",	SpeciesAngle::HarmonicForm,	677.80,	110.00 },
		{ this,	"CT",	"OS",	"CC",	SpeciesAngle::HarmonicForm,	694.50,	115.00 },
		{ this,	"HC*",	"CT",	"OS",	SpeciesAngle::HarmonicForm,	418.40,	109.50 },
		{ this,	"CT",	"CT",	"OS",	SpeciesAngle::HarmonicForm,	418.40,	109.50 },
		{ this,	"CT",	"CT",	"CC",	SpeciesAngle::HarmonicForm,	527.18,	111.10 },
		{ this,	"NA",	"CT",	"CC",	SpeciesAngle::HarmonicForm,	500.00,	112.7 },
		// ethers  OPLS-AA
//                { this,	"CT",	"CT",	"OE",	SpeciesAngle::HarmonicForm,	109.50,	418.40 },
//                { this,	"CT",	"OE",	"CT",	SpeciesAngle::HarmonicForm,	115.00,	694.50 },
//                { this,	"HC*",	"CT",	"OE",	SpeciesAngle::HarmonicForm,	109.50,	418.40 },
		// semifluorinated_alkanes  JPCA106:10116(2002)
		{ this,	"HC*",	"CT",	"CF",	SpeciesAngle::HarmonicForm,	313.80,	110.70 },
		{ this,	"F",	"CF",	"CT",	SpeciesAngle::HarmonicForm,	418.40,	109.50 },
		{ this,	"CT",	"CF",	"CF",	SpeciesAngle::HarmonicForm,	488.27,	112.70 },
		{ this,	"CT",	"CT",	"CF",	SpeciesAngle::HarmonicForm,	488.27,	112.70 },
		// pyridinium  OPLS-AA:Theochem424:145(1998),JCSPerkin2:2365(1999)
		{ this,	"CA",	"CA",	"CA",	SpeciesAngle::HarmonicForm,	527.2,	120.0 },
		{ this,	"CA",	"CA",	"HA*",	SpeciesAngle::HarmonicForm,	292.9,	120.0 },
		{ this,	"CA",	"CA",	"NA",	SpeciesAngle::HarmonicForm,	585.8,	120.0 },
		{ this,	"CA",	"NA",	"CA",	SpeciesAngle::HarmonicForm,	585.8,	120.4 },
		{ this,	"CA",	"NA",	"CT",	SpeciesAngle::HarmonicForm,	585.8,	119.8 },
		{ this,	"NA",	"CA",	"HA*",	SpeciesAngle::HarmonicForm,	292.9,	120.0 },
		{ this,	"NA",	"CT",	"HC*",	SpeciesAngle::HarmonicForm,	292.9,	109.5 },
		// ammonium(also_pyrrolidinium)  OPLS-AA:JACS121:4827(1999),AMBER
		{ this,	"N3",	"CT",	"CT",	SpeciesAngle::HarmonicForm,	669.4,	109.5 },
		{ this,	"CT",	"N3",	"CT",	SpeciesAngle::HarmonicForm,	418.4,	109.5 },
		{ this,	"HC*",	"CT",	"N3",	SpeciesAngle::HarmonicForm,	209.2,	109.5 },
		{ this,	"H",	"N3",	"CT",	SpeciesAngle::HarmonicForm,	418.4,	109.5 },
		// phosphosphonium  OPLS-AA:JPCB110:19586(2006)
		{ this,	"CT",	"P3",	"CT",	SpeciesAngle::HarmonicForm,	607.8,	109.5 },
		{ this,	"HC*",	"CT",	"P3",	SpeciesAngle::HarmonicForm,	389.9,	110.1 },
		{ this,	"CT",	"CT",	"P3",	SpeciesAngle::HarmonicForm,	509.1,	115.2 },
		// tetrafluoroborate_anion  ""
		{ this,	"FB",	"B",	"FB",	SpeciesAngle::HarmonicForm,	669.5,	109.5 },
		// hexafluorophosphate  JCSPerkin2:2365(1999)
		{ this,	"F",	"P",	"F",	SpeciesAngle::HarmonicForm,	1165,	90.0 },
		// nitrate  JPCB108:2038(2004)
		{ this,	"ON",	"NN",	"ON",	SpeciesAngle::HarmonicForm,	1011,	120.0 },
		// triflate_and_bistriflylamide  JPCB108:16893(2004)
//                { this,	"FBT",	"CBT",	"FBT",	SpeciesAngle::HarmonicForm,	107.1,	781 },
		{ this,	"F",	"CF",	"SBT",	SpeciesAngle::HarmonicForm,	694,	111.7 },
		{ this,	"OBT",	"SBT",	"OBT",	SpeciesAngle::HarmonicForm,	969,	118.5 },
		{ this,	"CF",	"SBT",	"OBT",	SpeciesAngle::HarmonicForm,	870,	102.6 },
		{ this,	"NBT",	"SBT",	"OBT",	SpeciesAngle::HarmonicForm,	789,	113.6 },
		{ this,	"NBT",	"SBT",	"CF",	SpeciesAngle::HarmonicForm,	764,	103.5 },
		{ this,	"SBT",	"NBT",	"SBT",	SpeciesAngle::HarmonicForm,	671,	125.6 },
		// bis(fluorosulfonyl)amide  ""
		{ this,	"F",	"SBT",	"OBT",	SpeciesAngle::HarmonicForm,	1077,	104.1 },
		{ this,	"F",	"SBT",	"NBT",	SpeciesAngle::HarmonicForm,	902,	103.0 },
		// longer_perfluoroalkanerosulfonylamides  ""
		{ this,	"SBT",	"CF",	"CF",	SpeciesAngle::HarmonicForm,	418.4,	115.9 },
		// perfluoroalkanes  OPLS-AA:JPCA105:4118(2001)
		{ this,	"F",	"CF",	"F",	SpeciesAngle::HarmonicForm,	644.34,	109.10 },
		{ this,	"F",	"CF",	"CF",	SpeciesAngle::HarmonicForm,	418.40,	109.50 },
		{ this,	"CF",	"CF",	"CF",	SpeciesAngle::HarmonicForm,	488.27,	112.70 },
		// dicyanamide  JPCB110:19586(2006)
		{ this,	"CZA",	"N3A",	"CZA",	SpeciesAngle::HarmonicForm,	362,	118.5 },
		{ this,	"N3A",	"CZA",	"NZA",	SpeciesAngle::HarmonicForm,	425,	175.2 },
		// alkylsulfates  JPCB112:5039(2008)
		{ this,	"OS4",	"SO",	"OS4",	SpeciesAngle::HarmonicForm,	969.00,	114.00 },
		{ this,	"OC4",	"SO",	"OS4",	SpeciesAngle::HarmonicForm,	1239.62,	103.50 },
		{ this,	"CT",	"OC4",	"SO",	SpeciesAngle::HarmonicForm,	300.45,	116.60 },
		{ this,	"HC*",	"CT",	"OC4",	SpeciesAngle::HarmonicForm,	488.68,	109.70 },
		{ this,	"CT",	"CT",	"OC4",	SpeciesAngle::HarmonicForm,	765.57,	107.80 },
		// alkylsulfonates  JPCB112:5039(2008)
		{ this,	"OS3",	"SO",	"OS3",	SpeciesAngle::HarmonicForm,	969.00,	114.00 },
		{ this,	"CT",	"SO",	"OS3",	SpeciesAngle::HarmonicForm,	870.00,	104.50 },
		{ this,	"HC*",	"CT",	"SO",	SpeciesAngle::HarmonicForm,	390.30,	107.30 },
		{ this,	"CT",	"CT",	"SO",	SpeciesAngle::HarmonicForm,	583.00,	113.30 },
		// acetate  OPLS-AA:JPCB(2004)
		{ this,	"O2",	"CO2",	"O2",	SpeciesAngle::HarmonicForm,	669.4,	126.00 },
		{ this,	"CT",	"CO2",	"O2",	SpeciesAngle::HarmonicForm,	585.8,	117.00 },
		{ this,	"HC*",	"CT",	"CO2",	SpeciesAngle::HarmonicForm,	298.9,	109.50 }
	};

	static ForcefieldTorsionTerm torsionTerms[] =
	{
		//      i       j       k      l        Type (CosineForm)              v1           v2          v3        v4
		// alkanes  OPLS-AA
		{ this,	"HC*",	"CT",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.3305,	0.0000 },
		{ this,	"CT",	"CT",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.5313,	0.0000 },
		{ this,	"CT",	"CT",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	7.2800,	-0.6569,	1.1673,	0.0000 },
		// dialkylimidazolium  JPCB108:2038(2004)
		{ this,	"CW",	"NA",	"CR",	"NA",	SpeciesTorsion::Cos3Form,	0.0000,	19.4600,	0.0000,	0.0000 },
		{ this,	"CW",	"NA",	"CR",	"HA*",	SpeciesTorsion::Cos3Form,	0.0000,	19.4600,	0.0000,	0.0000 },
		{ this,	"CT",	"NA",	"CR",	"NA",	SpeciesTorsion::Cos3Form,	0.0000,	19.4600,	0.0000,	0.0000 },
		{ this,	"CT",	"NA",	"CR",	"HA*",	SpeciesTorsion::Cos3Form,	0.0000,	19.4600,	0.0000,	0.0000 },
		{ this,	"CR",	"NA",	"CW",	"CW",	SpeciesTorsion::Cos3Form,	0.0000,	12.5500,	0.0000,	0.0000 },
		{ this,	"CR",	"NA",	"CW",	"HA*",	SpeciesTorsion::Cos3Form,	0.0000,	12.5500,	0.0000,	0.0000 },
		{ this,	"CT",	"NA",	"CW",	"CW",	SpeciesTorsion::Cos3Form,	0.0000,	12.5500,	0.0000,	0.0000 },
		{ this,	"CT",	"NA",	"CW",	"HA*",	SpeciesTorsion::Cos3Form,	0.0000,	12.5500,	0.0000,	0.0000 },
		{ this,	"NA",	"CW",	"CW",	"NA",	SpeciesTorsion::Cos3Form,	0.0000,	44.9800,	0.0000,	0.0000 },
		{ this,	"NA",	"CW",	"CW",	"HA*",	SpeciesTorsion::Cos3Form,	0.0000,	44.9800,	0.0000,	0.0000 },
		{ this,	"HA*",	"CW",	"CW",	"HA*",	SpeciesTorsion::Cos3Form,	0.0000,	44.9800,	0.0000,	0.0000 },
		{ this,	"CW",	"NA",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.5190,	0.0000 },
		{ this,	"CR",	"NA",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
		{ this,	"CW",	"NA",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	-7.1535,	6.1064,	0.7939,	0.0000 },
		{ this,	"CR",	"NA",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	-5.2691,	0.0000,	0.0000,	0.0000 },
		{ this,	"NA",	"CT",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	-7.4797,	3.1642,	-1.2026,	0.0000 },
		{ this,	"NA",	"CT",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.3670,	0.0000 },
		// alkylimidazolium_ring  AMBER
		{ this,	"HA*",	"NA",	"CR",	"NA",	SpeciesTorsion::Cos3Form,	0.0000,	19.4600,	0.0000,	0.0000 },
		{ this,	"HA*",	"NA",	"CR",	"HA*",	SpeciesTorsion::Cos3Form,	0.0000,	19.4600,	0.0000,	0.0000 },
		{ this,	"HA*",	"NA",	"CW",	"CW",	SpeciesTorsion::Cos3Form,	0.0000,	12.5500,	0.0000,	0.0000 },
		{ this,	"HA*",	"NA",	"CW",	"HA*",	SpeciesTorsion::Cos3Form,	0.0000,	12.5500,	0.0000,	0.0000 },
		// dialkylimethylmidazolium  JPCB112:5039(2008)
		{ this,	"CW",	"NA",	"CR",	"CT",	SpeciesTorsion::Cos3Form,	0.0000,	19.4600,	0.0000,	0.0000 },
		{ this,	"CT",	"NA",	"CR",	"CT",	SpeciesTorsion::Cos3Form,	0.0000,	19.4600,	0.0000,	0.0000 },
		{ this,	"NA",	"CR",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
		// alkenylimidazolium  (OPLS-AA)
		{ this,	"CT",	"CT",	"CM",	"HCM",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.3305,	0.0000 },
		{ this,	"CM",	"CT",	"CT",	"H",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.5313,	0.0000 },
		{ this,	"HCM",	"CM",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.3305,	0.0000 },
		{ this,	"CT",	"CM",	"CM",	"HCM",	SpeciesTorsion::Cos3Form,	0.0000,	58.5760,	0.0000,	0.0000 },
		{ this,	"HCM",	"CM",	"CM",	"HCM",	SpeciesTorsion::Cos3Form,	0.0000,	58.5760,	0.0000,	0.0000 },
		{ this,	"HCM",	"CM",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.3305,	0.0000 },
		{ this,	"CM",	"CM",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	-1.5565,	0.0000 },
		{ this,	"CM",	"CT",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	5.4392,	-0.2092,	0.8368,	0.0000 },
		{ this,	"HCM",	"CM",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.33051,	0.0000 },
		{ this,	"CM",	"CM",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	-1.5565,	0.0000 },
		{ this,	"CM",	"CT",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.5313,	0.0000 },
		{ this,	"CM",	"CM",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	1.4477,	1.6945,	-3.7823,	0.0000 },
		{ this,	"CM",	"CT",	"CT",	"NA",	SpeciesTorsion::Cos3Form,	-6.5464,	-0.9714,	4.2205,	0.0000 },
		// alcohols  OPLS-AA:JACS118:11225(1996);JPC100:18010(1996),AMBER98(OCCO)
		{ this,	"HC*",	"CT",	"OH",	"HO",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.8828,	0.0000 },
		{ this,	"CT",	"CT",	"OH",	"HO",	SpeciesTorsion::Cos3Form,	-1.4895,	-0.7280,	2.0585,	0.0000 },
		{ this,	"HC*",	"CT",	"CT",	"OH",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.9581,	0.0000 },
		{ this,	"OH",	"CT",	"CT",	"OH",	SpeciesTorsion::Cos3Form,	0.0000,	-9.8324,	1.2050,	0.0000 },
		// alkoxy_imidazolium  ""
		{ this,	"NA",	"CT",	"CT",	"OH",	SpeciesTorsion::Cos3Form,	-3.5787,	-1.6564,	4.9154,	0.0000 },
		// esters  OPLS-AA:JCC22:1340(2001)
		{ this,	"CT",	"CT",	"CC",	"OC",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
		{ this,	"HC*",	"CT",	"CC",	"OC",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
		{ this,	"CT",	"CT",	"CC",	"OS",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	-1.4853,	0.0000 },
		{ this,	"HC*",	"CT",	"CC",	"OS",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.5523,	0.0000 },
		{ this,	"CT",	"CC",	"OS",	"CT",	SpeciesTorsion::Cos3Form,	19.5351,	21.4388,	0.0000,	0.0000 },
		{ this,	"OC",	"CC",	"OS",	"CT",	SpeciesTorsion::Cos3Form,	0.0000,	21.4388,	0.0000,	0.0000 },
		{ this,	"CC",	"OS",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.8284,	0.0000 },
		{ this,	"CC",	"OS",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	-5.1045,	-0.5272,	1.7656,	0.0000 },
		{ this,	"HC*",	"CT",	"CT",	"CC",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	-1.4477,	0.0000 },
		{ this,	"CT",	"CT",	"CT",	"CC",	SpeciesTorsion::Cos3Form,	-7.1002,	-1.9079,	2.4476,	0.0000 },
		{ this,	"HC*",	"CT",	"CT",	"OS",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.6067,	0.0000 },
		{ this,	"CT",	"CT",	"CT",	"OS",	SpeciesTorsion::Cos3Form,	11.4307,	-0.9581,	2.0292,	0.0000 },
		{ this,	"CW",	"NA",	"CT",	"CC",	SpeciesTorsion::Cos3Form,	4.1274,	0.0000,	0.0000,	0.0000 },
		{ this,	"CR",	"NA",	"CT",	"CC",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
		{ this,	"NA",	"CT",	"CC",	"OC",	SpeciesTorsion::Cos3Form,	-9.1642,	14.1359,	1.0771,	0.0000 },
		{ this,	"NA",	"CT",	"CC",	"OS",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
		// ethers  OPLS-AA
//                { this,	"HC*",	"CT",	"CT",	"OE",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
//                { this,	"CT",	"CT",	"CT",	"OE",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
//                { this,	"OS",	"CT",	"CT",	"OE",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
//                { this,	"OE",	"CT",	"CT",	"OE",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
//                { this,	"HC*",	"CT",	"OE",	"CT",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
//                { this,	"CT",	"CT",	"OE",	"CT",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
		// semifluorinated_alkanes  JPCA:106:10116(2002)
		{ this,	"HC*",	"CT",	"CF",	"F",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.2134,	0.0000 },
		{ this,	"F",	"CF",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.2134,	0.0000 },
		{ this,	"F",	"CF",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.9372,	0.0000 },
		{ this,	"F",	"CF",	"CF",	"CT",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	2.7656,	0.0000 },
		{ this,	"HC*",	"CT",	"CT",	"CF",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.5565,	0.0000 },
		{ this,	"HC*",	"CT",	"CF",	"CF",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.7573,	0.0000 },
		{ this,	"CF",	"CT",	"CT",	"NA",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
		// pyridinium  AMBER(cycle)JPCB110:19586(2006)
		{ this,	"CA",	"CA",	"CA",	"CA",	SpeciesTorsion::Cos3Form,	0.0000,	30.3340,	0.0000,	0.0000 },
		{ this,	"NA",	"CA",	"CA",	"CA",	SpeciesTorsion::Cos3Form,	0.0000,	30.3340,	0.0000,	0.0000 },
		{ this,	"HA*",	"CA",	"CA",	"CA",	SpeciesTorsion::Cos3Form,	0.0000,	30.3340,	0.0000,	0.0000 },
		{ this,	"HA*",	"CA",	"CA",	"NA",	SpeciesTorsion::Cos3Form,	0.0000,	30.3340,	0.0000,	0.0000 },
		{ this,	"HA*",	"CA",	"CA",	"HA*",	SpeciesTorsion::Cos3Form,	0.0000,	30.3340,	0.0000,	0.0000 },
		{ this,	"CA",	"NA",	"CA",	"HA*",	SpeciesTorsion::Cos3Form,	0.0000,	12.5520,	0.0000,	0.0000 },
		{ this,	"CT",	"NA",	"CA",	"HA*",	SpeciesTorsion::Cos3Form,	0.0000,	12.5520,	0.0000,	0.0000 },
		{ this,	"CA",	"NA",	"CA",	"CA",	SpeciesTorsion::Cos3Form,	0.0000,	12.5520,	0.0000,	0.0000 },
		{ this,	"CT",	"NA",	"CA",	"CA",	SpeciesTorsion::Cos3Form,	0.0000,	12.5520,	0.0000,	0.0000 },
		{ this,	"HC*",	"CT",	"NA",	"CA",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
		// ammonium(also_pyrrolidinium)  OPLS-AA:JACS121:4827(1999)
		{ this,	"HC*",	"CT",	"CT",	"N3",	SpeciesTorsion::Cos3Form,	-4.2384,	-2.9665,	1.9790,	0.0000 },
		{ this,	"CT",	"CT",	"CT",	"N3",	SpeciesTorsion::Cos3Form,	10.0081,	-2.8200,	2.3012,	0.0000 },
		{ this,	"CT",	"N3",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	1.7405,	-0.5356,	2.9079,	0.0000 },
		{ this,	"HC*",	"CT",	"N3",	"CT",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	2.3430,	0.0000 },
		{ this,	"H",	"N3",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
		{ this,	"H",	"N3",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
		// Cholinium:  unpublished(Ferid+Agilio)
		{ this,	"OH",	"CT",	"CT",	"N3",	SpeciesTorsion::Cos3Form,	-44.0515,	-5.4349,	0.0000,	0.0000 },
		// phosphonium  OPLS-AA
		{ this,	"CT",	"P3",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.9270,	0.0000 },
		{ this,	"CT",	"P3",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.1330,	0.0000 },
		{ this,	"P3",	"CT",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.4650,	0.0000 },
		{ this,	"P3",	"CT",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	-3.2480,	0.9880,	-0.7150,	0.0000 },
		// triflate_and_bistriflylamide  JPCB108:16893(2004)
		{ this,	"OBT",	"SBT",	"CF",	"F",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.4510,	0.0000 },
		{ this,	"NBT",	"SBT",	"CF",	"F",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.3220,	0.0000 },
		{ this,	"OBT",	"SBT",	"NBT",	"SBT",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	-0.0150,	0.0000 },
		{ this,	"SBT",	"NBT",	"SBT",	"CF",	SpeciesTorsion::Cos3Form,	2.7730,	-10.4200,	-3.1950,	0.0000 },
		// longer_perfluoroalkanerosulfonylamides  ""
		{ this,	"SBT",	"CF",	"CF",	"F",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.4530,	0.0000 },
		{ this,	"OBT",	"SBT",	"CF",	"CF",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	-0.7400,	0.0000 },
		{ this,	"NBT",	"SBT",	"CF",	"CF",	SpeciesTorsion::Cos3Form,	-3.0940,	0.0000,	0.0000,	0.0000 },
		// perfluoroalkanes  OPLS-AAJPCA105:4118(2001)
		{ this,	"F",	"CF",	"CF",	"F",	SpeciesTorsion::Cos3Form,	-10.4600,	0.0000,	1.0460,	0.0000 },
		{ this,	"F",	"CF",	"CF",	"CF",	SpeciesTorsion::Cos3Form,	1.2552,	0.0000,	1.6736,	0.0000 },
		// dicyanamide  JPCB110:19586(2006)
		{ this,	"NZA",	"CZA",	"N3A",	"CZA",	SpeciesTorsion::Cos3Form,	4.0800,	0.0000,	0.0000,	0.0000 },
		// alkylsulfates  JPCB112:5039(2008)
		{ this,	"OS4",	"SO",	"OC4",	"CT",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	2.4815,	0.0000 },
		{ this,	"SO",	"OC4",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.6858,	0.0000 },
		{ this,	"SO",	"OC4",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	-6.0142,	-3.1133,	1.4941,	0.0000 },
		{ this,	"OC4",	"CT",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	2.0698,	0.0000 },
		{ this,	"OC4",	"CT",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	4.3893,	-1.8273,	2.9705,	0.0000 },
		// alkylsulfonates  JPCB112:5039(2008)
		{ this,	"OS3",	"SO",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.6250,	0.0000 },
		{ this,	"OS3",	"SO",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.3938,	0.0000 },
		{ this,	"SO",	"CT",	"CT",	"HC*",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	1.3797,	0.0000 },
		{ this,	"SO",	"CT",	"CT",	"CT",	SpeciesTorsion::Cos3Form,	-16.1000,	-2.0046,	0.7674,	0.0000 },
		// acetate  OPLS-AA:JPCB(2004)
		{ this,	"HC*",	"CT",	"CO2",	"O2",	SpeciesTorsion::Cos3Form,	0.0000,	0.0000,	0.0000,	0.0000 },
		// semifluorinated_alkanes  JPCA:106:10116(2002)
		{ this,	"CF",	"CF",	"CF",	"CT",	SpeciesTorsion::Cos4Form,	14.6750,	-0.9179,	-2.8980,	-2.0174 },
		{ this,	"CF",	"CF",	"CT",	"CT",	SpeciesTorsion::Cos4Form,	-0.5881,	-0.7642,	-0.3170,	-0.3179 },
		// pyridinium  AMBER(cycle)JPCB110:19586(2006)
		{ this,	"CT",	"CT",	"NA",	"CA",	SpeciesTorsion::Cos4Form,	0.0000,	1.1120,	0.0000,	0.6900 },
		// Cholinium:  unpublished(Ferid+Agilio)
//                { this,	"OH",	"CT",	"CT",	"N3",	SpeciesTorsion::Cos4Form,	-44.0515,	-5.0148,	0.0000,	-3.1510 },
		// bis(fluorosulfonyl)amide  ""
		{ this,	"F",	"SBT",	"NBT",	"SBT",	SpeciesTorsion::Cos4Form,	13.4514,	-12.3373,	-8.4874,	-2.8654 },
		// longer_perfluoroalkanerosulfonylamides  ""
		{ this,	"SBT",	"CF",	"CF",	"CF",	SpeciesTorsion::Cos4Form,	50.0900,	0.0000,	-4.6260,	-4.0080 },
		// perfluoroalkanes  OPLS-AA:JPCA105:4118(2001)
		{ this,	"CF",	"CF",	"CF",	"CF",	SpeciesTorsion::Cos4Form,	27.7064,	3.9664,	-5.8074,	-8.8617 }
	};

	static ForcefieldImproperTerm improperTerms[] =
	{
		//      i       j       k      l        Type (CosineForm)              v1           v2          v3     v4
		// improper  pyridinium
//                { this,	"CA",	"CA",	"CA",	"HA*",	SpeciesImproper::Cos3,	0.0000,	9.2000,	0.0000,	0.0000 },
//                { this,	"CA",	"CA",	"NA",	"HA*",	SpeciesImproper::Cos3,	0.0000,	9.2000,	0.0000,	0.0000 },
//                { this,	"CA",	"CA",	"NA",	"CT",	SpeciesImproper::Cos3,	0.0000,	8.3700,	0.0000,	0.0000 },
//                { this,	"CR",	"CW",	"NA",	"HA*",	SpeciesImproper::Cos3,	0.0000,	8.3700,	0.0000,	0.0000 },
//                { this,	"CR",	"CW",	"NA",	"C1",	SpeciesImproper::Cos3,	0.0000,	8.3700,	0.0000,	0.0000 },
//                { this,	"NA",	"NA",	"CR",	"HA*",	SpeciesImproper::Cos3,	0.0000,	9.2000,	0.0000,	0.0000 },
//                { this,	"NA",	"CW",	"CW",	"HA*",	SpeciesImproper::Cos3,	0.0000,	9.2000,	0.0000,	0.0000 }
	};
}

Forcefield_CLDP2010::~Forcefield_CLDP2010()
{
}

/*
 * Definition
 */

// Return name of Forcefield
const char* Forcefield_CLDP2010::name() const
{
        return "Molecular force field for ionic liquids (version 2010/09/16)";
}

// Return description for Forcefield
const char* Forcefield_CLDP2010::description() const
{
        return "Molecular force field for ionic liquids (version 2010/09/16); This forcefield is the work of K. Shimizu, A. Podgorsek, F. Hammami, L. Gontrani, J. Canongia Lopes, and A. Padua; Questions regarding the forcefield itself to: agilio.padua@univ-bpclermont.fr; Latest version donwloadable at: http://therm10.univ-bpclermont.fr/~apadua/";
}

// Return short-range interaction style for AtomTypes
Forcefield::ShortRangeType Forcefield_CLDP2010::shortRangeType() const
{
        return Forcefield::LennardJonesType;
}
