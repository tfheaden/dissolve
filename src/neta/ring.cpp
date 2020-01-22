/*
	*** NETA Ring Node
	*** src/neta/ring.cpp
	Copyright T. Youngs 2019-2020

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

#include "neta/ring.h"
#include "data/ffatomtype.h"
#include "classes/speciesatom.h"
#include "templates/dynamicarray.h"
#include "templates/refdatalist.h"

// Constructor
NETARingNode::NETARingNode(NETADefinition* parent) : NETANode(parent, NETANode::RingNode)
{
	repeatCount_ = 1;
	repeatCountOperator_ = NETANode::EqualTo;
	sizeValue_ = -1;
	sizeValueOperator_ = NETANode::EqualTo;
	noRings_ = false;
}

// Destructor
NETARingNode::~NETARingNode()
{
}

/*
 * Modifiers
 */

// Return enum options for NETARingModifiers
EnumOptions<NETARingNode::NETARingModifier> NETARingNode::modifiers()
{
	static EnumOptionsList ModifierOptions = EnumOptionsList() <<
		EnumOption(SizeModifier,		"size") <<
		EnumOption(RepeatRingModifier,		"n");
	
	static EnumOptions<NETARingNode::NETARingModifier> options("RingModifier", ModifierOptions);

	return options;
}

// Return whether the specified modifier is valid for this node
bool NETARingNode::isValidModifier(const char* s) const
{
	return (modifiers().isValid(s));
}

// Set value and comparator for specified modifier
bool NETARingNode::setModifier(const char* modifier, ComparisonOperator op, int value)
{
	// Check that the supplied index is valid
	if (!modifiers().isValid(modifier)) return Messenger::error("Invalid modifier '%s' passed to NETARingNode.\n", modifier);

	switch (modifiers().enumeration(modifier))
	{
		case (NETARingNode::SizeModifier):
			sizeValue_ = value;
			sizeValueOperator_ = op;
			break;
		case (NETARingNode::RepeatRingModifier):
			repeatCount_ = value;
			repeatCountOperator_ = op;
			break;
		default:
			return Messenger::error("Don't know how to handle modifier '%s' in ring node.\n", modifier);
	}

	return true;
}

/*
 * Flags
 */


// Return enum options for NETARingFlags
EnumOptions<NETARingNode::NETARingFlag> NETARingNode::flags()
{
	static EnumOptionsList FlagOptions = EnumOptionsList() <<
		EnumOption(NoRingsFlag,			"none");

	static EnumOptions<NETARingNode::NETARingFlag> options("RingFlag", FlagOptions);

	return options;
}

// Return whether the specified flag is valid for this node
bool NETARingNode::isValidFlag(const char* s) const
{
	return (flags().isValid(s));
}

// Set specified flag
bool NETARingNode::setFlag(const char* flag, bool state)
{
	// Check that the supplied index is valid
	if (!flags().isValid(flag)) return Messenger::error("Invalid flag '%s' passed to NETARingNode.\n", flag);

	switch (flags().enumeration(flag))
	{
		case (NETARingNode::NoRingsFlag):
			noRings_ = true;
			break;
		default:
			return Messenger::error("Don't know how to handle flag '%s' in connection node.\n", flag);
	}

	return true;
}

/*
 * Scoring
 */

// Locate rings in which the specified atom is involved
void NETARingNode::findRings(const SpeciesAtom* currentAtom, List<SpeciesRing>& rings, PointerArray<const SpeciesAtom>& path, const int minSize, const int maxSize) const
{
	// Check whether the path is already at the maximum size - if so, return immediately.
	if (path.nItems() == maxSize) return;

	// Add the current atom to the path
	path.append(currentAtom);

	// Loop over bonds to the atom
	const SpeciesAtom* j;
	SpeciesRing* ring;
	const PointerArray<SpeciesBond>& bonds = currentAtom->bonds();
	for (int n=0; n<bonds.nItems(); ++n)
	{
		/*
		 * Get the partner atom and compare to first atom in the current path.
		 * If it is the currentAtom then we have found a cyclic route back to the originating atom.
		 * If not, check whether the atom is already elsewhere in the path - if so, continue with the next bond.
		 */
		j = bonds.at(n)->partner(currentAtom);
		if ((path.nItems() >= minSize) && (j == path.at(0)))
		{
			// Special case - if NotEqualTo was specified as the comparison operator, check that against the maximum size
			if ((sizeValueOperator_ == NETANode::NotEqualTo) && (path.nItems() == maxSize)) continue;

			// Add new ring
			ring = rings.add();
			ring->setAtoms(path);

			// Continue with the next bond, unless the "None" modifier has been set
			if (noRings_) break;
			else continue;
		}
		else if (path.sniatnoc(j)) continue;

		// The current atom j is not in the path, so recurse
		findRings(j, rings, path, minSize, maxSize);
	}

	// Remove current atom from the path
	path.removeLast();
}

// Evaluate the node and return its score
int NETARingNode::score(const SpeciesAtom* i, RefList<const SpeciesAtom>& matchPath) const
{
// 	printf("I AM THE RING - matchPath size = %i:\n", matchPath.nItems());
// 	RefListIterator<const SpeciesAtom> matchIterator(matchPath);
// 	while (const SpeciesAtom* iii = matchIterator.iterate()) printf("   -- %p %i %s\n", iii, iii->userIndex(), iii->element()->symbol());
// 	printf("SITTING ON SPECIESATOM %i (%s)\n", i->userIndex(), i->element()->symbol());

	// Generate array of rings of specified size that the atom 'i' is present in
	List<SpeciesRing> rings;
	PointerArray<const SpeciesAtom> ringPath;
	if (sizeValue_ == -1) findRings(i, rings, ringPath, 3, 6);
	else if (sizeValueOperator_ == NETANode::EqualTo) findRings(i, rings, ringPath, sizeValue_, sizeValue_);
	else if (sizeValueOperator_ == NETANode::LessThan) findRings(i, rings, ringPath, 3, sizeValue_-1);
	else if (sizeValueOperator_ == NETANode::LessThanEqualTo) findRings(i, rings, ringPath, 3, sizeValue_);
	else if (sizeValueOperator_ == NETANode::GreaterThan) findRings(i, rings, ringPath, sizeValue_+1, 99);
	else if (sizeValueOperator_ == NETANode::GreaterThanEqualTo) findRings(i, rings, ringPath, sizeValue_, 99);
	else findRings(i, rings, ringPath, 3, 99);

	// Was the "None" modifier given?
	if (noRings_) return (rings.nItems() == 0 ? 1 : NETANode::NoMatch);

	// Prune rings for duplicates
	ListIterator<SpeciesRing> ringIterator(rings);
	while (SpeciesRing* ring = ringIterator.iterate())
	{
		// Check this ring against others in the list - if we find a duplicate, we can remove it and then move on with the next ring.
		for (SpeciesRing* other = ring->next(); other != NULL; other = other->next())
		{
			if ((*ring) == (*other))
			{
				rings.remove(other);
				break;
			}
		}
	}
// 	ringIterator.restart();
// 	while (SpeciesRing* ring = ringIterator.iterate()) ring->print();

	// Loop over rings
	int nMatches = 0, totalScore = 0, nodeScore;
	ringIterator.restart();
	while (SpeciesRing* ring = ringIterator.iterate())
	{
		// Copy the atoms in the ring into an array we can modify
		PointerArray<const SpeciesAtom> ringAtoms = ring->atoms();

		// Check through atoms in the ring - either in order or not - to see if the ring matches
		if (false)
		{
			// Ordered search
		}
		else
		{
			// Disordered search - try to match the branch definition against this ring, in any order (provide all atoms in the ring at once)
			RefList<const SpeciesAtom> ringAtoms;
			for (int n=0; n<ring->size(); ++n) ringAtoms.append(ring->atom(n));

			const SpeciesAtom* matchedAtom;
			ListIterator<NETANode> branchIterator(branch_);
			while (NETANode* node = branchIterator.iterate())
			{
				nodeScore = node->score(NULL, ringAtoms);
				if (nodeScore == NETANode::NoMatch) break;

				// Match found
				totalScore += nodeScore;
			}

			// If we didn't find a match for the ring, exit the loop now
			if (nodeScore == NETANode::NoMatch) break;

			++nMatches;

			// Don't match more than we need to - check the repeatCount
			if (compareValues(nMatches, repeatCountOperator_, repeatCount_)) break;
		}
	}

	// Did we find the required number of ring matches?
	if (!compareValues(nMatches, repeatCountOperator_, repeatCount_)) return NETANode::NoMatch;

	return totalScore;
}
