/*
	*** Renderable Group
	*** src/gui/viewer/render/renderablegroupmanager.cpp
	Copyright T. Youngs 2013-2019

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

#include "gui/viewer/render/renderablegroupmanager.h"
#include "gui/viewer/render/renderable.h"
#include "base/sysfunc.h"

// Constructor
RenderableGroupManager::RenderableGroupManager()
{
	stockColourUsageCount_.initialise(ColourDefinition::nStockColours);
	verticalShift_ = NoVerticalShift;
}

// Create named group, or return existing group by the same name
RenderableGroup* RenderableGroupManager::createGroup(const char* name)
{
	// Does a group with this name already exist?
	RenderableGroup* renderableGroup = group(name);
	if (!renderableGroup)
	{
		// No existing group, so must add a new one
		// First, find the StockColour with the lowest usage count
		int lowestId = 0;
		for (int colourId = 0; colourId < ColourDefinition::nStockColours; ++colourId)
		{
			if (stockColourUsageCount_[colourId] < stockColourUsageCount_[lowestId]) lowestId = colourId;
		}

		renderableGroup = new RenderableGroup(name, (ColourDefinition::StockColour) lowestId);
		groups_.own(renderableGroup);
		++stockColourUsageCount_[lowestId];

		// Need to update vertical shifts
		setRenderableGroupShifts();
	}

	return renderableGroup;
}

// Add Renderable to its specified group, creating / associating as necessary
RenderableGroup* RenderableGroupManager::addToGroup(Renderable* renderable, const char* groupName)
{
	// Check to see if the Renderable is already associated to a group...
	if (renderable->group())
	{
		if (DissolveSys::sameString(renderable->group()->name(), groupName))
		{
			Messenger::print("Renderable '%s' already associated to group '%s'...\n", renderable->name(), renderable->group()->name());
			return renderable->group();
		}

		// Remove it from the current group
		renderable->group()->removeRenderable(renderable);
	}

	// Create / retrieve the group
	RenderableGroup* renderableGroup = createGroup(groupName);

	// Add unique Renderable reference to the group
	renderableGroup->associateRenderable(renderable);
	renderable->setGroup(renderableGroup);

	return renderableGroup;
}

// Return named group, if it exists
RenderableGroup* RenderableGroupManager::group(const char* name)
{
	for (RenderableGroup* group = groups_.first(); group != NULL; group = group->next) if (DissolveSys::sameString(group->name(), name)) return group;
	return NULL;
}

// Return group for specified Renderable, if one has been assigned
RenderableGroup* RenderableGroupManager::group(Renderable* renderable)
{
	for (RenderableGroup* group = groups_.first(); group != NULL; group = group->next) if (group->usedByRenderable(renderable)) return group;
	return NULL;
}

// Remove Renderable from its specified group
void RenderableGroupManager::removeFromGroup(Renderable* renderable)
{
	// If no group is currently set in the Renderable, nothing more to do here
	RenderableGroup* renderableGroup = renderable->group();
	if (!renderableGroup) return;

	// Remove the Renderable from the group and nullify its pointer
	renderableGroup->removeRenderable(renderable);

	// If the group is now empty, we can delete it
	if (renderableGroup->isEmpty())
	{
		--stockColourUsageCount_[renderableGroup->stockColour()];
		groups_.remove(renderableGroup);
	}
}

// Vertical shifts
double VerticalShiftAmounts[] = { 0.0, 0.5, 1.0, 2.0 };

// Set vertical shifts for current RenderableGroups
void RenderableGroupManager::setRenderableGroupShifts()
{
	// Loop over RenderableGroups
	double verticalShift = 0.0;
	for (RenderableGroup* group = groups_.first(); group != NULL; group = group->next)
	{
		group->setVerticalShift(verticalShift_ > 0, verticalShift);

		// Increase shift amount for the next group
		verticalShift += VerticalShiftAmounts[verticalShift_];
	}
}

// Return colour definition for specified Renderable
const ColourDefinition& RenderableGroupManager::colourDefinition(Renderable* renderable)
{
	RenderableGroup* renderableGroup = group(renderable);
	return (renderableGroup ? renderableGroup->colour() : renderable->colour());
}

// Cycle vertical shift applied to RenderableGroups
void RenderableGroupManager::cycleVerticalShifts()
{
	verticalShift_ = (RenderableGroupManager::VerticalShift) ((verticalShift_+1)%nVerticalShifts);

	setRenderableGroupShifts();
}

// Set vertical shift applied to RenderableGroups
void RenderableGroupManager::setVerticalShift(VerticalShift shiftType)
{
	verticalShift_ = shiftType;

	setRenderableGroupShifts();
}

// Return current vertical shift type
RenderableGroupManager::VerticalShift RenderableGroupManager::verticalShift() const
{
	return verticalShift_;
}

// Remove all vertical shifts from RenderableGroups
void RenderableGroupManager::removeVerticalShifts()
{
	verticalShift_ = NoVerticalShift;

	setRenderableGroupShifts();
}
