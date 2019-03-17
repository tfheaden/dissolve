/*
	*** Species Widget - Functions 
	*** src/gui/viewer/specieswidget_funcs.cpp
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

#include "gui/viewer/specieswidget.h"
#include "gui/widgets/elementselector.hui"
#include "classes/empiricalformula.h"
#include "classes/species.h"

// Constructor
SpeciesWidget::SpeciesWidget(QWidget* parent) : QWidget(parent)
{
	// Set up our UI
	ui_.setupUi(this);

	// Create a button group for the interaction modes
	QButtonGroup* group = new QButtonGroup;
	group->addButton(ui_.InteractionViewButton);
	group->addButton(ui_.InteractionDrawButton);

	// Connect signals / slots
	connect(ui_.SpeciesView, SIGNAL(atomSelectionChanged()), this, SLOT(updateStatusBar()));
	connect(ui_.SpeciesView, SIGNAL(interactionModeChanged()), this, SLOT(updateStatusBar()));

	// Make sure our controls are consistent with the underlying viewer / data
	updateToolbar();
	updateStatusBar();
}

// Destructor
SpeciesWidget::~SpeciesWidget()
{
}

// Return contained SpeciesViewer
SpeciesViewer* SpeciesWidget::speciesViewer()
{
	return ui_.SpeciesView;
}

/*
 * Toolbar
 */

void SpeciesWidget::on_InteractionViewButton_clicked(bool checked)
{
	if (checked) speciesViewer()->setInteractionMode(SpeciesViewer::DefaultInteraction);
}

void SpeciesWidget::on_InteractionDrawButton_clicked(bool checked)
{
	if (checked) speciesViewer()->setInteractionMode(SpeciesViewer::DrawInteraction);
}

void SpeciesWidget::on_InteractionDrawElementButton_clicked(bool checked)
{
	// Select a new element for drawing
	bool ok;
	Element* newElement = ElementSelector::getElement(this, "Choose Element", "Select element to use for drawn atoms", speciesViewer()->drawElement(), &ok);
	if (!ok) return;

	speciesViewer()->setDrawElement(newElement);

	updateToolbar();
}

void SpeciesWidget::on_ViewResetButton_clicked(bool checked)
{
	speciesViewer()->view().showAllData();
	speciesViewer()->view().resetViewMatrix();

	speciesViewer()->postRedisplay();
}

void SpeciesWidget::on_ViewAxesVisibleButton_clicked(bool checked)
{
	speciesViewer()->setAxesVisible(checked);

	speciesViewer()->postRedisplay();
}

/*
 * Signals / Slots
 */

// Update toolbar to reflect current viewer state
void SpeciesWidget::updateToolbar()
{
	// Set current interaction mode
	switch (speciesViewer()->interactionMode())
	{
		case (SpeciesViewer::DefaultInteraction):
			ui_.InteractionViewButton->setChecked(true);
			break;
		case (SpeciesViewer::DrawInteraction):
			ui_.InteractionDrawButton->setChecked(true);
			break;
	}

	// Set drawing element symbol
	ui_.InteractionDrawElementButton->setText(speciesViewer()->drawElement()->symbol());

	// Set checkable buttons
	ui_.ViewAxesVisibleButton->setChecked(speciesViewer()->axesVisible());
}

// Update status bar
void SpeciesWidget::updateStatusBar()
{
	// Get displayed Species
	const Species* sp = speciesViewer()->species();

	// Set interaction mode text
	ui_.ModeLabel->setText(speciesViewer()->interactionModeText());

	// Set / update empirical formula for the Species and its current atom selection
	ui_.FormulaLabel->setText(sp ? EmpiricalFormula::formula(sp, true) : "--");
	ui_.SelectionLabel->setText(sp && (sp->nSelectedAtoms() > 0) ? EmpiricalFormula::formula(sp->selectedAtoms(), true) : "--");
}