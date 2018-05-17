/*
	*** Module Widget
	*** src/modules/neutronsq/gui/modulewidget.h
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

#ifndef DUQ_NEUTRONSQMODULEWIDGET_H
#define DUQ_NEUTRONSQMODULEWIDGET_H

#include "modules/neutronsq/gui/ui_modulewidget.h"
#include "gui/modulewidget.h"

// Forward Declarations
class DUQ;
class Module;
class NeutronSQModule;
class PartialSet;
class UChromaViewWidget;

// Module Widget
class NeutronSQModuleWidget : public ModuleWidget
{
	// All Qt declarations derived from QObject must include this macro
	Q_OBJECT

	private:
	// Associated Module
	NeutronSQModule* module_;
	// UChromaViews contained within this widget
	UChromaViewWidget* partialGRGraph_, *partialSQGraph_, *totalGRGraph_, *totalSQGraph_;
	// Reference to DUQ
	DUQ& dUQ_;

	public:
	// Constructor / Destructor
	NeutronSQModuleWidget(QWidget* parent, Module* module, DUQ& dUQ);
	~NeutronSQModuleWidget();
	// Main form declaration
	Ui::NeutronSQModuleWidget ui;
	// Update controls within widget
	void updateControls();
	// Disable sensitive controls within widget, ready for main code to run
	void disableSensitiveControls();
	// Enable sensitive controls within widget, ready for main code to run
	void enableSensitiveControls();


	/*
	 * ModuleWidget Implementations
	 */
	public:
	// Write widget state through specified LineParser
	bool writeState(LineParser& parser);
	// Read widget state through specified LineParser
	bool readState(LineParser& parser);


	/*
	 * Widgets / Functions
	 */
	private:
	// Set data targets in graphs
	void setGraphDataTargets(NeutronSQModule* module);
};

#endif
