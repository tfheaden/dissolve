// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2021 Team Dissolve and contributors

#include "modules/calculate_avgmol/avgmol.h"
#include "modules/calculate_avgmol/gui/modulewidget.h"

CalculateAvgMolModuleWidget::CalculateAvgMolModuleWidget(QWidget *parent, CalculateAvgMolModule *module)
    : ModuleWidget(parent), module_(module)
{
    // Set up user interface
    ui_.setupUi(this);

    ui_.SpeciesView->setSpecies(&module_->averageSpecies());

    refreshing_ = false;
}

/*
 * UI
 */

// Update controls within widget
void CalculateAvgMolModuleWidget::updateControls(int flags) {}

// Disable sensitive controls within widget
void CalculateAvgMolModuleWidget::disableSensitiveControls() {}

// Enable sensitive controls within widget
void CalculateAvgMolModuleWidget::enableSensitiveControls() {}

/*
 * Widgets / Functions
 */
