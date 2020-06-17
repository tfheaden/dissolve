/*
	*** Keyword Widget - Species RefList
	*** src/gui/keywordwidgets/speciesreflist.h
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

#ifndef DISSOLVE_KEYWORDWIDGET_SPECIESREFLIST_H
#define DISSOLVE_KEYWORDWIDGET_SPECIESREFLIST_H

#include "gui/keywordwidgets/base.h"
#include "gui/keywordwidgets/dropdown.h"
#include "gui/keywordwidgets/ui_speciesreflist.h"
#include "keywords/speciesreflist.h"
#include <QWidget>

// Forward Declarations
class Species;

class SpeciesRefListKeywordWidget : public KeywordDropDown, public KeywordWidgetBase
{
	// All Qt declarations must include this macro
	Q_OBJECT

      public:
	// Constructor
	SpeciesRefListKeywordWidget(QWidget *parent, KeywordBase *keyword, const CoreData &coreData);

	/*
	 * Keyword
	 */
      private:
	// Associated keyword
	SpeciesRefListKeyword *keyword_;

	/*
	 * Widgets
	 */
      private:
	// Main form declaration
	Ui::SpeciesRefListWidget ui_;

      private:
	// Selection list update function
	void updateSelectionRow(int row, Species *sp, bool createItem);

      private slots:
	// List item changed
	void itemChanged(QListWidgetItem *item);

      signals:
	// Keyword value changed
	void keywordValueChanged(int flags);

	/*
	 * Update
	 */
      public:
	// Update value displayed in widget
	void updateValue();
	// Update widget values data based on keyword data
	void updateWidgetValues(const CoreData &coreData);
	// Update keyword data based on widget values
	void updateKeywordData();
	// Update summary text
	void updateSummaryText();
};

#endif
