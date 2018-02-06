/*
	*** DropDown for Keyword Widget
	*** src/gui/keywordwidgets/dropdown.h
	Copyright T. Youngs 2012-2018

	This file is part of DUQ.

	DUQ is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	DUQ is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with DUQ.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DUQ_KEYWORDDROPDOWN_H
#define DUQ_KEYWORDDROPDOWN_H

#include "gui/keywordwidgets/ui_dropdown.h"
#include "gui/keywordwidgets/dropwidget.hui"
#include <QWidget>

// Forward Declarations
/* none */

class KeywordDropDown : public QWidget
{
	// All Qt declarations must include this macro
	Q_OBJECT

	public:
	// Constructor
	KeywordDropDown(QWidget* parent);
        // Main form declaration
        Ui::KeywordDropDownControlWidget ui;


	/*
	 * Drop Widget
	 */
	private:
	// Widget to display as the drop-down
	DropWidget dropWidget_;

	public:
	// Return the drop widget
	DropWidget* dropWidget();


	/*
	 * Signals / Slots
	 */
	private slots:
	void on_CallDropWidgetButton_clicked(bool checked);
	void dropWidgetHidden();


	/*
	 * Update
	 */
	protected:
	// Update keyword data based on widget values
	virtual void updateKeywordData() = 0;
};

#endif
