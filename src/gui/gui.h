/*
	*** dUQ Main Window
	*** src/gui/gui.h
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

#ifndef DUQ_MAINWINDOW_H
#define DUQ_MAINWINDOW_H

#include "gui/ui_gui.h"
#include "gui/widgets/subwindow.h"
#include "gui/outputhandler.hui"
#include "gui/thread.hui"
#include "gui/maintab.h"
#include "base/charstring.h"
#include "templates/list.h"

// Forward Declarations
class BrowserWidget;
class Configuration;
class DUQ;
class QMdiSubWindow;

class DUQWindow : public QMainWindow
{
	// All Qt declarations must include this macro
	Q_OBJECT

	public:
	// Constructor / Destructor
	DUQWindow(DUQ& duq);
	~DUQWindow();
	// Main form declaration
	Ui::DUQWindow ui;

	protected:
	void closeEvent(QCloseEvent* event);
	void resizeEvent(QResizeEvent* event);


	/*
	 * DUQ Reference
	 */
	private:
	// DUQ reference
	DUQ& duq_;
	// Whether any data has been modified in the GUI
	bool modified_;

	public slots:
	// Flag that data has been modified via the GUI
	void setModified();

	public:
	// Return reference to DUQ
	DUQ& duq();


	/*
	 * File
	 */
	public:
	// Open specified input file
	bool openFile(const char* inputFile, bool ignoreRestartFile, bool ignoreLayoutFile);


	/*
	 * Update Functions
	 */
	private:
	// Whether window is currently refreshing
	bool refreshing_;
	// Whether window has been shown
	bool shown_;
	// Output handler for messaging in GUI
	GUIOutputHandler outputHandler_;

	public slots:
	// Refresh all controls
	void updateControls();
	// Update status
	void updateStatus();
	// Link output handler in to the Messenger
	void addOutputHandler();


	/*
	 * Main Menu
	 */
	private slots:
	// Session
	void on_SessionOpenAction_triggered(bool checked);
	void on_SessionSaveAction_triggered(bool checked);
	void on_SessionQuitAction_triggered(bool checked);
	// Simulation
	// Workspace
	void on_WorkspaceAddNewAction_triggered(bool checked);
	void addWidgetToCurrentWorkspace(bool checked);
	void removeWidgetFromCurrentWorkspace(QString windowTitle);

	public:
	// Update menu items (after change in Modules etc.)
	void updateMenuItems();
	// Create Module menu with specified QAction as parent
	void createModuleMenu(QMenu* parent);


	/*
	 * Run Control
	 */
	public:
	// dUQ State Enum
	enum DUQState {
		RunningState,		/* dUQ is currently running in the GUI */
		StoppedState,		/* dUQ is currently stopped in the GUI */
		MonitoringState,	/* dUQ is running in the background, and we are monitoring it */
		nDUQStates
	};

	private:
	// Thread controller
	DUQThreadController threadController_;
	// Current state of dUQ
	DUQState duqState_;

	public:
	// Return current state of dUQ
	DUQState duqState() const;

	private slots: 
	void on_ControlSetUpButton_clicked(bool checked);
	void on_ControlRunButton_clicked(bool checked);
	void on_ControlStepButton_clicked(bool checked);
	void on_ControlStepFiveButton_clicked(bool checked);
	void on_ControlPauseButton_clicked(bool checked);
	void on_ControlReloadButton_clicked(bool checked);

	public slots:
	// Set widgets ready for the main code to be run
	void setWidgetsForRun();
	// Set widgets after the main code has been run
	void setWidgetsAfterRun();
	// All iterations requested are complete
	void iterationsComplete();

	signals:
	void iterate(int);
	void stopIterating();


	/*
	 * Tab Management
	 */
	private:
	// List of all displayed tabs
	List<MainTab> tabs_;

	private slots:
	void on_MainTabs_currentChanged(int index);
	void mainTabsDoubleClicked(int index);

	private:
	// Clear all tabs
	void clearAllTabs();
	// Add core tabs
	void addCoreTabs();
	// Add tab for specified Configuration target
	void addConfigurationTab(Configuration* cfg);
	// Add processing workspace
	void addProcessingTab();
	// Add on an empty workspace tab
	MainTab* addWorkspaceTab(const char* title);
	// Find tab with title specified
	MainTab* findTab(const char* title);
	// Find tab with specified page widget
	MainTab* findTab(QWidget* page);
	// Return current tab
	MainTab* currentTab();
	// Make specified tab the current one
	void setCurrentTab(MainTab* tab);
	// Make specified tab the current one (by index)
	void setCurrentTab(int tabIndex);

	public:
	// Add tab for specified Module target
	void addModuleTab(Module* module);

	public slots:
	// Remove tab containing the specified page widget, as it has been deleted
	void removeDeletedTab(QWidget* page);


	/*
	 * SubWidgets
	 */
	private:
	// Create named SubWidget
	SubWidget* createSubWidget(const char* widgetName, const char* title);

	public slots:
	// Shift ModuleReference up in its list
	void shiftModuleUp(void* moduleReference);
	// Shift ModuleReference down in its list
	void shiftModuleDown(void* moduleReference);
	// Remove Module(Reference)
	void removeModule(void* moduleReference);


	/*
	 * Window State
	 */
	private:
	// Filename containing current window layout
	CharString windowLayoutFilename_;

	public:
	// Save current window layout
	bool saveWindowLayout();
	// Load window state
	bool loadWindowLayout();
};

#endif
