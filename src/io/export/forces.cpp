// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2021 Team Dissolve and contributors

#include "io/export/forces.h"
#include "base/lineparser.h"
#include "base/sysfunc.h"
#include "classes/atomtype.h"
#include "classes/box.h"
#include "classes/configuration.h"
#include "classes/speciesatom.h"
#include "data/atomicmasses.h"

ForceExportFileFormat::ForceExportFileFormat(std::string_view filename, ForceExportFormat format)
    : FileAndFormat(filename, format)
{
}

/*
 * Format Access
 */

// Return enum options for ForceExportFormat
EnumOptions<ForceExportFileFormat::ForceExportFormat> ForceExportFileFormat::forceExportFormats()
{
    return EnumOptions<ForceExportFileFormat::ForceExportFormat>(
        "ForceExportFileFormat", {{ForceExportFileFormat::SimpleForces, "simple", "Simple Free-Formatted Forces"}});
}

// Return number of available formats
int ForceExportFileFormat::nFormats() const { return ForceExportFileFormat::nForceExportFormats; }

// Return format keyword for supplied index
std::string_view ForceExportFileFormat::formatKeyword(int id) const { return forceExportFormats().keywordByIndex(id); }

// Return description string for supplied index
std::string_view ForceExportFileFormat::formatDescription(int id) const { return forceExportFormats().descriptionByIndex(id); }

// Return current format as ForceExportFormat
ForceExportFileFormat::ForceExportFormat ForceExportFileFormat::forceFormat() const
{
    return (ForceExportFileFormat::ForceExportFormat)format_;
}

/*
 * Export Functions
 */

// Export simple forces
bool ForceExportFileFormat::exportSimple(LineParser &parser, const Array<double> &fx, const Array<double> &fy,
                                         const Array<double> &fz)
{
    if (!parser.writeLine("# Atom        FX            FY            FZ"))
        return false;

    if (!parser.writeLineF("{}\n", fx.nItems()))
        return false;

    for (auto n = 0; n < fx.nItems(); ++n)
        if (!parser.writeLineF("  {:10d}  {:15.8e}  {:15.8e}  {:15.8e}\n", n + 1, fx.at(n), fy.at(n), fz.at(n)))
            return false;

    return true;
}

// Export forces using current filename and format
bool ForceExportFileFormat::exportData(const Array<double> &fx, const Array<double> &fy, const Array<double> &fz)
{
    // Open the file
    LineParser parser;
    if (!parser.openOutput(filename_))
    {
        parser.closeFiles();
        return false;
    }

    // Write data
    auto result = false;
    if (forceFormat() == ForceExportFileFormat::SimpleForces)
        result = exportSimple(parser, fx, fy, fz);
    else
    {
        Messenger::error("Unrecognised force format.\nKnown formats are:\n");
        printAvailableFormats();
    }

    return result;
}
