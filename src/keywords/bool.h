// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2021 Team Dissolve and contributors

#pragma once

#include "keywords/data.h"

// Forward Declarations
/* none */

// Keyword with Bool Data
class BoolKeyword : public KeywordData<bool>
{
    public:
    BoolKeyword(bool value);
    ~BoolKeyword();

    /*
     * Arguments
     */
    public:
    // Return minimum number of arguments accepted
    int minArguments() const;
    // Return maximum number of arguments accepted
    int maxArguments() const;
    // Parse arguments from supplied LineParser, starting at given argument offset
    bool read(LineParser &parser, int startArg, CoreData &coreData);
    // Write keyword data to specified LineParser
    bool write(LineParser &parser, std::string_view keywordName, std::string_view prefix);

    /*
     * Conversion
     */
    public:
    // Return value (as bool)
    bool asBool();
    // Return value (as int)
    int asInt();
    // Return value (as double)
    double asDouble();
    // Return value (as string)
    std::string asString();
};
