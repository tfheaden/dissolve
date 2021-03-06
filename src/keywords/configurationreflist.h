// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2021 Team Dissolve and contributors

#pragma once

#include "keywords/data.h"
#include "templates/reflist.h"
#include <optional>

// Forward Declarations
class Configuration;

// Keyword with Configuration RefList Data
class ConfigurationRefListKeyword : public KeywordData<RefList<Configuration> &>
{
    public:
    ConfigurationRefListKeyword(RefList<Configuration> &references, int maxListSize);
    ~ConfigurationRefListKeyword();

    /*
     * Data
     */
    private:
    // Maximum number of Configurations to allow in the list (or -1 for any number)
    int maxListSize_;

    protected:
    // Determine whether current data is 'empty', and should be considered as 'not set'
    bool isDataEmpty() const;

    public:
    // Return maximum number of Configurations to allow in the list
    int maxListSize() const;

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
     * Object Management
     */
    protected:
    // Prune any references to the supplied Configuration in the contained data
    void removeReferencesTo(Configuration *cfg);
};
