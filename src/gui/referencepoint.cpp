// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2021 Team Dissolve and contributors

#include "gui/referencepoint.h"

ReferencePoint::ReferencePoint() : ListItem<ReferencePoint>() {}

ReferencePoint::~ReferencePoint() {}

// Set suffix for data items
void ReferencePoint::setSuffix(std::string_view suffix) { suffix_ = suffix; }

// Return suffix for data items
std::string_view ReferencePoint::suffix() const { return suffix_; }

// Set restart file from which the reference point data was read
void ReferencePoint::setRestartFile(std::string_view restartFile) { restartFile_ = restartFile; }

// Return restart file from which the reference point data was read
std::string_view ReferencePoint::restartFile() const { return restartFile_; }
