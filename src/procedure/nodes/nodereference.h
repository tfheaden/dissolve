// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2021 Team Dissolve and contributors

#pragma once

#include "procedure/nodes/node.h"
#include "templates/listitem.h"

// Forward Declarations
class CoreData;
class AnalyseModule;

// Procedure Node Reference
class ProcedureNodeReference : public ListItem<ProcedureNodeReference>
{
    public:
    ProcedureNodeReference(ProcedureNode *node = nullptr);
    ~ProcedureNodeReference();

    /*
     * Data
     */
    private:
    // Target node
    ProcedureNode *node_;
    // Parent AnalyseModule owning the node (if relevant)
    AnalyseModule *analyseModuleParent_;
    // Allowable node types (when reading / setting)
    bool allowedTypes_[ProcedureNode::nNodeTypes];

    public:
    // Return target node
    ProcedureNode *node();
    // Return node type
    ProcedureNode::NodeType type() const;
    // Add allowable node type
    void addAllowableNodeType(ProcedureNode::NodeType nt);
    // Allow all node types
    void setAllowAllNodeTypes();
    // Return if node pointer is NULL
    bool isNull() const;

    /*
     * Operators
     */
    public:
    void operator=(ProcedureNode *node);
    void operator=(const ProcedureNodeReference &nodeRef);

    /*
     * Read / Write
     */
    public:
    // Read structure from specified LineParser
    bool read(LineParser &parser, int startArg, CoreData &coreData, const Procedure *procedure);
    // Write structure to specified LineParser
    bool write(LineParser &parser, std::string_view prefix);
};
