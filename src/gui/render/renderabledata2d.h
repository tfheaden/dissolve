// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2021 Team Dissolve and contributors

#pragma once

#include "gui/render/renderable.h"
#include "math/data2d.h"

// Forward Declarations
class Axes;

// Renderable for Data2D
class RenderableData2D : public Renderable
{
    public:
    RenderableData2D(const Data2D *source, std::string_view objectTag);
    ~RenderableData2D();

    /*
     * Data
     */
    private:
    // Source data
    const Data2D *source_;

    private:
    // Return whether a valid data source is available (attempting to set it if not)
    bool validateDataSource();
    // Invalidate the current data source
    void invalidateDataSource();

    public:
    // Return version of data
    int dataVersion();

    /*
     * Transform / Limits
     */
    private:
    // Transformed data
    Data2D transformedData_;

    protected:
    // Transform data according to current settings
    void transformValues();
    // Return reference to transformed data
    const Data2D &transformedData();

    /*
     * Rendering Primitives
     */

    private:
    // Create line strip primitive
    void constructLine(const std::vector<double> &displayXAbscissa, const std::vector<double> &displayYAbscissa,
                       const Array2D<double> &displayValues, const Axes &axes, const ColourDefinition &colourDefinition);

    protected:
    // Recreate necessary primitives / primitive assemblies for the data
    void recreatePrimitives(const View &view, const ColourDefinition &colourDefinition);
    // Send primitives for rendering
    const void sendToGL(const double pixelScaling);

    /*
     * Style
     */
    public:
    // Display Styles enum
    enum Data2DDisplayStyle
    {
        LinesStyle,
        nData2DDisplayStyles
    };
    // Return EnumOptions for Data2DDisplayStyle
    static EnumOptions<Data2DDisplayStyle> data2DDisplayStyles();

    private:
    // Display style for the renderable
    Data2DDisplayStyle displayStyle_;

    public:
    // Set display style for renderable
    void setDisplayStyle(Data2DDisplayStyle displayStyle);
    // Return display style for the renderable
    Data2DDisplayStyle displayStyle() const;

    /*
     * Style I/O
     */
    public:
    // Data2DStyle Keywords Enum
    enum Data2DStyleKeyword
    {
        DisplayKeyword,  /* 'Display' - General display style for renderable */
        EndStyleKeyword, /* 'EndStyle' - End of Style block */
        nData2DStyleKeywords
    };
    // Return enum option info for RenderableKeyword
    static EnumOptions<RenderableData2D::Data2DStyleKeyword> data2DStyleKeywords();
    // Write style information
    bool writeStyleBlock(LineParser &parser, int indentLevel = 0) const;
    // Read style information
    bool readStyleBlock(LineParser &parser);
};
