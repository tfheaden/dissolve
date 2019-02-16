/*
	*** Rendering Primitive Instance
	*** src/gui/viewer/render/primitiveinstance.cpp
	Copyright T. Youngs 2013-2019

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

#include "gui/viewer/render/primitiveinstance.h"
#include <string.h>

// Static members
PrimitiveInstance::InstanceType PrimitiveInstance::globalInstanceType_ = PrimitiveInstance::VBOInstance;

// Constructor
PrimitiveInstance::PrimitiveInstance() : ListItem<PrimitiveInstance>()
{
	// Private variables
	context_ = NULL;
	type_ = PrimitiveInstance::ListInstance;
	listObject_ = 0;
	vboVertexObject_ = 0;
	vboIndexObject_ = 0;
}

// Return global instance type to use
PrimitiveInstance::InstanceType PrimitiveInstance::globalInstanceType()
{
	return globalInstanceType_;
}

// Set global instance type to use
void PrimitiveInstance::setGlobalInstanceType(PrimitiveInstance::InstanceType instanceType)
{
	globalInstanceType_ = instanceType;
}

// Return context to which primitive instance is associated
const QOpenGLContext* PrimitiveInstance::context()
{
	return context_;
}

// Return type of instance
PrimitiveInstance::InstanceType PrimitiveInstance::type() const
{
	return type_;
}

// Set display list data
void PrimitiveInstance::setDisplayList(const QOpenGLContext* context, GLuint listObject)
{
	context_ = context;
	type_ = PrimitiveInstance::ListInstance;
	listObject_ = listObject;
}

// Set vbo object data
void PrimitiveInstance::setVBO(const QOpenGLContext* context, GLuint vertexObject, GLuint indexObject)
{
	context_ = context;
	type_ = PrimitiveInstance::VBOInstance;
	vboVertexObject_ = vertexObject;
	vboIndexObject_ = indexObject;
}

// Return display list object for instance
GLuint PrimitiveInstance::listObject() const
{
	return listObject_;
}

// Return VBO ID of vertex array for instance
GLuint PrimitiveInstance::vboVertexObject() const
{
	return vboVertexObject_;
}

// Return VBO ID of index array for instance
GLuint PrimitiveInstance::vboIndexObject() const
{
	return vboIndexObject_;
}
