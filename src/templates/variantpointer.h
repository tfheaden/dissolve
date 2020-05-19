/*
    *** Variant Pointer
    *** src/templates/variantpointer.h
    Copyright T. Youngs 2013-2020

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

#pragma once

#include <QVariant>

// Simple class to convert between QVariant pointer (void*) and a custom class pointer
template <class A> class VariantPointer
{
    private:
    // Pointer to target class
    A *pointer_;

    public:
    VariantPointer(A *ptr) { pointer_ = ptr; }
    VariantPointer(QVariant variant) { pointer_ = (A *)variant.value<void *>(); }

    operator QVariant() { return QVariant::fromValue((void *)pointer_); }

    operator A *() { return dynamic_cast<A *>(pointer_); }
};
