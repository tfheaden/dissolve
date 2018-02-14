/*
	*** Simple Ordered Pointer w/Data Array Class
	*** src/templates/orderedpointerdataarray.h
	Copyright T. Youngs 2012-2018

	This file is part of dUQ.

	dUQ is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	dUQ is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with dUQ.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DUQ_ORDEREDPOINTERDATAARRAY_H
#define DUQ_ORDEREDPOINTERDATAARRAY_H

#include <stdlib.h>
#include <stdio.h>
#include "base/messenger.h"

// Forward Declarations
/* none */

/*
 * Simple Orderd Pointer w/Data Array Class
 * A contiguous pointer array class designed to be as lightweight as possible.
 * Each pointer may exist in the list only once.
 * Order of pointers, based on their increasing value, is maintained.
 */
template <class T> class OrderedPointerArray
{
	public:
	// Constructor
	OrderedPointerArray<T>()
	{
		nItems_ = 0;
		arraySize_ = 0;
		items_ = NULL;
	}
	// Destructor
	~OrderedPointerArray()
	{
		clear();
	}
	// Copy Constructor
	OrderedPointerArray<T>(const OrderedPointerArray<T>& source)
	{
		nItems_ = 0;
		arraySize_ = 0;
		items_ = NULL;
		(*this) = source;
	}
	// Assignment operator
	OrderedPointerArray<T>& operator=(const OrderedPointerArray<T>& source)
	{
		// Clear any current data in the list...
		clear();

		initialise(source.nItems());

		for (int n=0; n<source.nItems(); ++n) add(source.value(n));

		return *this;
	}
	// Element access operator
	T* operator[](int index)
	{
#ifdef CHECKS
		if ((index < 0) || (index >= nItems_))
		{
			Messenger::error("OrderedPointerArray<T>::operator[](%i) - Array index out of bounds (%i items in array).\n", index, nItems_);
			return NULL;
		}
#endif
		return items_[index];
	}


	/*
	 * Basic Data
	 */
	protected:
	// Array of items
	T** items_;
	// Size of items array
	int arraySize_;
	// Number of items in array
	int nItems_;

	public:
	// Clear the list
	void clear()
	{
		if (items_) delete[] items_;
		items_ = NULL;
		nItems_ = 0;
		arraySize_ = 0;
	}
	// Create empty list of specified size
	void initialise(int size)
	{
		clear();

		arraySize_ = size;

		if (size == 0) items_ = NULL;
		else items_ = new T*[arraySize_];

		for (int n=0; n<arraySize_; ++n) items_[n] = NULL;
	}
	// Returns the number of items in the list
	int nItems() const
	{
		return nItems_;
	}
	// Const pointer access
	T* pointer(int index) const
	{
#ifdef CHECKS
		if ((index < 0) || (index >= nItems_))
		{
			Messenger::error("OrderedPointerArray<T>::index(%i) - Array index out of bounds (%i items in array).\n", index, nItems_);
			return NULL;
		}
#endif

		return items_[index];
	}
	// Return pointer array
	T** items()
	{
		return items_;
	}


	/*
	 * Item Management
	 */
	public:
	// Add an item to the list
	void add(T* ptr)
	{
		// If there is no space for the new item we must resize the array
		if (nItems_ == arraySize_)
		{
			// Create a new array with some more item space
			T** newArray = new T*[arraySize_+2];

			// Copy the existing pointers over
			// TODO Add the new pointer at the same time as copying the existing ones
			for (int n=0; n<arraySize_; ++n) newArray[n] = items_[n];

			// Delete the old array
			if (items_) delete[] items_;

			// Set the new array
			items_ = newArray;
			arraySize_ += 2;
		}

		// Loop over pointers, searching for the first existing item that has a higher pointer address.
		// If we find the same pointer address in the list, return immediately.
		int insertAt;
		for (insertAt = 0; insertAt<nItems_; ++insertAt)
		{
			if (items_[insertAt] == ptr) return;
			if (items_[insertAt] > ptr) break;
		}
		if (insertAt < nItems_)
		{
			// Shuffle all pointers from the insertion point up in the list, towards higher indices
			for (int n=startIndex; n<nItems_; ++n) items_[n+1] = items_[n];
		}

		// Now put the new item into the 'insertAt' position
		items_[insertAt] = ptr;
		++nItems_;
	}
	// Remove an item from the array, leaving the remaining items contiguous in memory
	void remove(T* ptr)
	{
		// Step through the items until we find the specified pointer
		for (int n=0; n<nItems_; ++n)
		{
			if (items_[n] == ptr)
			{
				// Found it. If it is at the end of the list, just decrease nItems_.
				// If not, switch the current last item into this slot, and then decrease nItems_.
				if (n < (nItems_-1))
				{
					items_[n] = items_[nItems_-1];
					items_[nItems_-1] = NULL;
				}
				--nItems_;
				return;
			}
		}
		Messenger::print("OrderedPointerArray<T>::remove(%p) - Couldn't find pointer in array.\n", ptr);
	}
	// Return array index of pointer within the list
	int indexOf(T* ptr) const
	{
		for (int n=0; n<nItems_; ++n) if (items_[n] == ptr) return n;

		return -1;
	}
};

#endif