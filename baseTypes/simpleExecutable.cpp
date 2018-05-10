/*
** Copyright (C) 2018 Martin Brain
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
** simpleExecutable.cpp
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 07/08/14
**
** The most simple executable implementation of bit-vectors.
** Limited in the ranges it supports but fast and suitable for reasoning.
**
*/


#include "symfpu/baseTypes/simpleExecutable.h"

#include <assert.h>
#include <math.h>
#include <fenv.h>

namespace symfpu {
  namespace simpleExecutable {


    // This would all be much easier if C++ allowed partial specialisation of member templates...

    template <>
    bool bitVector<int64_t>::isRepresentable (const bitWidthType w, const int64_t v) {
      uint64_t shiftSafe = *((uint64_t *)(&v));
      uint64_t top = (shiftSafe >> w);
      uint64_t signbit = shiftSafe & 0x8000000000000000;
      int64_t stop = *((int64_t *)(&top));
      return (signbit) ? (stop == bitVector<int64_t>::nOnes(bitVector<int64_t>::maxWidth() - w)) : (stop == 0LL);
    }

    template <>
    bool bitVector<uint64_t>::isRepresentable (const bitWidthType w, const uint64_t v) {
      uint64_t top = (v >> w);
      return (top == 0);
    }

    
    template <>
    uint64_t bitVector<uint64_t>::makeRepresentable (const bitWidthType w, const uint64_t v) {
      return v & bitVector<uint64_t>::nOnes(w);
    }

    template <>
    int64_t bitVector<int64_t>::makeRepresentable (const bitWidthType w, const int64_t v) {
      if (v <= ((1LL << (w - 1)) - 1) && (-(1LL << (w - 1)) <= v)) {
	return v;
      } else {
	return 0;
      }
    }

    
    template <>
    bitVector<int64_t> bitVector<int64_t>::maxValue (const bitWidthType &w) {
      PRECONDITION(w != 1);
      return bitVector<int64_t>(w, (1ULL << (w - 1)) - 1);
    }
    
    template <>
    bitVector<uint64_t> bitVector<uint64_t>::maxValue (const bitWidthType &w) {
      PRECONDITION(w != 1);
      return bitVector<uint64_t>(w, (1ULL << w) - 1);
    }

    template <>
    bitVector<int64_t> bitVector<int64_t>::minValue (const bitWidthType &w) {
      PRECONDITION(w != 1);
      return bitVector<int64_t>(w, -(1ULL << (w - 1)));
    }
    
    template <>
    bitVector<uint64_t> bitVector<uint64_t>::minValue (const bitWidthType &w) {
      return bitVector<uint64_t>(w, 0);
    }
    

    template <>
    bitVector<int64_t> bitVector<int64_t>::operator- (void) const {
      return bitVector<int64_t>(this->width, -this->value);
    }

    // Used in addition
    template <>
    bitVector<uint64_t> bitVector<uint64_t>::operator- (void) const {
      return bitVector<uint64_t>(this->width,
				 bitVector<uint64_t>::makeRepresentable(this->width, (~this->value) + 1));
    }

    template <>
    bitVector<uint64_t> bitVector<uint64_t>::operator~ (void) const {
      return bitVector<uint64_t>(this->width,
				 bitVector<uint64_t>::makeRepresentable(this->width, ~this->value));
    }


    // This is wrong in the signed case as the sign bit it tracks and the sign bit in int64_t are in different places!
    static uint64_t stickyRightShift(const bool value, const bitWidthType width, const uint64_t left, const uint64_t right) {
      //uint64_t bitsInWord = bitVector<uint64_t>::maxWidth();
      uint64_t newValue = left;
      uint64_t stickyBit = 0;
      uint64_t signBit = left & (1ULL << (width - 1));
 
      if (right <= width)  {
	for (uint64_t i = 1; i <= width; i <<= 1) {
	  if (right & i) {
	    uint64_t iOnes = ((1ULL << i) - 1);
	    stickyBit |= ((newValue & iOnes) ? 1 : 0);
	    
	    // Sign extending shift
	    if (signBit) {
	      newValue = (newValue >> i) | (iOnes << (width - i));
	    } else {
	      newValue = (newValue >> i);
	    }
	    
	  }
	}
      } else {
	newValue = (signBit) ? 0xFFFFFFFFFFFFFFFFULL : 0x0;
	stickyBit = (left) ? 0x1 : 0x0;
      }

      return (value) ? newValue : stickyBit;
    }



    template <>
    bitVector<uint64_t> bitVector<uint64_t>::signExtendRightShift (const bitVector<uint64_t> &op) const {
      PRECONDITION(this->width == op.width);
      return bitVector<uint64_t>(this->width,
				 bitVector<uint64_t>::makeRepresentable(this->width, 
									stickyRightShift(true, this->width, this->value, op.value)));
    }

    template <>
    bitVector<int64_t> bitVector<int64_t>::signExtendRightShift (const bitVector<int64_t> &op) const {
      PRECONDITION(this->width == op.width);
      PRECONDITION(this->width < CHAR_BIT*sizeof(int64_t));

      int64_t newValue;
      
      if (this->value < 0) {
	newValue = -(-(this->value) >> op.value) + ((this->value & 0x1) ? -1 : 0); // Rounds away
      } else {
	newValue = this->value >> op.value;
      }
      
      return bitVector<int64_t>(this->width,
				bitVector<int64_t>::makeRepresentable(this->width, newValue));
    }

    template<>
    bitVector<uint64_t> bitVector<uint64_t>::modularLeftShift (const bitVector<uint64_t> &op) const {
      PRECONDITION(this->width == op.width);
      return bitVector<uint64_t>(this->width, 
				 bitVector<uint64_t>::makeRepresentable(this->width,
									(op.value >= this->width) ? 0ULL : this->value << op.value));
    }

    template<>
    bitVector<uint64_t> bitVector<uint64_t>::modularRightShift (const bitVector<uint64_t> &op) const {
      PRECONDITION(this->width == op.width);
      return bitVector<uint64_t>(this->width, 
				 bitVector<uint64_t>::makeRepresentable(this->width,
									(op.value >= this->width) ? 0ULL : this->value >> op.value));
    }


    template <>
    bitVector<uint64_t> bitVector<uint64_t>::modularNegate (void) const {
      return bitVector<uint64_t>(this->width, 
				 bitVector<uint64_t>::makeRepresentable(this->width, ~this->value + 1));
    }

    template <>
    bitVector<int64_t> bitVector<int64_t>::modularNegate (void) const {
      return bitVector<int64_t>(this->width, 
				bitVector<int64_t>::makeRepresentable(this->width, -this->value));
    }


    // Only instantiated for unsigned
    template <>
    bitVector<uint64_t> bitVector<uint64_t>::extract(bitWidthType upper, bitWidthType lower) const {
      PRECONDITION(this->width > upper);
      PRECONDITION(upper >= lower);
      
      bitWidthType newLength = (upper - lower) + 1;
      
      return bitVector<uint64_t>(newLength, 
				 bitVector<uint64_t>::makeRepresentable(newLength, (this->value >> lower)));
    }

    template <>
    bitVector<uint64_t> bitVector<uint64_t>::append(const bitVector<uint64_t> &op) const {
      PRECONDITION(this->width + op.width <= bitVector<uint64_t>::maxWidth());
      
      return bitVector<uint64_t>(this->width + op.width, this->value << op.width | op.value);
    }

    
    
    template <>
    bitVector<typename modifySignedness<uint64_t>::signedVersion> bitVector<uint64_t>::toSigned (void) const {
      return bitVector<int64_t>(this->width, *((int64_t *)&this->value));
    }

    template <>
    bitVector<typename modifySignedness<int64_t>::unsignedVersion> bitVector<int64_t>::toUnsigned (void) const {
      // Note we need to mask out the (sign extensions) of the negative part.
      return bitVector<uint64_t>(this->width, (*((uint64_t *)&this->value)) & bitVector<int64_t>::nOnes(this->width));
    }


    roundingMode traits::RNE (void) { return roundingMode(FE_TONEAREST); }
    roundingMode traits::RNA (void) { return roundingMode(23); }          // Could be better...
    roundingMode traits::RTP (void) { return roundingMode(FE_UPWARD); }
    roundingMode traits::RTN (void) { return roundingMode(FE_DOWNWARD); }
    roundingMode traits::RTZ (void) { return roundingMode(FE_TOWARDZERO); }

  }

  #if 0
  template <>
  simpleExecutable::traits::ubv orderEncode<simpleExecutable::traits, simpleExecutable::traits::ubv> (const simpleExecutable::traits::ubv &b) {
    return orderEncodeBitwise<simpleExecutable::traits, simpleExecutable::traits::ubv>(b);
  }
  #endif
  
}
