/* Part of SymFPU, see LICENSE for licensing information */
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
      if (w == bitVector<int64_t>::maxWidth()) { return true; }
      uint64_t shiftSafe = *((uint64_t *)(&v));
      uint64_t top = (shiftSafe >> w);
      uint64_t signbit = shiftSafe & 0x8000000000000000;
      int64_t stop = *((int64_t *)(&top));
      return (signbit) ? (stop == bitVector<int64_t>::nOnes(bitVector<int64_t>::maxWidth() - w)) : (stop == 0LL);
    }

    template <>
    bool bitVector<uint64_t>::isRepresentable (const bitWidthType w, const uint64_t v) {
      if (w == bitVector<uint64_t>::maxWidth()) { return true; }
      uint64_t top = (v >> w);
      return (top == 0);
    }

    
    template <>
    uint64_t bitVector<uint64_t>::makeRepresentable (const bitWidthType w, const uint64_t v) {
      return v & bitVector<uint64_t>::nOnes(w);
    }

    template <>
    int64_t bitVector<int64_t>::makeRepresentable (const bitWidthType w, const int64_t v) {
      // Mask to w bits then sign-extend.  Avoids overflow at w == maxWidth()
      // and matches the two's-complement wrap that modular operations expect.
      uint64_t mask(bitVector<int64_t>::nOnes(w));
      uint64_t bits(*((uint64_t *)(&v)) & mask);
      uint64_t signBit(1ULL << (w - 1));
      uint64_t extended((bits & signBit) ? (bits | ~mask) : bits);
      return *((int64_t *)(&extended));
    }

    
    template <>
    bitVector<int64_t> bitVector<int64_t>::maxValue (const bitWidthType &w) {
      PRECONDITION(w != 1);
      return bitVector<int64_t>(w, (1ULL << (w - 1)) - 1);
    }
    
    template <>
    bitVector<uint64_t> bitVector<uint64_t>::maxValue (const bitWidthType &w) {
      PRECONDITION(w != 1);
      return bitVector<uint64_t>(w, bitVector<uint64_t>::nOnes(w));
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
      // Negate in unsigned arithmetic to avoid -INT64_MIN signed overflow UB.
      uint64_t bits(~(*((uint64_t *)(&this->value))) + 1);
      return bitVector<int64_t>(this->width,
				bitVector<int64_t>::makeRepresentable(this->width, *((int64_t *)(&bits))));
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

      // Reuse the unsigned helper, treating the stored bit pattern as the
      // width-bit two's-complement value to be arithmetically shifted.
      uint64_t bits(*((uint64_t *)(&this->value)) & bitVector<int64_t>::nOnes(this->width));
      uint64_t shifted(stickyRightShift(true, this->width, bits, static_cast<uint64_t>(op.value)));
      return bitVector<int64_t>(this->width, *((int64_t *)(&shifted)));
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
      // Negate in unsigned arithmetic to avoid -INT64_MIN signed overflow UB.
      uint64_t bits(~(*((uint64_t *)(&this->value))) + 1);
      return bitVector<int64_t>(this->width,
				bitVector<int64_t>::makeRepresentable(this->width, *((int64_t *)(&bits))));
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
