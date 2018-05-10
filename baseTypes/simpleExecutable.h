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
** simpleExecutable.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 03/06/14
**
** The most simple executable implementation of bit-vectors.
** Limited in the ranges it supports but fast and suitable for reasoning.
**
** Unless otherwise stated, bit operations require all operands to have the
** same widths and the result will have the same width (SMT-LIB style,
** but see the 'expanding' instructions below).  Also underflow and
** overflow are considered errors (unlike SMT-LIB but see the 'modular' instructions
** below) although error checking is not perfect as overflows in the
** underlying data type can mask errors.
**
*/

#include "symfpu/utils/properties.h"
#include "symfpu/core/ite.h"
#include "symfpu/baseTypes/shared.h"

#include <assert.h>
#include <stdint.h>
#include <limits.h>

#ifndef SYMFPU_SIMPLE_EXECUTABLE
#define SYMFPU_SIMPLE_EXECUTABLE

namespace symfpu {
  namespace simpleExecutable {

    typedef symfpu::shared::bitWidthType bitWidthType;
    typedef symfpu::shared::executable_proposition proposition;
    typedef symfpu::shared::floatingPointTypeInfo floatingPointTypeInfo;
    
    // Forwards definitions
    class roundingMode;
    template <class T> class bitVector;

    
    // This is the class that is used as a template argument
    class traits {
    public :
      typedef bitWidthType bwt;
      typedef roundingMode rm;
      typedef floatingPointTypeInfo fpt;
      typedef proposition prop;
      typedef bitVector< int64_t> sbv;
      typedef bitVector<uint64_t> ubv;

      static roundingMode RNE(void);
      static roundingMode RNA(void);
      static roundingMode RTP(void);
      static roundingMode RTN(void);
      static roundingMode RTZ(void);

      // As prop == bool only one set of these is needed
      inline static void precondition(const bool b) { assert(b); return; }
      inline static void postcondition(const bool b) { assert(b); return; }
      inline static void invariant(const bool b) { assert(b); return; }

    };

    // To simplify the property macros
    typedef traits t;

    
    
    class roundingMode {
    private :
      int value;
      
    public :
      roundingMode (int v) : value(v) {}
      roundingMode (const roundingMode &old) : value(old.value) {}
      
      roundingMode & operator = (const roundingMode &op) {
	this->value = op.value;
	return (*this);
      }
      
      proposition operator == (const roundingMode &op) const {
	return proposition(this->value == op.value);
      }

      // Only for executable back-ends
      int getValue (void) const {
	return this->value;
      }
    };
    


    template <typename T> struct modifySignedness;

    template <> struct modifySignedness< int64_t> { 
      typedef uint64_t unsignedVersion; 
      typedef  int64_t   signedVersion; 
    };
    template <> struct modifySignedness<uint64_t> { 
      typedef uint64_t unsignedVersion; 
      typedef  int64_t   signedVersion; 
    };



    template <class T>
    class bitVector {
    protected :
      bitWidthType width;
      T value;

      static bitWidthType maxWidth (void) {
	return sizeof(T)*CHAR_BIT;
      }

      static T nOnes (bitWidthType n) {
	if (n == 0) {
	  return 0;
	} else {
	  // Not (1 << n) - 1 to avoid overflow for n = maxWidth()
	  bitWidthType shift = bitVector<T>::maxWidth() - n;
	  return ((~0ULL) << shift) >> shift;
	}
      }

      // Bit vectors should store values in 2's complement using the full width
      // (and thus may have significant bits outside the width).
      // Thus need to make sure the value stored is representable within
      // bit-vectors of the specified width.
      static bool isRepresentable (const bitWidthType w, const T v);

      // Modular operations need to reduce operations back to
      // something representable.
      static T makeRepresentable (const bitWidthType w, const T v);
      

    public :

      // Ideally should protect but is used in subnormal rounding
      bitVector (const bitWidthType w, const T v) : width(w), value(v)
      {
	PRECONDITION(width <= bitVector<T>::maxWidth()); 
	PRECONDITION(0 < width); 
	PRECONDITION(bitVector<T>::isRepresentable(w,v));
      }

      bitVector (const proposition &p) : width(1), value(p) {}
      
      bitVector (const bitVector<T> &old) : width(old.width), value(old.value) {}
      
      // Constructors
      //   non-det   on other instances but not this so that
      //             instantiation catches this




      bitWidthType getWidth (void) const {
	return this->width;
      }
      
      // Would it be better to not have this and only have copy?
      bitVector<T> & operator= (const bitVector<T> &op) {
	PRECONDITION(op.width == this->width);
	
	this->value = op.value;
	
	return (*this);
      }


      /*** Constant creation and test ***/

      static bitVector<T> one (const bitWidthType &w) { return bitVector<T>(w,1); }
      static bitVector<T> zero (const bitWidthType &w)  { return bitVector<T>(w,0); }
      static bitVector<T> allOnes (const bitWidthType &w)  { return bitVector<T>(w,bitVector<T>::nOnes(w)); }
      
      inline proposition isAllOnes() const {return proposition(((~this->value) & nOnes(this->width)) == 0);}
      inline proposition isAllZeros() const {return proposition(this->value == 0);}

      static bitVector<T> maxValue (const bitWidthType &w);
      static bitVector<T> minValue (const bitWidthType &w);

      /*** Operators ***/
      // Need to inline the operations where possible
      inline bitVector<T> operator << (const bitVector<T> &op) const {
	PRECONDITION(this->width == op.width);
	PRECONDITION(op.value >= 0U && op.value < (T)this->width);
	return bitVector<T>(this->width, this->value << op.value);
      }

      inline bitVector<T> operator >> (const bitVector<T> &op) const {
	PRECONDITION(this->width == op.width);
	PRECONDITION(op.value >= 0U && op.value < (T)this->width);
	return bitVector<T>(this->width, this->value >> op.value);
      }

      inline bitVector<T> operator | (const bitVector<T> &op) const {
	PRECONDITION(this->width == op.width);
	return bitVector<T>(this->width, this->value | op.value);
      }

      inline bitVector<T> operator & (const bitVector<T> &op) const {
	PRECONDITION(this->width == op.width);
	return bitVector<T>(this->width, this->value & op.value);
      }

      inline bitVector<T> operator + (const bitVector<T> &op) const {
	PRECONDITION(this->width == op.width);
	return bitVector<T>(this->width, this->value + op.value);
      }

      inline bitVector<T> operator - (const bitVector<T> &op) const {
	PRECONDITION(this->width == op.width);
	return bitVector<T>(this->width, this->value - op.value);
      }

      inline bitVector<T> operator * (const bitVector<T> &op) const {
	PRECONDITION(this->width == op.width);
	return bitVector<T>(this->width, this->value * op.value);
      }

      inline bitVector<T> operator / (const bitVector<T> &op) const {
	PRECONDITION(op.value != 0);
	PRECONDITION(this->width == op.width);
	return bitVector<T>(this->width, this->value / op.value);
      }

      inline bitVector<T> operator % (const bitVector<T> &op) const {
	PRECONDITION(op.value != 0);
	PRECONDITION(this->width == op.width);
	return bitVector<T>(this->width, this->value % op.value);
      }

      bitVector<T> operator - (void) const;
      bitVector<T> operator ~ (void) const;

      inline bitVector<T> increment () const {
	return bitVector<T>(this->width, this->value + 1);
      }

      inline bitVector<T> decrement () const {
	return bitVector<T>(this->width, this->value - 1);
      }

      bitVector<T> signExtendRightShift (const bitVector<T> &op) const;


      /*** Modular operations ***/
      bitVector<T> modularLeftShift (const bitVector<T> &op) const;
      bitVector<T> modularRightShift (const bitVector<T> &op) const;

      inline bitVector<T> modularIncrement () const {
	return bitVector<T>(this->width, 
			    bitVector<T>::makeRepresentable(this->width, this->value + 1));
      }
      
      inline bitVector<T> modularDecrement () const {
	return bitVector<T>(this->width, 
			    bitVector<T>::makeRepresentable(this->width, this->value - 1));
      }

      inline bitVector<T> modularAdd (const bitVector<T> &op) const {
	return bitVector<T>(this->width, 
			    bitVector<T>::makeRepresentable(this->width, 
							    this->value + op.value));
      }

      bitVector<T> modularNegate () const;




      /*** Comparisons ***/

      inline proposition operator == (const bitVector<T> &op) const {
	PRECONDITION(this->width == op.width);
	return proposition(this->value == op.value);
      }

      inline proposition operator <= (const bitVector<T> &op) const {
	PRECONDITION(this->width == op.width);
	return proposition(this->value <= op.value);
      }

      inline proposition operator >= (const bitVector<T> &op) const {
	PRECONDITION(this->width == op.width);
	return proposition(this->value >= op.value);
      }

      inline proposition operator < (const bitVector<T> &op) const {
	PRECONDITION(this->width == op.width);
	return proposition(this->value < op.value);
      }

      inline proposition operator > (const bitVector<T> &op) const {
	PRECONDITION(this->width == op.width);
	return proposition(this->value > op.value);
      }


      /*** Type conversion ***/

      bitVector<typename modifySignedness<T>::signedVersion> toSigned (void) const;
      bitVector<typename modifySignedness<T>::unsignedVersion> toUnsigned (void) const;



      /*** Bit hacks ***/

      inline bitVector<T> extend (bitWidthType extension) const {
	PRECONDITION(this->width + extension <= bitVector<T>::maxWidth());

	// No extension needed, even in the signed case as already correctly represented
	return bitVector<T>(this->width + extension, this->value);
      }

      inline bitVector<T> contract (bitWidthType reduction) const {
	PRECONDITION(this->width > reduction);

	return bitVector<T>(this->width - reduction, this->value);
      }

      inline bitVector<T> resize (bitWidthType newSize) const {
	return bitVector<T>(newSize, 
			    bitVector<T>::makeRepresentable(newSize, this->value));
      }

      inline bitVector<T> matchWidth (const bitVector<T> &op) const {
	PRECONDITION(this->width <= op.width);
	return this->extend(op.width - this->width);
      }



      bitVector<T> append(const bitVector<T> &op) const;

      // Inclusive of end points, thus if the same, extracts just one bit
      bitVector<T> extract(bitWidthType upper, bitWidthType lower) const;

      // Only meaningful for executable implementations
      T contents (void) const { return this->value; }

    };

  }


#define SEITEDFN(T) template <>						\
    struct ite<simpleExecutable::traits::prop, T> {			\
    static const T & iteOp (const simpleExecutable::traits::prop &cond,	\
			    const T &l,					\
			    const T &r) {				\
      if (cond) {							\
	return l;							\
      } else {								\
	return r;							\
      }									\
    }									\
  }

  SEITEDFN(simpleExecutable::traits::rm);
  SEITEDFN(simpleExecutable::traits::prop);

#undef SEITEDFN
  
#define SEITEDFNW(T) template <>					\
    struct ite<simpleExecutable::traits::prop, T> {			\
    static const T & iteOp (const simpleExecutable::traits::prop &cond,	\
			    const T &l,					\
			    const T &r) {				\
      assert(l.getWidth() == r.getWidth());				\
									\
      if (cond) {							\
	return l;							\
      } else {								\
	return r;							\
      }									\
    }									\
  }

  SEITEDFNW(simpleExecutable::traits::sbv);
  SEITEDFNW(simpleExecutable::traits::ubv);

#undef SEITEDFNW

}


// For testing only; bitwise implementation is way slower for software
#if 0
#include "../core/operations.h"

namespace symfpu {

  template <>
    simpleExecutable::traits::ubv orderEncode<simpleExecutable::traits, simpleExecutable::traits::ubv> (const simpleExecutable::traits::ubv &b);

}
#endif


#endif
