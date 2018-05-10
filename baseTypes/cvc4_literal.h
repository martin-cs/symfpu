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
** cvc4_literal.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 12/06/14
**
** There are two back-ends for CVC4, one literal and one symbolic.
** This is due to a quirk of the CVC4 infrastructure.  Nodes can
** contain literals thus anything that is used to implement a literal
** type cannot contain a node or a massive circular definition
** snarl-up happens.  Thus : two back ends.  Not all systems will want
** or need this.
**
** This is the literal back-end which uses CVC4 literal datatypes to
** perform arbitrary precision floating-point computations.  It is to
** implement the CVC4 literal floating-point type.
**
** It is intended to be included *in the middle* of util/floatingpoint.h
**
*/

// Symfpu headers
#include "symfpu/utils/properties.h"
#include "symfpu/utils/numberOfRoundingModes.h"
#include "symfpu/core/ite.h"

// CVC4 headers

// Do not include util/floatingpoint.h as we are in the middle of it!
//#include "util/floatingpoint.h"
#include "util/bitvector.h"


#ifndef SYMFPU_CVC4_LITERAL
#define SYMFPU_CVC4_LITERAL

#include <iostream>

namespace symfpu {
  namespace cvc4_literal {

    typedef unsigned bitWidthType;
    typedef bool proposition;
    typedef ::CVC4::RoundingMode roundingMode;
    typedef ::CVC4::FloatingPointSize floatingPointTypeInfo;

    // Forward declaration
    template <bool T> class bitVector;
    
    // Wrap up the types into one template parameter
    class traits {
    public :
      typedef bitWidthType bwt;
      typedef roundingMode rm;
      typedef floatingPointTypeInfo fpt;
      typedef proposition prop;
      typedef bitVector< true> sbv;
      typedef bitVector<false> ubv;

      static roundingMode RNE (void);
      static roundingMode RNA (void);
      static roundingMode RTP (void);
      static roundingMode RTN (void);
      static roundingMode RTZ (void);

      static void precondition(const prop &p);
      static void postcondition(const prop &p);
      static void invariant(const prop &p);
    };
    



    // Type function
    template <bool T> struct signedToLiteralType;

    template <> struct signedToLiteralType< true> {
      typedef int literalType;
    };
    template <> struct signedToLiteralType<false> {
      typedef  unsigned int literalType;
    };



    template <bool isSigned>
    class bitVector : public ::CVC4::BitVector {
    protected :
      typedef typename signedToLiteralType<isSigned>::literalType literalType;
      typedef ::CVC4::BitVector CVC4BV;

      friend bitVector<!isSigned>;    // To allow conversion between the types
      friend ite<proposition, bitVector<isSigned> >;   // For ITE


    public :
      bitVector (const bitWidthType w, const unsigned v) : CVC4BV(w,v) {}
      bitVector (const proposition &p) : CVC4BV(1,p ? 1U : 0U) {}
      bitVector (const bitVector<isSigned> &old) : CVC4BV(old) {}
      bitVector (const CVC4BV &old) : CVC4BV(old) {}


      bitWidthType getWidth (void) const {
	return getSize();
      }


      /*** Constant creation and test ***/
      
      static bitVector<isSigned> one (const bitWidthType &w);
      static bitVector<isSigned> zero (const bitWidthType &w);
      static bitVector<isSigned> allOnes (const bitWidthType &w);
      
      proposition isAllOnes() const;
      proposition isAllZeros() const;

      static bitVector<isSigned> maxValue (const bitWidthType &w);
      static bitVector<isSigned> minValue (const bitWidthType &w);

      
      /*** Operators ***/
      bitVector<isSigned> operator << (const bitVector<isSigned> &op) const;
      bitVector<isSigned> operator >> (const bitVector<isSigned> &op) const;


      /* Inherited but ...
       * *sigh* if we use the inherited version then it will return a
       * CVC4::BitVector which can be converted back to a
       * bitVector<isSigned> but isn't done automatically when working
       * out types for templates instantiation.  ITE is a particular
       * problem as expressions and constants no longer derive the
       * same type.  There aren't any good solutions in C++, we could
       * use CRTP but Liana wouldn't appreciate that, so the following
       * pointless wrapping functions are needed.
       */

      bitVector<isSigned> operator | (const bitVector<isSigned> &op) const;
      bitVector<isSigned> operator & (const bitVector<isSigned> &op) const;
      bitVector<isSigned> operator + (const bitVector<isSigned> &op) const;
      bitVector<isSigned> operator - (const bitVector<isSigned> &op) const;
      bitVector<isSigned> operator * (const bitVector<isSigned> &op) const;
      bitVector<isSigned> operator / (const bitVector<isSigned> &op) const;
      bitVector<isSigned> operator % (const bitVector<isSigned> &op) const;
      bitVector<isSigned> operator - (void) const;
      bitVector<isSigned> operator ~ (void) const;
      

      bitVector<isSigned> increment () const;
      bitVector<isSigned> decrement () const;
      bitVector<isSigned> signExtendRightShift (const bitVector<isSigned> &op) const;


      /*** Modular opertaions ***/
      // No overflow checking so these are the same as other operations
      bitVector<isSigned> modularLeftShift (const bitVector<isSigned> &op) const;
      bitVector<isSigned> modularRightShift (const bitVector<isSigned> &op) const;
      bitVector<isSigned> modularIncrement () const;
      bitVector<isSigned> modularDecrement () const;
      bitVector<isSigned> modularAdd (const bitVector<isSigned> &op) const;
      bitVector<isSigned> modularNegate () const;



      /*** Comparisons ***/

      /* Inherited ... ish ... */
      proposition operator == (const bitVector<isSigned> &op) const;
      proposition operator <= (const bitVector<isSigned> &op) const;
      proposition operator >= (const bitVector<isSigned> &op) const;
      proposition operator < (const bitVector<isSigned> &op) const;
      proposition operator > (const bitVector<isSigned> &op) const;


      /*** Type conversion ***/
      // CVC4 nodes make no distinction between signed and unsigned, thus ...
      bitVector<true> toSigned (void) const;
      bitVector<false> toUnsigned (void) const;


      /*** Bit hacks ***/

      bitVector<isSigned> extend (bitWidthType extension) const;
      bitVector<isSigned> contract (bitWidthType reduction) const;
      bitVector<isSigned> resize (bitWidthType newSize) const;
      bitVector<isSigned> matchWidth (const bitVector<isSigned> &op) const;
      bitVector<isSigned> append(const bitVector<isSigned> &op) const;

      // Inclusive of end points, thus if the same, extracts just one bit
      bitVector<isSigned> extract(bitWidthType upper, bitWidthType lower) const;

    };
  };

#define CVC4LITITEDFN(T) template <>					\
    struct ite<cvc4_literal::proposition, T> {				\
    static const T & iteOp (const cvc4_literal::proposition &cond,	\
			    const T &l,					\
			    const T &r) {				\
      if (cond) {							\
	return l;							\
      } else {								\
	return r;							\
      }									\
    }									\
  }

  CVC4LITITEDFN(cvc4_literal::traits::rm);
  CVC4LITITEDFN(cvc4_literal::traits::prop);
  CVC4LITITEDFN(cvc4_literal::traits::sbv);
  CVC4LITITEDFN(cvc4_literal::traits::ubv);

#undef CVC4LITITEDFN
  
};

#endif



