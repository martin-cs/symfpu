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
** cprover_exprt.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 27/06/17
**
** Implement bit-vectors by generating CPROVER bit-vectors (bvt)
**
*/

// Symfpu headers
#include "../utils/properties.h"
#include "../utils/numberOfRoundingModes.h"
#include "../core/ite.h"
#include "../core/operations.h"

// CPROVER headers
#include <util/arith_tools.h>
#include <util/std_expr.h>
#include <util/std_types.h>

// Common components
#include "cprover_common.h"

#ifndef SYMFPU_CPROVER_EXPRT
#define SYMFPU_CPROVER_EXPRT

namespace symfpu {
  namespace cprover_exprt {

    typedef symfpu::cprover_common::bitWidthType bitWidthType;
    typedef symfpu::cprover_common::floatingPointTypeInfo floatingPointTypeInfo;

    // Forward declarations
    class roundingMode;
    class proposition;
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

      // Literal invariants
      inline static void precondition (const bool b) { assert(b); return; }
      inline static void postcondition (const bool b) { assert(b); return; }
      inline static void invariant (const bool b) { assert(b); return; }

      // Symbolic invariants
      // TODO : Add to the solver
      inline static void precondition(const prop &p) { return; }
      inline static void postcondition(const prop &p) { return; }
      inline static void invariant(const prop &p) { return; }
     
    };

    // To simplify the property macros
    typedef traits t;

   
    // exprt version using bools
    class proposition : public exprt {
    protected :
      friend ite<proposition, proposition>;   // For ITE

      exprt simplify (const exprt &e) const;
      exprt tag (const exprt &e) const;
      
    public : 
      proposition (const exprt e) : exprt(tag(simplify(e))) {}
      proposition (bool v) : exprt(tag(v ?
				       static_cast<constant_exprt>(true_exprt()) :
				       static_cast<constant_exprt>(false_exprt()) )) {}
      proposition (const proposition &old) : exprt(old) {}


      proposition operator ! (void) const {
	return proposition(not_exprt(*this));
      }

      proposition operator && (const proposition &op) const {
	return proposition(and_exprt(*this, op));
      }

      proposition operator || (const proposition &op) const {
	return proposition(or_exprt(*this, op));
      }

      proposition operator == (const proposition &op) const {
	return proposition(equal_exprt(*this, op));
      }

      proposition operator ^ (const proposition &op) const {
	return proposition(equal_exprt(*this, not_exprt(op)));
      }
    };



    class roundingMode : public exprt {
    protected:
      friend ite<proposition, roundingMode>;   // For ITE

    public :
      roundingMode (const exprt &op) : exprt(op) {}
      roundingMode (const unsigned v) : exprt(from_integer(v, signedbv_typet(32))) {}
      roundingMode (const roundingMode &old) : exprt(old) {}

      proposition valid (void) const {
	// TODO : Improve...
	return true;
      }

      proposition operator == (const roundingMode &op) const {
	return proposition(equal_exprt(*this, op));
      }

    };



    template <bool isSigned>
    class bitVector : public exprt {
    protected :
      friend bitVector<!isSigned>;    // To allow conversion between the types
      friend ite<proposition, bitVector<isSigned> >;   // For ITE

      exprt simplify (const exprt &e) const;
      exprt tag (const exprt &e) const;
      
    public :
      bitVector (const bitWidthType w, const unsigned v) :
      exprt(tag((isSigned) ? from_integer(v, signedbv_typet(w)) : from_integer(v, unsignedbv_typet(w))))
	{}
      bitVector (const exprt &e) : exprt(tag(simplify(e))) {}
      bitVector (const proposition &p) : exprt(tag((isSigned) ? typecast_exprt(p, signedbv_typet(1)) : typecast_exprt(p, unsignedbv_typet(1)))) {}
      bitVector (const bitVector<isSigned> &old) : exprt(old) {}

      bitWidthType getWidth (void) const {
	if (isSigned) {
	  return to_signedbv_type(this->type()).get_width();
	} else {
	  return to_unsignedbv_type(this->type()).get_width();
	}
      }


      /*** Constant creation and test ***/
      
      static bitVector<isSigned> one (const bitWidthType &w) { return bitVector<isSigned>(w,1); }
      static bitVector<isSigned> zero (const bitWidthType &w)  { return bitVector<isSigned>(w,0); }
      static bitVector<isSigned> allOnes (const bitWidthType &w) { return ~zero(w); }

      
      inline proposition isAllOnes() const { return (*this) == allOnes(this->getWidth()); }
      inline proposition isAllZeros() const { return (*this) == zero(this->getWidth()); }

      static bitVector<isSigned> maxValue (const bitWidthType &w) {
	if (isSigned) {
	  return signedbv_typet(w).largest_expr();
	} else {
	  return unsignedbv_typet(w).largest_expr();
	}
      }
      
      static bitVector<isSigned> minValue (const bitWidthType &w) {
	if (isSigned) {
	  return signedbv_typet(w).smallest_expr();
	} else {
	  return unsignedbv_typet(w).smallest_expr();
	}
      }

      
      /*** Operators ***/
      inline bitVector<isSigned> operator << (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(shl_exprt(*this, op));
      }

      inline bitVector<isSigned> operator >> (const bitVector<isSigned> &op) const {
	if (isSigned) {
	  return bitVector<isSigned>(ashr_exprt(*this, op));
	} else {
	  return bitVector<isSigned>(lshr_exprt(*this, op));
	}
      }


      inline bitVector<isSigned> operator | (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(bitor_exprt(*this, op));
      }

      inline bitVector<isSigned> operator & (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(bitand_exprt(*this, op));
      }

      inline bitVector<isSigned> operator + (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(plus_exprt(*this, op));
      }

      inline bitVector<isSigned> operator - (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(minus_exprt(*this, op));
      }

      inline bitVector<isSigned> operator * (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(mult_exprt(*this, op));
      }
      
      inline bitVector<isSigned> operator / (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(div_exprt(*this, op));
      }
      
      inline bitVector<isSigned> operator % (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(mod_exprt(*this, op));
      }

      
      inline bitVector<isSigned> operator - (void) const {
	return bitVector<isSigned>(unary_minus_exprt(*this));
      }

      inline bitVector<isSigned> operator ~ (void) const {
	return bitVector<isSigned>(bitnot_exprt(*this));
      }

      inline bitVector<isSigned> increment () const {
	return bitVector<isSigned>(*this + bitVector<isSigned>::one(getWidth()));
      }

      inline bitVector<isSigned> decrement () const {
	return bitVector<isSigned>(*this - bitVector<isSigned>::one(getWidth()));
      }

      inline bitVector<isSigned> signExtendRightShift (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(ashr_exprt(*this, op));
      }



      /*** Modular opertaions ***/
      // No overflow checking so these are the same as other operations
      inline bitVector<isSigned> modularLeftShift (const bitVector<isSigned> &op) const {
	return *this << op;
      }
      
      inline bitVector<isSigned> modularRightShift (const bitVector<isSigned> &op) const {
	return *this >> op;
      }

      inline bitVector<isSigned> modularIncrement () const {
	return this->increment();
      }

      inline bitVector<isSigned> modularDecrement () const {
	return this->decrement();
      }

      inline bitVector<isSigned> modularAdd (const bitVector<isSigned> &op) const {
	return *this + op;
      }

      inline bitVector<isSigned> modularNegate () const {
	return -(*this);
      }




      /*** Comparisons ***/

      inline proposition operator == (const bitVector<isSigned> &op) const {
	return proposition(equal_exprt(*this, op));
      }

      inline proposition operator <= (const bitVector<isSigned> &op) const {
	return proposition(binary_relation_exprt(*this, ID_le, op));
      }

      inline proposition operator >= (const bitVector<isSigned> &op) const {
	return proposition(binary_relation_exprt(*this, ID_ge, op));
      }

      inline proposition operator < (const bitVector<isSigned> &op) const {
	return proposition(binary_relation_exprt(*this, ID_lt, op));
      }

      inline proposition operator > (const bitVector<isSigned> &op) const {
	return proposition(binary_relation_exprt(*this, ID_gt, op));
      }

      /*** Type conversion ***/
      // CPROVER bvts make no distinction between signed and unsigned, thus ...
      bitVector<true> toSigned (void) const {
	return bitVector<true>(typecast_exprt(*this, signedbv_typet(getWidth())));
      }
      bitVector<false> toUnsigned (void) const {
	return bitVector<false>(typecast_exprt(*this, unsignedbv_typet(getWidth())));
      }



      /*** Bit hacks ***/

      inline bitVector<isSigned> extend (bitWidthType extension) const {
	if (isSigned) {
	  return bitVector<isSigned>(typecast_exprt(*this, signedbv_typet(getWidth() + extension)));
	} else {
	  return bitVector<isSigned>(typecast_exprt(*this, unsignedbv_typet(getWidth() + extension)));
	}
      }

      inline bitVector<isSigned> contract (bitWidthType reduction) const {
	if (isSigned) {
	  return bitVector<isSigned>(typecast_exprt(*this, signedbv_typet(getWidth() - reduction)));
	} else {
	  return bitVector<isSigned>(typecast_exprt(*this, unsignedbv_typet(getWidth() - reduction)));
	}
      }

      inline bitVector<isSigned> resize (bitWidthType newSize) const {
	bitWidthType width = this->getWidth();

	if (newSize > width) {
	  return this->extend(newSize - width);
	} else if (newSize < width) {
	  return this->contract(width - newSize);
	} else {
	  return *this;
	}
      }

      inline bitVector<isSigned> matchWidth (const bitVector<isSigned> &op) const {
	PRECONDITION(this->getWidth() <= op.getWidth());
	return this->extend(op.getWidth() - this->getWidth());
      }


      bitVector<isSigned> append(const bitVector<isSigned> &op) const {
	if (isSigned) {
	  return bitVector<isSigned>(concatenation_exprt(*this, op, signedbv_typet(this->getWidth() + op.getWidth())));
	} else {
	  return bitVector<isSigned>(concatenation_exprt(*this, op, unsignedbv_typet(this->getWidth() + op.getWidth())));
	}
      }

      // Inclusive of end points, thus if the same, extracts just one bit
      bitVector<isSigned> extract(bitWidthType upper, bitWidthType lower) const {
	PRECONDITION(upper >= lower);
	exprt u(bitVector<isSigned>(getWidth(),upper));
	exprt l(bitVector<isSigned>(getWidth(),lower));
	bitWidthType extractedWidth = upper - lower + 1;

	if (isSigned) {
	  return bitVector<isSigned>(extractbits_exprt(*this, u, l, signedbv_typet(extractedWidth)));
	} else {
	  return bitVector<isSigned>(extractbits_exprt(*this, u, l, unsignedbv_typet(extractedWidth)));
	}
      }
    };

  };

  
  template <>
  struct ite<cprover_exprt::proposition, cprover_exprt::proposition> { 
    static const cprover_exprt::proposition iteOp (const cprover_exprt::proposition &cond,
						 const cprover_exprt::proposition &l,
						 const cprover_exprt::proposition &r) {
      return cprover_exprt::proposition(if_exprt(cond, l, r));
    }
  };

#define CPROVERBVTITEDFN(T) template <>					\
    struct ite<cprover_exprt::proposition, T> {				\
    static const T iteOp (const cprover_exprt::proposition &cond,	\
			  const T &l,					\
			  const T &r) {					\
      assert(l.type() == r.type());					\
      return T(if_exprt(cond,l,r));					\
    }									\
  };

  CPROVERBVTITEDFN(cprover_exprt::traits::rm);
  CPROVERBVTITEDFN(cprover_exprt::traits::sbv);
  CPROVERBVTITEDFN(cprover_exprt::traits::ubv);

#undef CPROVERBVTITEDFN

  // Use improved version of some of the operations
  #ifndef SYMFPU_WORDLEVEL_ENCODINGS
  template <>
  cprover_exprt::traits::ubv orderEncode<cprover_exprt::traits, cprover_exprt::traits::ubv> (const cprover_exprt::traits::ubv &b);

  template <>
  stickyRightShiftResult<cprover_exprt::traits> stickyRightShift (const cprover_exprt::traits::ubv &input, const cprover_exprt::traits::ubv &shiftAmount);
  #endif
  
  template <>
  void probabilityAnnotation<cprover_exprt::traits, cprover_exprt::traits::prop> (const cprover_exprt::traits::prop &p, const probability &pr);

  
};

#endif
