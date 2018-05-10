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
** cprover_bvt.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 12/04/17
**
** Implement bit-vectors by generating CPROVER bit-vectors (bvt)
**
*/

// Symfpu headers
#include "../utils/properties.h"
#include "../utils/numberOfRoundingModes.h"
#include "../core/ite.h"

// CPROVER headers
#include <solvers/prop/literal.h>
#include <solvers/flattening/bv_utils.h>

// Common components
#include "cprover_common.h"

#ifndef SYMFPU_CPROVER_BVT
#define SYMFPU_CPROVER_BVT

namespace symfpu {
  namespace cprover_bvt {

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

    // TODO : Fix this hack
    extern bv_utilst *solver;
   
    // (_ BitVec 1) version
    class proposition : public literalt {
    protected :
      friend ite<proposition, proposition>;   // For ITE

    public : 
      proposition (const literalt l) : literalt(l) {}
      proposition (bool v) : literalt() {
	if (v) {
	  make_true();
	} else {
	  make_false();
	}
      }
      proposition (const proposition &old) : literalt(old) {}


      proposition operator ! (void) const {
	return proposition(this->literalt::operator!());
      }

      proposition operator && (const proposition &op) const {
	return proposition(solver->prop.land(*this,op));
      }

      proposition operator || (const proposition &op) const {
	return proposition(solver->prop.lor(*this,op));
      }

      proposition operator == (const proposition &op) const {
	return proposition(solver->prop.lequal(*this,op));
      }

      proposition operator ^ (const proposition &op) const {
	return proposition(solver->prop.lxor(*this,op));
      }

    };



    class roundingMode : public bvt {
    protected:
      friend ite<proposition, roundingMode>;   // For ITE

    public :
      roundingMode (const bvt &op) : bvt(op) {}
      roundingMode (const unsigned v) : bvt(solver->build_constant(v,32)) {}
      roundingMode (const roundingMode &old) : bvt(old) {}

      proposition valid (void) const {
	// TODO : Improve...
	return true;
      }

      proposition operator == (const roundingMode &op) const {
	return proposition(solver->equal(*this,op));
      }

    };



    template <bool isSigned>
    class bitVector : public bvt {
    protected :
      friend bitVector<!isSigned>;    // To allow conversion between the types
      friend ite<proposition, bitVector<isSigned> >;   // For ITE

    public :
      bitVector (const bitWidthType w, const unsigned v) : bvt(solver->build_constant(v,w)) {}
      bitVector (const bvt &old) : bvt(old) {}
      bitVector (const proposition &p) { push_back(p); }
      bitVector (const bitVector<isSigned> &old) : bvt(old) {}

      bitWidthType getWidth (void) const {
	return size();
      }


      /*** Constant creation and test ***/
      
      static bitVector<isSigned> one (const bitWidthType &w) { return bitVector<isSigned>(w,1); }
      static bitVector<isSigned> zero (const bitWidthType &w)  { return bitVector<isSigned>(w,0); }
      static bitVector<isSigned> allOnes (const bitWidthType &w) { return ~zero(w); }

      
      inline proposition isAllOnes() const { return proposition(solver->is_all_ones(*this)); }
      inline proposition isAllZeros() const { return proposition(solver->is_zero(*this)); }

      static bitVector<isSigned> maxValue (const bitWidthType &w) {
	if (isSigned) {
	  return zero(1).append(allOnes(w-1));
	} else {
	  return allOnes(w);
	}
      }
      
      static bitVector<isSigned> minValue (const bitWidthType &w) {
	if (isSigned) {
	  return allOnes(w);
	} else {
	  return zero(w);
	}
      }

      
      /*** Operators ***/
      inline bitVector<isSigned> operator << (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(solver->shift(*this, bv_utilst::shiftt::LEFT, op));
      }

      inline bitVector<isSigned> operator >> (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(solver->shift(*this, isSigned ? bv_utilst::shiftt::ARIGHT : bv_utilst::shiftt::LRIGHT, op));
      }


      inline bitVector<isSigned> operator | (const bitVector<isSigned> &op) const {
	bitVector<isSigned> bv(*this);

	unsigned width = bv.size();
	for(std::size_t i=0; i<width; i++)
	  bv[i]=solver->prop.lor(bv[i], op[i]);

	return bv;
      }

      inline bitVector<isSigned> operator & (const bitVector<isSigned> &op) const {
	bitVector<isSigned> bv(*this);

	unsigned width = bv.size();
	for(std::size_t i=0; i<width; i++)
	  bv[i]=solver->prop.land(bv[i], op[i]);

	return bv;
      }

      inline bitVector<isSigned> operator + (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(solver->add(*this, op));
      }

      inline bitVector<isSigned> operator - (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(solver->sub(*this, op));
      }

      inline bitVector<isSigned> operator * (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(isSigned ?
				   solver->signed_multiplier(*this, op) :
				   solver->unsigned_multiplier(*this, op));
      }
      
      inline bitVector<isSigned> operator / (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(solver->divider(*this, op, isSigned ? bv_utilst::representationt::SIGNED : bv_utilst::representationt::UNSIGNED));
      }
      
      inline bitVector<isSigned> operator % (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(solver->remainder(*this, op, isSigned ? bv_utilst::representationt::SIGNED : bv_utilst::representationt::UNSIGNED));
      }

      
      inline bitVector<isSigned> operator - (void) const {
	return bitVector<isSigned>(solver->negate(*this));
      }

      inline bitVector<isSigned> operator ~ (void) const {
	bitVector<isSigned> bv(*this);

	unsigned width = bv.size();
	for(std::size_t i=0; i<width; i++)
	  bv[i]=!bv[i];

	return bv;
      }

      inline bitVector<isSigned> increment () const {
	return bitVector<isSigned>(solver->inc(*this));
      }

      inline bitVector<isSigned> decrement () const {
	return *this - bitVector<isSigned>::one(getWidth());
      }

      inline bitVector<isSigned> signExtendRightShift (const bitVector<isSigned> &op) const {
	return bitVector<isSigned>(solver->shift(*this, bv_utilst::shiftt::ARIGHT, op));
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
	return proposition(solver->equal(*this, op));
      }

      inline proposition operator <= (const bitVector<isSigned> &op) const {
	return proposition(solver->rel(*this, ID_le, op, isSigned ? bv_utilst::representationt::SIGNED : bv_utilst::representationt::UNSIGNED));
      }

      inline proposition operator >= (const bitVector<isSigned> &op) const {
	return proposition(solver->rel(*this, ID_ge, op, isSigned ? bv_utilst::representationt::SIGNED : bv_utilst::representationt::UNSIGNED));
      }

      inline proposition operator < (const bitVector<isSigned> &op) const {
	return proposition(solver->rel(*this, ID_lt, op, isSigned ? bv_utilst::representationt::SIGNED : bv_utilst::representationt::UNSIGNED));
      }

      inline proposition operator > (const bitVector<isSigned> &op) const {
	return proposition(solver->rel(*this, ID_gt, op, isSigned ? bv_utilst::representationt::SIGNED : bv_utilst::representationt::UNSIGNED));
      }

      /*** Type conversion ***/
      // CPROVER bvts make no distinction between signed and unsigned, thus ...
      bitVector<true> toSigned (void) const {
	return bitVector<true>(*this);
      }
      bitVector<false> toUnsigned (void) const {
	return bitVector<false>(*this);
      }



      /*** Bit hacks ***/

      inline bitVector<isSigned> extend (bitWidthType extension) const {
	return bitVector<isSigned>(solver->extension(*this, this->size() + extension, isSigned ? bv_utilst::representationt::SIGNED : bv_utilst::representationt::UNSIGNED));	
      }

      inline bitVector<isSigned> contract (bitWidthType reduction) const {
	bitVector<isSigned> bv(*this);
	bv.bvt::resize(bv.size() - reduction);
	return bv;
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
	bitVector<isSigned> bv(op);

	bv.insert(bv.end(), this->begin(), this->end());

	return bv;
      }

      // Inclusive of end points, thus if the same, extracts just one bit
      bitVector<isSigned> extract(bitWidthType upper, bitWidthType lower) const {
	PRECONDITION(upper >= lower);

	bitVector<isSigned> bv(upper-lower + 1, 0);

	for (bitWidthType i = lower; i <= upper; ++i) {
	  bv[i - lower] = (*this)[i];
	}
	
	return bv;
      }
    };

  };

  
  template <>
  struct ite<cprover_bvt::proposition, cprover_bvt::proposition> {
    static const cprover_bvt::proposition iteOp (const cprover_bvt::proposition &cond,
						 const cprover_bvt::proposition &l,
						 const cprover_bvt::proposition &r) {
      return cprover_bvt::proposition(symfpu::cprover_bvt::solver->prop.lselect(cond,l,r));
    }
  };

#define CPROVERBVTITEDFN(T) template <>					\
    struct ite<cprover_bvt::proposition, T> {				\
    static const T iteOp (const cprover_bvt::proposition &cond,	\
			  const T &l,					\
			  const T &r) {					\
      assert(l.size() == r.size());					\
      return T(symfpu::cprover_bvt::solver->select(cond,l,r));		\
    }									\
  };

  CPROVERBVTITEDFN(cprover_bvt::traits::rm);
  CPROVERBVTITEDFN(cprover_bvt::traits::sbv);
  CPROVERBVTITEDFN(cprover_bvt::traits::ubv);

#undef CPROVERBVTITEDFN
  
};

#endif
