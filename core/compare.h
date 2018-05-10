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
** compare.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 25/08/14
**
** Comparison between floating-point numbers
**
*/

#include "symfpu/core/unpackedFloat.h"
#include "symfpu/core/ite.h"

#ifndef SYMFPU_COMPARE
#define SYMFPU_COMPARE

namespace symfpu {

  // SMT-LIB equality
  template <class t>
    typename t::prop smtlibEqual (const typename t::fpt &format, 
				  const unpackedFloat<t> &left,
				  const unpackedFloat<t> &right) {
    typedef typename t::prop prop;
    
    PRECONDITION(left.valid(format));
    PRECONDITION(right.valid(format));

    // Relies on a number of properties of the unpacked format
    // particularly the use of default exponents, significands and signs

    prop flagsEqual((left.getNaN() == right.getNaN()) &&
		    (left.getInf() == right.getInf()) &&
		    (left.getZero() == right.getZero()) &&
		    (left.getSign() == right.getSign()));
    
    prop flagsAndExponent(flagsEqual && left.getExponent() == right.getExponent());

    // Avoid comparing (and thus instantiating) the significand unless necessary
    probabilityAnnotation<t,prop>(flagsAndExponent, UNLIKELY);

    prop res(ITE(flagsAndExponent,
		 left.getSignificand() == right.getSignificand(),
		 prop(false)));

    return res;
  }
  
  // IEEE-754 Equality (not actually an equivalence relation but ...)
  template <class t>
    typename t::prop ieee754Equal (const typename t::fpt &format, 
				   const unpackedFloat<t> &left,
				   const unpackedFloat<t> &right) {

    typedef typename t::prop prop;

    PRECONDITION(left.valid(format));
    PRECONDITION(right.valid(format));

    prop neitherNan(!left.getNaN() && !right.getNaN());   // All comparison with NaN are false
    prop bothZero(left.getZero() && right.getZero());  // Both zeros are equal
    prop neitherZero(!left.getZero() && !right.getZero());

    prop flagsAndExponent(neitherNan &&
			  (bothZero || (neitherZero &&
					(left.getInf() == right.getInf() && 
					 left.getSign() == right.getSign() &&
					 left.getExponent() == right.getExponent()))));

    // Avoid comparing (and thus instantiating) the significand unless necessary
    probabilityAnnotation<t,prop>(flagsAndExponent, UNLIKELY);
    
    prop res(ITE(flagsAndExponent, left.getSignificand() == right.getSignificand(), prop(false)));

    return res;
  }
  

  // Share the common comparison code between functions
  // equality == true if the equal case returns true
  // IEEE-754 semantics for ordering with NaN
  // (i.e. unordered with everything, not even equal to itself)
  template <class t>
    typename t::prop ordering (const typename t::fpt &format, 
			       const unpackedFloat<t> &left,
			       const unpackedFloat<t> &right,
			       const typename t::prop equality) {

    typedef typename t::prop prop;

    PRECONDITION(left.valid(format));
    PRECONDITION(right.valid(format));

    // All comparison with NaN are false
    prop neitherNaN(!left.getNaN() && !right.getNaN());
    
    // Either is an infinity (wrong in the case of NaN but will be corrected)
    prop infCase( (left.isNegativeInf() && ITE(equality, prop(true), !right.isNegativeInf()) ) ||
		  (right.isPositiveInf() && ITE(equality, prop(true), !left.isPositiveInf()) ) ||
		  (ITE(equality, left.getInf() && right.getInf() && left.getSign() == right.getSign(), prop(false))) );


    // Either is a zero (wrong in the case of NaN but will be corrected)
    prop zeroCase( ( left.getZero() && !right.getZero() && !right.getSign()) ||
		   (right.getZero() && !left.getZero()  &&  left.getSign()) ||
		   (ITE(equality, left.getZero() && right.getZero(), prop(false))) );


    // Normal and subnormal case
    prop normalOrSubnormal(!left.getNaN()  && !right.getNaN() &&
			   !left.getInf()  && !right.getInf() &&
			   !left.getZero() && !right.getZero());

    prop negativeLessThanPositive(normalOrSubnormal && left.getSign() && !right.getSign());

    prop exponentNeeded(normalOrSubnormal && left.getSign() == right.getSign());
    probabilityAnnotation<t>(exponentNeeded, UNLIKELY);
    
    prop positiveCase(!left.getSign() && !right.getSign() &&
		      left.getExponent() < right.getExponent());
    prop negativeCase( left.getSign() &&  right.getSign() &&
		      left.getExponent() > right.getExponent());

    
    prop exponentEqual(left.getExponent() == right.getExponent());
    
    prop significandNeeded(exponentNeeded && exponentEqual);
    probabilityAnnotation<t>(significandNeeded, VERYUNLIKELY);

    prop positiveExEqCase(!left.getSign() && !right.getSign() &&
			  left.getSignificand() < right.getSignificand());
    prop negativeExEqCase( left.getSign() &&  right.getSign() &&
			   left.getSignificand() > right.getSignificand());

    prop positiveExEqCaseEq(!left.getSign() && !right.getSign() &&
			    left.getSignificand() <= right.getSignificand());
    prop negativeExEqCaseEq( left.getSign() &&  right.getSign() &&
			     left.getSignificand() >= right.getSignificand());

    return ITE(!normalOrSubnormal,
	       neitherNaN && (infCase || zeroCase),
	       ITE(!exponentNeeded,
		   negativeLessThanPositive,
		   ITE(!significandNeeded,
		       positiveCase || negativeCase,
		       ITE(equality,
			   positiveExEqCaseEq || negativeExEqCaseEq,
			   positiveExEqCase || negativeExEqCase))));
  }


  template <class t>
    typename t::prop lessThan (const typename t::fpt &format, 
			       const unpackedFloat<t> &left,
			       const unpackedFloat<t> &right) {
    PRECONDITION(left.valid(format));
    PRECONDITION(right.valid(format));

    typedef typename t::prop prop;

    return ordering(format, left, right, prop(false));
  }

  
  template <class t>
    typename t::prop lessThanOrEqual (const typename t::fpt &format, 
				      const unpackedFloat<t> &left,
				      const unpackedFloat<t> &right) {
    PRECONDITION(left.valid(format));
    PRECONDITION(right.valid(format));

    typedef typename t::prop prop;

    return ordering(format, left, right, prop(true));
  }


  // Note that IEEE-754 says that max(+0,-0) = +/-0 and max(-0,+0) = +/- 0
  template <class t>
  unpackedFloat<t> max (const typename t::fpt &format, 
			const unpackedFloat<t> &left,
			const unpackedFloat<t> &right,
			const typename t::prop &zeroCase) {
    return ITE(left.getNaN() || ordering(format, left, right, zeroCase),
	       right,
	       left);
  }

  // Note that IEEE-754 says that min(+0,-0) = +/-0 and min(-0,+0) = +/- 0
  // this will always return the left one.
  template <class t>
  unpackedFloat<t> min (const typename t::fpt &format, 
			const unpackedFloat<t> &left,
			const unpackedFloat<t> &right,
			const typename t::prop &zeroCase) {
    return ITE(right.getNaN() || ordering(format, left, right, zeroCase),
	       left,
	       right);
  }


  
  template <class t>
    typename t::prop originalLessThan (const typename t::fpt &format, 
				       const unpackedFloat<t> &left,
				       const unpackedFloat<t> &right) {

    typedef typename t::prop prop;

    PRECONDITION(left.valid(format));
    PRECONDITION(right.valid(format));

    // Optimisation : merge < and ==


    // All comparison with NaN are false
    prop neitherNan(!left.getNaN() && !right.getNaN());

    // Infinities are bigger than everything but themself
    prop eitherInf(left.getInf() || right.getInf());
    prop infCase(( left.isNegativeInf() && !right.isNegativeInf()) ||
		 (!left.isPositiveInf() &&  right.isPositiveInf()));


    // Both zero are equal
    prop eitherZero(left.getZero() || right.getZero());
    prop zeroCase(( left.getZero() && !right.getZero() && !right.getSign()) ||
		  (!left.getZero() &&   left.getSign() &&  right.getZero()));


    // Normal and subnormal

    prop negativeLessThanPositive(left.getSign() && !right.getSign());  // - < +
    prop positiveCase(!left.getSign() && !right.getSign() &&
		      ((left.getExponent() < right.getExponent()) ||
		       (left.getExponent() == right.getExponent() && 
			left.getSignificand() < right.getSignificand())));

    prop negativeCase(left.getSign() && right.getSign() &&
		      ((left.getExponent() > right.getExponent()) ||
		       (left.getExponent() == right.getExponent() && 
			left.getSignificand() > right.getSignificand())));
		 

    return neitherNan &&
      ITE(eitherInf,
	  infCase,
	  ITE(eitherZero,
	      zeroCase,
	      negativeLessThanPositive || positiveCase || negativeCase));
  }
  
  // Optimised combination of the two
  template <class t>
    typename t::prop originalLessThanOrEqual (const typename t::fpt &format, 
					      const unpackedFloat<t> &left,
					      const unpackedFloat<t> &right) {

    typedef typename t::prop prop;

    PRECONDITION(left.valid(format));
    PRECONDITION(right.valid(format));

    // Optimisation : merge < and ==


    // All comparison with NaN are false
    prop neitherNan(!left.getNaN() && !right.getNaN());

    // Infinities are bigger than everything but themself
    prop eitherInf(left.getInf() || right.getInf());
    prop infCase( (left.getInf() && right.getInf() && left.getSign() == right.getSign()) || 
		  left.isNegativeInf() ||
		  right.isPositiveInf());


    // Both zero are equal
    prop eitherZero(left.getZero() || right.getZero());
    prop zeroCase((left.getZero() && right.getZero()) ||
		  ( left.getZero() && !right.getSign()) ||
		  ( left.getSign() &&  right.getZero()));


    // Normal and subnormal

    prop negativeLessThanPositive(left.getSign() && !right.getSign());  // - < +
    prop positiveCase(!left.getSign() && !right.getSign() &&
		      ((left.getExponent() < right.getExponent()) ||
		       (left.getExponent() == right.getExponent() && 
			left.getSignificand() <= right.getSignificand())));

    prop negativeCase(left.getSign() && right.getSign() &&
		      ((left.getExponent() > right.getExponent()) ||
		       (left.getExponent() == right.getExponent() && 
			left.getSignificand() >= right.getSignificand())));
		 

    return neitherNan &&
      ITE(eitherInf,
	  infCase,
	  ITE(eitherZero,
	      zeroCase,
	      negativeLessThanPositive || positiveCase || negativeCase));
  }

}

#endif
