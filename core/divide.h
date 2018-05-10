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
** divide.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 04/02/16
**
** Division of arbitrary precision floats
**
*/

#include "symfpu/core/unpackedFloat.h"
#include "symfpu/core/ite.h"
#include "symfpu/core/rounder.h"
#include "symfpu/core/operations.h"

#ifndef SYMFPU_DIVIDE
#define SYMFPU_DIVIDE

namespace symfpu {

template <class t>
  unpackedFloat<t> addDivideSpecialCases (const typename t::fpt &format,
					  const unpackedFloat<t> &left,
					  const unpackedFloat<t> &right,
					  const typename t::prop &sign,
					  const unpackedFloat<t> &divideResult) {
  typedef typename t::prop prop;

  prop eitherArgumentNaN(left.getNaN() || right.getNaN());
  prop generateNaN((left.getInf() && right.getInf()) ||
		   (left.getZero() && right.getZero()));
  
  prop isNaN(eitherArgumentNaN || generateNaN);

  prop isInf((!left.getZero() && right.getZero()) ||
	     (left.getInf() && !right.getInf()));

  prop isZero((!left.getInf() && right.getInf()) ||
	      (left.getZero() && !right.getZero()));

  return ITE(isNaN,
	     unpackedFloat<t>::makeNaN(format),
	     ITE(isInf,
		 unpackedFloat<t>::makeInf(format, sign),
		 ITE(isZero,
		     unpackedFloat<t>::makeZero(format, sign),
		     divideResult)));
 }


 template <class t>
  unpackedFloat<t> arithmeticDivide (const typename t::fpt &format,
				       const unpackedFloat<t> &left,
				       const unpackedFloat<t> &right) {
  typedef typename t::bwt bwt;
  typedef typename t::prop prop;
  typedef typename t::ubv ubv;
  typedef typename t::sbv sbv;
  typedef typename t::fpt fpt;

  PRECONDITION(left.valid(format));
  PRECONDITION(right.valid(format));

  // Compute sign
  prop divideSign(left.getSign() ^ right.getSign());

  // Subtract up exponents
  sbv exponentDiff(expandingSubtract<t>(left.getExponent(),right.getExponent()));
  // Optimisation : do this late and use the increment as a carry in

  sbv min(unpackedFloat<t>::minSubnormalExponent(format));
  sbv max(unpackedFloat<t>::maxNormalExponent(format));
  INVARIANT(expandingSubtract<t>(min,max) <= exponentDiff);
  INVARIANT(exponentDiff <= expandingSubtract<t>(max, min));
  // Optimisation : use the if-then-lazy-else to avoid dividing for underflow and overflow
  //                subnormal / greater-than-2^sigwidth does not need to be evaluated


  // Divide the significands
  // We need significandWidth() + 1 bits in the result but the top one may cancel, so add two bits
  ubv extendedNumerator(left.getSignificand().append(ubv::zero(2)));
  ubv extendedDenominator(right.getSignificand().append(ubv::zero(2)));

  resultWithRemainderBit<t> divided(fixedPointDivide<t>(extendedNumerator, extendedDenominator));
  

  bwt resWidth(divided.result.getWidth());
  ubv topBit(divided.result.extract(resWidth - 1, resWidth - 1));
  ubv nextBit(divided.result.extract(resWidth - 2, resWidth - 2));

  // Alignment of inputs means at least one of the two MSB is 1
  //  i.e. [1,2) / [1,2) = [0.5,2)
  // Top bit is set by the first round of the divide and thus is 50/50 1 or 0
  prop topBitSet(topBit.isAllOnes());
  INVARIANT(topBitSet || nextBit.isAllOnes());
  INVARIANT(topBitSet == (left.getSignificand() >= right.getSignificand()));
  

  // Re-align
  sbv alignedExponent(conditionalDecrement<t>(!topBitSet, exponentDiff)); // Will not overflow as previously expanded
  ubv alignedSignificand(conditionalLeftShiftOne<t>(!topBitSet, divided.result)); // Will not loose information

  // Create the sticky bit, it is important that this is after alignment
  ubv finishedSignificand(alignedSignificand | ubv(divided.remainderBit).extend(resWidth - 1));
  
  // Put back together
  unpackedFloat<t> divideResult(divideSign, alignedExponent.extend(1), finishedSignificand);

  // A brief word about formats.
  // You might think that the extend above is unnecessary : it is from a overflow point of view.
  // It's needed so that it is a valid number with exponentWidth() + 2.
  // +1 is sufficient in almost all cases.  However:
  //    very large normal / very small subnormal
  // can have an exponent greater than very large normal * 2 ( + 1)
  // because the exponent range is asymmetric with more subnormal than normal.
  
  fpt extendedFormat(format.exponentWidth() + 2, format.significandWidth() + 2);
  POSTCONDITION(divideResult.valid(extendedFormat));

  return divideResult;
 }


// Put it all together...
template <class t>
  unpackedFloat<t> divide (const typename t::fpt &format,
			   const typename t::rm &roundingMode,
			   const unpackedFloat<t> &left,
			   const unpackedFloat<t> &right) {
  //typedef typename t::bwt bwt;
  //typedef typename t::prop prop;
  //typedef typename t::ubv ubv;
  //typedef typename t::sbv sbv;

  PRECONDITION(left.valid(format));
  PRECONDITION(right.valid(format));

  unpackedFloat<t> divideResult(arithmeticDivide(format, left, right));
  
  unpackedFloat<t> roundedDivideResult(rounder(format, roundingMode, divideResult));
  
  unpackedFloat<t> result(addDivideSpecialCases(format, left, right, roundedDivideResult.getSign(), roundedDivideResult));

  POSTCONDITION(result.valid(format));

  return result;
 }


}

#endif
