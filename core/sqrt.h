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
** sqrt.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 05/02/16
**
** Square root of arbitrary precision floats
**
*/

#include "symfpu/core/unpackedFloat.h"
#include "symfpu/core/ite.h"
#include "symfpu/core/rounder.h"
#include "symfpu/core/operations.h"

#ifndef SYMFPU_SQRT
#define SYMFPU_SQRT

namespace symfpu {

template <class t>
  unpackedFloat<t> addSqrtSpecialCases (const typename t::fpt &format,
					  const unpackedFloat<t> &uf,
					  const typename t::prop &sign,
					  const unpackedFloat<t> &sqrtResult) {
  typedef typename t::prop prop;

  prop generateNaN(uf.getSign() && !uf.getZero());
  prop isNaN(uf.getNaN() || generateNaN);

  prop isInf(uf.getInf() && !uf.getSign());

  prop isZero(uf.getZero());

  return ITE(isNaN,
	     unpackedFloat<t>::makeNaN(format),
	     ITE(isInf,
		 unpackedFloat<t>::makeInf(format, prop(false)),
		 ITE(isZero,
		     unpackedFloat<t>::makeZero(format, sign),
		     sqrtResult)));
 }


 template <class t>
  unpackedFloat<t> arithmeticSqrt (const typename t::fpt &format,
				       const unpackedFloat<t> &uf) {
  typedef typename t::bwt bwt;
  typedef typename t::prop prop;
  typedef typename t::ubv ubv;
  typedef typename t::sbv sbv;
  typedef typename t::fpt fpt;

  PRECONDITION(uf.valid(format));

  // Compute sign
  prop sqrtSign(uf.getSign());

  // Divide the exponent by 2
  sbv exponent(uf.getExponent());
  bwt exponentWidth(exponent.getWidth());
  prop exponentEven((exponent & sbv::one(exponentWidth)).isAllZeros());
  #if 0
  sbv exponentHalved(conditionalDecrement<t,sbv,prop>((exponent < sbv::zero(exponentWidth)) && !exponentEven,
						      exponent.signExtendRightShift(sbv::one(exponentWidth))));
  #endif
  sbv exponentHalved(exponent.signExtendRightShift(sbv::one(exponentWidth)));
  // Right shift rounds down for positive, and away for negative  (-5 >>> 1 == -3)
  //  sqrt(1.s * 2^{-(2n + 1)}) = sqrt(1.s * 2^{-2n - 2 + 1)})
  //                            = sqrt(1.s * 2^{-2(n + 1)} * 2) 
  //                            = sqrt(1.s * 2) * 2^{-(n + 1)}
  // Optimisation : improve the encoding of this operation
  
  // Sqrt the significands
  //  extend to allow alignment, pad so result has a guard bit
  ubv alignedSignificand(conditionalLeftShiftOne<t,ubv,prop>(!exponentEven, uf.getSignificand().extend(1).append(ubv::zero(1))));

  resultWithRemainderBit<t> sqrtd(fixedPointSqrt<t>(alignedSignificand));
  

  bwt resWidth(sqrtd.result.getWidth());
  ubv topBit(sqrtd.result.extract(resWidth - 1, resWidth - 1));
  ubv guardBit(sqrtd.result.extract(0,0));

  // Alignment of inputs means it is the range [1,4) so the result is in [1,2)
  // Also, the square root cannot be exactly between two numbers
  INVARIANT(topBit.isAllOnes());
  INVARIANT(IMPLIES(guardBit.isAllOnes(), sqrtd.remainderBit));
  // This also implies that no alignment of the exponent is needed

  ubv finishedSignificand(sqrtd.result.append(ubv(sqrtd.remainderBit)));

  unpackedFloat<t> sqrtResult(sqrtSign, exponentHalved, finishedSignificand);

  
  fpt extendedFormat(format.exponentWidth(), format.significandWidth() + 2);
  // format.exponentWidth() - 1 should also be true but requires shrinking the exponent and
  // then increasing it in the rounder
  POSTCONDITION(sqrtResult.valid(extendedFormat));

  return sqrtResult;
 }


// Put it all together...
template <class t>
  unpackedFloat<t> sqrt (const typename t::fpt &format,
			   const typename t::rm &roundingMode,
			   const unpackedFloat<t> &uf) {
  //typedef typename t::bwt bwt;
    typedef typename t::prop prop;
  //typedef typename t::ubv ubv;
  //typedef typename t::sbv sbv;

  PRECONDITION(uf.valid(format));

  unpackedFloat<t> sqrtResult(arithmeticSqrt(format, uf));

  // Exponent is divided by two, thus it can't overflow, underflow or generate a subnormal number.
  // The last one is quite subtle but you can show that the largest number generatable
  // by arithmeticSqrt is 111...111:0:1 with the last two as the guard and sticky bits.
  // Round up (when the sign is positive) and round down (when the sign is negative --
  // the result will be computed but then discarded) are the only cases when this can increment the significand.
  customRounderInfo<t> cri(prop(true), prop(true), prop(false), prop(true),
			   !((roundingMode == t::RTP() && !sqrtResult.getSign()) ||
			     (roundingMode == t::RTN() &&  sqrtResult.getSign())));
  unpackedFloat<t> roundedSqrtResult(customRounder(format, roundingMode, sqrtResult, cri));
  
  unpackedFloat<t> result(addSqrtSpecialCases(format, uf, roundedSqrtResult.getSign(), roundedSqrtResult));

  POSTCONDITION(result.valid(format));

  return result;
 }


}

#endif
