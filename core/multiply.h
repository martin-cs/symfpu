/*
** Copyright (C) 2018 Martin Brain
**
** See the file LICENSE for licensing information.
*/

/*
** multiply.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 25/08/14
**
** Multiplication of arbitrary precision floats
**
*/

#include "symfpu/core/unpackedFloat.h"
#include "symfpu/core/ite.h"
#include "symfpu/core/rounder.h"
#include "symfpu/core/operations.h"

#ifndef SYMFPU_MULTIPLY
#define SYMFPU_MULTIPLY

namespace symfpu {

  // sign == multiplyResult.getSign() normally but not for FMA, thus an argument is needed
template <class t>
  unpackedFloat<t> addMultiplySpecialCases (const typename t::fpt &format,
					    const unpackedFloat<t> &left,
					    const unpackedFloat<t> &right,
					    const typename t::prop &sign,
					    const unpackedFloat<t> &multiplyResult) {
  typedef typename t::prop prop;

  prop eitherArgumentNan(left.getNaN() || right.getNaN());
  prop generateNan((left.getInf() && right.getZero()) ||
		   (left.getZero() && right.getInf()));
  prop isNan(eitherArgumentNan || generateNan);

  prop isInf(left.getInf() || right.getInf());

  prop isZero(left.getZero() || right.getZero());

  return ITE(isNan,
	     unpackedFloat<t>::makeNaN(format),
	     ITE(isInf,
		 unpackedFloat<t>::makeInf(format, sign),
		 ITE(isZero,
		     unpackedFloat<t>::makeZero(format, sign),
		     multiplyResult)));
 }

template <class t>
  unpackedFloat<t> arithmeticMultiply (const typename t::fpt &format,
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
  prop multiplySign(left.getSign() ^ right.getSign());

  // Multiply the significands
  ubv significandProduct(expandingMultiply<t>(left.getSignificand(), right.getSignificand()));
  // Optimisation : low bits are not needed apart from the guard and sticky bits
  // Optimisation : top bits accurately predict whether re-alignment is needed

  bwt spWidth(significandProduct.getWidth());
  ubv topBit(significandProduct.extract(spWidth - 1, spWidth - 1));
  ubv nextBit(significandProduct.extract(spWidth - 2, spWidth - 2));

  // Alignment of inputs means at least one of the two MSB is 1
  //  i.e. [1,2) * [1,2) = [1,4)
  // topBitSet is the likely case
  prop topBitSet(topBit.isAllOnes());
  INVARIANT(topBitSet || nextBit.isAllOnes());
  probabilityAnnotation<t>(topBitSet, LIKELY);

  // Re-align
  ubv alignedSignificand(conditionalLeftShiftOne<t>(!topBitSet, significandProduct)); // Will not loose information

  // Add up exponents
  #if 0
  sbv exponentSum(expandingAdd<t>(left.getExponent(),right.getExponent()));
  sbv min(unpackedFloat<t>::minSubnormalExponent(format));
  sbv max(unpackedFloat<t>::maxNormalExponent(format));
  INVARIANT(expandingAdd<t>(min,min) <= exponentSum);
  INVARIANT(exponentSum <= expandingAdd<t>(max, max));
  // Optimisation : use the if-then-lazy-else to avoid multiplying for underflow and overflow
  //                subnormal * subnormal does not need to be evaluated
  //                may be best done in the rounder along with underflow
  #endif
  
  sbv alignedExponent(expandingAddWithCarryIn<t>(left.getExponent(),right.getExponent(), topBitSet));

  
  // Put back together
  fpt extendedFormat(format.exponentWidth() + 1, format.significandWidth() * 2);
  unpackedFloat<t> multiplyResult(extendedFormat, multiplySign, alignedExponent, alignedSignificand);

  POSTCONDITION(multiplyResult.valid(extendedFormat));

  return multiplyResult;
 }


// Put it all together...
template <class t>
  unpackedFloat<t> multiply (const typename t::fpt &format,
			     const typename t::rm &roundingMode,
			     const unpackedFloat<t> &left,
			     const unpackedFloat<t> &right) {
  //typedef typename t::bwt bwt;
  //typedef typename t::prop prop;
  //typedef typename t::ubv ubv;
  //typedef typename t::sbv sbv;

  PRECONDITION(left.valid(format));
  PRECONDITION(right.valid(format));

  unpackedFloat<t> multiplyResult(arithmeticMultiply(format, left, right));
  
  unpackedFloat<t> roundedMultiplyResult(rounder(format, roundingMode, multiplyResult));
  
  unpackedFloat<t> result(addMultiplySpecialCases(format, left, right, roundedMultiplyResult.getSign(), roundedMultiplyResult));

  POSTCONDITION(result.valid(format));

  return result;
 }


}

#endif

