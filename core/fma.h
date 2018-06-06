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
** fma.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 20/05/15
**
** Fused multiply and add :
**  fma(R,A,B,C) = round(R, A * B + C)
**
*/


#include "symfpu/core/unpackedFloat.h"
#include "symfpu/core/ite.h"
#include "symfpu/core/rounder.h"
#include "symfpu/core/multiply.h"
#include "symfpu/core/convert.h"
#include "symfpu/core/add.h"


#ifndef SYMFPU_FMA
#define SYMFPU_FMA

namespace symfpu {

 template <class t>
   unpackedFloat<t> fma (const typename t::fpt &format,
			 const typename t::rm &roundingMode,
			 const unpackedFloat<t> &leftMultiply,
			 const unpackedFloat<t> &rightMultiply,
			 const unpackedFloat<t> &addArgument) {
   
   //   typedef typename t::bwt bwt;
   typedef typename t::prop prop;
   //typedef typename t::ubv ubv;
   //typedef typename t::sbv sbv;
   typedef typename t::fpt fpt;
  
   PRECONDITION(leftMultiply.valid(format));
   PRECONDITION(rightMultiply.valid(format));
   PRECONDITION(addArgument.valid(format));

   /* First multiply */
   unpackedFloat<t> arithmeticMultiplyResult(arithmeticMultiply(format, leftMultiply, rightMultiply));
   
   fpt extendedFormat(format.exponentWidth() + 1, format.significandWidth() * 2);
   INVARIANT(arithmeticMultiplyResult.valid(extendedFormat));

   

   /* Then add */
   
   // Rounding mode doesn't matter as this is a strict extension
   unpackedFloat<t> extendedAddArgument(convertFloatToFloat(format, extendedFormat, t::RTZ(), addArgument));

   prop knownInCorrectOrder(false);
   exponentCompareInfo<t> ec(addExponentCompare<t>(arithmeticMultiplyResult.getExponent().getWidth() + 1,
						   arithmeticMultiplyResult.getSignificand().getWidth(),
						   arithmeticMultiplyResult.getExponent(),
						   extendedAddArgument.getExponent(),
						   knownInCorrectOrder));

   unpackedFloat<t> additionResult(arithmeticAdd(extendedFormat, roundingMode, arithmeticMultiplyResult, extendedAddArgument, prop(true), knownInCorrectOrder, ec).uf);
   // Custom rounder flags are ignored as they are not applicable in this case

   fpt evenMoreExtendedFormat(extendedFormat.exponentWidth() + 1, extendedFormat.significandWidth() + 2);
   INVARIANT(additionResult.valid(evenMoreExtendedFormat));


   /* Then round */
   
   unpackedFloat<t> roundedResult(rounder(format, roundingMode, additionResult));
   INVARIANT(roundedResult.valid(format));
   
   // This result is correct as long as neither of multiplyResult or extendedAddArgument is
   // 0, Inf or NaN.  Note that roundedResult may be zero from cancelation or underflow
   // or infinity due to rounding. If it is, it has the correct sign.



   /* Finally, the special cases */
   
   // One disadvantage to having a flag for zero and default exponents and significands for zero
   // that are not (min, 0) is that the x + (+/-)0 case has to be handled by the addition special cases.
   // This means that you need the value of x, rounded to the correct format.
   // arithmeticMultiplyResult is in extended format, thus we have to use a second rounder just for this case.
   // It is not zero, inf or NaN so it only matters when addArgument is zero when it would be returned.
   unpackedFloat<t> roundedMultiplyResult(rounder(format, roundingMode, arithmeticMultiplyResult));

   unpackedFloat<t> fullMultiplyResult(addMultiplySpecialCases(format, leftMultiply, rightMultiply, roundedMultiplyResult.getSign(), roundedMultiplyResult));

   
   // We need the flags from the multiply special cases, determined on the arithemtic result,
   // i.e. handling special values and not the underflow / overflow of the result.
   // But we will use roundedMultiplyResult instead of the value so ...
   unpackedFloat<t> dummyZero(unpackedFloat<t>::makeZero(format, prop(false)));
   unpackedFloat<t> dummyValue(dummyZero.getSign(), dummyZero.getExponent(), dummyZero.getSignificand());

   unpackedFloat<t> multiplyResultWithSpecialCases(addMultiplySpecialCases(format, leftMultiply, rightMultiply, arithmeticMultiplyResult.getSign(), dummyValue));

   
   unpackedFloat<t> result(addAdditionSpecialCasesWithID(format,
							 roundingMode,
							 multiplyResultWithSpecialCases,
							 fullMultiplyResult, // for the identity case
							 addArgument,
							 roundedResult,
							 prop(true)));
   
   POSTCONDITION(result.valid(format));
   
   return result;
 }

/*
 * BUGS : 
 * 1. sign of zero different for exact 0 and underflow
 * 2. large * -large  + inf  = inf  not NaN
 * 3. rounder decision bugs : one looks like an issue with too-eager overflow,
 *    one looks like a misplaced decision on highest subnormal exponent
 */
 
 template <class t>
   unpackedFloat<t> fmaBroken (const typename t::fpt &format,
			 const typename t::rm &roundingMode,
			 const unpackedFloat<t> &leftMultiply,
			 const unpackedFloat<t> &rightMultiply,
			 const unpackedFloat<t> &addArgument) {
   
   //   typedef typename t::bwt bwt;
   typedef typename t::prop prop;
   //typedef typename t::ubv ubv;
   //typedef typename t::sbv sbv;
   typedef typename t::fpt fpt;
  
   PRECONDITION(leftMultiply.valid(format));
   PRECONDITION(rightMultiply.valid(format));
   PRECONDITION(addArgument.valid(format));

   unpackedFloat<t> multiplyResult(arithmeticMultiply(format, leftMultiply, rightMultiply));
   
   fpt extendedFormat(format.exponentWidth() + 1, format.significandWidth() * 2);
   INVARIANT(multiplyResult.valid(extendedFormat));

   // Rounding mode doesn't matter as this is a strict extension
   unpackedFloat<t> extendedAddArgument(convertFloatToFloat(format, extendedFormat, t::RTZ(), addArgument));

   unpackedFloat<t> additionResult(arithmeticAdd(extendedFormat, roundingMode, multiplyResult, extendedAddArgument, prop(true), prop(false)).uf);
   // Custom rounder flags are ignored as they are not applicable in this case
   
   unpackedFloat<t> roundedResult(rounder(format, roundingMode, additionResult));
     
   // Note that multiplyResult.getSign() != roundedResult.getSign() in rare cases
   // the multiply special cases use the sign for zeros and infinities, thus the sign of the
   // result of the multiplication is needed (i.e. the xor of the sign of left and right multiply)
   // (-small, +inf, large) should trigger this as the desired result is -inf
   // but roundedResult.getSign() is positive.
   unpackedFloat<t> roundedResultWithMultiplyCases(addMultiplySpecialCases(format,
									   leftMultiply,
									   rightMultiply,
									   multiplyResult.getSign(),
									   roundedResult));

   
   // One disadvantage to having a flag for zero and default exponents and significands for zero
   // that are not (min, 0) is that the x + 0 case has to be handled by the addition special cases.
   // This means that you need the value of x, rounded to the correct format.
   // multiplyResult is in extended format, thus we have to use a second rounder just for this case.
   // It is not zero, inf or NaN so it only matters when addArgument is zero when it would be returned.

   unpackedFloat<t> roundedMultiplyResult(rounder(format, roundingMode, multiplyResult));
   // Optimisation : Try ITE before rounding so that only one rounder is needed


   // To make matters more awkward, we also need to apply the multiplicative special cases so that
   // (x*0) + y is correctly handled by the addition special cases.  Without applying the
   // multiplicative ones, (x*0) would not be correctly flagged as 0.
   unpackedFloat<t> roundedMultiplyResultWithMultiplyCases(addMultiplySpecialCases(format,
										   leftMultiply,
										   rightMultiply,
										   multiplyResult.getSign(),
										   roundedMultiplyResult));
   // Optimisation : consolidate the special cases and verify against this
   
   unpackedFloat<t> result(addAdditionSpecialCases(format,
						   roundingMode,
						   roundedMultiplyResultWithMultiplyCases,
						   addArgument,
						   roundedResultWithMultiplyCases,
						   prop(true)));
   
   POSTCONDITION(result.valid(format));
   
   return result;
 }


}

#endif
