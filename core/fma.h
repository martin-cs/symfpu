/* Part of SymFPU, see LICENSE for licensing information */
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
   
   typedef typename t::bwt bwt;
   typedef typename t::prop prop;
   //typedef typename t::ubv ubv;
   typedef typename t::sbv sbv;
   typedef typename t::fpt fpt;
  
   PRECONDITION(leftMultiply.valid(format));
   PRECONDITION(rightMultiply.valid(format));
   PRECONDITION(addArgument.valid(format));

   /* First multiply */
   unpackedFloat<t> arithmeticMultiplyResult(arithmeticMultiply(format, leftMultiply, rightMultiply));

   sbv min(unpackedFloat<t>::minSubnormalExponent(format));
   sbv max(unpackedFloat<t>::maxNormalExponent(format));
   sbv multiplyResultExponentUpperBound(expandingAddWithCarryIn<t>(max, max, true));  // + 1 for renormalisation of the top bit
   sbv multiplyResultExponentLowerBound(expandingAddWithCarryIn<t>(min, min, false));
   INVARIANT(arithmeticMultiplyResult.wellFormed(multiplyResultExponentLowerBound, multiplyResultExponentUpperBound));

   // To meet the preconditions of other methods, we need to find a format
   // which is able to represent the value we have.
   // It is tempting to say this should be:
   //  fpt extendedFormat(format.exponentWidth() + 1, format.significandWidth() * 2);
   // And for many formats, this will work.
   // However the amount that that exponentWidth adds is not necessarily
   // the same at (e,s) and (e+1,2*s).  So we have to do this which
   // may result in a slightly wider exponent than we need.

   fpt extendedFormat(arithmeticMultiplyResult.getExponent().getWidth(), format.significandWidth() * 2);

   bwt currentExponentWidth = arithmeticMultiplyResult.getExponent().getWidth();
   bwt targetExponentWidth = unpackedFloat<t>::exponentWidth(extendedFormat);
   INVARIANT(targetExponentWidth >= currentExponentWidth);
   bwt extension = targetExponentWidth - currentExponentWidth;
   
   unpackedFloat<t> formattedArithmeticMultiplyResult(arithmeticMultiplyResult.getSign(), arithmeticMultiplyResult.getExponent().extend(extension), arithmeticMultiplyResult.getSignificand());   
   INVARIANT(formattedArithmeticMultiplyResult.valid(extendedFormat));

   

   /* Then add */
   
   // Rounding mode doesn't matter as this is a strict extension
   unpackedFloat<t> extendedAddArgument(convertFloatToFloat(format, extendedFormat, t::RTZ(), addArgument));

   prop knownInCorrectOrder(false);
   exponentCompareInfo<t> ec(addExponentCompare<t>(formattedArithmeticMultiplyResult.getExponent().getWidth() + 1,
						   formattedArithmeticMultiplyResult.getSignificand().getWidth(),
						   formattedArithmeticMultiplyResult.getExponent(),
						   extendedAddArgument.getExponent(),
						   knownInCorrectOrder));

   unpackedFloat<t> additionResult(arithmeticAdd(extendedFormat, roundingMode, formattedArithmeticMultiplyResult, extendedAddArgument, prop(true), knownInCorrectOrder, ec).uf);
   // Custom rounder flags are ignored as they are not applicable in this case

   // The invariants on this are tighter than you might think
   // In most formats the range of the multiply dominates
   //  if x = max exponent in the input format
   //     y = min exponent in the input format
   //  then the exponent of the product is in [2y,2x+1]
   //  so the exponent of the result is in [min(2y,max(2y,y)+a), max(2x+1,x)+b] allowing for 0's
   //  where a,b \in [-(p-1),1] but depend on the alignment of the two inputs
   //  for a,b to be negative you need that the exponents are equal or 1 apart
   //  for a,b to be +1 you need that the exponent of the lower one is within the length of the longest possible sets of leading 1's in a product.
   // This is a conservative choice of invariant
   INVARIANT(additionResult.wellFormed(multiplyResultExponentLowerBound.matchWidth(additionResult.getExponent()), multiplyResultExponentUpperBound.matchWidth(additionResult.getExponent())));


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
   // formattedArithmeticMultiplyResult is in extended format, thus we have to use a second rounder just for this case.
   // It is not zero, inf or NaN so it only matters when addArgument is zero when it would be returned.
   unpackedFloat<t> roundedMultiplyResult(rounder(format, roundingMode, formattedArithmeticMultiplyResult));

   unpackedFloat<t> fullMultiplyResult(addMultiplySpecialCases(format, leftMultiply, rightMultiply, roundedMultiplyResult.getSign(), roundedMultiplyResult));

   
   // We need the flags from the multiply special cases, determined on the arithemtic result,
   // i.e. handling special values and not the underflow / overflow of the result.
   // But we will use roundedMultiplyResult instead of the value so ...
   unpackedFloat<t> dummyZero(unpackedFloat<t>::makeZero(format, prop(false)));
   unpackedFloat<t> dummyValue(dummyZero.getSign(), dummyZero.getExponent(), dummyZero.getSignificand());

   unpackedFloat<t> multiplyResultWithSpecialCases(addMultiplySpecialCases(format, leftMultiply, rightMultiply, formattedArithmeticMultiplyResult.getSign(), dummyValue));

   
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
