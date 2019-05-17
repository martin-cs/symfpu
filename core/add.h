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
** add.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 01/09/14
**
** Addition of arbitrary precision floats
**
** The current design is based on a two-path adder but it may be useful to use
** more paths.  There are five cases that are of interest:
**  1. effective add / very far
**     -> set the sticky bit only
**  2. effective add / far or near
**     -> align and add, realign down if needed
**  3. effective sub / very far
**     -> decrement, re-normalise and set sticky bits or (dependent on rounding-mode, skip entirely)
**  4. effective sub / far
**     -> align and subtract, realign up if needed
**  5. effective sub / near
**     -> align, subtract and normalise up
**
*/


/*
** Ideas
**  Optimisation : Collar the exponent difference, convert add to twice the width and thus unify the paths and simplify the shifting.
**  Enable the absolute max underapproximation
*/


#include "symfpu/core/unpackedFloat.h"
#include "symfpu/core/ite.h"
#include "symfpu/core/rounder.h"
#include "symfpu/core/sign.h"
#include "symfpu/core/operations.h"


#ifndef SYMFPU_ADD
#define SYMFPU_ADD

namespace symfpu {

  // There are a number of variants on how this should be done.
  // This is the implementation that handles all of them.
  // Below are restricted versions for particular special cases.
template <class t>
  unpackedFloat<t> addAdditionSpecialCasesComplete (const typename t::fpt &format,
						    const typename t::rm &roundingMode,
						    const unpackedFloat<t> &left,
						    const unpackedFloat<t> &leftID,
						    const unpackedFloat<t> &right,
						    const typename t::prop &returnLeft,
						    const typename t::prop &returnRight,
						    const unpackedFloat<t> &additionResult,
						    const typename t::prop &isAdd) {

  typedef typename t::prop prop;

  // NaN
  prop eitherArgumentNan(left.getNaN() || right.getNaN());
  prop bothInfinity(left.getInf() && right.getInf());
  prop signsMatch(left.getSign() == right.getSign());
  //prop compatableSigns(ITE(isAdd, signsMatch, !signsMatch));
  prop compatableSigns(isAdd ^ !signsMatch);

  prop generatesNaN(eitherArgumentNan || (bothInfinity && !compatableSigns));


  // Inf
  prop generatesInf((bothInfinity && compatableSigns) ||
		    ( left.getInf() && !right.getInf()) ||
		    (!left.getInf() &&  right.getInf()));

  prop signOfInf(ITE(left.getInf(), left.getSign(), prop(isAdd ^ !right.getSign())));

  
  // Zero
  prop bothZero(left.getZero() && right.getZero());
  prop flipRightSign(!isAdd ^ right.getSign());
  prop signOfZero(ITE((roundingMode == t::RTN()),
		      left.getSign() || flipRightSign,
		      left.getSign() && flipRightSign));

  prop  idLeft(!left.getZero() &&  right.getZero());
  prop idRight( left.getZero() && !right.getZero());

  // At most one of idLeft, idRight, generatesNaN, generatesInf and bothZero is true.
  // If used in addition additionResult is guaranteed to not be NaN.

  // Subtle trick : as the input to this will have been rounded it will have
  // an ITE with the default values "on top", thus doing the special cases
  // first (inner) rather than last (outer) allows them to be compacted better
  return ITE(idRight || returnRight,
	     ITE(isAdd,
		 right,
		 negate(format, right)),
	     ITE(idLeft || returnLeft,
		 leftID,
		 ITE(generatesNaN,
		     unpackedFloat<t>::makeNaN(format),
		     ITE(generatesInf,
			 unpackedFloat<t>::makeInf(format, signOfInf),
			 ITE(bothZero,
			     unpackedFloat<t>::makeZero(format, signOfZero),
			     additionResult)))));
 }


  // leftID is the value returned in the idLeft case (i.e. when left is not a
  // special number and right is zero).  This is needed by FMA as the flags
  // for left and leftID are computed differently and need to be handled differently.
template <class t>
  unpackedFloat<t> addAdditionSpecialCasesWithID (const typename t::fpt &format,
						  const typename t::rm &roundingMode,
						  const unpackedFloat<t> &left,
						  const unpackedFloat<t> &leftID,
						  const unpackedFloat<t> &right,
						  const unpackedFloat<t> &additionResult,
						  const typename t::prop &isAdd) {
  return addAdditionSpecialCasesComplete<t>(format, roundingMode, left, leftID,
					    right, typename t::prop(false), typename t::prop(false),
					    additionResult, isAdd);
  }


  // This is the usual case; use this one!
  template <class t>
  unpackedFloat<t> addAdditionSpecialCases (const typename t::fpt &format,
					    const typename t::rm &roundingMode,
					    const unpackedFloat<t> &left,
					    const unpackedFloat<t> &right,
					    const unpackedFloat<t> &additionResult,
					    const typename t::prop &isAdd) {
    return addAdditionSpecialCasesComplete<t>(format, roundingMode, left, left,
					      right, typename t::prop(false), typename t::prop(false),
					      additionResult, isAdd);
  }


  // As above but allows the (very) far path to be accelerated
  template <class t>
  unpackedFloat<t> addAdditionSpecialCasesWithBypass (const typename t::fpt &format,
						      const typename t::rm &roundingMode,
						      const unpackedFloat<t> &left,
						      const unpackedFloat<t> &right,
						      const typename t::prop &returnLeft,
						      const typename t::prop &returnRight,
						      const unpackedFloat<t> &additionResult,
						      const typename t::prop &isAdd) {
    return addAdditionSpecialCasesComplete<t>(format, roundingMode, left, left,
					      right, returnLeft, returnRight,
					      additionResult, isAdd);
  }



  /* Computes the normal / subnormal case only.
   * This allows multiple versions of the first phase to be used and
   * the first phase to be used for other things (e.g. FMA).
   */

template <class t>
struct exponentCompareInfo {
  typedef typename t::sbv sbv;
  typedef typename t::prop prop;

  prop leftIsMax;
  sbv maxExponent;
  sbv absoluteExponentDifference;
  prop diffIsZero;
  prop diffIsOne;
  prop diffIsGreaterThanPrecision;
  prop diffIsTwoToPrecision;
  prop diffIsGreaterThanPrecisionPlusOne;

  exponentCompareInfo(const prop &lil, const sbv &me, const sbv &aed,
		      const prop &diz, const prop &dio, const prop &digtp, const prop &dittp, const prop &disgtppo) :
    leftIsMax(lil), maxExponent(me), absoluteExponentDifference(aed),
      diffIsZero(diz), diffIsOne(dio), diffIsGreaterThanPrecision(digtp), diffIsTwoToPrecision(dittp), diffIsGreaterThanPrecisionPlusOne(disgtppo) {}

  exponentCompareInfo(const exponentCompareInfo<t> &old) :
    leftIsMax(old.leftIsMax), maxExponent(old.maxExponent), absoluteExponentDifference(old.absoluteExponentDifference),
    diffIsZero(old.diffIsZero), diffIsOne(old.diffIsOne),
    diffIsGreaterThanPrecision(old.diffIsGreaterThanPrecision),
    diffIsTwoToPrecision(old.diffIsTwoToPrecision),
    diffIsGreaterThanPrecisionPlusOne(old.diffIsGreaterThanPrecisionPlusOne) {}
};

template <class t>
  exponentCompareInfo<t> addExponentCompare(const typename t::bwt exponentWidth,
					    const typename t::bwt significandWidth,
					    const typename t::sbv &leftExponent,
					    const typename t::sbv &rightExponent,
					    const typename t::prop &knownInCorrectOrder) {
  PRECONDITION( leftExponent.getWidth() + 1 == exponentWidth);
  PRECONDITION(rightExponent.getWidth() + 1 == exponentWidth);
  
   typedef typename t::prop prop;
   typedef typename t::sbv sbv;

#if 0
   // The obvious implementation
   
   // Compute exponent distance
   sbv maxExponent(max<t,sbv>(leftExponent.extend(1), rightExponent.extend(1)));
   sbv minExponent(min<t,sbv>(leftExponent.extend(1), rightExponent.extend(1)));
   sbv absoluteExponentDifference(maxExponent - minExponent);

   prop leftIsMax(knownInCorrectOrder || leftExponent.extend(1) == maxExponent);
#else
   
   // A better implementation(?)
   sbv exponentDifference(leftExponent.extend(1) - rightExponent.extend(1));

   prop signBit(exponentDifference.toUnsigned().extract(exponentWidth - 1, exponentWidth - 1).isAllOnes());
   prop leftIsMax(knownInCorrectOrder || !signBit);

   sbv maxExponent(ITE(leftIsMax, leftExponent.extend(1), rightExponent.extend(1)));
   sbv absoluteExponentDifference(ITE(leftIsMax, exponentDifference, exponentDifference.modularNegate()));  // Largest negative value not obtainable so negate is safe
#endif

   INVARIANT(sbv::zero(exponentWidth) <= absoluteExponentDifference);

   // Optimisation : compact these comparisons at the bit-level
   prop diffIsZero(absoluteExponentDifference == sbv::zero(exponentWidth));
   prop diffIsOne(absoluteExponentDifference == sbv::one(exponentWidth));
   prop diffIsGreaterThanPrecision(sbv(exponentWidth, significandWidth) < absoluteExponentDifference);  // Assumes this is representable
   prop diffIsTwoToPrecision(!diffIsZero && !diffIsOne && !diffIsGreaterThanPrecision);
   prop diffIsGreaterThanPrecisionPlusOne(sbv(exponentWidth, significandWidth + 1) < absoluteExponentDifference);
  
   probabilityAnnotation<t,prop>(diffIsZero, UNLIKELY);
   probabilityAnnotation<t,prop>(diffIsOne, UNLIKELY);
   probabilityAnnotation<t,prop>(diffIsGreaterThanPrecision, LIKELY);  // In proving if not execution
   probabilityAnnotation<t,prop>(diffIsGreaterThanPrecisionPlusOne, LIKELY);  // In proving if not execution


   return exponentCompareInfo<t>(leftIsMax,
				 maxExponent, absoluteExponentDifference,
				 diffIsZero, diffIsOne, diffIsGreaterThanPrecision, diffIsTwoToPrecision, diffIsGreaterThanPrecisionPlusOne);
 }
  
  

// Note that the arithmetic part of add needs the rounding mode.
// This is an oddity due to the way that the sign of zero is generated.

 template <class t>
 struct floatWithCustomRounderInfo {
   unpackedFloat<t> uf;
   customRounderInfo<t> known;

 floatWithCustomRounderInfo(const unpackedFloat<t> &_uf, const customRounderInfo<t> &_known) : uf(_uf), known(_known) {}
 floatWithCustomRounderInfo(const floatWithCustomRounderInfo<t> &old) : uf(old.uf), known(old.known) {}
 };

 template <class t>
   floatWithCustomRounderInfo<t> arithmeticAdd (const typename t::fpt &format,
						const typename t::rm &roundingMode,
						const unpackedFloat<t> &left,
						const unpackedFloat<t> &right,
						const typename t::prop &isAdd,
						const typename t::prop &knownInCorrectOrder,
						const exponentCompareInfo<t> &ec) {
   
   typedef typename t::bwt bwt;
   typedef typename t::fpt fpt;
   typedef typename t::prop prop;
   typedef typename t::ubv ubv;
   typedef typename t::sbv sbv;
   
   PRECONDITION(left.valid(format));
   PRECONDITION(right.valid(format));

   // Work out if an effective subtraction
   prop effectiveAdd((left.getSign() ^ right.getSign()) ^ isAdd);
   
   bwt exponentWidth(left.getExponent().getWidth() + 1);
   bwt significandWidth(left.getSignificand().getWidth());
   
   /* Exponent difference and effective add implies a large amount about the output exponent and flags
   ** R denotes that this is possible via rounding up and incrementing the exponent
   **
   ** Case       A. max(l,r) + 1,    B. max(l,r)    C. max(l,r) - 1      D. max(l,r) - k       E. zero
   ** Eff. Add      Y                   Y
   **  diff = 0     Y, sticky 0 
   **  diff = 1     Y, sticky 0, R      Y, sticky 0
   **  diff : [2,p] decreasing prob., R Y
   **  diff > p     R                   Y
   **
   ** Eff. Sub                          Y              Y                    Y, exact              Y, exact
   **  diff = 0                                        Y, exact             prob. drop with k     low prob.
   **  diff = 1                         Y, sticky 0    Y, exact             prob. drop with k
   **  diff : [2,p]                     Y, R           decreasing prob.
   **  diff > p                         Y, R           low prob.
   **
   */
   // Optimisation : add 'bypass' invariants that link the final exponent to the input using this.
   
   
   // Rounder flags
   prop noOverflow(!effectiveAdd);
   prop noUnderflow(true);
   //prop exact(); // Need to see if it cancels, see below
   prop subnormalExact(true);
   prop noSignificandOverflow((effectiveAdd && ec.diffIsZero) ||
			      (!effectiveAdd && (ec.diffIsZero || ec.diffIsOne)));

   prop stickyBitIsZero(ec.diffIsZero || ec.diffIsOne);


   // Work out ordering
   prop leftLarger(knownInCorrectOrder ||
		   (ec.leftIsMax &&
		    ITE(!ec.diffIsZero,
			prop(true),
			left.getSignificand() >= right.getSignificand())));
   // Optimisation : can we avoid this comparison completely and allow the result to be negative?
   // This may be hard as a compare is cheaper than a negate after, particularly as there has to be an ITE here

   // Extend the significands to give room for carry plus guard and sticky bits
   ubv lsig((ITE(leftLarger, left.getSignificand(), right.getSignificand())).extend(1).append(ubv::zero(2)));
   ubv ssig((ITE(leftLarger, right.getSignificand(), left.getSignificand())).extend(1).append(ubv::zero(2)));

   prop resultSign(ITE(leftLarger,
		       left.getSign(),
		       prop(!isAdd ^ right.getSign())));

   // Extended so no info lost, negate before shift so that sign-extension works
   ubv negatedSmaller(conditionalNegate<t,ubv,prop>(!effectiveAdd, ssig));

   ubv shiftAmount(ec.absoluteExponentDifference.toUnsigned() // Safe as >= 0
		   .resize(negatedSmaller.getWidth()));  // Safe as long as the significand has more bits than the exponent
   INVARIANT(exponentWidth <= significandWidth);


   // Shift the smaller significand
   stickyRightShiftResult<t> shifted(stickyRightShift<t>(negatedSmaller, shiftAmount));
   
   ubv negatedAlignedSmaller(ITE(ec.diffIsGreaterThanPrecisionPlusOne, // Fast path the common case, +1 to avoid issues with the guard bit
				 ITE(effectiveAdd,
				      ubv::zero(negatedSmaller.getWidth()),
				     ~ubv::zero(negatedSmaller.getWidth())),				     
				 shifted.signExtendedResult));
   ubv shiftedStickyBit(ITE(ec.diffIsGreaterThanPrecision,
			    ubv::one(negatedSmaller.getWidth()),
			    shifted.stickyBit));  // Have to separate otherwise align up may convert it to the guard bit
   
   
   // Sum and re-align
   ubv sum(lsig.modularAdd(negatedAlignedSmaller));

   bwt sumWidth(sum.getWidth());
   ubv topBit(sum.extract(sumWidth - 1, sumWidth - 1));
   ubv alignedBit(sum.extract(sumWidth - 2, sumWidth - 2));
   ubv lowerBit(sum.extract(sumWidth - 3, sumWidth - 3));

   
   prop overflow(!(topBit.isAllZeros()));
   prop cancel(topBit.isAllZeros() && alignedBit.isAllZeros());
   prop minorCancel(cancel && lowerBit.isAllOnes());
   prop majorCancel(cancel && lowerBit.isAllZeros());
   prop fullCancel(majorCancel && sum.isAllZeros());
   
   probabilityAnnotation<t,prop>(overflow, UNLIKELY);
   probabilityAnnotation<t,prop>(cancel, UNLIKELY);
   probabilityAnnotation<t,prop>(minorCancel, UNLIKELY);
   probabilityAnnotation<t,prop>(majorCancel, VERYUNLIKELY);
   probabilityAnnotation<t,prop>(fullCancel, VERYUNLIKELY);

   INVARIANT(IMPLIES(effectiveAdd && ec.diffIsZero, overflow));
   INVARIANT(IMPLIES(overflow, effectiveAdd && (!ec.diffIsGreaterThanPrecision))); // That case can only overflow by rounding
   INVARIANT(IMPLIES(cancel, !effectiveAdd));
   INVARIANT(IMPLIES(majorCancel, ec.diffIsZero || ec.diffIsOne));
   
   probabilityAnnotation<t,prop>(overflow && ec.diffIsTwoToPrecision, UNLIKELY);
   probabilityAnnotation<t,prop>(cancel && ec.diffIsTwoToPrecision, UNLIKELY);
   probabilityAnnotation<t,prop>(cancel && ec.diffIsGreaterThanPrecision, VERYUNLIKELY);
   
   prop exact(cancel && (ec.diffIsZero || ec.diffIsOne)); // For completeness
   
   ubv alignedSum(conditionalLeftShiftOne<t,ubv,prop>(minorCancel,
						      conditionalRightShiftOne<t,ubv,prop>(overflow, sum)));

   sbv exponentCorrectionTerm(ITE(minorCancel,
				  -sbv::one(exponentWidth),
				  ITE(overflow,
				      sbv::one(exponentWidth),
				      sbv::zero(exponentWidth))));
   
   sbv correctedExponent(ec.maxExponent + exponentCorrectionTerm); // Safe due to extension

   // Watch closely...
   ubv stickyBit(ITE(stickyBitIsZero || majorCancel,
		     ubv::zero(alignedSum.getWidth()),
		     (shiftedStickyBit | ITE(!overflow, ubv::zero(1), sum.extract(0,0)).extend(alignedSum.getWidth() - 1))));
   
   
   // Put it back together
   unpackedFloat<t> sumResult(resultSign, correctedExponent, (alignedSum | stickyBit).contract(1));

   // We return something in an extended format
   //  *. One extra exponent bit to deal with the 'overflow' case
   //  *. Two extra significand bits for the guard and sticky bits
   fpt extendedFormat(format.exponentWidth() + 1, format.significandWidth() + 2);
   
   // Deal with the major cancellation case
   // It would be nice to use normaliseUpDetectZero but the sign
   // of the zero depends on the rounding mode.
   unpackedFloat<t> additionResult(ITE(fullCancel,
				       unpackedFloat<t>::makeZero(extendedFormat, roundingMode == t::RTN()),
				       ITE(majorCancel,
					   sumResult.normaliseUp(extendedFormat),
					   sumResult)));
   
   // Some thought is required here to convince yourself that 
   // there will be no subnormal values that violate this.
   // See 'all subnormals generated by addition are exact'
   // and the extended exponent.
   POSTCONDITION(additionResult.valid(extendedFormat));
   
   return floatWithCustomRounderInfo<t>(additionResult, customRounderInfo<t>(noOverflow, noUnderflow, exact, subnormalExact, noSignificandOverflow));   
 }
 
 template <class t>
   unpackedFloat<t> dualPathArithmeticAdd (const typename t::fpt &format,
					   const typename t::rm &roundingMode,
					   const unpackedFloat<t> &left,
					   const unpackedFloat<t> &right,
					   const typename t::prop &isAdd) {
   
   typedef typename t::bwt bwt;
   typedef typename t::fpt fpt;
   typedef typename t::prop prop;
   typedef typename t::ubv ubv;
   typedef typename t::sbv sbv;
   
   PRECONDITION(left.valid(format));
   PRECONDITION(right.valid(format));

   // We return something in an extended format
   //  *. One extra exponent bit to deal with the 'overflow' case
   //  *. Two extra significand bits for the guard and sticky bits
   fpt extendedFormat(format.exponentWidth() + 1, format.significandWidth() + 2);


   // Compute exponent difference and swap the two arguments if needed
   sbv initialExponentDifference(expandingSubtract<t>(left.getExponent(), right.getExponent()));
   bwt edWidth(initialExponentDifference.getWidth());
   sbv edWidthZero(sbv::zero(edWidth));
   prop orderingCorrect( initialExponentDifference >  edWidthZero ||
			(initialExponentDifference == edWidthZero &&
			 left.getSignificand() >= right.getSignificand()));

   unpackedFloat<t> larger(ITE(orderingCorrect, left, right));
   unpackedFloat<t> smaller(ITE(orderingCorrect, right, left));
   sbv exponentDifference(ITE(orderingCorrect,
			      initialExponentDifference,
			      -initialExponentDifference));

   prop resultSign(ITE(orderingCorrect,
		       left.getSign(),
		       prop(!isAdd ^ right.getSign())));


   // Work out if an effective subtraction
   prop effectiveAdd(larger.getSign() ^ smaller.getSign() ^ isAdd);

   
   // Extend the significands to give room for carry plus guard and sticky bits
   ubv lsig( larger.getSignificand().extend(1).append(ubv::zero(2)));
   ubv ssig(smaller.getSignificand().extend(1).append(ubv::zero(2)));


   // This is a two-path adder, so determine which of the two paths to use 
   // The near path is only needed for things that can cancel more than one bit
   prop farPath(exponentDifference > sbv::one(edWidth) || effectiveAdd);


   // Far path : Align
   ubv negatedSmaller(ITE(effectiveAdd, ssig, ssig.modularNegate())); // Extended so no info lost
                                                                      // Negate before shift so that sign-extension works

   sbv significandWidth(edWidth, lsig.getWidth());
   prop noOverlap(exponentDifference > significandWidth);

   ubv shiftAmount(exponentDifference.toUnsigned() // Safe as >= 0
		   .resize(ssig.getWidth()));      // This looses information but the case in which it does is handles by noOverlap


   ubv negatedAlignedSmaller(negatedSmaller.signExtendRightShift(shiftAmount));
   ubv shiftedStickyBit(rightShiftStickyBit<t>(negatedSmaller, shiftAmount));  // Have to separate otherwise align up may convert it to the guard bit

   // Far path : Sum and re-align
   ubv sum(lsig.modularAdd(negatedAlignedSmaller));

   bwt sumWidth(sum.getWidth());
   ubv topBit(sum.extract(sumWidth - 1, sumWidth - 1));
   ubv centerBit(sum.extract(sumWidth - 2, sumWidth - 2));

   prop noOverflow(topBit.isAllZeros()); // Only correct if effectiveAdd is set
   prop noCancel(centerBit.isAllOnes());


   
   // TODO : Add invariants

   ubv alignedSum(ITE(effectiveAdd,
		      ITE(noOverflow,
			  sum,
			  (sum >> ubv::one(sumWidth)) | (sum & ubv::one(sumWidth))),  // Cheap sticky right shift
		      ITE(noCancel,
			  sum,
			  sum.modularLeftShift(ubv::one(sumWidth))))); // In the case when this looses data, the result is not used

   sbv extendedLargerExponent(larger.getExponent().extend(1));  // So that increment and decrement don't overflow
   sbv correctedExponent(ITE(effectiveAdd,
			     ITE(noOverflow,
				 extendedLargerExponent,
				 extendedLargerExponent.increment()),
			     ITE(noCancel,
				 extendedLargerExponent,
				 extendedLargerExponent.decrement())));

   // Far path : Construct result
   unpackedFloat<t> farPathResult(resultSign, correctedExponent, (alignedSum | shiftedStickyBit).contract(1));




   // Near path : Align
   prop exponentDifferenceAllZeros(exponentDifference.isAllZeros());
   ubv nearAlignedSmaller(ITE(exponentDifferenceAllZeros, ssig, ssig >> ubv::one(ssig.getWidth())));


   // Near path : Sum and realign
   ubv nearSum(lsig - nearAlignedSmaller);
   // Optimisation : the two paths can be merged up to here to give a pseudo-two path encoding

   prop fullCancel(nearSum.isAllZeros());
   prop nearNoCancel(nearSum.extract(sumWidth - 2, sumWidth - 2).isAllOnes());

   ubv choppedNearSum(nearSum.extract(sumWidth - 3,1)); // In the case this is used, cut bits are all 0 
   unpackedFloat<t> cancellation(resultSign, 
				 larger.getExponent().decrement(),
				 choppedNearSum);


   // Near path : Construct result
   unpackedFloat<t> nearPathResult(resultSign, extendedLargerExponent, nearSum.contract(1));



   // Bring the paths together
   // Optimisation : fix the noOverlap / very far path for directed rounding modes
   unpackedFloat<t> additionResult(ITE(farPath,
				       /* ITE(noOverlap,
					   ITE((isAdd || orderingCorrect),
					       larger,
					       negate(format, larger)),
					       farPathResult), */
				       farPathResult,
				       ITE(fullCancel,
					   unpackedFloat<t>::makeZero(extendedFormat, roundingMode == t::RTN()),
					   ITE(nearNoCancel,
					       nearPathResult,
					       cancellation.normaliseUp(format).extend(1,2)))));
   
   // Some thought is required here to convince yourself that 
   // there will be no subnormal values that violate this.
   // See 'all subnormals generated by addition are exact'
   // and the extended exponent.
   POSTCONDITION(additionResult.valid(extendedFormat));
   
   return additionResult;
 }




 template <class t>
   unpackedFloat<t> dualPathAdd (const typename t::fpt &format,
				 const typename t::rm &roundingMode,
				 const unpackedFloat<t> &left,
				 const unpackedFloat<t> &right,
				 const typename t::prop &isAdd) {
   
   PRECONDITION(left.valid(format));
   PRECONDITION(right.valid(format));

   unpackedFloat<t> additionResult(dualPathArithmeticAdd(format, roundingMode, left, right, isAdd));

   unpackedFloat<t> roundedAdditionResult(rounder(format, roundingMode, additionResult));

   unpackedFloat<t> result(addAdditionSpecialCases(format, roundingMode, left, right, roundedAdditionResult, isAdd));
   
   POSTCONDITION(result.valid(format));
   
   return result;
 }

template <class t>
   unpackedFloat<t> add (const typename t::fpt &format,
			 const typename t::rm &roundingMode,
			 const unpackedFloat<t> &left,
			 const unpackedFloat<t> &right,
			 const typename t::prop &isAdd) {

   //typedef typename t::bwt bwt;
   typedef typename t::prop prop;
   //typedef typename t::ubv ubv;
   //typedef typename t::sbv sbv;

   PRECONDITION(left.valid(format));
   PRECONDITION(right.valid(format));

   // Optimisation : add a flag which assumes that left and right are in the correct order
   prop knownInCorrectOrder(false);

   exponentCompareInfo<t> ec(addExponentCompare<t>(left.getExponent().getWidth() + 1, left.getSignificand().getWidth(),
						   left.getExponent(), right.getExponent(), knownInCorrectOrder));

   floatWithCustomRounderInfo<t> additionResult(arithmeticAdd(format, roundingMode, left, right, isAdd, knownInCorrectOrder, ec));

   unpackedFloat<t> roundedAdditionResult(customRounder(format, roundingMode, additionResult.uf, additionResult.known));

   unpackedFloat<t> result(addAdditionSpecialCases(format, roundingMode, left, right, roundedAdditionResult, isAdd));

   POSTCONDITION(result.valid(format));

   return result;
 }

 template <class t>
   unpackedFloat<t> addWithBypass (const typename t::fpt &format,
				   const typename t::rm &roundingMode,
				   const unpackedFloat<t> &left,
				   const unpackedFloat<t> &right,
				   const typename t::prop &isAdd) {
   
   //typedef typename t::bwt bwt;
   typedef typename t::prop prop;
   //typedef typename t::ubv ubv;
   //typedef typename t::sbv sbv;
   
   PRECONDITION(left.valid(format));
   PRECONDITION(right.valid(format));

   // Optimisation : add a flag which assumes that left and right are in the correct order
   prop knownInCorrectOrder(false);


   exponentCompareInfo<t> ec(addExponentCompare<t>(left.getExponent().getWidth() + 1, left.getSignificand().getWidth(),
						   left.getExponent(), right.getExponent(), knownInCorrectOrder));

   floatWithCustomRounderInfo<t> additionResult(arithmeticAdd(format, roundingMode, left, right, isAdd, knownInCorrectOrder, ec));

   unpackedFloat<t> roundedAdditionResult(customRounder(format, roundingMode, additionResult.uf, additionResult.known));

   // In the "very far path", i.e. exponent difference is greater than significand length + 1
   // addition becomes max(left,right) or max(left,right)+/-1 ULP.  This is rare in execution
   // (to the point of being a software quality issue) but common in theorem proving.
   // Given we have to have cases for "return left" and "return right" to handle zeros,
   // we might as well make use of these cases to handle when addition behaves like max...
   // Note that this is possible but more complex with just diffIsGreaterThanPrecision.

   prop enableBypass(ec.diffIsGreaterThanPrecisionPlusOne &&
		      !left.getNaN() &&  !left.getInf() &&  !left.getZero() && // Handle this as special cases
		     !right.getNaN() && !right.getInf() && !right.getZero());


   // Duplication but easier to recompute than to pass
   prop effectiveAdd((left.getSign() ^ right.getSign()) ^ isAdd);
   prop resultSign(ITE((knownInCorrectOrder || ec.leftIsMax),  // CAUTION : only true in the enableBypass case!
		       left.getSign(),
		       prop(!isAdd ^ right.getSign())));

   prop significandEven(true);  // This is an optimisation that assumes only RNE uses this bit
                                // This needs to be changed to implement things like roundToOdd
                                // or for the diffIsGreaterThanPrecision case.
   prop farRoundUp(roundingDecision<t>(roundingMode, resultSign, significandEven, !effectiveAdd, prop(true), prop(false)));

   // Returns left or right unchanged if adding and rounded down or subtracting and rounded up
   prop roundInCorrectDirection(effectiveAdd ^ farRoundUp);

   prop  returnLeft(enableBypass &&  ec.leftIsMax && roundInCorrectDirection);
   prop returnRight(enableBypass && !ec.leftIsMax && roundInCorrectDirection);

   unpackedFloat<t> result(addAdditionSpecialCasesWithBypass(format, roundingMode, left, right, returnLeft, returnRight, roundedAdditionResult, isAdd));
   
   POSTCONDITION(result.valid(format));
   
   return result;
 }


 // True if and only if adding these would result in a catastrophic cancellation
 // I.E. if the addition cancells out cancelAmount or more MSBs leaving only LSBs
 template <class t>
   typename t::prop isCatastrophicCancellation (const typename t::fpt &format,
						const unpackedFloat<t> &left,
						const unpackedFloat<t> &right,
						const typename t::bwt &cancelAmount,
						const typename t::prop &isAdd) {

   typedef typename t::bwt bwt;
   typedef typename t::prop prop;
   typedef typename t::ubv ubv;
   //typedef typename t::sbv sbv;

   PRECONDITION(left.valid(format));
   PRECONDITION(right.valid(format));
   PRECONDITION(cancelAmount >= 2);   // cancel 0 is not meaningful
                                      // cancel 1 is common on subtract and arguably not an error
   PRECONDITION(cancelAmount <= format.significandWidth());  // can't cancel more than you have

   // 1. It has to be an effective subtraction
   // Duplication but easier to recompute than to pass
   prop effectiveAdd((left.getSign() ^ right.getSign()) ^ isAdd);

   // 2. It must be either normal or subnormal numbers
   prop leftSpecial(left.getNaN() || left.getInf() || left.getZero());
   prop rightSpecial(right.getNaN() || right.getInf() || right.getZero());

   // 3.A. exponents are equal and so are the leading cancelAmount bits
   // 3.B. exponent diff is one and the smaller one is 11111, larger one is 10000

   // Optimisation : add a flag which assumes that left and right are in the correct order
   prop knownInCorrectOrder(false);
   exponentCompareInfo<t> ec(addExponentCompare<t>(left.getExponent().getWidth() + 1, left.getSignificand().getWidth(),
						   left.getExponent(), right.getExponent(), knownInCorrectOrder));

   // Can ignore the MSB of the significand as by invariants this is always 1
   bwt significandWidth(format.significandWidth());
   bwt topBit(significandWidth - 2);
   bwt bottomBit(significandWidth - cancelAmount);

   ubv leftExtract(left.getSignificand.extract(topBit, bottomBit));
   ubv rightExtract(right.getSignificand.extract(topBit, bottomBit));

   prop result(ITE(!effectiveAdd && !leftSpecial && !rightSpecial,
		   ITE(ec.diffIsZero,
		       leftExtract == rightExtract,
		       ITE(ec.diffIsOne,
			   ITE(ec.leftIsMax,
			        leftExtract.isAllZeros() && rightExtract.isAllOnes(),
			       rightExtract.isAllZeros() &&  leftExtract.isAllOnes()),
			   false)),
		   false));

   return result;
 }

}

#endif

