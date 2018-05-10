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
** rounder.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 26/08/14
**
** Rounding arbitrary length unpacked floats back to their correct length.
**
*/

/*
 * Rounder Ideas
 *
 *  - get the rounder to handle things with 1 or 2 possible leading
 *    0's rather than normalise before and after rounding.
 *    Subnormal rounding should ideally be temporarily denormalised to
 *    avoid the expensive variable position rounding.
 *
 *  1. Take exponent and add a 0 at the start to deal with overflow.
 *     (may not be needed if the exponent has already been extended due
 *      to subnormals, addition, etc.)
 *  2. For each of the possible positions of the binary point, produce
 *     an incremented version and round and sticky bits.
 *     (given the alignment is known, may be able to just pick how much
 *      to increment by and then increment once.)
 *  3. The leading zeros of the original exponent pick which set of
 *     round and sticky to use.  Thus pick which of the incremented /
 *     normal to use and align apart from 1.
 *  4. Finally align, most of the information needed for this will be
 *     already known.
 *  (Catch -- this can introduce diamonds)
 *
 * - can fold the increment from the rounder into the circuitry / 
 *   last round of accumulation.
 *   In fact, can take this further.  a + b + 1 takes the same
 *   circuitry as a + b.  a * b + 1 shares most of the circuit with
 *   a * b and neither of them require a sign extension.  Extending
 *   this to a + b + k and a * b + k for positional rounding is an 
 *   open question.
 *
 * - add a 'non-deterministic rounding' mode for underapproximation.
 *
 * - add 'round-to-odd' and 'round-away-from-zero' mode
 *
 * - add 'flush subnormals to zero' option
 *
 * - Rather than increment and re-align, take all but the top bit of the
 *   significand, concatinate on to the exponent and then increment.
 *   This is effectively using the ordering trick for packed rounding.
 *
 * - The sticky bit should be a flag rather than concatinated on to the
 *   number to make arithmetic reasoning easier.
 *   (Conversely ... having it as bits is better for semi-normalised things)
 *
 * - Specialised rounders
 *    Additive
 *     Only one overflow increment is possible, even for a + b + 1
 *     When a subnormal number if generated, it is exact and thus you
 *      don't need to think about positional rounding for the
 *      subnormals.
 *     No need to check for underflow to zero for similar reasons.
 *     The rounder for the near path can be even more specialised as 
 *      the sticky bit is always 0 and it can't overflow.
 *     Can only overflow if it is an effective add.
 *
 *    Multiplicative
 *     May be able to increment without overflow.
 *      It only seems to be possible for subnormals.
 *     Don't need to extend the exponent to increment when normalising.
 *     Can use injection bits during the multiplication reduction tree,
 *      the problem with the carry up can be fixed afterwards with an
 *      increment (add for subnormal) or by spotting the carry up during
 *      the reduction.
 *
 *    FMA
 *     More like multiply.
 */


#include "symfpu/core/operations.h"
#include "symfpu/core/unpackedFloat.h"

#ifndef SYMFPU_ROUNDER
#define SYMFPU_ROUNDER

namespace symfpu {

  // The final reconstruction of the rounded result
  // Handles the overflow and underflow conditions
  template <class t>
  unpackedFloat<t> rounderSpecialCases (const typename t::fpt &format,
					const typename t::rm &roundingMode,
					const unpackedFloat<t> &roundedResult,
					const typename t::prop &overflow,
					const typename t::prop &underflow,
					const typename t::prop &isZero)
  {
    typedef typename t::prop prop;
    typedef typename t::ubv ubv;

    /*** Underflow and overflow ***/
    
    // On overflow either return inf or max
    prop returnInf(roundingMode == t::RNE() || 
		   roundingMode == t::RNA() ||
		   (roundingMode == t::RTP() && !roundedResult.getSign()) ||
		   (roundingMode == t::RTN() &&  roundedResult.getSign()));
    probabilityAnnotation<t>(returnInf, LIKELY);  // Inf is more likely than max in most application scenarios
    
    // On underflow either return 0 or minimum subnormal
    prop returnZero(roundingMode == t::RNE() || 
		    roundingMode == t::RNA() ||
		    roundingMode == t::RTZ() ||
		    (roundingMode == t::RTP() &&  roundedResult.getSign()) ||
		    (roundingMode == t::RTN() && !roundedResult.getSign()));
    probabilityAnnotation<t>(returnZero, LIKELY);   // 0 is more likely than min in most application scenarios


    
    /*** Reconstruct ***/
    unpackedFloat<t> inf(unpackedFloat<t>::makeInf(format, roundedResult.getSign()));
    unpackedFloat<t> max(roundedResult.getSign(), unpackedFloat<t>::maxNormalExponent(format), ubv::allOnes(unpackedFloat<t>::significandWidth(format)));
    unpackedFloat<t> min(roundedResult.getSign(), unpackedFloat<t>::minSubnormalExponent(format), unpackedFloat<t>::leadingOne(unpackedFloat<t>::significandWidth(format)));
    unpackedFloat<t> zero(unpackedFloat<t>::makeZero(format, roundedResult.getSign()));
    
    unpackedFloat<t> result(ITE(isZero, 
				zero,
				ITE(underflow,
				    ITE(returnZero, zero, min),
				    ITE(overflow,
					ITE(returnInf, inf, max),
					roundedResult))));
    return result;
  }

  
  // Decide whether to round up or not
  template <class t>
    typename t::prop roundingDecision (const typename t::rm &roundingMode,
				     const typename t::prop &sign,
				     const typename t::prop &significandEven,
				     const typename t::prop &guardBit,
				     const typename t::prop &stickyBit,
				     const typename t::prop &knownRoundDown) {
    typedef typename t::prop prop;

    prop roundUpRNE(roundingMode == t::RNE() && guardBit && (stickyBit || !significandEven));
    prop roundUpRNA(roundingMode == t::RNA() && guardBit);
    prop roundUpRTP(roundingMode == t::RTP() && !sign && (guardBit || stickyBit));
    prop roundUpRTN(roundingMode == t::RTN() &&  sign && (guardBit || stickyBit));
    prop roundUpRTZ(roundingMode == t::RTZ() && prop(false));
    prop roundUp(!knownRoundDown &&
		 (roundUpRNE || roundUpRNA || roundUpRTP || roundUpRTN || roundUpRTZ));

    return roundUp;
  }



  template <class t>
  struct significandRounderResult{
    typename t::ubv significand;
    typename t::prop incrementExponent;
      
    significandRounderResult(const typename t::ubv &sig, const typename t::prop &inc) :
      significand(sig), incrementExponent(inc) {}
      
    significandRounderResult(const significandRounderResult &old) :
      significand(old.significand), incrementExponent(old.incrementExponent) {}
  };

  // Handles rounding the significand to a fixed width
  // If knownRoundDown is true should simplify to just extract
  // Not quite the same as either rounder so can't quite be refactored
  template <class t>
  significandRounderResult<t> fixedPositionRound(const typename t::rm &roundingMode,
						 const typename t::prop &sign,
						 const typename t::ubv &significand,
						 const typename t::bwt &targetWidth,
						 const typename t::prop &knownLeadingOne,
						 const typename t::prop &knownRoundDown) {
    typedef typename t::bwt bwt;
    typedef typename t::prop prop;
    typedef typename t::ubv ubv;

    bwt sigWidth(significand.getWidth());
    PRECONDITION(sigWidth >= targetWidth + 2);
    // Extract
    ubv extractedSignificand(significand.extract(sigWidth - 1, sigWidth - targetWidth).extend(1)); // extended to catch the overflow

    prop significandEven(extractedSignificand.extract(0,0).isAllZeros());

    
    // Normal guard and sticky bits
    bwt guardBitPosition(sigWidth - (targetWidth + 1));
    prop guardBit(significand.extract(guardBitPosition, guardBitPosition).isAllOnes());

    prop stickyBit(!significand.extract(guardBitPosition - 1,0).isAllZeros());

    // Rounding decision
    prop roundUp(roundingDecision<t>(roundingMode, sign, significandEven,
				     guardBit, stickyBit, knownRoundDown));

    // Conditional increment
    ubv roundedSignificand(conditionalIncrement<t>(roundUp, extractedSignificand));

    ubv overflowBit(roundedSignificand.extract(targetWidth, targetWidth) & ubv(roundUp));
    ubv carryUpMask((overflowBit | ubv(knownLeadingOne)).append(ubv::zero(targetWidth - 1)));   // Cheaper than conditional shift
    
    // Build result
    significandRounderResult<t> result(roundedSignificand.extract(targetWidth-1,0) | carryUpMask,
				    overflowBit.isAllOnes());

    return result;
  }

  
  // Handles rounding the significand to a fixed width
  // If knownRoundDown is true should simplify to just mask
  // Not quite the same as either rounder so can't quite be refactored
  template <class t>
  significandRounderResult<t> variablePositionRound(const typename t::rm &roundingMode,
						    const typename t::prop &sign,
						    const typename t::ubv &significand,
						    const typename t::ubv &roundPosition,
						    const typename t::prop &knownLeadingOne,
						    const typename t::prop &knownRoundDown) {
    typedef typename t::bwt bwt;
    typedef typename t::prop prop;
    typedef typename t::ubv ubv;

    bwt sigWidth(significand.getWidth());
    
    // Set up significand
    // Round-up-from-sticky bit and overflow bit at MSB, (fall-back) guard and sticky bits at LSB
    ubv expandedSignificand(significand.extend(2).append(ubv::zero(2)));
    bwt exsigWidth(expandedSignificand.getWidth());

    
    // Identify the increment, guard and sticky bits
    ubv incrementLocation(ubv(exsigWidth, (0x1U << 2U)) << roundPosition.matchWidth(expandedSignificand));
    ubv guardLocation(incrementLocation >> ubv::one(exsigWidth));
    ubv stickyLocations(guardLocation.decrement());

    prop significandEven((incrementLocation & expandedSignificand).isAllZeros());
    prop guardBit(!(guardLocation & expandedSignificand).isAllZeros());
    prop stickyBit(!(stickyLocations & expandedSignificand).isAllZeros());

    // Rounding decision
    prop roundUp(roundingDecision<t>(roundingMode, sign, significandEven,
				     guardBit, stickyBit, knownRoundDown));

    // Conditional increment
    ubv roundedSignificand(expandedSignificand + ITE(roundUp,
						     incrementLocation,
						     ubv::zero(exsigWidth)));
    
    // Mask out rounded bits and extract
    ubv maskedRoundedSignificand(roundedSignificand & (~(stickyLocations << ubv::one(exsigWidth)))); // LSB is wrong but gets cut

    ubv roundUpFromSticky(roundedSignificand.extract(exsigWidth - 1, exsigWidth - 1));  // Only true when rounding up and whole significand is sticky
    ubv overflowBit(roundedSignificand.extract(exsigWidth - 2, exsigWidth - 2));
    ubv maskTrigger((roundUpFromSticky | overflowBit) & ubv(roundUp));
    ubv carryUpMask((maskTrigger | ubv(knownLeadingOne)).append(ubv::zero(sigWidth - 1)));   // Cheaper than conditional shift
    
    // Build result
    significandRounderResult<t> result(maskedRoundedSignificand.extract(sigWidth + 1, 2) | carryUpMask,
				       maskTrigger.isAllOnes());

    return result;
  }

  
  
  // Allows various of the key branches in the rounder to be fixed / removed
  // using information that is known from the operation in which the rounder
  // is used. Setting these to false gives the usual rounder.
  template <class t>
  struct customRounderInfo {
    typedef typename t::prop prop;

    prop noOverflow;
    prop noUnderflow;
    prop exact;                 // Significand does not need to be changed
    prop subnormalExact;        // If the value is subnormal then it is exact
    prop noSignificandOverflow; // Incrementing the significand will not cause overflow
    //prop stickyIsZero;        // Can be fixed in the number directly
    
    customRounderInfo (const prop &noO, const prop &noU,
		       const prop &e, const prop &sE,
		       const prop &nSO) :
      noOverflow(noO), noUnderflow(noU), exact(e),
      subnormalExact(sE), noSignificandOverflow(nSO) {}
  };
  
template <class t>
  unpackedFloat<t> customRounder (const typename t::fpt &format,
				  const typename t::rm &roundingMode,
				  const unpackedFloat<t> &uf,
				  const customRounderInfo<t> &known) {

  typedef typename t::bwt bwt;
  typedef typename t::prop prop;
  typedef typename t::ubv ubv;
  typedef typename t::sbv sbv;

  //PRECONDITION(uf.valid(format));
  // Not a precondition because
  //  1. Exponent and significand may be extended.
  //  2. Their values may also be outside of the correct range.
  //
  // However some thing do hold:
  //  1. Leading bit of the significant is 1   (if you want to get a meaningful answer)
  //     (if not, a cancellation on the near path of add can cause
  //      this, but the result will not be used, so it can be incorrect)
  ubv psig(uf.getSignificand());
  bwt sigWidth(psig.getWidth());
  //PRECONDITION(psig.extract(sigWidth-1, sigWidth-1).isAllOnes());
  ubv sig(psig | unpackedFloat<t>::leadingOne(sigWidth));

  //  2. Must have round and sticky bits
  bwt targetSignificandWidth(unpackedFloat<t>::significandWidth(format));
  PRECONDITION(sigWidth >= targetSignificandWidth + 2);

  //  3. Must have at least enough exponent bits
  sbv exp(uf.getExponent());
  bwt expWidth(exp.getWidth());
  bwt targetExponentWidth(unpackedFloat<t>::exponentWidth(format));
  PRECONDITION(expWidth >= targetExponentWidth);

  // Also, do not round special values.
  // Note that this relies on these having default values and 
  // the code that calls the rounder constructing uf from parts.
  //PRECONDITION(!uf.getNaN());
  //PRECONDITION(!uf.getInf());
  //PRECONDITION(!uf.getZero());
  // The safe choice of default values means this should work OK



  /*** Early underflow and overflow detection ***/
  bwt exponentExtension(expWidth - targetExponentWidth);
  prop earlyOverflow(exp > unpackedFloat<t>::maxNormalExponent(format).extend(exponentExtension));
  prop earlyUnderflow(exp < unpackedFloat<t>::minSubnormalExponent(format).extend(exponentExtension).decrement());
  // Optimisation : if the precondition on sigWidth and targetSignificandWidth is removed
  //                then can change to:
  //                   exponent >= minSubnormalExponent - 1
  //                && sigWidth > targetSignificandBits
  probabilityAnnotation<t>(earlyOverflow, UNLIKELY);  // (over,under)flows are generally rare events
  probabilityAnnotation<t>(earlyUnderflow, UNLIKELY);

  prop potentialLateOverflow(exp == unpackedFloat<t>::maxNormalExponent(format).extend(exponentExtension));
  prop potentialLateUnderflow(exp == unpackedFloat<t>::minSubnormalExponent(format).extend(exponentExtension).decrement());
  probabilityAnnotation<t>(potentialLateOverflow, VERYUNLIKELY);
  probabilityAnnotation<t>(potentialLateUnderflow, VERYUNLIKELY);
  


  /*** Normal or subnormal rounding? ***/
  prop normalRoundingRange(exp >= unpackedFloat<t>::minNormalExponent(format).extend(exponentExtension));
  probabilityAnnotation<t>(normalRoundingRange, LIKELY);
  prop normalRounding(normalRoundingRange || known.subnormalExact);

  

  /*** Round to correct significand. ***/
  ubv extractedSignificand(sig.extract(sigWidth - 1, sigWidth - targetSignificandWidth).extend(1)); // extended to catch the overflow

  // Normal guard and sticky bits
  bwt guardBitPosition(sigWidth - (targetSignificandWidth + 1));
  prop guardBit(sig.extract(guardBitPosition, guardBitPosition).isAllOnes());

  prop stickyBit(!sig.extract(guardBitPosition - 1,0).isAllZeros());
  

  // For subnormals, locating the guard and stick bits is a bit more involved
  //sbv subnormalAmount(uf.getSubnormalAmount(format)); // Catch is, uf isn't in the given format, so this doesn't work
  sbv subnormalAmount(expandingSubtract<t>(unpackedFloat<t>::minNormalExponent(format).matchWidth(exp),exp));
  INVARIANT((subnormalAmount < sbv(expWidth + 1, sigWidth - 1)) || earlyUnderflow);
  // Note that this is negative if normal, giving a full subnormal mask
  // but the result will be ignored (see the next invariant)

  ubv subnormalShiftPrepared(subnormalAmount.toUnsigned().matchWidth(extractedSignificand));

  // Compute masks
  ubv subnormalMask(orderEncode<t>(subnormalShiftPrepared)); // Invariant implies this if all ones, it will not be used
  ubv subnormalStickyMask(subnormalMask >> ubv::one(targetSignificandWidth + 1)); // +1 as the exponent is extended

  // Apply
  ubv subnormalMaskedSignificand(extractedSignificand & (~subnormalMask));
  ubv subnormalMaskRemoved(extractedSignificand & subnormalMask);
  // Optimisation : remove the masking with a single orderEncodeBitwise style construct
  
  prop subnormalGuardBit(!(subnormalMaskRemoved & (~subnormalStickyMask)).isAllZeros());
  prop subnormalStickyBit(guardBit || stickyBit || 
			  !((subnormalMaskRemoved & subnormalStickyMask).isAllZeros()));


  ubv subnormalIncrementAmount((subnormalMask.modularLeftShift(ubv::one(targetSignificandWidth + 1))) & ~subnormalMask); // The only case when this looses info is earlyUnderflow
  INVARIANT(IMPLIES(subnormalIncrementAmount.isAllZeros(), earlyUnderflow || normalRounding));
  

  // Have to choose the right one dependent on rounding mode
  prop choosenGuardBit(ITE(normalRounding, guardBit, subnormalGuardBit));
  prop choosenStickyBit(ITE(normalRounding, stickyBit, subnormalStickyBit));
  
  prop significandEven(ITE(normalRounding,
			   extractedSignificand.extract(0,0).isAllZeros(),
			   ((extractedSignificand & subnormalIncrementAmount).isAllZeros())));
  prop roundUp(roundingDecision<t>(roundingMode, uf.getSign(), significandEven,
				   choosenGuardBit, choosenStickyBit,
				   known.exact || (known.subnormalExact && !normalRoundingRange)));


  // Perform the increment as needed
  ubv leadingOne(unpackedFloat<t>::leadingOne(targetSignificandWidth));
  // Not actually true, consider minSubnormalExponent - 1 : not an early underfow and empty significand
  //INVARIANT(!(subnormalMaskedSignificand & leadingOne).isAllZeros() ||
  //          earlyUnderflow); // This one really matters, it means only the early underflow path is wrong
  

  // Convert the round up flag to a mask
  ubv normalRoundUpAmount(ubv(roundUp).matchWidth(extractedSignificand));
  ubv subnormalRoundUpMask(ubv(roundUp).append(ubv::zero(targetSignificandWidth)).signExtendRightShift(ubv(targetSignificandWidth + 1, targetSignificandWidth)));
  ubv subnormalRoundUpAmount(subnormalRoundUpMask & subnormalIncrementAmount);
  
  ubv rawRoundedSignificand((ITE(normalRounding,
				 extractedSignificand,
				 subnormalMaskedSignificand)
			     +
			     ITE(normalRounding,
				 normalRoundUpAmount,
				 subnormalRoundUpAmount)));
  
  // We might have lost the leading one, if so, re-add and note that we need to increment the significand
  prop significandOverflow(rawRoundedSignificand.extract(targetSignificandWidth, targetSignificandWidth).isAllOnes());
  INVARIANT(IMPLIES(significandOverflow, roundUp));
  
  ubv extractedRoundedSignificand(rawRoundedSignificand.extract(targetSignificandWidth - 1, 0));
  ubv roundedSignificand(extractedRoundedSignificand | leadingOne);
  INVARIANT(IMPLIES(significandOverflow, extractedRoundedSignificand.isAllZeros()));

  


  /*** Round to correct exponent. ***/

  // The extend is almost certainly unnecessary (see specialised rounders)
  sbv extendedExponent(exp.extend(1));

  prop incrementExponentNeeded(roundUp && significandOverflow);  // The roundUp is implied but kept for signal forwarding
  probabilityAnnotation<t>(incrementExponentNeeded, VERYUNLIKELY);
  prop incrementExponent(!known.noSignificandOverflow && incrementExponentNeeded);
  INVARIANT(IMPLIES(known.noSignificandOverflow, !incrementExponentNeeded));
  
  sbv correctedExponent(conditionalIncrement<t>(incrementExponent, extendedExponent));

  // Track overflows and underflows
  sbv maxNormal(unpackedFloat<t>::maxNormalExponent(format).matchWidth(correctedExponent));
  sbv minSubnormal(unpackedFloat<t>::minSubnormalExponent(format).matchWidth(correctedExponent));
  
  sbv correctedExponentInRange(collar<t>(correctedExponent, minSubnormal, maxNormal));

  
  // This can over and underflow but these values will not be used
  bwt currentExponentWidth(correctedExponentInRange.getWidth());
  sbv roundedExponent(correctedExponentInRange.contract(currentExponentWidth - targetExponentWidth));


  
  /*** Finish ***/

  prop computedOverflow(potentialLateOverflow && incrementExponentNeeded);
  prop computedUnderflow(potentialLateUnderflow && !incrementExponentNeeded);
  probabilityAnnotation<t>(computedOverflow, UNLIKELY);
  probabilityAnnotation<t>(computedUnderflow, UNLIKELY);

  prop lateOverflow(!earlyOverflow && computedOverflow);
  prop lateUnderflow(!earlyUnderflow && computedUnderflow);
  probabilityAnnotation<t>(lateOverflow, VERYUNLIKELY);
  probabilityAnnotation<t>(lateUnderflow, VERYUNLIKELY);

  
  // So that ITE abstraction works...
  prop overflow(!known.noOverflow && ITE(lateOverflow, prop(true), earlyOverflow));
  prop underflow(!known.noUnderflow && ITE(lateUnderflow, prop(true), earlyUnderflow));
  
  unpackedFloat<t> roundedResult(uf.getSign(), roundedExponent, roundedSignificand);
  unpackedFloat<t> result(rounderSpecialCases<t>(format, roundingMode, roundedResult,
						 overflow, underflow, uf.getZero()));
					      
  POSTCONDITION(result.valid(format));

  return result;
 }

template <class t>
unpackedFloat<t> originalRounder (const typename t::fpt &format,
				  const typename t::rm &roundingMode,
				  const unpackedFloat<t> &uf) {

  typedef typename t::bwt bwt;
  typedef typename t::prop prop;
  typedef typename t::ubv ubv;
  typedef typename t::sbv sbv;

  //PRECONDITION(uf.valid(format));
  // Not a precondition because
  //  1. Exponent and significand may be extended.
  //  2. Their values may also be outside of the correct range.
  //
  // However some thing do hold:
  //  1. Leading bit of the significant is 1   (if you want to get a meaningful answer)
  //     (if not, a cancellation on the near path of add can cause
  //      this, but the result will not be used, so it can be incorrect)
  ubv psig(uf.getSignificand());
  bwt sigWidth(psig.getWidth());
  //PRECONDITION(psig.extract(sigWidth-1, sigWidth-1).isAllOnes());
  ubv sig(psig | unpackedFloat<t>::leadingOne(sigWidth));

  //  2. Must have round and sticky bits
  bwt targetSignificandWidth(unpackedFloat<t>::significandWidth(format));
  PRECONDITION(sigWidth >= targetSignificandWidth + 2);

  //  3. Must have at least enough exponent bits
  sbv exp(uf.getExponent());
  bwt expWidth(exp.getWidth());
  bwt targetExponentWidth(unpackedFloat<t>::exponentWidth(format));
  PRECONDITION(expWidth >= targetExponentWidth);

  // Also, do not round special values.
  // Note that this relies on these having default values and 
  // the code that calls the rounder constructing uf from parts.
  PRECONDITION(!uf.getNaN());
  PRECONDITION(!uf.getInf());
  //PRECONDITION(!uf.getZero());   // The safe choice of default values means this should work OK



  /*** Early underflow and overflow detection ***/
  bwt exponentExtension(expWidth - targetExponentWidth);
  prop earlyOverflow(exp > unpackedFloat<t>::maxNormalExponent(format).extend(exponentExtension));
  prop earlyUnderflow(exp < unpackedFloat<t>::minSubnormalExponent(format).extend(exponentExtension).decrement());
  // Optimisation : if the precondition on sigWidth and targetSignificandWidth is removed
  //                then can change to:
  //                   exponent >= minSubnormalExponent - 1
  //                && sigWidth > targetSignificandBits



  /*** Normal or subnormal rounding? ***/
  prop normalRounding(exp >= unpackedFloat<t>::minNormalExponent(format).extend(exponentExtension));
  probabilityAnnotation<t>(normalRounding, LIKELY);
  

  /*** Round to correct significand. ***/
  ubv extractedSignificand(sig.extract(sigWidth - 1, sigWidth - targetSignificandWidth));

  bwt guardBitPosition(sigWidth - (targetSignificandWidth + 1));
  prop guardBit(sig.extract(guardBitPosition, guardBitPosition).isAllOnes());

  prop stickyBit(!sig.extract(guardBitPosition - 1,0).isAllZeros());
  

  // For subnormals, locating the guard and stick bits is a bit more involved
  sbv subnormalAmount(unpackedFloat<t>::maxSubnormalExponent(format).extend(exponentExtension) - exp);
  prop belowLimit(subnormalAmount <= sbv::zero(expWidth));    // Not subnormal
  prop aboveLimit(subnormalAmount >= sbv(expWidth, targetSignificandWidth));  // Will underflow
  sbv subnormalShift(ITE((belowLimit || aboveLimit), sbv::zero(expWidth), subnormalAmount));
  // Optimisation : collar
  
  ubv subnormalShiftPrepared(subnormalShift.toUnsigned().extend(targetSignificandWidth - expWidth));
  ubv guardLocation(ubv::one(targetSignificandWidth) << subnormalShiftPrepared);
  ubv stickyMask(guardLocation.decrement());


  prop subnormalGuardBit(!(extractedSignificand & guardLocation).isAllZeros());
  prop subnormalStickyBit(guardBit || stickyBit || 
			  !((extractedSignificand & stickyMask).isAllZeros()));




  // Can overflow but is handled
  ubv incrementedSignificand(extractedSignificand.modularIncrement());
  prop incrementedSignificandOverflow(incrementedSignificand.isAllZeros());
  // Optimisation : conditional increment
  // Optimisation : use top bit of significand to increment the exponent

  ubv correctedIncrementedSignificand(ITE(!incrementedSignificandOverflow,
					  incrementedSignificand,
					  unpackedFloat<t>::leadingOne(unpackedFloat<t>::significandWidth(format))));


  ubv incrementAmount(guardLocation.modularLeftShift(ubv::one(guardLocation.getWidth()))); // Overflows (safely) in the case of rounding up to the least subnormal.
  ubv mask(guardLocation | stickyMask);
  ubv maskedSignificand(extractedSignificand & ~mask);

  ubv subnormalIncrementedSignificand(maskedSignificand.modularAdd(incrementAmount));
  prop subnormalIncrementedSignificandOverflow(subnormalIncrementedSignificand.isAllZeros());
  ubv subnomalCorrectedIncrementedSignificand(ITE(!subnormalIncrementedSignificandOverflow,
						  subnormalIncrementedSignificand,
						  unpackedFloat<t>::leadingOne(unpackedFloat<t>::significandWidth(format))));
 


  // Have to choose the right one dependent on rounding mode
  prop choosenGuardBit(ITE(normalRounding, guardBit, subnormalGuardBit));
  prop choosenStickyBit(ITE(normalRounding, stickyBit, subnormalStickyBit));
  
  prop significandEven(ITE(normalRounding,
			   extractedSignificand.extract(0,0).isAllZeros(),
			   ((extractedSignificand & incrementAmount).isAllZeros())));
  prop roundUp(roundingDecision<t>(roundingMode, uf.getSign(), significandEven,
				   choosenGuardBit, choosenStickyBit, prop(false)));
  
  ubv roundedSignificand(ITE(normalRounding,
			     ITE((!roundUp), extractedSignificand, correctedIncrementedSignificand),
			     ITE((!roundUp), maskedSignificand, subnomalCorrectedIncrementedSignificand)));



  /*** Round to correct exponent. ***/

  // The extend is almost certainly unnecessary (see specialised rounders)
  sbv extendedExponent(exp.extend(1));

  prop incrementExponent(ITE(normalRounding,
			     incrementedSignificandOverflow,
			     subnormalIncrementedSignificandOverflow)
			 && roundUp);
  probabilityAnnotation<t>(incrementExponent, VERYUNLIKELY);
  
  sbv correctedExponent(conditionalIncrement<t>(incrementExponent, extendedExponent));

  // This can over and underflow but these values will not be used
  bwt currentExponentWidth(correctedExponent.getWidth());
  sbv roundedExponent(correctedExponent.contract(currentExponentWidth - targetExponentWidth));


  /*** Finish ***/
  prop computedOverflow(correctedExponent > unpackedFloat<t>::maxNormalExponent(format).extend(currentExponentWidth - targetExponentWidth));
  prop computedUnderflow(correctedExponent < unpackedFloat<t>::minSubnormalExponent(format).extend(currentExponentWidth - targetExponentWidth));

  // So that ITE abstraction works...
  prop overflow(ITE(earlyOverflow, prop(true), computedOverflow));
  prop underflow(ITE(earlyUnderflow, prop(true), computedUnderflow));

  unpackedFloat<t> roundedResult(uf.getSign(), roundedExponent, roundedSignificand);
  unpackedFloat<t> result(rounderSpecialCases<t>(format, roundingMode, roundedResult,
						 overflow, underflow, uf.getZero()));

  POSTCONDITION(result.valid(format));

  return result;
 }


template <class t>
  unpackedFloat<t> rounder (const typename t::fpt &format,
			    const typename t::rm &roundingMode,
			    const unpackedFloat<t> &uf) {
  typedef typename t::prop prop;
  customRounderInfo<t> cri(prop(false), prop(false), prop(false), prop(false), prop(false));  // Default is to know nothing

  #ifdef USE_ORIGINAL_ROUNDER
  return originalRounder(format, roundingMode, uf);  // Allow old versions to be compared
  #else
  return customRounder(format, roundingMode, uf, cri);
  #endif
 }


}

#endif
