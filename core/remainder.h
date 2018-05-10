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
** remainder.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 14/12/16
**
** Computing the IEEE-754 remainder of arbitrary precision floats
**
*/

#include "symfpu/core/unpackedFloat.h"
#include "symfpu/core/ite.h"
#include "symfpu/core/rounder.h"
#include "symfpu/core/operations.h"
#include "symfpu/core/add.h"
#include "symfpu/core/sign.h"


#ifndef SYMFPU_REMAINDER
#define SYMFPU_REMAINDER

namespace symfpu {

template <class t>
  unpackedFloat<t> addRemainderSpecialCases (const typename t::fpt &format,
					    const unpackedFloat<t> &left,
					    const unpackedFloat<t> &right,
					    const unpackedFloat<t> &remainderResult) {
  typedef typename t::prop prop;

  prop eitherArgumentNan(left.getNaN() || right.getNaN());
  prop generateNan(left.getInf() || right.getZero());
  prop isNan(eitherArgumentNan || generateNan);

  prop passThrough((!(left.getInf() || left.getNaN()) && right.getInf()) ||
		   left.getZero());

  return ITE(isNan,
	     unpackedFloat<t>::makeNaN(format),
	     ITE(passThrough,
		 left,
		 remainderResult));
 }

/* Let left = x*2^e, right = y*2^f, x \in [1,2), y \in [1,2)
 * x/y \in (0.5,2)   x > y  x/y \in (1,2)   x < y (0.5,1)
 *
 *  rem =  x*2^e     - (y*2^f * int((x*2^e) / (y*2^f)))
 *      =  x*2^e     - (y*2^f * int((x/y) * 2^{e-f}))
 *      = (x*2^{e-f} - (y     * int((x/y) * 2^{e-f}))) * 2^f
 *
 *
 * If e - f >  0
 *      = (x*2^{e-f} - (y     * int((x*2^{e-f})/y)) * 2^f
 *
 *
 * If e - f == 0
 *      = (x         - (y     * int((x/y)          ))) * 2^f
 *      = ITTE(x ?= y,
 *             (x - (y * int[guard=1,sticky=1])) * 2^f
 *             (x -  y) * 2^f,
 *             ...)
 *      = ITTE(x ?= y,
 *             (x - (y * int[guard=1,sticky=1])) * 2^f
 *             left - right,
 *             ...)
 *
 *
 * If e - f == -1
 *      = (x*2^{-1}  - (y     * int((x/y) * 2^{-1 }))) * 2^f
 *      = ITTE(x ?= y,
 *             (x*2^{-1}  - (y * int[guard=0,sticky=1])) * 2^f
 *             (x*2^{-1}  - (y * int[guard=1,sticky=0])) * 2^f
 *             (x*2^{-1}  - (y * int[guard=1,sticky=1])) * 2^f
 *
 * If e - f <= -2
 *      = (x*2^{e-f}  - (y     * int[guard=0,sticky=1])) * 2^f
 *      = ITE(int[guard=0,sticky=1],
 *            x*2^e  - y*2^f,
 *            left)
 *      = ITE(int[guard=0,sticky=1],
 *            left  - right,
 *            left)
 *
 */

// Divide for max(e - f, 0) cycles
// The equal case, if you divide you use to extract the even bit of n, also save the rem.
// Then one more cycle for the guard bit
// Use that remainder to work out the sticky bit
// Round and either subtract or not from saved rem
// Output at 2^f
 
template <class t>
  unpackedFloat<t> arithmeticRemainder (const typename t::fpt &format,
					const typename t::rm &roundingMode,
					const unpackedFloat<t> &left,
					const unpackedFloat<t> &right) {
  typedef typename t::bwt bwt;
  typedef typename t::prop prop;
  typedef typename t::ubv ubv;
  typedef typename t::sbv sbv;
  //typedef typename t::fpt fpt;

  PRECONDITION(left.valid(format));
  PRECONDITION(right.valid(format));

  // Compute sign
  prop remainderSign(left.getSign());


  // Compute exponent difference
  sbv exponentDifference(expandingSubtract<t>(left.getExponent(), right.getExponent()));
  bwt edWidth(exponentDifference.getWidth());
  
  // Extend for divide steps
  ubv lsig(left.getSignificand().extend(1));
  ubv rsig(right.getSignificand().extend(1));

  
  ubv first(divideStep<t>(lsig,rsig).result);
  ubv *running = new ubv(first); // To avoid running out of stack space loop with a pointer

  bwt maxDifference = unpackedFloat<t>::maximumExponentDifference(format);
  for (bwt i = maxDifference - 1; i > 0; i--) {
    prop needPrevious(exponentDifference > sbv(edWidth, i));
    probabilityAnnotation<t>(needPrevious, (i > (maxDifference / 2)) ? VERYUNLIKELY : UNLIKELY);
    
    ubv r(ITE(needPrevious, *running, lsig));
    delete running;  // We assume the value / reference has been transfered to ITE
    running = new ubv(divideStep<t>(r, rsig).result);
  }

  // The zero exponent difference case is a little different
  // as we need the result bit for the even flag
  // and the actual result for the final
  prop lsbRoundActive(exponentDifference > -sbv::one(edWidth));  // i.e. >= 0
  
  prop needPrevious(exponentDifference > sbv::zero(edWidth));
  probabilityAnnotation<t>(needPrevious, UNLIKELY);
    
  ubv r0(ITE(needPrevious, *running, lsig));
  delete running;
  resultWithRemainderBit<t> dsr(divideStep<t>(r0, rsig));

  prop integerEven(!lsbRoundActive || !dsr.remainderBit);  // Note negation of guardBit

  
  // The same to get the guard flag
  prop guardRoundActive(exponentDifference > -sbv(edWidth,2));  // i.e. >= -1

  ubv rm1(ITE(lsbRoundActive, dsr.result, lsig));
  resultWithRemainderBit<t> dsrg(divideStep<t>(rm1, rsig));

  prop guardBit(guardRoundActive && dsrg.remainderBit);
  
  prop stickyBit(!ITE(guardRoundActive,
		      dsrg.result,
		      lsig).isAllZeros());
		 

  // The base result if lsbRoundActive
  unpackedFloat<t> reconstruct(remainderSign,
			       right.getExponent(),
			       dsr.result.extract(lsig.getWidth() - 1,1)); // dsr shifts right as last action so is safe

  
  probabilityAnnotation<t>(needPrevious, UNLIKELY);   // Perhaps stretching it a bit but good for approximation
  unpackedFloat<t> candidateResult(ITE(lsbRoundActive,
				       reconstruct.normaliseUpDetectZero(format),
				       left));

  // The final subtract is a little different as previous ones were
  // guaranteed to be positive
  // TODO : This could be improved as these don't need special cases, etc.

  // From the rounding of the big integer multiple
  prop bonusSubtract(roundingDecision<t>(roundingMode,
					 remainderSign,
					 integerEven,
					 guardBit,
					 stickyBit,
					 prop(false)));
  probabilityAnnotation<t>(bonusSubtract, UNLIKELY); // Again, more like 50/50

  // The big integer has sign left.getSign() ^ right.getSign() so we subtract something of left.getSign().
  // For the integer part we handle this by working with absolutes (ignoring the sign) and
  // adding it back in at the end.
  // However for the correction for the rounded part we need to take it into account
  unpackedFloat<t> signCorrectedRight(right, left.getSign());
  unpackedFloat<t> remainderResult(ITE(bonusSubtract,
				       add<t>(format,
					      roundingMode,
					      candidateResult,
					      signCorrectedRight,
					      false),
				       candidateResult));
  
  // TODO : fast path if first.isAllZeros()
  
  POSTCONDITION(remainderResult.valid(format));

  return remainderResult;
 }


// Put it all together...
template <class t>
  unpackedFloat<t> remainderWithRounding (const typename t::fpt &format,
					  const typename t::rm &roundingMode,
					  const unpackedFloat<t> &left,
					  const unpackedFloat<t> &right) {
  //typedef typename t::bwt bwt;
  //typedef typename t::prop prop;
  //typedef typename t::ubv ubv;
  //typedef typename t::sbv sbv;

  PRECONDITION(left.valid(format));
  PRECONDITION(right.valid(format));

  unpackedFloat<t> remainderResult(arithmeticRemainder(format, roundingMode, left, right));
  
  //unpackedFloat<t> result(addRemainderSpecialCases(format, left, right, roundedRemainderResult));
  unpackedFloat<t> result(addRemainderSpecialCases(format, left, right, remainderResult));

  POSTCONDITION(result.valid(format));

  return result;
 }

// IEEE-754 remainder always uses round to nearest, ties to even
template <class t>
  unpackedFloat<t> remainder (const typename t::fpt &format,
			      const unpackedFloat<t> &left,
			      const unpackedFloat<t> &right) {

  return remainderWithRounding<t>(format, t::RNE(), left, right);
 }


}

#endif

