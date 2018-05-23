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
** operations.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 15/11/14
**
** A number of compound operations on bit-vectors.  These are default
** implementations to reduce the amount of code that is needed for a
** back-end (and the risk of making a mistake).  Back-ends can provide
** their own implementations of these if they can handle them in a
** smarter way than the default.
**
*/

#include <cassert>
#include <map>

#include "../utils/common.h"


#ifndef SYMFPU_OPERATIONS
#define SYMFPU_OPERATIONS

namespace symfpu {

  /*** Expanding operations ***/
  template <class t, class bv>
  bv expandingAdd (const bv &op1, const bv &op2) {
    PRECONDITION(op1.getWidth() == op2.getWidth());

    bv x(op1.extend(1));
    bv y(op2.extend(1));
    
    return x + y;
  }

  template <class t, class bv, class prop>
    bv expandingAddWithCarryIn (const bv &op1, const bv &op2, const prop &cin) {
    PRECONDITION(op1.getWidth() == op2.getWidth());

    bv x(op1.extend(1));
    bv y(op2.extend(1));
    
    bv sum(x + y);

    typename t::bwt w(sum.getWidth());
    bv carry(ITE(cin, bv::one(w), bv::zero(w)));
    bv res(sum.modularAdd(carry)); // Modular is safe due to the extension
                                   // (2^n - 1) + (2^n - 1) + 1 == 2^(n+1) - 1
                                   // -(2^n) + -(2^n) + 1 > 2^(n+1)
    return res;
  }

  template <class t, class bv>
  bv expandingSubtract (const bv &op1, const bv &op2) {
    PRECONDITION(op1.getWidth() == op2.getWidth());

    bv x(op1.extend(1));
    bv y(op2.extend(1));
    
    return x - y;
  }

  template <class t, class bv>
  bv expandingMultiply (const bv &op1, const bv &op2) {
    typename t::bwt width = op1.getWidth();
    PRECONDITION(width == op2.getWidth());
    
    bv x(op1.extend(width));
    bv y(op2.extend(width));

    return x * y;
  }


  /*** Conditional Operations ***/
  template <class t, class bv, class prop>
  bv conditionalIncrement (const prop &p, const bv &b) {
    PRECONDITION(IMPLIES(p, b <  bv::maxValue(b.getWidth())));

    typename t::bwt w(b.getWidth());
    bv inc(ITE(p, bv::one(w), bv::zero(w)));
    
    return b + inc;
  }

  template <class t, class bv, class prop>
  bv conditionalDecrement (const prop &p, const bv &b) {
    PRECONDITION(IMPLIES(p, bv::minValue(b.getWidth()) < b));

    typename t::bwt w(b.getWidth());
    bv inc(ITE(p, bv::one(w), bv::zero(w)));

    return b - inc;
  }

  template <class t, class bv, class prop>
  bv conditionalLeftShiftOne (const prop &p, const bv &b) {
    typename t::bwt w(b.getWidth());
    PRECONDITION(IMPLIES(p, (b.extract(w - 1, w - 1).isAllZeros())));

    bv shifted(b.modularLeftShift(bv::one(w)));
    return bv(ITE(p, shifted, b));
  }

  template <class t, class bv, class prop>
  bv conditionalRightShiftOne (const prop &p, const bv &b) {
    typename t::bwt w(b.getWidth());
    // PRECONDITION(IMPLIES(p, (b.extract(0, 0).isAllZeros())));  // Adder uses and compensates for this case.

    bv shifted(b.modularRightShift(bv::one(w)));
    return bv(ITE(p, shifted, b));
  }

  template <class t, class bv, class prop>
  bv conditionalNegate (const prop &p, const bv &b) {
    typename t::bwt w(b.getWidth());
    PRECONDITION(w >= 2);
    PRECONDITION(IMPLIES(p, !(b.extract(w - 1, w - 1).isAllOnes() &&
			      b.extract(w - 2,     0).isAllZeros())));
    
    return bv(ITE(p, -b, b));
  }

  template <class t, class bv>
  bv abs(const bv &b) {
    return conditionalNegate<t, bv, typename t::prop>(b < bv::zero(b.getWidth()), b);
  }

  

  /*** Probability Annotations ***/
  enum probability {
    VERYLIKELY = 100,
    LIKELY = 50,
    NEUTRAL = 0,
    UNLIKELY = -50,
    VERYUNLIKELY = -100
  };

  template <class t, class prop>
  void probabilityAnnotation (const prop &, const probability &) {
    // Back-ends can make use of this information if they want
    return;
  }


  /*** Max and min ***/
  template <class t, class bv>
  bv max (const bv &op1, const bv &op2) {
    return ITE(op1 <= op2, op2, op1);
  }

  template <class t, class bv>
  bv min (const bv &op1, const bv &op2) {
    return ITE(op1 <= op2, op1, op2);
  }

  template <class t, class bv>
  bv collar(const bv &op, const bv &lower, const bv &upper) {
    return ITE(op < lower,
	       lower,
	       ITE(upper < op,
		   upper,
		   op));
  }
  

  /*** Unary/Binary operations ***/
  template <class t, class bv, class prop, class bwt>
  bv countLeadingZerosRec (const bv &op, const bwt position, const prop &allPreceedingZeros) {
    typename t::bwt w(op.getWidth());
    
    PRECONDITION(0 <= position && position < w);
    
    bv bit(op.extract(position, position));
    
    prop isLeadingOne(allPreceedingZeros && (bit.isAllOnes()));
    prop continuingZero(allPreceedingZeros && (bit.isAllZeros()));

    if (position == 0) {
      return ITE(isLeadingOne, bv(w, w - 1), bv(w, w));
    } else {
      return ITE(isLeadingOne,
		 bv(w, w - (position + 1)),
		 countLeadingZerosRec<t>(op, position - 1, continuingZero));
    }
  }
  
  template <class t, class bv>
  bv countLeadingZeros (const bv &op) {
    typedef typename t::bwt bwt;
    typedef typename t::prop prop;    
    bwt w(op.getWidth());

    return countLeadingZerosRec<t>(op, w - 1, prop(true));
  }
  
  // This is sort of the opposite of count trailing 1's (a.k.a. clz(reverse(not(x))) )
  template <class t, class bv>
  bv orderEncode (const bv &op) {
    typename t::bwt w(op.getWidth());
    
    //PRECONDITION(bv::zero(w) <= op && op <= bv(w, w)); // Not needed as using modular shift

    bv tmp((bv::one(w + 1).modularLeftShift(op.resize(w + 1))).modularDecrement().extract(w-1,0));
    return tmp;
  }


  // The comparison we need to do is
  //   op.extract(relevantBits - 1, 0) == bv(relevantBits, position + 1)
  // bits can be shared between instances of this.
  // The following is a dynamic programming style solution.
  // HOWEVER : it may not actually give a net saving depending on your application
  template <class t, class bv>
  struct fragmentMap {
    typedef typename t::bwt bwt;
    typedef typename t::prop prop;
    
    typedef std::pair<bwt, bwt> fragment;
    typedef std::map<fragment, prop> map;

  protected :
    const bv &op;
    map m;

    prop getComparitorRec(bwt length, bwt value)
    {
      PRECONDITION(length > 0);
      PRECONDITION(bitsToRepresent(value) <= length);
      typename map::const_iterator it = m.find(std::make_pair(length, value));
      
      if (it != m.end()) {
	return it->second;
      } else {
	bwt leadingBit = bwt(1) << (length - 1);
	prop leadingBitIsOne(op.extract(length - 1, length - 1).isAllOnes());
	prop correctComparison((value & leadingBit) ?
			       leadingBitIsOne : !leadingBitIsOne);
	prop *step = NULL;
	
	if (length == 1) {
	  step = new prop(correctComparison);
	} else {
	  prop rec(getComparitorRec(length - 1, value & (~leadingBit)));
	  step = new prop(correctComparison && rec);
	}

	prop res(*step);
	delete step;
	
	m.insert(std::make_pair(std::make_pair(length, value), res));
	return res;
      }
    }
    
  public :
    fragmentMap(const bv &_op) : op(_op) {}

    prop getComparitor(bwt length, bwt value)
    {
      PRECONDITION(length > 0);
      PRECONDITION(bitsToRepresent(value) <= length);

      prop res(getComparitorRec(length, value));
      
      POSTCONDITION(bv(res) == (op.extract(length - 1, 0) == bv(length, value)));
      return res;
    }    
  };
 
  
  // A more compact, bitwise implementation of orderEncode for SAT encoding
  // Intended to be used in specialisation of orderEncode
  template <class t, class bv>
  bv orderEncodeBitwise (const bv &op) {
    typedef typename t::bwt bwt;
    bwt w(op.getWidth());

    fragmentMap<t,bv> m(op);

    // If op is too large, then set everything to 1
    bv outOfRange(op >= bv(w, w));
    
    // SAND to fill in the remaining bits
    bv * working = new bv(outOfRange);
    for (bwt i = w; i > 0; --i) {
      bwt position = i - 1;    // Position in the output bitvectors
      bwt relevantBits = bitsToRepresent(position + 1);
      INVARIANT(relevantBits > 0);
      
      //bv activateBit(m.getComparitor(relevantBits, position + 1));  // No more compact and slower
      bv activateBit(op.extract(relevantBits - 1, 0) == bv(relevantBits, position + 1));
      bv nextBit(working->extract(0,0) | activateBit);
      
      bv * tmp = working;
      working = new bv(working->append(nextBit));
      delete tmp;
    }

    bv output(working->extract(w - 1,0));
    delete working;

    POSTCONDITION(output == (bv::one(w + 1).modularLeftShift(op.resize(w + 1))).modularDecrement().extract(w-1,0));
    
    return output;
  }

  
  /*** Custom shifts ***/
  // 1 if and only if the right shift moves at least one 1 out of the word
  template <class t, class bv>
  bv rightShiftStickyBit (const bv &op, const bv &shift) {
    bv stickyBit(ITE((orderEncode<t>(shift) & op).isAllZeros(),
		     bv::zero(op.getWidth()),
		     bv::one(op.getWidth())));
    
    return stickyBit;
  }


  // It is easier to compute along with the shift
  template <class t>
  struct stickyRightShiftResult {
    typedef typename t::ubv ubv;

    ubv signExtendedResult;
    ubv stickyBit;

    stickyRightShiftResult(const ubv &ser, const ubv &sb) : signExtendedResult(ser), stickyBit(sb) {}
  stickyRightShiftResult(const stickyRightShiftResult &old) : signExtendedResult(old.signExtendedResult), stickyBit(old.stickyBit) {}
  };

  
  template <class t>
  stickyRightShiftResult<t> stickyRightShift (const typename t::ubv &input, const typename t::ubv &shiftAmount) {    
    stickyRightShiftResult<t> res(input.signExtendRightShift(shiftAmount), rightShiftStickyBit<t>(input, shiftAmount));

    return res;
  }

  
  template <class t>
  stickyRightShiftResult<t> stickyRightShiftBitwise (const typename t::ubv &input, const typename t::ubv &shiftAmount) {
    typedef typename t::bwt bwt;
    typedef typename t::prop prop;
    typedef typename t::ubv ubv;

    bwt width(input.getWidth());
    bwt startingPosition(positionOfLeadingOne(width));
    INVARIANT(0 < startingPosition && startingPosition < width);
    
    // Catch the out of bounds case
    PRECONDITION(shiftAmount.getWidth() == width);
    prop fullShift(shiftAmount >= ubv(width, width));
    // Note the shiftAmount is treated as unsigned...

    ubv *working = new ubv(input);
    prop *stickyBit = new prop(ITE(fullShift, !input.isAllZeros(), prop(false)));
    
    for (bwt i = startingPosition + 1; i > 0; --i)
    {
      bwt shiftAmountPosition = i - 1;

      prop shiftEnabled = fullShift || shiftAmount.extract(shiftAmountPosition, shiftAmountPosition).isAllOnes();

      prop stickyAccumulate(shiftEnabled && !(working->extract((1ULL << shiftAmountPosition) - 1, 0).isAllZeros())); // Would this loose data?
      
      prop * tmp = stickyBit;
      stickyBit = new prop(*stickyBit || stickyAccumulate);
      delete tmp;


      // Note the slightly unexpected sign extension
      ubv shifted(working->signExtendRightShift(ubv::one(width) << ubv(width, shiftAmountPosition)));
      
      ubv * workingTmp = working;
      working = new ubv(ITE(shiftEnabled, shifted, *working));
      delete workingTmp;
    }
    
    stickyRightShiftResult<t> res(*working, ubv(*stickyBit).extend(width - 1));

    delete working;
    delete stickyBit;

    POSTCONDITION(res.signExtendedResult == input.signExtendRightShift(shiftAmount));
    POSTCONDITION(res.stickyBit == rightShiftStickyBit<t>(input, shiftAmount));

    return res;
  }
  
  
  
  template <class t>
  struct normaliseShiftResult {
    typedef typename t::ubv ubv;
    typedef typename t::prop prop;
    
    ubv normalised;
    ubv shiftAmount;
    prop isZero;
    
  normaliseShiftResult(const ubv &n, const ubv &s, const prop &z) : normalised(n), shiftAmount(s), isZero(z) {}
  normaliseShiftResult(const normaliseShiftResult<t> &old) :  normalised(old.normalised), shiftAmount(old.shiftAmount), isZero(old.isZero) {}
  };

  template <class t>
  normaliseShiftResult<t> normaliseShift (const typename t::ubv input) {
    typedef typename t::bwt bwt;
    typedef typename t::prop prop;
    typedef typename t::ubv ubv;

    bwt width(input.getWidth());
    bwt startingMask(previousPowerOfTwo(width));
    INVARIANT(startingMask < width);
    
    // Catch the zero case
    prop zeroCase(input.isAllZeros());

    ubv *working = new ubv(input);
    ubv *shiftAmount = NULL;
    prop *deactivateShifts = new prop(zeroCase);
    
    for (bwt i = startingMask; i > 0; i >>= 1) {
      prop newDeactivateShifts = *deactivateShifts || working->extract(width-1,width-1).isAllOnes();
      delete deactivateShifts;
      deactivateShifts = new prop(newDeactivateShifts);
      
      ubv mask(ubv::allOnes(i).append(ubv::zero(width - i)));
      prop shiftNeeded(!(*deactivateShifts) && (mask & *working).isAllZeros());

      // Modular is safe because of the mask comparison
      ubv shifted(ITE(shiftNeeded, working->modularLeftShift(ubv(width, i)), *working));
      delete working;
      working = new ubv(shifted);

      if (shiftAmount == NULL) {
	shiftAmount = new ubv(shiftNeeded);
      } else {
	ubv newShiftAmount = shiftAmount->append(ubv(shiftNeeded));
	delete shiftAmount;
	shiftAmount = new ubv(newShiftAmount);
      }
    }

    normaliseShiftResult<t> res(*working, *shiftAmount, zeroCase);

    delete deactivateShifts;
    delete working;
    delete shiftAmount;

    POSTCONDITION(res.normalised.extract(width-1,width-1).isAllZeros() == res.isZero);
    POSTCONDITION(IMPLIES(res.isZero, res.shiftAmount.isAllZeros()));

    bwt shiftAmountWidth(res.shiftAmount.getWidth());
    bwt widthBits(bitsToRepresent(width));
    POSTCONDITION(shiftAmountWidth == widthBits ||
		  shiftAmountWidth == widthBits - 1); // If width is an exact power of 2
    ubv widthBV(widthBits, width);
    POSTCONDITION(res.shiftAmount.matchWidth(widthBV) < widthBV);

    return res;
  }


  /*** Dividers ***/
  template <class t>
  struct resultWithRemainderBit {
    typedef typename t::ubv ubv;
    typedef typename t::prop prop;
    
    ubv result;
    prop remainderBit;
    
  resultWithRemainderBit(const ubv &o, const prop &r) : result(o), remainderBit(r) {}
  resultWithRemainderBit(const resultWithRemainderBit<t> &old) : result(old.result), remainderBit(old.remainderBit) {}
  };
  
  // x and y are fixed-point numbers in the range [1,2)
  // Compute o \in [0.5,2), r \in [0,\delta) such that:  x = o*y + r
  // Return (o, r != 0)
  template <class t>
  resultWithRemainderBit<t> fixedPointDivide (const typename t::ubv &x, const typename t::ubv &y) {
    typename t::bwt w(x.getWidth());

    // Same width and both have MSB ones
    PRECONDITION(y.getWidth() == w);
    PRECONDITION(x.extract(w - 1, w - 1).isAllOnes());
    PRECONDITION(y.extract(w - 1, w - 1).isAllOnes());

    typedef typename t::ubv ubv;
    
    // Not the best way of doing this but pretty universal
    ubv ex(x.append(ubv::zero(w - 1)));
    ubv ey(y.extend(w - 1));

    ubv div(ex / ey);
    ubv rem(ex % ey);

    return resultWithRemainderBit<t>(div.extract(w - 1, 0), !(rem.isAllZeros()));
  }

  
  // x is a fixed-point number in the range [1,4) with 2/p bits
  // Compute o \in [1,sqrt(2)), r \in [0,o*2 + 1) such that x = o*o + r with 1/p bits
  // Return (o, r != 0)
  template <class t>
  resultWithRemainderBit<t> fixedPointSqrt (const typename t::ubv &x) {
    typedef typename t::bwt bwt;
    typedef typename t::ubv ubv;
    typedef typename t::prop prop;

    // The default algorithm given here isn't a great one
    // However it is simple and has a very simple termination criteria.
    // Plus most symbolic back ends will prefer
    // o = nondet(), r = nondet(), assert(r < o*2 + 1), assert(x = o*o + r)

    bwt inputWidth(x.getWidth());
    bwt outputWidth(inputWidth - 1);

    // To compare against, we need to pad x to 2/2p
    ubv xcomp(x.append(ubv::zero(inputWidth - 2))); 

    // Start at 1
    ubv working(ubv::one(outputWidth) << ubv(outputWidth, outputWidth - 1));

    bwt location;
    for (location = outputWidth - 1; location > 0; --location) { // Offset by 1 for easy termination
      ubv shift(ubv(outputWidth, location - 1));
      
      ubv candidate(working | (ubv::one(outputWidth) << shift));

      prop addBit(expandingMultiply<t, ubv>(candidate, candidate) <= xcomp);

      working = working | (ubv(addBit).extend(outputWidth - 1) << shift);
    }
    
    return resultWithRemainderBit<t>(working, !(expandingMultiply<t, ubv>(working, working) == xcomp));
  }

  // One step of a divider
  // Here the "remainder bit" is actual the result bit and
  // The result is the remainder
  template <class t>
  resultWithRemainderBit<t> divideStep (const typename t::ubv &x, const typename t::ubv &y) {
    typedef typename t::bwt bwt;
    typedef typename t::ubv ubv;
    typedef typename t::prop prop;

    bwt xWidth(x.getWidth());
    bwt yWidth(y.getWidth());

    PRECONDITION(xWidth == yWidth);
    PRECONDITION(yWidth >= 2);
    PRECONDITION(y.extract(yWidth - 2, yWidth - 2).isAllOnes());  // Assume y is aligned
    
    prop canSubtract(x >= y);
    ubv sub(x.modularAdd(y.modularNegate())); // TODO : modular subtract or better
    ubv step(ITE(canSubtract, sub, x));

    return resultWithRemainderBit<t>(step << ubv::one(xWidth), canSubtract);
  }
}

#endif
