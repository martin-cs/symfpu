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
** convert.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 03/02/14
**
** Conversion from unpacked floats in one format to another.
**
*/

#include "symfpu/core/unpackedFloat.h"
#include "symfpu/core/rounder.h"

#ifndef SYMFPU_CONVERT
#define SYMFPU_CONVERT

namespace symfpu {

template <class t>
unpackedFloat<t> convertFloatToFloat (const typename t::fpt &sourceFormat,
				      const typename t::fpt &targetFormat,
				      const typename t::rm &roundingMode,
				      const unpackedFloat<t> &input) {

  PRECONDITION(input.valid(sourceFormat));

  typedef typename t::bwt bwt;
  //typedef typename t::prop prop;
  //typedef typename t::ubv ubv;
  //typedef typename t::sbv sbv;

  // increased includes equality
  bool exponentIncreased = unpackedFloat<t>::exponentWidth(sourceFormat) <= unpackedFloat<t>::exponentWidth(targetFormat);
  bool significandIncreased = unpackedFloat<t>::significandWidth(sourceFormat) <= unpackedFloat<t>::significandWidth(targetFormat);

  bwt expExtension = (exponentIncreased) ? unpackedFloat<t>::exponentWidth(targetFormat) - unpackedFloat<t>::exponentWidth(sourceFormat) : 0;
  bwt sigExtension = (significandIncreased) ? unpackedFloat<t>::significandWidth(targetFormat) - unpackedFloat<t>::significandWidth(sourceFormat) : 0;

  unpackedFloat<t> extended(input.extend(expExtension, sigExtension));

  // Format sizes are literal so it is safe to branch on them
  if (exponentIncreased && significandIncreased) {
    // Fast path strict promotions

    POSTCONDITION(extended.valid(targetFormat));
    
    return extended;

  } else {

    unpackedFloat<t> rounded(rounder(targetFormat, roundingMode, extended));

    unpackedFloat<t> result(ITE(input.getNaN(),
				unpackedFloat<t>::makeNaN(targetFormat),
				ITE(input.getInf(),
				    unpackedFloat<t>::makeInf(targetFormat, input.getSign()),
				    ITE(input.getZero(),
					unpackedFloat<t>::makeZero(targetFormat, input.getSign()),
					rounded))));
    
    POSTCONDITION(result.valid(targetFormat));
    
    return result;
  }
}


template <class t>
unpackedFloat<t> roundToIntegral (const typename t::fpt &format,
				  const typename t::rm &roundingMode,
				  const unpackedFloat<t> &input) {

  PRECONDITION(input.valid(format));

  typedef typename t::bwt bwt;
  typedef typename t::prop prop;
  typedef typename t::ubv ubv;
  typedef typename t::sbv sbv;

  sbv exponent(input.getExponent());
  bwt exponentWidth(exponent.getWidth());
  
  sbv packedSigWidth(exponentWidth, format.packedSignificandWidth());
  sbv unpackedSigWidth(exponentWidth, format.significandWidth());
  
  // Fast path for things that must be integral
  prop isIntegral(exponent >= packedSigWidth);
  prop isSpecial(input.getNaN() || input.getInf() || input.getZero());
  prop isID(isIntegral || isSpecial);
  probabilityAnnotation<t>(isID, LIKELY);
  // TODO : fast path the cases that don't round up

  
  // Otherwise, compute rounding location
  sbv initialRoundingPoint(expandingSubtract<t>(packedSigWidth,exponent));  // Expansion only needed in obscure formats
  sbv roundingPoint(collar<t>(initialRoundingPoint,
			      sbv::zero(exponentWidth + 1),
			      unpackedSigWidth.extend(1).increment()));

  // Round
  ubv significand(input.getSignificand());
  significandRounderResult<t> roundedResult(variablePositionRound<t>(roundingMode, input.getSign(), significand,
								     roundingPoint.toUnsigned().matchWidth(significand),
								     prop(false), // TODO : Could actually be exponent >= 0
								     isID));      // The fast-path case so just deactives some code

  // Reconstruct
  // Note this is not in a valid form if significand is all zeros
  // The max is necessary to catch cases when we round up to one from very small numbers
  // The rounder ensures these are zero if they don't round up
  unpackedFloat<t> reconstructed(input.getSign(),
				 max<t>(conditionalIncrement<t>(roundedResult.incrementExponent, exponent),
					sbv::zero(exponentWidth)),
				 roundedResult.significand);
					    
  
  unpackedFloat<t> result(ITE(isID,
			      input,
			      ITE(roundedResult.significand.isAllZeros(),
				  unpackedFloat<t>::makeZero(format, input.getSign()),
				  reconstructed)));

  POSTCONDITION(result.valid(format));
  
  return result;

  
  #if 0
  // increased includes equality
  bool exponentIncreased = unpackedFloat<t>::exponentWidth(sourceFormat) <= unpackedFloat<t>::exponentWidth(targetFormat);
  bool significandIncreased = unpackedFloat<t>::significandWidth(sourceFormat) <= unpackedFloat<t>::significandWidth(targetFormat);

  bwt expExtension = (exponentIncreased) ? unpackedFloat<t>::exponentWidth(targetFormat) - unpackedFloat<t>::exponentWidth(sourceFormat) : 0;
  bwt sigExtension = (significandIncreased) ? unpackedFloat<t>::significandWidth(targetFormat) - unpackedFloat<t>::significandWidth(sourceFormat) : 0;

  unpackedFloat<t> extended(input.extend(expExtension, sigExtension));

  // Format sizes are literal so it is safe to branch on them
  if (exponentIncreased && significandIncreased) {
    // Fast path strict promotions

    POSTCONDITION(extended.valid(targetFormat));
    
    return extended;

  } else {

    unpackedFloat<t> rounded(rounder(targetFormat, roundingMode, extended));

    unpackedFloat<t> result(ITE(input.getNaN(),
				unpackedFloat<t>::makeNaN(targetFormat),
				ITE(input.getInf(),
				    unpackedFloat<t>::makeInf(targetFormat, input.getSign()),
				    ITE(input.getZero(),
					unpackedFloat<t>::makeZero(targetFormat, input.getSign()),
					rounded))));
    
    POSTCONDITION(result.valid(targetFormat));
    
    return result;
  }
  #endif

}


template <class t>
  unpackedFloat<t> convertUBVToFloat (const typename t::fpt &targetFormat,
				      const typename t::rm &roundingMode,
				      const typename t::ubv &input,
				      const typename t::bwt &decimalPointPosition = 0) {
  
  typedef typename t::bwt bwt;
  typedef typename t::prop prop;
  typedef typename t::sbv sbv;
  typedef typename t::fpt fpt;

  bwt inputWidth(input.getWidth());

  PRECONDITION(decimalPointPosition <= inputWidth);
  
  // Devise an appropriate format 
  bwt initialExponentWidth(bitsToRepresent<bwt>(inputWidth) + 1); // +1 as unsigned -> signed
  fpt initialFormat(initialExponentWidth, inputWidth);
  bwt actualExponentWidth(unpackedFloat<t>::exponentWidth(initialFormat));

  // Build
  unpackedFloat<t> initial(prop(false), sbv(actualExponentWidth, (inputWidth - 1) - decimalPointPosition), input);  // inputWidth - 1 as we want one bit above the decimal point
  
  // Normalise
  unpackedFloat<t> normalised(initial.normaliseUpDetectZero(initialFormat));

  // Round (the conversion will catch the cases where no rounding is needed)
  return convertFloatToFloat(initialFormat, targetFormat, roundingMode, normalised);
 }

 
template <class t>
  unpackedFloat<t> convertSBVToFloat (const typename t::fpt &targetFormat,
				      const typename t::rm &roundingMode,
				      const typename t::sbv &input,
				      const typename t::bwt &decimalPointPosition = 0) {
  typedef typename t::bwt bwt;
  typedef typename t::prop prop;
  typedef typename t::sbv sbv;
  typedef typename t::fpt fpt;

  bwt inputWidth(input.getWidth());

  PRECONDITION(decimalPointPosition <= inputWidth);
  
  // Devise an appropriate format 
  bwt initialExponentWidth(bitsToRepresent<bwt>(inputWidth) + 1); // +1 as unsigned -> signed
  fpt initialFormat(initialExponentWidth, inputWidth + 1);        // +1 as signed -> unsigned
  bwt actualExponentWidth(unpackedFloat<t>::exponentWidth(initialFormat));

  // Work out the sign
  prop negative(input < sbv::zero(inputWidth));

  // Build
  unpackedFloat<t> initial(negative, sbv(actualExponentWidth, inputWidth - decimalPointPosition), (abs<t,sbv>(input.extend(1))).toUnsigned());
  
  // Normalise
  unpackedFloat<t> normalised(initial.normaliseUpDetectZero(initialFormat));

  // Round (the conversion will catch the cases where no rounding is needed)
  return convertFloatToFloat(initialFormat, targetFormat, roundingMode, normalised);
 }


 // Common conversion code for both convert to sgined and to unsigned.
 // Note that the results will be junk if it is not in bounds, etc.
 // convertFloatToUBV and convertFloatToSBV handle all of that logic.
 template <class t>
   significandRounderResult<t> convertFloatToBV (const typename t::fpt &format,
						 const typename t::rm &roundingMode,
						 const unpackedFloat<t> &input,
						 const typename t::bwt &targetWidth,
						 const typename t::bwt &decimalPointPosition) {
   
   typedef typename t::bwt bwt;
   typedef typename t::prop prop;
   typedef typename t::ubv ubv;
   typedef typename t::sbv sbv;


   PRECONDITION(decimalPointPosition < targetWidth);


   // TODO : fast path the RTZ / don't need to round case

   bwt maxShift(targetWidth + 1); // + 1 as we have to shift over the guard bit
   bwt maxShiftBits(bitsToRepresent(maxShift) + 1); // +1 as we want it to be signed

   bwt exponentWidth(input.getExponent().getWidth());
   bwt workingExponentWidth((exponentWidth >= maxShiftBits) ?
			    exponentWidth : maxShiftBits);

   sbv maxShiftAmount(workingExponentWidth, maxShift);
   sbv exponent(input.getExponent().matchWidth(maxShiftAmount));


   // Optimisation : compact the significand in the case targetWidth < significantWidth
   ubv inputSignificand(input.getSignificand());
   bwt inputSignificandWidth(inputSignificand.getWidth());
   ubv *working = NULL;
   if (targetWidth + 2 < inputSignificandWidth) {

     ubv dataAndGuard(inputSignificand.extract(inputSignificandWidth - 1, (inputSignificandWidth - targetWidth) - 1));
     prop sticky(!inputSignificand.extract((inputSignificandWidth - targetWidth) - 2, 0).isAllZeros());

     working = new ubv(dataAndGuard.append(ubv(sticky)));
   } else {
     working = new ubv(inputSignificand);
   }
   ubv significand(*working);
   delete working;
   bwt significandWidth(significand.getWidth());

   // Handle zero
   ubv zerodSignificand(significand &
			ITE(input.getZero(), ubv::zero(significandWidth), ubv::allOnes(significandWidth)));
   ubv expandedSignificand(zerodSignificand.extend(maxShift)); // Start with the significand in the sticky position.
                                                               // targetWidth +1 is for the guard bit
   
   // Align
   sbv shiftAmount(collar<t>(expandingAdd<t>(exponent,
					     sbv(workingExponentWidth, decimalPointPosition + 2)),  // +1 to guard, +1 to LSB
			     sbv::zero(workingExponentWidth + 1),
			     maxShiftAmount.extend(1)));
   ubv convertedShiftAmount(shiftAmount.resize(bitsToRepresent(maxShift) + 1 /* +1 for sign bit, safe due to collar */
					       ).toUnsigned().matchWidth(expandedSignificand));
   ubv aligned(expandedSignificand << convertedShiftAmount); // Safe by collar


   // Fixed position round
   significandRounderResult<t> rounded(fixedPositionRound<t>(roundingMode, input.getSign(),
							     aligned, targetWidth,
							     prop(false), prop(false)));

   return rounded;
 }

 // A more compact version for round to zero
 // Only handles normal, subnormal and zero cases, overflow of targetWidth will give junk.
 // Inf, NaN, and overflow must be handled by the caller.
 template <class t>
   significandRounderResult<t> convertFloatToBVRTZ (const typename t::fpt &format,
						    const unpackedFloat<t> &input,
						    const typename t::bwt &targetWidth,
						    const typename t::bwt &decimalPointPosition) {
   typedef typename t::bwt bwt;
   typedef typename t::prop prop;
   typedef typename t::ubv ubv;
   typedef typename t::sbv sbv;

   PRECONDITION(targetWidth > 0);
   PRECONDITION(decimalPointPosition < targetWidth);

   // A working significand of the right length
   ubv significand(input.getSignificand());
   bwt significandWidth(significand.getWidth());

   ubv significantSignificand(significand.extract(significandWidth - 1,
						  ((targetWidth < significandWidth) ? significandWidth - targetWidth : 0)));
   bwt ssWidth(significantSignificand.getWidth());


   // Handle zero and fractional cases
   sbv exponent(input.getExponent());
   bwt exponentWidth(exponent.getWidth());

   prop fraction(input.getExponent() < sbv::zero(exponentWidth));
   ubv zerodSignificand(significantSignificand &
			ITE(input.getZero() || fraction, ubv::zero(ssWidth), ubv::allOnes(ssWidth)));

   ubv expandedSignificand(zerodSignificand.extend(targetWidth - 1)); // Start with the significand in the LSB of output


   // Prepare exponent
   bwt maxShift(targetWidth - 1); // - 1 as we are already at LSB
   bwt maxShiftBits(bitsToRepresent(maxShift)); // Don't care about it being signed

   ubv convertedExponent(exponent.toUnsigned());
   bwt topExtractedBit(((maxShiftBits >  (exponentWidth - 1)) ? exponentWidth : maxShiftBits) - 1);

   ubv shiftBits(convertedExponent.extract(topExtractedBit, 0));
   ubv shiftOperand(shiftBits.matchWidth(expandedSignificand));

   // Align
   ubv shifted(expandedSignificand.modularLeftShift(shiftOperand));
   bwt shiftedWidth(shifted.getWidth());

   // Extract
   ubv result(shifted.extract(shiftedWidth - 1, shiftedWidth - targetWidth));

   return significandRounderResult<t>(result, prop(false));
 }



 // Decimal point position in the bit in the output on the left hand side of the decimal point
 // I.E. if it is positive then it is converting to a fix-point number
 template <class t>
   typename t::ubv convertFloatToUBV (const typename t::fpt &format,
				      const typename t::rm &roundingMode,
				      const unpackedFloat<t> &input,
				      const typename t::bwt &targetWidth,
				      const typename t::ubv &undefValue,
				      const typename t::bwt &decimalPointPosition = 0) {

   typedef typename t::bwt bwt;
   typedef typename t::prop prop;
   typedef typename t::ubv ubv;
   typedef typename t::sbv sbv;


   PRECONDITION(decimalPointPosition < targetWidth);


   // Invalid cases
   prop specialValue(input.getInf() || input.getNaN());

   bwt maxExponentValue(targetWidth);
   bwt maxExponentBits(bitsToRepresent(maxExponentValue) + 1);

   bwt exponentWidth(input.getExponent().getWidth());
   bwt workingExponentWidth((exponentWidth >= maxExponentBits) ?
			    exponentWidth : maxExponentBits);

   sbv maxExponent(workingExponentWidth, maxExponentValue);
   sbv exponent(input.getExponent().matchWidth(maxExponent));

   prop tooLarge(exponent >= maxExponent);

   prop tooNegative(input.getSign() &&
		    !input.getZero() &&  // Zero is handled elsewhere
		    sbv::zero(workingExponentWidth) <= exponent);  // Can't round to 0
   
   prop earlyUndefinedResult(specialValue || tooLarge || tooNegative);
   probabilityAnnotation<t>(earlyUndefinedResult, LIKELY); // Convertable values are rare


   // Fixed position round
   significandRounderResult<t> rounded(convertFloatToBV(format, roundingMode, input,
							targetWidth, decimalPointPosition));

   // TODO : fast path negative by converting exp==0 into guard and exp < 0 into sticky
   
   // Put the result together
   prop undefinedResult(earlyUndefinedResult ||
			rounded.incrementExponent ||    // Overflow
			(input.getSign() && !rounded.significand.isAllZeros()));  // Negative case
   
   ubv result(ITE(undefinedResult,
		  undefValue,
		  rounded.significand));

   return result;
 }

  // Decimal point position in the bit in the output on the left hand side of the decimal point
  // I.E. if it is positive then it is converting to a fix-point number
  template <class t>
    typename t::sbv convertFloatToSBV (const typename t::fpt &format,
				       const typename t::rm &roundingMode,
				       const unpackedFloat<t> &input,
				       const typename t::bwt &targetWidth,
				       const typename t::sbv &undefValue,
				       const typename t::bwt &decimalPointPosition = 0) {

   typedef typename t::bwt bwt;
   typedef typename t::prop prop;
   //typedef typename t::ubv ubv;
   typedef typename t::sbv sbv;


   PRECONDITION(decimalPointPosition < targetWidth);


   // Invalid cases
   prop specialValue(input.getInf() || input.getNaN());

   bwt maxExponentValue(targetWidth);
   bwt maxExponentBits(bitsToRepresent(maxExponentValue) + 1);

   bwt exponentWidth(input.getExponent().getWidth());
   bwt workingExponentWidth((exponentWidth >= maxExponentBits) ?
			    exponentWidth : maxExponentBits);

   sbv maxExponent(workingExponentWidth, maxExponentValue);
   sbv exponent(input.getExponent().matchWidth(maxExponent));

   prop tooLarge(exponent >= maxExponent);

   prop earlyUndefinedResult(specialValue || tooLarge);
   probabilityAnnotation<t>(earlyUndefinedResult, LIKELY); // Convertable values are rare


   // Fixed position round
   // (It is tempting to think that this could be done with targetWidth - 1 bits but that
   // missed the case of things like -128.05 -> int8_t)
   significandRounderResult<t> rounded(convertFloatToBV(format, roundingMode, input,
							targetWidth, decimalPointPosition));

   // Put the result together
   bwt roundSigWidth(rounded.significand.getWidth());
   prop undefinedResult(earlyUndefinedResult ||
			rounded.incrementExponent ||    // Definite Overflow
			(rounded.significand.extract(roundSigWidth - 1,
						     roundSigWidth - 1).isAllOnes() &&
			 !(input.getSign() && rounded.significand.extract(roundSigWidth - 2, 0).isAllZeros()))); // -2^{n-1} is the only safe "overflow" case

   
   sbv result(ITE(undefinedResult,
		  undefValue,
		  conditionalNegate<t,sbv,prop>(input.getSign(), rounded.significand.toSigned())));

   return result;
 }

}

#endif
