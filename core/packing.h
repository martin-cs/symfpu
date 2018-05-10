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
** packing.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 03/06/14
**
** Algorithms for converting from a packed float (a bit vector) into
** the working, unpacked float format.  Should be instantiated with a
** traits structure from one of the baseTypes/ implementations.
**
*/

#include "symfpu/core/unpackedFloat.h"

#ifndef SYMFPU_PACKING
#define SYMFPU_PACKING

namespace symfpu {

  template<class t>
    unpackedFloat<t> unpack (const typename t::fpt format, const typename t::ubv &packedFloat) {
    typedef typename t::bwt bwt;
    typedef typename t::prop prop;
    typedef typename t::ubv ubv;
    typedef typename t::sbv sbv;

    bwt pWidth = format.packedWidth();
    bwt exWidth = format.packedExponentWidth();
    bwt sigWidth = format.packedSignificandWidth();

    PRECONDITION(packedFloat.getWidth() == pWidth);

    // Extract
    ubv packedSignificand(packedFloat.extract(sigWidth - 1, 0));
    ubv packedExponent(packedFloat.extract(sigWidth + exWidth - 1, sigWidth));
    prop sign(packedFloat.extract(pWidth - 1, sigWidth + exWidth).isAllOnes());

    // Prepare the normal and subnormal cases
    bwt unpackedExWidth = unpackedFloat<t>::exponentWidth(format);
    bwt unpackedSigWidth = unpackedFloat<t>::significandWidth(format);

    INVARIANT(unpackedExWidth > exWidth); // To avoid overflow
    sbv biasedExponent(packedExponent.extend(unpackedExWidth - exWidth).toSigned() - unpackedFloat<t>::bias(format));
    // Optimisation : both the normal and subnormal paths subtract a constant
    // from the exponent as the last step (bias and minNormalExponent respectively)

    ubv significandWithLeadingZero(packedSignificand.extend(unpackedSigWidth - sigWidth));
    ubv significandWithLeadingOne(unpackedFloat<t>::leadingOne(unpackedFloat<t>::significandWidth(format)) | significandWithLeadingZero);

    unpackedFloat<t> ufNormal(sign, biasedExponent, significandWithLeadingOne);
    unpackedFloat<t> ufSubnormalBase(sign, unpackedFloat<t>::minNormalExponent(format), significandWithLeadingZero);

    // Optimisation : right shift the significand by one if subnormal
    //                and then set the carry in when you add to the exponent.
    //                May be sufficient to just assert that it has a leading zero

    
    // Analyse
    prop zeroExponent(packedExponent.isAllZeros());
    prop onesExponent(packedExponent.isAllOnes());
    prop zeroSignificand(significandWithLeadingZero.isAllZeros()); // Shared with normaliseUp

    // Identify the cases
    prop isZero(zeroExponent && zeroSignificand);
    prop isSubnormal(zeroExponent && !zeroSignificand);
    prop isNormal(!zeroExponent && !onesExponent);
    prop isInf(onesExponent && zeroSignificand);
    prop isNaN(onesExponent && !zeroSignificand);

    INVARIANT(isZero || isSubnormal || isNormal || isInf || isNaN);

    probabilityAnnotation<t,prop>(isSubnormal, UNLIKELY);
    
    // Splice together
    unpackedFloat<t> uf(ITE(isNaN,
			    unpackedFloat<t>::makeNaN(format),
			    ITE(isInf,
				unpackedFloat<t>::makeInf(format, sign),
				ITE(isZero,
				    unpackedFloat<t>::makeZero(format, sign),
				    ITE(!isSubnormal,
					ufNormal,
					ufSubnormalBase.normaliseUp(format) )))));

    POSTCONDITION(uf.valid(format));

    return uf;
  }


  template<class t>
    typename t::ubv pack (const typename t::fpt &format, const unpackedFloat<t> &uf) {
    typedef typename t::bwt bwt;
    typedef typename t::prop prop;
    typedef typename t::ubv ubv;
    typedef typename t::sbv sbv;

    PRECONDITION(uf.valid(format));

    // Sign
    ubv packedSign(uf.getSign());

    // Exponent
    bwt packedExWidth = format.packedExponentWidth();

    prop inNormalRange(uf.inNormalRange(format, prop(true)));
    INVARIANT(inNormalRange || uf.inSubnormalRange(format, prop(true)));     // Default values ensure this.
    //prop inSubnormalRange(uf.inSubnormalRange(format));        // Allowing this optimisation
    prop inSubnormalRange(!inNormalRange);

    probabilityAnnotation<t,prop>(inNormalRange, LIKELY);
    probabilityAnnotation<t,prop>(inSubnormalRange, UNLIKELY);
    
    sbv biasedExp(uf.getExponent() + unpackedFloat<t>::bias(format));
    // Will be correct for normal values only, subnormals may still be negative.
    ubv packedBiasedExp(biasedExp.toUnsigned().extract(packedExWidth - 1,0));

    ubv maxExp(ubv::allOnes(packedExWidth));
    ubv minExp(ubv::zero(packedExWidth));

    prop hasMaxExp(uf.getNaN() || uf.getInf());
    prop hasMinExp(uf.getZero() || inSubnormalRange);
    prop hasFixedExp(hasMaxExp || hasMinExp);
    
    ubv packedExp(ITE(hasFixedExp,
		      ITE(hasMaxExp, maxExp, minExp),
		      packedBiasedExp));
    

    // Significand
    bwt packedSigWidth = format.packedSignificandWidth();
    ubv unpackedSignificand(uf.getSignificand());

    INVARIANT(packedSigWidth == unpackedSignificand.getWidth() - 1);
    ubv dropLeadingOne(unpackedSignificand.extract(packedSigWidth - 1,0));
    ubv correctedSubnormal((unpackedSignificand >> (uf.getSubnormalAmount(format).toUnsigned().matchWidth(unpackedSignificand))).extract(packedSigWidth - 1,0));

    prop hasFixedSignificand(uf.getNaN() || uf.getInf() || uf.getZero());
    
    ubv packedSig(ITE(hasFixedSignificand,
		      ITE(uf.getNaN(),
			  unpackedFloat<t>::nanPattern(packedSigWidth),
			  ubv::zero(packedSigWidth)),
		      ITE(inNormalRange,
			  dropLeadingOne,
			  correctedSubnormal)));


    // Finish up
    ubv packed(packedSign.append(packedExp).append(packedSig));

    POSTCONDITION(packed.getWidth() == format.packedWidth());

    return packed;
  }

}

#endif
