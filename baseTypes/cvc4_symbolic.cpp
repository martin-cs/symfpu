/* Part of SymFPU, see LICENSE for licensing information */
/*
** cvc4_symbolic.cpp
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 05/01/15
**
** Non-templated functions for linking.
** It's best not to ask how this is built...
**
*/

#include "symfpu/baseTypes/cvc4_symbolic.h"

// CVC4 headers
#include "theory/fp/fp_converter.h"

namespace symfpu {
  namespace cvc4_symbolic {

    roundingMode traits::RNE (void) { return roundingMode(0x01); };
    roundingMode traits::RNA (void) { return roundingMode(0x02); };
    roundingMode traits::RTP (void) { return roundingMode(0x04); };
    roundingMode traits::RTN (void) { return roundingMode(0x08); };
    roundingMode traits::RTZ (void) { return roundingMode(0x10); };


    template <>
    bitVector<true> bitVector<true>::maxValue (const bitWidthType &w) {
      bitVector<true> leadingZero(bitVector<true>::zero(1));
      bitVector<true> base(bitVector<true>::allOnes(w-1));
      
      return bitVector<true>(::CVC4::NodeManager::currentNM()->mkNode(::CVC4::kind::BITVECTOR_CONCAT, leadingZero.node, base.node));
    }

    template <>
    bitVector<false> bitVector<false>::maxValue (const bitWidthType &w) {
      return bitVector<false>::allOnes(w);
    }

    template <>
    bitVector<true> bitVector<true>::minValue (const bitWidthType &w) {
      bitVector<true> leadingOne(bitVector<true>::one(1));
      bitVector<true> base(bitVector<true>::zero(w-1));
      
      return bitVector<true>(::CVC4::NodeManager::currentNM()->mkNode(::CVC4::kind::BITVECTOR_CONCAT, leadingOne.node, base.node));
    }

    template <>
    bitVector<false> bitVector<false>::minValue (const bitWidthType &w) {
      return bitVector<false>::zero(w);
    }    
    
  };

  template <>
  cvc4_symbolic::traits::ubv orderEncode<cvc4_symbolic::traits, cvc4_symbolic::traits::ubv> (const cvc4_symbolic::traits::ubv &b) {
    return orderEncodeBitwise<cvc4_symbolic::traits, cvc4_symbolic::traits::ubv>(b);
  }

  template <>
  stickyRightShiftResult<cvc4_symbolic::traits> stickyRightShift (const cvc4_symbolic::traits::ubv &input, const cvc4_symbolic::traits::ubv &shiftAmount) {
    return stickyRightShiftBitwise<cvc4_symbolic::traits>(input, shiftAmount);
  }

  template <>
  void probabilityAnnotation<cvc4_symbolic::traits, cvc4_symbolic::traits::prop> (const cvc4_symbolic::traits::prop &p, const probability &pr) {
    ::CVC4::theory::fp::FpConverter::currentConverter->registerProbabilityAnnotation(p.getNode(), pr);
    return;
  }

};
