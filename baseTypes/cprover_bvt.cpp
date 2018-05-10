/* Part of SymFPU, see LICENSE for licensing information */
#include "symfpu/baseTypes/cprover_bvt.h"

namespace symfpu {
  namespace cprover_bvt {
    bv_utilst *solver = NULL;

    #define BITS 32

    roundingMode traits::RNE (void) {
      return roundingMode(solver->build_constant(0/*ieee_floatt::ROUND_TO_EVEN*/, BITS));
    }
    
    roundingMode traits::RNA (void) {
      return roundingMode(solver->build_constant(4, BITS));
    }
    
    roundingMode traits::RTP (void) {
      return roundingMode(solver->build_constant(2/*ieee_floatt::ROUND_TO_PLUS_INF*/, BITS));
    }
    
    roundingMode traits::RTN (void) {
      return roundingMode(solver->build_constant(1/*ieee_floatt::ROUND_TO_MINUS_INF*/, BITS));
    }

    roundingMode traits::RTZ (void) {
      return roundingMode(solver->build_constant(3/*ieee_floatt::ROUND_TO_ZERO*/, BITS));
    }

  }
}
