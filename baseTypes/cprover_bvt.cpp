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
