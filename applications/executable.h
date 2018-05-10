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

// DEPRECIATED : USE implementation.h instead

/*
** executable.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 06/08/14
**
** A set of basic functions constructed from the simple executable
** implementation of bit-vectors.
**
*/

#include <limits.h>
#include <math.h>
#include <fenv.h>

#include "symfpu/core/unpackedFloat.h"
#include "symfpu/core/packing.h"
#include "symfpu/core/sign.h"
#include "symfpu/core/classify.h"
#include "symfpu/core/compare.h"
#include "symfpu/core/multiply.h"
#include "symfpu/core/add.h"



#ifndef SYMFPU_EXECUTABLE
#define SYMFPU_EXECUTABLE



template <class execBV, class execFloat, class traits>
  class executableTests {

 public :
    // Wrapped in a struct to make type scoping easier
    // and to save on typenames.
    // Object is stateless.

    typedef typename traits::rm rm;
    typedef typename traits::bwt bwt;
    typedef typename traits::fpt fpt;
    typedef typename traits::ubv ubv;
    typedef typename traits::prop prop;
    typedef symfpu::unpackedFloat<traits> uf;

    static bwt bitsInExecBV () {
      return sizeof(execBV) * CHAR_BIT;
    }



    static execBV unpackPack (const fpt &format, const execBV bv) {
      ubv packed(bitsInExecBV(),bv);
      
      uf unpacked(symfpu::unpack<traits>(format, packed));
      
      ubv repacked(symfpu::pack<traits>(format, unpacked));
      
      return repacked.contents();
    }

    static execBV unpackPackReference (const fpt &, execBV bv) {
      return bv;
    }



    static execBV negate (const fpt &format, execBV bv) {
      ubv packed(bitsInExecBV(),bv);
 
      uf unpacked(symfpu::unpack<traits>(format, packed));

      uf negated(symfpu::negate<traits>(format, unpacked));

      ubv repacked(symfpu::pack<traits>(format, negated));

      return repacked.contents();
    }

    static execBV negateReference (const fpt &, execBV bv) {
      execFloat f = *((execFloat *)&bv);

      f = -f;

      return *((execBV *)&f);
    }



    static execBV absolute (const fpt &format, execBV bv) {
      ubv packed(bitsInExecBV(),bv);
 
      uf unpacked(symfpu::unpack<traits>(format, packed));

      uf abs(symfpu::absolute<traits>(format, unpacked));

      ubv repacked(symfpu::pack<traits>(format, abs));

      return repacked.contents();
    }

    static execBV absoluteReference (const fpt &, execBV bv) {
      execFloat f = *((execFloat *)&bv);

      f = nativeFunctions<execFloat>::abs(f);

      return *((execBV *)&f);
    }


    static bool isNormal (const fpt &format, execBV bv) {
      ubv packed(bitsInExecBV(),bv);
 
      uf unpacked(symfpu::unpack<traits>(format, packed));

      prop result(symfpu::isNormal<traits>(format, unpacked));

      return result;  
    }

    static bool isNormalReference (const fpt &, execBV bv) {
      execFloat f = *((execFloat *)&bv);

      return nativeFunctions<execFloat>::isNormal(f);
    }



    static bool isSubnormal (const fpt &format, execBV bv) {
      ubv packed(bitsInExecBV(),bv);
 
      uf unpacked(symfpu::unpack<traits>(format, packed));

      prop result(symfpu::isSubnormal<traits>(format, unpacked));

      return result;
    }

    static bool isSubnormalReference (const fpt &, execBV bv) {
      execFloat f = *((execFloat *)&bv);

      return nativeFunctions<execFloat>::isSubnormal(f);
    }




    static bool isZero (const fpt &format, execBV bv) {
      ubv packed(bitsInExecBV(),bv);
 
      uf unpacked(symfpu::unpack<traits>(format, packed));

      prop result(symfpu::isZero<traits>(format, unpacked));

      return result;
    }

    static bool isZeroReference (const fpt &, execBV bv) {
      execFloat f = *((execFloat *)&bv);

      return nativeFunctions<execFloat>::isZero(f);
    }



    static bool isInfinite (const fpt &format, execBV bv) {
      ubv packed(bitsInExecBV(),bv);
 
      uf unpacked(symfpu::unpack<traits>(format, packed));

      prop result(symfpu::isInfinite<traits>(format, unpacked));

      return result;
    }

    static bool isInfiniteReference (const fpt &, execBV bv) {
      execFloat f = *((execFloat *)&bv);

      return nativeFunctions<execFloat>::isInf(f);
    }



    static bool isNaN (const fpt &format, execBV bv) {
      ubv packed(bitsInExecBV(),bv);
 
      uf unpacked(symfpu::unpack<traits>(format, packed));

      prop result(symfpu::isNaN<traits>(format, unpacked));

      return result;
    }

    static bool isNaNReference (const fpt &, execBV bv) {
      execFloat f = *((execFloat *)&bv);

      return nativeFunctions<execFloat>::isNaN(f);
    }




    static bool isPositive (const fpt &format, execBV bv) {
      ubv packed(bitsInExecBV(),bv);
 
      uf unpacked(symfpu::unpack<traits>(format, packed));

      prop result(symfpu::isPositive<traits>(format, unpacked));

      return result;
    }

    static bool isPositiveReference (const fpt &, execBV bv) {
      execFloat f = *((execFloat *)&bv);

      return nativeFunctions<execFloat>::isPositive(f);
    }



    static bool isNegative (const fpt &format, execBV bv) {
      ubv packed(bitsInExecBV(),bv);
 
      uf unpacked(symfpu::unpack<traits>(format, packed));

      prop result(symfpu::isNegative<traits>(format, unpacked));

      return result;
    }

    static bool isNegativeReference (const fpt &, execBV bv) {
      execFloat f = *((execFloat *)&bv);

      return nativeFunctions<execFloat>::isNegative(f);
    }


    static bool smtlibEqual (const fpt &format, execBV bv1, execBV bv2) {
      ubv packed1(bitsInExecBV(),bv1);
      ubv packed2(bitsInExecBV(),bv2);
 
      uf unpacked1(symfpu::unpack<traits>(format, packed1));
      uf unpacked2(symfpu::unpack<traits>(format, packed2));

      prop result(symfpu::smtlibEqual<traits>(format, unpacked1, unpacked2));

      return result;
    }

    static bool smtlibEqualReference (const fpt &, execBV bv1, execBV bv2) {
      execFloat f = *((execFloat *)(&bv1));
      execFloat g = *((execFloat *)(&bv2));

      return (bv1 == bv2) || (nativeFunctions<execFloat>::isNaN(f) && nativeFunctions<execFloat>::isNaN(g));
    }



    static bool ieee754Equal (const fpt &format, execBV bv1, execBV bv2) {
      ubv packed1(bitsInExecBV(),bv1);
      ubv packed2(bitsInExecBV(),bv2);
 
      uf unpacked1(symfpu::unpack<traits>(format, packed1));
      uf unpacked2(symfpu::unpack<traits>(format, packed2));

      prop result(symfpu::ieee754Equal<traits>(format, unpacked1, unpacked2));

      return result;
    }


    static bool ieee754EqualReference (const fpt &, execBV bv1, execBV bv2) {
      execFloat f = *((execFloat *)(&bv1));
      execFloat g = *((execFloat *)(&bv2));

      return (f == g);
    }



    static bool lessThan (const fpt &format, execBV bv1, execBV bv2) {
      ubv packed1(bitsInExecBV(),bv1);
      ubv packed2(bitsInExecBV(),bv2);
 
      uf unpacked1(symfpu::unpack<traits>(format, packed1));
      uf unpacked2(symfpu::unpack<traits>(format, packed2));

      prop result(symfpu::lessThan<traits>(format, unpacked1, unpacked2));

      return result;
    }

    static bool lessThanReference (const fpt &, execBV bv1, execBV bv2) {
      execFloat f = *((execFloat *)(&bv1));
      execFloat g = *((execFloat *)(&bv2));

      return (f < g);
    }





    static bool lessThanOrEqual  (const fpt &format, execBV bv1, execBV bv2) {
      ubv packed1(bitsInExecBV(),bv1);
      ubv packed2(bitsInExecBV(),bv2);
 
      uf unpacked1(symfpu::unpack<traits>(format, packed1));
      uf unpacked2(symfpu::unpack<traits>(format, packed2));

      prop result(symfpu::lessThanOrEqual<traits>(format, unpacked1, unpacked2));

      return result;
    }


    static bool lessThanOrEqualReference  (const fpt &, execBV bv1, execBV bv2) {
      execFloat f = *((execFloat *)(&bv1));
      execFloat g = *((execFloat *)(&bv2));

      return (f <= g);
    }




    static execBV multiply (const fpt &format, const rm &roundingMode, execBV bv1, execBV bv2) {
      ubv packed1(bitsInExecBV(),bv1);
      ubv packed2(bitsInExecBV(),bv2);
 
      uf unpacked1(symfpu::unpack<traits>(format, packed1));
      uf unpacked2(symfpu::unpack<traits>(format, packed2));

      uf multiplied(symfpu::multiply<traits>(format, roundingMode, unpacked1, unpacked2));
 
      ubv repacked(symfpu::pack<traits>(format, multiplied));

      return repacked.contents();
    }


    static execBV multiplyReference (const fpt &, const rm &roundingMode, execBV bv1, execBV bv2) {
      execFloat f = *((execFloat *)(&bv1));
      execFloat g = *((execFloat *)(&bv2));

      fesetround(roundingMode.getValue());

      execFloat h = f * g;

      return *((execBV *)&h);
    }

    static execBV add (const fpt &format, const rm &roundingMode, execBV bv1, execBV bv2) {
      ubv packed1(bitsInExecBV(),bv1);
      ubv packed2(bitsInExecBV(),bv2);
 
      uf unpacked1(symfpu::unpack<traits>(format, packed1));
      uf unpacked2(symfpu::unpack<traits>(format, packed2));

      uf added(symfpu::add<traits>(format, roundingMode, unpacked1, unpacked2, prop(true)));
 
      ubv repacked(symfpu::pack<traits>(format, added));

      return repacked.contents();
    }


    static execBV addReference (const fpt &, const rm &roundingMode, execBV bv1, execBV bv2) {
      execFloat f = *((execFloat *)(&bv1));
      execFloat g = *((execFloat *)(&bv2));

      fesetround(roundingMode.getValue());

      execFloat h = f + g;

      return *((execBV *)&h);
    }


    static execBV sub (const fpt &format, const rm &roundingMode, execBV bv1, execBV bv2) {
      ubv packed1(bitsInExecBV(),bv1);
      ubv packed2(bitsInExecBV(),bv2);
 
      uf unpacked1(symfpu::unpack<traits>(format, packed1));
      uf unpacked2(symfpu::unpack<traits>(format, packed2));

      uf added(symfpu::add<traits>(format, roundingMode, unpacked1, unpacked2, prop(false)));
 
      ubv repacked(symfpu::pack<traits>(format, added));

      return repacked.contents();
    }


    static execBV subReference (const fpt &, const rm &roundingMode, execBV bv1, execBV bv2) {
      execFloat f = *((execFloat *)(&bv1));
      execFloat g = *((execFloat *)(&bv2));

      fesetround(roundingMode.getValue());

      execFloat h = f - g;

      return *((execBV *)&h);
    }

    // The SMT-LIB notion of equality
    //bool compareFloat (execBV bv1, execBV bv2);

  };




#endif
