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
** implementations.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 27/01/16
**
** Several different implementations of the floating-point operations.
** Object should only have static functions so that function pointers 
** can be used.  This means thee significand and exponent width should
** be fixed by the templating.
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
#include "symfpu/core/divide.h"
#include "symfpu/core/sqrt.h"
#include "symfpu/core/fma.h"
#include "symfpu/core/remainder.h"

#ifndef SYMFPU_IMPLEMENTATIONS
#define SYMFPU_IMPLEMENTATIONS


template <class execFloat>
class nativeFunctions {
 public :
  static execFloat abs (execFloat f);
  static execFloat max (execFloat f, execFloat g);
  static execFloat min (execFloat f, execFloat g);
  static execFloat sqrt (execFloat f);
  static execFloat rti (execFloat f);
  static execFloat fma (execFloat f, execFloat g, execFloat h);
  static execFloat rem (execFloat f, execFloat g);


  
  // These will work for any inbuilt float-type, specialisations need for others
  // Name changes because they are macros...

  static int fpClassify (execFloat f) { return fpclassify(f); }
  static int isNormal (execFloat f) { return isnormal(f); }
  static int isSubnormal (execFloat f) { return fpclassify(f) == FP_SUBNORMAL; }
  static int isZero (execFloat f);
  static int isInf (execFloat f) { return isinf(f); }
  static int isNaN (execFloat f) { return isnan(f); }
  static int isPositive (execFloat f) { return !isnan(f) && !signbit(f); }
  static int isNegative (execFloat f) { return !isnan(f) &&  signbit(f); }

  
};

// Specialisations for the obvious execFloat types
// To use other types as the reference, more instantiations will be needed

template <>
float nativeFunctions<float>::abs (float f) {
  return fabsf(f);
}

template <>
double nativeFunctions<double>::abs (double f) {
  return fabs(f);
}

template <>
long double nativeFunctions<long double>::abs (long double f) {
  return fabsl(f);
}



template <>
float nativeFunctions<float>::min (float f, float g) {
  return fminf(f,g);
}

template <>
double nativeFunctions<double>::min (double f, double g) {
  return fmin(f,g);
}

template <>
long double nativeFunctions<long double>::min (long double f, long double g) {
  return fminl(f,g);
}



template <>
float nativeFunctions<float>::max (float f, float g) {
  return fmaxf(f,g);
}

template <>
double nativeFunctions<double>::max (double f, double g) {
  return fmax(f,g);
}

template <>
long double nativeFunctions<long double>::max (long double f, long double g) {
  return fmaxl(f,g);
}



template <>
float nativeFunctions<float>::sqrt (float f) {
  return sqrtf(f);
}

template <>
double nativeFunctions<double>::sqrt (double f) {
  return sqrt(f);
}

template <>
long double nativeFunctions<long double>::sqrt (long double f) {
  return sqrtl(f);
}

template <>
float nativeFunctions<float>::rti (float f) {
  switch (fegetround()) {
  case FE_TONEAREST : return rintf(f); break;
  case FE_UPWARD : return ceilf(f); break;
  case FE_DOWNWARD : return floorf(f); break;
  case FE_TOWARDZERO : return truncf(f); break;
    //  case RNA : return roundf(f); break;
  default : assert(0); break;
  }
  return 0.0f;
}

template <>
double nativeFunctions<double>::rti (double f) {
  switch (fegetround()) {
  case FE_TONEAREST : return rint(f); break;
  case FE_UPWARD : return ceil(f); break;
  case FE_DOWNWARD : return floor(f); break;
  case FE_TOWARDZERO : return trunc(f); break;
    //  case RNA : return round(f); break;
  default : assert(0); break;
  }
  return 0.0;
}

template <>
long double nativeFunctions<long double>::rti (long double f) {
  switch (fegetround()) {
  case FE_TONEAREST : return rintl(f); break;
  case FE_UPWARD : return ceill(f); break;
  case FE_DOWNWARD : return floorl(f); break;
  case FE_TOWARDZERO : return truncl(f); break;
    //  case RNA : return roundl(f); break;
  default : assert(0); break;
  }
  return 0.0f;
}


// 1000000 tests
// libc fma 1739 bugs / 94 not sign of zero
// double 1861 bugs all not sign of zero
// float 1861 bugs all not sign of zero

template <>
float nativeFunctions<float>::fma (float f, float g, float h) {

  // vfmadd132ss
  return __builtin_fmaf(f,g,h);
  
  /*
    // Needs -mfma4 which is not the fma attribute in /proc/cpuinfo
  typedef float v4sf __attribute__ ((vector_size (16)));

  v4sf vf = {f,f,f,f};
  v4sf vg = {g,g,g,g};
  v4sf vh = {h,h,h,h};

  v4sf result = __builtin_ia32_vfmaddps (vf, vg, vh);

  return result[0];
  */

  /*
  float mult = f * g;
  return mult + h;
  */  

  /*
  double df = f;
  double dg = g;
  double dh = h;

  double mult = f * g;
  return (float)(mult + h);
  */
  
  //return fmaf(f,g,h);
}

template <>
double nativeFunctions<double>::fma (double f, double g, double h) {
  return fma(f,g,h);
}


template <>
long double nativeFunctions<long double>::fma (long double f, long double g, long double h) {
  return fmal(f,g,h);
}



template <>
int nativeFunctions<float>::isZero (float f) {
  return f == 0.0f;
}

template <>
int nativeFunctions<double>::isZero (double f) {
  return f == 0.0;
}

template <>
int nativeFunctions<long double>::isZero (long double f) {
  return f == 0.0l;
}


template <>
float nativeFunctions<float>::rem (float f, float g) {
  return remainderf(f,g);
}

template <>
double nativeFunctions<double>::rem (double f, double g) {
  return remainder(f,g);
}

template <>
long double nativeFunctions<long double>::rem (long double f, long double g) {
  return remainderl(f,g);
}




template <class execBV, class execFloat>
class native {
 public :

  static void setRoundingMode (const int roundingMode) {
    fesetround(roundingMode);
  }

  static execBV unpackPack (execBV bv) {
    return bv;
  }
  
  static execBV negate (execBV bv) {
    execFloat f = *((execFloat *)&bv);
    
    f = -f;
    
    return *((execBV *)&f);
  }
  
  static execBV absolute (execBV bv) {
    execFloat f = *((execFloat *)&bv);
    
    f = nativeFunctions<execFloat>::abs(f);
    
    return *((execBV *)&f);
  }
  
  static execBV sqrt (execBV bv) {
    execFloat f = *((execFloat *)&bv);
    
    f = nativeFunctions<execFloat>::sqrt(f);
    
    return *((execBV *)&f);
  }
  
  static execBV rti (execBV bv) {
    execFloat f = *((execFloat *)&bv);
    
    f = nativeFunctions<execFloat>::rti(f);
    
    return *((execBV *)&f);
  }
  
  static bool isNormal (execBV bv) {
    execFloat f = *((execFloat *)&bv);
    
    return nativeFunctions<execFloat>::isNormal(f);
  }
  
  static bool isSubnormal (execBV bv) {
    execFloat f = *((execFloat *)&bv);
    
    return nativeFunctions<execFloat>::isSubnormal(f);
  }
  
  static bool isZero (execBV bv) {
    execFloat f = *((execFloat *)&bv);
    
    return nativeFunctions<execFloat>::isZero(f);
  }
  
  static bool isInfinite (execBV bv) {
    execFloat f = *((execFloat *)&bv);
    
    return nativeFunctions<execFloat>::isInf(f);
  }
  
  static bool isNaN (execBV bv) {
    execFloat f = *((execFloat *)&bv);
    
    return nativeFunctions<execFloat>::isNaN(f);
  }
  
  static bool isPositive (execBV bv) {
    execFloat f = *((execFloat *)&bv);
    
    return nativeFunctions<execFloat>::isPositive(f);
  }
  
  static bool isNegative (execBV bv) {
    execFloat f = *((execFloat *)&bv);
    
    return nativeFunctions<execFloat>::isNegative(f);
  }
  
  static bool smtlibEqual (execBV bv1, execBV bv2) {
    execFloat f = *((execFloat *)(&bv1));
    execFloat g = *((execFloat *)(&bv2));
    
    return (bv1 == bv2) || (nativeFunctions<execFloat>::isNaN(f) && nativeFunctions<execFloat>::isNaN(g));
  }
  
  static bool ieee754Equal (execBV bv1, execBV bv2) {
    execFloat f = *((execFloat *)(&bv1));
    execFloat g = *((execFloat *)(&bv2));
    
    return (f == g);
  }
  
  static bool lessThan (execBV bv1, execBV bv2) {
    execFloat f = *((execFloat *)(&bv1));
    execFloat g = *((execFloat *)(&bv2));
    
    return (f < g);
  }
  
  static bool lessThanOrEqual  (execBV bv1, execBV bv2) {
    execFloat f = *((execFloat *)(&bv1));
    execFloat g = *((execFloat *)(&bv2));
    
    return (f <= g);
  }
  
  static execBV multiply (execBV bv1, execBV bv2) {
    execFloat f = *((execFloat *)(&bv1));
    execFloat g = *((execFloat *)(&bv2));
      
    execFloat h = f * g;
    
    return *((execBV *)&h);
  }

  static execBV add (execBV bv1, execBV bv2) {
    execFloat f = *((execFloat *)(&bv1));
    execFloat g = *((execFloat *)(&bv2));
    
    execFloat h = f + g;
    
    return *((execBV *)&h);
  }
  
  static execBV sub (execBV bv1, execBV bv2) {
    execFloat f = *((execFloat *)(&bv1));
    execFloat g = *((execFloat *)(&bv2));
    
    execFloat h = f - g;
    
    return *((execBV *)&h);
  }

  static execBV div (execBV bv1, execBV bv2) {
    execFloat f = *((execFloat *)(&bv1));
    execFloat g = *((execFloat *)(&bv2));
    
    execFloat h = f / g;
    
    return *((execBV *)&h);
  }

  static execBV max (execBV bv1, execBV bv2) {
    execFloat f = *((execFloat *)(&bv1));
    execFloat g = *((execFloat *)(&bv2));
    
    execFloat h = nativeFunctions<execFloat>::max(f,g);
    
    return *((execBV *)&h);
  }

  static execBV min (execBV bv1, execBV bv2) {
    execFloat f = *((execFloat *)(&bv1));
    execFloat g = *((execFloat *)(&bv2));
    
    execFloat h = nativeFunctions<execFloat>::min(f,g);
    
    return *((execBV *)&h);
  }

  static execBV fma (execBV bv1, execBV bv2, execBV bv3) {
    execFloat f = *((execFloat *)(&bv1));
    execFloat g = *((execFloat *)(&bv2));
    execFloat h = *((execFloat *)(&bv3));
    
    execFloat p = nativeFunctions<execFloat>::fma(f,g,h);
    
    return *((execBV *)&p);
  }

  static execBV rem (execBV bv1, execBV bv2) {
    execFloat f = *((execFloat *)(&bv1));
    execFloat g = *((execFloat *)(&bv2));
    
    execFloat h = nativeFunctions<execFloat>::rem(f,g);
    
    return *((execBV *)&h);
  }
  
};





template <class execBV, class traits>
class sympfuImplementation {

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

 protected :
  static rm * mode;
  static fpt * format;
    
 public :
  
  static void setRoundingMode (const int roundingMode) {
    if (mode != NULL) {
      delete mode;
    }

    switch (roundingMode) {
    case FE_TONEAREST :
      mode = new rm(traits::RNE());
      break;
    case FE_UPWARD :
      mode = new rm(traits::RTP());
      break;
    case FE_DOWNWARD :
      mode = new rm(traits::RTN());
      break;
    case FE_TOWARDZERO :
      mode = new rm(traits::RTZ());
      break;
      /* // Disabled until a suitable reference implementation is available
    case ??? :
      mode = new traits::rm(traits::RNA());
      break;
      */
    default :
      assert(0);
      break;
    }
  }

  static void setFormat (const fpt &newFormat) {
    if (format != NULL) {
      delete format;
    }
    format = new fpt(newFormat);
    return;
  }

  static void destroyFormat() {
    delete format;
    return;
  }
  
  static execBV unpackPack (const execBV bv) {
    ubv packed(bitsInExecBV(),bv);
    
    uf unpacked(symfpu::unpack<traits>(*format, packed));
    
    ubv repacked(symfpu::pack<traits>(*format, unpacked));
    
    return repacked.contents();
  }

  static execBV negate (execBV bv) {
    ubv packed(bitsInExecBV(),bv);
    
    uf unpacked(symfpu::unpack<traits>(*format, packed));
    
    uf negated(symfpu::negate<traits>(*format, unpacked));
    
    ubv repacked(symfpu::pack<traits>(*format, negated));
    
    return repacked.contents();
  }

  static execBV absolute (execBV bv) {
    ubv packed(bitsInExecBV(),bv);
 
    uf unpacked(symfpu::unpack<traits>(*format, packed));
    
    uf abs(symfpu::absolute<traits>(*format, unpacked));
    
    ubv repacked(symfpu::pack<traits>(*format, abs));
    
    return repacked.contents();
  }
  
  static execBV sqrt (execBV bv) {
    ubv packed(bitsInExecBV(),bv);
 
    uf unpacked(symfpu::unpack<traits>(*format, packed));
    
    uf sqrt(symfpu::sqrt<traits>(*format, *mode, unpacked));
    
    ubv repacked(symfpu::pack<traits>(*format, sqrt));
    
    return repacked.contents();
  }

  static execBV rti (execBV bv) {
    ubv packed(bitsInExecBV(),bv);
 
    uf unpacked(symfpu::unpack<traits>(*format, packed));
    
    uf rti(symfpu::roundToIntegral<traits>(*format, *mode, unpacked));
    
    ubv repacked(symfpu::pack<traits>(*format, rti));
    
    return repacked.contents();
  }

  static bool isNormal (execBV bv) {
    ubv packed(bitsInExecBV(),bv);
    
    uf unpacked(symfpu::unpack<traits>(*format, packed));
    
    prop result(symfpu::isNormal<traits>(*format, unpacked));
    
    return result;  
  }

  static bool isSubnormal (execBV bv) {
    ubv packed(bitsInExecBV(),bv);
    
    uf unpacked(symfpu::unpack<traits>(*format, packed));
    
    prop result(symfpu::isSubnormal<traits>(*format, unpacked));
    
    return result;
  }

  static bool isZero (execBV bv) {
    ubv packed(bitsInExecBV(),bv);
    
    uf unpacked(symfpu::unpack<traits>(*format, packed));
    
    prop result(symfpu::isZero<traits>(*format, unpacked));
    
    return result;
  }

  static bool isInfinite (execBV bv) {
    ubv packed(bitsInExecBV(),bv);
    
    uf unpacked(symfpu::unpack<traits>(*format, packed));
    
    prop result(symfpu::isInfinite<traits>(*format, unpacked));
    
    return result;
  }

  static bool isNaN (execBV bv) {
    ubv packed(bitsInExecBV(),bv);
    
    uf unpacked(symfpu::unpack<traits>(*format, packed));
    
    prop result(symfpu::isNaN<traits>(*format, unpacked));
    
    return result;
  }

  static bool isPositive (execBV bv) {
    ubv packed(bitsInExecBV(),bv);
    
    uf unpacked(symfpu::unpack<traits>(*format, packed));
    
    prop result(symfpu::isPositive<traits>(*format, unpacked));
    
    return result;
  }

  static bool isNegative (execBV bv) {
    ubv packed(bitsInExecBV(),bv);
    
    uf unpacked(symfpu::unpack<traits>(*format, packed));
    
    prop result(symfpu::isNegative<traits>(*format, unpacked));
    
    return result;
  }

  static bool smtlibEqual (execBV bv1, execBV bv2) {
    ubv packed1(bitsInExecBV(),bv1);
    ubv packed2(bitsInExecBV(),bv2);
    
    uf unpacked1(symfpu::unpack<traits>(*format, packed1));
    uf unpacked2(symfpu::unpack<traits>(*format, packed2));
    
    prop result(symfpu::smtlibEqual<traits>(*format, unpacked1, unpacked2));
    
    return result;
  }

  static bool ieee754Equal (execBV bv1, execBV bv2) {
    ubv packed1(bitsInExecBV(),bv1);
    ubv packed2(bitsInExecBV(),bv2);
    
    uf unpacked1(symfpu::unpack<traits>(*format, packed1));
    uf unpacked2(symfpu::unpack<traits>(*format, packed2));
    
    prop result(symfpu::ieee754Equal<traits>(*format, unpacked1, unpacked2));
    
    return result;
  }

  static bool lessThan (execBV bv1, execBV bv2) {
    ubv packed1(bitsInExecBV(),bv1);
    ubv packed2(bitsInExecBV(),bv2);
    
    uf unpacked1(symfpu::unpack<traits>(*format, packed1));
    uf unpacked2(symfpu::unpack<traits>(*format, packed2));
    
    prop result(symfpu::lessThan<traits>(*format, unpacked1, unpacked2));
    
    return result;
  }

  static bool lessThanOrEqual  (execBV bv1, execBV bv2) {
    ubv packed1(bitsInExecBV(),bv1);
    ubv packed2(bitsInExecBV(),bv2);
    
    uf unpacked1(symfpu::unpack<traits>(*format, packed1));
    uf unpacked2(symfpu::unpack<traits>(*format, packed2));
    
    prop result(symfpu::lessThanOrEqual<traits>(*format, unpacked1, unpacked2));
    
    return result;
  }

  static execBV multiply (execBV bv1, execBV bv2) {
    ubv packed1(bitsInExecBV(),bv1);
    ubv packed2(bitsInExecBV(),bv2);
    
    uf unpacked1(symfpu::unpack<traits>(*format, packed1));
    uf unpacked2(symfpu::unpack<traits>(*format, packed2));
    
    uf multiplied(symfpu::multiply<traits>(*format, *mode, unpacked1, unpacked2));
    
    ubv repacked(symfpu::pack<traits>(*format, multiplied));
    
    return repacked.contents();
  }

  static execBV add (execBV bv1, execBV bv2) {
    ubv packed1(bitsInExecBV(),bv1);
    ubv packed2(bitsInExecBV(),bv2);
    
    uf unpacked1(symfpu::unpack<traits>(*format, packed1));
    uf unpacked2(symfpu::unpack<traits>(*format, packed2));
    
    uf added(symfpu::add<traits>(*format, *mode, unpacked1, unpacked2, prop(true)));
    
    ubv repacked(symfpu::pack<traits>(*format, added));
    
    return repacked.contents();
  }

  static execBV sub (execBV bv1, execBV bv2) {
    ubv packed1(bitsInExecBV(),bv1);
    ubv packed2(bitsInExecBV(),bv2);
    
    uf unpacked1(symfpu::unpack<traits>(*format, packed1));
    uf unpacked2(symfpu::unpack<traits>(*format, packed2));
    
    uf added(symfpu::add<traits>(*format, *mode, unpacked1, unpacked2, prop(false)));
    
    ubv repacked(symfpu::pack<traits>(*format, added));
    
    return repacked.contents();
  }

  static execBV div (execBV bv1, execBV bv2) {
    ubv packed1(bitsInExecBV(),bv1);
    ubv packed2(bitsInExecBV(),bv2);
    
    uf unpacked1(symfpu::unpack<traits>(*format, packed1));
    uf unpacked2(symfpu::unpack<traits>(*format, packed2));
    
    uf added(symfpu::divide<traits>(*format, *mode, unpacked1, unpacked2));
    
    ubv repacked(symfpu::pack<traits>(*format, added));
    
    return repacked.contents();
  }

  #define INTELSSEMAXSTYLE true
  #define INTELSSEMINSTYLE false
  
  static execBV max (execBV bv1, execBV bv2) {
    ubv packed1(bitsInExecBV(),bv1);
    ubv packed2(bitsInExecBV(),bv2);
    
    uf unpacked1(symfpu::unpack<traits>(*format, packed1));
    uf unpacked2(symfpu::unpack<traits>(*format, packed2));
    
    uf max(symfpu::max<traits>(*format, unpacked1, unpacked2, INTELSSEMAXSTYLE));
    
    ubv repacked(symfpu::pack<traits>(*format, max));
    
    return repacked.contents();
  }

  static execBV min (execBV bv1, execBV bv2) {
    ubv packed1(bitsInExecBV(),bv1);
    ubv packed2(bitsInExecBV(),bv2);
    
    uf unpacked1(symfpu::unpack<traits>(*format, packed1));
    uf unpacked2(symfpu::unpack<traits>(*format, packed2));
    
    uf min(symfpu::min<traits>(*format, unpacked1, unpacked2, INTELSSEMINSTYLE));
    
    ubv repacked(symfpu::pack<traits>(*format, min));
    
    return repacked.contents();
  }

  static execBV fma (execBV bv1, execBV bv2, execBV bv3) {
    ubv packed1(bitsInExecBV(),bv1);
    ubv packed2(bitsInExecBV(),bv2);
    ubv packed3(bitsInExecBV(),bv3);
    
    uf unpacked1(symfpu::unpack<traits>(*format, packed1));
    uf unpacked2(symfpu::unpack<traits>(*format, packed2));
    uf unpacked3(symfpu::unpack<traits>(*format, packed3));
    
    uf fma(symfpu::fma<traits>(*format, *mode, unpacked1, unpacked2, unpacked3));
    
    ubv repacked(symfpu::pack<traits>(*format, fma));
    
    return repacked.contents();
  }

  static execBV rem (execBV bv1, execBV bv2) {
    ubv packed1(bitsInExecBV(),bv1);
    ubv packed2(bitsInExecBV(),bv2);
    
    uf unpacked1(symfpu::unpack<traits>(*format, packed1));
    uf unpacked2(symfpu::unpack<traits>(*format, packed2));
    
    uf min(symfpu::remainder<traits>(*format, unpacked1, unpacked2));
    
    ubv repacked(symfpu::pack<traits>(*format, min));
    
    return repacked.contents();
  }

  
  // The SMT-LIB notion of equality
  //bool compareFloat (execBV bv1, execBV bv2);

};

template <class execBV, class traits>
typename sympfuImplementation<execBV, traits>::rm * sympfuImplementation<execBV, traits>::mode = NULL;

template <class execBV, class traits>
typename sympfuImplementation<execBV, traits>::fpt * sympfuImplementation<execBV, traits>::format = NULL;


#endif
