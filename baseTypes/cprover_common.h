/* Part of SymFPU, see LICENSE for licensing information */
/*
** cprover_common.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 24/06/17
**
** Common CPROVER components.
**
*/

// CPROVER headers
#include <util/ieee_float.h>

#ifndef SYMFPU_CPROVER_COMMON
#define SYMFPU_CPROVER_COMMON

namespace symfpu {
  namespace cprover_common {

    typedef unsigned bitWidthType;


    class floatingPointTypeInfo : public ieee_float_spect {
    protected :
      //friend ite<proposition, floatingPointTypeInfo>;   // For ITE

    public :
      // CPROVER's format's significand doesn't include the hidden bit
    floatingPointTypeInfo(const ieee_float_spect &old) : ieee_float_spect(old) {}
    floatingPointTypeInfo(const floatbv_typet &t) : ieee_float_spect(t) {}
    floatingPointTypeInfo(unsigned exp, unsigned sig) : ieee_float_spect(sig - 1, exp) {}
      floatingPointTypeInfo(const floatingPointTypeInfo &old) : ieee_float_spect(old) {}
      
      bitWidthType exponentWidth(void) const    { return this->e; }
      bitWidthType significandWidth(void) const { return this->f + 1; }
      
      bitWidthType packedWidth(void) const            { return this->e + this->f + 1; }
      bitWidthType packedExponentWidth(void) const    { return this->e; }
      bitWidthType packedSignificandWidth(void) const { return this->f; }
    };

  };  
};

#endif
