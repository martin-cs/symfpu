/* Part of SymFPU, see LICENSE for licensing information */
/*
** shared.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 11/04/17
**
** Things that are used by multiple back-ends.
**
*/

#include <cstdint>
#include <cassert>

#ifndef SYMFPU_SHARED
#define SYMFPU_SHARED

namespace symfpu {
  namespace shared {

    
    // Must be able to contain the number of bit used in the bit-vector type to avoid overflow
    typedef uint64_t bitWidthType;

    // We can use bools for propositions
    typedef bool executable_proposition;

    // In SMT-LIB style -- significand includes hidden bit
    class floatingPointTypeInfo {
    private :
      bitWidthType exponentBits;
      bitWidthType significandBits;
      
    public :
      floatingPointTypeInfo (bitWidthType eb, bitWidthType sb) : exponentBits(eb), significandBits(sb) {
	assert(eb > 1);  // Not precondition as we don't have a traits class to use
	assert(sb > 1);
      }
      
      floatingPointTypeInfo (const floatingPointTypeInfo &old) : 
      exponentBits(old.exponentBits), significandBits(old.significandBits) {}
      
      floatingPointTypeInfo & operator= (const floatingPointTypeInfo &old) {
	this->exponentBits = old.exponentBits;
	this->significandBits = old.significandBits;
	
	return *this;
      }

      bitWidthType exponentWidth(void) const    { return this->exponentBits; }
      bitWidthType significandWidth(void) const { return this->significandBits; }

      
      bitWidthType packedWidth(void) const            { return this->exponentBits + this->significandBits; }
      bitWidthType packedExponentWidth(void) const    { return this->exponentBits; }
      bitWidthType packedSignificandWidth(void) const {	return this->significandBits - 1; }

      
    };
    

  }
}

#endif
