/*
** Copyright (C) 2018 Martin Brain
**
** See the file LICENSE for licensing information.
*/

/*
** common.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 05/08/14
**
** Commonly used utility functions.
**
*/

#include <assert.h>
#include <stdint.h>

#include "symfpu/utils/properties.h"

#ifndef SYMFPU_COMMON
#define SYMFPU_COMMON

namespace symfpu {
  template <class T>
  T previousPowerOfTwo (T x) {
    assert(x > 1);
    //PRECONDITION(x > 1);

    T current = 1;
    T next = current << 1;
    
    while (next < x) {
      current = next;
      next <<= 1;
    }

    return current;
  }
  
  template <class T>
  T leftmostBit (T x) {
    assert(x > 1);
    //PRECONDITION(x > 1);

    T current = 1;
    T next = current << 1;
    
    while (next <= x) {
      current = next;
      next <<= 1;
    }

    return current;
  }

  
  // The number of bits required to represent a number
  //  == the position of the leading 0 + 1
  //  == ceil(log_2(value + 1))
  template <class T>
  T bitsToRepresent (const T value) {
    T i = 0; 
    //unsigned T working = *((unsigned T)&value);   // Implementation defined for signed types
    T working = value;

    while (working != 0) {
      ++i;
      working >>= 1;
    }

    return i;
  }

  template <class T>
  T positionOfLeadingOne (const T value) {
    //PRECONDITION(value != 0);
    assert(value != 0);

    T i = 0;
    //unsigned T working = *((unsigned T)&value);   // Implementation defined for signed types
    T working = value;

    while (working != 0) {
      ++i;
      working >>= 1;
    }

    return i - 1;
  }
}

#endif
