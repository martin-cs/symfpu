SymFPU : The Symbolic Floating Point Unit
=========================================

SymFPU is an implementation of the SMT-LIB / IEEE-754 operations in
terms of bit-vector operations.  It is templated in terms of the
bit-vectors, propositions, floating-point formats and rounding mode
types used.  This allow the same code to be executed as an arbitrary
precision "SoftFloat" library (although it's performance would not be
good) or to be used to build symbolic representation of floating-point
operations suitable for use in "bit-blasting" SMT solvers (you could
also generate circuits from them but again, performance will likely
not be good).

A considerable amount of effort has gone in to checking that these
encodings are correct and so please do report any discrepancies you
see.

The library is Free Software licensed under the GPL V3.  If this poses
a particular challenge for your application, please contact the author.


A Quick Start
-------------

1. Create a "back-end", a class with the following members:

```
class traits {
  public :
  // The base types to use.
  // Must implement the SMT-LIB-like interfaces used by other back-ends
  typedef YourBitWidthType bwt;
  typedef YourRoundingMode rm;
  typedef YourFloatingPointTypeInfo fpt;
  typedef YourProposition prop;
  typedef YourSignedBitVector sbv;
  typedef YourUnsignedBitVector ubv;

  // Return an instance of each rounding mode.
  static rm RNE(void);
  static rm RNA(void);
  static rm RTP(void);
  static rm RTN(void);
  static rm RTZ(void);

  // Handle various invariants.
  // These can be empty to start with.
  static void precondition(const bool b) { assert(b); }
  static void postcondition(const bool b) { assert(b); }
  static void invariant(const bool b) { assert(b); }
  static void precondition(const prop &p) {}
  static void postcondition(const prop &p) {}
  static void invariant(const prop &p) {}
};
```

If you are stuck; start from one of the existing (symbolic or literal)
ones.  `baseTypes/shared.h` gives suitable implementations of `fpt` and
`bwt`, plus if it is executable, you can use bool for prop.


2. Make sure symfpu is on your path and include the following headers:

```
#include "symfpu/core/unpackedFloat.h"
#include "symfpu/core/packing.h"
#include "symfpu/core/add.h"
```

3. To generate a (32-bit) unsigned bit-vector containing the single
precision addition of two other 32-bit bit-vectors, use the following code:

```
fpt format(8,24);
ubv packed1 = ...  // Must be 32-bit.
ubv packed2 = ...  // Must be 32-bit.

uf unpacked1(symfpu::unpack<traits>(format, packed1));
uf unpacked2(symfpu::unpack<traits>(format, packed2));
    
uf added(symfpu::add<traits>(format, traits::RNE(), unpacked1, unpacked2, prop(true)));
    
ubv repacked(symfpu::pack<traits>(format, added));
```

See `applications/implementation.h` for examples of other operations
(although, really, it is pretty similar).

