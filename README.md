SymFPU : The Symbolic Floating Point Unit
=========================================

SymFPU is an implementation of the SMT-LIB theory of (IEEE-754)
floating-point in terms of bit-vector operations.  It is templated in
terms of the bit-vectors, propositions, floating-point formats and
rounding mode types used.  By providing different implementations of
these templates it is possible to:

- Run SymFPU as an arbitrary precision software floating-point
  library.  This is done for testing and consistent constant
  evaluation in some tools.  Note that it will not perform as well as
  a dedicated software floating-point library such as SoftFloat or
  MPFR.

- Run SymFPU "symbolically" to generate a representation of the
  bit-vector operations used in the floating-point operations.  By
  providing different implementations of the templates it is possible
  to produce the "symbolic" output in a variety of different formats.
  This is normally used for "bit-blasting" in SMT solvers.
  Alternatively you could use it to generate code in SMT-LIB, C or
  even Verilog if you wanted (SymFPU is not designed to produce good
  hardware so YMMV).

A considerable amount of effort has gone in to checking that these
encodings are correct and so please do report any discrepancies you
see.

The library is available under your choice of two licenses:
- The Free Software Foundations' GPL V3.
- The "3 clause" BSD license.
If this poses particular challenge for your application, please
contact the authors.



A Quick Start
-------------

To integrate SymFPU with your system, you will need to:

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

