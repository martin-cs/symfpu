SymFPU Internal Design
======================

This document explains the design ideas behind SymFPU so that you can
read (and hopefully make sense of!) the code.  It assumes that you are
familiar with some of the ideas of IEEE-754:

- The five classes of number; normal, subnormal, zeros, infinities and
  NaN.
- Binary representation of floating-point numbers including details
  such as the biased representation of exponents and the "hidden bit".
- How the bits of normal and subnormal numbers are interpreted as
  rational numbers.
- The handling of special values in operations (for example,
  +Inf + Normal -> +Inf).
- The detail of how rounding mode handle inexact, overflow and
  underflow cases.
- How catastrophic cancellation happens and its consequences.

If you don't feel confident with these then the Handbook of
Floating-point or the SMT-LIB Theory of Floating-point would be a good
starting place.

Many of the high-level ideas in this document are also covered in:
"Building Better Bit-Blasting for Floating-Point Problems", TACAS'19

The first file to read is baseTypes/simpleExecutable.h .  This
provides a simple, executable implementation of the bit-vector types
that SymFPU uses a building blocks.  The class traits is the most
important as it gives an implementation of the interface that SymFPU
uses.  If you want to integrate SymFPU to use your own bit-vector
implementations (either executable or symbolic), this is the interface
you will need to implement.

The most basic approach to implementing floating-point operations is
as four separate stages:

UNPACK : Unpack the IEEE-754 representation of the float into an
 internal format.
OPERATE : Perform the operation on the unpacked format to generate an
 exact, but often larger result in the internal format.
ROUND : Round the result back to the desired size of internal format.
PACK : Pack the internal format back into the IEEE-754 format.

In practice many implementations blur the boundaries between these
steps or even skip some of them, as we will see below.

Although IEEE-754 formats are fixed, designers have greater freedom
with their internal format.  The simplest of these could be to just
split the IEEE-754 formatted bit-vector into its three components
(sign, exponent and significand).  However the design space is
actually much bigger, including:

- Flags for infinities, NaN, zeros, etc.
- Inclusion of the hidden bit implied for normal numbers.
- Removing the bias from the exponent.
- How subnormal numbers are handled.

When designing the internal format, there is a trade-off.  The more
elaborate it is, the more work unpack and pack steps have to do but it
may be possible to save effort on operate and round.

If the output of an operation goes straight into the next one and pack
followed by unpack leaves the internal representation unchanged, then
it is possible to skip the intermediate packing completely and keep
the intermediate value in the internal format.  Furthermore, if the
operation feeds into a comparison or classification operation, there
may not even be a need for a final pack step.  Also if the original
input numbers can be generated in the internal format then the only
place that packing and unpacking will be needed is when explicitly
converting to/from IEEE-754 format.

SymFPU assumes that explicit conversions will be relatively rare and
so uses a more complex internal format, with more elaborate pack and
unpack steps.  This simplifies a number of operations at very little
cost as pack and unpack are rarely used.  Key features of the internal
format (unpackedFloat) in SymFPU:

- It has flags for tracking when infinities, NaN and zero are
  represented.  This allows for fast "short-cut" handling of these
  cases.
- A default exponent and significand value are used for these special
  cases which reduces the number of cases that the operation code
  needs to handle.
- The exponent is unbiased so it can be handled as a normal signed
  bit-vector.  It is also extended slightly to help handling of
  subnormals.  The calculation of how much to extend the number is
  slightly complex.
- All subnormal numbers are represented in a normalised form.  Their
  significant is shifted up and exponent is reduced by the
  corresponding amount.
- The hidden bit is given explicitly and is always 1, even for
  subnormals after normalisation.

The key benefit of this representation is that in almost all cases,
the position of the leading 1 in the significand is known exactly or
almost exactly.  This means that we can skip searches for the leading
1 and expensive re-normalisation during almost all operations.  The
main exception is catastrophic cancellation but a re-normalisation
step is always needed in this case.  There is a slight cost to the
format as it makes rounding slightly more complex if the result is
subnormal, but this is a relatively rare case.

core/packing.h gives two functions which perform the unpack and pack
steps.  These illustrate the steps needed to convert to and from the
internal format.

core/unpackedFloat.h gives the unpackedFloat class which is used to
represent a float in the internal format.  unpackedFloat::valid is a
key method that checks (or generates the conditions that describe)
when a unpackedFloat is in the right format.  This is used for
creating variables in the internal format (by assuming they are valid)
and is a pre and post condition of operations.  One thing to note is
that the significand of subnormal values must have a number of
trailing 0's corresponding to what was past the LSB before normalisation.

core/sign.h gives the operations for taking the absolute value
(setting the sign bit to zero) and negation (flipping the sign bit).
These are good examples of the simplest possible operations.

core/classify.h gives the classification functions.  These mostly make
use of the flags in unpackedFloat.  The exception is testing for
normal and subnormal numbers.  Although it would be trivial to
generate a subnormal flag on unpacking, it is not naturally generated
by other operations in the same way as infinities and NaN.  So it
would require additional, possibly unnecessary checks to maintain.

core/compare.h gives ordering and equality operations.  These make
more use of the exponent and significand of the unpackedFloats and are
a good example of how these can be used.

core/multiply.h is the simplest of the core operations.  It is split
into three steps, computing the arithmetic result, rounding and then
adding the special cases.  This roughly corresponds to the operate and
round stages in the original model.  Special cases are added last so
they can "short cut" the whole computation.

rounder

divide

operations

add

fma

sqrt

remainder