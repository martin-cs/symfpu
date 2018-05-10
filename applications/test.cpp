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
** test.cpp
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 06/08/14
**
** The main loop of the test program
**
*/

#include <math.h>
#include <ieee754.h>
#include <fenv.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "symfpu/baseTypes/simpleExecutable.h"

#include "symfpu/applications/implementations.h"


/*** Test Vector Generation ***/

#define NUMBER_OF_FLOAT_TESTS 124
static float floatTestValue [NUMBER_OF_FLOAT_TESTS] = {
  0x0p+0f, -0x0p+0f,                        // Zeros
  0x1p+0f, -0x1p+0f,                        // Ones
  INFINITY, -INFINITY, NAN, -NAN,           // Special values
  M_E, M_LOG2E, M_LOG10E, M_LN2, M_LN10,    // Mathematical constants
  M_PI, M_PI_2, M_PI_4, M_1_PI, M_2_PI,
  M_2_SQRTPI, M_SQRT2, M_SQRT1_2,
   0x1.fffffcp-127,  0x1.000004p-127,  0x1p-127,  0x1p-148,  0x1p-149,  // Subnormals
  -0x1.fffffcp-127, -0x1.000004p-127, -0x1p-127, -0x1p-148, -0x1p-149,
   0x1.ffffffcp-103,  0x1.ffffffcp-99,  0x1.ffffffcp-75,      // Normals
   0x1.ffffffcp-51,   0x1.ffffffcp-27,  0x1.ffffffcp-03,
   0x1.ffffffcp+51,   0x1.ffffffcp+27,  0x1.ffffffcp+03,
   0x1.ffffffcp+103,  0x1.ffffffcp+99,  0x1.ffffffcp+75,
  -0x1.ffffffcp-103, -0x1.ffffffcp-99, -0x1.ffffffcp-75,
  -0x1.ffffffcp-51,  -0x1.ffffffcp-27, -0x1.ffffffcp-03,
  -0x1.ffffffcp+51,  -0x1.ffffffcp+27, -0x1.ffffffcp+03,
  -0x1.ffffffcp+103, -0x1.ffffffcp+99, -0x1.ffffffcp+75,
   0x1.fffffep+127,                                            // From the CBMC regression tests
   0x1.4p+4,
   0x1.fffffep-105f,
   0x1.0p-23,
   0x1.0p-126f,
   0x1.fffffep-3f,
   0x1.fffffep+123f,
   0x1.000002p+124f,
   0x1.0p+124f,
   0x1.fffffep-3f,
   0x1p+124f,
   0x1p-126f,
   0x1.000000p-63f,
   0x1.fffffep-64f,
   0x1.084c64p-63f,
   0x1.efec9ep-64f,
   0x1.47e8c2p-63f,
   0x1.8fb86cp-64f,
   0x1.1fcf1cp-63f,
   0x1.c769c0p-64f,
   0x1.b1fffcp-63f,
   0x1.2e025ep-64f,
   0x1.000000p-62f,
   0x1.fffffep-65f,
   0x1.000000p-61f,
   0x1.fffffep-66f,
   0x1.000000p-50f,
   0x1.fffffep-77f,
   0x1.000000p-30f,
   0x1.fffffep-97f,
   0x1.000000p-10f,
   0x1.fffffep-117f,
   0x1.000000p-1f,
   0x1.fffffep-126f,
   0x1.9e0c22p-101f,
  -0x1.3c9014p-50f,
   0x1.8p-24,
   0x1.fffffep-1f,
   0x1.000002p+0f,
   0x1.7ffffep-24f,
   0x1.800002p-24f,
   0x1.800000p-24f,
  -0x1.0p-127f,
  -0x1.6b890ep+29,
   0x1.6b890ep+29,
   0x1.000002p+25f,
   0x1.000004p+25f,
   0x1.0017ecp+22f,            // Divide!
   0x1.7ffff8p+21f,
   0x1.557542p+0f,
   0x1p+120f,                  // Distributivity
  -0x1p+120f,
   0x1.46p+7f,                 // e^{pi * sqrt{163}} = 640320^3 + 744
   0x1.38a8p+19f,
   0x1.74p+9f,
   0x1.79999ap+5f,             // 47.2
   0x1.99999ap-4f,             // For the patriots..
   0x1.5p+5f,
   0x1.7p+4f,
   0x1.fffffep+24,             // To test carry on increment
   0x1.fffffcp+24,
   0x1.0p-1,                   // Half for a laugh
   0x1.000002p-75f,            // To test rounding on multiply
   0x1.0p-75f,
   0x1.8p+0f,                  // Carry in to top bit of fraction when half is added
   0x1.fffffep+125f,           // Hunt a specific bug
   0x1.fffffep+126f,
   0x1.8p+1f,
   0x1.000004p+125f
};

float getTestValue (uint64_t index) {
  if (index < NUMBER_OF_FLOAT_TESTS) {
    return floatTestValue[index];
  } else {
    index -= NUMBER_OF_FLOAT_TESTS;

    ieee754_float c;

    // The low bits give the sign and exponent so that a wide range is covered quickly
    c.ieee.negative = index & 0x1;
    index >>= 1;

    c.ieee.exponent = ((index & 0xFE) >> 1) | ((index & 0x1) << 7);
    index >>= 8;

    assert(index < (1 << 23));

    // Use the next bit to optionally negates the word
    //index = (index & 0x1) ? ~(index >> 1) : (index >> 1);

    // Even bits affect the MSBs of the significand, odd bits affect the LSBs of the significand
    uint32_t lsb = 0x0;
    uint32_t msb = 0x0;


    for (uint32_t i = 0; i < 23; ++i) {
      if (i & 0x1) {
        lsb |= (((index >> i) & 0x1) << (i >> 1));
      } else {
        msb <<= 1;
        msb |= ((index >> i) & 0x1);
      }
    }

    assert((lsb & (msb << 11)) == 0);

    c.ieee.mantissa = lsb | (msb << 11);

    return c.f;
  }
}




/*** Types ***/

// We are testing the 'simple executable' back-end
typedef symfpu::simpleExecutable::traits traits;
typedef traits::fpt fpt;

// We are also testing it only for single precision
traits::fpt singlePrecisionFormatObject(8,24);

// Thus we are using the functions in ...
typedef sympfuImplementation<uint32_t, traits> singlePrecisionExecutableSymfpu;

// Which we then compare against
typedef native<uint32_t, float> singlePrecisionHardware;



/*** Output helpers ***/

#define BUFFERLENGTH 256

FILE * openOutputFile (const char *string, const char *name, const char *roundingMode, const uint64_t testNumber) {
  char filename[BUFFERLENGTH];

  snprintf(filename, BUFFERLENGTH, string, name, roundingMode, testNumber);

  FILE *out = fopen(filename, "w");
  if (out == NULL) { 
    fprintf(stderr,"Couldn't open %s\n", filename);
    perror("Can't open file for writing");
    exit(1);
  }

  return out;
}

FILE * startOutputC (const char *name, const char *roundingMode, const uint64_t testNumber) {
  FILE * out = openOutputFile("testC-%s-%s-%d.c", name, roundingMode, testNumber);

  // TODO : version and date
  fprintf(out, "// Test case created by symfpu for operation %s, rounding mode %s, test %lx\n\n", name, roundingMode, testNumber);
  fprintf(out, "#include <assert.h>\n");
  fprintf(out, "#include <math.h>\n");
  fprintf(out, "#include <fenv.h>\n\n");

  fprintf(out, "int compare (float ref, float computed) {\n\n");

  fprintf(out, "int isrefnan = isnan(ref);\n");
  fprintf(out, "int iscomputednan = isnan(computed);\n");
  fprintf(out, "int equal = (ref == computed);\n");
  fprintf(out, "int signref = signbit(ref);\n");
  fprintf(out, "int signcomp = signbit(computed);\n");

  fprintf(out, "return ((isrefnan && iscomputednan) || \n");
  fprintf(out, "        (equal && ((signref == 0 && signcomp == 0) || \n");
  fprintf(out, "                   (signref != 0 && signcomp != 0))));\n");
  fprintf(out, "}\n\n");

  fprintf(out, "int main (void) {\n");

  return out;
}

void finishOutputC (FILE *out) {
  fprintf(out,"return 1;\n");
  fprintf(out, "}\n\n");
  fclose(out);
  return;
}

void printFloatC (FILE *out, float f) {
  if (isnan(f)) {
    fprintf(out, "NAN");
  } else if (isinf(f)) {
    fprintf(out, "%cINFINITY", (signbit(f) == 0) ? ' ' : '-');
  } else {
    fprintf(out, "%af", f);
  }

  return;
}

void printFloatSMT (FILE *out, uint32_t f) {
  fprintf(out, "(fp (_ bv%d 1) (_ bv%d 8) (_ bv%d 23))",
	  (f & 0x80000000) >> 31,
	  (f & 0x7F800000) >> 23,
	  (f & 0x007FFFFF) >> 0);
  return;
}

FILE * startOutputSMT (const char *name, const char *roundingMode, const uint64_t testNumber) {
  FILE * out = openOutputFile("testSMT-%s-%s-%d.smt2", name, roundingMode, testNumber);

  // TODO : version and date
  fprintf(out, "(set-logic ALL_SUPPORTED)\n");
  fprintf(out, "; Should be SAT\n");

  return out;
}

void finishOutputSMT (FILE *out) {

  fprintf(out, "(check-sat)\n");

  fclose(out);
  return;
}




/*** Test Execution ***/


typedef uint32_t (*unaryFunctionTFP) (uint32_t);

template <unaryFunctionTFP test, unaryFunctionTFP ref>
void unaryFunctionTest (const int verbose, const uint64_t start, const uint64_t end) {
  uint64_t i;

  for (i = start; i < end; ++i) {
      float f = getTestValue(i);
      
      uint32_t input = *((uint32_t *)(&f));
     
      uint32_t reference = ref(input);
      uint32_t computed = test(input);
      
      if (verbose || !singlePrecisionHardware::smtlibEqual(computed, reference)) {
	fprintf(stdout,"vector[%d] ", (uint32_t)i);
	fprintf(stdout,"input = 0x%x, computed = 0x%x, real = 0x%x\n", input, computed, reference);
	fflush(stdout);
      }

      if ((i & 0xFFFF) == 0) {
	fprintf(stdout,".");
	fflush(stdout);
      }
  }

  return;
}

template <unaryFunctionTFP ref>
void unaryFunctionPrintC (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *cPrintString, const char *) {
  uint64_t i;
  FILE * out;

  for (i = start; i < end; ++i) {
      float f = getTestValue(i);
      uint32_t input = *((uint32_t *)(&f));
     
      uint32_t reference = ref(input);
      float fref = *((float *)(&reference));

      out = startOutputC(name, "NA", i);

      fprintf(out, "float f = ");
      printFloatC(out, f);
      fprintf(out,";\n");

      fprintf(out, "float ref = ");
      printFloatC(out, fref);
      fprintf(out,";\n");

      fprintf(out, "float computed = %s;\n", cPrintString);
      fprintf(out, "assert(compare(ref, computed));\n");

      finishOutputC(out);

  }

  return;
}

template <unaryFunctionTFP ref>
void unaryFunctionPrintSMT (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *SMTPrintString, const char *) {
  uint64_t i;
  FILE * out;

  for (i = start; i < end; ++i) {
      float f = getTestValue(i);
      uint32_t input = *((uint32_t *)(&f));
     
      uint32_t reference = ref(input);
      //float fref = *((float *)(&reference));

      out = startOutputSMT(name, "NA", i);

      fprintf(out, "(define-fun f () Float32 ");
      printFloatSMT(out, input);
      fprintf(out, ")\n");

      fprintf(out, "(define-fun ref () Float32 ");
      printFloatSMT(out, reference);
      fprintf(out, ")\n");

      fprintf(out, "(define-fun result () Float32 %s )\n", SMTPrintString);
      fprintf(out, "(assert (= ref result))\n");

      finishOutputSMT(out);
  }

  return;
}

typedef uint32_t (*unaryRoundedFunctionTFP) (uint32_t);

template <unaryRoundedFunctionTFP test, unaryRoundedFunctionTFP ref>
void unaryRoundedFunctionTest (const int verbose, const uint64_t start, const uint64_t end) {
  uint64_t i;

  for (i = start; i < end; ++i) {
      float f = getTestValue(i);
      
      uint32_t input = *((uint32_t *)(&f));
     
      uint32_t reference = ref(input);
      uint32_t computed = test(input);
      
      if (verbose || !singlePrecisionHardware::smtlibEqual(computed, reference)) {
	fprintf(stdout,"vector[%d] ", (uint32_t)i);
	fprintf(stdout,"input = 0x%x, computed = 0x%x, real = 0x%x\n", input, computed, reference);
	fflush(stdout);
      }

      if ((i & 0xFFFF) == 0) {
	fprintf(stdout,".");
	fflush(stdout);
      }
  }

  return;
}

template <unaryRoundedFunctionTFP ref>
void unaryRoundedFunctionPrintC (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *cPrintString, const char *roundingModeString) {
  uint64_t i;
  FILE * out;

  for (i = start; i < end; ++i) {
      float f = getTestValue(i);
      uint32_t input = *((uint32_t *)(&f));
     
      uint32_t reference = ref(input);
      float fref = *((float *)(&reference));

      out = startOutputC(name, roundingModeString, i);

      fprintf(out, "float f = ");
      printFloatC(out, f);
      fprintf(out,";\n");

      fprintf(out, "float ref = ");
      printFloatC(out, fref);
      fprintf(out,";\n");
      
      fprintf(out, "fesetround(%s);\n", roundingModeString);
      fprintf(out, "float computed = %s;\n", cPrintString);
      fprintf(out, "assert(compare(ref, computed));\n");

      finishOutputC(out);

  }

  return;
}

template <unaryRoundedFunctionTFP ref>
void unaryRoundedFunctionPrintSMT (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *SMTPrintString, const char *roundingModeString) {
  uint64_t i;
  FILE * out;

  for (i = start; i < end; ++i) {
      float f = getTestValue(i);
      uint32_t input = *((uint32_t *)(&f));
     
      uint32_t reference = ref(input);
      //float fref = *((float *)(&reference));

      out = startOutputSMT(name, roundingModeString, i);

      fprintf(out, "(define-fun f () Float32 ");
      printFloatSMT(out, input);
      fprintf(out, ")\n");

      fprintf(out, "(define-fun ref () Float32 ");
      printFloatSMT(out, reference);
      fprintf(out, ")\n");

      fprintf(out, "(define-fun rm () RoundingMode %s )\n", roundingModeString);

      fprintf(out, "(define-fun result () Float32 %s )\n", SMTPrintString);
      fprintf(out, "(assert (= ref result))\n");

      finishOutputSMT(out);
  }

  return;
}




typedef bool (*unaryPredicateTFP) (uint32_t);

template <unaryPredicateTFP test, unaryPredicateTFP ref>
void unaryPredicateTest (const int verbose, const uint64_t start, const uint64_t end) {
  uint64_t i;

  for (i = start; i < end; ++i) {
      float f = getTestValue(i);
      
      uint32_t input = *((uint32_t *)(&f));
     
      bool reference = ref(input);
      bool computed = test(input);
      
      if (verbose || !(computed == reference)) {
	fprintf(stdout,"vector[%d] ", (uint32_t)i);
	fprintf(stdout,"input = 0x%x, computed = %d, real = %d\n", input, computed, reference);
	fflush(stdout);
      }

      if ((i & 0xFFFF) == 0) {
	fprintf(stdout,".");
	fflush(stdout);
      }
  }

  return;
}


template <unaryPredicateTFP ref>
void unaryPredicatePrintC (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *cPrintString, const char *) {
  uint64_t i;
  FILE *out;

  for (i = start; i < end; ++i) {
      float f = getTestValue(i);
      
      uint32_t input = *((uint32_t *)(&f));
     
      bool reference = ref(input);

      out = startOutputC(name, "NA", i);

      fprintf(out, "float f = ");
      printFloatC(out, f);
      fprintf(out,";\n");

      fprintf(out, "assert(%c(%s));\n", (reference) ? ' ' : '!', cPrintString);

      finishOutputC(out);
  }

  return;}

template <unaryPredicateTFP ref>
void unaryPredicatePrintSMT (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *SMTPrintString, const char *) {
  uint64_t i;
  FILE *out;

  for (i = start; i < end; ++i) {
      float f = getTestValue(i);
      
      uint32_t input = *((uint32_t *)(&f));
     
      bool reference = ref(input);

      out = startOutputSMT(name, "NA", i);

      fprintf(out, "(define-fun f () Float32 ");
      printFloatSMT(out, input);
      fprintf(out, ")\n");

      fprintf(out, "(define-fun ref () Bool ");
      if (reference) {
	fprintf(out, "true");
      } else {
	fprintf(out, "false");
      }
      fprintf(out, ")\n");

      fprintf(out, "(define-fun result () Bool %s )\n", SMTPrintString);

      fprintf(out, "(assert (= ref result))\n");


      finishOutputSMT(out);
  }

  return;
}


uint64_t splitRight (uint64_t input) {
  uint64_t output = 0;

  for (uint64_t i = 0; i <= 63; i +=2) {
    output |= (input & (1ULL << i)) ? 1ULL << (i >> 1ULL) : 0;
  }

  return output;
}

uint64_t splitLeft (uint64_t input) {
  uint64_t output = 0;

  for (uint64_t i = 1; i <= 63; i +=2) {
    output |= (input & (1ULL << i)) ? 1ULL << (i >> 1ULL) : 0;
  }

  return output;
}


typedef bool (*binaryPredicateTFP) (uint32_t, uint32_t);

template <binaryPredicateTFP test, binaryPredicateTFP ref>
void binaryPredicateTest (const int verbose, const uint64_t start, const uint64_t end) {
  uint64_t i;

  for (i = start; i < end; ++i) {
    uint64_t right = splitRight(i);
    uint64_t left = splitLeft(i);

    float f = getTestValue(right);
    float g = getTestValue(left);
    
    uint32_t input1 = *((uint32_t *)(&f));
    uint32_t input2 = *((uint32_t *)(&g));
    
    bool reference = ref(input1, input2);
    bool computed = test(input1, input2);
    
    if (verbose || !(computed == reference)) {
      fprintf(stdout,"vector[%d -> (%d,%d)] ", (uint32_t)i, (uint32_t)right, (uint32_t)left);
      fprintf(stdout,"input1 = 0x%x, input2 = 0x%x, computed = %d, real = %d\n", input1, input2, computed, reference);
      fflush(stdout);
    }
    
    if ((i & 0xFFFF) == 0) {
      fprintf(stdout,".");
      fflush(stdout);
    }
  }
  
  return;
}

template <binaryPredicateTFP ref>
void binaryPredicatePrintC (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *cPrintString, const char *) {
  uint64_t i;
  FILE *out;

  for (i = start; i < end; ++i) {
    uint64_t right = splitRight(i);
    uint64_t left = splitLeft(i);

    float f = getTestValue(right);
    float g = getTestValue(left);
    
    uint32_t input1 = *((uint32_t *)(&f));
    uint32_t input2 = *((uint32_t *)(&g));
    
    bool reference = ref(input1, input2);

    out = startOutputC(name, "NA", i);

    fprintf(out, "float f = ");
    printFloatC(out, f);
    fprintf(out,";\n");

    fprintf(out, "float g = ");
    printFloatC(out, g);
    fprintf(out,";\n");

    fprintf(out, "assert(%c(%s));\n", (reference) ? ' ' : '!', cPrintString);

    finishOutputC(out);
  }
  
  return;
}

template <binaryPredicateTFP ref>
void binaryPredicatePrintSMT (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *SMTPrintString, const char *) {
  uint64_t i;
  FILE *out;

  for (i = start; i < end; ++i) {
    uint64_t right = splitRight(i);
    uint64_t left = splitLeft(i);

    float f = getTestValue(right);
    float g = getTestValue(left);
    
    uint32_t input1 = *((uint32_t *)(&f));
    uint32_t input2 = *((uint32_t *)(&g));
    
    bool reference = ref(input1, input2);

    out = startOutputSMT(name, "NA", i);

    fprintf(out, "(define-fun f () Float32 ");
    printFloatSMT(out, input1);
    fprintf(out, ")\n");

    fprintf(out, "(define-fun g () Float32 ");
    printFloatSMT(out, input2);
    fprintf(out, ")\n");
    
    fprintf(out, "(define-fun ref () Bool ");
    if (reference) {
      fprintf(out, "true");
    } else {
      fprintf(out, "false");
    }
    fprintf(out, ")\n");
    
    fprintf(out, "(define-fun result () Bool %s )\n", SMTPrintString);
    
    fprintf(out, "(assert (= ref result))\n");
    


    finishOutputSMT(out);
  }
  
  return;
}








typedef uint32_t (*binaryFunctionTFP) (uint32_t, uint32_t);

template <binaryFunctionTFP test, binaryFunctionTFP ref>
void binaryFunctionTest (const int verbose, const uint64_t start, const uint64_t end) {
  uint64_t i;

  for (i = start; i < end; ++i) {
    uint64_t right = splitRight(i);
    uint64_t left = splitLeft(i);

    float f = getTestValue(right);
    float g = getTestValue(left);
    
    uint32_t input1 = *((uint32_t *)(&f));
    uint32_t input2 = *((uint32_t *)(&g));
    
    uint32_t reference = ref(input1, input2);
    uint32_t computed = test(input1, input2);

    if (verbose || !singlePrecisionHardware::smtlibEqual(computed, reference)) {
      fprintf(stdout,"vector[%d -> (%d,%d)] ", (uint32_t)i, (uint32_t)right, (uint32_t)left);
      fprintf(stdout,"input1 = 0x%x, input2 = 0x%x, computed = 0x%x, real = 0x%x\n", input1, input2, computed, reference);
      fflush(stdout);
    }

    if ((i & 0xFFFF) == 0) {
      fprintf(stdout,".");
      fflush(stdout);
    }
  }
  
  return;
}

template <binaryFunctionTFP ref>
void binaryFunctionPrintC (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *cPrintString, const char */*roundingModeString*/) {
  uint64_t i;
  FILE *out;

  for (i = start; i < end; ++i) {
    uint64_t right = splitRight(i);
    uint64_t left = splitLeft(i);

    float f = getTestValue(right);
    float g = getTestValue(left);
    
    uint32_t input1 = *((uint32_t *)(&f));
    uint32_t input2 = *((uint32_t *)(&g));
    
    uint32_t reference = ref(input1, input2);
    float fref = *((float *)(&reference));

    out = startOutputC(name, "NA", i);

    fprintf(out, "float f = ");
    printFloatC(out, f);
    fprintf(out,";\n");

    fprintf(out, "float g = ");
    printFloatC(out, g);
    fprintf(out,";\n");

    fprintf(out, "float ref = ");
    printFloatC(out, fref);
    fprintf(out,";\n");

    fprintf(out, "float computed = %s;\n", cPrintString);
    fprintf(out, "assert(compare(ref, computed));\n");

    finishOutputC(out);

  }
  
  return;
}

template <binaryFunctionTFP ref>
void binaryFunctionPrintSMT (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *SMTPrintString, const char */*roundingModeString*/) {
  uint64_t i;
  FILE *out;

  for (i = start; i < end; ++i) {
    uint64_t right = splitRight(i);
    uint64_t left = splitLeft(i);

    float f = getTestValue(right);
    float g = getTestValue(left);
    
    uint32_t input1 = *((uint32_t *)(&f));
    uint32_t input2 = *((uint32_t *)(&g));
    
    uint32_t reference = ref(input1, input2);
    //uint32_t fref = *((float *)(&reference));

    out = startOutputSMT(name, "NA", i);

    fprintf(out, "(define-fun f () Float32 ");
    printFloatSMT(out, input1);
    fprintf(out, ")\n");

    fprintf(out, "(define-fun g () Float32 ");
    printFloatSMT(out, input2);
    fprintf(out, ")\n");
    
    fprintf(out, "(define-fun ref () Float32 ");
    printFloatSMT(out, reference);
    fprintf(out, ")\n");

    fprintf(out, "(define-fun result () Float32 %s )\n", SMTPrintString);
    
    fprintf(out, "(assert (= ref result))\n");
    

    finishOutputSMT(out);

  }
  
  return;
}




typedef uint32_t (*binaryRoundedFunctionTFP) (uint32_t, uint32_t);

template <binaryRoundedFunctionTFP test, binaryRoundedFunctionTFP ref>
void binaryRoundedFunctionTest (const int verbose, const uint64_t start, const uint64_t end) {
  uint64_t i;

  for (i = start; i < end; ++i) {
    uint64_t right = splitRight(i);
    uint64_t left = splitLeft(i);

    float f = getTestValue(right);
    float g = getTestValue(left);
    
    uint32_t input1 = *((uint32_t *)(&f));
    uint32_t input2 = *((uint32_t *)(&g));
    
    uint32_t reference = ref(input1, input2);
    uint32_t computed = test(input1, input2);

    if (verbose || !singlePrecisionHardware::smtlibEqual(computed, reference)) {
      fprintf(stdout,"vector[%d -> (%d,%d)] ", (uint32_t)i, (uint32_t)right, (uint32_t)left);
      fprintf(stdout,"input1 = 0x%x, input2 = 0x%x, computed = 0x%x, real = 0x%x\n", input1, input2, computed, reference);
      fflush(stdout);
    }

    if ((i & 0xFFFF) == 0) {
      fprintf(stdout,".");
      fflush(stdout);
    }
  }
  
  return;
}

template <binaryRoundedFunctionTFP ref>
void binaryRoundedFunctionPrintC (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *cPrintString, const char *roundingModeString) {
  uint64_t i;
  FILE *out;

  for (i = start; i < end; ++i) {
    uint64_t right = splitRight(i);
    uint64_t left = splitLeft(i);

    float f = getTestValue(right);
    float g = getTestValue(left);
    
    uint32_t input1 = *((uint32_t *)(&f));
    uint32_t input2 = *((uint32_t *)(&g));
    
    uint32_t reference = ref(input1, input2);
    float fref = *((float *)(&reference));

    out = startOutputC(name, roundingModeString, i);

    fprintf(out, "float f = ");
    printFloatC(out, f);
    fprintf(out,";\n");

    fprintf(out, "float g = ");
    printFloatC(out, g);
    fprintf(out,";\n");

    fprintf(out, "float ref = ");
    printFloatC(out, fref);
    fprintf(out,";\n");

    fprintf(out, "fesetround(%s);\n", roundingModeString);
    fprintf(out, "float computed = %s;\n", cPrintString);
    fprintf(out, "assert(compare(ref, computed));\n");

    finishOutputC(out);

  }
  
  return;
}

template <binaryRoundedFunctionTFP ref>
void binaryRoundedFunctionPrintSMT (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *SMTPrintString, const char *roundingModeString) {
  uint64_t i;
  FILE *out;

  for (i = start; i < end; ++i) {
    uint64_t right = splitRight(i);
    uint64_t left = splitLeft(i);

    float f = getTestValue(right);
    float g = getTestValue(left);
    
    uint32_t input1 = *((uint32_t *)(&f));
    uint32_t input2 = *((uint32_t *)(&g));
    
    uint32_t reference = ref(input1, input2);
    //uint32_t fref = *((float *)(&reference));

    out = startOutputSMT(name, roundingModeString, i);

    fprintf(out, "(define-fun f () Float32 ");
    printFloatSMT(out, input1);
    fprintf(out, ")\n");

    fprintf(out, "(define-fun g () Float32 ");
    printFloatSMT(out, input2);
    fprintf(out, ")\n");
    
    fprintf(out, "(define-fun ref () Float32 ");
    printFloatSMT(out, reference);
    fprintf(out, ")\n");

    fprintf(out, "(define-fun rm () RoundingMode %s )\n", roundingModeString);
        
    fprintf(out, "(define-fun result () Float32 %s )\n", SMTPrintString);
    
    fprintf(out, "(assert (= ref result))\n");
    

    finishOutputSMT(out);

  }
  
  return;
}







uint64_t splitOneOfThree (uint64_t input) {
  uint64_t output = 0;

  for (uint64_t i = 0; i <= 63; i +=3) {
    output |= (input & (1ULL << i)) ? 1ULL << (i >> 1ULL) : 0;
  }

  return output;
}

uint64_t splitTwoOfThree (uint64_t input) {
  uint64_t output = 0;

  for (uint64_t i = 1; i <= 63; i +=3) {
    output |= (input & (1ULL << i)) ? 1ULL << (i >> 1ULL) : 0;
  }

  return output;
}

uint64_t splitThreeOfThree (uint64_t input) {
  uint64_t output = 0;

  for (uint64_t i = 2; i <= 63; i +=3) {
    output |= (input & (1ULL << i)) ? 1ULL << (i >> 1ULL) : 0;
  }

  return output;
}

typedef uint32_t (*ternaryRoundedFunctionTFP) (uint32_t, uint32_t, uint32_t);

template <ternaryRoundedFunctionTFP test, ternaryRoundedFunctionTFP ref>
void ternaryRoundedFunctionTest (const int verbose, const uint64_t start, const uint64_t end) {
  uint64_t i;

  for (i = start; i < end; ++i) {
    uint64_t right = splitOneOfThree(i);
    uint64_t middle = splitTwoOfThree(i);
    uint64_t left = splitThreeOfThree(i);

    float f = getTestValue(right);
    float g = getTestValue(middle);
    float h = getTestValue(left);
    
    uint32_t input1 = *((uint32_t *)(&f));
    uint32_t input2 = *((uint32_t *)(&g));
    uint32_t input3 = *((uint32_t *)(&h));
    
    uint32_t reference = ref(input1, input2, input3);
    uint32_t computed = test(input1, input2, input3);

    if (verbose || !singlePrecisionHardware::smtlibEqual(computed, reference)) {
      fprintf(stdout,"vector[%d -> (%d,%d,%d)] ", (uint32_t)i, (uint32_t)right, (uint32_t)middle, (uint32_t)left);
      fprintf(stdout,"input1 = 0x%x, input2 = 0x%x, input3 = 0x%x, computed = 0x%x, real = 0x%x\n", input1, input2, input3, computed, reference);
      fflush(stdout);
    }

    if ((i & 0xFFFF) == 0) {
      fprintf(stdout,".");
      fflush(stdout);
    }
  }
  
  return;
}

template <ternaryRoundedFunctionTFP ref>
void ternaryRoundedFunctionPrintC (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *cPrintString, const char *roundingModeString) {
  uint64_t i;
  FILE *out;

  for (i = start; i < end; ++i) {
    uint64_t right = splitOneOfThree(i);
    uint64_t middle = splitTwoOfThree(i);
    uint64_t left = splitThreeOfThree(i);

    float f = getTestValue(right);
    float g = getTestValue(middle);
    float h = getTestValue(left);
    
    uint32_t input1 = *((uint32_t *)(&f));
    uint32_t input2 = *((uint32_t *)(&g));
    uint32_t input3 = *((uint32_t *)(&h));
    
    uint32_t reference = ref(input1, input2, input3);
    float fref = *((float *)(&reference));

    out = startOutputC(name, roundingModeString, i);

    fprintf(out, "float f = ");
    printFloatC(out, f);
    fprintf(out,";\n");

    fprintf(out, "float g = ");
    printFloatC(out, g);
    fprintf(out,";\n");

    fprintf(out, "float h = ");
    printFloatC(out, g);
    fprintf(out,";\n");

    fprintf(out, "float ref = ");
    printFloatC(out, fref);
    fprintf(out,";\n");

    fprintf(out, "fesetround(%s);\n", roundingModeString);
    fprintf(out, "float computed = %s;\n", cPrintString);
    fprintf(out, "assert(compare(ref, computed));\n");

    finishOutputC(out);
  }
  
  return;
}

template <ternaryRoundedFunctionTFP ref>
void ternaryRoundedFunctionPrintSMT (const int /*verbose*/, const uint64_t start, const uint64_t end, const char *name, const char *SMTPrintString, const char *roundingModeString) {
  uint64_t i;
  FILE *out;

  for (i = start; i < end; ++i) {
    uint64_t right = splitOneOfThree(i);
    uint64_t middle = splitTwoOfThree(i);
    uint64_t left = splitThreeOfThree(i);

    float f = getTestValue(right);
    float g = getTestValue(middle);
    float h = getTestValue(left);
    
    uint32_t input1 = *((uint32_t *)(&f));
    uint32_t input2 = *((uint32_t *)(&g));
    uint32_t input3 = *((uint32_t *)(&h));
    
    uint32_t reference = ref(input1, input2, input3);
    // float fref = *((float *)(&reference));

    out = startOutputSMT(name, roundingModeString, i);

    fprintf(out, "(define-fun f () Float32 ");
    printFloatSMT(out, input1);
    fprintf(out, ")\n");

    fprintf(out, "(define-fun g () Float32 ");
    printFloatSMT(out, input2);
    fprintf(out, ")\n");

    fprintf(out, "(define-fun h () Float32 ");
    printFloatSMT(out, input2);
    fprintf(out, ")\n");

    fprintf(out, "(define-fun ref () Float32 ");
    printFloatSMT(out, reference);
    fprintf(out, ")\n");

    fprintf(out, "(define-fun rm () RoundingMode %s )\n", roundingModeString);
        
    fprintf(out, "(define-fun result () Float32 %s )\n", SMTPrintString);
    
    fprintf(out, "(assert (= ref result))\n");


    finishOutputSMT(out);
  }
  
  return;
}



typedef void (*testFunction) (const int, const uint64_t, const uint64_t);
typedef void (*printFunction) (const int, const uint64_t, const uint64_t, const char *, const char *, const char *);

struct testStruct {
  int enable;
  const int usesRounding;
  const char *name;
  const testFunction run;
  const printFunction printC;
  const printFunction printSMT;
  const char *cPrintString;
  const char *SMTPrintString;
};

struct roundingModeTestStruct {
  int enable;
  const char *name;
  const int value;
  const char *cPrintString;
};





/*** Application ***/



#define INCREMENT 0xFFFFFF

#define TEST 0
#define PRINTC 1
#define PRINTSMT 2

// Save on typing!
#define INST(T,F) T##Test<singlePrecisionExecutableSymfpu::F, singlePrecisionHardware::F>, T##PrintC<singlePrecisionHardware::F>, T##PrintSMT<singlePrecisionHardware::F>

int main (int argc, char **argv) {
  struct testStruct tests[] = {
    {0,0,         "unpackPack", INST( unaryFunction, unpackPack),      "f",        "f"},
    {0,0,             "negate", INST( unaryFunction, negate),          "-f",       "(fp.neg f)"},
    {0,0,           "absolute", INST( unaryFunction, absolute),        "fabsf(f)", "(fp.abs f)"},
    {0,0,           "isNormal", INST(unaryPredicate, isNormal),        "isnormal(f)", "(fp.isNormal f)"},
    {0,0,        "isSubnormal", INST(unaryPredicate, isSubnormal),     "fpclassify(f) == FP_SUBNORMAL", "(fp.isSubnormal f)"},
    {0,0,             "isZero", INST(unaryPredicate, isZero),          "(f) == 0.0f", "(fp.isZero f)"},
    {0,0,         "isInfinite", INST(unaryPredicate, isInfinite),      "isinf(f)", "(fp.isInfinite f)"},
    {0,0,              "isNaN", INST(unaryPredicate, isNaN),           "isnan(f)", "(fp.isNaN f)"},
    {0,0,         "isPositive", INST(unaryPredicate, isPositive),      "!isnan(f) && signbit(f) == 0", "(fp.isPositive f)"},
    {0,0,         "isNegative", INST(unaryPredicate, isNegative),      "!isnan(f) && signbit(f) != 0", "(fp.isNegative f)"},
    {0,0,      "SMT-LIB_equal", INST(binaryPredicate, smtlibEqual),    "compare(f,g)", "(= f g)"},
    {0,0,       "IEE754_equal", INST(binaryPredicate, ieee754Equal),    "f == g",  "(fp.eq f g)"},
    {0,0,          "less_than", INST(binaryPredicate, lessThan),         "f < g",  "(fp.lt f g)"},
    {0,0, "less_than_or_equal", INST(binaryPredicate, lessThanOrEqual), "f <= g",  "(fp.leq f g)"},
    {0,1,           "multiply", INST(binaryRoundedFunction, multiply),  "f * g",  "(fp.mul rm f g)"},
    {0,1,                "add", INST(binaryRoundedFunction, add),       "f + g",  "(fp.add rm f g)"},
    {0,1,           "subtract", INST(binaryRoundedFunction, sub),       "f - g",  "(fp.sub rm f g)"},
    {0,1,             "divide", INST(binaryRoundedFunction, div),       "f / g",  "(fp.div rm f g)"},
    {0,0,                "max", INST(binaryFunction, max),              "fmaxf(f,g)",  "(fp.max f g)"},
    {0,0,                "min", INST(binaryFunction, min),              "fminf(f,g)",  "(fp.min f g)"},
    {0,1,               "sqrt", INST(unaryRoundedFunction, sqrt),       "sqrtf(f)",  "(fp.sqrt rm f)"},
    {0,1,  "round_to_integral", INST(unaryRoundedFunction, rti),        "(fegetround()==FE_TONEAREST) ? rintf(f) : (fegetround()==FE_UPWARD) ? ceilf(f) : (fegetround()==FE_DOWNWARD) ? floorf(f) : truncf(f)",  "(fp.roundToIntegral rm f)"},
    {0,1,                "fma", INST(ternaryRoundedFunction, fma),      "fmaf(f,g)",  "(fp.fma rm f g h)"},
    {0,0,          "remainder", INST(binaryFunction, rem),              "remainderf(f,g)",  "(fp.remainder f g)"},
    {0,0,                 NULL, NULL, NULL, NULL,                           NULL,  NULL}
  };

  struct roundingModeTestStruct roundingModeTests[] = {
    {0, "RNE",  FE_TONEAREST, "FE_TONEAREST"},
    {0, "RTP",     FE_UPWARD, "FE_UPWARD"},
    {0, "RTN",   FE_DOWNWARD, "FE_DOWNWARD"},
    {0, "RTZ", FE_TOWARDZERO, "FE_TOWARDZERO"},
  //{0, "RNA",        FE_???, "FE_???"},   // Disabled until a suitable reference is available
    {0,  NULL,             0, NULL}
  };

  uint64_t start = 0;
  uint64_t end = INCREMENT;
  int verbose = 0;
  int help = 0;
  int enableAllTests = 0;
  int enableAllRoundingModes = 0;
  int continuous = 0;
  int action = TEST;

  struct option options[] = {
    {         "verbose",        no_argument,                          &verbose,  1 },
    {            "help",        no_argument,                             &help,  1 },

    {           "start",  required_argument,                              NULL, 's'},
    {             "end",  required_argument,                              NULL, 'e'},
    {   "specialValues",        no_argument,                              NULL, 't'},
    {      "continuous",        no_argument,                       &continuous,  1 },

    {        "allTests",        no_argument,                   &enableAllTests,  1 },
    {"allRoundingModes",        no_argument,           &enableAllRoundingModes,  1 },

    {          "printC",        no_argument,                           &action,  PRINTC },
    {        "printSMT",        no_argument,                           &action,  PRINTSMT },

    {      "unpackPack",        no_argument,                &(tests[0].enable),  1 },
    {          "negate",        no_argument,                &(tests[1].enable),  1 },
    {        "absolute",        no_argument,                &(tests[2].enable),  1 },
    {        "isNormal",        no_argument,                &(tests[3].enable),  1 },
    {     "isSubnormal",        no_argument,                &(tests[4].enable),  1 },
    {          "isZero",        no_argument,                &(tests[5].enable),  1 },
    {      "isInfinite",        no_argument,                &(tests[6].enable),  1 },
    {           "isNaN",        no_argument,                &(tests[7].enable),  1 },
    {      "isPositive",        no_argument,                &(tests[8].enable),  1 },
    {      "isNegative",        no_argument,                &(tests[9].enable),  1 },
    {     "smtlibEqual",        no_argument,               &(tests[10].enable),  1 },
    {    "ieee754Equal",        no_argument,               &(tests[11].enable),  1 },
    {        "lessThan",        no_argument,               &(tests[12].enable),  1 },
    { "lessThanOrEqual",        no_argument,               &(tests[13].enable),  1 },
    {        "multiply",        no_argument,               &(tests[14].enable),  1 },
    {             "add",        no_argument,               &(tests[15].enable),  1 },
    {        "subtract",        no_argument,               &(tests[16].enable),  1 },
    {          "divide",        no_argument,               &(tests[17].enable),  1 },
    {             "max",        no_argument,               &(tests[18].enable),  1 },
    {             "min",        no_argument,               &(tests[19].enable),  1 },
    {            "sqrt",        no_argument,               &(tests[20].enable),  1 },
    {             "rti",        no_argument,               &(tests[21].enable),  1 },
    {             "fma",        no_argument,               &(tests[22].enable),  1 },
    {       "remainder",        no_argument,               &(tests[23].enable),  1 },
    {             "rne",        no_argument,    &(roundingModeTests[0].enable),  1 },
    {             "rtp",        no_argument,    &(roundingModeTests[1].enable),  1 },
    {             "rtn",        no_argument,    &(roundingModeTests[2].enable),  1 },
    {             "rtz",        no_argument,    &(roundingModeTests[3].enable),  1 },
  //{             "rna",        no_argument,    &(roundingModeTests[4].enable),  1 },
    {             "RNE",        no_argument,    &(roundingModeTests[0].enable),  1 },
    {             "RTP",        no_argument,    &(roundingModeTests[1].enable),  1 },
    {             "RTN",        no_argument,    &(roundingModeTests[2].enable),  1 },
    {             "RTZ",        no_argument,    &(roundingModeTests[3].enable),  1 },
  //{             "RNA",        no_argument,    &(roundingModeTests[4].enable),  1 },
    {              NULL,                  0,                              NULL,  0 }


  }; 

  bool parseOptions = true;
  while (parseOptions) {
    int currentOption = 0;
    int response = getopt_long(argc, argv, "vs:e:t", options, &currentOption);

    switch(response) {
    case 'v' :
      verbose = 1;
      break;

    case 's' :
      start = strtoull(optarg,NULL,0);
      break;

    case 'e' :
      end = strtoull(optarg,NULL,0);
      break;

    case 't' :
      end = NUMBER_OF_FLOAT_TESTS  * NUMBER_OF_FLOAT_TESTS; // TODO : split doesn't work like this so this isn't quite right
      break;

    case 0 :    /* Flag set */
      break;

    case -1 :   /* End of options */
      parseOptions = false;
      break;

    case '?' :  /* Unknown option */
      fprintf(stderr,"Unknown option : \"%s\"\n", optarg);
      return 1;
      break;

    default :   /* Unrecognised getopt response */
      fprintf(stderr,"Unrecognised getopt response \'%d\' from \"%s\"\n",response, optarg);
      return 1;
      break;
    }
  }

  /* Help text */
  if (help) {
    fprintf(stderr, "TODO: Super helpful message here!\n");

    fprintf(stderr, "Options : \n");
    struct option *opt = &(options[0]);
    while (opt->name != NULL) {
      switch (opt->has_arg) {
      case no_argument : fprintf(stderr, "\t--%s\n", opt->name); break;
      case required_argument : fprintf(stderr, "\t--%s\t  argument\n", opt->name); break;
      case optional_argument : fprintf(stderr, "\t--%s\t[ argument ]\n", opt->name); break;
      default :	fprintf(stderr, "\t--%s\t???\n", opt->name); break;
      }

      ++opt;
    }

    return 0;
  }

  int roundingModeSet = 0;
  uint64_t j = 0;
  while (roundingModeTests[j].name != NULL) {
    roundingModeSet |= roundingModeTests[j].enable;
    ++j;
  }
  if (!roundingModeSet && !enableAllRoundingModes) {
    // Set default rounding mode of RNE
    roundingModeTests[0].enable = 1;
  }


  singlePrecisionExecutableSymfpu::setFormat(singlePrecisionFormatObject);
  
 top :

  /* Run the tests */
  int i = 0;
  while (tests[i].name != NULL) {
    if (enableAllTests || tests[i].enable) {
      if (tests[i].usesRounding) {
	
	int j = 0;
	while (roundingModeTests[j].name != NULL) {
	  if (enableAllRoundingModes || roundingModeTests[j].enable) {
	    fprintf(stdout, "Running test for %s %s : ", tests[i].name, roundingModeTests[j].name);
	    fflush(stdout);

	    singlePrecisionExecutableSymfpu::setRoundingMode(roundingModeTests[j].value);
	    singlePrecisionHardware::setRoundingMode(roundingModeTests[j].value);
	    
	    switch (action) {
	    case TEST :
	      tests[i].run(verbose, start, end);
	      break;
	      
	    case PRINTC :
	      tests[i].printC(verbose, start, end, tests[i].name, tests[i].cPrintString, roundingModeTests[j].cPrintString);
	      break;
	      
	    case PRINTSMT :
	      tests[i].printSMT(verbose, start, end, tests[i].name, tests[i].SMTPrintString, roundingModeTests[j].name);
	      break;
	      
	    default :
	      assert(0);
	      break;
	    }

	    fprintf(stdout, "\n");
	    fflush(stdout);
	  }
	  ++j;
	}


      } else {
	fprintf(stdout, "Running test for %s : ", tests[i].name);
	fflush(stdout);

	switch (action) {
	case TEST :
	  tests[i].run(verbose, start, end);
	  break;
	  
	case PRINTC :
	  tests[i].printC(verbose, start, end, tests[i].name, tests[i].cPrintString, NULL);
	  break;
	  
	case PRINTSMT :
	  tests[i].printSMT(verbose, start, end, tests[i].name, tests[i].SMTPrintString, NULL);
	  break;
	  
	default :
	  assert(0);
	  break;
	}

	fprintf(stdout, "\n");
	fflush(stdout);

      }
      
    }
    i++;

  }

  if (continuous) {
    uint64_t oldEnd = end;
    end += (end - start);
    start = oldEnd;
    goto top;
  }

  singlePrecisionExecutableSymfpu::destroyFormat();

  return 1;
}
