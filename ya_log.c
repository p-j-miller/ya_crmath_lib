/* ya_log.c
  ======

   Uses basically the same algorithms as the log function in llvm-project-23.1.1 (\llvm-project-23.1.1\libc\src\__support\FPUtil & \llvm-project-23.1.1\libc\src\__support\math ).
   Note the llvm original is in C++, whereas this is in C.
   Where practical, variable names (and many comments) have been kept the same in this code as the llvm code, to make them easier to compare.
   It leverages the double-double maths library ya_double_double.   
   
   This code assumes double FMA is available (either as a single instruction [ideally], or emulated in software) 
 
   The "accurate" evaluation is done using using functions in cr_log() - in the test program these are only used 0.001% of the time so it was not worthwhile converting the accurate algorithms
   used in the llvm implementations(which might or might not be faster, but will make no difference to the average timing and negligible difference to the worse case).

   Checking for special cases (inf, NAN, etc) is done in the same was as cr_log.c, as that solution was found to be faster than the approach used in llvm (& supports CORE_MATH_SUPPORT_ERRNO).
   
 	According to  https://members.loria.fr/PZimmermann/papers/accuracy.pdf the llvm double log function has a maximum error of 0.500ulp,
 	and this confirmed at https://libc.llvm.org/headers/math/index.html#higher-math-functions (<Func> (double) column).
	This implementation also achieves this as confirmed by the test program in this directory (main.c).
	Note that this test program only checks the default rounding mode (round to even), the other rounding modes have not been tested.  
	
   This C code was created by Peter Miller 17-9-2026   
*/
/*
This is licensed under the following standard MIT license:

----------------------------------------------------------------------
Copyright © 2026 Peter Miller.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
----------------------------------------------------------------------
*/ 
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <inttypes.h>
#include <ctype.h>
#define YA_DD_LO_FIRST	 /* we want lo first in DoubleDouble to match llm ordering (required for initialisations) */
#include "../ya_double_double/ya_double_double.h" // double-double routines here assume we have FMA - note this will be implemented in software if not available as an instruction in the processor (available for AVX2)
#include "ya_crmath_config.h" /* get configuration - like CORE_MATH_SUPPORT_ERRNO */
#ifdef CORE_MATH_SUPPORT_ERRNO
 #include <errno.h>
#endif

#include "ya_crmath.h"

/* code below cannot be compiled with -Ofast as this makes the compiler break some C rules that we need, so make sure of this here */
/* we also need -msse2 and -mfpmath=sse to actually use the sse instructions for float and double maths */
/* there seems to be no way to duplicate "-fexcess-precision=standard" using a pragma - so that must be present on the command line [see comments at head of this file that suggest "-fexcess-precision=standard" is not required any more ] */
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)) || defined(__clang__)
 #if !defined(__BORLANDC__)
  #pragma GCC push_options
  #pragma GCC optimize ("-O3") /* cannot use Ofast, normally -O3 is OK. Note macro expansion does not work here ! */
  #pragma GCC target("avx2") /*  -mavx2 [for fma] */
 #endif
 // based on  https://jdebp.uk/FGA/predefined-macros-processor.html "__i386__" is set by GCC,Clang,Intel which is good enough as the outer #if limits us to gcc and clang
 #ifdef __i386__
   #pragma GCC target("sse2,fpmath=sse,avx2,fma") /* -msse2 and -mfpmath=sse, -mavx2 [the later for fma] */
 #endif 
#endif

// PMi code
typedef union {double f; uint64_t u;} b64u64_u;
// characteristics of ieee format doubles
#define DOUBLE_MANTISSA_BITS 52
#define DOUBLE_EXPONENT_BITS 11
#define DOUBLE_EXPONENT_BIAS 1023
#define DOUBLE_IMPLICIT_BIT ((uint64_t)1<< DOUBLE_MANTISSA_BITS) 
#ifdef __GNUC__  /* gcc or clang - not events most be very common/uncommon for this to give performance gains, > 1 in 10,000 is recommended. Use in if statements */
 #define ya_LIKELY(b) __builtin_expect((!!(b)),1)
 #define ya_UNLIKELY(b) __builtin_expect((!!(b)),0)
#else
 #define ya_LIKELY(b) (!!(b)) /* (!!(b)) appears to be faster than just (b) with gcc 16.1.0 and AMD Ryzen 5 processor */
 #define ya_UNLIKELY(b) (!!(b))
#endif

#define multiply_add(a,b,c) ( ((a)*(b))+(c) ) /* this does not need to be a "real" FMA from an accuracy viewpoint , but compiler can use FMA instruction if its faster */

// log.h

// A simple upper bound for the error of e_x * log(2) - log(r).
static const double HI_ERR = 0x1.0p-85;
// log(2) generated by Sollya with:
// > a = 2^-43 * nearestint(2^43*log(2));
// LSB = 2^-43 is chosen so that e_x * LOG_2_HI is exact for -1075 < e_x < 1024.
static const double LOG_2_HI = 0x1.62e42fefa38p-1; // LSB = 2^-43
// > b = round(log10(2) - a, D, RN);
static const double LOG_2_LO = 0x1.ef35793c7673p-45; // LSB = 2^-97

// Extra errors from P is from using x^2 to reduce evaluation latency.
static const double P_ERR = 0x1.0p-50;

// Range reduction constants for logarithms.
// r(0) = 1, r(127) = 0.5
// r(k) = 2^-8 * ceil(2^8 * (1 - 2^-8) / (1 + k*2^-7))
// The constants are chosen so that v = fma(r, m_x, -1) is exact in single
// precision, and -2^-8 <= v < 2^-7.
static const double RD[128] = {
    0x1p0,     0x1.fcp-1, 0x1.f8p-1, 0x1.f4p-1, 0x1.fp-1,  0x1.ecp-1, 0x1.e8p-1,
    0x1.e4p-1, 0x1.ep-1,  0x1.dep-1, 0x1.dap-1, 0x1.d6p-1, 0x1.d4p-1, 0x1.dp-1,
    0x1.ccp-1, 0x1.cap-1, 0x1.c6p-1, 0x1.c4p-1, 0x1.cp-1,  0x1.bep-1, 0x1.bap-1,
    0x1.b8p-1, 0x1.b4p-1, 0x1.b2p-1, 0x1.aep-1, 0x1.acp-1, 0x1.a8p-1, 0x1.a6p-1,
    0x1.a4p-1, 0x1.ap-1,  0x1.9ep-1, 0x1.9cp-1, 0x1.98p-1, 0x1.96p-1, 0x1.94p-1,
    0x1.92p-1, 0x1.9p-1,  0x1.8cp-1, 0x1.8ap-1, 0x1.88p-1, 0x1.86p-1, 0x1.84p-1,
    0x1.8p-1,  0x1.7ep-1, 0x1.7cp-1, 0x1.7ap-1, 0x1.78p-1, 0x1.76p-1, 0x1.74p-1,
    0x1.72p-1, 0x1.7p-1,  0x1.6ep-1, 0x1.6cp-1, 0x1.6ap-1, 0x1.68p-1, 0x1.66p-1,
    0x1.64p-1, 0x1.62p-1, 0x1.6p-1,  0x1.5ep-1, 0x1.5cp-1, 0x1.5ap-1, 0x1.58p-1,
    0x1.56p-1, 0x1.54p-1, 0x1.54p-1, 0x1.52p-1, 0x1.5p-1,  0x1.4ep-1, 0x1.4cp-1,
    0x1.4ap-1, 0x1.4ap-1, 0x1.48p-1, 0x1.46p-1, 0x1.44p-1, 0x1.42p-1, 0x1.4p-1,
    0x1.4p-1,  0x1.3ep-1, 0x1.3cp-1, 0x1.3ap-1, 0x1.3ap-1, 0x1.38p-1, 0x1.36p-1,
    0x1.34p-1, 0x1.34p-1, 0x1.32p-1, 0x1.3p-1,  0x1.3p-1,  0x1.2ep-1, 0x1.2cp-1,
    0x1.2cp-1, 0x1.2ap-1, 0x1.28p-1, 0x1.28p-1, 0x1.26p-1, 0x1.24p-1, 0x1.24p-1,
    0x1.22p-1, 0x1.2p-1,  0x1.2p-1,  0x1.1ep-1, 0x1.1cp-1, 0x1.1cp-1, 0x1.1ap-1,
    0x1.1ap-1, 0x1.18p-1, 0x1.16p-1, 0x1.16p-1, 0x1.14p-1, 0x1.14p-1, 0x1.12p-1,
    0x1.1p-1,  0x1.1p-1,  0x1.0ep-1, 0x1.0ep-1, 0x1.0cp-1, 0x1.0cp-1, 0x1.0ap-1,
    0x1.0ap-1, 0x1.08p-1, 0x1.08p-1, 0x1.06p-1, 0x1.06p-1, 0x1.04p-1, 0x1.04p-1,
    0x1.02p-1, 0x1.0p-1};

// Generated by Sollya with:
// for i from 0 to 127 do {
//     r = 2^-8 * ceil( 2^8 * (1 - 2^(-8)) / (1 + i*2^-7) );
//     b = nearestint(log(r)*2^43) * 2^-43;
//     c = round(log(r) - b, D, RN);
//     print("{", -c, ",", -b, "},");
//   };
// We replace LOG_R[0] with log10(1.0) == 0.0
static const DoubleDouble LOG_R_DD[128] = {
    {0.0, 0.0},
    {-0x1.0c76b999d2be8p-46, 0x1.010157589p-7},
    {-0x1.3dc5b06e2f7d2p-45, 0x1.0205658938p-6},
    {-0x1.aa0ba325a0c34p-45, 0x1.8492528c9p-6},
    {0x1.111c05cf1d753p-47, 0x1.0415d89e74p-5},
    {-0x1.c167375bdfd28p-45, 0x1.466aed42ep-5},
    {-0x1.97995d05a267dp-46, 0x1.894aa149fcp-5},
    {-0x1.a68f247d82807p-46, 0x1.ccb73cdddcp-5},
    {-0x1.7e5dd7009902cp-46, 0x1.08598b59e4p-4},
    {-0x1.5325d560d9e9bp-45, 0x1.1973bd1466p-4},
    {0x1.cc85ea5db4ed7p-45, 0x1.3bdf5a7d1ep-4},
    {-0x1.c69063c5d1d1ep-45, 0x1.5e95a4d97ap-4},
    {0x1.c1e8da99ded32p-49, 0x1.700d30aeacp-4},
    {0x1.3115c3abd47dap-45, 0x1.9335e5d594p-4},
    {-0x1.390802bf768e5p-46, 0x1.b6ac88dad6p-4},
    {0x1.646d1c65aacd3p-45, 0x1.c885801bc4p-4},
    {-0x1.dc068afe645ep-45, 0x1.ec739830a2p-4},
    {-0x1.534d64fa10afdp-45, 0x1.fe89139dbep-4},
    {0x1.1ef78ce2d07f2p-45, 0x1.1178e8227ep-3},
    {0x1.ca78e44389934p-45, 0x1.1aa2b7e23fp-3},
    {0x1.39d6ccb81b4a1p-47, 0x1.2d1610c868p-3},
    {0x1.62fa8234b7289p-51, 0x1.365fcb0159p-3},
    {0x1.5837954fdb678p-45, 0x1.4913d8333bp-3},
    {0x1.633e8e5697dc7p-45, 0x1.527e5e4a1bp-3},
    {0x1.9cf8b2c3c2e78p-46, 0x1.6574ebe8c1p-3},
    {-0x1.5118de59c21e1p-45, 0x1.6f0128b757p-3},
    {0x1.e0ddb9a631e83p-46, 0x1.823c16551ap-3},
    {-0x1.73d54aae92cd1p-47, 0x1.8beafeb39p-3},
    {0x1.7f22858a0ff6fp-47, 0x1.95a5adcf7p-3},
    {-0x1.8724350562169p-45, 0x1.a93ed3c8aep-3},
    {-0x1.c358d4eace1aap-47, 0x1.b31d8575bdp-3},
    {-0x1.d4bc4595412b6p-45, 0x1.bd087383bep-3},
    {-0x1.84a7e75b6f6e4p-47, 0x1.d1037f2656p-3},
    {-0x1.aff2af715b035p-45, 0x1.db13db0d49p-3},
    {0x1.212276041f43p-51, 0x1.e530effe71p-3},
    {-0x1.a211565bb8e11p-51, 0x1.ef5ade4ddp-3},
    {0x1.bcbecca0cdf3p-46, 0x1.f991c6cb3bp-3},
    {0x1.89cdb16ed4e91p-48, 0x1.07138604d58p-2},
    {0x1.7188b163ceae9p-45, 0x1.0c42d67616p-2},
    {-0x1.c210e63a5f01cp-45, 0x1.1178e8227e8p-2},
    {0x1.b9acdf7a51681p-45, 0x1.16b5ccbacf8p-2},
    {0x1.ca6ed5147bdb7p-45, 0x1.1bf99635a68p-2},
    {0x1.c93c1df5bb3b6p-45, 0x1.269621134d8p-2},
    {0x1.a9cfa4a5004f4p-45, 0x1.2bef07cdc9p-2},
    {-0x1.8e27ad3213cb8p-45, 0x1.314f1e1d36p-2},
    {0x1.16ecdb0f177c8p-46, 0x1.36b6776be1p-2},
    {0x1.83b54b606bd5cp-46, 0x1.3c25277333p-2},
    {0x1.8e436ec90e09dp-47, 0x1.419b423d5e8p-2},
    {-0x1.f27ce0967d675p-45, 0x1.4718dc271c8p-2},
    {-0x1.e20891b0ad8a4p-45, 0x1.4c9e09e173p-2},
    {0x1.ebe708164c759p-45, 0x1.522ae0738ap-2},
    {0x1.fadedee5d40efp-46, 0x1.57bf753c8dp-2},
    {-0x1.a0b2a08a465dcp-47, 0x1.5d5bddf596p-2},
    {-0x1.db623e731aep-45, 0x1.630030b3abp-2},
    {0x1.0a0d32756ebap-45, 0x1.68ac83e9c68p-2},
    {0x1.721657c222d87p-46, 0x1.6e60ee6af18p-2},
    {0x1.d8b0949dc60b3p-45, 0x1.741d876c678p-2},
    {0x1.9ec7d2efd1778p-45, 0x1.79e26687cf8p-2},
    {-0x1.72090c812566ap-45, 0x1.7fafa3bd818p-2},
    {0x1.fd56f3333778ap-45, 0x1.85855776dc8p-2},
    {-0x1.05ae1e5e7047p-45, 0x1.8b639a88b3p-2},
    {-0x1.766b52ee6307dp-46, 0x1.914a8635bf8p-2},
    {-0x1.52313a502d9fp-46, 0x1.973a3431358p-2},
    {-0x1.6279e10d0c0bp-45, 0x1.9d32bea15fp-2},
    {0x1.3c6457f9d79f5p-45, 0x1.a33440224f8p-2},
    {0x1.3c6457f9d79f5p-45, 0x1.a33440224f8p-2},
    {0x1.e36f2bea77a5dp-46, 0x1.a93ed3c8ad8p-2},
    {-0x1.17cc552774458p-45, 0x1.af5295248dp-2},
    {0x1.095252d841995p-46, 0x1.b56fa044628p-2},
    {0x1.7d85bf40a666dp-45, 0x1.bb9611b80ep-2},
    {0x1.cec807fe8e18p-45, 0x1.c1c60693fap-2},
    {0x1.cec807fe8e18p-45, 0x1.c1c60693fap-2},
    {-0x1.9b6ddc15249aep-45, 0x1.c7ff9c74558p-2},
    {-0x1.797c33ec7a6bp-47, 0x1.ce42f180648p-2},
    {0x1.35bafe9a767a8p-45, 0x1.d490246def8p-2},
    {-0x1.ea42d60dc616ap-46, 0x1.dae75484c98p-2},
    {-0x1.326b207322938p-46, 0x1.e148a1a2728p-2},
    {-0x1.326b207322938p-46, 0x1.e148a1a2728p-2},
    {-0x1.465505372bd08p-45, 0x1.e7b42c3ddbp-2},
    {0x1.f27f45a470251p-45, 0x1.ee2a156b41p-2},
    {0x1.2cde56f014a8bp-46, 0x1.f4aa7ee0318p-2},
    {0x1.2cde56f014a8bp-46, 0x1.f4aa7ee0318p-2},
    {0x1.085fa3c164935p-47, 0x1.fb358af7a48p-2},
    {-0x1.53ba3b1727b1cp-47, 0x1.00e5ae5b208p-1},
    {-0x1.4c45fe79539ep-47, 0x1.04360be7604p-1},
    {-0x1.4c45fe79539ep-47, 0x1.04360be7604p-1},
    {0x1.6812241edf5fdp-45, 0x1.078bf0533c4p-1},
    {0x1.f486b887e7e27p-46, 0x1.0ae76e2d054p-1},
    {0x1.f486b887e7e27p-46, 0x1.0ae76e2d054p-1},
    {0x1.c299807801742p-46, 0x1.0e4898611ccp-1},
    {-0x1.58647bb9ddcb2p-45, 0x1.11af823c75cp-1},
    {-0x1.58647bb9ddcb2p-45, 0x1.11af823c75cp-1},
    {-0x1.edd97a293ae49p-45, 0x1.151c3f6f298p-1},
    {0x1.4cc4ef8ab465p-46, 0x1.188ee40f23cp-1},
    {0x1.4cc4ef8ab465p-46, 0x1.188ee40f23cp-1},
    {0x1.cacdeed70e667p-51, 0x1.1c07849ae6p-1},
    {-0x1.a7242c9fe81d3p-45, 0x1.1f8635fc618p-1},
    {-0x1.a7242c9fe81d3p-45, 0x1.1f8635fc618p-1},
    {0x1.2fc066e48667bp-46, 0x1.230b0d8bebcp-1},
    {-0x1.b61f10522625p-47, 0x1.269621134dcp-1},
    {-0x1.b61f10522625p-47, 0x1.269621134dcp-1},
    {0x1.06d2be797882dp-45, 0x1.2a2786d0ecp-1},
    {-0x1.7a6e507b9dc11p-46, 0x1.2dbf557b0ep-1},
    {-0x1.7a6e507b9dc11p-46, 0x1.2dbf557b0ep-1},
    {-0x1.74e93c5a0ed9cp-45, 0x1.315da443408p-1},
    {-0x1.74e93c5a0ed9cp-45, 0x1.315da443408p-1},
    {0x1.0b83f9527e6acp-46, 0x1.35028ad9d8cp-1},
    {-0x1.18b7abb5569a4p-45, 0x1.38ae2171978p-1},
    {-0x1.18b7abb5569a4p-45, 0x1.38ae2171978p-1},
    {-0x1.2b7367cfe13c2p-47, 0x1.3c6080c36cp-1},
    {-0x1.2b7367cfe13c2p-47, 0x1.3c6080c36cp-1},
    {-0x1.6ce7930f0c74cp-45, 0x1.4019c2125ccp-1},
    {-0x1.d984f481051f7p-48, 0x1.43d9ff2f924p-1},
    {-0x1.d984f481051f7p-48, 0x1.43d9ff2f924p-1},
    {-0x1.2cb6af94d60aap-45, 0x1.47a1527e8a4p-1},
    {-0x1.2cb6af94d60aap-45, 0x1.47a1527e8a4p-1},
    {0x1.f7115ed4c541cp-49, 0x1.4b6fd6f970cp-1},
    {0x1.f7115ed4c541cp-49, 0x1.4b6fd6f970cp-1},
    {-0x1.e6c516d93b8fbp-45, 0x1.4f45a835a5p-1},
    {-0x1.e6c516d93b8fbp-45, 0x1.4f45a835a5p-1},
    {0x1.5ccc45d257531p-47, 0x1.5322e268678p-1},
    {0x1.5ccc45d257531p-47, 0x1.5322e268678p-1},
    {0x1.9980bff3303ddp-47, 0x1.5707a26bb8cp-1},
    {0x1.9980bff3303ddp-47, 0x1.5707a26bb8cp-1},
    {0x1.dfa63ac10c9fbp-45, 0x1.5af405c3648p-1},
    {0x1.dfa63ac10c9fbp-45, 0x1.5af405c3648p-1},
    {0x1.202380cda46bep-45, 0x1.5ee82aa2418p-1},
    {0.0, 0.0},
};

// Minimax polynomial for (log(1 + x) - x)/x^2, generated by sollya with:
// > P = fpminimax((log(1 + x) - x)/x^2, 5, [|D...|], [-2^-8, 2^-7]);
static const double LOG_COEFFS[6] = {
    -0x1.fffffffffffffp-2, 0x1.5555555554a9bp-2,  -0x1.0000000094567p-2,
    0x1.99999dcc9823cp-3,  -0x1.55550ac2e537ap-3, 0x1.21a02c4e624d7p-3};
    

double ya_log(double x) 
{
  int x_e = -DOUBLE_EXPONENT_BIAS;
 // handle special cases as per cr_log()
  b64u64_u v = {.f = x};
  uint64_t x_u = v.u;
  if(x_u==(uint64_t)DOUBLE_EXPONENT_BIAS<<52)  // if x==1 return 0 as log(1)=0 . 1 has unbiased exp=0 and mantissa =0 [+hidden bit] 
  	return 0;  	
  int e = (v.u >> 52) - 0x3ff;
  if (e >= 0x400 || e == -0x3ff) /* x <= 0 or NaN/Inf or subnormal */  
  {   
    const b64u64_u minf = {.u = 0xfffull << 52};
    if (e == 0x400 || (e == 0xc00 && x != minf.f)) /* +Inf or NaN */
      return x + x;
    if (x <= 0.0) // x<=0
    {
      /* f(x<0) is NaN, f(+/-0) is -Inf and raises DivByZero */
      if (x<0) { //  x<0
#ifdef CORE_MATH_SUPPORT_ERRNO
        errno = EDOM;
#endif
        return NAN;
      }
      else {
#ifdef CORE_MATH_SUPPORT_ERRNO
        errno = ERANGE; // pole error
#endif
        return 1.0 / -0.0;// -INF
      }
    }
    if (e == -0x3ff) /* subnormal */
    {
      v.f *= 0x1p52;
      x_e -= 52;
      x_u = v.u;
    }
  }
  // end of special case handler
  /* now x > 0 */


  // log(x) = log(2^x_e * x_m)
  //        = x_e * log(2) + log(x_m)

  // Range reduction for log(x_m):
  // For each x_m, we would like to find r such that:
  //   -2^-8 <= r * x_m - 1 < 2^-7
  int shifted = x_u >> 45;
  int index = shifted & 0x7F;
  double r = RD[index];

  // Add unbiased exponent. Add an extra 1 if the 8 leading fractional bits are
  // all 1's.
  x_e += (int)((x_u + (1ULL << 45)) >> 52);
  double e_x = x_e;

  // hi is exact
  double hi = multiply_add(e_x, LOG_2_HI, LOG_R_DD[index].hi);
  // lo errors ~ e_x * LSB(LOG_2_LO) + LSB(LOG_R[index].lo) + rounding err
  //           <= 2 * (e_x * LSB(LOG_2_LO) + LSB(LOG_R[index].lo))
  double lo = multiply_add(e_x, LOG_2_LO, LOG_R_DD[index].lo);

  // Set m = 1.mantissa.
  uint64_t x_m = (x_u & 0x000FFFFFFFFFFFFFULL) | 0x3FF0000000000000ULL;
  v.u=x_m;
  double m = v.f;

  double u, u_sq;
  DoubleDouble r1;

  // Perform exact range reduction
  u = __builtin_fma(r, m, -1.0); // exact

  // Exact sum:
  //   r1.hi + r1.lo = e_x * log(2)_hi - log(r)_hi + u
  r1 = exact_add(hi, u);

  // Error of u_sq = ulp(u^2);
  u_sq = u * u;
  // Total error is bounded by ~ C * ulp(u^2).
  // Degree-7 minimax polynomial
  double p0 = multiply_add(u, LOG_COEFFS[1], LOG_COEFFS[0]);
  double p1 = multiply_add(u, LOG_COEFFS[3], LOG_COEFFS[2]);
  double p2 = multiply_add(u, LOG_COEFFS[5], LOG_COEFFS[4]);
  // Evaluate polynomial using Horner's Scheme:
  // With polyeval(x, a_0, a_1, ..., a_n) = a_n * x^n + ... + a_1 * x + a_0, we
  // evaluated it as:  a_0 + x * (a_1 + x * ( ... (a_(n-1) + x * a_n) ... ) ) ).
  // We will use FMA instructions if available.
  // Example: to evaluate x^3 + 2*x^2 + 3*x + 4, call
  //   polyeval( x, 4.0, 3.0, 2.0, 1.0 )
  double p=  multiply_add(u_sq,multiply_add(u_sq,multiply_add(u_sq,p2,p1),p0),lo+r1.lo);
  // double p = polyeval(u_sq, lo + r1.lo, p0, p1, p2);

  // Technically error of r1.lo is bounded by:
  //    hi*ulp(log(2)_lo) + C*ulp(u^2)
  // To simplify the error computation a bit, we replace |hi|*ulp(log(2)_lo)
  // with the upper bound: 2^11 * ulp(log(2)_lo) = 2^-85.
  // Total error is bounded by ~ C * ulp(u^2) + 2^-85.
  double err = multiply_add(u_sq, P_ERR, HI_ERR);

  // Lower bound from the result
  double left = r1.hi + (p - err);
  // Upper bound from the result
  double right = r1.hi + (p + err);

  // Ziv's test if fast pass is accurate enough.
  if (left == right)
    return left;
  double _cr_log_accurate(double x); // in cr_log.c
  return _cr_log_accurate(x);
}

