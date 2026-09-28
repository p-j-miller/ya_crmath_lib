/* ya_sin_cos.c
   ========

   Uses basically the same algorithms as the sin & cos functions in llvm-project-23.1.1 (\llvm-project-23.1.1\libc\src\__support\FPUtil & \llvm-project-23.1.1\libc\src\__support\math ).
   Note the llvm original is in C++, whereas this is in C.
   Where practical, variable names (and many comments) have been kept the same in this code as the llvm code, to make them easier to compare.
   It leverages the double-double maths library ya_double_double.
   
   This code assumes double FMA is available (either as a single instruction [ideally], or emulated in software) 
 
   The "accurate" evaluation is done using functions in cr_sin() and cr_cos() - in the test program these are only used 0.001% of the time so it was not worthwhile converting the accurate algorithms
   used in the llvm implementations(which might or might not be faster, but will make no difference to the average timing and negligible difference to the worse case).

   Checking for special cases (inf, NAN, etc) is done in the same was as cr_sin.c/cr_cos.c, as that solution supports CORE_MATH_SUPPORT_ERRNO.

 	According to  https://members.loria.fr/PZimmermann/papers/accuracy.pdf the llvm double sin/cos functions have a maximum error of 0.500ulp,
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
// #define INSTRUMENT /* if defined add print_ll_sin_counters() & print_ll_cos_counters()*/

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



#ifdef INSTRUMENT
#include <stdio.h>
 // for ll_sin(x) 
unsigned int scount_all=0,scount1=0,scount2=0,scount3=0,scount4=0,scount5=0,scount6=0;
void print_ll_sin_counters(void)
 {printf("ll_sin.c: count_all=%u count1=%u count2=%u count3=%u count4=%u count5=%u count6=%u\n",scount_all,scount1,scount2,scount3,scount4,scount5,scount6);
  printf("ll_sin.c: count_all=%.4f%% count1=%.4f%% count2=%.4f%% count3=%.4f%% count4=%.4f%% count5=%.4f%% count6=%.4f%%\n",100.0*scount_all/(double)scount_all,
  100.0*scount1/(double)scount_all,100.0*scount2/(double)scount_all,100.0*scount3/(double)scount_all,100.0*scount4/(double)scount_all,100.0*scount5/(double)scount_all,100.0*scount6/(double)scount_all);
 }
 // for ll_cos(x)
unsigned int ccount_all=0,ccount1=0,ccount2=0,ccount3=0,ccount4=0,ccount5=0,ccount6=0;
void print_ll_cos_counters(void)
 {printf("ll_cos.c: count_all=%u count1=%u count2=%u count3=%u count4=%u count5=%u count6=%u\n",ccount_all,ccount1,ccount2,ccount3,ccount4,ccount5,ccount6);
  printf("ll_cos.c: count_all=%.4f%% count1=%.4f%% count2=%.4f%% count3=%.4f%% count4=%.4f%% count5=%.4f%% count6=%.4f%%\n",100.0*ccount_all/(double)ccount_all,
  100.0*ccount1/(double)ccount_all,100.0*ccount2/(double)ccount_all,100.0*ccount3/(double)ccount_all,100.0*ccount4/(double)ccount_all,100.0*ccount5/(double)ccount_all,100.0*ccount6/(double)ccount_all);
 
 } 
#endif

#define multiply_add(a,b,c) ( ((a)*(b))+(c) ) /* this does not need to be a "real" FMA from an accuracy viewpoint , but compiler can use FMA instruction if its faster */

// from range_reduction_double_common.h

#define FAST_PASS_EXPONENT 16

// For 2^-7 < |x| < 2^16, return k and u such that:
//   k = round(x * 128/pi)
//   x mod pi/128 = x - k * pi/128 ~ u.hi + u.lo
// Error bound:
//   |(x - k * pi/128) - (u_hi + u_lo)| <= max(ulp(ulp(u_hi)), 2^-119)
//                                      <= 2^-111.
inline unsigned range_reduction_small(double x, DoubleDouble *u) 
{
  // Values of -pi/128 used for inputs with absolute value <= 2^16.
  // The first 3 parts are generated with (53 - 21 = 32)-bit precision, so that
  // the product k * MPI_OVER_128[i] is exact.
  // Generated by Sollya with:
  // > display = hexadecimal!;
  // > a = round(pi/128, 32, RN);
  // > b = round(pi/128 - a, 32, RN);
  // > c = round(pi/128 - a - b, D, RN);
  // > print(-a, ",", -b, ",", -c);
  const double MPI_OVER_128[3] = {-0x1.921fb544p-6, -0x1.0b4611a6p-40,
                                      -0x1.3198a2e037073p-75};
  const double ONE_TWENTY_EIGHT_OVER_PI_D = 0x1.45f306dc9c883p5;
  double prod_hi = x * ONE_TWENTY_EIGHT_OVER_PI_D;
  double kd = __builtin_rint(prod_hi);

  // Let y = x - k * (pi/128)
  // Then |y| < pi / 256
  // With extra rounding errors, we can bound |y| < 1.6 * 2^-7.
  double y_hi = multiply_add(kd, MPI_OVER_128[0], x); // Exact
  // |u.hi| < 1.6*2^-7
  u->hi = multiply_add(kd, MPI_OVER_128[1], y_hi);
  double u0 = y_hi - u->hi; // Exact
  // |u.lo| <= max(ulp(u.hi), |kd * MPI_OVER_128[2]|)
  double u1 = multiply_add(kd, MPI_OVER_128[1], u0); // Exact
  u->lo = multiply_add(kd, MPI_OVER_128[2], u1);
  // Error bound:
  // |x - k * pi/128| - (u.hi + u.lo) <= ulp(u.lo)
  //                                  <= ulp(max(ulp(u.hi), kd*MPI_OVER_128[2]))
  //                                  <= 2^(-7 - 104) = 2^-111.

  return kd;
}

// Digits of 2^(16*i) / pi, generated by Sollya with:
// > procedure ulp(x, n) { return 2^(floor(log2(abs(x))) - n); };
// > for i from 0 to 63 do {
//     if i < 3 then { pi_inv = 0.25 + 2^(16*(i - 3)) / pi; }
//     else { pi_inv = 2^(16*(i-3)) / pi; };
//     pn = nearestint(pi_inv);
//     pi_frac = pi_inv - pn;
//     a = round(pi_frac, 51, RN);
//     b = round(pi_frac - a, 51, RN);
//     c = round(pi_frac - a - b, 51, RN);
//     d = round(pi_frac - a - b - c, D, RN);
//     print("{", 2^7 * a, ",", 2^7 * b, ",", 2^7 * c, ",", 2^7 * d, "},");
// };
//
// Notice that for [0..2] the leading bit of 2^(16*(i - 3)) / pi is very small,
// so we add 0.25 so that the conditions for the algorithms are still satisfied,
// and one of those conditions guarantees that ulp(0.25 * x_reduced) >= 2, and
// will safely be discarded.

static const double ONE_TWENTY_EIGHT_OVER_PI[64][4] = {
    {0x1.0000000000014p5, 0x1.7cc1b727220a8p-49, 0x1.4fe13abe8fa9cp-101,
     -0x1.911f924eb5336p-153},
    {0x1.0000000145f3p5, 0x1.b727220a94fep-49, 0x1.3abe8fa9a6eep-101,
     0x1.b6c52b3278872p-155},
    {0x1.000145f306dc8p5, 0x1.c882a53f84ebp-47, -0x1.70565911f925p-101,
     0x1.4acc9e21c821p-153},
    {0x1.45f306dc9c884p5, -0x1.5ac07b1505c14p-47, -0x1.96447e493ad4cp-99,
     -0x1.b0ef1bef806bap-152},
    {-0x1.f246c6efab58p4, -0x1.ec5417056591p-49, -0x1.f924eb53361ep-101,
     0x1.c820ff28b1d5fp-153},
    {0x1.391054a7f09d4p4, 0x1.f47d4d377036cp-48, 0x1.8a5664f10e41p-100,
     0x1.fe5163abdebbcp-154},
    {0x1.529fc2757d1f4p2, 0x1.34ddc0db62958p-50, 0x1.93c439041fe5p-102,
     0x1.63abdebbc561bp-154},
    {-0x1.ec5417056591p-1, -0x1.f924eb53361ep-53, 0x1.c820ff28b1d6p-105,
     -0x1.0a21d4f246dc9p-157},
    {-0x1.505c1596447e4p5, -0x1.275a99b0ef1cp-48, 0x1.07f9458eaf7bp-100,
     -0x1.0ea79236e4717p-152},
    {-0x1.596447e493ad4p1, -0x1.9b0ef1bef806cp-52, 0x1.63abdebbc561cp-106,
     -0x1.1b7238b7b645ap-159},
    {0x1.bb81b6c52b328p5, -0x1.de37df00d74e4p-49, 0x1.5ef5de2b0db94p-101,
     -0x1.c8e2ded9169p-153},
    {0x1.b6c52b3278874p5, -0x1.f7c035d38a844p-47, 0x1.778ac36e48dc8p-99,
     -0x1.6f6c8b47fe6dbp-152},
    {0x1.2b3278872084p5, -0x1.ae9c5421443a8p-50, -0x1.e48db91c5bdb4p-102,
     0x1.d2e006492eea1p-154},
    {-0x1.8778df7c035d4p5, 0x1.d5ef5de2b0db8p-49, 0x1.2371d2126e97p-101,
     0x1.924bba8274648p-160},
    {-0x1.bef806ba71508p4, -0x1.443a9e48db91cp-50, -0x1.6f6c8b47fe6dcp-104,
     0x1.77504e8c90e7fp-157},
    {-0x1.ae9c5421443a8p-2, -0x1.e48db91c5bdb4p-54, 0x1.d2e006492eeap-106,
     0x1.3a32439fc3bd6p-159},
    {-0x1.38a84288753c8p5, -0x1.1b7238b7b645cp-47, 0x1.c00c925dd413cp-99,
     -0x1.cdbc603c429c7p-151},
    {-0x1.0a21d4f246dc8p3, -0x1.c5bdb22d1ff9cp-50, 0x1.25dd413a32438p-103,
     0x1.fc3bd63962535p-155},
    {-0x1.d4f246dc8e2ep3, 0x1.26e9700324978p-49, -0x1.5f62e6de301e4p-102,
     0x1.eb1cb129a73efp-154},
    {-0x1.236e4716f6c8cp4, 0x1.700324977505p-49, -0x1.736f180f10a7p-101,
     -0x1.a76b2c608bbeep-153},
    {0x1.b8e909374b8p4, 0x1.924bba8274648p-48, 0x1.cfe1deb1cb128p-102,
     0x1.a73ee88235f53p-154},
    {0x1.09374b801924cp4, -0x1.15f62e6de302p-50, 0x1.deb1cb129a74p-102,
     -0x1.177dca0ad144cp-154},
    {-0x1.68ffcdb688afcp3, 0x1.d1921cfe1debp-50, 0x1.cb129a73ee884p-102,
     -0x1.ca0ad144bb7b1p-154},
    {0x1.924bba8274648p0, 0x1.cfe1deb1cb128p-54, 0x1.a73ee88235f54p-106,
     -0x1.144bb7b16639p-158},
    {-0x1.a22bec5cdbc6p5, -0x1.e214e34ed658cp-50, -0x1.177dca0ad144cp-106,
     0x1.213a671c09ad1p-160},
    {0x1.3a32439fc3bd8p1, -0x1.c69dacb1822fp-51, 0x1.1afa975da2428p-105,
     -0x1.6638fd94ba082p-158},
    {-0x1.b78c0788538d4p4, 0x1.29a73ee88236p-50, -0x1.5a28976f62cc8p-103,
     0x1.c09ad17df904ep-156},
    {0x1.fc3bd63962534p5, 0x1.cfba208d7d4bcp-48, -0x1.12edec598e3f8p-100,
     0x1.ad17df904e647p-152},
    {-0x1.4e34ed658c118p2, 0x1.046bea5d7689p-51, 0x1.3a671c09ad17cp-104,
     0x1.f904e64758e61p-156},
    {0x1.62534e7dd1048p5, -0x1.415a28976f62cp-47, -0x1.8e3f652e8207p-100,
     0x1.3991d63983534p-154},
    {-0x1.63045df7282b4p4, -0x1.44bb7b16638fcp-50, -0x1.94ba081bec67p-102,
     0x1.d639835339f4ap-154},
    {0x1.d1046bea5d768p5, 0x1.213a671c09adp-48, 0x1.7df904e64759p-100,
     -0x1.9f2b3182d8defp-152},
    {0x1.afa975da24274p3, 0x1.9c7026b45f7e4p-50, 0x1.3991d63983534p-106,
     -0x1.82d8dee81d108p-160},
    {-0x1.a28976f62cc7p5, -0x1.fb29741037d8cp-47, -0x1.b8a719f2b3184p-100,
     0x1.272117e2ef7e5p-152},
    {-0x1.76f62cc71fb28p5, -0x1.741037d8cdc54p-47, 0x1.cc1a99cfa4e44p-101,
     -0x1.d03a21036be27p-153},
    {0x1.d338e04d68bfp5, -0x1.bec66e29c67ccp-50, 0x1.339f49c845f8cp-102,
     -0x1.081b5f13801dap-156},
    {0x1.c09ad17df905p4, -0x1.9b8a719f2b318p-48, -0x1.6c6f740e8840cp-103,
     -0x1.af89c00ed0004p-155},
    {0x1.68befc827323cp5, -0x1.38cf9598c16c8p-47, 0x1.08bf177bf2508p-99,
     -0x1.3801da00087eap-152},
    {-0x1.037d8cdc538dp5, 0x1.a99cfa4e422fcp-49, 0x1.77bf250763ffp-103,
     0x1.2fffbc0b301fep-155},
    {-0x1.8cdc538cf9598p5, -0x1.82d8dee81d108p-48, -0x1.b5f13801dap-104,
     -0x1.0fd33f8086877p-157},
    {-0x1.4e33e566305bp3, -0x1.bdd03a21036cp-49, 0x1.d8ffc4bffef04p-101,
     -0x1.33f80868773a5p-153},
    {-0x1.f2b3182d8dee8p4, -0x1.d1081b5f138p-52, -0x1.da00087e99fcp-104,
     -0x1.0d0ee74a5f593p-158},
    {-0x1.8c16c6f740e88p5, -0x1.036be27003b4p-49, -0x1.0fd33f8086878p-109,
     0x1.8b5a0a6d1f6d3p-162},
    {0x1.3908bf177bf24p5, 0x1.0763ff12fffbcp-47, 0x1.6603fbcbc462cp-104,
     0x1.6829b47db4dap-156},
    {0x1.7e2ef7e4a0ec8p4, -0x1.da00087e99fcp-56, -0x1.0d0ee74a5f594p-110,
     0x1.1f6d367ecf27dp-162},
    {-0x1.081b5f13801dcp4, 0x1.fff7816603fbcp-48, 0x1.788c5ad05369p-101,
     -0x1.25930261b069fp-155},
    {-0x1.af89c00ed0004p5, -0x1.fa67f010d0ee8p-50, 0x1.6b414da3eda6cp-103,
     0x1.fb3c9f2c26dd4p-156},
    {-0x1.c00ed00043f4cp5, -0x1.fc04343b9d298p-48, 0x1.4da3eda6cfdap-103,
     -0x1.b069ec9161738p-155},
    {0x1.2fffbc0b301fcp5, 0x1.e5e2316b414dcp-47, -0x1.c125930261b08p-99,
     0x1.6136e9e8c7ecdp-151},
    {-0x1.0fd33f8086878p3, 0x1.8b5a0a6d1f6d4p-50, -0x1.30261b069ec9p-103,
     -0x1.61738132c3403p-155},
    {-0x1.9fc04343b9d28p4, -0x1.7d64b824b2604p-48, -0x1.86c1a7b24585cp-101,
     -0x1.c09961a015d29p-154},
    {-0x1.0d0ee74a5f594p2, 0x1.1f6d367ecf27cp-50, 0x1.6136e9e8c7eccp-103,
     0x1.3cbfd45aea4f7p-155},
    {-0x1.dce94beb25c14p5, 0x1.a6cfd9e4f9614p-47, -0x1.22c2e70265868p-100,
     -0x1.5d28ad8453814p-158},
    {-0x1.4beb25c12593p5, -0x1.30d834f648b0cp-50, 0x1.8fd9a797fa8b4p-104,
     0x1.d49eeb1faf97cp-156},
    {0x1.b47db4d9fb3c8p4, 0x1.f2c26dd3d18fcp-48, 0x1.9a797fa8b5d48p-100,
     0x1.eeb1faf97c5edp-152},
    {-0x1.25930261b06ap5, 0x1.36e9e8c7ecd3cp-47, 0x1.7fa8b5d49eebp-100,
     0x1.faf97c5ecf41dp-152},
    {0x1.fb3c9f2c26dd4p4, -0x1.738132c3402bcp-51, 0x1.aea4f758fd7ccp-103,
     -0x1.d0985f18c10ebp-159},
    {-0x1.b069ec9161738p5, -0x1.32c3402ba515cp-51, 0x1.eeb1faf97c5ecp-104,
     0x1.e839cfbc52949p-157},
    {-0x1.ec9161738132cp5, -0x1.a015d28ad8454p-50, 0x1.faf97c5ecf41cp-104,
     0x1.cfbc529497536p-157},
    {-0x1.61738132c3404p5, 0x1.45aea4f758fd8p-47, -0x1.a0e84c2f8c608p-102,
     -0x1.d6b5b45650128p-156},
    {0x1.fb34f2ff516bcp3, -0x1.6c229c0a0d074p-49, -0x1.30be31821d6b4p-104,
     -0x1.b4565012813b8p-156},
    {0x1.3cbfd45aea4f8p5, -0x1.4e050683a130cp-48, 0x1.ce7de294a4ba8p-104,
     0x1.afed7ec47e357p-156},
    {-0x1.5d28ad8453814p2, -0x1.a0e84c2f8c608p-54, -0x1.d6b5b45650128p-108,
     -0x1.3b81ca8bdea7fp-164},
    {-0x1.15b08a702834p5, -0x1.d0985f18c10ecp-47, 0x1.4a4ba9afed7ecp-100,
     0x1.1f8d5d0856033p-154},
};

// New PMi definitions
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


// For large range |x| >= 2^16, we perform the range reduction computations as:
//   u = x - k * pi/128 = (pi/128) * (x * (128/pi) - k).
// We use the exponent of x to find 4 double-chunks of 128/pi:
// c_hi, c_mid, c_lo, c_lo_2 such that:
//   1) ulp(round(x * c_hi, D, RN)) >= 2^8 = 256,
//   2) If x * c_hi = ph_hi + ph_lo and x * c_mid = pm_hi + pm_lo, then
//        min(ulp(ph_lo), ulp(pm_hi)) >= 2^-53.
// This will allow us to drop the high part ph_hi and the addition:
//   (ph_lo + pm_hi) mod 1
// can be exactly representable in a double precision.
// This will allow us to do split the computations as:
//   (x * 256/pi) ~ x * (c_hi + c_mid + c_lo + c_lo_2)    (mod 256)
//                ~ (ph_lo + pm_hi) + (pm_lo + x * c_lo) + x * c_lo_2.
// Then,
//   round(x * 128/pi) = round(ph_lo + pm_hi)    (mod 256)
// And the high part of fractional part of (x * 128/pi) can simply be:
//   {x * 128/pi}_hi = {ph_lo + pm_hi}.
// To prevent overflow when x is very large, we simply scale up
// (c_hi, c_mid, c_lo, c_lo_2) by a fixed power of 2 (based on the index) and
// scale down x by the same amount.

// from range_reduction_double_fma.h -> here we must use FMA
static inline unsigned LargeRangeReduction_fast(double x, DoubleDouble *u) 
{
  //int x_e_m62 = xbits.get_biased_exponent() - (FPBits::EXP_BIAS + 62);
  b64u64_u t = {.f = x};
  int e = (t.u >> DOUBLE_MANTISSA_BITS) & 0x7ff; /* get biased exponent */
  int x_e_m62=e-(DOUBLE_EXPONENT_BIAS+62);
  unsigned idx=((x_e_m62 >> 4) + 3);
  // Scale x down by 2^(-(16 * (idx - 3))
  // xbits.set_biased_exponent((x_e_m62 & 15) + FPBits::EXP_BIAS + 62);
  t.u &= ~(((uint64_t)0x7ff) <<DOUBLE_MANTISSA_BITS); // blank out existing exponent
  t.u |= (((uint64_t)((x_e_m62 & 15) + DOUBLE_EXPONENT_BIAS + 62))&(uint64_t)0x7ff)<<DOUBLE_MANTISSA_BITS;// set exponent
  // 2^62 <= |x_reduced| < 2^(62 + 16) = 2^78
  //double x_reduced = xbits.get_val();
  double x_reduced=t.f;
  // x * c_hi = ph.hi + ph.lo exactly.
  DoubleDouble ph = exact_mult(x_reduced, ONE_TWENTY_EIGHT_OVER_PI[idx][0]);
  // x * c_mid = pm.hi + pm.lo exactly.
  DoubleDouble pm = exact_mult(x_reduced, ONE_TWENTY_EIGHT_OVER_PI[idx][1]);
  // x * c_lo = pl.hi + pl.lo exactly.
  DoubleDouble pl = exact_mult(x_reduced, ONE_TWENTY_EIGHT_OVER_PI[idx][2]);
  // Extract integral parts and fractional parts of (ph.lo + pm.hi).
  double sum_hi = ph.lo + pm.hi;
  //double kd = fputil::nearest_integer(sum_hi);
  double kd =  __builtin_rint(sum_hi);// rint() - but should convert to a single instruction with sse4.2 or above
  // x * 128/pi mod 1 ~ y_hi + y_mid + y_lo
  double y_hi = (ph.lo - kd) + pm.hi; // Exact
  DoubleDouble y_mid = exact_add(pm.lo, pl.hi);
  double y_lo = pl.lo;

  // y_l = x * c_lo_2 + pl.lo
  double y_l = __builtin_fma(x_reduced, ONE_TWENTY_EIGHT_OVER_PI[idx][3], y_lo);
  DoubleDouble y = exact_add(y_hi, y_mid.hi);
  y.lo += (y_mid.lo + y_l);

  // Digits of pi/128, generated by Sollya with:
  // > a = round(pi/128, D, RN);
  // > b = round(pi/128 - a, D, RN);
  const DoubleDouble PI_OVER_128_DD = {0x1.1a62633145c07p-60,
                                           0x1.921fb54442d18p-6};

  // Error bound: with {a} denote the fractional part of a, i.e.:
  //   {a} = a - round(a)
  // Then,
  //   | {x * 128/pi} - (y_hi + y_lo) | <=  ulp(ulp(y_hi)) <= 2^-105
  //   | {x mod pi/128} - (u.hi + u.lo) | < 2 * 2^-6 * 2^-105 = 2^-110
  DoubleDouble r_u = mult_dd_dd(y, PI_OVER_128_DD);
  *u=r_u;
  return kd; // note kd is a double and we return an unsigned
}

// Lookup table for sin(k * pi / 128) with k = 0, ..., 255.
// Table is generated with Sollya as follow:
// > display = hexadecimal;
// > for k from 0 to 255 do {
//     a = D(sin(k * pi/128)); };
//     b = D(sin(k * pi/128) - a);
//     print("{", b, ",", a, "},");
//   };
static const DoubleDouble SIN_K_PI_OVER_128[] = {
    {0, 0},
    {-0x1.b1d63091a013p-64, 0x1.92155f7a3667ep-6},
    {-0x1.912bd0d569a9p-61, 0x1.91f65f10dd814p-5},
    {-0x1.9a088a8bf6b2cp-59, 0x1.2d52092ce19f6p-4},
    {-0x1.e2718d26ed688p-60, 0x1.917a6bc29b42cp-4},
    {0x1.a2704729ae56dp-59, 0x1.f564e56a9730ep-4},
    {0x1.13000a89a11ep-58, 0x1.2c8106e8e613ap-3},
    {0x1.531ff779ddac6p-57, 0x1.5e214448b3fc6p-3},
    {-0x1.26d19b9ff8d82p-57, 0x1.8f8b83c69a60bp-3},
    {-0x1.af1439e521935p-62, 0x1.c0b826a7e4f63p-3},
    {-0x1.42deef11da2c4p-57, 0x1.f19f97b215f1bp-3},
    {0x1.824c20ab7aa9ap-56, 0x1.111d262b1f677p-2},
    {-0x1.5d28da2c4612dp-56, 0x1.294062ed59f06p-2},
    {0x1.0c97c4afa2518p-56, 0x1.4135c94176601p-2},
    {-0x1.efdc0d58cf62p-62, 0x1.58f9a75ab1fddp-2},
    {-0x1.44b19e0864c5dp-56, 0x1.7088530fa459fp-2},
    {-0x1.72cedd3d5a61p-57, 0x1.87de2a6aea963p-2},
    {0x1.6da81290bdbabp-57, 0x1.9ef7943a8ed8ap-2},
    {0x1.5b362cb974183p-57, 0x1.b5d1009e15ccp-2},
    {0x1.6850e59c37f8fp-58, 0x1.cc66e9931c45ep-2},
    {0x1.e0d891d3c6841p-58, 0x1.e2b5d3806f63bp-2},
    {-0x1.2ec1fc1b776b8p-60, 0x1.f8ba4dbf89abap-2},
    {-0x1.a5a014347406cp-55, 0x1.073879922ffeep-1},
    {-0x1.ef23b69abe4f1p-55, 0x1.11eb3541b4b23p-1},
    {0x1.b25dd267f66p-55, 0x1.1c73b39ae68c8p-1},
    {-0x1.5da743ef3770cp-55, 0x1.26d054cdd12dfp-1},
    {-0x1.efcc626f74a6fp-57, 0x1.30ff7fce17035p-1},
    {0x1.e3e25e3954964p-56, 0x1.3affa292050b9p-1},
    {0x1.8076a2cfdc6b3p-57, 0x1.44cf325091dd6p-1},
    {0x1.3c293edceb327p-57, 0x1.4e6cabbe3e5e9p-1},
    {-0x1.75720992bfbb2p-55, 0x1.57d69348cecap-1},
    {-0x1.251b352ff2a37p-56, 0x1.610b7551d2cdfp-1},
    {-0x1.bdd3413b26456p-55, 0x1.6a09e667f3bcdp-1},
    {0x1.0d4ef0f1d915cp-55, 0x1.72d0837efff96p-1},
    {-0x1.0f537acdf0ad7p-56, 0x1.7b5df226aafafp-1},
    {-0x1.6f420f8ea3475p-56, 0x1.83b0e0bff976ep-1},
    {-0x1.2c5e12ed1336dp-55, 0x1.8bc806b151741p-1},
    {0x1.3d419a920df0bp-55, 0x1.93a22499263fbp-1},
    {-0x1.30ee286712474p-55, 0x1.9b3e047f38741p-1},
    {-0x1.128bb015df175p-56, 0x1.a29a7a0462782p-1},
    {0x1.9f630e8b6dac8p-60, 0x1.a9b66290ea1a3p-1},
    {-0x1.926da300ffccep-55, 0x1.b090a581502p-1},
    {-0x1.bc69f324e6d61p-55, 0x1.b728345196e3ep-1},
    {-0x1.825a732ac700ap-55, 0x1.bd7c0ac6f952ap-1},
    {-0x1.6e0b1757c8d07p-56, 0x1.c38b2f180bdb1p-1},
    {-0x1.2fb761e946603p-58, 0x1.c954b213411f5p-1},
    {-0x1.e7b6bb5ab58aep-58, 0x1.ced7af43cc773p-1},
    {-0x1.4ef5295d25af2p-55, 0x1.d4134d14dc93ap-1},
    {0x1.457e610231ac2p-56, 0x1.d906bcf328d46p-1},
    {0x1.83c37c6107db3p-55, 0x1.ddb13b6ccc23cp-1},
    {-0x1.014c76c126527p-55, 0x1.e212104f686e5p-1},
    {-0x1.16b56f2847754p-57, 0x1.e6288ec48e112p-1},
    {0x1.760b1e2e3f81ep-55, 0x1.e9f4156c62ddap-1},
    {0x1.e82c791f59cc2p-56, 0x1.ed740e7684963p-1},
    {0x1.52c7adc6b4989p-56, 0x1.f0a7efb9230d7p-1},
    {-0x1.d7bafb51f72e6p-56, 0x1.f38f3ac64e589p-1},
    {0x1.562172a361fd3p-56, 0x1.f6297cff75cbp-1},
    {0x1.ab256778ffcb6p-56, 0x1.f8764fa714ba9p-1},
    {-0x1.7a0a8ca13571fp-55, 0x1.fa7557f08a517p-1},
    {0x1.1ec8668ecaceep-55, 0x1.fc26470e19fd3p-1},
    {-0x1.87df6378811c7p-55, 0x1.fd88da3d12526p-1},
    {0x1.521ecd0c67e35p-57, 0x1.fe9cdad01883ap-1},
    {-0x1.c57bc2e24aa15p-57, 0x1.ff621e3796d7ep-1},
    {-0x1.1354d4556e4cbp-55, 0x1.ffd886084cd0dp-1},
    {0, 1},
    {-0x1.1354d4556e4cbp-55, 0x1.ffd886084cd0dp-1},
    {-0x1.c57bc2e24aa15p-57, 0x1.ff621e3796d7ep-1},
    {0x1.521ecd0c67e35p-57, 0x1.fe9cdad01883ap-1},
    {-0x1.87df6378811c7p-55, 0x1.fd88da3d12526p-1},
    {0x1.1ec8668ecaceep-55, 0x1.fc26470e19fd3p-1},
    {-0x1.7a0a8ca13571fp-55, 0x1.fa7557f08a517p-1},
    {0x1.ab256778ffcb6p-56, 0x1.f8764fa714ba9p-1},
    {0x1.562172a361fd3p-56, 0x1.f6297cff75cbp-1},
    {-0x1.d7bafb51f72e6p-56, 0x1.f38f3ac64e589p-1},
    {0x1.52c7adc6b4989p-56, 0x1.f0a7efb9230d7p-1},
    {0x1.e82c791f59cc2p-56, 0x1.ed740e7684963p-1},
    {0x1.760b1e2e3f81ep-55, 0x1.e9f4156c62ddap-1},
    {-0x1.16b56f2847754p-57, 0x1.e6288ec48e112p-1},
    {-0x1.014c76c126527p-55, 0x1.e212104f686e5p-1},
    {0x1.83c37c6107db3p-55, 0x1.ddb13b6ccc23cp-1},
    {0x1.457e610231ac2p-56, 0x1.d906bcf328d46p-1},
    {-0x1.4ef5295d25af2p-55, 0x1.d4134d14dc93ap-1},
    {-0x1.e7b6bb5ab58aep-58, 0x1.ced7af43cc773p-1},
    {-0x1.2fb761e946603p-58, 0x1.c954b213411f5p-1},
    {-0x1.6e0b1757c8d07p-56, 0x1.c38b2f180bdb1p-1},
    {-0x1.825a732ac700ap-55, 0x1.bd7c0ac6f952ap-1},
    {-0x1.bc69f324e6d61p-55, 0x1.b728345196e3ep-1},
    {-0x1.926da300ffccep-55, 0x1.b090a581502p-1},
    {0x1.9f630e8b6dac8p-60, 0x1.a9b66290ea1a3p-1},
    {-0x1.128bb015df175p-56, 0x1.a29a7a0462782p-1},
    {-0x1.30ee286712474p-55, 0x1.9b3e047f38741p-1},
    {0x1.3d419a920df0bp-55, 0x1.93a22499263fbp-1},
    {-0x1.2c5e12ed1336dp-55, 0x1.8bc806b151741p-1},
    {-0x1.6f420f8ea3475p-56, 0x1.83b0e0bff976ep-1},
    {-0x1.0f537acdf0ad7p-56, 0x1.7b5df226aafafp-1},
    {0x1.0d4ef0f1d915cp-55, 0x1.72d0837efff96p-1},
    {-0x1.bdd3413b26456p-55, 0x1.6a09e667f3bcdp-1},
    {-0x1.251b352ff2a37p-56, 0x1.610b7551d2cdfp-1},
    {-0x1.75720992bfbb2p-55, 0x1.57d69348cecap-1},
    {0x1.3c293edceb327p-57, 0x1.4e6cabbe3e5e9p-1},
    {0x1.8076a2cfdc6b3p-57, 0x1.44cf325091dd6p-1},
    {0x1.e3e25e3954964p-56, 0x1.3affa292050b9p-1},
    {-0x1.efcc626f74a6fp-57, 0x1.30ff7fce17035p-1},
    {-0x1.5da743ef3770cp-55, 0x1.26d054cdd12dfp-1},
    {0x1.b25dd267f66p-55, 0x1.1c73b39ae68c8p-1},
    {-0x1.ef23b69abe4f1p-55, 0x1.11eb3541b4b23p-1},
    {-0x1.a5a014347406cp-55, 0x1.073879922ffeep-1},
    {-0x1.2ec1fc1b776b8p-60, 0x1.f8ba4dbf89abap-2},
    {0x1.e0d891d3c6841p-58, 0x1.e2b5d3806f63bp-2},
    {0x1.6850e59c37f8fp-58, 0x1.cc66e9931c45ep-2},
    {0x1.5b362cb974183p-57, 0x1.b5d1009e15ccp-2},
    {0x1.6da81290bdbabp-57, 0x1.9ef7943a8ed8ap-2},
    {-0x1.72cedd3d5a61p-57, 0x1.87de2a6aea963p-2},
    {-0x1.44b19e0864c5dp-56, 0x1.7088530fa459fp-2},
    {-0x1.efdc0d58cf62p-62, 0x1.58f9a75ab1fddp-2},
    {0x1.0c97c4afa2518p-56, 0x1.4135c94176601p-2},
    {-0x1.5d28da2c4612dp-56, 0x1.294062ed59f06p-2},
    {0x1.824c20ab7aa9ap-56, 0x1.111d262b1f677p-2},
    {-0x1.42deef11da2c4p-57, 0x1.f19f97b215f1bp-3},
    {-0x1.af1439e521935p-62, 0x1.c0b826a7e4f63p-3},
    {-0x1.26d19b9ff8d82p-57, 0x1.8f8b83c69a60bp-3},
    {0x1.531ff779ddac6p-57, 0x1.5e214448b3fc6p-3},
    {0x1.13000a89a11ep-58, 0x1.2c8106e8e613ap-3},
    {0x1.a2704729ae56dp-59, 0x1.f564e56a9730ep-4},
    {-0x1.e2718d26ed688p-60, 0x1.917a6bc29b42cp-4},
    {-0x1.9a088a8bf6b2cp-59, 0x1.2d52092ce19f6p-4},
    {-0x1.912bd0d569a9p-61, 0x1.91f65f10dd814p-5},
    {-0x1.b1d63091a013p-64, 0x1.92155f7a3667ep-6},
    {0, 0},
    {0x1.b1d63091a013p-64, -0x1.92155f7a3667ep-6},
    {0x1.912bd0d569a9p-61, -0x1.91f65f10dd814p-5},
    {0x1.9a088a8bf6b2cp-59, -0x1.2d52092ce19f6p-4},
    {0x1.e2718d26ed688p-60, -0x1.917a6bc29b42cp-4},
    {-0x1.a2704729ae56dp-59, -0x1.f564e56a9730ep-4},
    {-0x1.13000a89a11ep-58, -0x1.2c8106e8e613ap-3},
    {-0x1.531ff779ddac6p-57, -0x1.5e214448b3fc6p-3},
    {0x1.26d19b9ff8d82p-57, -0x1.8f8b83c69a60bp-3},
    {0x1.af1439e521935p-62, -0x1.c0b826a7e4f63p-3},
    {0x1.42deef11da2c4p-57, -0x1.f19f97b215f1bp-3},
    {-0x1.824c20ab7aa9ap-56, -0x1.111d262b1f677p-2},
    {0x1.5d28da2c4612dp-56, -0x1.294062ed59f06p-2},
    {-0x1.0c97c4afa2518p-56, -0x1.4135c94176601p-2},
    {0x1.efdc0d58cf62p-62, -0x1.58f9a75ab1fddp-2},
    {0x1.44b19e0864c5dp-56, -0x1.7088530fa459fp-2},
    {0x1.72cedd3d5a61p-57, -0x1.87de2a6aea963p-2},
    {-0x1.6da81290bdbabp-57, -0x1.9ef7943a8ed8ap-2},
    {-0x1.5b362cb974183p-57, -0x1.b5d1009e15ccp-2},
    {-0x1.6850e59c37f8fp-58, -0x1.cc66e9931c45ep-2},
    {-0x1.e0d891d3c6841p-58, -0x1.e2b5d3806f63bp-2},
    {0x1.2ec1fc1b776b8p-60, -0x1.f8ba4dbf89abap-2},
    {0x1.a5a014347406cp-55, -0x1.073879922ffeep-1},
    {0x1.ef23b69abe4f1p-55, -0x1.11eb3541b4b23p-1},
    {-0x1.b25dd267f66p-55, -0x1.1c73b39ae68c8p-1},
    {0x1.5da743ef3770cp-55, -0x1.26d054cdd12dfp-1},
    {0x1.efcc626f74a6fp-57, -0x1.30ff7fce17035p-1},
    {-0x1.e3e25e3954964p-56, -0x1.3affa292050b9p-1},
    {-0x1.8076a2cfdc6b3p-57, -0x1.44cf325091dd6p-1},
    {-0x1.3c293edceb327p-57, -0x1.4e6cabbe3e5e9p-1},
    {0x1.75720992bfbb2p-55, -0x1.57d69348cecap-1},
    {0x1.251b352ff2a37p-56, -0x1.610b7551d2cdfp-1},
    {0x1.bdd3413b26456p-55, -0x1.6a09e667f3bcdp-1},
    {-0x1.0d4ef0f1d915cp-55, -0x1.72d0837efff96p-1},
    {0x1.0f537acdf0ad7p-56, -0x1.7b5df226aafafp-1},
    {0x1.6f420f8ea3475p-56, -0x1.83b0e0bff976ep-1},
    {0x1.2c5e12ed1336dp-55, -0x1.8bc806b151741p-1},
    {-0x1.3d419a920df0bp-55, -0x1.93a22499263fbp-1},
    {0x1.30ee286712474p-55, -0x1.9b3e047f38741p-1},
    {0x1.128bb015df175p-56, -0x1.a29a7a0462782p-1},
    {-0x1.9f630e8b6dac8p-60, -0x1.a9b66290ea1a3p-1},
    {0x1.926da300ffccep-55, -0x1.b090a581502p-1},
    {0x1.bc69f324e6d61p-55, -0x1.b728345196e3ep-1},
    {0x1.825a732ac700ap-55, -0x1.bd7c0ac6f952ap-1},
    {0x1.6e0b1757c8d07p-56, -0x1.c38b2f180bdb1p-1},
    {0x1.2fb761e946603p-58, -0x1.c954b213411f5p-1},
    {0x1.e7b6bb5ab58aep-58, -0x1.ced7af43cc773p-1},
    {0x1.4ef5295d25af2p-55, -0x1.d4134d14dc93ap-1},
    {-0x1.457e610231ac2p-56, -0x1.d906bcf328d46p-1},
    {-0x1.83c37c6107db3p-55, -0x1.ddb13b6ccc23cp-1},
    {0x1.014c76c126527p-55, -0x1.e212104f686e5p-1},
    {0x1.16b56f2847754p-57, -0x1.e6288ec48e112p-1},
    {-0x1.760b1e2e3f81ep-55, -0x1.e9f4156c62ddap-1},
    {-0x1.e82c791f59cc2p-56, -0x1.ed740e7684963p-1},
    {-0x1.52c7adc6b4989p-56, -0x1.f0a7efb9230d7p-1},
    {0x1.d7bafb51f72e6p-56, -0x1.f38f3ac64e589p-1},
    {-0x1.562172a361fd3p-56, -0x1.f6297cff75cbp-1},
    {-0x1.ab256778ffcb6p-56, -0x1.f8764fa714ba9p-1},
    {0x1.7a0a8ca13571fp-55, -0x1.fa7557f08a517p-1},
    {-0x1.1ec8668ecaceep-55, -0x1.fc26470e19fd3p-1},
    {0x1.87df6378811c7p-55, -0x1.fd88da3d12526p-1},
    {-0x1.521ecd0c67e35p-57, -0x1.fe9cdad01883ap-1},
    {0x1.c57bc2e24aa15p-57, -0x1.ff621e3796d7ep-1},
    {0x1.1354d4556e4cbp-55, -0x1.ffd886084cd0dp-1},
    {0, -1},
    {0x1.1354d4556e4cbp-55, -0x1.ffd886084cd0dp-1},
    {0x1.c57bc2e24aa15p-57, -0x1.ff621e3796d7ep-1},
    {-0x1.521ecd0c67e35p-57, -0x1.fe9cdad01883ap-1},
    {0x1.87df6378811c7p-55, -0x1.fd88da3d12526p-1},
    {-0x1.1ec8668ecaceep-55, -0x1.fc26470e19fd3p-1},
    {0x1.7a0a8ca13571fp-55, -0x1.fa7557f08a517p-1},
    {-0x1.ab256778ffcb6p-56, -0x1.f8764fa714ba9p-1},
    {-0x1.562172a361fd3p-56, -0x1.f6297cff75cbp-1},
    {0x1.d7bafb51f72e6p-56, -0x1.f38f3ac64e589p-1},
    {-0x1.52c7adc6b4989p-56, -0x1.f0a7efb9230d7p-1},
    {-0x1.e82c791f59cc2p-56, -0x1.ed740e7684963p-1},
    {-0x1.760b1e2e3f81ep-55, -0x1.e9f4156c62ddap-1},
    {0x1.16b56f2847754p-57, -0x1.e6288ec48e112p-1},
    {0x1.014c76c126527p-55, -0x1.e212104f686e5p-1},
    {-0x1.83c37c6107db3p-55, -0x1.ddb13b6ccc23cp-1},
    {-0x1.457e610231ac2p-56, -0x1.d906bcf328d46p-1},
    {0x1.4ef5295d25af2p-55, -0x1.d4134d14dc93ap-1},
    {0x1.e7b6bb5ab58aep-58, -0x1.ced7af43cc773p-1},
    {0x1.2fb761e946603p-58, -0x1.c954b213411f5p-1},
    {0x1.6e0b1757c8d07p-56, -0x1.c38b2f180bdb1p-1},
    {0x1.825a732ac700ap-55, -0x1.bd7c0ac6f952ap-1},
    {0x1.bc69f324e6d61p-55, -0x1.b728345196e3ep-1},
    {0x1.926da300ffccep-55, -0x1.b090a581502p-1},
    {-0x1.9f630e8b6dac8p-60, -0x1.a9b66290ea1a3p-1},
    {0x1.128bb015df175p-56, -0x1.a29a7a0462782p-1},
    {0x1.30ee286712474p-55, -0x1.9b3e047f38741p-1},
    {-0x1.3d419a920df0bp-55, -0x1.93a22499263fbp-1},
    {0x1.2c5e12ed1336dp-55, -0x1.8bc806b151741p-1},
    {0x1.6f420f8ea3475p-56, -0x1.83b0e0bff976ep-1},
    {0x1.0f537acdf0ad7p-56, -0x1.7b5df226aafafp-1},
    {-0x1.0d4ef0f1d915cp-55, -0x1.72d0837efff96p-1},
    {0x1.bdd3413b26456p-55, -0x1.6a09e667f3bcdp-1},
    {0x1.251b352ff2a37p-56, -0x1.610b7551d2cdfp-1},
    {0x1.75720992bfbb2p-55, -0x1.57d69348cecap-1},
    {-0x1.3c293edceb327p-57, -0x1.4e6cabbe3e5e9p-1},
    {-0x1.8076a2cfdc6b3p-57, -0x1.44cf325091dd6p-1},
    {-0x1.e3e25e3954964p-56, -0x1.3affa292050b9p-1},
    {0x1.efcc626f74a6fp-57, -0x1.30ff7fce17035p-1},
    {0x1.5da743ef3770cp-55, -0x1.26d054cdd12dfp-1},
    {-0x1.b25dd267f66p-55, -0x1.1c73b39ae68c8p-1},
    {0x1.ef23b69abe4f1p-55, -0x1.11eb3541b4b23p-1},
    {0x1.a5a014347406cp-55, -0x1.073879922ffeep-1},
    {0x1.2ec1fc1b776b8p-60, -0x1.f8ba4dbf89abap-2},
    {-0x1.e0d891d3c6841p-58, -0x1.e2b5d3806f63bp-2},
    {-0x1.6850e59c37f8fp-58, -0x1.cc66e9931c45ep-2},
    {-0x1.5b362cb974183p-57, -0x1.b5d1009e15ccp-2},
    {-0x1.6da81290bdbabp-57, -0x1.9ef7943a8ed8ap-2},
    {0x1.72cedd3d5a61p-57, -0x1.87de2a6aea963p-2},
    {0x1.44b19e0864c5dp-56, -0x1.7088530fa459fp-2},
    {0x1.efdc0d58cf62p-62, -0x1.58f9a75ab1fddp-2},
    {-0x1.0c97c4afa2518p-56, -0x1.4135c94176601p-2},
    {0x1.5d28da2c4612dp-56, -0x1.294062ed59f06p-2},
    {-0x1.824c20ab7aa9ap-56, -0x1.111d262b1f677p-2},
    {0x1.42deef11da2c4p-57, -0x1.f19f97b215f1bp-3},
    {0x1.af1439e521935p-62, -0x1.c0b826a7e4f63p-3},
    {0x1.26d19b9ff8d82p-57, -0x1.8f8b83c69a60bp-3},
    {-0x1.531ff779ddac6p-57, -0x1.5e214448b3fc6p-3},
    {-0x1.13000a89a11ep-58, -0x1.2c8106e8e613ap-3},
    {-0x1.a2704729ae56dp-59, -0x1.f564e56a9730ep-4},
    {0x1.e2718d26ed688p-60, -0x1.917a6bc29b42cp-4},
    {0x1.9a088a8bf6b2cp-59, -0x1.2d52092ce19f6p-4},
    {0x1.912bd0d569a9p-61, -0x1.91f65f10dd814p-5},
    {0x1.b1d63091a013p-64, -0x1.92155f7a3667ep-6},
};

// from sincos_eval.h
static inline double sincos_eval(const DoubleDouble u, DoubleDouble *sin_u,
                               DoubleDouble *cos_u) 
{
  // Evaluate sin(y) = sin(x - k * (pi/128))
  // We use the degree-7 Taylor approximation:
  //   sin(y) ~ y - y^3/3! + y^5/5! - y^7/7!
  // Then the error is bounded by:
  //   |sin(y) - (y - y^3/3! + y^5/5! - y^7/7!)| < |y|^9/9! < 2^-54/9! < 2^-72.
  // For y ~ u_hi + u_lo, fully expanding the polynomial and drop any terms
  // < ulp(u_hi^3) gives us:
  //   y - y^3/3! + y^5/5! - y^7/7! = ...
  // ~ u_hi + u_hi^3 * (-1/6 + u_hi^2 * (1/120 - u_hi^2 * 1/5040)) +
  //        + u_lo (1 + u_hi^2 * (-1/2 + u_hi^2 / 24))
  double u_hi_sq = u.hi * u.hi; // Error < ulp(u_hi^2) < 2^(-6 - 52) = 2^-58.
  // p1 ~ 1/120 + u_hi^2 / 5040.
  double p1 = multiply_add(u_hi_sq, -0x1.a01a01a01a01ap-13,
                                   0x1.1111111111111p-7);
  // q1 ~ -1/2 + u_hi^2 / 24.
  double q1 = multiply_add(u_hi_sq, 0x1.5555555555555p-5, -0x1.0p-1);
  double u_hi_3 = u_hi_sq * u.hi;
  // p2 ~ -1/6 + u_hi^2 (1/120 - u_hi^2 * 1/5040)
  double p2 = multiply_add(u_hi_sq, p1, -0x1.5555555555555p-3);
  // q2 ~ 1 + u_hi^2 (-1/2 + u_hi^2 / 24)
  double q2 = multiply_add(u_hi_sq, q1, 1.0);
  double sin_lo = multiply_add(u_hi_3, p2, u.lo * q2);
  // Overall, |sin(y) - (u_hi + sin_lo)| < 2*ulp(u_hi^3) < 2^-69.

  // Evaluate cos(y) = cos(x - k * (pi/128))
  // We use the degree-8 Taylor approximation:
  //   cos(y) ~ 1 - y^2/2 + y^4/4! - y^6/6! + y^8/8!
  // Then the error is bounded by:
  //   |cos(y) - (...)| < |y|^10/10! < 2^-81
  // For y ~ u_hi + u_lo, fully expanding the polynomial and drop any terms
  // < ulp(u_hi^3) gives us:
  //   1 - y^2/2 + y^4/4! - y^6/6! + y^8/8! = ...
  // ~ 1 - u_hi^2/2 + u_hi^4(1/24 + u_hi^2 (-1/720 + u_hi^2/40320)) +
  //     + u_hi u_lo (-1 + u_hi^2/6)
  // We compute 1 - u_hi^2 accurately:
  //   v_hi + v_lo ~ 1 - u_hi^2/2
  // with error <= 2^-105.
  double u_hi_neg_half = (-0.5) * u.hi;
  DoubleDouble v;
  // these next 2 must be FMA
  v.hi = __builtin_fma(u.hi, u_hi_neg_half, 1.0);
  v.lo = 1.0 - v.hi; // Exact
  v.lo = __builtin_fma(u.hi, u_hi_neg_half, v.lo);

  // r1 ~ -1/720 + u_hi^2 / 40320
  double r1 = multiply_add(u_hi_sq, 0x1.a01a01a01a01ap-16,
                                   -0x1.6c16c16c16c17p-10);
  // s1 ~ -1 + u_hi^2 / 6
  double s1 = multiply_add(u_hi_sq, 0x1.5555555555555p-3, -1.0);
  double u_hi_4 = u_hi_sq * u_hi_sq;
  double u_hi_u_lo = u.hi * u.lo;
  // r2 ~ 1/24 + u_hi^2 (-1/720 + u_hi^2 / 40320)
  double r2 = multiply_add(u_hi_sq, r1, 0x1.5555555555555p-5);
  // s2 ~ v_lo + u_hi * u_lo * (-1 + u_hi^2 / 6)
  double s2 = multiply_add(u_hi_u_lo, s1, v.lo);
  double cos_lo = multiply_add(u_hi_4, r2, s2);
  // Overall, |cos(y) - (v_hi + cos_lo)| < 2*ulp(u_hi^4) < 2^-75.

  DoubleDouble r_sin_u = exact_add(u.hi, sin_lo);
  DoubleDouble r_cos_u = exact_add(v.hi, cos_lo);
  *sin_u=r_sin_u;
  *cos_u=r_cos_u;

  return multiply_add(__builtin_fabs(u_hi_3),0x1.0p-51, 0x1.0p-105);
}


// from sin.h

double ya_sin(double x) {
  b64u64_u t = {.f = x};
  uint16_t x_e = (t.u >> DOUBLE_MANTISSA_BITS) & 0x7ff; /* get biased exponent */

  DoubleDouble y;
  unsigned k = 0;
#ifdef INSTRUMENT		
  scount_all++; // total count
#endif  
  // |x| < 2^16
  if (ya_LIKELY(x_e < DOUBLE_EXPONENT_BIAS + FAST_PASS_EXPONENT)) {
    // |x| < 2^-4
    if (ya_UNLIKELY(x_e < DOUBLE_EXPONENT_BIAS - 4)) {
      // |x| < 2^-26, |sin(x) - x| < ulp(x)/2.
      if (ya_UNLIKELY(x_e < DOUBLE_EXPONENT_BIAS - 26)) {
        if (ya_UNLIKELY(x == 0.0))  // Signed zeros.
          return x + x; // Make sure it works with FTZ/DAZ.
	   // Taylor expansion of sin(x) is x - x^3/6 around zero
		// for x=-0, fma (x, -0x1p-54, x) returns +0
		/* We have underflow when 0 < |x| < 2^-1022 or when |x| = 2^-1022
		   and rounding towards zero. */
		double res = __builtin_fma (x, -0x1p-54, x);// must be FMA
#ifdef CORE_MATH_SUPPORT_ERRNO
		if (__builtin_fabs (x) < 0x1p-1022 || __builtin_fabs (res) < 0x1p-1022)
		  errno = ERANGE; // underflow
#endif
		return res;
      }
      // No range reduction needed.

      // Use degree-9 polynomial approximation:
      //   sin(x) ~ x + a1 * x^3 + a2 * x^5 + a3 * x^7 + a4 * x^9
      //          ~ x + x^3 * Q(x^2).
      // > P = fpminimax(sin(x)/x, [|0, 2, 4, 6, 8|], [|1, D...|], [0, 2^-4]);
      // > dirtyinfnorm((sin(x) - x*P)/sin(x), [-2^-4, 2^-4]);
      // 0x1.3c2e...p-69
      // > P;
      const double COEFFS[] = {-0x1.5555555555555p-3, 0x1.111111110f491p-7,
                                   -0x1.a01a00e16af3ep-13,
                                   0x1.71c24233f1bafp-19};
      double x_sq = x * x;
      double c0 =  multiply_add(x_sq, COEFFS[1], COEFFS[0]);
      double c1 =  multiply_add(x_sq, COEFFS[3], COEFFS[2]);
      double x4 = x_sq * x_sq;
      double x3 = x * x_sq;
      double r_lo =  multiply_add(x4, c1, c0) * x3;

      // Overall errors <= 2 * ulp(x^3/6) + |x| * 2^-68.
      double err =  multiply_add(x_sq, 0x1.0p-53, 0x1.0p-68);
      double r_lo_u =  multiply_add(x, err, r_lo);
      double r_lo_l =  multiply_add(-x, err, r_lo);
      double r_upper = x + r_lo_u;
      double r_lower = x + r_lo_l;

      if (ya_LIKELY(r_upper == r_lower))
		{
#ifdef INSTRUMENT		
		 scount1++;
#endif		
         return r_upper;
		}        

      //k = range_reduction_small(x, y);
      //return sin_accurate(x, x_e, k, range_reduction_large);
#ifdef INSTRUMENT		
	  scount2++;
#endif      
	  double _cr_sin_accurate(double x);// in cr_sin.c
      return _cr_sin_accurate(x);
    } else {
      // Small range reduction.
#ifdef INSTRUMENT		
	  scount3++;
#endif      
      k = range_reduction_small(x, &y);
    }
  } else {
  if (__builtin_expect (x_e == 0x7ff, 0)) /* NaN, +Inf and -Inf. */
    {
      if ((t.u << 1) == 0x7ffull<<53){ // +/-Inf
#ifdef CORE_MATH_SUPPORT_ERRNO
        errno = EDOM;
#endif
        return x - x; // raises invalid
      }
      return x + x;
    }	  

    // Large range reduction.
#ifdef INSTRUMENT		
	scount4++;
#endif    
    k = LargeRangeReduction_fast(x, &y);
  }

  DoubleDouble sin_y, cos_y;

  double err=sincos_eval(y, &sin_y, &cos_y);

  // Look up sin(k * pi/128) and cos(k * pi/128)
  // Fast look up version, but needs 256-entry table.
  // cos(k * pi/128) = sin(k * pi/128 + pi/2) = sin((k + 64) * pi/128).
  DoubleDouble sin_k = SIN_K_PI_OVER_128[k & 255];
  DoubleDouble cos_k = SIN_K_PI_OVER_128[(k + 64) & 255];

  // After range reduction, k = round(x * 128 / pi) and y = x - k * (pi / 128).
  // So k is an integer and -pi / 256 <= y <= pi / 256.
  // Then sin(x) = sin((k * pi/128 + y)
  //             = sin(y) * cos(k*pi/128) + cos(y) * sin(k*pi/128)
  DoubleDouble sin_k_cos_y = mult_dd_dd(cos_y, sin_k);
  DoubleDouble cos_k_sin_y = mult_dd_dd(sin_y, cos_k);
  // When k != 0 mod 128,
  //   |sin( k * pi/128 )| > pi/128 - epsilon > |y| >= |sin(y)|,
  // and cos(y) > 1 - pi/128.  So we can use Fast2Sum for the addition:
  //   sin(y) * cos(k*pi/128) + cos(y) * sin(k*pi/128).
  DoubleDouble rr = exact_add(sin_k_cos_y.hi, cos_k_sin_y.hi);
  rr.lo += sin_k_cos_y.lo + cos_k_sin_y.lo;


  // Accurate test and pass for correctly rounded implementation.

  double rlp = rr.lo + err;
  double rlm = rr.lo - err;

  double r_upper = rr.hi + rlp; // (rr.lo + ERR);
  double r_lower = rr.hi + rlm; // (rr.lo - ERR);

  // Ziv's rounding test.
  if (ya_LIKELY(r_upper == r_lower))
  	{
#ifdef INSTRUMENT		
	 scount5++;
#endif
     return r_upper;
	}
  //return sin_accurate(x, x_e, k, range_reduction_large);
#ifdef INSTRUMENT		
  scount6++;
#endif  
  double _cr_sin_accurate (double x);// in cr_sin.c
  return _cr_sin_accurate(x);
}

// from cos.h

double ya_cos(double x) {
  b64u64_u t = {.f = x};
  uint16_t x_e = (t.u >> DOUBLE_MANTISSA_BITS) & 0x7ff; /* get biased exponent */

  DoubleDouble y;
  unsigned k = 0;
#ifdef INSTRUMENT		
		 ccount_all++; // total count
#endif  
  // |x| < 2^16
  if (ya_LIKELY(x_e < DOUBLE_EXPONENT_BIAS + FAST_PASS_EXPONENT)) {
    // |x| < 2^-4
    if (ya_UNLIKELY(x_e < DOUBLE_EXPONENT_BIAS - 4)) {
      // |x| < 2^-27
      if (ya_UNLIKELY(x_e < DOUBLE_EXPONENT_BIAS - 27)) {
        // Signed zeros.
        if (ya_UNLIKELY(x == 0.0))
          return 1.0; 
		// For |x| < 2^-27, |cos(x) - 1| < |x|^2/2 < 2^-54 = ulp(1 - 2^-53)/2.
        return 1.0;
      }
      // No range reduction needed.

	   // Use degree-8 polynomial approximation:
      //   cos(x) ~ 1 + a1 * x^2 + a2 * x^4 + a3 * x^6 + a4 * x^8
      //          ~ 1 + x^2 * Q(x^2).
      // > P = fpminimax(cos(x), [|0, 2, 4, 6, 8|], [|1, D...|], [0, 2^-4]);
      // > dirtyinfnorm(cos(x) - P, [-2^-4, 2^-4]);
      // 0x1.3cfe...p-70
      // > P;
      const double COEFFS[] = {-0x1p-1, 0x1.5555555555262p-5,
                                   -0x1.6c16c1508bff1p-10,
                                   0x1.a00ffd769159ap-16};
      double x_sq = x * x;
      double c0 =  multiply_add(x_sq, COEFFS[1], COEFFS[0]);
      double c1 =  multiply_add(x_sq, COEFFS[3], COEFFS[2]);
      double x4 = x_sq * x_sq;

      double r_lo =  multiply_add(x4, c1, c0) * x_sq;

      // Overall errors <= ulp(x^2/2) + 2^-69.
      double err = multiply_add(x_sq, 0x1.0p-53, 0x1.0p-69);
      double r_lo_u = r_lo + err;
      double r_lo_l = r_lo - err;
      double r_upper = 1.0 + r_lo_u;
      double r_lower = 1.0 + r_lo_l;

      if (ya_LIKELY(r_upper == r_lower))
		{
#ifdef INSTRUMENT		
		 ccount1++;
#endif		
         return r_upper;
		}        

      //k = range_reduction_small(x, y);
      //return cos_accurate(x, x_e, k, range_reduction_large);
#ifdef INSTRUMENT		
	  ccount2++;
#endif  
	  double _cr_cos_accurate(double x) ; // in cr_cos.c   
	  t.u &= 0x7fffffffffffffff;// cos(-x)=cos(x) , as _cr_cos_accurate(x) needs positive x 
      return _cr_cos_accurate(t.f);
    } else {
      // Small range reduction.
#ifdef INSTRUMENT		
	  ccount3++;
#endif      
      k = range_reduction_small(x, &y);
    }
  } else {
  if (__builtin_expect (x_e == 0x7ff, 0)) /* NaN, +Inf and -Inf. */
    {
#ifdef CORE_MATH_SUPPORT_ERRNO
      if ((t.u << 1) == 0x7ffull<<53) // Inf
		{errno = EDOM;
		 return 0.0 / 0.0; // raise invalid flag
		}
#else
      if ((t.u << 1) == 0x7ffull<<53) // Inf
        return 0.0 / 0.0; // raise invalid flag
#endif
      return x + x; // return qNaN
    }	  

    // Large range reduction.
#ifdef INSTRUMENT		
	ccount4++;
#endif    
    k = LargeRangeReduction_fast(x, &y);
  }

  DoubleDouble sin_y, cos_y;

  double err=sincos_eval(y, &sin_y, &cos_y);

  // Look up sin(k * pi/128) and cos(k * pi/128)
  // Fast look up version, but needs 256-entry table.
   // -sin(k * pi/128) = sin((k + 128) * pi/128)
  // cos(k * pi/128) = sin(k * pi/128 + pi/2) = sin((k + 64) * pi/128).
  DoubleDouble msin_k = SIN_K_PI_OVER_128[(k + 128) & 255];
  DoubleDouble cos_k = SIN_K_PI_OVER_128[(k + 64) & 255];

  // After range reduction, k = round(x * 128 / pi) and y = x - k * (pi / 128).
  // So k is an integer and -pi / 256 <= y <= pi / 256.
  // Then cos(x) = cos((k * pi/128 + y)
  //             = cos(y) * cos(k*pi/128) - sin(y) * sin(k*pi/128)
  DoubleDouble cos_k_cos_y  = mult_dd_dd(cos_y, cos_k);
  DoubleDouble msin_k_sin_y = mult_dd_dd(sin_y, msin_k);
  // When k != 64 mod 128,
  //   |cos( k * pi/128 )| > pi/128 - epsilon > |y| >= |sin(y)|,
  // and cos(y) > 1 - pi/128.  So we can use Fast2Sum for the subtraction:
  //   cos(y) * cos(k*pi/128) - sin(y) * sin(k*pi/128).
  DoubleDouble rr = exact_add(cos_k_cos_y.hi, msin_k_sin_y.hi);
  rr.lo += msin_k_sin_y.lo + cos_k_cos_y.lo;  

  // Accurate test and pass for correctly rounded implementation.

  double rlp = rr.lo + err;
  double rlm = rr.lo - err;

  double r_upper = rr.hi + rlp; // (rr.lo + ERR);
  double r_lower = rr.hi + rlm; // (rr.lo - ERR);

  // Ziv's rounding test.
  if (ya_LIKELY(r_upper == r_lower))
  	{
#ifdef INSTRUMENT		
	 ccount5++;
#endif
     return r_upper;
	}
  //return cos_accurate(x, x_e, k, range_reduction_large);
#ifdef INSTRUMENT		
  ccount6++;
#endif  
  double _cr_cos_accurate(double x) ; // in cr_cos.c   
  t.u &= 0x7fffffffffffffff;// cos(-x)=cos(x) , as _cr_cos_accurate(x) needs positive x 
  return _cr_cos_accurate(t.f);  
}