/* ya_exp.c
  ======

   Uses basically the same algorithms as the exp function in llvm-project-23.1.1 (\llvm-project-23.1.1\libc\src\__support\FPUtil & \llvm-project-23.1.1\libc\src\__support\math ).
   Note the llvm original is in C++, whereas this is in C.
   Where practical, variable names (and many comments) have been kept the same in this code as the llvm code, to make them easier to compare.
   It leverages the double-double maths library ya_double_double.   
   
   This code assumes double FMA is available (either as a single instruction [ideally], or emulated in software) 
 
   The "accurate" evaluation is done using using functions in cr_exp() - in the test program these are only used 0.001% of the time so it was not worthwhile converting the accurate algorithms
   used in the llvm implementations(which might or might not be faster, but will make no difference to the average timing and negligible difference to the worse case).

   Checking for special cases (inf, NAN, etc) is done in the same was as cr_exp.c, as that solution was found to be faster than the approach used in llvm (& supports CORE_MATH_SUPPORT_ERRNO).
   
 	According to  https://members.loria.fr/PZimmermann/papers/accuracy.pdf the llvm double exp function has a maximum error of 0.500ulp,
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

// #define USE_INTERNAL_DD_POLY /* if defined use internal double-double polynomial, otherwise just call cr_exp() directly after initial double poly approximation */
							 /* Warning: using USE_INTERNAL_DD_POLY is slower than using just calling cr_exp , at least on an Intel u3-10100 processor, with a 64 bit compilation  */

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

#if 0 /* debug */
 #include <stdio.h>
 #include "../my_printf/my_printf.h"
 #define dprintf(...) my_printf(__VA_ARGS__)
#else
 #define dprintf(x,...) /* nothing */
#endif

struct _TripleDouble
{	
	/* Note lo is first to match llvm - so smallest magnitude value is 1st in initialisers  */
	double lo; 
	double mid; 
	double hi;
	
};
typedef struct _TripleDouble TripleDouble;/* 3 doubles */

// Lookup table for 2^(k * 2^-6) with k = 0..63.
// Generated by Sollya with:
// > display=hexadecimal;
// > prec = 500;
// > for i from 0 to 63 do {
//     a = 2^(i * 2^-6);
//     b = round(a, D, RN);
//     c = round(a - b, D, RN);
//     d = round(a - b - c, D, RN);
//     print("{", d, ",", c, ",", b, "},");
//   };
static const TripleDouble EXP2_MID1[64] = {
    {0, 0, 0x1p0},
    {-0x1.9085b0a3d74d5p-110, -0x1.19083535b085dp-56, 0x1.02c9a3e778061p0},
    {0x1.05ff94f8d257ep-110, 0x1.d73e2a475b465p-55, 0x1.059b0d3158574p0},
    {0x1.15820d96b414fp-111, 0x1.186be4bb284ffp-57, 0x1.0874518759bc8p0},
    {-0x1.67c9bd6ebf74cp-108, 0x1.8a62e4adc610bp-54, 0x1.0b5586cf9890fp0},
    {-0x1.5aa76994e9ddbp-113, 0x1.03a1727c57b53p-59, 0x1.0e3ec32d3d1a2p0},
    {0x1.9d58b988f562dp-109, -0x1.6c51039449b3ap-54, 0x1.11301d0125b51p0},
    {-0x1.2fe7bb4c76416p-108, -0x1.32fbf9af1369ep-54, 0x1.1429aaea92dep0},
    {0x1.4f2406aa13ffp-109, -0x1.19041b9d78a76p-55, 0x1.172b83c7d517bp0},
    {0x1.ad36183926ae8p-111, 0x1.e5b4c7b4968e4p-55, 0x1.1a35beb6fcb75p0},
    {0x1.ea62d0881b918p-110, 0x1.e016e00a2643cp-54, 0x1.1d4873168b9aap0},
    {-0x1.781dbc16f1ea4p-111, 0x1.dc775814a8495p-55, 0x1.2063b88628cd6p0},
    {-0x1.4d89f9af532ep-109, 0x1.9b07eb6c70573p-54, 0x1.2387a6e756238p0},
    {0x1.277393a461b77p-110, 0x1.2bd339940e9d9p-55, 0x1.26b4565e27cddp0},
    {0x1.de5448560469p-111, 0x1.612e8afad1255p-55, 0x1.29e9df51fdee1p0},
    {-0x1.ee9d8f8cb9307p-110, 0x1.0024754db41d5p-54, 0x1.2d285a6e4030bp0},
    {0x1.7b7b2f09cd0d9p-110, 0x1.6f46ad23182e4p-55, 0x1.306fe0a31b715p0},
    {-0x1.406a2ea6cfc6bp-108, 0x1.32721843659a6p-54, 0x1.33c08b26416ffp0},
    {0x1.87e3e12516bfap-108, -0x1.63aeabf42eae2p-54, 0x1.371a7373aa9cbp0},
    {0x1.9b0b1ff17c296p-111, -0x1.5e436d661f5e3p-56, 0x1.3a7db34e59ff7p0},
    {-0x1.808ba68fa8fb7p-109, 0x1.ada0911f09ebcp-55, 0x1.3dea64c123422p0},
    {-0x1.32b43eafc6518p-114, -0x1.ef3691c309278p-58, 0x1.4160a21f72e2ap0},
    {-0x1.0ac312de3d922p-114, 0x1.89b7a04ef80dp-59, 0x1.44e086061892dp0},
    {0x1.e1eebae743acp-111, 0x1.3c1a3b69062fp-56, 0x1.486a2b5c13cdp0},
    {0x1.c06c7745c2b39p-113, 0x1.d4397afec42e2p-56, 0x1.4bfdad5362a27p0},
    {-0x1.1aa1fd7b685cdp-112, -0x1.4b309d25957e3p-54, 0x1.4f9b2769d2ca7p0},
    {0x1.fa733951f214cp-111, -0x1.07abe1db13cadp-55, 0x1.5342b569d4f82p0},
    {-0x1.ff86852a613ffp-111, 0x1.9bb2c011d93adp-54, 0x1.56f4736b527dap0},
    {-0x1.744ee506fdafep-109, 0x1.6324c054647adp-54, 0x1.5ab07dd485429p0},
    {-0x1.95f9ab75fa7d6p-108, 0x1.ba6f93080e65ep-54, 0x1.5e76f15ad2148p0},
    {0x1.5d8e757cfb991p-111, -0x1.383c17e40b497p-54, 0x1.6247eb03a5585p0},
    {0x1.4a337f4dc0a3bp-108, -0x1.bb60987591c34p-54, 0x1.6623882552225p0},
    {0x1.57d3e3adec175p-108, -0x1.bdd3413b26456p-54, 0x1.6a09e667f3bcdp0},
    {0x1.a59f88abbe778p-115, -0x1.bbe3a683c88abp-57, 0x1.6dfb23c651a2fp0},
    {-0x1.269796953a4c3p-109, -0x1.16e4786887a99p-55, 0x1.71f75e8ec5f74p0},
    {-0x1.8f8e7fa19e5e8p-108, -0x1.0245957316dd3p-54, 0x1.75feb564267c9p0},
    {-0x1.4217a932d10d4p-113, -0x1.41577ee04992fp-55, 0x1.7a11473eb0187p0},
    {0x1.70a1427f8fcdfp-112, 0x1.05d02ba15797ep-56, 0x1.7e2f336cf4e62p0},
    {0x1.0f6ad65cbbac1p-112, -0x1.d4c1dd41532d8p-54, 0x1.82589994cce13p0},
    {-0x1.f16f65181d921p-109, -0x1.fc6f89bd4f6bap-54, 0x1.868d99b4492edp0},
    {-0x1.30644a7836333p-110, 0x1.6e9f156864b27p-54, 0x1.8ace5422aa0dbp0},
    {0x1.3bf26d2b85163p-114, 0x1.5cc13a2e3976cp-55, 0x1.8f1ae99157736p0},
    {0x1.697e257ac0db2p-111, -0x1.75fc781b57ebcp-57, 0x1.93737b0cdc5e5p0},
    {0x1.7edb9d7144b6fp-108, -0x1.d185b7c1b85d1p-54, 0x1.97d829fde4e5p0},
    {0x1.6376b7943085cp-110, 0x1.c7c46b071f2bep-56, 0x1.9c49182a3f09p0},
    {0x1.354084551b4fbp-109, -0x1.359495d1cd533p-54, 0x1.a0c667b5de565p0},
    {-0x1.bfd7adfd63f48p-111, -0x1.d2f6edb8d41e1p-54, 0x1.a5503b23e255dp0},
    {0x1.8b16ae39e8cb9p-109, 0x1.0fac90ef7fd31p-54, 0x1.a9e6b5579fdbfp0},
    {0x1.a7fbc3ae675eap-108, 0x1.7a1cd345dcc81p-54, 0x1.ae89f995ad3adp0},
    {0x1.2babc0edda4d9p-111, -0x1.2805e3084d708p-57, 0x1.b33a2b84f15fbp0},
    {0x1.aa64481e1ab72p-111, -0x1.5584f7e54ac3bp-56, 0x1.b7f76f2fb5e47p0},
    {0x1.9a164050e1258p-109, 0x1.23dd07a2d9e84p-55, 0x1.bcc1e904bc1d2p0},
    {0x1.99e51125928dap-110, 0x1.11065895048ddp-55, 0x1.c199bdd85529cp0},
    {-0x1.fc44c329d5cb2p-109, 0x1.2884dff483cadp-54, 0x1.c67f12e57d14bp0},
    {0x1.d8765566b032ep-110, 0x1.503cbd1e949dbp-56, 0x1.cb720dcef9069p0},
    {-0x1.e7044039da0f6p-108, -0x1.cbc3743797a9cp-54, 0x1.d072d4a07897cp0},
    {-0x1.ab053b05531fcp-111, 0x1.2ed02d75b3707p-55, 0x1.d5818dcfba487p0},
    {0x1.7f6246f0ec615p-108, 0x1.c2300696db532p-54, 0x1.da9e603db3285p0},
    {0x1.b7225a944efd6p-108, -0x1.1a5cd4f184b5cp-54, 0x1.dfc97337b9b5fp0},
    {0x1.1e92cb3c2d278p-109, 0x1.39e8980a9cc8fp-55, 0x1.e502ee78b3ff6p0},
    {-0x1.fc0f242bbf3dep-109, -0x1.e9c23179c2893p-54, 0x1.ea4afa2a490dap0},
    {0x1.f6dd5d229ff69p-108, 0x1.dc7f486a4b6bp-54, 0x1.efa1bee615a27p0},
    {-0x1.4019bffc80ef3p-110, 0x1.9d3e12dd8a18bp-54, 0x1.f50765b6e454p0},
    {0x1.dc060c36f7651p-112, 0x1.74853f3a5931ep-55, 0x1.fa7c1819e90d8p0},
};

// Lookup table for 2^(k * 2^-12) with k = 0..63.
// Generated by Sollya with:
// > display=hexadecimal;
// > prec = 500;
// > for i from 0 to 63 do {
//     a = 2^(i * 2^-12);
//     b = round(a, D, RN);
//     c = round(a - b, D, RN);
//     d = round(a - b - c, D, RN);
//     print("{", d, ",", c, ",", b, "},");
//   };
static const TripleDouble EXP2_MID2[64] = {
    {0, 0, 0x1p0},
    {0x1.39726694630e3p-108, 0x1.ae8e38c59c72ap-54, 0x1.000b175effdc7p0},
    {0x1.e5e06ddd31156p-112, -0x1.7b5d0d58ea8f4p-58, 0x1.00162f3904052p0},
    {0x1.5a0768b51f609p-111, 0x1.4115cb6b16a8ep-54, 0x1.0021478e11ce6p0},
    {0x1.d008403605217p-111, -0x1.d7c96f201bb2fp-55, 0x1.002c605e2e8cfp0},
    {0x1.89bc16f765708p-109, 0x1.84711d4c35e9fp-54, 0x1.003779a95f959p0},
    {-0x1.4535b7f8c1e2dp-109, -0x1.0484245243777p-55, 0x1.0042936faa3d8p0},
    {-0x1.8ba92f6b25456p-108, -0x1.4b237da2025f9p-54, 0x1.004dadb113dap0},
    {-0x1.30c72e81f4294p-113, -0x1.5e00e62d6b30dp-56, 0x1.0058c86da1c0ap0},
    {-0x1.34a5384e6f0b9p-110, 0x1.a1d6cedbb9481p-54, 0x1.0063e3a559473p0},
    {0x1.f8d0580865d2ep-108, -0x1.4acf197a00142p-54, 0x1.006eff583fc3dp0},
    {-0x1.002bcb3ae9a99p-111, -0x1.eaf2ea42391a5p-57, 0x1.007a1b865a8cap0},
    {0x1.c3c5aedee9851p-111, 0x1.da93f90835f75p-56, 0x1.0085382faef83p0},
    {0x1.7217851d1ec6ep-109, -0x1.6a79084ab093cp-55, 0x1.00905554425d4p0},
    {-0x1.80cbca335a7c3p-110, 0x1.86364f8fbe8f8p-54, 0x1.009b72f41a12bp0},
    {-0x1.706bd4eb22595p-110, -0x1.82e8e14e3110ep-55, 0x1.00a6910f3b6fdp0},
    {-0x1.b55dd523f3c08p-111, -0x1.4f6b2a7609f71p-55, 0x1.00b1afa5abcbfp0},
    {0x1.90a1e207cced1p-110, -0x1.e1a258ea8f71bp-56, 0x1.00bcceb7707ecp0},
    {0x1.78d0472db37c5p-110, 0x1.4362ca5bc26f1p-56, 0x1.00c7ee448ee02p0},
    {-0x1.bcd4db3cb52fep-109, 0x1.095a56c919d02p-54, 0x1.00d30e4d0c483p0},
    {-0x1.cf1b131575ec2p-112, -0x1.406ac4e81a645p-57, 0x1.00de2ed0ee0f5p0},
    {-0x1.6aaa1fa7ff913p-112, 0x1.b5a6902767e09p-54, 0x1.00e94fd0398ep0},
    {0x1.68f236dff3218p-110, -0x1.91b2060859321p-54, 0x1.00f4714af41d3p0},
    {-0x1.e8bb58067e60ap-109, 0x1.427068ab22306p-55, 0x1.00ff93412315cp0},
    {0x1.d4cd5e1d71fdfp-108, 0x1.c1d0660524e08p-54, 0x1.010ab5b2cbd11p0},
    {0x1.e4ecf350ebe88p-108, -0x1.e7bdfb3204be8p-54, 0x1.0115d89ff3a8bp0},
    {0x1.6a2aa2c89c4f8p-109, 0x1.843aa8b9cbbc6p-55, 0x1.0120fc089ff63p0},
    {0x1.1ca368a20ed05p-110, -0x1.34104ee7edae9p-56, 0x1.012c1fecd613bp0},
    {0x1.edb1095d925cfp-114, -0x1.2b6aeb6176892p-56, 0x1.0137444c9b5b5p0},
    {-0x1.488c78eded75fp-111, 0x1.a8cd33b8a1bb3p-56, 0x1.01426927f5278p0},
    {-0x1.7480f5ea1b3c9p-113, 0x1.2edc08e5da99ap-56, 0x1.014d8e7ee8d2fp0},
    {-0x1.ae45989a04dd5p-111, 0x1.57ba2dc7e0c73p-55, 0x1.0158b4517bb88p0},
    {0x1.bf48007d80987p-109, 0x1.b61299ab8cdb7p-54, 0x1.0163da9fb3335p0},
    {0x1.1aa91a059292cp-109, -0x1.90565902c5f44p-54, 0x1.016f0169949edp0},
    {0x1.b6663292855f5p-110, 0x1.70fc41c5c2d53p-55, 0x1.017a28af25567p0},
    {0x1.e7fbca6793d94p-108, 0x1.4b9a6e145d76cp-54, 0x1.018550706ab62p0},
    {-0x1.5b9f5c7de3b93p-110, -0x1.008eff5142bf9p-56, 0x1.019078ad6a19fp0},
    {0x1.4638bf2f6acabp-110, -0x1.77669f033c7dep-54, 0x1.019ba16628de2p0},
    {-0x1.ab237b9a069c5p-109, -0x1.09bb78eeead0ap-54, 0x1.01a6ca9aac5f3p0},
    {0x1.3ab358be97cefp-108, 0x1.371231477ece5p-54, 0x1.01b1f44af9f9ep0},
    {-0x1.4027b2294bb64p-110, 0x1.5e7626621eb5bp-56, 0x1.01bd1e77170b4p0},
    {0x1.656394426c99p-111, -0x1.bc72b100828a5p-54, 0x1.01c8491f08f08p0},
    {0x1.bf9785189bdd8p-111, -0x1.ce39cbbab8bbep-57, 0x1.01d37442d507p0},
    {0x1.7c12f86114fe3p-109, 0x1.16996709da2e2p-55, 0x1.01de9fe280ac8p0},
    {-0x1.653d5d24b5d28p-109, -0x1.c11f5239bf535p-55, 0x1.01e9cbfe113efp0},
    {0x1.04a0cdc1d86d7p-109, 0x1.e1d4eb5edc6b3p-55, 0x1.01f4f8958c1c6p0},
    {0x1.c678c46149782p-109, -0x1.afb99946ee3fp-54, 0x1.020025a8f6a35p0},
    {0x1.48524e1e9df7p-108, -0x1.8f06d8a148a32p-54, 0x1.020b533856324p0},
    {0x1.9953ea727ff0bp-109, -0x1.2bf310fc54eb6p-55, 0x1.02168143b0281p0},
    {-0x1.ccfbbec22d28ep-108, -0x1.c95a035eb4175p-54, 0x1.0221afcb09e3ep0},
    {0x1.9e2bb6e181de1p-108, -0x1.491793e46834dp-54, 0x1.022cdece68c4fp0},
    {0x1.f17609ae29308p-110, -0x1.3e8d0d9c49091p-56, 0x1.02380e4dd22adp0},
    {-0x1.c7dc2c476bfb8p-110, -0x1.314aa16278aa3p-54, 0x1.02433e494b755p0},
    {-0x1.fab994971d4a3p-109, 0x1.48daf888e9651p-55, 0x1.024e6ec0da046p0},
    {0x1.848b62cbdd0afp-109, 0x1.56dc8046821f4p-55, 0x1.02599fb483385p0},
    {-0x1.bf603ba715d0cp-109, 0x1.45b42356b9d47p-54, 0x1.0264d1244c719p0},
    {0x1.89434e751e1aap-110, -0x1.082ef51b61d7ep-56, 0x1.027003103b10ep0},
    {-0x1.03b54fd64e8acp-110, 0x1.2106ed0920a34p-56, 0x1.027b357854772p0},
    {0x1.7785ea0acc486p-109, -0x1.fd4cf26ea5d0fp-54, 0x1.0286685c9e059p0},
    {-0x1.ce447fdb35ff9p-109, -0x1.09f8775e78084p-54, 0x1.02919bbd1d1d8p0},
    {0x1.5b884aab5642ap-112, 0x1.64cbba902ca27p-58, 0x1.029ccf99d720ap0},
    {-0x1.cfb3e46d7c1cp-108, 0x1.4383ef231d207p-54, 0x1.02a803f2d170dp0},
    {-0x1.0d40cee4b81afp-112, 0x1.4a47a505b3a47p-54, 0x1.02b338c811703p0},
    {0x1.6ae7d36d7c1f7p-109, 0x1.e47120223467fp-54, 0x1.02be6e199c811p0},
};


// log2(e)
static const double LOG2_E = 0x1.71547652b82fep+0;

// Error bounds:
// Errors when using double precision.
static const double EXP_ERR_D = 0x1.8p-63;

#ifdef USE_INTERNAL_DD_POLY
// Errors when using double-double precision.
static const double EXP_ERR_DD = 0x1.0p-99;
#endif

// -2^-12 * log(2)
// > a = -2^-12 * log(2);
// > b = round(a, 30, RN);
// > c = round(a - b, 30, RN);
// > d = round(a - b - c, D, RN);
// Errors < 1.5 * 2^-133
static const double MLOG_2_EXP2_M12_HI = -0x1.62e42ffp-13;
static const double MLOG_2_EXP2_M12_MID = 0x1.718432a1b0e26p-47;

#ifdef USE_INTERNAL_DD_POLY
static const double MLOG_2_EXP2_M12_MID_30 = 0x1.718432ap-47;
static const double MLOG_2_EXP2_M12_LO = 0x1.b0e2633fe0685p-79;
#endif

// Polynomial approximations with double precision:
// Return expm1(dx) / x ~ 1 + dx / 2 + dx^2 / 6 + dx^3 / 24.
// For |dx| < 2^-13 + 2^-30:
//   | output - expm1(dx) / dx | < 2^-51.
static double poly_approx_d(double dx) {
  // dx^2
  double dx2 = dx * dx;
  // c0 = 1 + dx / 2
  double c0 = multiply_add(dx, 0.5, 1.0);
  // c1 = 1/6 + dx / 24
  double c1 = multiply_add(dx, 0x1.5555555555555p-5, 0x1.5555555555555p-3);
  // p = dx^2 * c1 + c0 = 1 + dx / 2 + dx^2 / 6 + dx^3 / 24
  double p = multiply_add(dx2, c1, c0);
  return p;
}

#ifdef USE_INTERNAL_DD_POLY
// Polynomial approximation with double-double precision:
// Return exp(dx) ~ 1 + dx + dx^2 / 2 + ... + dx^6 / 720
// For |dx| < 2^-13 + 2^-30:
//   | output - exp(dx) | < 2^-101
static DoubleDouble poly_approx_dd(const DoubleDouble dx)
 {
  // Taylor polynomial.
  const DoubleDouble COEFFS[] = {
      {0, 0x1p0},                                      // 1
      {0, 0x1p0},                                      // 1
      {0, 0x1p-1},                                     // 1/2
      {0x1.5555555555555p-57, 0x1.5555555555555p-3},   // 1/6
      {0x1.5555555555555p-59, 0x1.5555555555555p-5},   // 1/24
      {0x1.1111111111111p-63, 0x1.1111111111111p-7},   // 1/120
      {-0x1.f49f49f49f49fp-65, 0x1.6c16c16c16c17p-10}, // 1/720
  };
  // Evaluate polynomial using Horner's Scheme:
  // With polyeval(x, a_0, a_1, ..., a_n) = a_n * x^n + ... + a_1 * x + a_0, we
  // evaluated it as:  a_0 + x * (a_1 + x * ( ... (a_(n-1) + x * a_n) ... ) ) ).
  // We will use FMA instructions if available.
  // Example: to evaluate x^3 + 2*x^2 + 3*x + 4, call
  //   polyeval( x, 4.0, 3.0, 2.0, 1.0 )
#define d_mult_add(a,b,c) add_dd_dd(c,mult_dd_dd(a,b))
  //DoubleDouble p = polyeval(dx, COEFFS[0], COEFFS[1], COEFFS[2],
  //                                  COEFFS[3], COEFFS[4], COEFFS[5], COEFFS[6]);
  DoubleDouble p=COEFFS[6];
  for(int i=5;i>=0;--i)
  	p=d_mult_add(dx,p,COEFFS[i]);
  return p;
 }


// Compute exp(x) with double-double precision.
static DoubleDouble exp_double_double(double x, double kd,
                                           const DoubleDouble exp_mid) {
  // Recalculate dx:
  //   dx = x - k * 2^-12 * log(2)
  double t1 = multiply_add(kd, MLOG_2_EXP2_M12_HI, x); // exact
  double t2 = kd * MLOG_2_EXP2_M12_MID_30;                     // exact
  double t3 = kd * MLOG_2_EXP2_M12_LO;                         // Error < 2^-130

  DoubleDouble dx = exact_add(t1, t2);
  dx.lo += t3;

  // Degree-6 Taylor polynomial approximation in double-double precision.
  // | p - exp(x) | < 2^-100.
  DoubleDouble p = poly_approx_dd(dx);

  // Error bounds: 2^-99.
  DoubleDouble r = mult_dd_dd(exp_mid, p);

  return r;
}
#endif


double ya_exp(double x) 
{
  // code from cr_exp.c for special cases (this is faster than the approach used in llvm's exp.h
  b64u64_u ix = {.f = x};
  uint64_t aix = ix.u & (~(uint64_t)0>>1);
  // exp(x) rounds to 1 to nearest for |x| <= 0x1p-54
  if(__builtin_expect(aix <= 0x3c90000000000000ull, 0)) // |x| <= 0x1p-54
    return 1.0 + x;
  if(__builtin_expect(aix>=0x40862e42fefa39f0ull, 0)){ // |x| >= 0x1.62e42fefa39fp+9
    if(aix>0x7ff0000000000000ull) return x + x; // nan
    if(aix==0x7ff0000000000000ull){ // |x| = inf
      if(ix.u>>63)
	return 0.0; // x = -inf
      else
	return x; // x = inf
    }
    if(!(ix.u>>63)){ // x >= 0x1.62e42fefa39fp+9
#ifdef CORE_MATH_SUPPORT_ERRNO
      errno = ERANGE;
#endif
      volatile double z = 0x1p1023;
      return z*z;
    }
    if (aix>=0x40874910d52d3052ull) { // x <= -0x1.74910d52d3052p+9
#ifdef CORE_MATH_SUPPORT_ERRNO
      errno = ERANGE; // underflow
#endif
      volatile double z = 0x1p-1022;
      return z * z;
    }
  }
   
  // Now log(2^-1075) <= x <= -2^-53 or 2^-53 <= x < log(2^1023 * (2 - 2^-52))

  // Range reduction:
  // Let x = log(2) * (hi + mid1 + mid2) + lo
  // in which:
  //   hi is an integer
  //   mid1 * 2^6 is an integer
  //   mid2 * 2^12 is an integer
  // then:
  //   exp(x) = 2^hi * 2^(mid1) * 2^(mid2) * exp(lo).
  // With this formula:
  //   - multiplying by 2^hi is exact and cheap, simply by adding the exponent
  //     field.
  //   - 2^(mid1) and 2^(mid2) are stored in 2 x 64-element tables.
  //   - exp(lo) ~ 1 + lo + a0 * lo^2 + ...
  //
  // They can be defined by:
  //   hi + mid1 + mid2 = 2^(-12) * round(2^12 * log_2(e) * x)
  // If we store L2E = round(log2(e), D, RN), then:
  //   log2(e) - L2E ~ 1.5 * 2^(-56)
  // So the errors when computing in double precision is:
  //   | x * 2^12 * log_2(e) - D(x * 2^12 * L2E) | <=
  //  <= | x * 2^12 * log_2(e) - x * 2^12 * L2E | +
  //     + | x * 2^12 * L2E - D(x * 2^12 * L2E) |
  //  <= 2^12 * ( |x| * 1.5 * 2^-56 + eps(x))  for RN
  //     2^12 * ( |x| * 1.5 * 2^-56 + 2*eps(x)) for other rounding modes.
  // So if:
  //   hi + mid1 + mid2 = 2^(-12) * round(x * 2^12 * L2E) is computed entirely
  // in double precision, the reduced argument:
  //   lo = x - log(2) * (hi + mid1 + mid2) is bounded by:
  //   |lo| <= 2^-13 + (|x| * 1.5 * 2^-56 + 2*eps(x))
  //         < 2^-13 + (1.5 * 2^9 * 1.5 * 2^-56 + 2*2^(9 - 52))
  //         < 2^-13 + 2^-41
  //

  // The following trick computes the round(x * L2E) more efficiently
  // than using the rounding instructions, with the trade-off for less accuracy,
  // and hence a slightly larger range for the reduced argument `lo`.
  //
  // To be precise, since |x| < |log(2^-1075)| < 1.5 * 2^9,
  //   |x * 2^12 * L2E| < 1.5 * 2^9 * 1.5 < 2^23,
  // So we can fit the rounded result round(x * 2^12 * L2E) in int32_t.
  // Thus, the goal is to be able to use an additional addition and fixed width
  // shift to get an int32_t representing round(x * 2^12 * L2E).
  //
  // Assuming int32_t using 2-complement representation, since the mantissa part
  // of a double precision is unsigned with the leading bit hidden, if we add an
  // extra constant C = 2^e1 + 2^e2 with e1 > e2 >= 2^25 to the product, the
  // part that are < 2^e2 in resulted mantissa of (x*2^12*L2E + C) can be
  // considered as a proper 2-complement representations of x*2^12*L2E.
  //
  // One small problem with this approach is that the sum (x*2^12*L2E + C) in
  // double precision is rounded to the least significant bit of the dominant
  // factor C.  In order to minimize the rounding errors from this addition, we
  // want to minimize e1.  Another constraint that we want is that after
  // shifting the mantissa so that the least significant bit of int32_t
  // corresponds to the unit bit of (x*2^12*L2E), the sign is correct without
  // any adjustment.  So combining these 2 requirements, we can choose
  //   C = 2^33 + 2^32, so that the sign bit corresponds to 2^31 bit, and hence
  // after right shifting the mantissa, the resulting int32_t has correct sign.
  // With this choice of C, the number of mantissa bits we need to shift to the
  // right is: 52 - 33 = 19.
  //
  // Moreover, since the integer right shifts are equivalent to rounding down,
  // we can add an extra 0.5 so that it will become round-to-nearest, tie-to-
  // +infinity.  So in particular, we can compute:
  //   hmm = x * 2^12 * L2E + C,
  // where C = 2^33 + 2^32 + 2^-1, then if
  //   k = int32_t(lower 51 bits of double(x * 2^12 * L2E + C) >> 19),
  // the reduced argument:
  //   lo = x - log(2) * 2^-12 * k is bounded by:
  //   |lo| <= 2^-13 + 2^-41 + 2^-12*2^-19
  //         = 2^-13 + 2^-31 + 2^-41.
  //
  // Finally, notice that k only uses the mantissa of x * 2^12 * L2E, so the
  // exponent 2^12 is not needed.  So we can simply define
  //   C = 2^(33 - 12) + 2^(32 - 12) + 2^(-13 - 12), and
  //   k = int32_t(lower 51 bits of double(x * L2E + C) >> 19).

  // Rounding errors <= 2^-31 + 2^-41.
  double tmp = multiply_add(x, LOG2_E, 0x1.800000004p21);
  ix.f=tmp;
  int32_t k = ix.u >> 19;// as noted above if we shift as a 64 bit, then we end up with a 2's complement signed number!
  double kd = k;

  uint32_t idx1 = (k >> 6) & 0x3f;
  uint32_t idx2 = k & 0x3f;
  int hi = k >> 12;

  bool denorm = (hi <= -1022);
  dprintf("tmp=%.2f, k=%d (0x%08x) kd=%.2f idx1=%u idx2=%u hi=%d denorm=%s\n",tmp,k,k,kd,idx1,idx2,hi,denorm?"True":"False");
  DoubleDouble exp_mid1={EXP2_MID1[idx1].mid, EXP2_MID1[idx1].hi};
  DoubleDouble exp_mid2={EXP2_MID2[idx2].mid, EXP2_MID2[idx2].hi};

  DoubleDouble exp_mid = mult_dd_dd(exp_mid1, exp_mid2);

  // |x - (hi + mid1 + mid2) * log(2) - dx| < 2^11 * eps(M_LOG_2_EXP2_M12.lo)
  //                                        = 2^11 * 2^-13 * 2^-52
  //                                        = 2^-54.
  // |dx| < 2^-13 + 2^-30.
  double lo_h = multiply_add(kd, MLOG_2_EXP2_M12_HI, x); // exact
  double dx = multiply_add(kd, MLOG_2_EXP2_M12_MID, lo_h);

  // We use the degree-4 Taylor polynomial to approximate exp(lo):
  //   exp(lo) ~ 1 + lo + lo^2 / 2 + lo^3 / 6 + lo^4 / 24 = 1 + lo * P(lo)
  // So that the errors are bounded by:
  //   |P(lo) - expm1(lo)/lo| < |lo|^4 / 64 < 2^(-13 * 4) / 64 = 2^-58
  // Let P_ be an evaluation of P where all intermediate computations are in
  // double precision.  Using either Horner's or Estrin's schemes, the evaluated
  // errors can be bounded by:
  //      |P_(dx) - P(dx)| < 2^-51
  //   => |dx * P_(dx) - expm1(lo) | < 1.5 * 2^-64
  //   => 2^(mid1 + mid2) * |dx * P_(dx) - expm1(lo)| < 1.5 * 2^-63.
  // Since we approximate
  //   2^(mid1 + mid2) ~ exp_mid.hi + exp_mid.lo,
  // We use the expression:
  //    (exp_mid.hi + exp_mid.lo) * (1 + dx * P_(dx)) ~
  //  ~ exp_mid.hi + (exp_mid.hi * dx * P_(dx) + exp_mid.lo)
  // with errors bounded by 1.5 * 2^-63.

  double mid_lo = dx * exp_mid.hi;

  // Approximate expm1(dx)/dx ~ 1 + dx / 2 + dx^2 / 6 + dx^3 / 24.
  double p = poly_approx_d(dx);

  double lo = multiply_add(p, mid_lo, exp_mid.lo);

/*
  if (ya_UNLIKELY(denorm)) {
    if (auto r = ziv_test_denorm(hi, exp_mid.hi, lo, EXP_ERR_D);
        ya_LIKELY(r.has_value()))
      return r.value();
  } else {
*/
// Rounding tests for 2^hi * (mid + lo) when the output might be denormal. We
// assume further that 1 <= mid < 2, mid + lo < 2, and |lo| << mid.
// Notice that, if 0 < x < 2^-1022,
//   double(2^-1022 + x) - 2^-1022 = double(x).
// So if we scale x up by 2^1022, we can use
//   double(1.0 + 2^1022 * x) - 1.0 to test how x is rounded in denormal range.

//ziv_test_denorm(int hi, double mid, double lo, double err) {
  if(ya_UNLIKELY(denorm)) 
  	 {double mid=exp_mid.hi,err=EXP_ERR_D; // hi & lo are already correctly named
	
	  // Scaling factor = 1/(min normal number) = 2^1022
	  int64_t exp_hi =
	      (int64_t)(hi + 1022) * (1LL << DOUBLE_MANTISSA_BITS);	      
	  //double mid_hi = cpp::bit_cast<double>(exp_hi + cpp::bit_cast<int64_t>(mid));
	  ix.f=mid;
	  ix.u+=exp_hi;
	  double mid_hi=ix.f;
	  //double lo_scaled =
	  //    (lo != 0.0) ? cpp::bit_cast<double>(exp_hi + cpp::bit_cast<int64_t>(lo))
	  //                : 0.0;
	  ix.f=lo;
	  ix.u+=exp_hi;	  
	  double lo_scaled = (lo != 0.0) ? ix.f : 0.0;	                  
	
	  double extra_factor = 0.0;
	  uint64_t scale_down = 0x3FE0000000000000; // 1022 in the exponent field.
	
	  // Result is denormal if (mid_hi + lo_scale < 1.0).
	  if ((1.0 - mid_hi) > lo_scaled) 
	   {
	    // Extra rounding step is needed, which adds more rounding errors.
	    err += 0x1.0p-52;
	    extra_factor = 1.0;
	    scale_down = 0x3FF0000000000000; // 1023 in the exponent field.
	   }
	
	  // By adding 1.0, the results will have similar rounding points as denormal
	  // outputs.

	  //double err_scaled =
	  //    cpp::bit_cast<double>(exp_hi + cpp::bit_cast<int64_t>(err));
	  ix.f=err;
	  ix.u+=exp_hi;
	  double err_scaled=ix.f;
	  double lo_u = lo_scaled + err_scaled;
	  double lo_l = lo_scaled - err_scaled;
	
	  double upper = extra_factor + (mid_hi + lo_u);
	  double lower = extra_factor + (mid_hi + lo_l);
	
	  if (ya_LIKELY(upper == lower)) 
		 {
	      //return cpp::bit_cast<double>(cpp::bit_cast<uint64_t>(upper) - scale_down);
	      ix.f=upper;
	      ix.u-=scale_down;
	      return ix.f;
	     }
     } 
  else // if(ya_UNLIKELY(denorm))
  	{
	  double upper = exp_mid.hi + (lo + EXP_ERR_D);
	  double lower = exp_mid.hi + (lo - EXP_ERR_D);
	
	  if (ya_LIKELY(upper == lower))
	   {
	    // to multiply by 2^hi, a fast way is to simply add hi to the exponent
	    // field.
	    ix.f=upper;   
	    int64_t exp_hi = (int64_t)(hi) * (1LL << DOUBLE_MANTISSA_BITS);
	    ix.u+=exp_hi; // add hi to exponent
	    return ix.f;
	   }
	}
#ifdef USE_INTERNAL_DD_POLY	
  if(ya_LIKELY(!denorm))
  	{// do double-double case only for normalised numbers, as that's the common case
  	 DoubleDouble r_dd = exp_double_double(x, kd, exp_mid);
     double upper_dd = r_dd.hi + (r_dd.lo + EXP_ERR_DD);
     double lower_dd = r_dd.hi + (r_dd.lo - EXP_ERR_DD);
     if (ya_LIKELY(upper_dd == lower_dd)) 
	   {
	 	ix.f=upper_dd;   
	    int64_t exp_hi = (int64_t)(hi) * (1LL << DOUBLE_MANTISSA_BITS);
	    ix.u+=exp_hi; // add hi to exponent
	    return ix.f;
	   } 
	} 	 
#endif	
  double _cr_as_exp_accurate(double x); // in cr_exp.c 
  return _cr_as_exp_accurate(x); // accurate calculation 

}


