/* simple test program for ya_math_lib
   Written by Peter Miller, 1/9/2026
 
 The concept for this is https://members.loria.fr/PZimmermann/papers/accuracy.pdf which gives the accuracy of functions in a large number of maths libaries 
 
  While this works OK, and shows all max errors are 0.5ulp (as expected), note while the tests for all functions are almost identical, several tests have been "tweaked":
  	- sinq/cosq (rather than sinl/cosl) are used as references for sin/cos as sinl/cosl gave significant errors
    - range/increment in for loop in "searching for larger errors" varies depending on function being tested  (but all "special values" are always checked)
    	- in all cases over 66 Million values are checked for each function (for sin,cos,atan2 its over 69 Million)
    
 Most values in the table of special values are from https://gitlab.inria.fr/zimmerma/math_accuracy/-/blob/master/binary64/check_sample.c?ref_type=heads 
 
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <quadmath.h> /* see https://gcc.gnu.org/onlinedocs/libquadmath/quadmath_005fsnprintf.html#quadmath_005fsnprintf - also needs quadmath library linking in */
#include "../my_printf/my_printf.h"
#include "ya_crmath.h"


/* Peter Miller: code below cannot be compiled with -Ofast as this makes the compiler break some C rules that we need, so make sure of this here */
/* we also need -msse2 and -mfpmath=sse to actually use the sse instructions for float and double maths */
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)) || defined(__clang__)
 #if !defined(__BORLANDC__)
  #pragma GCC push_options
  #pragma GCC optimize ("-O3") /* cannot use Ofast, normally -O3 is OK. Note macro expansion does not work here ! */
 #endif
 // based on  https://jdebp.uk/FGA/predefined-macros-processor.html "__i386__" is set by GCC,Clang,Intel which is good enough as the outer #if limits us to gcc and clang
 #ifdef __i386__
   #pragma GCC target("sse2,fpmath=sse") /* -msse2 and -mfpmath=sse */
 #endif 
#endif


#define STRINGIZER(x) # x
#define TO_STRING(x) STRINGIZER(x)
#define nos_elements_in(x) (sizeof(x)/(sizeof(x[0]))) /* number of elements in x , max index is 1 less than this as we index 0... */

double special_vals[]=
	{NAN,-NAN,nan(""),-nan(""),nan("0"),-nan("0"),nan("1"),-nan("1"),INFINITY,-INFINITY,
	 -DBL_MAX, -DBL_MIN,-__DBL_DENORM_MIN__, -0.0,0.0,__DBL_DENORM_MIN__,DBL_MIN,DBL_MAX,
	 -1,-2,-4,-8,-16,-256,
	 1,2,4,8,16,256,
	 // values from here on from https://gitlab.inria.fr/zimmerma/math_accuracy/-/blob/master/binary64/check_sample.c?ref_type=heads on 3/9/2026 from its array worst
    /* acos */
    0x1.dffffb3488a4p-1,    /* glibc */
    0x1.6c05eb219ec46p-1,   /* icc */
    0x1.35b03e336a82bp-1,   /* AMD LibM */
    -0x1.0068b067c6feep-1,  /* Newlib */
    -0x1.0068b067c6feep-1,  /* OpenLibm */
    -0x1.0068b067c6feep-1,  /* Musl */
    -0x1.8d313198a2e03p-53, /* Apple 1.05154 */
    -0x1.cb3b399d747f3p-55, /* llvm */
    0x1.35a0de2b038fep-1,   /* libmvec 2.30611 */
    0x1.ffc00159839aep-1,   /* libmvec 2.05784 */
    0x1.35b9bac9f42c6p-1,   /* libmvec 1.82629 */
    -0x1.01f97a9d82a5dp-1,  /* msvc-x64-acos 0.931491 */
    -0x1.0068b067c6feep-1,  /* freebsd 0.929812 */
    0x1.251869c3f7881p-1,   /* ArmPL 1.51061 */

    /* acosh */
    0x1.c29124ea23a1cp+950, /* GNU libc */
    0x1.01825ca7da7e5p+0,   /* icx 0.508764 */
    0x1.209fae707a0edp+0,   /* AMD LibM */
    0x1.0001ff6afc4bap+0,   /* Newlib */
    0x1.0001ff6afc4bap+0,   /* OpenLibm */
    -0x1.34e729fd08086p+21, /* Musl Inf */
    0x1.00007fb3703ddp+0,   /* Apple 2.24182 */
    0,                      /* llvm */
    0x1.001ffe4f6d239p+0,   /* libmvec 3.31361 */
    0x1.007ff5e6aae25p+0,   /* libmvec 3.28513 */
    0x1.0007ffe4f42f8p+0,   /* libmvec 1.98412 */
    0x1.0007fd4307b75p+0,   /* msvc-x64-acosh 3.21636 */
    0x1.0001ff6afc4bap+0,   /* freebsd 2.24122 */
    0x1.071334daf83adp+0,   /* ArmPL 2.65022 */

    /* asin */
    -0x1.0000045b2c904p-3,  /* glibc */
    0x1.6c042a6378102p-1,   /* icc 19.1.3.304 */
    -0x1.00d44cccfa99p-1,   /* AMD LibM */
    -0x1.004d1c5a9400bp-1,  /* Newlib */
    -0x1.004d1c5a9400bp-1,  /* OpenLibm */
    -0x1.004d1c5a9400bp-1,  /* Musl 1.2.1 */
    0x1.eaeca101513d2p-2,   /* Apple 0.710076 */
    -0x1.7137449123ef6p-26, /* llvm */
    0x1.000c80481e7f9p-1,   /* libmvec 3.48486 */
    0x1.000fb59dbb7ffp-1,   /* libmvec 2.95635 */
    -0x1.0312655c1d169p-1,  /* libmvec 2.69006 */
    0x1.003d8c089658bp-1,   /* msvc-x64-asin 0.959391 */
    -0x1.004d1c5a9400bp-1,  /* freebsd 0.980823 */
    0x1.0479b37d95e5cp-1,   /* ArmPL 2.68343 */

    /* asinh */
    -0x1.7137449123ef7p-26, /* GNU libc */
    -0x1.000276d9cf31ap-4,  /* icx 0.506163 */
    0x1.005ae8d126f7ep+0,   /* AMD LibM */
    -0x1.02657ff36d5f3p-2,  /* Newlib */
    -0x1.02657ff36d5f3p-2,  /* OpenLibm */
    -0x1.0240f2bdb3f25p-2,  /* Musl 1.2.1 */
    -0x1.fdefd03df4cd7p-3,  /* Apple 1.57940 */
    0,                      /* llvm */
    -0x1.fff13ade9918p-12,  /* libmvec 3.58616 */
    0x1.fffdfee9d0656p-12,  /* libmvec 3.58591 */
    -0x1.fff14d29165f4p-8,  /* libmvec 1.52036 */
    -0x1.00212bb59c31ep-4,  /* msvc-x64-asinh 2.04895 */
    -0x1.02657ff36d5f3p-2,  /* freebsd 1.91356 */
    -0x1.000eeed78380ap+0,  /* ArmPL 2.03455 */

    /* atan */
    0x1.f9004c4fef9eap-4,   /* glibc */
    -0x1.ffff8020d3d1dp-7,  /* icc */
    -0x1.60370d15396b7p-1,  /* AMD LibM 0.862182 */
    0x1.62ff6a1682c25p-1,   /* Newlib */
    0x1.62ff6a1682c25p-1,   /* OpenLibm */
    0x1.62ff6a1682c25p-1,   /* Musl */
    0x1.01e0be37af68fp+1,   /* Apple 0.875241 */
    -0x1.085b3fdb46cbp+0,   /* llvm */
    0x1.000a5e848c0dp-3,    /* libmvec 2.64340 */
    0x1.0029e0e2db7dp-3,    /* libmvec 2.64176 */
    -0x1.0010aea41501p-3,   /* libmvec 2.64024 */
    -0x1.60370d15396b7p-1,  /* msvc-x64-atan 0.862182 */
    0x1.62ff6a1682c25p-1,   /* freebsd 0.860246 */
    0x1.032b4811f3dc5p+0,   /* ArmPL 2.23722 */

    /* atanh */
    0x1.a6a58d55e307bp-26,  /* GNU libc */
    -0x1.e2cfb2667f17ep-9,  /* icc */
    -0x1.d8fb311a52173p-2,  /* AMD LibM */
    -0x1.f97fabc0650c4p-4,  /* Newlib */
    -0x1.f97fabc0650c4p-4,  /* OpenLibm */
    -0x1.f8a404597baf4p-4,  /* Musl 1.2.1 */
    0x1.ffd834a270fp-10,    /* Apple 2.00021 */
    0,                      /* llvm */
    -0x1.fff0caf4c48dp-13,  /* libmvec 3.58751 */
    -0x1.ffe2abaa5690dp-13, /* libmvec 3.58538 */
    0x1.85cb7cc1e1318p-6,   /* libmvec 1.51047 */
    -0x1.ffbe8dd88527fp-9,  /* msvc-x64-atanh 2.49418 */
    -0x1.f97fabc0650c4p-4,  /* freebsd 1.80512 */
    -0x1.e7c1f36602014p-4,  /* ArmPL 2.9904 */

    /* cbrt */
    0x1.7a337e1ba1ec2p-257,  /* GNU libc */
    -0x1.f7af4893d1d51p-616, /* icc */
    -0x1.0c8125f4f1c39p+998, /* AMD LibM 0.502920 */
    -0x1.00ddafe7d9deep-885, /* Newlib */
    -0x1.13a5ccd87c9bbp+1008,/* OpenLibm */
    -0x1.13a5ccd87c9bbp+1008,/* Musl */
    0x1.fed1fe9f1122dp+11,   /* Apple 0.728127 */
    0x1.a202bfc89ddffp-1,    /* llvm 0.5 */
    -0x0.bd54cbc41f0b9p-1022,/* libmvec 3.42039 */
    0x0.bdf2e4b035cc5p-1022, /* libmvec 3.40348 */
    0x1.477fc84889eabp-511,  /* libmvec 1.83354 */
    -0x1.55cd285f321f6p-64,  /* msvc-x64-cbrt 1.85081 */
    -0x1.13a5ccd87c9bbp+1008,/* freebsd 0.66715 */
    0x1.fffad1cec59fep-332,  /* ArmPL 1.78544 */

    /* cos */
    -0x1.7120161c92674p+0,   /* GNU libc */
    -0x1.d19ebc5567dcdp+311, /* icc */
    0x1.9eb7148f354d6p+20,   /* AMD LibM */
    -0x1.4ae182c1ab422p+21,  /* Newlib */
    -0x1.34e729fd08086p+21,  /* OpenLibm */
    -0x1.34e729fd08086p+21,  /* Musl */
    0x1.2f29eb4e99fa2p+7,    /* Apple 0.947350 */
    0x1.16e534ee3658p-4,     /* llvm 0.5 */
    0x1.852e715836081p+18,   /* libmvec 3.84216 */
    -0x1.f5ec1ef4d1fb8p+3,   /* libmvec 3.66518 */
    -0x1.9a4f79002782p-6,    /* libmvec 3.65239 */
    -0x1.9200634d4471fp-1,   /* msvc-x64-cos 0.896141 */
    -0x1.34e729fd08086p+21,  /* freebsd 0.833135 */
    -0x1.7120161c92674p+0,   /* ArmPL 0.515628 */

    /* cosh */
    0x1p-26,   /* GNU libc */
    -0x1.5a364e6b98134p+9,   /* icc */
    0x1.ff76fb3f476d5p+0,    /* AMD LibM */
    0x1.633cc2ae1c934p+9,    /* Newlib */
    -0x1.6310ab92794a8p+9,   /* OpenLibm */
    -0x1.502bf5ad80729p+0,   /* Musl */
    -0x1.62dabd4848dc4p-2,   /* Apple 0.522287 */
    0,                       /* llvm */
    -0x1.633c654fee2bap+9,   /* libmvec 1.92222 */
    -0x1.633c654fee2bap+9,   /* libmvec 1.92222 */
    -0x1.2b3088f4a6e98p+4,   /* libmvec 2.02425 */
    0x1.ed991002dc978p+0,    /* msvc-x64-cosh 1.41338 */
    -0x1.6310ab92794a8p+9,   /* freebsd 1.46619 */
    -0x1.628af341989dap+9,   /* ArmPL 1.92568 */

    /* erf */
    -0x1.6e9b204d2167dp-539, /* GNU libc */
    -0x1.223fe3e7458a3p+0,   /* icx 0.541205 */
    0x1.11642f2eab9edp+0,    /* AMD LibM 0.999489 */
    -0x1.c57541b55c8ebp-16,  /* Newlib */
    -0x1.c57541b55c8ebp-16,  /* OpenLibm */
    -0x1.c57541b55c8ebp-16,  /* Musl 1.2.1 */
    -0x1.e057e7a0e494cp-2,   /* Apple 6.40351 */
    0,                       /* llvm */
    0x1.0000b7af4dcdp-8,     /* libmvec 2.54483 */
    0x1.00005abf94234p-8,    /* libmvec 2.54487 */
    0x1.00001d2920fb7p-8,    /* libmvec 2.54487 */
    -0x1.bf3c82e5f58dfp-8,   /* msvc-x64-erf 3.31306 */
    -0x1.c57541b55c8ebp-16,  /* freebsd 1.01967 */
    0x1.0000187085e56p-8,    /* ArmPL 2.2867 */

    /* erfc */
    -0x1.c5bf891b4ef6ap-54,  /* GNU libc */
    0x1.a8cf81fa2a28fp+4,    /* icx */
    0x1.3cf8ce875afeap+0,    /* AMD LibM 4.48819 */
    0x1.5182d8799b84bp+0,    /* Newlib */
    0x1.5182d8799b84bp+0,    /* OpenLibm */
    0x1.527f4fb0d9331p+0,    /* Musl 1.2.1 */
    0x1.bba14dc3507ccp+1,    /* Apple 10.6637 */
    0,                       /* llvm */
    /* GLIBC LIBMVEC has the medium accuracy Intel SVML functions
       (high accuracy: < 1 ulp, low accuracy: many ulps) */
    0x1.78afff6f3044cp+4,    /* libmvec 2.21867 */
    0x1.78affead86a26p+4,    /* libmvec 2.20423 */
    0x1.78afff9d452cp+4,     /* libmvec 2.20537 */
    0x1.f5c5a6822a341p+1,    /* msvc-x64-erfc 6.03199 */
    0x1.5182d8799b84bp+0,    /* freebsd 4.07682 */
    0x1.46cffdf330b13p+4,    /* ArmPL 1.70623 */

    /* exp */
    -0x1.49f33ad2c1c58p+9,   /* GNU libc 0.510359 */
    -0x1.9df9aecc00001p-29,  /* GNU libc RNDZ 1.24399 */
    0x1.50f5a68c91c18p-1,    /* GNU libc RNDU 1.01597 */
    -0x1.9df9aecc00001p-29,  /* GNU libc RNDD 1.24399 */
    0x1.fce66609f7428p+5,    /* icc 19.1.3.304 */
    0x1.b97dc8345c55p+5,     /* AMD LibM */
    0x1.2e8f20cf3cbe7p+8,    /* Newlib */
    0x1.2e8f20cf3cbe7p+8,    /* OpenLibm */
    -0x1.18209ecd19a8cp+6,   /* Musl */
    -0x1.4133f4fd79c1cp-13,  /* Apple 0.520417 */
    0x1p-53,                 /* llvm 0.5 */
    -0x1.205968aae119fp-8,   /* libmvec 3.20379 */
    -0x1.2059763f8882bp-8,   /* libmvec 3.20362 */
    -0x1.205968a73d4abp-8,   /* libmvec 3.20361 */
    -0x1.ffffffa000001p-27,  /* msvc-x64-exp 1.5 */
    0x1.2e8f20cf3cbe7p+8,    /* freebsd 0.948331 */
    -0x1.49f33ad2c1c58p+9,   /* ArmPL 0.510359 */

    /* exp10 */
    -0x1.585fadeab79d4p-7,   /* GNU libc */
    -0x1.5cd9d94d49a85p+1,   /* icc */
    -0x1.285d82b75258fp+2,   /* AMD LibM */
    0x1.ce7ef793d4b0ap-2,    /* Newlib */
    0,                       /* OpenLibm */
    -0x1.fe8c27141c94ap+3,   /* Musl */
    -0x1.c37443e446523p-16,  /* Apple 0.520238 */
    0x1.bcb7b1526e50dp-55,   /* llvm 0.5 */
    0x1.33f4082f47b74p+8,    /* libmvec 2.00015 */
    0x1.33f4082f47b74p+8,    /* libmvec 2.00015 */
    0x1.33f4082f47b74p+8,    /* libmvec 2.00015 */
    0,                       /* msvc_x64-exp10 0.0 */
    0,                        /* freebsd: NA */
    -0x1.5acf96f42165bp-7,   /* ArmPL 0.509874 */

    /* exp2 */
    -0x1.1a4ce073ea908p-5,   /* GNU libc */
    -0x1.8002f29666b99p-6,   /* icx */
    0x1.fffbff4152bafp+9,    /* AMD LibM 1.04695 */
    -0x1.ff95ecb4e6331p-2,   /* Newlib */
    -0x1.ff1eb5acee46bp+9,   /* OpenLibm */
    -0x1.1a4ce073ea908p-5,   /* Musl */
    -0x1.b3d9b47ad1b2fp-13,  /* Apple 0.520087 */
    0x1.71547652b82fep-53,   /* llvm */
    -0x1.4c31866f6d3bbp-6,   /* libmvec 1.64293 */
    -0x1.4c3c931a5de98p-6,   /* libmvec 1.64097 */
    -0x1.8000e569a5545p-3,   /* libmvec 1.05024 */
    -0x1.e383ca1fd4462p+9,   /* msvc-x64-exp2 0.505 */
    -0x1.ff1eb5acee46bp+9,   /* freebsd 0.750695 */
    -0x1.f7087fb1cf9e8p+9,   /* ArmPL 0.508998 */

    /* expm1 */
    0x1.63e063d87938p-2,     /* GNU libc 0.912463 */
    -0x1.62fe464c64f65p-8,   /* icc */
    0x1.9a5796c0b73b8p-2,    /* AMD LibM 0.669268 */
    0x1.62ff47a01658fp-2,    /* Newlib */
    0x1.62ff47a01658fp-2,    /* OpenLibm */
    0x1.62ff47a01658fp-2,    /* Musl */
    0x1.e7f93188565ecp-5,    /* Apple 0.705059 */
    0x1p-52,                 /* llvm 0.5 */
    0x1.86f57e8de4a5p-9,     /* libmvec 2.96513 */
    0x1.856b41d86994cp-9,    /* libmvec 2.96463 */
    0x1.c3b7c858f0575p-6,    /* libmvec 2.11697 */
    -0x1.62d7c116d2e32p-1,   /* msvc-x64-expm1 3.05702 */
    0x1.62ff47a01658fp-2,    /* freebsd 0.908146 */
    0x1.633f993a730c9p-2,    /* ArmPL 2.17324 */

    /* j0 */
    0x1.33d152e971b4p+1,     /* GNU libc */
    0x1.aff859518c846p+7,    /* icc */
    0,                       /* AMD LibM */
    0x1.45f3067a0f4b2p+847,  /* Newlib */
    0x1.33d152e971b4p+1,     /* OpenLibm */
    -0x1.33d152e971b4p+1,    /* Musl 1.2.1 */
    0x1.6148f5b2c2e45p+2,    /* Apple 3.82247e+14 */
    0,                       /* llvm */
    0,                       /* libmvec */
    0x1.782b7a20df6d4p+66,   /* msvc-x64-_j0 2.41665e+33 */
    0x1.33d152e971b4p+1,     /* freebsd 4.51e+14 */
    0x1.6148f5b2c2e45p+2,    /* ArmPL 3.82247e+14 */

    /* j1 */
    -0x1.ea75575af6f09p+1,   /* GNU libc */
    -0x1.67b5541c7d8b7p+7,   /* icc 19.1.3.304 */
    0,                       /* AMD LibM */
    0x1.45f3066f80258p+325,  /* Newlib */
    -0x1.ea75575af6f09p+1,   /* OpenLibm */
    0x1.ea75575af6f09p+1,    /* Musl 1.2.1 */
    -0x1.ea75575af6f09p+1,   /* Apple 1.10e+15 */
    0,                       /* llvm */
    0,                       /* libmvec */
    0x1.99caa5236feeap+75,   /* msvc-x64-_j1 8.18612e+32 */
    -0x1.ea75575af6f09p+1,   /* freebsd 1.10e+15 */
    -0x1.ea75575af6f09p+1,   /* ArmPL 8.95204e+14 */

    /* lgamma */
    0x0.00000002f6806p-1022, /* GNU libc */
    -0x1.3f62c60e23b31p+2,   /* icc 19.1.3.304 */
    0,                       /* AMD LibM */
    -0x1.3a7fc9600f86cp+1,   /* RedHat Newlib 4.0.0 */
    -0x1.3a7fc9600f86cp+1,   /* OpenLibm */
    -0x1.3a7fc9600f86cp+1,   /* Musl 1.2.1 */
    -0x1.bffcbf76b86fp+2,    /* Apple 2.32471e+16 */
    0,                       /* llvm */
    0,                       /* libmvec */
    -0x1.e383ca1fd4462p+9,   /* msvc-x64-lgamma Inf */
    -0x1.3a7fc9600f86cp+1,   /* freebsd 4.45e+15 */
    -0x1.f60bdba970b3cp+1,   /* ArmPL 9.52506 */

    /* log */
    0x1.1211bef8f68e9p+0,    /* GNU libc */
    0x1.008000db2e8bep+0,    /* icc */
    0x1.108058fed54b3p+0,    /* AMD LibM 0.610401 */
    0x1.48ae5a67204f5p+0,    /* Newlib */
    0x1.48ae5a67204f5p+0,    /* OpenLibm */
    0x1.dc0b586f2b26p-1,     /* Musl */
    0x1.490af72a25a81p-1,    /* Apple 0.507847 */
    0x1.5b6e7e4e96f86p+2,    /* llvm */
    0x1.00e000c7fa1c3p+0,    /* libmvec 1.58569 */
    0x1.002001ffaa4ap+0,     /* libmvec 1.58883 */
    0x1.001f01ac83b3p+0,     /* libmvec 1.59111 */
    0x1.0ffc349469a2fp+0,    /* msvc-x64-log 0.576920 */
    0x1.48ae5a67204f5p+0,    /* freebsd 0.945306 */
    0x1.1211bef8f68e9p+0,    /* ArmPL 0.519023 */

    /* log10 */
    0x1.de02157073b31p-1,    /* GNU libc */
    0x1.feda7b62c1033p-1,    /* icc 19.1.3.304 */
    0x1.10fdf4211fd45p+0,    /* AMD LibM */
    0x1.55535a0140a21p+0,    /* Newlib */
    0x1.553e1cb579ee9p+0,    /* OpenLibm */
    0x1.553e1cb579ee9p+0,    /* Musl */
    0x1.2501ee5628b08p-1,    /* Apple 0.513076 */
    0x1.e12d66744ff81p+429,  /* llvm 0.5 */
    0x1.00201204555c8p+0,    /* libmvec 2.09444 */
    0x1.001fffbd3f495p+0,    /* libmvec 2.09480 */
    0x1.f03ebdaea826bp-1,    /* libmvec 1.95902 */
    0x1.e005e1e891807p-1,    /* msvc-x64-log10 0.632134 */
    0x1.553e1cb579ee9p+0,    /* freebsd 0.813227 */
    0x1.de02157073b31p-1,    /* ArmPL 1.61828 */

    /* log1p */
    -0x1.2bf183e0344b2p-2,   /* GNU libc */
    0x1.000aee2a2757fp-9,    /* icc */
    0x1.e0013fd35cbbp-4,     /* AMD LibM */
    -0x1.2bf1de6b04a8ap-2,   /* Newlib */
    -0x1.2bf1de6b04a8ap-2,   /* OpenLibm */
    -0x1.2bf32aaf122e2p-2,   /* Musl */
    -0x1.ffffff3fffffdp-28,  /* Apple 0.666667 */
    0x1p-53,                 /* llvm */
    0x1.000bcdec306p-11,     /* libmvec 2.58927 */
    0x1.fff86f9b9acp-12,     /* libmvec 2.58923 */
    0x1.075745181aabp-6,     /* libmvec 1.94684 */
    -0x1.8000000000000p-53,  /* msvc-x64-log1p 1.4375 */
    -0x1.2bf1de6b04a8ap-2,   /* freebsd 0.89558 */
    -0x1.2e496d25897ecp-2,   /* ArmPL 1.73667 */

    /* log2 */
    0x1.1406d79e1b574p+0,    /* GNU libc */
    0x1.fea11b9451f9bp-1,    /* icx 0.509 */
    0x1.0b541b6746bd1p+0,    /* AMD LibM */
    0x1.68d778f076021p+0,    /* Newlib */
    0x1.67eaf07ce24d1p+0,    /* OpenLibm */
    0x1.0b53197bd66c8p+0,    /* Musl */
    0x1.6d0092bcf066dp-1,    /* Apple 0.514489 */
    0x1.b4ebe40c95a01p+0,    /* llvm */
    0x1.002000d8e91c5p+0,    /* libmvec 2.08932 */
    0x1.002003e5a80e3p+0,    /* libmvec 2.08921 */
    0x1.ede4ac763282bp-1,    /* libmvec 1.86313 */
    0x1.160732376ad7fp+0,    /* msvc-x64-log2 0.811615 */
    0x1.67eaf07ce24d1p+0,    /* freebsd 0.920898 */
    0x1.0b53f741fb8c4p+0,    /* ArmPL 0.553243 */

    /* sin */
    -0x1.f8b791cafcdefp+4,   /* GNU libc */
    -0x1.0e16eb809a35dp+944, /* icc */
    0x1.9eb7148f354d6p+21,   /* AMD LibM */
    -0x1.842d8ec8f752fp+21,  /* Newlib */
    0x1.4d84db080b9fdp+21,   /* OpenLibm */
    0x1.4d84db080b9fdp+21,   /* Musl */
    -0x1.07e4c92b5349dp+4,   /* Apple 0.943598 */
    0x1.1031fb019618p-8,     /* llvm 0.5 */
    0x1.a5a68e24971a3p+20,   /* libmvec 3.83437 */
    0x1.9977bea4253f1p+0,    /* libmvec 3.48842 */
    -0x1.99631ed67b43fp+0,   /* libmvec 3.48873 */
    -0x1.11b624b546894p+9,   /* msvc-x64-sin 0.798290 */
    0x1.4d84db080b9fdp+21,   /* freebsd 0.83076 */
    -0x1.f8b791cafcdefp+4,   /* ArmPL 0.515623 */

    /* sinh */
    0x1.2b705c15abeb9p+5,    /* GNU libc 1.93369e+19 */
    -0x1.adc135eb544c1p-2,   /* icc */
    0x1.1feb2a79f307p+3,     /* AMD LibM */
    -0x1.633cae1335f26p+9,   /* Newlib */
    -0x1.63324af2fb5b7p-1,   /* OpenLibm */
    -0x1.63324af2fb5b7p-1,   /* Musl */
    0x1.d7131e11fc6b3p-2,    /* Apple 0.538334 */
    0,                       /* llvm */
    -0x1.c5c9440e9422dp-9,   /* libmvec 2.39492 */
    -0x1.633c654fee2bap+9,   /* libmvec 1.92222 */
    -0x1.633c654fee2bap+9,   /* libmvec 1.92222 */
    0x1.1fd8c752e39f5p+3,    /* msvc-x64-sinh 1.40707 */
    -0x1.63324af2fb5b7p-1,   /* freebsd 1.87405 */
    0x1.9fa32b1149d35p-2,    /* ArmPL 2.58013 */

    /* sqrt */
    0x1.fffffffffffffp-1,    /* GNU libc */
    0x1.fffffffffffffp-1,    /* icc 19.1.3.304 */
    0x1.fffffffffffffp-1,    /* AMD LibM */
    0x1.fffffffffffffp-1,    /* Newlib */
    0x1.fffffffffffffp-1,    /* OpenLibm */
    0x1.fffffffffffffp-1,    /* Musl 1.2.1 */
    0x1.fffffffffffffp-1,    /* Apple */
    0x1.fffffffffffffp-1,    /* llvm */
    0x1.fffffffffffffp-1,    /* libmvec */
    0x1.fffffffffffffp-1,    /* msvc-x64-sqrt */
    0x1.fffffffffffffp-1,    /* freebsd */
    0x1.fffffffffffffp-1,    /* ArmPL 0.5 */

    /* tan */
    -0x1.317cd745dd37cp+9,   /* glibc */
    0x1.49adfd996a81dp+18,   /* icc */
    0x1.93c05c9ed3cbcp+18,   /* AMD LibM */
    0x1.3f9605aaeb51bp+21,   /* Newlib */
    0x1.3f9605aaeb51bp+21,   /* OpenLibm */
    0x1.3f9605aaeb51bp+21,   /* Musl 1.2.2 */
    -0x1.a81d98fc58537p+6,   /* Apple 3.52535 */
    0x1.d12ed0af1a27fp-27,   /* llvm 0.5 */
    0x1.72e90b4651593p+15,   /* libmvec 3.96274 */
    0x1.3fab696843fbfp+8,    /* libmvec 3.53263 */
    -0x1.780c9aeca907cp+17,  /* libmvec 3.98992 */
    -0x1.4d7c8b8320237p+11,  /* msvc-x64-tan 1.31973 */
    0x1.3f9605aaeb51bp+21,   /* freebsd 1.01166 */
    -0x1.317cd745dd37cp+9,   /* ArmPL 0.618367 */

    /* tanh */
    0x1.d12ed0af1a27fp-27,   /* GNU libc CR */
    0x1.002629fd74484p+0,    /* icc */
    -0x1.fde5bd2769a01p-1,   /* AMD LibM */
    -0x1.e134557098e37p-3,   /* Newlib */
    -0x1.e134557098e37p-3,   /* OpenLibm */
    -0x1.e134557098e37p-3,   /* Musl */
    0x1.00cf9f273d84p+1,     /* Apple 0.612229 */
    0,                       /* llvm */
    -0x1.000a02c5a8b47p-2,   /* libmvec 2.17016 */
    -0x1.00010c3967f16p-2,   /* libmvec 2.13884 */
    -0x1.001bf41f56582p-1,   /* libmvec 1.19944 */
    0x1.fefc7206bc932p-1,    /* msvc-x64-tanh 1.32827 */
    -0x1.e134557098e37p-3,   /* freebsd 2.21499 */
    -0x1.c416448380debp-3,   /* ArmPL 2.75884 */

    /* tgamma */
    0x1.8b9710d0adf83p-702,  /* GNU libc */
    -0x1.3e0001ad3bee3p+6,   /* icc */
    0,                       /* AMD LibM */
    -0x1.535175475cc8dp+7,   /* Newlib */
    -0x1.540b170c4e65ep+7,   /* OpenLibm */
    -0x1.fc4b534c8eccp+2,    /* Musl */
    -0x1.5456e56919a19p+7,   /* Apple 1026.83 */
    0,                       /* llvm */
    0,                       /* libmvec */
    -0x1.5c00000003c15p+7,   /* msvc-x64-tgamma 9.00718e+15 */
    -0x1.547cf3ddaaddap+7,   /* freebsd 1026.47 */
    -0x1.c18caecc00f7bp+2,   /* ArmPL 8.67183 */

    /* y0 */
    0x1.c982eb8d417eap-1,    /* GNU libc */
    0x1.4cdee58a47eddp-31,   /* icc */
    0,                       /* AMD LibM */
    0x1.c982eb8d417eap-1,    /* Newlib */
    0x1.c982eb8d417eap-1,    /* OpenLibm */
    0x1.c982eb8d417eap-1,    /* Musl 1.2.1 */
    0x1.c982eb8d417eap-1,    /* Apple 1.42e+15 */
    0,                       /* llvm */
    0,                       /* libmvec */
    0x1.99caa5236feeap+75,   /* msvc-x64-_y0 8.18612e+32 */
    0x1.c982eb8d417eap-1,    /* freebsd 1.42e+15 */
    0x1.c982eb8d417eap-1,    /* ArmPL 3.48447e+15 */

    /* y1 */
    0x1.193bed4dff243p+1,    /* GNU libc */
    0x1.c513c569fe78ep+0,    /* icc */
    0,                       /* AMD LibM */
    0x1.193bed4dff243p+1,    /* Newlib */
    0x1.193bed4dff243p+1,    /* OpenLibm */
    0x1.193bed4dff243p+1,    /* Musl 1.2.1 */
    0x1.193bed4dff243p+1,    /* Apple 5.56e+15 */
    0,                       /* llvm */
    0,                       /* libmvec */
    0x1.782b7a20df6d4p+66,   /* msvc-x64-_y1 2.41665e+33 */
    0x1.193bed4dff243p+1,    /* freebsd 5.56e+15 */
    0x1.193bed4dff243p+1,    /* ArmPL 4.70515e+15 */

    /* acospi */
    0x1.da0ee5e5d4298p-1,    /* GNU libc 1.50743 */
    0x1.6a1efc80ddbb1p-1,    /* icc 0.540236 */
    0,                       /* AMD LibM */
    0,                       /* Newlib */
    0,                       /* OpenLibm */
    0,                       /* Musl */
    0,                       /* Apple */
    0,                       /* llvm */
    0,                       /* libmvec */
    0,                       /* msvc-x64-acospi */
    0,                       /* freebsd */
    0,                       /* ArmPL */

    /* asinpi */
    -0x1.8805060cb885cp-3,   /* GNU libc 1.49879 */
    0x1.921fabd2a2b25p-750,  /* icc 0.520834 */
    0,                       /* AMD LibM */
    0,                       /* Newlib */
    0,                       /* OpenLibm */
    0,                       /* Musl */
    0,                       /* Apple */
    0,                       /* llvm */
    0,                       /* libmvec */
    0,                       /* msvc-x64-asinpi */
    0,                       /* freebsd */
    0,                       /* ArmPL */

    /* atanpi */
    0x1.9601a785ff35dp-3,    /* GNU libc 1.50040 */
    0x1.a7fb417ac774ap-2,    /* icx 0.941040 */
    0,                       /* AMD LibM */
    0,                       /* Newlib */
    0,                       /* OpenLibm */
    0,                       /* Musl */
    0,                       /* Apple */
    0,                       /* llvm */
    0,                       /* libmvec */
    0,                       /* msvc-x64-atanpi */
    0,                       /* freebsd */
    0,                       /* ArmPL */

    /* cospi */
    -0x1.1a0a2fa299b92p+6,   /* GNU libc 1.85020 */
    -0x1.f01d619f61e5bp-8,   /* icx 0.502370 */
    -0x1.fffffffffffffp-28,  /* AMD LibM 2.46741 */
    0,                       /* Newlib */
    0,                       /* OpenLibm */
    0,                       /* Musl */
    -0x1.554bb66ede184p-1,   /* Apple */
    0,                       /* llvm */
    0,                       /* libmvec */
    0,                       /* msvc-x64-cospi */
    -0x1.fe3bb5207682dp-3,   /* freebsd 0.811621 */
    0x1.5a33cc258p-22,       /* ArmPL 3.17244 */

    /* exp10m1 */
    0x1.e880c5bafbd41p-2,    /* GNU libc 3.53871 */
    0,                       /* icx */
    0,                       /* AMD LibM */
    0,                       /* Newlib */
    0,                       /* OpenLibm */
    0,                       /* Musl */
    0,                       /* Apple */
    0,                       /* llvm */
    0,                       /* libmvec */
    0,                       /* msvc-x64-exp10m1 */
    0,                       /* freebsd */
    0,                       /* ArmPL */

    /* exp2m1 */
    0x1.d047583a6c6dp-1,     /* GNU libc 1.92510 */
    0,                       /* icx */
    0,                       /* AMD LibM */
    0,                       /* Newlib */
    0,                       /* OpenLibm */
    0,                       /* Musl */
    0,                       /* Apple */
    0,                       /* llvm */
    0,                       /* libmvec */
    0,                       /* msvc-x64-exp2m1 */
    0,                       /* freebsd */
    0,                       /* ArmPL */

    /* log10p1 */
    -0x1.4c2971893052fp-1,   /* GNU libc 1.94190 */
    0,                       /* icx */
    0,                       /* AMD LibM */
    0,                       /* Newlib */
    0,                       /* OpenLibm */
    0,                       /* Musl */
    0,                       /* Apple */
    0,                       /* llvm */
    0,                       /* libmvec */
    0,                       /* msvc-x64-log10p1 */
    0,                       /* freebsd */
    0,                       /* ArmPL */

    /* log2p1 */
    0x1.a7b725780ff2cp-2,    /* GNU libc 1.87575 */
    0,                       /* icx */
    0,                       /* AMD LibM */
    0,                       /* Newlib */
    0,                       /* OpenLibm */
    0,                       /* Musl */
    0,                       /* Apple */
    0,                       /* llvm */
    0,                       /* libmvec */
    0,                       /* msvc-x64-log2p1 */
    0,                       /* freebsd */
    0,                       /* ArmPL */

    /* sinpi */
    -0x1.45f3e53e1d707p-7,   /* GNU libc 1.85066 */
    -0x0.07e17f27e5086p-1022,/* icx 0.531251 */
    -0x1.ffffed20ca4dep-14,  /* AMD LibM 2.55084 */
    0,                       /* Newlib */
    0,                       /* OpenLibm */
    0,                       /* Musl */
    0x1.54c657059ae8dp-3,    /* Apple */
    0,                       /* llvm */
    0,                       /* libmvec */
    0,                       /* msvc-x64-sinpi */
    0x1.400ad6f5900efp+0,    /* freebsd 0.811311 */
    -0x1.00000d362cd4dp-1,   /* ArmPL 3.16284 */

    /* tanpi */
    -0x1.fae7d0ef22d4ep-2,   /* GNU libc 2.84904 */
    0x1.49b79692667bp+46,    /* icc 1.00 */
    0x1.ffffce32e4a7cp-15,   /* AMD LibM 2.55083 */
    0,                       /* Newlib */
    0,                       /* OpenLibm */
    0,                       /* Musl */
    -0x1.fe9eb8d802c61p-3,   /* Apple */
    0,                       /* llvm */
    0,                       /* libmvec */
    0,                       /* msvc-x64-tanpi */
    0x1.c9542a6e8f18fp+0,    /* freebsd 0.968891 */
    0x1.68bf660eea82p-2,     /* ArmPL 2.22117 */

    /* rsqrt */
    0x1.000080005p-620,       /* GNU libc 1.499990463176802 */
    0x1.0000000000002p-926, /* GNU libc RNDZ 2.0 */
	0x1.0000007ffffffp-62, /* GNU libc RNDU 2.0 */
    0x1.0000000000002p-748,  /* GNU libc RNDD 2.0 */
    0x1.00018f5913816p-458,  /* icc 0.500733 */
    0,                       /* AMD LibM */
    0,                       /* Newlib */
    0,                       /* OpenLibm */
    0,                       /* Musl */
    0,                       /* Apple */
    0,                       /* llvm */
    0,                       /* libmvec */
    0,                       /* msvc-x64-rsqrt */
    0,                       /* freebsd */
    0,                       /* ArmPL */
 // following vallues from check_sample.c array extra
  /* don't remove the following values: they should give an error > 0.5 for
     glibc asin after commit f67f9c9 */
  0x1.fcd5742999ab8p-1,
  -0x1.ee2b43286db75p-1,
  -0x1.f692ba202abcp-4,
  -0x1.9915e876fc062p-1,
  -0x1.fd7d13f1663afp-1, /* 5.0000005122065272e-01 for asin after f67f9c9 */
  0x1.16c08b622e36p-1, /* 4.9999999999919309e-01 for asin with glibc-2.32 */
  /* same for glibc acos */
  0x1.f63845056f35ep-1,
  0x1.fffff3634acd6p-1, /* glibc-2.32 acos 0.5000000044534775 */
  -0x1.cb3b399d747f3p-55,  /* crlibm bug acos_rn */
  0x1.1a62633145c07p-54,   /* crlibm bug acos_ru */
  0x1.800010834fed5p-1,    /* crlibm bug acospi_rn */
  -0x1.80052ddfdd07cp-1,   /* crlibm bug asin_rn */
  0x1p+0,                  /* crlibm bug atanpi_ru */
  0x1.499c0c2050962p-997,  /* crlibm bug sinpi_rn */
  -0x1.010e23e83c2a7p-997, /* crlibm bug sinpi_rz, sinpi_ru and tanpi_ru */
  -0x1.08b9e4381190dp-999, /* crlibm bug sinpi_rd */
  0x1.dd113d1bb494bp-998,  /* crlibm bug tanpi_rn */
  -0x1.ffffab1933322p-32,  /* crlibm bug tanpi_rz */
  0x1.fffbbb0468123p-32,   /* crlibm bug tanpi_rd */
  +0.0,                    /* msvc _y1 yields NaN */
  -0.0,
  0x1p+0, /* 10^0 */
  0x1.4p+3, /* 10^1 */
  0x1.9p+6, /* 10^2 */
  0x1.f4p+9, /* 10^3 */
  0x1.388p+13, /* 10^4 */
  0x1.86ap+16, /* 10^5 */
  0x1.e848p+19, /* 10^6 */
  0x1.312dp+23, /* 10^7 */
  0x1.7d784p+26, /* 10^8 */
  0x1.dcd65p+29, /* 10^9 */
  0x1.2a05f2p+33, /* 10^10 */
  0x1.74876e8p+36, /* 10^11 */
  0x1.d1a94a2p+39, /* 10^12 */
  0x1.2309ce54p+43, /* 10^13 */
  0x1.6bcc41e9p+46, /* 10^14 */
  0x1.c6bf52634p+49, /* 10^15 */
  0x1.1c37937e08p+53, /* 10^16 */
  0x1.6345785d8ap+56, /* 10^17 */
  0x1.bc16d674ec8p+59, /* 10^18 */
  0x1.158e460913dp+63, /* 10^19 */
  0x1.5af1d78b58c4p+66, /* 10^20 */
  0x1.b1ae4d6e2ef5p+69, /* 10^21 */
  0x1.0f0cf064dd592p+73, /* 10^22 */
  /* the following are inputs near zeros of Bessel functions j0,j1,y0,y1
     from Table 3 of https://www.cl.cam.ac.uk/~jrh13/papers/bessel.pdf */
  0x1.782b7a20df6d4p+66, /* j0 and y1 */
  0x1.66bd5424e5655p+89, /* j0 and y1 */
  0x1.07f285bb70bfp+38,  /* j0 */
  0x1.99caa5236feeap+75, /* j1 and y0 */
  0x1.e0664dbedfec5p+80, /* j0 and y1 */
  0x1.dbabcc5913d2fp+26, /* y0 */
  0x1.1a209b98a791fp+68, /* j1 and y0 */
  0x1.f19e5d71b26bap+84, /* j0 and y1 */
  0x1.b8fd7a4646112p+34, /* j1 and y0 */
  0x1.51eae0a757998p+10, /* j1 */
  0x1.4997d3ab50ac4p+30, /* j0 and y1 */
  0x1.097db7f70d377p+22, /* j1 and y0 */
  0x1.d63658a917489p+68, /* j0 and y1 */
  0x1.56a4aa740a5a7p+52, /* j0 and y1 */
  0x1.2ace52fd26fc2p+42, /* j0 and y1 */
  0x1.737d16f5eb092p+36, /* j1 and y0 */
  0x1.6f25c6bc40e4bp+17, /* j0 and y1 */
  0x1.2a32976bff5p+26,   /* j1 and y0 */
  0x1.ad01ea3bcbd83p+32, /* j0 and y1 */
  0x1.dfd42a8864ad7p+26, /* j1 and y0 */
  0x1.908de75d3884fp+7,  /* j0 */
  0x1.fc015e6663696p+12, /* j1 and y0 */
  0x1.cd4732748519dp+23, /* j1 and y0 */
  0x1.c982eb8d417eap-1,  /* y0 */
  0x1.0192a33615d8cp+0,  /* acosh icx */
  0x1.0413747aa0c3fp+0, /* llvm 21.1.8 atan rndu 1.9842 */
  // the following values are from https://sourceware.org/pipermail/libc-alpha/2026-July/178756.html
    103993.0,
    104348.0,
    208341.0,
    260515.0,
    312689.0,
    573204.0,
    833719.0,
    1146408.0,
    1980127.0,
    2292816.0,
    3126535.0,
    4272943.0,
    4846147.0,
    5419351.0,
    6565759.0,
    9692294.0,
    10265498.0,
    10838702.0,
    37362253.0,
    42781604.0,
    74724506.0,
    80143857.0,
    122925461.0,
    165707065.0,
    171126416.0,
    245850922.0,
    411557987.0,
    534483448.0,
    657408909.0,
    1068966896.0,
    1480524883.0,
    2137933792.0,
    2549491779.0,
    3083975227.0,
    6167950454.0,
    8717442233.0,
    12335900908.0,
    14885392687.0,
    17969367914.0,
    21053343141.0,
    35938735828.0,
    42106686282.0,
    56992078969.0,
    63160029423.0,
    78045422110.0,
    84213372564.0,
    99098765251.0,
    105266715705.0,
    120152108392.0,
    126320058846.0,
    231586774551.0,
    442120205961.0,
    881156436695.0,
    902209779836.0,
    1720206187108.0,
    1783366216531.0,
    2685575996367.0,
    3587785776203.0,
    5371151992734.0,
    8958937768937.0,
    17917875537874.0,
    65398140378926.0,
    74357078147863.0,
    139755218526789.0,
    145126370519523.0,
    148714156295726.0,
    214112296674652.0,
    279510437053578.0,
    288469374822515.0,
    428224593349304.0,
    567979811876093.0,
    856449186698608.0,
    1284673780047912.0,
    3137327371971917.0,
    5706674932067741.0,
    5920787228742393.0,
    6027843377079719.0,
    6134899525417045.0,
    12269799050834090.0
 
	 
	};

int main(int argc, char *argv[]) 
{double x,y_cr,me_x,m_rel_e_x,m_rel_e_y;// m_rel_e_y as a double is correct, as we use nextafter() on this to get the value for 1upl (as a double)
 long double y,e,me,me_y,rel_e,m_rel_e,m_rel_a; // compare against long double function
 double yr;
 int errs=0;
 
 // check sqrt
#define ref_sqrt_LD(x) __builtin_sqrtl(x) /* function thats taken to be "exact" - either sqrtl(x) or __builtin_sqrtl(x) - must take a long double argument (x) and return a long double */
#define ref_sqrt_DBL(x) __builtin_sqrt(x) /* function thats taken to be "exact" - either sqrt(x) or __builtin_sqrt(x) - must take a ldouble argument (x) and return a double - only used for NAM,INF */
#define sqrt_under_test(x) cr_sqrt(x) /* function under test - normally here that is cr_sqrt(x), but can also test library sqrt - must take a double argument (x) and return a double */ 
 
 // first check is just used to initialise variables to sensible values, actual value of x used is not critical, but it must be a valid argument for the function under test  
 x=0x1.fffffffffffffp-1;
 printf("Sizeof double=%u size of long double=%u\n",(unsigned)sizeof(double),(unsigned)sizeof(long double));
 printf("Complied with C standard version=%u on GCC %s\n\n",(unsigned int)__STDC_VERSION__,__VERSION__);
 printf("Testing %s vs %s : \n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)));
 y=ref_sqrt_LD((long double)x);
 y_cr=sqrt_under_test(x);
 e=fabsl((long double)y_cr-y); // error
 me=e;// use this value as "baseline" for wider test below
 me_x=x;
 me_y=y;
 rel_e=fabsl(e/y);// rel error [ relative to correct result]
 m_rel_e=rel_e;
 m_rel_e_x=x;
 m_rel_e_y=y;
 m_rel_a=e;	
 rel_e=fabsl(me/y);  
 printf("Checking values in special_vals[]...\n");
 for(int i=0;i<nos_elements_in(special_vals);++i)
	{
	 x=special_vals[i];
	 yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;	  	
	 	  	}
	 	}	 	
	
	}
 printf("searching for larger errors...\n");	
 /* now check a range of values 1009 is prime */
 for(x=0;x<65536.5;x+=1.0/1009.0)
 	{yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;		  	
	 	  	}
	 	}
 	}
 rel_e=fabsl(me/me_y);
 printf("%s vs %s max abs difference found was %g (=%.3f bits) at %.20g\n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)),(double)me,(double)(-log2l(rel_e)),me_x);
 printf(" max relative error found was %g (=%.3f bits) at %.20g\n",(double)m_rel_e,(double)(-log2l(m_rel_e)),m_rel_e_x); 
 printf(" So accuracy of approximation is %.3f bits\n",exp2((double)DBL_MANT_DIG-(double)(-log2l(m_rel_e))));
 // calculate ulp according to "accuracy" paper
 double nexty=nextafter(m_rel_e_y,DBL_MAX);// next floating point number up
 printf("\n1ulp at x=%.20g, y=%.20g is %g\n",m_rel_e_x,m_rel_e_y,nexty-m_rel_e_y);
 printf("So accuracy of approximation vs %s is %.3f ulp\n",TO_STRING(ref_sqrt_LD(x)),(double)(m_rel_a/(long double)(nexty-m_rel_e_y )));
 if(errs!=0) printf("Warning: %d errors found in NAN/INF handling\n",errs);

 // now test log(x)
#undef ref_sqrt_LD 
#undef ref_sqrt_DBL
#undef sqrt_under_test 
#define ref_sqrt_LD(x) __builtin_logl(x) /* function thats taken to be "exact" - either sqrtl(x) or __builtin_sqrtl(x) - must take a long double argument (x) and return a long double */
#define ref_sqrt_DBL(x) __builtin_log(x) /* function thats taken to be "exact" - either sqrt(x) or __builtin_sqrt(x) - must take a ldouble argument (x) and return a double - only used for NAM,INF */
#define sqrt_under_test(x) cr_log(x) /* function under test - normally here that is cr_sqrt(x), but can also test library sqrt - must take a double argument (x) and return a double */ 
 errs=0;
 // first check is just used to initialise variables to sensible values, actual value of x used is not critical, but it must be a valid argument for the function under test 
 x=0x1.fffffffffffffp-1;
 printf("\nTesting %s vs %s : \n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)));
 y=ref_sqrt_LD((long double)x);
 y_cr=sqrt_under_test(x);
 e=fabsl((long double)y_cr-y); // error
 me=e;// use this value as "baseline" for wider test below
 me_x=x;
 me_y=y;
 rel_e=fabsl(e/y);// rel error [ relative to correct result]
 m_rel_e=rel_e;
 m_rel_e_x=x;
 m_rel_e_y=y;
 m_rel_a=e;	
 rel_e=fabsl(me/y);
 printf("Checking values in special_vals[]...\n");
 for(int i=0;i<nos_elements_in(special_vals);++i)
	{
	 x=special_vals[i];
	 yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;		  	
	 	  	}
	 	}	 	
	
	}
 printf("searching for larger errors...\n");	
 /* now check a range of values 1009 is prime */
 for(x=0;x<65536.5;x+=1.0/1009.0)
 	{yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;		  	
	 	  	}
	 	}
 	}
 rel_e=fabsl(me/me_y);
 printf("%s vs %s max abs difference found was %g (=%.3f bits) at %.20g\n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)),(double)me,(double)(-log2l(rel_e)),me_x);
 printf(" max relative error found was %g (=%.3f bits) at %.20g\n",(double)m_rel_e,(double)(-log2l(m_rel_e)),m_rel_e_x); 
 printf(" So accuracy of approximation is %.3f bits\n",exp2((double)DBL_MANT_DIG-(double)(-log2l(m_rel_e))));
 // calculate ulp according to "accuracy" paper
 nexty=nextafter(m_rel_e_y,DBL_MAX);// next floating point number up
 printf("\n1ulp at x=%.20g, y=%.20g is %g\n",m_rel_e_x,m_rel_e_y,nexty-m_rel_e_y);
 printf("So accuracy of approximation vs %s is %.3f ulp\n",TO_STRING(ref_sqrt_LD(x)),(double)(m_rel_a/(long double)(nexty-m_rel_e_y )));
 if(errs!=0) printf("Warning: %d errors found in NAN/INF handling\n",errs); 

 // Now test exp(x)
#undef ref_sqrt_LD 
#undef ref_sqrt_DBL
#undef sqrt_under_test 
#define ref_sqrt_LD(x) __builtin_expl(x) /* function thats taken to be "exact" - either sqrtl(x) or __builtin_sqrtl(x) - must take a long double argument (x) and return a long double */
#define ref_sqrt_DBL(x) __builtin_exp(x) /* function thats taken to be "exact" - either sqrt(x) or __builtin_sqrt(x) - must take a ldouble argument (x) and return a double - only used for NAM,INF */
#define sqrt_under_test(x) cr_exp(x) /* function under test - normally here that is cr_sqrt(x), but can also test library sqrt - must take a double argument (x) and return a double */ 
 errs=0;
 // first check is just used to initialise variables to sensible values, actual value of x used is not critical, but it must be a valid argument for the function under test  
 x=0x1.fffffffffffffp-1;
 printf("\nTesting %s vs %s : \n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)));
 y=ref_sqrt_LD((long double)x);
 y_cr=sqrt_under_test(x);
 e=fabsl((long double)y_cr-y); // error
 me=e;// use this value as "baseline" for wider test below
 me_x=x;
 me_y=y;
 rel_e=fabsl(e/y);// rel error [ relative to correct result]
 m_rel_e=rel_e;
 m_rel_e_x=x;
 m_rel_e_y=y;
 m_rel_a=e;	
 rel_e=fabsl(me/y);
 printf("Checking values in special_vals[]...\n");
 for(int i=0;i<nos_elements_in(special_vals);++i)
	{
	 x=special_vals[i];
	 yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y;
			 m_rel_a=e;	 	  	
	 	  	}
	 	}	 	
	
	}
 printf("searching for larger errors...\n");	
 /* now check a range of values 1009 is prime */
 for(x=0;x<65536.5;x+=1.0/1009.0)
 	{yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;		  	
	 	  	}
	 	}
 	}
 rel_e=fabsl(me/me_y);
 printf("%s vs %s max abs difference found was %g (=%.3f bits) at %.20g\n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)),(double)me,(double)(-log2l(rel_e)),me_x);
 printf(" max relative error found was %g (=%.3f bits) at %.20g\n",(double)m_rel_e,(double)(-log2l(m_rel_e)),m_rel_e_x); 
 printf(" So accuracy of approximation is %.3f bits\n",exp2((double)DBL_MANT_DIG-(double)(-log2l(m_rel_e))));
 // calculate ulp according to "accuracy" paper
 nexty=nextafter(m_rel_e_y,DBL_MAX);// next floating point number up
 printf("\n1ulp at x=%.20g, y=%.20g is %g\n",m_rel_e_x,m_rel_e_y,nexty-m_rel_e_y);
 printf("So accuracy of approximation vs %s is %.3f ulp\n",TO_STRING(ref_sqrt_LD(x)),(double)(m_rel_a/(long double)(nexty-m_rel_e_y )));
 if(errs!=0) printf("Warning: %d errors found in NAN/INF handling\n",errs); 


 // check sin(x) vs f128 sin(x) - had to do this as using sinl() gave big errors!
#undef ref_sqrt_LD 
#undef ref_sqrt_DBL
#undef sqrt_under_test 
#define ref_sqrt_LD(x) sinq(x) /* function thats taken to be "exact" - either sqrtl(x) or __builtin_sqrtl(x) - must take a long double argument (x) and return a long double */
//#define ref_sqrt_LD(x) sinl(x) /* gives big max. error */
#define ref_sqrt_DBL(x) __builtin_sin(x) /* function thats taken to be "exact" - either sqrt(x) or __builtin_sqrt(x) - must take a ldouble argument (x) and return a double - only used for NAM,INF */
#define sqrt_under_test(x) cr_sin(x) /* function under test - normally here that is cr_sqrt(x), but can also test library sqrt - must take a double argument (x) and return a double */ 
 errs=0;
 __float128 y128,e128,m_rel_a128  ;
 // first check is just used to initialise variables to sensible values, actual value of x used is not critical, but it must be a valid argument for the function under test  
 x=0x1.fffffffffffffp-1;
 printf("\nTesting %s vs %s : \n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)));
 y128=ref_sqrt_LD((__float128)x);
 y=(long double)y128;
 y_cr=sqrt_under_test(x);
 e128=fabsq((__float128)y_cr-y); // error
 e=e128;
 me=e;// use this value as "baseline" for wider test below
 me_x=x;
 me_y=y;
 rel_e=fabsl(e/y);// rel error [ relative to correct result]
 m_rel_e=rel_e;
 m_rel_e_x=x;
 m_rel_e_y=y;
 m_rel_a=e;	
 m_rel_a128=e128;
 rel_e=fabsl(me/y);
 printf("Checking values in special_vals[]...\n");
 for(int i=0;i<nos_elements_in(special_vals);++i)
	{
	 x=special_vals[i];
	 yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 if(x>1000.0 || x<-1000.0) continue; // x is in radians so large values are "silly"
	 y128=ref_sqrt_LD((__float128)x);
	 y=(long double)y128;
	 e128=fabsq((__float128)y_cr-y); // error
	 e=e128;
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 	
			 m_rel_a=e;	  
			 m_rel_a128=e128;	
	 	  	}
	 	}	 	
	
	}
 printf("searching for larger errors...\n");	
 /* now check a range of values 1062599 is prime 65.5365 is a little above 20*PI */
 for(x=0;x<65.5365;x+=1.0/1062599.0)
 	{yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y128=ref_sqrt_LD((__float128)x);
	 y=(long double)y128;
	 e128=fabsq((__float128)y_cr-y); // error
	 e=e128;
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 	
			 m_rel_a=e;	  
			 m_rel_a128=e128;	
	 	  	}
	 	}
 	}
 rel_e=fabsl(me/me_y);
 printf("%s vs %s max abs difference found was %g (=%.3f bits) at %.20g\n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)),(double)me,(double)(-log2l(rel_e)),me_x);
 my_printf(" max relative error found was %g (=%.3f bits) at %.20g where abs error was %.20Lg and expected y was %.20g\n",(double)m_rel_e,(double)(-log2l(m_rel_e)),m_rel_e_x,(long double)m_rel_a128,m_rel_e_y); 
 printf(" So accuracy of approximation is %.3f bits\n",exp2((double)DBL_MANT_DIG-(double)(-log2l(m_rel_e))));
 // calculate ulp according to "accuracy" paper
 nexty=nextafter(m_rel_e_y,DBL_MAX);// next floating point number up
 printf("\n1ulp at x=%.20g, y=%.20g is %g\n",m_rel_e_x,m_rel_e_y,nexty-m_rel_e_y);
 printf("So accuracy of approximation vs %s is %.3f ulp\n",TO_STRING(ref_sqrt_LD(x)),(double)(m_rel_a/(long double)(nexty-m_rel_e_y )));
 if(errs!=0) printf("Warning: %d errors found in NAN/INF handling\n",errs);

 // check cos(x) vs f128 cos(x) - had to do this as using cosl() gave big errors!
#undef ref_sqrt_LD 
#undef ref_sqrt_DBL
#undef sqrt_under_test 
#define ref_sqrt_LD(x) cosq(x) /* function thats taken to be "exact" - either sqrtl(x) or __builtin_sqrtl(x) - must take a long double argument (x) and return a long double */
#define ref_sqrt_DBL(x) __builtin_cos(x) /* function thats taken to be "exact" - either sqrt(x) or __builtin_sqrt(x) - must take a ldouble argument (x) and return a double - only used for NAM,INF */
#define sqrt_under_test(x) cr_cos(x) /* function under test - normally here that is cr_sqrt(x), but can also test library sqrt - must take a double argument (x) and return a double */ 
 errs=0;
 // first check is just used to initialise variables to sensible values, actual value of x used is not critical, but it must be a valid argument for the function under test  
 x=0x1.fffffffffffffp-1;
 printf("\nTesting %s vs %s : \n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)));
 y128=ref_sqrt_LD((__float128)x);
 y=(long double)y128;
 y_cr=sqrt_under_test(x);
 e128=fabsq((__float128)y_cr-y); // error
 e=e128;
 me=e;// use this value as "baseline" for wider test below
 me_x=x;
 me_y=y;
 rel_e=fabsl(e/y);// rel error [ relative to correct result]
 m_rel_e=rel_e;
 m_rel_e_x=x;
 m_rel_e_y=y;
 m_rel_a=e;	
 m_rel_a128=e128;
 rel_e=fabsl(me/y); 
 printf("Checking values in special_vals[]...\n");
 for(int i=0;i<nos_elements_in(special_vals);++i)
	{
	 x=special_vals[i];
	 yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 if(x>1000.0 || x<-1000.0) continue; // x is in radians so large values are "silly"
	 y128=ref_sqrt_LD((__float128)x);
	 y=(long double)y128;
	 e128=fabsq((__float128)y_cr-y); // error
	 e=e128;
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 	
			 m_rel_a=e;	  
			 m_rel_a128=e128;	
	 	  	}
	 	}	 	
	
	}
 printf("searching for larger errors...\n");	
 /* now check a range of values 1062599 is prime 65.5365 is a little above 20*PI */
 for(x=0;x<65.5365;x+=1.0/1062599.0)
 	{yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y128=ref_sqrt_LD((__float128)x);
	 y=(long double)y128;
	 e128=fabsq((__float128)y_cr-y); // error
	 e=e128;
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 	
			 m_rel_a=e;	  
			 m_rel_a128=e128;	
	 	  	}
	 	}
 	}
 rel_e=fabsl(me/me_y);
 printf("%s vs %s max abs difference found was %g (=%.3f bits) at %.20g\n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)),(double)me,(double)(-log2l(rel_e)),me_x);
 my_printf(" max relative error found was %g (=%.3f bits) at %.20g where abs error was %.20Lg and expected y was %.20g\n",(double)m_rel_e,(double)(-log2l(m_rel_e)),m_rel_e_x,(long double)m_rel_a128,m_rel_e_y); 
 printf(" So accuracy of approximation is %.3f bits\n",exp2((double)DBL_MANT_DIG-(double)(-log2l(m_rel_e))));
 // calculate ulp according to "accuracy" paper
 nexty=nextafter(m_rel_e_y,DBL_MAX);// next floating point number up
 printf("\n1ulp at x=%.20g, y=%.20g is %g\n",m_rel_e_x,m_rel_e_y,nexty-m_rel_e_y);
 printf("So accuracy of approximation vs %s is %.3f ulp\n",TO_STRING(ref_sqrt_LD(x)),(double)(m_rel_a/(long double)(nexty-m_rel_e_y )));
 if(errs!=0) printf("Warning: %d errors found in NAN/INF handling\n",errs);


 // check atan2(y,x)
#undef ref_sqrt_LD 
#undef ref_sqrt_DBL
#undef sqrt_under_test
// in the line below __builtin_sin(x) & cos as arguments are correct as we want to compare atan2() with identical argument
#define ref_sqrt_LD(x) __builtin_atan2l((long double)__builtin_sin(x),(long double)__builtin_cos(x)) /* function thats taken to be "exact" - either sqrtl(x) or __builtin_sqrtl(x) - must take a long double argument (x) and return a long double */
#define ref_sqrt_DBL(x) __builtin_atan2(__builtin_sin(x),__builtin_cos(x)) /* function thats taken to be "exact" - either sqrt(x) or __builtin_sqrt(x) - must take a ldouble argument (x) and return a double - only used for NAM,INF */
#define sqrt_under_test(x) cr_atan2(__builtin_sin(x),__builtin_cos(x)) /* function under test - normally here that is cr_sqrt(x), but can also test library sqrt - must take a double argument (x) and return a double */ 
 errs=0;
 // first check is just used to initialise variables to sensible values, actual value of x used is not critical, but it must be a valid argument for the function under test  
 x=0x1.fffffffffffffp-1;
 printf("\nTesting %s vs %s : \n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)));
 y=ref_sqrt_LD((long double)x);
 y_cr=sqrt_under_test(x);
 e=fabsl((long double)y_cr-y); // error
 me=e;// use this value as "baseline" for wider test below
 me_x=x;
 me_y=y;
 rel_e=fabsl(e/y);// rel error [ relative to correct result]
 m_rel_e=rel_e;
 m_rel_e_x=x;
 m_rel_e_y=y;
 m_rel_a=e;	
 rel_e=fabsl(me/y);
 printf("Checking values in special_vals[]...\n");
 for(int i=0;i<nos_elements_in(special_vals);++i)
	{
	 x=special_vals[i];
	 yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 if(x>1000.0 || x<-1000.0) continue; // x is in radians so large values are "silly"
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;		  	
	 	  	}
	 	}	 	
	
	}
 printf("searching for larger errors...\n");	
 /* now check a range of values 1062599 is prime 65.5365 is a little above 20*PI */
 for(x=0;x<65.5365;x+=1.0/1062599.0)
 	{yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;		  	
	 	  	}
	 	}
 	}
 rel_e=fabsl(me/me_y);
 printf("%s vs %s max abs difference found was %g (=%.3f bits) at %.20g\n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)),(double)me,(double)(-log2l(rel_e)),me_x);
 printf(" max relative error found was %g (=%.3f bits) at %.20g\n",(double)m_rel_e,(double)(-log2l(m_rel_e)),m_rel_e_x); 
 printf(" So accuracy of approximation is %.3f bits\n",exp2((double)DBL_MANT_DIG-(double)(-log2l(m_rel_e))));
 // calculate ulp according to "accuracy" paper
 nexty=nextafter(m_rel_e_y,DBL_MAX);// next floating point number up
 printf("\n1ulp at x=%.20g, y=%.20g is %g\n",m_rel_e_x,m_rel_e_y,nexty-m_rel_e_y);
 printf("So accuracy of approximation vs %s is %.3f ulp\n",TO_STRING(ref_sqrt_LD(x)),(double)(m_rel_a/(long double)(nexty-m_rel_e_y )));
 if(errs!=0) printf("Warning: %d errors found in NAN/INF handling\n",errs); 
 
 // check pow(x,y) for x^3
#undef ref_sqrt_LD 
#undef ref_sqrt_DBL
#undef sqrt_under_test
#define ref_sqrt_LD(x) ((long double)x*(long double)x*(long double)x) /* __builtin_powl(x,3) function thats taken to be "exact" - either sqrtl(x) or __builtin_sqrtl(x) - must take a long double argument (x) and return a long double */
#define ref_sqrt_DBL(x) __builtin_pow(x,3) /* function thats taken to be "exact" - either sqrt(x) or __builtin_sqrt(x) - must take a ldouble argument (x) and return a double - only used for NAM,INF */
#define sqrt_under_test(x) cr_pow(x,3) /* function under test - normally here that is cr_sqrt(x), but can also test library sqrt - must take a double argument (x) and return a double */ 
 errs=0;
 // first check is just used to initialise variables to sensible values, actual value of x used is not critical, but it must be a valid argument for the function under test  
 x=0x1.fffffffffffffp-1;
 printf("\nTesting %s vs %s : \n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)));
 y=ref_sqrt_LD((long double)x);
 y_cr=sqrt_under_test(x);
 e=fabsl((long double)y_cr-y); // error
 me=e;// use this value as "baseline" for wider test below
 me_x=x;
 me_y=y;
 rel_e=fabsl(e/y);// rel error [ relative to correct result]
 m_rel_e=rel_e;
 m_rel_e_x=x;
 m_rel_e_y=y;
 m_rel_a=e;	
 rel_e=fabsl(me/y);
 printf("Checking values in special_vals[]...\n");
 for(int i=0;i<nos_elements_in(special_vals);++i)
	{
	 x=special_vals[i];
	 yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;		  	
	 	  	}
	 	}	 	
	
	}
 printf("searching for larger errors...\n");	
 /* now check a range of values 1009 is prime */
 for(x=0;x<65536.5;x+=1.0/1009.0)
 	{yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;		  	
	 	  	}
	 	}
 	}
 rel_e=fabsl(me/me_y);
 printf("%s vs %s max abs difference found was %g (=%.3f bits) at %.20g\n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)),(double)me,(double)(-log2l(rel_e)),me_x);
 printf(" max relative error found was %g (=%.3f bits) at %.20g\n",(double)m_rel_e,(double)(-log2l(m_rel_e)),m_rel_e_x); 
 printf(" So accuracy of approximation is %.3f bits\n",exp2((double)DBL_MANT_DIG-(double)(-log2l(m_rel_e))));
 // calculate ulp according to "accuracy" paper
 nexty=nextafter(m_rel_e_y,DBL_MAX);// next floating point number up
 printf("\n1ulp at x=%.20g, y=%.20g is %g\n",m_rel_e_x,m_rel_e_y,nexty-m_rel_e_y);
 printf("So accuracy of approximation vs %s is %.3f ulp\n",TO_STRING(ref_sqrt_LD(x)),(double)(m_rel_a/(long double)(nexty-m_rel_e_y )));
 if(errs!=0) printf("Warning: %d errors found in NAN/INF handling\n",errs); 

 // check pow(x,y) for x^0.5 (sqrt(x)
#undef ref_sqrt_LD 
#undef ref_sqrt_DBL
#undef sqrt_under_test
#define ref_sqrt_LD(x) sqrtl((long double)x) /* __builtin_powl(x,0.5) function thats taken to be "exact" - either sqrtl(x) or __builtin_sqrtl(x) - must take a long double argument (x) and return a long double */
#define ref_sqrt_DBL(x) __builtin_pow(x,0.5) /* function thats taken to be "exact" - either sqrt(x) or __builtin_sqrt(x) - must take a ldouble argument (x) and return a double - only used for NAM,INF */
#define sqrt_under_test(x) cr_pow(x,0.5) /* function under test - normally here that is cr_sqrt(x), but can also test library sqrt - must take a double argument (x) and return a double */ 
 errs=0;
 // first check is just used to initialise variables to sensible values, actual value of x used is not critical, but it must be a valid argument for the function under test  
 x=0x1.fffffffffffffp-1;
 printf("\nTesting %s vs %s : \n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)));
 y=ref_sqrt_LD((long double)x);
 y_cr=sqrt_under_test(x);
 e=fabsl((long double)y_cr-y); // error
 me=e;// use this value as "baseline" for wider test below
 me_x=x;
 me_y=y;
 rel_e=fabsl(e/y);// rel error [ relative to correct result]
 m_rel_e=rel_e;
 m_rel_e_x=x;
 m_rel_e_y=y;
 m_rel_a=e;	
 rel_e=fabsl(me/y); 
 printf("Checking values in special_vals[]...\n");
 for(int i=0;i<nos_elements_in(special_vals);++i)
	{
	 x=special_vals[i];
	 yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;		  	
	 	  	}
	 	}	 	
	
	}
 printf("searching for larger errors...\n");	
 /* now check a range of values 1009 is prime */
 for(x=0;x<65536.5;x+=1.0/1009.0)
 	{yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;		  	
	 	  	}
	 	}
 	}
 rel_e=fabsl(me/me_y);
 printf("%s vs %s max abs difference found was %g (=%.3f bits) at %.20g\n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)),(double)me,(double)(-log2l(rel_e)),me_x);
 printf(" max relative error found was %g (=%.3f bits) at %.20g\n",(double)m_rel_e,(double)(-log2l(m_rel_e)),m_rel_e_x); 
 printf(" So accuracy of approximation is %.3f bits\n",exp2((double)DBL_MANT_DIG-(double)(-log2l(m_rel_e))));
 // calculate ulp according to "accuracy" paper
 nexty=nextafter(m_rel_e_y,DBL_MAX);// next floating point number up
 printf("\n1ulp at x=%.20g, y=%.20g is %g\n",m_rel_e_x,m_rel_e_y,nexty-m_rel_e_y);
 printf("So accuracy of approximation vs %s is %.3f ulp\n",TO_STRING(ref_sqrt_LD(x)),(double)(m_rel_a/(long double)(nexty-m_rel_e_y )));
 if(errs!=0) printf("Warning: %d errors found in NAN/INF handling\n",errs); 
 
 // check pow(x,y) for x^-0.5 (1/sqrt(x))
#undef ref_sqrt_LD 
#undef ref_sqrt_DBL
#undef sqrt_under_test
#define ref_sqrt_LD(x) sqrtl(1.0L/(long double)x) /* __builtin_powl(x,-0.5) function thats taken to be "exact" - either sqrtl(x) or __builtin_sqrtl(x) - must take a long double argument (x) and return a long double */
#define ref_sqrt_DBL(x) __builtin_pow(x,-0.5) /* function thats taken to be "exact" - either sqrt(x) or __builtin_sqrt(x) - must take a ldouble argument (x) and return a double - only used for NAM,INF */
#define sqrt_under_test(x) cr_pow(x,-0.5) /* function under test - normally here that is cr_sqrt(x), but can also test library sqrt - must take a double argument (x) and return a double */ 
 errs=0;
 // first check is just used to initialise variables to sensible values, actual value of x used is not critical, but it must be a valid argument for the function under test  
 x=0x1.fffffffffffffp-1;
 printf("\nTesting %s vs %s : \n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)));
 y=ref_sqrt_LD((long double)x);
 y_cr=sqrt_under_test(x);
 e=fabsl((long double)y_cr-y); // error
 me=e;// use this value as "baseline" for wider test below
 me_x=x;
 me_y=y;
 rel_e=fabsl(e/y);// rel error [ relative to correct result]
 m_rel_e=rel_e;
 m_rel_e_x=x;
 m_rel_e_y=y;
 m_rel_a=e;	
 rel_e=fabsl(me/y); 
 printf("Checking values in special_vals[]...\n");
 for(int i=0;i<nos_elements_in(special_vals);++i)
	{
	 x=special_vals[i];
	 yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;		  	
	 	  	}
	 	}	 	
	
	}
 printf("searching for larger errors...\n");	
 /* now check a range of values 1009 is prime */
 for(x=0;x<65536.5;x+=1.0/1009.0)
 	{yr=ref_sqrt_DBL(x);
	 y_cr=sqrt_under_test(x);
	 if(isnan(yr) || isinf(yr))
	 	{
		 if(isnan(y_cr)!=isnan(yr) && isinf(y_cr)!=isinf(yr))
	 		{++errs;
	 		 printf("Error: %s returns %g but should be (from %s) %g\n",TO_STRING(sqrt_under_test(x)),y_cr,TO_STRING(ref_sqrt_DBL(x)),yr);
	 		}
	 	 continue;// allows indenting below to match next loop below
	 	}
	 y=ref_sqrt_LD((long double)x);
 	 e=fabsl((long double)y_cr-y); // abs error
 	 if(e>me)
 	 	{me=e;
 	 	 me_x=x;
 	 	 me_y=y;
 	 	}
 	  if(fabs(yr)>=DBL_MIN) /* we need to avoid dividing by denorms as that reduces the dynamic range and we end up with a smaller result in ulp terms */
 	  	{
	 	 rel_e=fabsl(e/y); // rel error [ relative to correct result]
	 	 if(rel_e>m_rel_e)
	 	  	{
			 m_rel_e=rel_e;
			 m_rel_e_x=x; 
			 m_rel_e_y=y; 
			 m_rel_a=e;		  	
	 	  	}
	 	}
 	}
 rel_e=fabsl(me/me_y);
 printf("%s vs %s max abs difference found was %g (=%.3f bits) at %.20g\n",TO_STRING(sqrt_under_test(x)),TO_STRING(ref_sqrt_LD(x)),(double)me,(double)(-log2l(rel_e)),me_x);
 printf(" max relative error found was %g (=%.3f bits) at %.20g\n",(double)m_rel_e,(double)(-log2l(m_rel_e)),m_rel_e_x); 
 printf(" So accuracy of approximation is %.3f bits\n",exp2((double)DBL_MANT_DIG-(double)(-log2l(m_rel_e))));
 // calculate ulp according to "accuracy" paper
 nexty=nextafter(m_rel_e_y,DBL_MAX);// next floating point number up
 printf("\n1ulp at x=%.20g, y=%.20g is %g\n",m_rel_e_x,m_rel_e_y,nexty-m_rel_e_y);
 printf("So accuracy of approximation vs %s is %.3f ulp\n",TO_STRING(ref_sqrt_LD(x)),(double)(m_rel_a/(long double)(nexty-m_rel_e_y )));
 if(errs!=0) printf("Warning: %d errors found in NAN/INF handling\n",errs);  


}
