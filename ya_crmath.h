/* Header file for maths functions in ya_math_lib
 Written by Peter Miller, 1/9/2026
 
 Note if YA_MATH_LIB_REPLACE  is defined cr_xxx and ya_xxx routines will replace standard C functions (via #defines), making it easy to swap existing code to use these functions
 e.g use:
   #define YA_MATH_LIB_REPLACE  // specify we want to replace standard C functions MUST come before #include on line below to be effective
   #include "../ya_math_lib/ya_math.h" // for cc_xxx maths functions 
  
Version 1.1 
	adds "derived" functions ya_tan(), ya_asin(), ya_acos() which are not correctly rounded and have errors > 0.500ulp (very similar to the errors with the UCRT), but will give consistent results between different compilers/operating systems / processors
	ya_sin(),ya_cos(),ya_exp() and ya_log() functions added which use a hybrid of the algorithms in llvm maths functions and the CORE MATHS functions to give faster execution speed  (and less sensitivity to the use of the FMA instruction except in ll_exp) with identical accuracy
 
  These functions optionally support errno - this can be configured in ya_crmath_config.h (by default they are configured not to use errno)
  
Execution times (relative to "built in double" function) 1v1b [ so values >1* take longer]. 
Note "no fma" means use of hardware fma instruction was disabled, and it was emulated in software/

Function no fma	-mfma	comments
 cr_sqrt		1.00* 	1.00	(using __builtin_sqrt(x) )
 cr_sqrt		4.01* 	4.02	(using C implementation in cr_sqrt.c )
 cr_log			2.85*	1.30
 cr_exp			1.78*	1.21
 cr_sin			6.25*	3.45
 cr_cos			5.77*	3.27
 cr_atan2		1.35*	1.20
 cr_power	2.58*	1.21	power(x,3)
					2.94*	1.37	power(x,0.5)
					2.65*	1.19	power(x,-0.5)
				
ya_xxx functions:
 ya_log			1.52*	1.19 
 ya_exp		1.90*	1.20
 ya_sin			2.63*	1.91
 ya_cos		2.73*	1.96
				
derived functions:			error without FMA 	error using FMA		error using UCRT function
 ya_tan			3.54 *	2.43   	1.212ulp					1.212ulp				1.32 ulp
 ya_asin		3.65*	2.95		1.099 ulp  				0.869 ulp				0.960 ulp
 ya_acos 		3.51*	2.83		1.311 ulp  				1.195 ulp				0.932 ulp
 sin/cos		0.54*	0.41		1.7e9 ulp using recurrence below  ("sum" in test program is accurate, this error is near x=2pi where sin(x) is 0 so the size of 1ulp is 2.65-23, the abs error is only  4.6e-14 )

 
 Numerical Recipes says "If your programs runtime is dominated by evaluating trigonometric functions you are probably doing something wrong"!
 It then offers the following recurrence to step sin/cos by a constant increment d from the previous value x (if initial x=0, sin(0)=0 and cos(0)=1)
 
	 cos(x+d)=cos(x)-(a*cos(x)+b*sin(x))
	 sin(x+d)=sin(x)-(a*sin(x)-b*cos(x))
 
 where a and b are precomputed as:
	 a=2*(sin(d/2))^2
	 b=sin(d)
 => Numerical Recipes says the reason for this implementation is that a & b do not lose significance if the increment d is small.
 
 
with gcc  :
	 -msse4.2
	-mfpmath=sse
	-mfma
	-mavx2    (this is not essential, but if fma is present then avx2 is also present so we might as well let the compiler use it! )
	-O3
	-fno-math-errno
	-DNDEBUG
	
with C++ Builder 13 CE 64-bit Modern (Clang 20) project options:
	C++ Compiler - General compilation
		Instruction set : AVX2   (sets -target-feature +avx2)
	Advanced - floating point
		require maths functions to indicate errors : false
	Advanced - Other options
		Additional options to pass to the compiler : -ffp-contract=fast -target-feature +fma
	Optimizations
		Generate most optimised code (-O3) : true
		
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
#ifndef _ya_math_lib_h
#define _ya_math_lib_h

 #ifdef __cplusplus
 extern "C" {
 #endif 
	double cr_atan2(double y0, double x0);
	double cr_cos(double x);
	double cr_exp(double x);
	double cr_log(double x);
	double cr_pow(double x, double y);
	double cr_sin(double x);
#if defined(__has_builtin) && __has_builtin(__builtin_sqrt)  && defined(__SSE2__) /*  is builtin for gcc and clang - we also need at least SSE2 as that introduced the SQRTSD instruction that gives an accurate sqrt of a double */
 #define cr_sqrt(x) __builtin_sqrt(x) /* https://members.loria.fr/PZimmermann/papers/accuracy.pdf says built in functions on gcc and clang are accurate to 0.500 ulp and are correctly rounded [ and they are significantly faster than C code in cr_sqrt.c as they typically compile to 1 assembly instruction ] */
 // note that for gcc -fno-math-errno  also gives a significant speedup where errno is not needed.
 // These functions optionally support errno - this can be configured in ya_crmath_config.h
#else /* using C code */	
	double cr_sqrt(double x);	
#endif	

// functions based on llvm algorithms (actually they are a hybrid of the llvm algorithms and the cr_xxx algorithms above)
// These are also accurate to 0.500ulp and are faster than the corresponding cr_xxx functions above in the test program with an i3-10100 and winlibs gcc 16.2.0 -m64 -mfma 
// According to "accuracy.pdf", the llvm atan2() and power() functions have an error >0.500 ulp so were not considered for use here
  double ya_sin(double x);
  double ya_cos(double x);
  double ya_log(double x);
  double ya_exp(double x);
 
// "derived" functions 
// tan=sin/cos by definition 
// max error found by test program is 1.21 ulp
// according to the paper "accuracy.pdf", the UCRT's tan(x)  is 1.32ulp
inline static double ya_tan(double x)
{
 return ya_sin(x)/ya_cos(x); // use ya_sin/cos as they are significantly faster than cr_sin/cos
}

 // asin(x) = atan2 (x, sqrt ((1.0 + x) * (1.0 - x)))
 // (1+x)*(1-x) is more accurate than 1-x*x.
 // fma(x,-x,1.0) is slightly more accurate, but a little slower when fma is a single instruction, and quite a lot slower if fma is a function.
 // note asin(x) only valid for -1<=x<=1
// max error found by test program without fma is 1.10 ulp
// max error found by test program with fma is 0.87 ulp
// according to the paper "accuracy.pdf", the UCRT's asin(x)  is 0.960ulp
inline static double ya_asin(double x)
{
 //return cr_atan2(x,cr_sqrt((1.0 + x) * (1.0 - x)));
 return cr_atan2(x,cr_sqrt(__builtin_fma(x,-x,1.0)) );
}

 // acos(x) = atan2 (sqrt ((1.0 + x) * (1.0 - x)), x)
 // (1+x)*(1-x) is more accurate than 1-x*x.
 // fma(x,-x,1.0) is slightly more accurate, but a little slower when fma is a single instruction, and quite a lot slower if fma is a function.
 // note acos(x) only valid for -1<=x<=1
 // max error found by test program without fma is 1.31 ulp
 // max error found by test program with fma is 1.20  ulp
 // according to the paper "accuracy.pdf", the UCRT's acos(x)  is 0.932ulp
inline static double ya_acos(double x)
{
 //return cr_atan2(cr_sqrt((1.0 + x) * (1.0 - x)),x);
 return cr_atan2(cr_sqrt(__builtin_fma(x,-x,1.0)),x);
}
 
 #ifdef YA_CRMATH_LIB_REPLACE  /* if defined cr_xxx,ll_xxx and ya_xxx routines will replace standard C functions */
  #if 1 /* use the ya_xxx functions instead of cr_xxx functions , this should give a small speedup while giving exactly the same values */
  #define sin(x) ya_sin(x)
  #define cos(x) ya_cos(x)
  #define log(x) ya_log(x)
  #define exp(x) ya_exp(x)
#else /* use cr_xxx functions, slightly slower, but should give identical results */
  #define sin(x) cr_sin(x)	
  #define cos(x) cr_cos(x)
  #define log(x) cr_log(x)  
  #define exp(x) cr_exp(x)  
#endif  
  #define atan2(y,x) cr_atan2(y,x)
  #define pow(x,y) cr_pow(x,y)
  #define  sqrt(x) cr_sqrt(x)
  // note the next few functions give an error that is worse than 0.500 ulp, (see above but <=1.31ulp so still very small and almost identical to the <=1.32 with the UCRT)  but their results are repeatable
  #define tan(x) ya_tan(x)
  #define asin(x) ya_asin(x)
  #define acos(x) ya_acos(x)
 #endif
 
 #ifdef __cplusplus
    }
 #endif

#endif // ifndef _ya_math_lib_h