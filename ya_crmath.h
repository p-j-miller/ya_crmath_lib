/* Header file for maths functions in ya_math_lib
 Written by Peter Miller, 1/9/2026
 
 Note if YA_MATH_LIB_REPLACE  is defined cr_xxx routines will replace standard C functions (via #defines), making it easy to swap existing code to use these functions
 e.g use:
   #define YA_MATH_LIB_REPLACE  // specify we want to replace standard C functions MUST come before #include on line below to be effective
   #include "../ya_math_lib/ya_math.h" // for cc_xxx maths functions 
   
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
	double cr_sqrt(double x);	
	
 #ifdef YA_CRMATH_LIB_REPLACE  /* if defined cr_xxx routines will replace standard C functions */
  #define atan2(y,x) cr_atan2(y,x)
  #define cos(x) cr_cos(x)
  #define exp(x) cr_exp(x)
  #define log(x) cr_log(x)
  #define pow(x,y) cr_pow(x,y)
  #define sin(x) cr_sin(x)
  #define  sqrt(x) cr_sqrt(x)
 #endif
 
 #ifdef __cplusplus
    }
 #endif

#endif // ifndef _ya_math_lib_h