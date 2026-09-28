# ya_crmath_lib
ya_crmath_lib is a portable C maths library that should always give the same results independent of the C compiler / processor / operating system.
Functions are named cr_* which is a syntax proposed for a future C standard, "cr_" standing for "correct rounding".
Note this library contains a small subset of all the C maths functions - just enough for wmawk2.c initially (sin,cos,atan2,exp,log,sqrt).

It was created as the wmawk2 test program showed up differences between the windows 32 and 64 bit executables with sin() giving slightly different results.
According to https://blog.r-project.org/2026/01/06/debugging-sensitivity-to-c-math-library-and-mingw-w64-v12/  
	“The change in mingw-w64 v12 switched about 90 math functions from internal implementations to UCRT,
	 the Windows C runtime provided by Microsoft and shipped with the OS.
	 The goal of the change was improving performance and reducing maintenance cost of mingw-w64.”
 It also says “One should keep in mind, as reported in one of the earlier posts, that UCRT results are OS-specific, they differ between Windows versions.”.  
There is more detail in https://blog.r-project.org/2025/04/24/sensitivity-to-c-math-library-and-mingw-w64-v12/  and https://members.loria.fr/PZimmermann/papers/accuracy.pdf 

The last reference (accuracy.pdf) says  "The square root function seems to be correctly rounded for all (tested) libraries, as required by IEEE 754",
I selected the sqrt function from MUSL-1.2.6 library (one of the ones tested in "accuracy.pdf") as it does its internal calculations using integer maths, rather than floating point maths which many of the others I looked at did,
on the basis that cpu's tend to be faster on integer code.  For example the sqrt function in crlibm (see e.g.  https://github.com/taschini/crlibm/blob/master/docs/latex/sqrt.tex ) uses triple-double floating point maths for its final iteration.
The sqrt function at https://gitlab.com/gtd-gmbh/libmcs/-/blob/core-math/libm/mathd/sqrtd.c?ref_type=heads  (the SUN Microsystems implementation) also uses integer maths, but
 it generates the sqrt one bit at a time rather than MUSL which uses an iterative approach with a quadratic convergence rate (from an initial reasonably good estimate).

"accuracy.pdf" also provides the definition of "ulp" that is used here as a measure of accuracy, with 0.500 ulp meaning the results are as accurate as possible and larger values meaning lower accuracy.

cr_sqrt.c is from musl-1.2.6  - which is MIT licensed see below.

Other functions are from the CORE-MATH project ( https://core-math.gitlabpages.inria.fr/ ) which is also under the MIT license.

The CORE-MATH projects mission is to provide on-the-shelf high performance open-source mathematical functions with correct rounding that can be integrated into current mathematical libraries.
This library does not provide a binary64 sqrt function, but does provide sin,cos,atan2,exp,log.
Note that I actually used the code from core-math-binary64, a cross-platform library with binary64 methods from CORE-MATH designed by Bob Burger at https://github.com/indigobio/core-math-binary64 .
This says "The implementation is a snapshot of the src/binary64 directory of the CORE-MATH project at commit 4acb80ff2b6fc614023a1ada2110bf7625d3cf8c from 2026-07-29, with the build and testing infrastructure removed. The function sources are otherwise unchanged."

The author of this library (ya_crmath_lib, Peter Miller) added some portability improvements to these functions (in particular for Intel-32 bit targets), and for the sqrt function provided a significant speedup for compilers that provide 128 bit integer variables.
He also provided the test code.
 
Overall its believed these functions will operate correctly on any processor/C compiler which implements IEEE 754 format double precision floating point (binary64).

Version 1v1
===========
This version adds a number of functions (ya_exp(), ya_log(), ya_sin(), ya_cos()) that give identical results to the cr_xxx functions of release 1v0, but execute faster. 
Also, as "accuracy.pdf" shows the built in square root for most compilers is accurate to 0.500ulp, this is used by default for sqrt as most processors have a built in square root instruction (Intel processors have this with SSE2) which is significantly faster than the software implementation in cr_sqrt.c (1.45 secs vs 5.86 secs with an Intel i3-10100).

It should be noted that these new functions have only been tested with the default rounding mode ("round to even").

The new functions are a hybrid of the CORE-MATH algorithms and the LLVM algorithms in the libc maths library from  llvm-project-23.1.1  (\llvm-project-23.1.1\libc\src\__support\FPUtil & \llvm-project-23.1.1\libc\src\__support\math https://libc.llvm.org/headers/math/index.html#higher-math-functions ), and only exist because of the excellent work of both the CORE-MATH and LLVM teams.

These new functions also leverage the double-double basic maths functions in https://github.com/p-j-miller/ya_double_double

Execution times (secs) within the supplied test program for the complete set of tests for each function are (using i3-10100 processor, with 1v0 being the cr_xxx time and 1v1 being the ya_xxx time):
|function| 1v0 | 1v1 |
|--------|-----|-----|
|sqrt	 | 5.86| 1.45|
|log	 | 9.03| 4.03|
|exp	 | 7.37| 4.37|
|sin	 |28.71| 8.90|
|cos	 |26.51| 8.85|

When 1v1 is used within wmawk2 ( https://github.com/p-j-miller/wmawk2 ) a single call to sqrt() in an awk script (using i3-10100 processor) takes ~ 3ns, while a call to exp() takes ~ 5ns and sin()/cos() take ~ 15ns. 

1v1 also adds a number of derived functions that leverage the built in functions to give reproducible results, but have errors > 0.500 ulp. These are:
~~~
 double ya_tan(double x);
 double ya_asin(double x);
 double ya_acos(double x)
~~~
The errors of these are <=1.31ulp (for comparison the same functions in the UCRT are <=1.32ulp)

The test program also tests two recurrences, one that provides sin(n*x) and cos(n*x) and the other just sin(n*x), with n increasing by one each step. These are significantly faster than using ya_sin/ya_cos but while they give reproducible results, they have errors significantly more than 0.500 ulp for very small values of x when n*x is very close to an integer multiple of pi (i.e. when sin(x)=0). The absolute error in this case is still small.

USE
===
At the start of files using these functions add (you may need to edit the relative path depending on your layout of directories):
~~~
#define YA_CRMATH_LIB_REPLACE  /* if defined cr_xxx/ya_xxx routines will replace standard C functions */
#include "../ya_crmath_lib/ya_crmath.h
~~~
If YA_CRMATH_LIB_REPLACE is not defined then these functions must be called as for example cr_cos(x), if it is defined then cos(x) will also use the cr_cos(x) function (which makes it simpler to add these functions into existing code).

As of 1v1, by default using YA_CRMATH_LIB_REPLACE will give the ya_xxx functions (as they give the same results as the cr_xxx functions, just faster). A "#if" in ya_crmath.h allows you to quickly revert to the cr_xxx functions if required, and both cr_xxx and ya_xxx functions can be called directly.

You will also need to compile and link to all the cr_*.c and ya_*.c files in the directory ya_crmath_lib (see the example compilation in the INSTALLATION section below).  
You must compile the cr_xxx functions even if you only use the ya_xxx functions, as the ya_xxx functions leverage capability within the cr_xxx.c files.

INSTALLATION 
============
Copy the files into a directory called ya_crmath_lib

A makefile and dev-c++ project file (ya_maths.dev) for the test program (main.c) is included.

A suitable gcc command line (for Windows or Linux) for the test program is:
~~~
gcc -O3 main.c cr_*.c ya_*.c ../my_printf/my_printf.c ../hr_timer/hr_timer.c -lm -lquadmath
~~~
This should give no warnings or errors, but will give 3 notes like:
~~~
cr_exp.c:71:11: note: '#pragma message: roundeven_finite(x) using C function'
   71 |   #pragma message( "roundeven_finite(x) using C function")
      |           ^~~~~~~
In file included from cr_pow.c:88:
cr_pow.h:76:11: note: '#pragma message: roundeven_finite(x) using C function'
   76 |   #pragma message( "roundeven_finite(x) using C function")
      |           ^~~~~~~
cr_sqrt.c:111:11: note: '#pragma message: mul64() using __int128'
  111 |   #pragma message( "mul64() using __int128")
      |           ^~~~~~~
~~~
These notes will not appear if NDEBUG is defined (-DNDEBUG as 1st argument to gcc command), they tell the user what options within the code are being used.

my_printf is available at https://github.com/p-j-miller/my_printf and is only needed for the test program (this is used as the Windows gcc port (mingw64) printf does not correctly print out long doubles). 

hr_timer is available at https://github.com/p-j-miller/hr_timer and again this is only needed for the test program.

ya_double_double.h is available at https://github.com/p-j-miller/ya_double_double

Note that the test program must be linked with -lquadmath (this is also not required in normal use of ya_crmath_lib).

The test program is run as "a" under windows, or "./a.out" under Linux. The tests typically take around 11 minutes to execute, with the tests for cr_cos and cr_sin the slowest.

The test programs checks each function with over 66 Million values.
Each function is checked against a higher accuracy version - either long double or __float128.

For example the expected output from the cr_sqrt() function tests is:
~~~
Testing cr_sqrt(x) vs __builtin_sqrtl(x) :
Checking values in special_vals[]...
searching for larger errors...
cr_sqrt(x) vs __builtin_sqrtl(x) max abs difference found was 7.44283e+137 (=54.000 bits) at 1.7976931348623157081e+308
 max relative error found was 1.11022e-16 (=53.000 bits) at 1.0000000000000144329
 So accuracy of approximation is 1.000 bits

1ulp at x=1.0000000000000144329, y=1.0000000000000071054 is 2.22045e-16
So accuracy of approximation vs __builtin_sqrtl(x) is 0.500 ulp
~~~
The tests are passed if the final line for each function is a line like that below showing an accuracy of 0.500 ulp :
~~~
So accuracy of approximation vs __builtin_sqrtl(x) is 0.500 ulp
~~~
The test program does require that long double is different from double (typically sizeof double=8 and sizeof long double = 12 or 16 - these are printed at the start of the test program), and that the compiler supports __float128's. The author mainly uses the winlibs gcc compiler on Windows ( https://winlibs.com/ ) and gcc on Linux and has tested this code using gcc versions between 13.3.0 and 16.2.0 .

MUSL license:
=============
musl as a whole is licensed under the following standard MIT license:

----------------------------------------------------------------------
Copyright © 2005-2020 Rich Felker, et al.

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

