# ya_crmath_lib
ya_crmath_lib is a portable C maths library that should always give the same results independent of the C compiler / operating system.
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
on the basis that cpu's tend to be faster on integer code.  In particular the sqrt function in crlibm (see e.g.  https://github.com/taschini/crlibm/blob/master/docs/latex/sqrt.tex which uses triple-double fp. maths for the final iteration).
The sqrt function at https://gitlab.com/gtd-gmbh/libmcs/-/blob/core-math/libm/mathd/sqrtd.c?ref_type=heads  (the SUN Microsystems implementation) also uses integer maths, but
 it generates the sqrt one bit at a time rather than MUSL which uses an iterative approach with a quadratic convergence rate (from an initial reasonably good estimate).


cr_sqrt.c is from musl-1.2.6  - which is MIT licensed see below.

Other functions are from the CORE-MATH project ( https://core-math.gitlabpages.inria.fr/ ) which is also under the MIT license.
The CORE-MATH projects mission is to provide on-the-shelf high performance open-source mathematical functions with correct rounding that can be integrated into current mathematical libraries.
This library does not provide a binary64 sqrt function, but does provide sin,cos,atan2,exp,log.
Note that I actually used the code from core-math-binary64, a cross-platform library with binary64 methods from CORE-MATH designed by Bob Burger at https://github.com/indigobio/core-math-binary64 .
This says "The implementation is a snapshot of the src/binary64 directory of the CORE-MATH project at commit 4acb80ff2b6fc614023a1ada2110bf7625d3cf8c from 2026-07-29, with the build and testing infrastructure removed. The function sources are otherwise unchanged."

The author of this library (ya_crmath_lib, Peter Miller) added some portability improvements to these functions (in particular for Intel-32 bit targets), and for the sqrt function provided a significant speedup for compilers that provide 128 bit integer variables.
He also provided the test code.
 
Overall its believed these functions will operate correctly on any processor/C compiler which implements IEEE 754 format double precision floating point (binary64).

MSUL license:
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

