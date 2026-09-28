/* general settings (#defines) for ya_crmath_lib functions

  Created by Peter Miller, 11/9/2026
*/
  
 // #define CORE_MATH_SUPPORT_ERRNO /* if defined then support errno in core-maths functions - note cr_sqrt does not support this at present */
			// also note that supporting errno makes a significant difference to the compiled speed with sqty - so the use of the gcc arg -fmo-math-errno or #pragma GCC optimize ("O3,no-math-errno") is recommended unless errno is really required
			
  #define SQRT_FENV_SUPPORT 1 /*  For cr_sqrt() only : 1 means support to handle rounding modes and inexact exception, 0 means only basic rounding [ which is faster] */
  