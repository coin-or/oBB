/* Copyright (c) 2007-2014 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#include <math.h>
#include <float.h>

#include "nlopt-internal.h"

/*************************************************************************/

int nlopt_isinf(double x) {
     return fabs(x) >= HUGE_VAL * 0.99
#ifdef HAVE_ISINF
	  || isinf(x)
#endif
	  ;
}

int nlopt_isfinite(double x) {
    return fabs(x) <= DBL_MAX;
}

int nlopt_istiny(double x)
{
#if defined(HAVE_FPCLASSIFY)
    return x == 0.0 || fpclassify(x) == FP_SUBNORMAL;
#elif defined(_WIN32)
    if (x == 0.0)
        return 1;
    else {
        int c = _fpclass(x);
        return c == _FPCLASS_ND || c == _FPCLASS_PD;
    }
#else
    return fabs(x) < 2.2250738585072014e-308; /* assume IEEE 754 double */
#endif
}

int nlopt_isnan(double x)
{
#if defined(HAVE_ISNAN)
    return isnan(x);
#elif defined(_WIN32)
    return _isnan(x);
#else
    return x != x; /* might fail with aggressive optimization */
#endif
}

/*************************************************************************/

void NLOPT_STDCALL nlopt_version(int *major, int *minor, int *bugfix)
{
     *major = 2;
     *minor = 5;
     *bugfix = 0;
}

/*************************************************************************/

static const char nlopt_algorithm_names[NLOPT_NUM_ALGORITHMS][256] = {
     "Sequential Quadratic Programming (SQP) (local, derivative)"
};

const char * NLOPT_STDCALL nlopt_algorithm_name(nlopt_algorithm a)
{
     if (((int) a) < 0 || a >= NLOPT_NUM_ALGORITHMS) return "UNKNOWN";
     return nlopt_algorithm_names[a];
}
