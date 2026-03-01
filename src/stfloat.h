#pragma once
#include <math.h>

#ifdef STARTRACK_USE_FLOAT

#define stfloat float
#define stsqrt(x)    sqrtf(x)
#define stsin(x)     sinf(x)
#define stcos(x)     cosf(x)
#define sttan(x)     tanf(x)
#define stpow(x, n)  powf(x, n)
#define stacos(x)    acosf(x)
#define stasin(x)    asinf(x)
#define statan2(x,y) atan2f(x, y)
#define stexp(x)     expf(x)
#define stlog(x)     logf(x)
#define stfabs(x)    fabsf(x)
#define stceil(x)    ceilf(x)

#else

#define stfloat double
#define stsqrt(x)    sqrt(x)
#define stsin(x)     sin(x)
#define stcos(x)     cos(x)
#define sttan(x)     tan(x)
#define stpow(x, n)  pow(x, n)
#define stacos(x)    acos(x)
#define stasin(x)    asin(x)
#define statan2(x,y) atan2(x, y)
#define stexp(x)     exp(x)
#define stlog(x)     log(x)
#define stfabs(x)    fabs(x)
#define stceil(x)    ceil(x)

#endif

