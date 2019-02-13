#ifndef NDF_H
#define NDF_H

#include <cmath>

#if !defined(M_PI)
#define M_PI (4 * atan(1.0))
#endif


//A handy approximation for the error function and its inverse... by Sergei Winitzki.
template <typename T>
inline double erfinv2(T x)
{
   double tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0 : 1.0;

   x = (1 - x)*(1 + x);        // x = 1 - x*x;
   lnx = log(x);

   tt1 = 2/(M_PI*0.147) + 0.5 * lnx;
   tt2 = 1/(0.147) * lnx;

   return(sgn*sqrtf(-tt1 + sqrt(tt1*tt1 - tt2)));
}

/* compute natural logarithm with a maximum error of 0.85089 ulp */
template <typename T>
inline float logf1(T a)
{
    float i, m, r, s, t;
    int e;

    m = frexpf (a, &e);
    if (m < 0.666666667f) { // 0x1.555556p-1
        m = m + m;
        e = e - 1;
    }
    i = (float)e;
    /* m in [2/3, 4/3] */
    m = m - 1.0f;
    s = m * m;
    /* Compute log1p(m) for m in [-1/3, 1/3] */
    r =             -0.130310059f;  // -0x1.0ae000p-3
    t =              0.140869141f;  //  0x1.208000p-3
    r = fmaf (r, s, -0.121484190f); // -0x1.f19968p-4
    t = fmaf (t, s,  0.139814854f); //  0x1.1e5740p-3
    r = fmaf (r, s, -0.166846052f); // -0x1.55b362p-3
    t = fmaf (t, s,  0.200120345f); //  0x1.99d8b2p-3
    r = fmaf (r, s, -0.249996200f); // -0x1.fffe02p-3
    r = fmaf (t, m, r);
    r = fmaf (r, m,  0.333331972f); //  0x1.5554fap-2
    r = fmaf (r, m, -0.500000000f); // -0x1.000000p-1
    r = fmaf (r, s, m);
    r = fmaf (i,  0.693147182f, r); //  0x1.62e430p-1 // log(2)
    if (!((a > 0.0f) && (a <= 3.40282346e+38f))) { // 0x1.fffffep+127
        r = a + a;  // silence NaNs if necessary
        if (a  < 0.0f) r = ( 0.0f / 0.0f); //  NaN
        if (a == 0.0f) r = (-1.0f / 0.0f); // -Inf
    }
    return r;
}

/* compute inverse error functions with maximum error of 2.35793 ulp */
template <typename T>
inline double erfinv1(T a)
{
    float p, r, t;
    t = fmaf (a, 0.0 - a, 1.0);
    t = logf1 (t);
    if (fabsf(t) > 6.125) { // maximum ulp error = 2.35793
        p =              3.03697567e-10f; //  0x1.4deb44p-32 
        p = fmaf (p, t,  2.93243101e-8f); //  0x1.f7c9aep-26 
        p = fmaf (p, t,  1.22150334e-6f); //  0x1.47e512p-20 
        p = fmaf (p, t,  2.84108955e-5f); //  0x1.dca7dep-16 
        p = fmaf (p, t,  3.93552968e-4f); //  0x1.9cab92p-12 
        p = fmaf (p, t,  3.02698812e-3f); //  0x1.8cc0dep-9 
        p = fmaf (p, t,  4.83185798e-3f); //  0x1.3ca920p-8 
        p = fmaf (p, t, -2.64646143e-1f); // -0x1.0eff66p-2 
        p = fmaf (p, t,  8.40016484e-1f); //  0x1.ae16a4p-1 
    } else { // maximum ulp error = 2.35456
        p =              5.43877832e-9f;  //  0x1.75c000p-28 
        p = fmaf (p, t,  1.43286059e-7f); //  0x1.33b458p-23 
        p = fmaf (p, t,  1.22775396e-6f); //  0x1.49929cp-20 
        p = fmaf (p, t,  1.12962631e-7f); //  0x1.e52bbap-24 
        p = fmaf (p, t, -5.61531961e-5f); // -0x1.d70c12p-15 
        p = fmaf (p, t, -1.47697705e-4f); // -0x1.35be9ap-13 
        p = fmaf (p, t,  2.31468701e-3f); //  0x1.2f6402p-9 
        p = fmaf (p, t,  1.15392562e-2f); //  0x1.7a1e4cp-7 
        p = fmaf (p, t, -2.32015476e-1f); // -0x1.db2aeep-3 
        p = fmaf (p, t,  8.86226892e-1f); //  0x1.c5bf88p-1 
    }
    r = a * p;
    return (double)r;
}


template <typename T>
inline double dnorm(T x, double mean = 0.0, double sd = 1.0, bool _log = false)
{
    double density = exp( - pow(x - mean, 2.0) / 2.0*sd*sd ) / sqrt(2.0*M_PI*sd*sd);
    if (_log)
	return log( density );
    return density;
}

template <typename T>
inline double pnorm(T x, double mean = 0.0, double sd = 1.0, bool _lowtail = true, bool _logp = false)
{
    double distrib = (1 + erf( (x - mean) / (sd * sqrt(2.0)) )) / 2.0;
    if (!_lowtail)
	distrib = 1.0 - distrib;
    if (_logp)
	return log( distrib );
    return distrib;
}

template <typename T>
inline double qnorm(T p, double mean = 0.0, double sd = 1.0, bool _lowtail = true, bool _logp = false)
{
    double quantile;
    if (_logp)
	p = exp(p);
    quantile = sd * sqrt(2.0) * erfinv1(2.0*p - 1);
    if (_lowtail)
	quantile = mean + quantile;
    else
	quantile = mean - quantile;
    return quantile;
}



#endif //NDF_H
