#ifndef NDF_H
#define NDF_H

#include <cmath>
#include <cfloat>
#include <limits>

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
        density = - pow(x - mean, 2.0) / 2.0*sd*sd - log(sqrt(2.0*M_PI*sd*sd));
    return density;
}

template <typename T>
inline double pnorm(T x, double mean = 0.0, double sd = 1.0, bool _lowtail = true, bool _logp = false)
{
    long double distrib = (1 + erfl( (long double)(x - mean) / (sd * sqrt(2.0)) )) / 2.0;
    if (!_lowtail)
	distrib = 1.0 - distrib;
    if (_logp)
    {
        std::cout << log(1- (1 + erfl( (long double)(4 - mean) / (sd * sqrt(2.0)) )) / 2.0 ) << "\n";
        std::cout << log(1- 0.5 * erfcl(-(4-mean)/(sd*sqrt(2))) ) << "\n";

        std::cout << log(1- (1 + erfl( (long double)(5 - mean) / (sd * sqrt(2.0)) )) / 2.0 ) << "\n";
        std::cout << log(1- 0.5 * erfcl(-(5-mean)/(sd*sqrt(2))) ) << "\n";

        std::cout << log(1- (1 + erfl( (long double)(10 - mean) / (sd * sqrt(2.0)) )) / 2.0 ) << "\n";
        std::cout << log(1- 0.5 * erfcl(-(10-mean)/(sd*sqrt(2))) ) << "\n";


        std::cout << "***" << "\n";


        return log( distrib );
    }
    return distrib;
}

template <typename T>
inline double qnorm(T p, double mean = 0.0, double sd = 1.0, bool _lowtail = true, bool _logp = false)
{
    double quantile;
    if (_logp)
	p = exp(p);
    quantile = sd * sqrt(2.0) * erfinv2(2.0*p - 1);
    if (_lowtail)
	quantile = mean + quantile;
    else
	quantile = mean - quantile;
    return quantile;
}

inline void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p)
{
/* i_tail in {0,1,2} means: "lower", "upper", or "both" :
   if(lower) return  *cum := P[X <= x]
   if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
*/
    const static double a[5] = {
	2.2352520354606839287,
	161.02823106855587881,
	1067.6894854603709582,
	18154.981253343561249,
	0.065682337918207449113
    };
    const static double b[4] = {
	47.20258190468824187,
	976.09855173777669322,
	10260.932208618978205,
	45507.789335026729956
    };
    const static double c[9] = {
	0.39894151208813466764,
	8.8831497943883759412,
	93.506656132177855979,
	597.27027639480026226,
	2494.5375852903726711,
	6848.1904505362823326,
	11602.651437647350124,
	9842.7148383839780218,
	1.0765576773720192317e-8
    };
    const static double d[8] = {
	22.266688044328115691,
	235.38790178262499861,
	1519.377599407554805,
	6485.558298266760755,
	18615.571640885098091,
	34900.952721145977266,
	38912.003286093271411,
	19685.429676859990727
    };
    const static double p[6] = {
	0.21589853405795699,
	0.1274011611602473639,
	0.022235277870649807,
	0.001421619193227893466,
	2.9112874951168792e-5,
	0.02307344176494017303
    };
    const static double q[5] = {
	1.28426009614491121,
	0.468238212480865118,
	0.0659881378689285515,
	0.00378239633202758244,
	7.29751555083966205e-5
    };

    double xden, xnum, temp, del, eps, xsq, y;
    int i, lower, upper;
    /* Consider changing these : */
    eps = std::numeric_limits<double>::epsilon() * 0.5;

    /* i_tail in {0,1,2} =^= {lower, upper, both} */
    lower = i_tail != 1;
    upper = i_tail != 0;

    y = fabs(x);
    if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
	if (y > eps) {
	    xsq = x * x;
	    xnum = a[4] * xsq;
	    xden = xsq;
	    for (i = 0; i < 3; ++i) {
		xnum = (xnum + a[i]) * xsq;
		xden = (xden + b[i]) * xsq;
	    }
	} else xnum = xden = 0.0;

	temp = x * (xnum + a[3]) / (xden + b[3]);
	if(lower)  *cum = 0.5 + temp;
	if(upper) *ccum = 0.5 - temp;
	if(log_p) {
	    if(lower)  *cum = log(*cum);
	    if(upper) *ccum = log(*ccum);
	}
    }
    else if (y <= 5.656854249492380195206754896838) { //M_SQRT_32

	/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

	xnum = c[8] * y;
	xden = y;
	for (i = 0; i < 7; ++i) {
	    xnum = (xnum + c[i]) * y;
	    xden = (xden + d[i]) * y;
	}
	temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)							\
	xsq = trunc(X * 16) / 16;				\
	del = (X - xsq) * (X + xsq);					\
	if(log_p) {							\
	    *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);	\
	    if((lower && x > 0.) || (upper && x <= 0.))			\
		  *ccum = log1p(-exp(-xsq * xsq * 0.5) *		\
				exp(-del * 0.5) * temp);		\
	}								\
	else {								\
	    *cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;	\
	    *ccum = 1.0 - *cum;						\
	}

#define swap_tail						\
	if (x > 0.) {/* swap  ccum <--> cum */			\
	    temp = *cum; if(lower) *cum = *ccum; *ccum = temp;	\
	}

	do_del(y);
	swap_tail;
    }

/* else	  |x| > sqrt(32) = 5.657 :
 * the next two case differentiations were really for lower=T, log=F
 * Particularly	 *not*	for  log_p !
 * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
 *
 * Note that we do want symmetry(0), lower/upper -> hence use y
 */
    else if((log_p && y < 1e170) /* avoid underflow below */
	/*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
	 * Then, make use of  Abramowitz & Stegun, 26.2.13, something like
	 xsq = x*x;
	 if(xsq * DBL_EPSILON < 1.)
	    del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
	 else
	    del = 0.;
	 *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
	 *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./
 	 swap_tail;
	 [Yes, but xsq might be infinite.]
	*/
	    || (lower && -37.5193 < x  &&  x < 8.2924)
	    || (upper && -8.2924  < x  &&  x < 37.5193)
	) {

	/* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
	xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
	xnum = p[5] * xsq;
	xden = xsq;
	for (i = 0; i < 4; ++i) {
	    xnum = (xnum + p[i]) * xsq;
	    xden = (xden + q[i]) * xsq;
	}
	temp = xsq * (xnum + p[4]) / (xden + q[4]);
	temp = (0.398942280401432677939946059934 - temp) / y; // M_1_SQRT_2PI

	do_del(x);
	swap_tail;
    } else { /* large x such that probs are 0 or 1 */
	if(x > 0) {	*cum = (log_p ? 0. : 1.); *ccum = (log_p ? -INFINITY : 0.);	}
	else {	        *cum = (log_p ? -INFINITY : 0.); *ccum = (log_p ? 0. : 1.);	}
    }
    return;
}

inline double pnormR(double x, double mu=0.0, double sigma=1.0, int lower_tail=false, int log_p=false)
{
    double p, cp;

    /* Note: The structure of these checks has been carefully thought through.
     * For example, if x == mu and sigma == 0, we get the correct answer 1.
     */
    p = (x - mu) / sigma;
    x = p;
    pnorm_both(x, &p, &cp, (lower_tail ? 0 : 1), log_p);
    return(lower_tail ? p : cp);
}

inline double qnormR(double p, double mu=0.0, double sigma=1.0, int lower_tail=false, int log_p=false)
{
    double p_, q, r, val;
    p_ = (log_p ? (lower_tail ? exp(p) : - std::expm1(p)) : (lower_tail ? (p) : (0.5 - (p) + 0.5)));/* real lower_tail prob. p */
    q = p_ - 0.5;

    if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
	val =
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary */

	/* r = min(p, 1-p) < 0.075 */
	if (q > 0)
	    r = (log_p ? (lower_tail ? -std::expm1(p) : exp(p)) : (lower_tail ? (0.5 - (p) + 0.5) : (p)));/* 1-p */
	else
	    r = p_;/* = R_DT_Iv(p) ^=  p */

	r = sqrt(- ((log_p &&
		     ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
		    p : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
        }
        else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
        }

	if(q < 0.0)
	    val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return mu + sigma * val;
}

#endif //NDF_H
