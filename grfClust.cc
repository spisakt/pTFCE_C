#include <iostream>
#include <cmath>
#include "gsl_integration_wrapper.h"
#include "ndf.h"
#include "mathutil.h"

#include "grfClust.h"

#if !defined(M_PI)
#define M_PI (4 * atan(1.0))
#endif

#define GSL_INTEGRATION_INTERVAL_LIMIT 8
#define GSL_INTEGRATION_UPPER std::numeric_limits<double>::infinity()
//#define GSL_INTEGRATION_UPPER 6.0
#define GSL_EPS_ABS 0.1
#define GSL_EPS_REL 0.1

double dvox(double h)
{
    return dnorm(h);
}


double pvox(double h)
{
    return pnormR(h, 0.0, 1.0, false, false);
}


double Es(double h, double V, double Rd)
{
    double ret = log(V) + pnormR(h, 0.0, 1.0, false, true);
    double h2 = h*h;

    if (h >= 1.1)
        ret = ret - ( log(Rd) + log(h2-1) - h2/2 - 2*log(2*M_PI) );

if (isnan2(exp(ret))) { std::cout << "nan:ES" << std::endl; exit; }

    return exp(ret);
}


double dcl(double h, void * p)
{
    struct dcl_params* params = (struct dcl_params*) p;
    double V = params->V;
    double Rd      = params->Rd;
    double c       = params->c;
    double ZestThr = params->ZestThr;

    if (h < ZestThr)
        return 0.0;

    double lambda = pow( (Es(h, V, Rd) / tgamma(2.5)), -2.0/3.0 );
    double dclust = lambda * exp( -lambda * pow(c, 2.0/3.0) );

if (isnan2(dclust)) { std::cout << "nan:dcl" << std::endl; exit; }

    return dclust;
}


double dclust( double h, void * p )
{
    struct dcl_params* params = (struct dcl_params*) p;
    double ZestThr = params->ZestThr;

    //double result;

    //result = quad( &dcl, p, {ZestThr, GSL_INTEGRATION_UPPER}, GSL_EPS_ABS, GSL_EPS_REL, GSL_INTEGRATION_INTERVAL_LIMIT);

    return dcl(h, params); // / result;
}


double dvox_dclust( double h, void * p )
{

if (isnan2(dvox(h) * dclust(h, p))) { std::cout << "nan:dvox_dclust" << std::endl; exit; }

    return dvox(h) * dclust(h, p);
}


double dvox_clust(double h, void * p)
{
    struct dcl_params* params = (struct dcl_params*) p;
    double ZestThr = params->ZestThr;

    double result;

    if ( std::isnan(dclust(h, params)) )
        return dvox(h);

    result = quad( &dvox_dclust, p, {ZestThr, GSL_INTEGRATION_UPPER}, GSL_EPS_ABS, GSL_EPS_REL, GSL_INTEGRATION_INTERVAL_LIMIT);

if (isnan2(result) || result == 0.0) { std::cout << "nan:dvox_clust " << result << std::endl; exit;}

    return dvox_dclust(h, params) / result;
}


double pvox_clust(double actH, void * p)
{
    struct dcl_params* params = (struct dcl_params*) p;
    double ZestThr = params->ZestThr;

    double result;

    if ( actH < ZestThr )  // ZestThr: GRF theory might not apply at low thresholds
        return pvox(actH);
    if ( std::isnan(dvox_clust(actH, params)) )
        return exp(-745);

    result = quad( &dvox_clust, p, {actH, GSL_INTEGRATION_UPPER}, GSL_EPS_ABS, GSL_EPS_REL, GSL_INTEGRATION_INTERVAL_LIMIT);

if (isnan2(result) || result == 0.0) { std::cout << "nan:pvox_clust " << result << std::endl; exit;}

    return result;
}
