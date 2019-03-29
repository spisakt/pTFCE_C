#include <iostream>
#include <cmath>

#include <functional>
#include <memory>
#include <utility>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

template < typename F >
class gsl_quad
{
  F f;
  void * p;
  int limit;
  std::unique_ptr < gsl_integration_workspace,
                    std::function < void(gsl_integration_workspace*) >
                    > workspace;

  /*static double gsl_wrapper(double x, void * p)
  {
    gsl_quad * t = reinterpret_cast<gsl_quad*>(p);
    return t->f(x);
  }*/

public:
  gsl_quad(F f, void * p, int limit)
    : f(f)
    , p(p)
    , limit(limit)
    , workspace(gsl_integration_workspace_alloc(limit), gsl_integration_workspace_free)
  {}

  double integrate(double min, double max, double epsabs, double epsrel)
  {
    gsl_function gsl_f;
    //gsl_f.function = &gsl_wrapper;
    gsl_f.function = this->f;
    gsl_f.params = this->p;

    double result, error;
    if ( !std::isinf(min) && !std::isinf(max) )
    {
      gsl_integration_qags ( &gsl_f, min, max,
                             epsabs, epsrel, limit,
                             workspace.get(), &result, &error );
    }
    else if ( std::isinf(min) && !std::isinf(max) )
    {
      gsl_integration_qagil( &gsl_f, max,
                             epsabs, epsrel, limit,
                             workspace.get(), &result, &error );
    }
    else if ( !std::isinf(min) && std::isinf(max) )
    {
      gsl_integration_qagiu( &gsl_f, min,
                             epsabs, epsrel, limit,
                             workspace.get(), &result, &error );
    }
    else
    {
      gsl_integration_qagi ( &gsl_f,
                             epsabs, epsrel, limit,
                             workspace.get(), &result, &error );
    }

    //std::cout << "i = " << result << "; a.err = " << error << std::endl;
    return result;
  }
};



template < typename F >
double quad(F func,
            void * params,
            std::pair<double,double> const& range,
            double epsabs = 0.01, double epsrel = 0.01,
            int limit = 50)
{
  return gsl_quad<F>(func, params, limit).integrate(range.first, range.second, epsabs, epsrel);
}

//Examples
//without parameters (p can be null pointer)
//result = quad([](double x, void*) { return x*x; }, p, {0,1});
//result = quad([](double x, void*) { return 1/(x*x); }, p, {1,INFINITY});
//result = quad([](double x, void*){return std::exp(-x*x);}, p, {-INFINITY, INFINITY});
//with parameters (dcl(F x, void * p))
//result = quad( &dcl, p, {thr, INFINITY}, 10e-7, 0.0, GSL_INTEGRATION_INTERVAL_LIMIT);
