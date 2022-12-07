#include <iostream>
#include <cmath>
#include "libprob/libprob.h" //FSL <6.0.5
//#include "cprob/libprob.h" //FSL >=6.0.5
#include "pvalutil.h"

#define MINP 1e-15

using namespace MISCMATHS;

double fwerp2z(double nresels, double p, bool twotailed, bool grf)
{
  double result;

  if (twotailed)
    p/=2;

  if (grf)
    p/=nresels*0.11694; /* (4ln2)^1.5 / (2pi)^2 */

  if (p<MINP) p=MINP;

  if (p>0.5)
    result = 0.0;
  else
    {
      if (grf)
	{
	  double l=2, u=100, z=0, pp;
	  while (u-l>0.001)
	    {
	      z=(u+l)/2;
	      pp=exp(-0.5*z*z)*(z*z-1);
	      if (pp<p) u=z;
	      else      l=z;
	    }
	  result = z;
	}
      else
	{
	  result = ndtri(1-p);
	}
    }

  if (false) std::cout << "fwerp2z(" << nresels << ", " << p << "): " << result << std::endl;
  return(result);
}


double fwerz2p(double nresels, double z, bool twotailed, bool grf)
{
  double p;

  if (twotailed)
    z = fabs(z);

  if (grf)
    {
      if (z<2)
        p = 1; /* Below z of 2 E(EC) becomes non-monotonic */
      else
        p = nresels * 0.11694 * exp(-0.5*z*z) * (z*z-1);  /* 0.11694 = (4ln2)^1.5 / (2pi)^2 */
    }
  else
    {
      p = 1 - ndtr(z);
    }

  if (twotailed)
    p*=2;

  p=std::min(p,1.0);

  if (false) std::cout << "fwerz2p(" << nresels << ", " << z << "): " << p << std::endl;
  return(p);
}


//g++ -o pvalutil.o -c pvalutil.cc -I$FSLDIR/extras/include -std=c++11 -lm
