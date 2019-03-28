#ifndef PTFCE_H
#define PTFCE_H

#include "newimage/newimageall.h"
#include "newmatio.h"
#include "ndf.h"

using namespace NEWIMAGE;

template<class T>
class pTFCE
{
  private:
    volume<T> img;		//Nitfi Z-score image to enhance
    volume<T> mask;		//Mask
    int Nh;			//Number of thresholds
    double logpmin;		//min threshold
    double logpmax;		//max threshold
    double dh;			//difference between thresholds
    volume4D<T> residual;	//4D residual data for better estimate of image smoothness
    double dof;			//Degrees of freedom (optional, but obligatory if residual is specified)
    double Rd;			//Resel count
    unsigned long V;		//Number of voxels in mask
    double resels;		//
    double dLh;			//
    unsigned int autosmooth;	//Smoothness estimation mode (0-manual; 1-Z image; 2-residual image)
    double ZestThr;		//Cluster-forming Z threshold below which P(h|c) is estimated as P(h), due to limitations of GRF theory (default: 1.3)
    bool _verbose;		//Print progress bar and diagnostic messages

    void estimateSmoothness();
    double aggregateLogp(NEWMAT::ColumnVector &pvals);

  public:
    pTFCE(volume<T> img, volume<T> mask) : img(img), mask(mask)
    {
	Nh = 100;
	dof = 0.0;
	logpmin = 0.0;
	logpmax = -1 * pnorm(img.max(mask), 0.0, 1.0, false, true);
	dh = (logpmax - logpmin) / (Nh-1);
	ZestThr = 1.3;
	resels = 0.0;
	autosmooth = 1;
	_verbose = false;
	estimateSmoothness();
    }

    pTFCE(volume<T> img, volume<T> mask, volume4D<T> residual, double dof): img(img), mask(mask), residual(residual), dof(dof)
    {
	Nh = 100;
	logpmin = 0.0;
	logpmax = -1 * pnorm(img.max(mask), 0.0, 1.0, false, true);
	dh = (logpmax - logpmin) / (Nh-1);
	ZestThr = 1.3;
	resels = 0.0;
	autosmooth = 2;
	_verbose = false;
	estimateSmoothness();
    }

    ~pTFCE();
    void destroy();

    volume<T> pTFCEimg;
    volume<T> logp_pTFCE;
    volume<T> Z_pTFCE;
    double number_of_resels;
    double fwer005_Z;

    int calculate();

    void setThresholdCount(int N);
    void setSmoothness(double Rd, unsigned long V, double resels);
    void printSmoothness();
    void setZestThr(double Z);
    void quiet();
    void verbose();
};



#endif //PTFCE_H
