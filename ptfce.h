#ifndef PTFCE_H
#define PTFCE_H

#include "newimage/newimageall.h"
#include "newmatio.h"
#include "ndf.h"

using namespace NEWIMAGE;

enum RPVEstimationMode {NORPV = 0, FSLRPV, SPMRPV};

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
    volume<T> RPV;              //Resels per Voxel image
    volume<T> FWHMimg;          //FWHM image
    bool adjustClusterSize;     //RFT cluster-size adjustment based on RPV
    unsigned int autosmooth;	//Smoothness estimation mode (0-manual; 1-Z image; 2-residual image) TODO proper mode selection
    int RPVMode;		//Local smoothness implementation mode (1-FSL*; 2-SPM)
    double ZestThr;		//Cluster-forming Z threshold below which P(h|c) is estimated as P(h), due to limitations of GRF theory (default: 1.3)
    bool _verbose;		//Print progress bar and diagnostic messages

    double aggregateLogp(NEWMAT::ColumnVector &pvals);

  public:
    pTFCE(volume<T> img, volume<T> mask) : img(img), mask(mask)
    {
	Nh = 100;
	dof = 0.0;
	logpmin = 0.0;
	logpmax = -1 * pnormR(img.max(mask), 0.0, 1.0, false, true);
	dh = (logpmax - logpmin) / (Nh-1);
	ZestThr = 1.3;
	Rd = 0.0;
	resels = 0.0;
	autosmooth = 1;
	RPVMode = 1;
	adjustClusterSize = false;
	_verbose = false;
    }

    pTFCE(volume<T> img, volume<T> mask, volume4D<T> residual, double dof): img(img), mask(mask), residual(residual), dof(dof)
    {
	Nh = 100;
	logpmin = 0.0;
	logpmax = -1 * pnormR(img.max(mask), 0.0, 1.0, false, true);
	dh = (logpmax - logpmin) / (Nh-1);
	ZestThr = 1.3;
	Rd = 0.0;
	resels = 0.0;
	autosmooth = 2;
	RPVMode = 1;
	adjustClusterSize = false;
	_verbose = false;
    }

    ~pTFCE();
    void destroy();

    volume<T> pTFCEimg;
    volume<T> logp_pTFCE;
    volume<T> Z_pTFCE;
    double number_of_resels;
    double fwer005_Z;

    void estimateSmoothness();
    int calculate();

    void setThresholdCount(int N);
    void setSmoothness(double Rd, unsigned long V, double resels);
    void setRPVEstimationMode(int mode);
    int  getRPVEstimationMode();
    void printSmoothness();
    void loadRPV(const string& filename);
    void saveRPV(const string& filename);
    void saveFWHM(const string& filename);
    void setRFTAdjust(bool a);
    int  getRFTAdjust();
    void setZestThr(double Z);
    void quiet();
    void verbose();
};



#endif //PTFCE_H
