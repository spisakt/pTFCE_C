#include <iostream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include "ptfce.h"
#include "newimage/newimageall.h" //volume, connected_components
#include "smoothest_ext.h" //smoothest
#include "grfClust.h" //pvox_clust
#include "pvalutil.h" //fwerp2z
#include "clustervalue.h" //projectClusterValues
#include "newmatio.h" //ColumnVector

using namespace NEWIMAGE;

template <class T>
void pTFCE<T>::destroy() {}

template <class T>
pTFCE<T>::~pTFCE()
{
    this->destroy();
}


template <class T>
void pTFCE<T>::estimateSmoothness()
{
    double FWHM[3];
    double FWHMmm[3];
    double sigmasq[3];

    if (_verbose) std::cout << "1. Smoothness estimation" << std::endl;
    if (_verbose) std::cout << "********************************" << std::endl;

    if (_verbose) std::cout << "RPV estimation mode = " << RPVMode << std::endl;

    switch (autosmooth)
    {
	case 0:
	{
	    if (_verbose) std::cout << "Using specified smoothness values" << std::endl;
	    break;
	}
	case 1:
	{
	    if (_verbose) std::cout << "Computing smoothness from Z-score image" << std::endl;
	    if (RPVMode == FSLRPV)
	    {
		smoothestVox(dLh, V, resels, FWHM, FWHMmm, sigmasq, RPV, FWHMimg,
	                  img, mask,
	                  100.0, _verbose);
	    }
	    if (RPVMode == SPMRPV)
	    {
		smoothest(dLh, V, resels, FWHM, FWHMmm, sigmasq,
	                  img, mask,
	                  100.0, _verbose);
		estimateRPV(RPV, FWHMimg,
	                  img, mask,
	                  100.0, _verbose);
	    }
	    if (RPVMode == NORPV)
	    {
		smoothest(dLh, V, resels, FWHM, FWHMmm, sigmasq,
	                  img, mask,
	                  100.0, _verbose);
	    }
	    Rd = dLh * V;
	    break;
	}
	case 2:
	{
	    if (_verbose) std::cout << "Computing smoothness from 4D residual image" << std::endl;
	    if (RPVMode == FSLRPV)
	    {
		smoothestVox(dLh, V, resels, FWHM, FWHMmm, sigmasq, RPV, FWHMimg,
	                  residual, mask,
	                  dof, _verbose);
	    }
	    if (RPVMode == SPMRPV)
	    {
		smoothest(dLh, V, resels, FWHM, FWHMmm, sigmasq,
	                  residual, mask,
	                  dof, _verbose);
		estimateRPV(RPV, FWHMimg,
	                  residual, mask,
	                  dof, _verbose);
	    }
	    if (RPVMode == NORPV)
	    {
		smoothest(dLh, V, resels, FWHM, FWHMmm, sigmasq,
	                  residual, mask,
	                  dof, _verbose);
	    }
	    Rd = dLh * V;
	    break;
	}
	default: break;
    }
}


template <class T>
double pTFCE<T>::aggregateLogp(NEWMAT::ColumnVector &pvals, bool printpvals)
{
    double nlogmin = -1.0 * log(std::numeric_limits<double>::min());
    double aggr=0.0;
    double tmp=0.0;

    for (int i = 1; i <= pvals.n_rows; ++i)
    {
	tmp = -1.0 * log(pvals(i));
	if (isinf(tmp) || isnan(tmp)) tmp = nlogmin;
	aggr+=tmp;
    }
    if (printpvals) cout << "diag: " << "aggr:" << " aggregated nlogp: " << aggr << endl;
    aggr = 0.5 * (sqrt(dh * (8.0 * aggr + dh)) - dh);

    return aggr;
}


template <class T>
int pTFCE<T>::calculate()
{
    // TODO check for NA/nan values
    // TODO check for (-)inf values

    /*volume4D<double> PVC;        //PVC.copyproperties(img);
    copybasicproperties(img, PVC); // now done during pTFCE class instantiation*/
    copyconvert(img, pTFCEimg);
    copyconvert(img, logp_pTFCE);
    copyconvert(img, Z_pTFCE);
    volume<T> thr;          copyconvert(img, thr);
    volume<int> ccc;        copyconvert(img, ccc);
    ////volume4D<int> CLUST;    CLUST.copyproperties(ccc);
    /////*volume4D<int> LABEL; //defined as member of pTFCE class*/
    ////LABEL.copyproperties(ccc);
    ////volumeClust<int> sizes;
    volumeClust<double> pvoxclust;

    img = mask_volume(img, mask);

    if (Rd == 0.0)
    {
	this->estimateSmoothness();
    }

    if (_verbose) std::cout << "2. Calculate pTFCE" << std::endl;
    if (_verbose) std::cout << "********************************" << std::endl;

    if (_verbose) std::cout << "* Performing pTFCE..." << std::endl;
    if (_verbose) std::cout << "../" << Nh << std::endl;

    for (unsigned int i = 0; i < Nh; ++i)
    {
	double p_thres = exp( -1 * ( logpmin + i * dh ) );
	double thres = qnormR(p_thres, 0.0, 1.0, false, false);
	NEWMAT::ColumnVector clustersizes;
	NEWMAT::ColumnVector clusterpvox;
	NEWMAT::ColumnVector clusterresels;

	if (_verbose && (i%10) == 0) std::cout << i+1 << " " << std::endl;
	if (logpmin == 0.0 && i == 0)
	    thres = std::numeric_limits<double>::lowest();
	    //thres = -9999999999;

	copyconvert(img, thr);
	thr.binarise(thres);
	ccc = connected_components(thr, clustersizes, 6);  //6, 18 or 26

	if (adjustClusterSize)
	{
	    volumeClust<T> RPC; copyconvert(ccc, RPC);
	    RPC.sumClusterRPVValues(RPV, clusterresels);
	}

	////copyconvert(ccc, sizes); //
	////sizes.projectClusterValues(clustersizes, 0.0f); //
	////CLUST.addvolume(sizes); //
	////if (_verbose)
	////{
	////    copyconvert(ccc, sizes);
	////    sizes.projectClusterValues(clustersizes, 0.0f);
	////    CLUST.addvolume(sizes);
	////    LABEL.addvolume(ccc);
	////}

	clusterpvox.ReSize(clustersizes.n_rows);
	for (int k = 1; k <= clusterpvox.n_rows; ++k)
	{
	    double pvc;
	    if (!adjustClusterSize)
	    {
		struct dcl_params params = {V, Rd, (double)clustersizes(k), ZestThr};
		pvc = pvox_clust(thres, &params);
	    }
	    else
	    {
		//clustersize given is resel must be converted to voxel units, hence the multiplications with "resels" (resel size)
		//mean of the RPV image * "resels" does not exactle equal 1 (it's ~ 0.993) TODO - check border voxels
		if (false && _verbose) cout << "*RESELS************************** " << resels << " threshold " << thres << " clustersizez(i) " << clustersizes(k) << " clusterresels(i) " << clusterresels(k) << " clusterresels(i) * resels " << clusterresels(k) * resels << endl;
		struct dcl_params params = {V, Rd, (double)clusterresels(k) * resels, ZestThr};
		pvc = pvox_clust(thres, &params);
	    }

	    clusterpvox(k) = pvc; //applying GRF approach
	    if (false && _verbose) std::cout << p_thres << " " << thres << " " << clustersizes(k) << " " << clusterpvox(k) << std::endl;
	}
	copyconvert(ccc, pvoxclust);
	pvoxclust.projectClusterValues(clusterpvox, 1.0f);

	PVC.addvolume(pvoxclust);
    }

    if (_verbose) std::cout << "...done" << std::endl;

    if (false && _verbose)
    {
	/*
	std::cout << "PVC x\t" << PVC.ColumnsX << " " << PVC.Xdim << std::endl;
	std::cout << "PVC y\t" << PVC.RowsY << " " << PVC.Ydim << std::endl;
	std::cout << "PVC z\t" << PVC.SlicesZ << " " << PVC.Zdim << std::endl;
	std::cout << "PVC t\t" << PVC.dim4 << " " << PVC.p_TR << std::endl;
	std::cout << "PVC min\t" << PVC.getDisplayMinimum() << std::endl;
	std::cout << "PVC max\t" << PVC.getDisplayMaximum() << std::endl;
	*/
    }

    for (int64_t x = 0; x < pTFCEimg.xsize(); ++x)
	for (int64_t y = 0; y < pTFCEimg.ysize(); ++y)
	    for (int64_t z = 0; z < pTFCEimg.zsize(); ++z)
	    {
		if (mask.value(x,y,z) == 0.0) continue;
		
		//aggregate logp values
		NEWMAT::ColumnVector allpvox;
		double aggr;
		allpvox = PVC.voxelts(x, y, z);

		aggr = aggregateLogp(allpvox);
		double p = exp(-aggr);
		if (p == 0.0) p = std::numeric_limits<double>::min();
		if (p == 1.0) p = 1.0 - std::numeric_limits<double>::epsilon();

		pTFCEimg.value(x,y,z)  = p;
		logp_pTFCE.value(x,y,z) = aggr;
		Z_pTFCE.value(x,y,z)  = qnormR(p, 0.0, 1.0, false, false);
	    }

    if (autosmooth || resels > 0)
    {
	number_of_resels = V / resels;
	fwer005_Z = fwerp2z(number_of_resels, 0.05);
    }
    else
    {
	if (resels == 0.0)
	{
	    if (_verbose) std::cout << "For GRF-based FWER correction, please specify resels, or use smoothness estimation based on the data, by not specifying Rd and V!" << std::endl;
	    number_of_resels = 0.0;
	    fwer005_Z = 0.0;
	}
    }

    if (_verbose) std::cout << "number_of_resels: " << number_of_resels << std::endl;
    if (_verbose) std::cout << "fwer_0.05_Z: " << fwer005_Z << std::endl;

    return 0;
}



template <class T>
void pTFCE<T>::setThresholdCount(int N)
{
    if (Nh > 1)
    {
	Nh = N;
	dh = (logpmax - logpmin) / (Nh-1);
    }
    else
    {
	std::cout << "Warning: Nh must be greater than 1" << std::endl;
    }
}

template <class T>
int pTFCE<T>::getThresholdCount()
{
    return this->Nh;
}

template <class T>
void pTFCE<T>::setSmoothness(double Rd, unsigned long V, double resels)
{
    this->Rd     = Rd;
    this->V      = V;
    this->resels = resels;
    autosmooth   = 0;
    if (Rd == 0.0 || V == 0) autosmooth = 1;
}

template <class T>
void pTFCE<T>::setRPVEstimationMode(int mode)
{
    this->RPVMode = mode;
}

template <class T>
int pTFCE<T>::getRPVEstimationMode()
{
    return this->RPVMode;
}

template <class T>
void pTFCE<T>::printSmoothness()
{
    std::cout << "global smoothness: V(" << V << "), Rd(" << Rd << "), dLh(" << dLh << "), resels(" << resels << ")" << std::endl;
}

template <class T>
void pTFCE<T>::loadRPV(const string& filename)
{
    if(_verbose) cerr << "Reading RPV image...." << endl;
    read_volume(this->RPV, filename);
    if(_verbose) print_volume_info(this->RPV,"RPV");

    this->RPVMode = NORPV;
    this->adjustClusterSize = true;
}

template <class T>
void pTFCE<T>::saveRPV(const string& filename)
{
    //TODO - check if RPV image is valid
    //save_volume( mask_volume(smooth(this->RPV, this->img.xdim()), this->mask), filename );
    save_volume( this->RPV, filename );
}

template <class T>
void pTFCE<T>::saveFWHM(const string& filename)
{
    //TODO - check if FWHM image is valid
    //save_volume( mask_volume(smooth(this->FWHMimg, this->img.xdim()), this->mask), filename );
    save_volume( this->FWHMimg, filename );
}

////template <class T>
////void pTFCE<T>::saveLABEL(const string& filename)
////{
////    //TODO - check if LABEL image is valid
////    save_volume( this->LABEL, filename );
////}

template <class T>
void pTFCE<T>::setRFTAdjust(bool a)
{
    adjustClusterSize = a;
}

template <class T>
int pTFCE<T>::getRFTAdjust( )
{
    return this->adjustClusterSize;
}

template <class T>
void pTFCE<T>::setZestThr(double Z)
{
    ZestThr = Z;
}

template <class T>
void pTFCE<T>::quiet()
{
    _verbose = false;
}

template <class T>
void pTFCE<T>::verbose()
{
    _verbose = true;
}

template class pTFCE<float>;
//template class pTFCE<double>;


