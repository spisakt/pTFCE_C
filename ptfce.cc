#include <iostream>
#include <cmath>
#include <limits>
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
	    smoothest(dLh, V, resels, FWHM, FWHMmm, sigmasq,
	              img, mask,
	              100.0, _verbose);
	    Rd = dLh * V;
	    break;
	}
	case 2:
	{
	    if (_verbose) std::cout << "Computing smoothness from 4D residual image" << std::endl;
	    smoothest(dLh, V, resels, FWHM, FWHMmm, sigmasq,
	              residual, mask,
	              dof, _verbose);
	    Rd = dLh * V;
	    break;
	}
	default: break;
    }
}


template <class T>
double pTFCE<T>::aggregateLogp(NEWMAT::ColumnVector &pvals)
{

    double aggr=0;
    for (int i = 1; i <= pvals.n_rows; ++i)
    {
	    pvals(i) = -1 * log(pvals(i));
	    if (isinf(pvals(i))) pvals(i) = 745;
	    aggr+=pvals(i);
    }
    //aggr = armawrap::Sum(pvals);
    aggr = 0.5 * (sqrt(dh * (8.0 * aggr + dh)) - dh);
    // return the logp
    //aggr = exp(-1 * aggr);
    return aggr;
}


template <class T>
int pTFCE<T>::calculate()
{
    // TODO check for NA/nan values
    // TODO check for (-)inf values

    volume4D<T> PVC;        PVC.copyproperties(img);
    copyconvert(img, pTFCEimg);
    copyconvert(img, logp_pTFCE);
    copyconvert(img, Z_pTFCE);
    volume<T> thr;          copyconvert(img, thr);
    volume<int> ccc;        copyconvert(img, ccc);
    volume4D<int> CLUST;    CLUST.copyproperties(ccc);
    volume4D<int> LABEL;    LABEL.copyproperties(ccc);
    volumeClust<int> sizes;
    volumeClust<T> pvoxclust;

    img = mask_volume(img, mask);

    if (_verbose) std::cout << "* Performing pTFCE..." << std::endl;

    for (unsigned int i = 0; i < Nh; ++i)
    {
	double p_thres = exp( -1 * ( logpmin + i * dh ) );
	double thres = qnormR(p_thres, 0.0, 1.0, false, false);
	NEWMAT::ColumnVector clustersizes;
	NEWMAT::ColumnVector clusterpvox;

	if (_verbose) std::cout << i << std::endl;
	if (logpmin == 0.0 && i == 0)
	    thres = std::numeric_limits<T>::lowest();

	copyconvert(img, thr);
	thr.binarise(thres);
	ccc = connected_components(thr, clustersizes, 6);  //6, 18 or 26

	if (_verbose)
	{
	    copyconvert(ccc, sizes);
	    sizes.projectClusterValues(clustersizes, 0.0f);
	    CLUST.addvolume(sizes);
	    LABEL.addvolume(ccc);
	}

	clusterpvox.ReSize(clustersizes.n_rows);
	for (int i = 1; i <= clusterpvox.n_rows; ++i)
	{
	    struct dcl_params params = {V, Rd, clustersizes(i), ZestThr};
	    double pvc = pvox_clust(thres, &params);

	    clusterpvox(i) = pvc; //applying GRF approach
	    if (_verbose) std::cout << p_thres << " " << thres << " " << clustersizes(i) << " " << clusterpvox(i) << std::endl;
	}
	copyconvert(ccc, pvoxclust);
	pvoxclust.projectClusterValues(clusterpvox, 1.0f);

	PVC.addvolume(pvoxclust);
    }

    if (_verbose)
    {
	save_volume(PVC, "testdata/PVC.nii.gz");
	save_volume(CLUST, "testdata/CLUST.nii.gz");
	save_volume(LABEL, "testdata/LABEL.nii.gz");
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
		//if (mask.value(x,y,z) > 0.0)
		//{
		    //aggregate logp values
		    NEWMAT::ColumnVector allpvox;
		    double aggr;
		    allpvox = PVC.voxelts(x, y, z);
		    aggr = aggregateLogp(allpvox);

		    logp_pTFCE.value(x,y,z) = aggr;
            //}
            //else
            //{
            //    pTFCEimg.value(x,y,z) = std::numeric_limits<T>::min();
            //}
            double p = exp(-aggr);
            pTFCEimg.value(x,y,z)  = p;
            if (pTFCEimg.value(x,y,z) == 0) pTFCEimg.value(x,y,z)=std::numeric_limits<T>::min();
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
	    std::cout << "For GRF-based FWER correction, please specify resels, or use smoothness estimation based on the data, by not specifying Rd and V!" << std::endl;
	    number_of_resels = 0.0;
	    fwer005_Z = 0.0;
	}
    }

    std::cout << "number_of_resels: " << number_of_resels << std::endl;
    std::cout << "fwer_0.05_Z: " << fwer005_Z << std::endl;

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
void pTFCE<T>::setSmoothness(double Rd, unsigned long V, double resels)
{
    this->Rd     = Rd;
    this->V      = V;
    this->resels = resels;
    autosmooth   = 0;
    if (Rd == 0.0 || V == 0) autosmooth = 1;
}

template <class T>
void pTFCE<T>::printSmoothness()
{
    std::cout << "smoothness: V(" << V << "), Rd(" << Rd << "), dLh(" << dLh << "), resels(" << resels << ")" << std::endl;
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



