#include <iostream>
#include <cmath>
#include "ndf.h"
#include "newimage/newimageall.h"
#include "smoothest_ext.h"
#include "pvalutil.h"
#include "grfClust.h"
#include "ptfce.h"

#if !defined(M_PI)
#define M_PI (4 * atan(1.0))
#endif

int main(int argc, char* argv[])
{


    //TEST general functions
    std::cout << std::endl << "TEST general functions" << std::endl;
    std::cout << "dens " << dnorm(1.5) << std::endl;
    std::cout << "dist " << pnorm(1.5) << std::endl;
    std::cout << "quant " << qnorm(0.05) << std::endl;
    std::cout << "quant uppertail " << qnorm(0.05, 0.0, 1.0, false, false) << std::endl;
    std::cout << "quant logp" << qnorm(-7.0, 0.0, 1.0, true, true) << std::endl;
    std::cout << "quant logp uppertail" << qnorm(-7.0, 0.0, 1.0, false, true) << std::endl;
    std::cout << "gamma(2.5) " << tgamma(2.5) << std::endl;
    std::cout << "log(0.0) " << log(0.0) << std::endl;
    std::cout << "qnorm(0.0) " << qnorm(0.0) << std::endl;


    //TEST pvox_clust, dvox_clust, dculst, dcl, Es
    std::cout << std::endl << "TEST pvox_clust, dvox_clust, dculst, dcl, Es" << std::endl;
    //gsl_set_error_handler_off();
    struct dcl_params params = {20, 60, 10, 1.3};
    double actH = 1.5;
    std::cout << "Es(" << actH << ", 20, 60) " << Es(actH, 20, 60) << std::endl;
    std::cout << "dcl(" << actH << ", 20, 60, 10, 1.3) R:0.02529095 Cpp:" << dcl( actH, &params ) << std::endl;
    std::cout << "dclust(" << actH << ", 20, 60, 10, 1.3) R:2.000079 Cpp:" << dclust( actH, &params ) << std::endl;
    std::cout << "dvox_clust(" << actH << ", 20, 60, 10, 1.3) R:1.955195 Cpp:" << dvox_clust( actH, &params ) << std::endl;
    std::cout << "pvox_clust(" << actH << ", 20, 60, 10, 1.3) R:0.2904333 Cpp:" << pvox_clust( actH, &params ) << std::endl;


    //TEST smoothest
    std::cout << std::endl << "TEST smoothest" << std::endl;
    double dLh, resels;
    double FWHM[3];
    double FWHMmm[3];
    double sigmasq[3];
    unsigned long mask_volume;
    NEWIMAGE::volume4D<float> R;
    NEWIMAGE::volume<float> mask;
    enum {X = 0, Y, Z};

    read_volume(mask, "testdata/mask.nii.gz");
    read_volume(R,    "testdata/zstat3.nii.gz");

    smoothest(dLh, mask_volume, resels, FWHM, FWHMmm, sigmasq,
              R, mask,
              100.0, false);

    std::cout << "DLH " << dLh << std::endl;
    std::cout << "VOLUME " << mask_volume << std::endl;
    std::cout << "RESELS " << resels << std::endl;
    std::cout << "FWHMvoxel " << FWHM[X] << " " <<  FWHM[Y] << " " << FWHM[Z] << std::endl;
    std::cout << "FWHMmm " << FWHMmm[X] <<  " " << FWHMmm[Y] << " " << FWHMmm[Z] << std::endl;
    std::cout << "sigmasq " << sigmasq[X] << " " << sigmasq[Y] << " " << sigmasq[Z] << std::endl;


    //TEST z2p p2z
    std::cout << std::endl << "TEST z2p p2z" << std::endl;
    double pez;
    pez = fwerp2z(resels);
    pez = fwerz2p(resels, pez);


    //TEST pTFCE
    std::cout << std::endl << "TEST pTFCE" << std::endl;
    pTFCE<float> enhance(R, mask);
    enhance.printSmoothness();
    enhance.saveRPV("testdata/RPV.nii.gz");
    enhance.setThresholdCount(12);
    enhance.calculate();
    save_volume(enhance.pTFCEimg, "testdata/pTFCE.nii.gz");
    save_volume(enhance.logp_pTFCE, "testdata/logpTFCE.nii.gz");
    save_volume(enhance.Z_pTFCE, "testdata/ZpTFCE.nii.gz");

    enhance.setRFTAdjust(true);
    enhance.calculate();
    save_volume(enhance.pTFCEimg, "testdata/pTFCE_rpv.nii.gz");
    save_volume(enhance.logp_pTFCE, "testdata/logpTFCE_rpv.nii.gz");
    save_volume(enhance.Z_pTFCE, "testdata/ZpTFCE_rpv.nii.gz");

    return 0;
}

//g++ -o pTFCE-test smoothest_ext.o fwerp2z.o pTFCE-test.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/
// -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack -lgsl -lgslcblas

