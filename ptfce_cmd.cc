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
    // Load data
    NEWIMAGE::volume4D<float> R;
    NEWIMAGE::volume<float> mask;
    read_volume(R,    argv[1]);
    read_volume(mask, argv[2]);

    // perform pTFCE
    pTFCE<float> enhance(R, mask);
    enhance.printSmoothness();
    enhance.setThresholdCount(100); //deafult in R package
    enhance.calculate();

    // save enhanced volume
    save_volume(enhance.pTFCEimg, "mytest/pTFCE.nii.gz");
    save_volume(enhance.logp_pTFCE, "mytest/logpTFCE.nii.gz");
    save_volume(enhance.Z_pTFCE, "mytest/ZpTFCE.nii.gz");
}