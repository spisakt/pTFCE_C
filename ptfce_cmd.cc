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

#include "utils/options.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace Utilities;
using namespace NEWIMAGE;

Option<bool> verbose(string("-v,--verbose"), false,
	     string("switch on diagnostic messages"),
	     false, no_argument);
Option<bool> help(string("-h,--help"), false,
	  string("display this message"),
	  false, no_argument);
Option<float> dof(string("-d,--dof"), 100.0,
	  string("number of degrees of freedom"),
	  false, requires_argument);
Option<int> nthr(string("-t,--thrs"), 100,
	  string("number of Z value thresholds"),
	  false, requires_argument);
Option<string> maskname(string("-m,--mask"), "mask",
	    string("brain mask volume"),
	    true, requires_argument);
Option<string> zstatname(string("-z,--zstat"), "zstat",
	     string("filename of zstat image (not with -d)"),
	     true, requires_argument);
Option<string> residname(string("-r,--res"), "res4d",
	     string("filename of `residual-fit' image (use -d)"),
	     false, requires_argument);
Option<string> rpvname(string("-l,--rpv"), "rpv",
	     string("local roughness image (RPV)"),
	     false, requires_argument);
Option<int> rpvmode(string("-s,--smmode"), 0,
	  string("local smoothness estimation mode 0-none, 1-fsl, 2-spm"),
	  false, requires_argument);


string title = "\
pTFCE ";

string examples = "\
\tpTFCE -d <number> -r <filename> -m <filename>\n\
\tpTFCE -z <filename> -m <filename>";

int main(int argc, char* argv[])
{

    int status = 0;

    OptionParser options(title, examples);
    options.add(zstatname);
    options.add(maskname);
    options.add(residname);
    options.add(dof);
    options.add(nthr);
    options.add(rpvmode);
    options.add(rpvname);
    options.add(verbose);
    options.add(help);
    options.parse_command_line(argc, argv);


    if(help.value())
    {
	options.usage();
	exit(EXIT_SUCCESS);
    }

    if(   (zstatname.unset()) ||
          (maskname.unset()) /*||
          (residname.set() && dof.unset())*/ )
    {
	options.usage();
	cerr << endl;
	cerr << "***************************************************************************" << endl;
	cerr << "You must specify a zstat image." << endl;
	cerr << "You must specify a mask volume image filename" << endl; 
	//cerr << "If you are processing a 4D residual image, then you must specify the degrees of freedom." << endl;
	cerr << "***************************************************************************" << endl;
	cerr << endl;
	exit(EXIT_FAILURE);
    }


    // Load data
    NEWIMAGE::volume<float> img;
    NEWIMAGE::volume4D<float> R;
    NEWIMAGE::volume<float> mask;
    float _dof = dof.value();

    if (verbose.value()) cout << "Loading data" << endl;
    if (verbose.value()) cout << "********************************" << endl;

    // Read the mask image (single volume)
    if(verbose.value()) cerr << "Reading mask...." << endl;
    read_volume(mask, maskname.value());
    if(verbose.value()) print_volume_info(mask,"mask");

    // Read the zstat image to boost (single volume)
    if(verbose.value()) cerr << "Reading image...." << endl;
    read_volume(img, zstatname.value());
    if(verbose.value()) print_volume_info(img,"Zstat");

    // Read the residual images for smoothness estimation (array of volumes)
    if(residname.set())
    {
        if(verbose.value()) cerr << "Reading residuals...." << endl;
        read_volume4D(R, residname.value());
        if(verbose.value()) print_volume_info(R,"residuals");
        if(dof.unset())
        {
            if(verbose.value()) cerr << "Dof unspecified. Using residual length as degrees of freedom: " << R.tsize() << endl;
            _dof = R.tsize();
        }
    }
    else
    {
        if(verbose.value()) cerr << "Residuals not found. Using Zstat for smoothness estimation...." << endl;;
        read_volume4D(R, zstatname.value());
        if(verbose.value()) print_volume_info(R,"Zstat");
    }

    if ( !samesize(R[0], mask) || !samesize(img, mask) )
    {
        cerr << "Mask and Data (residuals/zstat) volumes MUST be the same size!" << endl;
        exit(EXIT_FAILURE);
    }

    // Settings
    if (verbose.value()) cout << "Settings" << endl;
    if (verbose.value()) cout << "********************************" << endl;

    // Perform pTFCE
    pTFCE<float> enhance(img, mask, R, _dof);

    if (verbose.value()) enhance.verbose();

    enhance.setRPVEstimationMode(rpvmode.value());
    if (rpvmode.value() != 0)
        enhance.setRFTAdjust(true);
    else
        enhance.setRFTAdjust(false);
    if (rpvname.set())
        enhance.loadRPV( rpvname.value() );
    enhance.estimateSmoothness();
    if (verbose.value()) enhance.printSmoothness();

    enhance.setThresholdCount(nthr.value()); //deafult in R package


    if(true || verbose.value())
    {
	cout << "verbose = " << verbose.value() << endl;
	cout << "help = " << help.value() << endl;
	cout << "dof = " << dof.value() << endl;
	cout << "maskname = " << maskname.value() << endl;
	cout << "residname = " << residname.value() << endl;
	cout << "zstatname = " << zstatname.value() << endl;
	cout << "rpvname = " << rpvname.value() << endl;
	cout << "rpvmode = " << enhance.getRPVEstimationMode() << endl;
	cout << "rpvadjustment = " << enhance.getRFTAdjust() << endl;
	cout << "thresholdCount = " << nthr.value() << endl;
    }


    enhance.calculate();

    // Save enhanced volume TODO
    string fn_p, fn_lp, fn_z, fn_rpv, fn_fwhm;
    switch( enhance.getRPVEstimationMode() )
    {
	case NORPV:
	{
	    if ( ! enhance.getRFTAdjust() )
	    {
		fn_p  = "testdata/pTFCE_global.nii.gz";
		fn_lp = "testdata/logpTFCE_global.nii.gz";
		fn_z  = "testdata/ZpTFCE_global.nii.gz";
	    }
	    else
	    {
		fn_p  = "testdata/pTFCE_local.nii.gz";
		fn_lp = "testdata/logpTFCE_local.nii.gz";
		fn_z  = "testdata/ZpTFCE_local.nii.gz";
	    }
	} break;
	case FSLRPV:
	{
	    fn_p  = "testdata/pTFCE_fslrpv.nii.gz";
	    fn_lp = "testdata/logpTFCE_fslrpv.nii.gz";
	    fn_z  = "testdata/ZpTFCE_fslrpv.nii.gz";
	    fn_rpv  = "testdata/RPV_fslrpv.nii.gz";
	    fn_fwhm = "testdata/FWHM_fslrpv.nii.gz";
	} break;
	case SPMRPV:
	{
	    fn_p  = "testdata/pTFCE_spmrpv.nii.gz";
	    fn_lp = "testdata/logpTFCE_spmrpv.nii.gz";
	    fn_z  = "testdata/ZpTFCE_spmrpv.nii.gz";
	    fn_rpv  = "testdata/RPV_spmrpv.nii.gz";
	    fn_fwhm = "testdata/FWHM_spmrpv.nii.gz";
	} break;
    }

    save_volume(enhance.pTFCEimg, fn_p);
    save_volume(enhance.logp_pTFCE, fn_lp);
    save_volume(enhance.Z_pTFCE, fn_z);
    if ( enhance.getRPVEstimationMode() != 0 )
    {
	enhance.saveRPV(fn_rpv);
	enhance.saveFWHM(fn_fwhm);
    }

    return status;
}
