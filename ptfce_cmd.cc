#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>
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

static void appendLineToFile(string filepath, string line);
static void makeDirectoryPortable(string dirname);

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
Option<float> gRd(string("-x,--rd"), 0.0,
	  string("\tglobal smoothness: Rd"),
	  false, requires_argument);
Option<int> gVoxels(string("-y,--voxels"), 0,
	  string("global smoothness: Voxels"),
	  false, requires_argument);
Option<float> gResels(string("-w,--resels"), 0.0,
	  string("global smoothness: Resels"),
	  false, requires_argument);
Option<int> rpvmode(string("-s,--smmode"), 0,
	  string("local smoothness estimation mode 0-none, 1-fsl, 2-spm"),
	  false, requires_argument);
Option<string> operationtime(string("-o,--otime"), "otime",
	     string("filename to save operation time into"),
	     false, requires_argument);


string title = "\
pTFCE ";

string examples = "\
\tpTFCE -d <number> -r <filename> -m <filename>\n\
\tpTFCE -z <filename> -m <filename>";

int main(int argc, char* argv[])
{

    int status = 0;

    // 1. Parse options

    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> tSmoothest;
    std::chrono::duration<double> tPTFCE;

    OptionParser options(title, examples);
    options.add(zstatname);
    options.add(maskname);
    options.add(residname);
    options.add(dof);
    options.add(nthr);
    options.add(rpvmode);
    options.add(rpvname);
    options.add(gRd);
    options.add(gVoxels);
    options.add(gResels);
    options.add(operationtime);
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


    // 2. Load data
    NEWIMAGE::volume<float> img;
    NEWIMAGE::volume4D<float> R;
    NEWIMAGE::volume<float> mask;
    float _dof = dof.value();

    if (verbose.value()) cout << "Loading data" << endl;
    if (verbose.value()) cout << "********************************" << endl;

    // Read the mask image (single volume)
    if (verbose.value()) cerr << "Reading mask...." << endl;
    read_volume(mask, maskname.value());
    if (verbose.value()) print_volume_info(mask,"mask");

    // Read the zstat image to boost (single volume)
    if (verbose.value()) cerr << "Reading image...." << endl;
    read_volume(img, zstatname.value());
    if (verbose.value()) print_volume_info(img,"Zstat");

    // Read the residual images for smoothness estimation (array of volumes)
    if (residname.set())
    {
        if (verbose.value()) cerr << "Reading residuals...." << endl;
        read_volume4D(R, residname.value());
        if (verbose.value()) print_volume_info(R,"residuals");
        if (dof.unset())
        {
            if (verbose.value()) cerr << "Dof unspecified. Using residual length as degrees of freedom: " << R.tsize() << endl;
            _dof = R.tsize();
        }
    }
    else
    {
        if (verbose.value()) cerr << "Residuals not found. Using Zstat for smoothness estimation...." << endl;;
        read_volume4D(R, zstatname.value());
        if (verbose.value()) print_volume_info(R,"Zstat");
    }

    if ( !samesize(R[0], mask) || !samesize(img, mask) )
    {
        cerr << "Mask and Data (residuals/zstat) volumes MUST be the same size!" << endl;
        exit(EXIT_FAILURE);
    }


    // 3. Settings
    if (verbose.value()) cout << "Settings" << endl;
    if (verbose.value()) cout << "********************************" << endl;


    // 3.a) Define pTFCE
    pTFCE<float> enhance(img, mask, R, _dof);

    enhance.setThresholdCount(nthr.value());
    if (verbose.value()) enhance.verbose();


    // 3.b) Smoothness set-up or estimation
    if (gRd.set() && gVoxels.set() && gResels.set())
        enhance.setSmoothness( gRd.value(), gVoxels.value(), gResels.value() );
    // TODO - cover cases when global smoothness is set, but local calculations are needed
    enhance.setRPVEstimationMode(rpvmode.value());
    if (rpvmode.value() != 0)
        enhance.setRFTAdjust(true);
    else
        enhance.setRFTAdjust(false);
    if (rpvname.set())
        enhance.loadRPV( rpvname.value() );
    t1 = std::chrono::high_resolution_clock::now();
    enhance.estimateSmoothness();
    t2 = std::chrono::high_resolution_clock::now();
    tSmoothest = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    if (verbose.value()) enhance.printSmoothness();
    if (verbose.value()) cout << "Running time (smoothness): " << tSmoothest.count() << endl;


    // 3.c) Check settings
    if (true || verbose.value())
    {
	cout << "verbose = " << verbose.value() << endl;
	cout << "verbose = " << verbose.value() << endl;
	cout << "help = " << help.value() << endl;
	cout << "dof = " << dof.value() << endl;
	cout << "maskname = " << maskname.value() << endl;
	cout << "residname = " << (residname.set() ? residname.value() : "-") << endl;
	cout << "zstatname = " << (zstatname.set() ? zstatname.value() : "-") << endl;
	cout << "rpvname = " << (rpvname.set() ? rpvname.value() : "-") << endl;
	cout << "rpvmode = " << enhance.getRPVEstimationMode() << endl;
	cout << "rpvadjustment = " << enhance.getRFTAdjust() << endl;
	cout << "thresholdCount = " << nthr.value() << endl;
	cout << "timelogFile = " << (operationtime.set() ? operationtime.value() : "-") << endl;
    }


    // 4. Calculate pTFCE
    t1 = std::chrono::high_resolution_clock::now();
    enhance.calculate();
    t2 = std::chrono::high_resolution_clock::now();
    tPTFCE = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    if (verbose.value()) cout << "Running time (pTFCE): " << tPTFCE.count() << endl;


    // 5. Save results
    std::string thisDir ("./");
    std::string zstatName (zstatname.value());
    if ( zstatName.find_last_of("/\\") == std::string::npos ) zstatName = thisDir + zstatName;
    size_t zstatPathSep = zstatName.find_last_of("/\\");
    string zstatPath = zstatName.substr(0,zstatPathSep);
    string zstatFile = zstatName.substr(zstatPathSep);
    size_t zstatExtSep = zstatFile.find_last_of(".nii");
    string zstatFileBase = zstatFile.substr(0,zstatExtSep-4);
    if ( zstatPath == zstatFile ) zstatPath = ".";
    string outdir = zstatPath + "/" + "pTFCE";
    makeDirectoryPortable(outdir);

    string fn_p, fn_lp, fn_z, fn_rpv, fn_fwhm;
    switch( enhance.getRPVEstimationMode() )
    {
	case NORPV:
	{
	    if ( ! enhance.getRFTAdjust() )
	    {
		fn_p  = outdir + "/" + zstatFileBase + "_global_pTFCE.nii.gz";
		fn_lp = outdir + "/" + zstatFileBase + "_global_logpTFCE.nii.gz";
		fn_z  = outdir + "/" + zstatFileBase + "_global_ZpTFCE.nii.gz";
	    }
	    else
	    {
		fn_p  = outdir + "/" + zstatFileBase + "_local_pTFCE.nii.gz";
		fn_lp = outdir + "/" + zstatFileBase + "_local_logpTFCE.nii.gz";
		fn_z  = outdir + "/" + zstatFileBase + "_local_ZpTFCE.nii.gz";
	    }
	} break;
	case FSLRPV:
	{
	    fn_p  = outdir + "/" + zstatFileBase + "_fslrpv_pTFCE.nii.gz";
	    fn_lp = outdir + "/" + zstatFileBase + "_fslrpv_logpTFCE.nii.gz";
	    fn_z  = outdir + "/" + zstatFileBase + "_fslrpv_ZpTFCE.nii.gz";
	    fn_rpv  = outdir + "/" + zstatFileBase + "_fslrpv_RPV.nii.gz";
	    fn_fwhm = outdir + "/" + zstatFileBase + "_fslrpv_FWHM.nii.gz";
	} break;
	case SPMRPV:
	{
	    fn_p  = outdir + "/" + zstatFileBase + "_spmrpv_pTFCE.nii.gz";
	    fn_lp = outdir + "/" + zstatFileBase + "_spmrpv_logpTFCE.nii.gz";
	    fn_z  = outdir + "/" + zstatFileBase + "_spmrpv_ZpTFCE.nii.gz";
	    fn_rpv  = outdir + "/" + zstatFileBase + "_spmrpv_RPV.nii.gz";
	    fn_fwhm = outdir + "/" + zstatFileBase + "_spmrpv_FWHM.nii.gz";
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
    ////if ( verbose.value() )
    ////{
    ////    string fn_label = outdir + "/" + zstatFileBase + "_LABEL_4D.nii.gz";
    ////    enhance.saveLABEL(fn_label);
    ////}

    if ( operationtime.set() )
    {
	std::stringstream tLineStream;
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	tLineStream << std::put_time(&tm, "%d-%m-%Y_%H-%M-%S") << ","
	            << zstatname.value() << ","
	            << enhance.getRPVEstimationMode() << ","
	            << enhance.getThresholdCount() << ","
	            << tSmoothest.count() << ","
	            << tPTFCE.count();
	appendLineToFile(operationtime.value(), tLineStream.str());
    }

    return status;
}


static void appendLineToFile(string filepath, string line)
{
    std::ofstream file;
    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open(filepath, std::ios::out | std::ios::app);
    if (file.fail())
        throw std::ios_base::failure(std::strerror(errno));

    //make sure write fails with exception if something is wrong
    file.exceptions(file.exceptions() | std::ios::failbit | std::ifstream::badbit);

    file << line << std::endl;
}

static void makeDirectoryPortable(string dirname)
{
    mode_t nMode = 0733;
    int nError = 0;
#if defined(_WIN32)
	nError = _mkdir(dirname.c_str()); // can be used on Windows
#else
	nError = mkdir(dirname.c_str(),nMode); // can be used on non-Windows
#endif
    if (nError != 0) {
    // handle your error here
    }
}
