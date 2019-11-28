#include <iostream>
#include <cmath>
#include <string>
#include <map>

#include "smoothest_ext.h"

#include "utils/options.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace Utilities;
using namespace NEWIMAGE;

Option<bool> verbose(string("-V,--verbose"), false,
	     string("switch on diagnostic messages"),
	     false, no_argument);
Option<bool> help(string("-h,--help"), false,
	  string("display this message"),
	  false, no_argument);
Option<float> dof(string("-d,--dof"), 100.0,
	  string("number of degrees of freedom"),
	  true, requires_argument);
Option<int> local(string("-l,--rpv"), 0,
	  string("local smoothness estimation (0* - none; 1 - FSLRPV; 2 - SPMRPV)"),
	  true, requires_argument);
Option<string> maskname(string("-m,--mask"), "mask",
	    string("brain mask volume"),
	    true, requires_argument);
Option<string> zstatname(string("-z,--zstat"), "zstat",
	     string("filename of zstat image (not with -d)"),
	     true, requires_argument);
Option<string> residname(string("-r,--res"), "res4d",
	     string("filename of `residual-fit' image (use -d)"),
	     true, requires_argument);


string title = "\
smoothest \nCopyright(c) 2000-2002, University of Oxford (Dave Flitney and Mark Jenkinson)";

string examples = "\
\tsmoothest -d <number> -r <filename> -m <filename>\n\
\tsmoothest -z <filename> -m <filename>";

int main(int argc, char **argv) {

  int status;

  double dlh;
  unsigned long mask_volume;
  double resels;
  double FWHM_vox[3];
  double FWHM_mm[3];
  double sigmasq[3];
  volume<float> RPV;
  volume<float> FWHMimg;

  OptionParser options(title, examples);

  options.add(verbose);
  options.add(help);
  options.add(dof);
  options.add(local);
  options.add(maskname);
  options.add(residname);
  options.add(zstatname);

  options.parse_command_line(argc, argv);

//  if(verbose.value()) options.check_compulsory_arguments(true);

  if(help.value()) {
    options.usage();
    exit(EXIT_SUCCESS);
  }

  if( !((zstatname.set() && residname.unset()) || (residname.set() && zstatname.unset())) ||
      (zstatname.set() && (residname.set() || dof.set())) || 
      (residname.set() && dof.unset()) ||
      (maskname.unset()) ) {
    options.usage();
    cerr << endl;
    cerr << "***************************************************************************" << endl;
    cerr << "You must specify either a zstat image OR a 4d residual image." << endl;
    cerr << "If processing a zstat image then you should not set degrees of freedom." << endl; 
    cerr << "You must specify a mask volume image filename" << endl; 
    cerr << "You must specify the degrees of freedom for processing a 4d residual image." << endl; 
    cerr << "***************************************************************************" << endl;
    cerr << endl;
    exit(EXIT_FAILURE);    
  }

  if(verbose.value()) {
    cout << "verbose = " << verbose.value() << endl;
    cout << "help = " << help.value() << endl;
    cout << "dof = " << dof.value() << endl;
    cout << "maskname = " << maskname.value() << endl;
    cout << "residname = " << residname.value() << endl;
    cout << "zstatname = " << zstatname.value() << endl;
  }

  // Read the AVW mask image (single volume)
  if(verbose.value()) cerr << "Reading mask....";
  volume<float> mask;
  read_volume(mask,maskname.value());
  if(verbose.value()) cerr << "done" << endl;
  if (verbose.value()) print_volume_info(mask,"mask");

  if(verbose.value()) cerr << "Reading datafile....";
  string datafilename;
  if(residname.set()) {
    // Read the AVW residual images (array of volumes)
    datafilename = residname.value();
  } else {
    // Read the AVW zstat image (array of one volume)
    datafilename = zstatname.value();
  }
  if(verbose.value()) cerr << "done" << endl;


  volume4D<float> R;
  read_volume4D(R,datafilename);
  if (verbose.value()) print_volume_info(R,"Data (residuals/zstat)");

  if (!samesize(R[0],mask)) {
    cerr << "Mask and Data (residuals/zstat) volumes MUST be the same size!"
     << endl;
    exit(EXIT_FAILURE);
  }


  status = smoothest(dlh, mask_volume, resels, FWHM_vox, FWHM_mm, sigmasq,
	     R,
	     mask,
	     dof.value(), verbose.value());

  if (local.value() == 1)
  {
    smoothestVox(dlh, mask_volume, resels, FWHM_vox, FWHM_mm, sigmasq, RPV, FWHMimg,
                  R, mask,
                  dof.value(), verbose.value());
    save_volume(RPV, "testdata/smoothest_app_fslrpv.nii.gz");

  }
  if (local.value() == 2)
  {
    estimateRPV(RPV, FWHMimg,
                  R, mask,
                  dof.value(), verbose.value());
    save_volume(RPV, "testdata/smoothest_app_spmrpv.nii.gz");
  }


  return status;
}


//g++ -o smoothest smoothest_ext.o smoothest_app.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/
// -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack
