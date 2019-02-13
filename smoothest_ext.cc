#include <iostream>
#include <cmath>
#include <string>
#include <map>

#include "smoothest_ext.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1


using namespace NEWIMAGE;


namespace SMOOTHEST {

class Interpolate {
public:
  Interpolate() {
    lut[5]   = 1.5423138; lut[6]   = 1.3757105; lut[7]   = 1.2842680;
    lut[8]   = 1.2272151; lut[9]   = 1.1885232; lut[10]  = 1.1606988;
    lut[11]  = 1.1398000; lut[12]  = 1.1235677; lut[13]  = 1.1106196;
    lut[14]  = 1.1000651; lut[15]  = 1.0913060; lut[16]  = 1.0839261;
    lut[17]  = 1.0776276; lut[18]  = 1.0721920; lut[19]  = 1.0674553;
    lut[20]  = 1.0632924; lut[25]  = 1.0483053; lut[30]  = 1.0390117;
    lut[40]  = 1.0281339; lut[50]  = 1.0219834; lut[60]  = 1.0180339;
    lut[70]  = 1.0152850; lut[80]  = 1.0132621; lut[90]  = 1.0117115;
    lut[100] = 1.0104851; lut[150] = 1.0068808; lut[200] = 1.0051200;
    lut[300] = 1.0033865; lut[500] = 1.0020191;
  }

  inline float operator()(float v) {

    float retval = 0;

    if (v<6) return 1.1; // ?? no idea - steve ??

    map<int, float>::iterator i = lut.lower_bound(int(v));

    if(i != lut.end()) {
      if(i != lut.begin()) {
    map<int, float>::iterator j = i--;

    retval = (j->second - i->second)/(j->first - i->first)*(v - i->first)
      + j->second;
      } 
    } else {
      retval = 1.0321/v + 1;
    }
    
    retval = pow(retval, 0.5);

    return retval;
  }

private:
  map<int, float> lut;
};

Interpolate interpolate;

}

//////////////////////////////////////////////////////////////////////////////
// Standardise the residual field (assuming gaussianity)
unsigned long standardise(volume<float>& mask, 
	      volume4D<float>& R)
{
  unsigned long count = 0;
  int M=R.tsize();

  for (int z=mask.minz(); z<=mask.maxz(); z++) {
    for (int y=mask.miny(); y<=mask.maxy(); y++) {
      for (int x=mask.minx(); x<=mask.maxx(); x++) {

    if( mask(x,y,z) > 0.5) {
      
      count ++;
      
      if( M > 2 ) {
        
        // For each voxel 
        //    calculate mean and standard deviation...
        double Sx = 0.0, SSx = 0.0;
        
        for ( int t = 0; t < M; t++ ) {
          double R_it = R(x,y,z,t);
          
          Sx += R_it;
          SSx += Sqr(R_it);
        }
        
        double mean = Sx / M;
        double sdsq = (SSx - (Sqr(Sx) / M)) / (M - 1) ;
        
        if (sdsq<=0) {
          // trap for differences between mask and invalid data
          mask(x,y,z)=0;
          count--;
        } else {
          //    ... and use them to standardise to N(0, 1).
          for ( unsigned short t = 0; t < M; t++ ) {
	R(x,y,z,t) = (R(x,y,z,t) - mean) / sqrt(sdsq);
          }
        } 
      }
    }
      }
    }
  }  
  return count;
}


int smoothest(double &dLh, unsigned long &mask_volume, double &resels, double *FWHM, double *FWHMmm, double *sigmasq,
	      NEWIMAGE::volume4D<float>& R,
	      NEWIMAGE::volume<float>& mask,
	      double dof, bool verbose)
{
  if(verbose) cerr << "Standardising....";	//aranyics @ every case of verbose.value()
  mask_volume = standardise(mask, R);		//aranyics
  if(verbose) cerr << "done" << endl;

  if(verbose) cerr << "Masked-in voxels = " << mask_volume << endl;

  unsigned long N = 0;

  // MJ additions to make it cope with 2D images
  bool usez = true;
  if (R.zsize() <= 1) { usez = false; }
  if ((!usez) && verbose) {
    cout << "Using 2D image mode." << endl;
  }

  // Estimate the smoothness of the normalised residual field
  // see TR00DF1 for mathematical description of the algorithm.
  enum {X = 0, Y, Z};
  double SSminus[3] = {0, 0, 0}, S2[3] = {0, 0, 0};

  int zstart=1;
  if (!usez) zstart=0;
  for ( unsigned short z = zstart; z < R.zsize() ; z++ )
    for ( unsigned short y = 1; y < R.ysize() ; y++ )
      for ( unsigned short x = 1; x < R.xsize() ; x++ )
	// Sum over N
	if( (mask(x, y, z)>0.5) &&
	    (mask(x-1, y, z)>0.5) &&
	    (mask(x, y-1, z)>0.5) &&
	    ( (!usez) || (mask(x, y, z-1)>0.5) ) ) {
	
	  N++;
	
	  for ( unsigned short t = 0; t < R.tsize(); t++ ) {
	    // Sum over M
	    SSminus[X] += R(x, y, z, t) * R(x-1, y, z, t);
	    SSminus[Y] += R(x, y, z, t) * R(x, y-1, z, t);
	    if (usez) SSminus[Z] += R(x, y, z, t) * R(x, y, z-1, t);
	
	    S2[X] += 0.5 * (Sqr(R(x, y, z, t)) + Sqr(R(x-1, y, z, t)));
	    S2[Y] += 0.5 * (Sqr(R(x, y, z, t)) + Sqr(R(x, y-1, z, t)));
	    if (usez) S2[Z] += 0.5 * (Sqr(R(x, y, z, t)) + Sqr(R(x, y, z-1, t)));
	  }
	}

  double norm = 1.0/(double) N;
  double v = dof;	// v - degrees of freedom (nu) //aranyics
  if(R.tsize() > 1) {
    if(verbose) {
      cerr << "Non-edge voxels = " << N << endl;
      cerr << "(v - 2)/(v - 1) = " << (v - 2)/(v - 1) << endl;
    }
    norm = (v - 2) / ((v - 1) * N * R.tsize());
  }

//    SSminus[X] *= norm;
//    SSminus[Y] *= norm;
//    SSminus[Z] *= norm;

//    S2[X] *= norm;
//    S2[Y] *= norm;
//    S2[Z] *= norm;

  if(verbose) {
    cout << "SSminus[X] = " << SSminus[X] << ", SSminus[Y] = " << SSminus[Y] << ", SSminus[Z] = " << SSminus[Z]
	 << ", S2[X] = " << S2[X] << ", S2[Y] = " << S2[Y] << ", S2[Z] = " << S2[Z]
	 << endl;
  }

  // for extreme smoothness.
  if (SSminus[X]>=0.99999999*S2[X]) {
    SSminus[X]=0.99999*S2[X];
    cerr << "WARNING: Extreme smoothness detected in X - possibly biased"
	 << " global estimate." << endl; }
  if (SSminus[Y]>=0.99999999*S2[Y]) {
    SSminus[Y]=0.99999*S2[Y];
    cerr << "WARNING: Extreme smoothness detected in Y - possibly biased"
	 << " global estimate." << endl; }
  if (usez) {
    if (SSminus[Z]>=0.99999999*S2[Z]) {
      SSminus[Z]=0.99999*S2[Z];
      cerr << "WARNING: Extreme smoothness detected in Z - possibly biased"
	   << " global estimate." << endl; }
  }

  // Convert to sigma squared
  //double sigmasq[3]; //aranyics

  sigmasq[X] = -1.0 / (4 * log(fabs(SSminus[X]/S2[X])));
  sigmasq[Y] = -1.0 / (4 * log(fabs(SSminus[Y]/S2[Y])));
  if (usez) { sigmasq[Z] = -1.0 / (4 * log(fabs(SSminus[Z]/S2[Z]))); }
  else { sigmasq[Z]=0; }
  // the following is determininant of Lambda to the half
  //   i.e. dLh = | Lambda |^(1/2)
  // Furthermore, W_i = 1/(2.lambda_i) = sigma_i^2 =>.
  //   det(Lambda) = det( lambda_i ) = det ( (2 W_i)^-1 ) = (2^D det(W))^-1
  //   where D = number of dimensions (2 or 3)
  //double dLh; //aranyics
  if (usez) { dLh=pow(sigmasq[X]*sigmasq[Y]*sigmasq[Z], -0.5)*pow(8, -0.5); }
  else { dLh = pow(sigmasq[X]*sigmasq[Y], -0.5)*pow(4, -0.5); }

  if(verbose) {
    cout << "DLH " << dLh << " voxels^-3 before correcting for temporal DOF" << endl;
  }
  if(R.tsize() > 1) dLh *= SMOOTHEST::interpolate(v);

  // Convert to full width half maximum
  //double FWHM[3]; //aranyics
  FWHM[X] = sqrt(8 * log(2) * sigmasq[X]);
  FWHM[Y] = sqrt(8 * log(2) * sigmasq[Y]);
  if (usez) { FWHM[Z] = sqrt(8 * log(2) * sigmasq[Z]); }
  else { FWHM[Z]=0; }
  resels = FWHM[X] * FWHM[Y]; //aranyics
  if (usez) resels *= FWHM[Z];

  if(verbose) {
    cout << "FWHMx = " << FWHM[X] << " voxels, "
	 << "FWHMy = " << FWHM[Y] << " voxels";
    if (usez) cout << ", FWHMz = " << FWHM[Z] << " voxels";
    cout << endl;
  }

  //aranyics
  //double FWHMmm[3] = { FWHM[X]*R.xdim(), FWHM[Y]*R.ydim(), FWHM[Z]*R.zdim() };
  FWHMmm[X] = FWHM[X]*R.xdim();
  FWHMmm[Y] = FWHM[Y]*R.ydim();
  FWHMmm[Z] = FWHM[Z]*R.zdim();

  /*if(verbose)
  {
    cout << "FWHMx = " << FWHM[X] << " mm, "
	 << "FWHMy = " << FWHM[Y] << " mm";
    if (usez) cout << ", FWHMz = " << FWHM[Z] << " mm";
    cout << endl;
    cout << "DLH " << dLh << " voxels^-3" << endl;
    cout << "VOLUME " << mask_volume << " voxels" << endl;
    cout << "RESELS " << resels << " voxels per resel" << endl;
  }*/ //aranyics

  if(verbose)
  {
  cout << "DLH " << dLh << endl;
  cout << "VOLUME " << mask_volume << endl;
  cout << "RESELS " << resels << endl;
  cout << "FWHMvoxel " << FWHM[X] << " " <<  FWHM[Y] << " " << FWHM[Z] << endl;
  cout << "FWHMmm " << FWHMmm[X] <<  " " << FWHMmm[Y] << " " << FWHMmm[Z] << endl;
  cout << "sigmasq " << sigmasq[X] << " " << sigmasq[Y] << " " << sigmasq[Z] << endl; //aranyics
  }

  return EXIT_SUCCESS;
}

//g++ -o smoothest_ext.o -c smoothest_ext.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11
