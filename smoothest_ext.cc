#include <iostream>
#include <cmath>
#include <string>
#include <map>

#include "mathutil.h"
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

    for ( unsigned short t = 0; t < M; t++ ) {   // aranyics
	if( isnan2(R(x,y,z,t)) ) R(x,y,z,t) = 0.0;
    }

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


int smoothestVox(double &dLh, unsigned long &mask_volume, double &resels, double *FWHM, double *FWHMmm, double *sigmasq, NEWIMAGE::volume<float>& RPV, NEWIMAGE::volume<float>& FWHMimg,
              NEWIMAGE::volume4D<float>& R,
              NEWIMAGE::volume<float>& mask,
              double dof, bool verbose)
{
  if(verbose) cerr << "Standardising....";      //aranyics @ every case of verbose.value()
  mask_volume = standardise(mask, R);           //aranyics
  if(verbose) print_volume_info(R, "standardized");
  if(verbose) cerr << "done" << endl;

  if(verbose) cerr << "Masked-in voxels = " << mask_volume << endl;

  unsigned long N = 0;
  unsigned int D = 3;

  // MJ additions to make it cope with 2D images
  bool usez = true;
  if (R.zsize() <= 1) { usez = false; D = 2; }
  if ((!usez) && verbose) {
    cout << "Using 2D image mode." << endl;
  }

  // Creating volumes for RPV - aranyics
  RPV.reinitialize(mask, false); //FSL < 6.0.2
  //RPV.reinitialize(mask, TEMPLATE); // FSL <= 6.0.2
  RPV *= 0;
  copyconvert(RPV, FWHMimg);
  NEWIMAGE::volume<float> SS_X, SS_Y, SS_Z;
  NEWIMAGE::volume<float> S2X, S2Y, S2Z;
  NEWIMAGE::volume<float> SSQX, SSQY, SSQZ;
  NEWIMAGE::volume<float> FWHMX, FWHMY, FWHMZ;
  NEWIMAGE::volume<float> RPVX, RPVY, RPVZ;
  copyconvert(RPV,SS_X);
  copyconvert(RPV,SS_Y);
  copyconvert(RPV,SS_Z);
  copyconvert(RPV,S2X);
  copyconvert(RPV,S2Y);
  copyconvert(RPV,S2Z);
  copyconvert(RPV,SSQX);
  copyconvert(RPV,SSQY);
  copyconvert(RPV,SSQZ);
  copyconvert(RPV,FWHMX);
  copyconvert(RPV,FWHMY);
  copyconvert(RPV,FWHMZ);
  copyconvert(RPV,RPVX);
  copyconvert(RPV,RPVY);
  copyconvert(RPV,RPVZ);


  // Estimate the smoothness of the normalised residual field
  // see TR00DF1 for mathematical description of the algorithm.
  enum {X = 0, Y, Z};
  double SSminus[3] = {0, 0, 0}, S2[3] = {0, 0, 0};

  int zstart=1;
  if (!usez) zstart=0;
  for ( unsigned short z = zstart; z < R.zsize() ; z++ )
    for ( unsigned short y = 1; y < R.ysize() ; y++ )
      for ( unsigned short x = 1; x < R.xsize() ; x++ ) {
	// Sum over N
	// Losing 1 voxel on far edges on the whole volume resels computation
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

	if( (mask(x, y, z)>0.5 ) &&
	    (mask(x-1, y, z)>0.5) &&
	    (mask(x, y-1, z)>0.5) &&
	    ( (!usez) || (mask(x, y, z-1)>0.5) ) ) {
	  // TODO RPV image is translated by 1 voxel
	  for ( unsigned short t = 0; t < R.tsize(); t++ ) {
	    // Sum per voxels for RPV - aranyics
	    if (mask(x-1, y, z)>0.5) SS_X(x, y, z) += R(x, y, z, t) * R(x-1, y, z, t);
	    //if (mask(x+1, y, z)>0.5) SS_X(x, y, z) += R(x, y, z, t) * R(x+1, y, z, t);
	    if (mask(x, y-1, z)>0.5) SS_Y(x, y, z) += R(x, y, z, t) * R(x, y-1, z, t);
	    //if (mask(x, y+1, z)>0.5) SS_Y(x, y, z) += R(x, y, z, t) * R(x, y+1, z, t);
	    if (usez) {
	      if (mask(x, y, z-1)>0.5) SS_Z(x, y, z) += R(x, y, z, t) * R(x, y, z-1, t);
	      //if (mask(x, y, z+1)>0.5) SS_Z(x, y, z) += R(x, y, z, t) * R(x, y, z+1, t);
	    }
	    if (mask(x-1, y, z)>0.5) S2X(x, y, z) += 0.5 * (Sqr(R(x, y, z, t)) + Sqr(R(x-1, y, z, t)));
	    //if (mask(x+1, y, z)>0.5) S2X(x, y, z) += 0.5 * (Sqr(R(x, y, z, t)) + Sqr(R(x+1, y, z, t)));
	    if (mask(x, y-1, z)>0.5) S2Y(x, y, z) += 0.5 * (Sqr(R(x, y, z, t)) + Sqr(R(x, y-1, z, t)));
	    //if (mask(x, y+1, z)>0.5) S2Y(x, y, z) += 0.5 * (Sqr(R(x, y, z, t)) + Sqr(R(x, y+1, z, t)));
	    if (usez) {
	      if (mask(x, y, z-1)>0.5) S2Z(x, y, z) += 0.5 * (Sqr(R(x, y, z, t)) + Sqr(R(x, y, z-1, t)));
	      //if (mask(x, y, z+1)>0.5) S2Z(x, y, z) += 0.5 * (Sqr(R(x, y, z, t)) + Sqr(R(x, y, z+1, t)));
	    }
	  }
	  // for extreme smoothness
	  if (SS_X(x, y, z) >= 0.9999*S2X(x, y, z)) SS_X(x, y, z) = 0.99*S2X(x, y, z);
	  if (SS_Y(x, y, z) >= 0.9999*S2Y(x, y, z)) SS_Y(x, y, z) = 0.99*S2Y(x, y, z);
	  if (SS_Z(x, y, z) >= 0.9999*S2Z(x, y, z)) SS_Z(x, y, z) = 0.99*S2Z(x, y, z);
	}
      }

  if (false && verbose) {
    save_volume(RPV, "testdata/RPV_test.nii.gz");
    save_volume(SS_X, "testdata/SS_X_test.nii.gz");
    save_volume(SS_Y, "testdata/SS_Y_test.nii.gz");
    save_volume(SS_Z, "testdata/SS_Z_test.nii.gz");
    save_volume(S2X, "testdata/S2X_test.nii.gz");
    save_volume(S2Y, "testdata/S2Y_test.nii.gz");
    save_volume(S2Z, "testdata/S2Z_test.nii.gz");
  }

  double norm = 1.0/(double) N;
  double v = dof;       // v - degrees of freedom (nu) //aranyics
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

  if(false && verbose) {
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


  // Convert to sigma squared to FWHM for RPV - aranyics
  zstart=1;
  if (!usez) zstart=0;
  for ( unsigned short z = zstart; z < R.zsize() ; z++ )
    for ( unsigned short y = 1; y < R.ysize() ; y++ )
      for ( unsigned short x = 1; x < R.xsize() ; x++ )
      {
	if( (mask(x, y, z)>0.5) &&
	    (mask(x-1, y, z)>0.5) &&
	    (mask(x, y-1, z)>0.5) &&
	    ( (!usez) || (mask(x, y, z-1)>0.5) ) ) {

	    SSQX(x, y, z) = -1.0 / (4 * log(fabs(SS_X(x, y, z) / S2X(x, y, z))));
	    SSQY(x, y, z) = -1.0 / (4 * log(fabs(SS_Y(x, y, z) / S2Y(x, y, z))));
	    if (usez) SSQZ(x, y, z) = -1.0 / (4 * log(fabs(SS_Z(x, y, z) / S2Z(x, y, z))));

	    //TODO fwhm
	    /*FWHMX(x, y, z) = sqrt(8 * log(2) * SSQX(x, y, z));
	    FWHMY(x, y, z) = sqrt(8 * log(2) * SSQY(x, y, z));
	    if (usez) FWHMZ(x, y, z) = sqrt(8 * log(2) * SSQZ(x, y, z));

	    RPV(x, y, z) = (FWHMX(x, y, z)*R.xdim()) * (FWHMY(x, y, z)*R.ydim());
	    if (usez) RPV(x, y, z) = RPV(x, y, z) * (FWHMZ(x, y, z)*R.zdim());
	    //FWHMimg(x, y, z) = RPV(x, y, z);
	    RPV(x, y, z) = (R.xdim()*R.ydim()*R.zdim()) / RPV(x, y, z);
	    FWHMimg(x, y, z) = usez ? pow(RPV(x, y, z), -0.5) : pow(RPV(x, y, z), -0.33333333);*/

	    //TODO nan values at mask edges
	    RPVX(x, y, z) = pow( 4 * log( fabs(S2X(x, y, z) / SS_X(x, y, z)) ) , 0.5 ) * pow(8 * log(2), -0.5);
	    RPVY(x, y, z) = pow( 4 * log( fabs(S2Y(x, y, z) / SS_Y(x, y, z)) ) , 0.5 ) * pow(8 * log(2), -0.5);
	    RPVZ(x, y, z) = pow( 4 * log( fabs(S2Z(x, y, z) / SS_Z(x, y, z)) ) , 0.5 ) * pow(8 * log(2), -0.5);

	    RPV(x, y, z) = RPVX(x, y, z) * RPVY(x, y, z) * RPVZ(x, y, z);
	}
	    //if ( isnan2(RPV(x, y, z)) ) RPV(x, y, z) = 0.0;

	    FWHMimg(x, y, z) = pow(RPV(x, y, z), -(0.5));
      }

  if (verbose) {
    //save_volume(SSQX, "testdata/SSQX_test.nii.gz");
    //save_volume(SSQY, "testdata/SSQY_test.nii.gz");
    //save_volume(SSQZ, "testdata/SSQZ_test.nii.gz");
    //save_volume(FWHMX, "testdata/FWHMX_test.nii.gz");
    //save_volume(FWHMY, "testdata/FWHMY_test.nii.gz");
    //save_volume(FWHMZ, "testdata/FWHMZ_test.nii.gz");

    cout << "Residual " << R(R.xsize()/2-1, R.ysize()/2, R.zsize()/2, 0) << " "
                        << R(R.xsize()/2,   R.ysize()/2, R.zsize()/2, 0) << " "
                        << R(R.xsize()/2+1, R.ysize()/2, R.zsize()/2, 0) << endl;
    cout << "SS_X " << SS_X(R.xsize()/2, R.ysize()/2, R.zsize()/2) << endl;
    cout << "S2X " << S2X(R.xsize()/2, R.ysize()/2, R.zsize()/2) << endl;
    cout << "SSQX " << SSQX(R.xsize()/2, R.ysize()/2, R.zsize()/2) << endl;
    cout << "FWHMX " << FWHMX(R.xsize()/2, R.ysize()/2, R.zsize()/2) << endl;
    cout << "FWHMY " << FWHMY(R.xsize()/2, R.ysize()/2, R.zsize()/2) << endl;
    cout << "FWHMZ " << FWHMZ(R.xsize()/2, R.ysize()/2, R.zsize()/2) << endl;
    cout << "FWHM " << FWHMimg(R.xsize()/2, R.ysize()/2, R.zsize()/2) << endl;
    cout << "RPV " << RPV(R.xsize()/2, R.ysize()/2, R.zsize()/2) << endl;
    cout << "dims " << R.xdim() << " " << R.ydim() << " " << R.zdim() << endl;
  }


  if(false && verbose)
  {
    cout << "FWHMx = " << FWHM[X] << " mm, "
         << "FWHMy = " << FWHM[Y] << " mm";
    if (usez) cout << ", FWHMz = " << FWHM[Z] << " mm";
    cout << endl;
    cout << "DLH " << dLh << " voxels^-3" << endl;
    cout << "VOLUME " << mask_volume << " voxels" << endl;
    cout << "RESELS " << resels << " voxels per resel" << endl;
  } //aranyics

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



int estimateRPV(NEWIMAGE::volume<float>& RPV, NEWIMAGE::volume<float>& FWHMimg,
              NEWIMAGE::volume4D<float>& R,
              NEWIMAGE::volume<float>& mask,
              double dof, bool verbose)
{
  if(verbose) cerr << "Standardising....";      //aranyics @ every case of verbose.value()
  unsigned long mask_volume = standardise(mask, R);           //aranyics
  if(verbose) cerr << "done" << endl;

  //if(verbose) cerr << "Masked-in voxels = " << mask_volume << endl;

  double n = dof; //TODO
  unsigned short D = 3; // TODO
  if (D == 0) exit(EXIT_SUCCESS); //TODO
  unsigned long N = 0;

  enum {X = 0, Y, Z};
  double SSminus[3] = {0, 0, 0}, S2[3] = {0, 0, 0};

  // Creating volumes for RPV - aranyics
  //RPV.reinitialize(mask, false); // FSL < 6.0.2
  RPV.reinitialize(mask, TEMPLATE); // FSL >= 6.0.2
  RPV *= 0;
  copyconvert(RPV, FWHMimg);
  NEWIMAGE::volume<float> SSQ;
  NEWIMAGE::volume<float> Lxx, Lxy, Lyy, Lxz, Lyz, Lzz;
  copyconvert(RPV,SSQ);
  copyconvert(RPV,Lxx);
  copyconvert(RPV,Lxy);
  copyconvert(RPV,Lyy);
  copyconvert(RPV,Lxz);
  copyconvert(RPV,Lyz);
  copyconvert(RPV,Lzz);


  bool usez = true;
  int zstart=1;
  if (!usez) zstart=0;
  for ( unsigned short z = zstart; z < R.zsize() ; z++ )
    for ( unsigned short y = 1; y < R.ysize() ; y++ )
      for ( unsigned short x = 1; x < R.xsize() ; x++ )
      {
	// Sum over N
	N++;

	for ( unsigned short t = 0; t < R.tsize() ; t++ )
	    if( (mask(x, y, z)>0.5) ) {
		double dx = (x+1 < R.xsize()) ? R(x, y, z, t) - R(x+1, y, z, t) : 0.0;
		double dy = (y+1 < R.ysize()) ? R(x, y, z, t) - R(x, y+1, z, t) : 0.0;
		double dz = (z+1 < R.zsize()) ? R(x, y, z, t) - R(x, y, z+1, t) : 0.0;

		SSQ(x, 1, 1) += R(x, y, z) * R(x, y, z);

		if ( D >= 1 ) { Lxx(x, y, z) += dx * dx; };
		if ( D >= 2 ) { Lxy(x, y, z) += dx * dy; Lyy(x, y, z) += dy * dy; };
		if ( D >= 3 ) { Lxz(x, y, z) += dx * dz; Lyz(x, y, z) += dy * dz; Lzz(x, y, z) += dz * dz; };
	    }
      }

  double q = pow( (4 * log(2)), D );
  for ( unsigned short z = zstart; z < R.zsize() ; z++ )
    for ( unsigned short y = 1; y < R.ysize() ; y++ )
      for ( unsigned short x = 1; x < R.xsize() ; x++ )
      {
	if (mask(x, y, z) > 0.5)
	{
	  Lxx(x, y, z) = (Lxx(x, y, z) / R.tsize()) * (n/dof);
	  Lxy(x, y, z) = (Lxy(x, y, z) / R.tsize()) * (n/dof);
	  Lyy(x, y, z) = (Lyy(x, y, z) / R.tsize()) * (n/dof);
	  Lxz(x, y, z) = (Lxz(x, y, z) / R.tsize()) * (n/dof);
	  Lyz(x, y, z) = (Lyz(x, y, z) / R.tsize()) * (n/dof);
	  Lzz(x, y, z) = (Lzz(x, y, z) / R.tsize()) * (n/dof);

	  double reselVox = Lxx(x, y, z) * Lyy(x, y, z) * Lzz(x, y, z)     +
                            Lxy(x, y, z) * Lyz(x, y, z) * Lxz(x, y, z) * 2 -
                            Lxx(x, y, z) * Lyz(x, y, z) * Lyz(x, y, z)     -
                            Lxy(x, y, z) * Lxy(x, y, z) * Lzz(x, y, z)     -
                            Lxz(x, y, z) * Lyy(x, y, z) * Lxz(x, y, z);

	  if (reselVox < 0.0) reselVox = 0.0;

	  //RPV(x, y, z) = reselVox;
	  RPV(x, y, z) = sqrt( reselVox / q );
	}
	  if ( isnan2(RPV(x, y, z)) ) RPV(x, y, z) = 0.0;
	  FWHMimg(x, y, z) = pow(RPV(x, y, z), -(0.5));
      }

  if(false){
    int kk = 100;
  for ( unsigned short z = zstart; z < R.zsize() ; z++ )
    for ( unsigned short y = 1; y < R.ysize() ; y++ )
      for ( unsigned short x = 1; x < R.xsize() ; x++ )
      {
	if (mask(x,y,z) > 0.5 && kk > 0) {
	    kk--;
	    //cout << Lxx(x, y, z) << " "  << RPV(x, y, z) << endl;
	}
      }
  }


  return EXIT_SUCCESS;
}

//g++ -o smoothest_ext.o -c smoothest_ext.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11
