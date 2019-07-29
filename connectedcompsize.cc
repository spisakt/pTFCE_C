#include <iostream>
#include <string>
#include "clustervalue.h"
#include "newimage/newimageall.h"

using namespace NEWIMAGE;

int main(int argc, char* argv[])
{
  volume<float> invol;
  volume<int> outvol;
  volumeClust<int> outsizevol;
  NEWMAT::ColumnVector clustersizes;

  if (argc<2) {
    cerr << "Usage: " << argv[0] << " <in_volume> [outputvol [num_connect]]" << endl;
    return -1;
  }

  string inname = argv[1], outname, outsizename;
  if (argc>2) {
    outname = argv[2];
  } else {
    outname = inname + "_label";
    outsizename = inname + "_size";
  }

  int num_connect=26;
  if (argc>3) {
    num_connect=atoi(argv[3]);
    if ((num_connect!=6) && (num_connect!=18) && (num_connect!=26)) {
      cerr << "Can only have num_connect equal 4, 18 or 26" << endl;
      return 1;
    }
  }

  read_volume(invol, inname);
  outvol = connected_components(invol, clustersizes, num_connect);
  save_volume(outvol, outname);
  copyconvert(outvol, outsizevol);
  outsizevol.projectClusterValues(clustersizes, 0.0f);
  save_volume(outsizevol, outsizename);

  std::cout << "Clus\tVoxels" << std::endl;
  for (int i = 1; i <= clustersizes.n_rows; ++i)
  {
    std::cout << i << "\t" << clustersizes(i) << std::endl;
  }

  return 0;
}

//g++ -o connectedcompsize connectedcompsize.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I/site/fsl-6.0/include/ -I/site/fsl-6.0/extras/include/armawrap/armawrap/ -std=c++11
// -L/site/fsl-6.0/lib/ -L/site/fsl-6.0/extras/lib/ -lm -lutils -lnewimage -lmiscmaths -lNewNifti -lznz -lz -lblas -llapack
