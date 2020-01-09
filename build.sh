#!/bin/bash

# for release:
OPT="-Ofast -D_GLIBCXX_PARALLEL"
PIE="-no-pie"

echo "g++ -o mathutil.o -c mathutil.cc $PIE -std=c++11 -lm $OPT"
g++ -o mathutil.o -c mathutil.cc $PIE -std=c++11 -lm $OPT

echo "g++ -o pvalutil.o -c pvalutil.cc -I$FSLDIR/extras/include -std=c++11 -lm $OPT"
g++ -o pvalutil.o -c pvalutil.cc $PIE -I$FSLDIR/extras/include -std=c++11 -lm $OPT

echo "g++ -o grfClust.o -c grfClust.cc -std=c++11 $OPT"
g++ -o grfClust.o -c grfClust.cc $PIE -std=c++11 $OPT

echo "g++ -o smoothest_ext.o -c smoothest_ext.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11 $OPT"
g++ -o smoothest_ext.o -c smoothest_ext.cc $PIE -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11 $OPT

echo "g++ -o smoothest smoothest_ext.o smoothest_app.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack $OPT"
g++ -o smoothest mathutil.o smoothest_ext.o smoothest_app.cc $PIE -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack $OPT

echo "g++ -o connectedcompsize connectedcompsize.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/include/ -I$FSLDIR/extras/include/armawrap/armawrap/ -std=c++11 -L$FSLDIR/lib/ -L$FSLDIR/extras/lib/ -lm -lutils -lnewimage -lmiscmaths -lNewNifti -lznz -lz -lblas -llapack $OPT"
g++ -o connectedcompsize connectedcompsize.cc $PIE -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/include/ -I$FSLDIR/extras/include/armawrap/armawrap/ -std=c++11 -L$FSLDIR/lib/ -L$FSLDIR/extras/lib/ -lm -lutils -lnewimage -lmiscmaths -lNewNifti -lznz -lz -lblas -llapack $OPT

echo "g++ -o ptfce.o -c ptfce.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11 $OPT"
g++ -o ptfce.o -c ptfce.cc $PIE -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11 $OPT

echo "g++ -o pTFCE-test smoothest_ext.o pvalutil.o grfClust.o ptfce.o pTFCE-test.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/  -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack -lgsl $OPT"
g++ -o pTFCE-test mathutil.o smoothest_ext.o pvalutil.o grfClust.o ptfce.o pTFCE-test.cc $PIE -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/  -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack -lgsl -lgslcblas $OPT

echo "g++ -o ptfce smoothest_ext.o pvalutil.o grfClust.o ptfce.o ptfce_cmd.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/  -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack -lgsl $OPT"
g++ -o pTFCE mathutil.o smoothest_ext.o pvalutil.o grfClust.o ptfce.o ptfce_cmd.cc $PIE -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/  -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack -lgsl -lgslcblas $OPT
