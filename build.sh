#!/bin/bash

OPT="-O3"

echo "g++ -o grfClust.o -c grfClust.cc -std=c++11 $LINK"
g++ -o grfClust.o -c grfClust.cc -std=c++11 $LINK

echo "g++ -o pvalutil.o -c pvalutil.cc -I$FSLDIR/extras/include -std=c++11 -lm $LINK"
g++ -o pvalutil.o -c pvalutil.cc -I$FSLDIR/extras/include -std=c++11 -lm $LINK

echo "g++ -o smoothest_ext.o -c smoothest_ext.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11 $LINK"
g++ -o smoothest_ext.o -c smoothest_ext.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11 $LINK

echo "g++ -o smoothest smoothest_ext.o smoothest_app.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack $LINK"
g++ -o smoothest smoothest_ext.o smoothest_app.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack $LINK

echo "g++ -o connectedcompsize connectedcompsize.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/include/ -I$FSLDIR/extras/include/armawrap/armawrap/ -std=c++11 -L$FSLDIR/lib/ -L$FSLDIR/extras/lib/ -lm -lutils -lnewimage -lmiscmaths -lNewNifti -lznz -lz -lblas -llapack $LINK"
g++ -o connectedcompsize connectedcompsize.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/include/ -I$FSLDIR/extras/include/armawrap/armawrap/ -std=c++11 -L$FSLDIR/lib/ -L$FSLDIR/extras/lib/ -lm -lutils -lnewimage -lmiscmaths -lNewNifti -lznz -lz -lblas -llapack $LINK

echo "g++ -o ptfce.o -c ptfce.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11 $LINK"
g++ -o ptfce.o -c ptfce.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11 $LINK

echo "g++ -o pTFCE-test smoothest_ext.o pvalutil.o grfClust.o ptfce.o pTFCE-test.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/  -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack -lgsl $LINK"
g++ -o pTFCE-test smoothest_ext.o pvalutil.o grfClust.o ptfce.o pTFCE-test.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/  -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack -lgsl -lgslcblas $LINK

echo "g++ -o ptfce smoothest_ext.o pvalutil.o grfClust.o ptfce.o pTFCE_cmd.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/  -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack -lgsl $LINK"
g++ -o ptfce smoothest_ext.o pvalutil.o grfClust.o ptfce.o pTFCE_cmd.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/  -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack -lgsl -lgslcblas $LINK
