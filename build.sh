#!/bin/bash

#echo "g++ -o grfClust.o -c grfClust.cc -std=c++11"
g++ -o grfClust.o -c grfClust.cc -std=c++11

#echo "g++ -o pvalutil.o -c pvalutil.cc -I$FSLDIR/extras/include -std=c++11 -lm"
g++ -o pvalutil.o -c pvalutil.cc -I$FSLDIR/extras/include -std=c++11 -lm

#echo "g++ -o smoothest_ext.o -c smoothest_ext.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11"
g++ -o smoothest_ext.o -c smoothest_ext.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11

#echo "g++ -o smoothest smoothest_ext.o smoothest_app.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack"
g++ -o smoothest smoothest_ext.o smoothest_app.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack

#echo "g++ -o connectedcompsize connectedcompsize.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/include/ -I$FSLDIR/extras/include/armawrap/armawrap/ -std=c++11 -L$FSLDIR/lib/ -L$FSLDIR/extras/lib/ -lm -lutils -lnewimage -lmiscmaths -lNewNifti -lznz -lz -lblas -llapack"
g++ -o connectedcompsize connectedcompsize.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/include/ -I$FSLDIR/extras/include/armawrap/armawrap/ -std=c++11 -L$FSLDIR/lib/ -L$FSLDIR/extras/lib/ -lm -lutils -lnewimage -lmiscmaths -lNewNifti -lznz -lz -lblas -llapack

#echo "g++ -o ptfce.o -c ptfce.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11"
g++ -o ptfce.o -c ptfce.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/ -std=c++11

#echo "g++ -o pTFCE-test smoothest_ext.o pvalutil.o grfClust.o ptfce.o pTFCE-test.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/  -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack -lgsl"
g++ -o pTFCE-test smoothest_ext.o pvalutil.o grfClust.o ptfce.o pTFCE-test.cc -D_GLIBCXX_USE_CXX11_ABI=0 -I$FSLDIR/extras/include/armawrap/armawrap/ -I$FSLDIR/extras/include/ -I$FSLDIR/include/  -L$FSLDIR/extras/lib/ -L$FSLDIR/lib/ -std=c++11 -lm -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lNewNifti -lznz -lnewmat -lz -lblas -llapack -lgsl -lgslcblas
