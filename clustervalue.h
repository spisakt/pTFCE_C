#ifndef CLUSTERSIZE_H
#define CLUSTERSIZE_H

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"

#include "mathutil.h"

namespace NEWIMAGE{

template <class T>
class volumeClust: public volume<T>
{
  public:
    void projectClusterValues(ColumnVector &clustervalues, T background);
    void sumClusterRPVValues(volume<T> &RPV, ColumnVector &sums);
};


template <class T>
void volumeClust<T>::projectClusterValues(ColumnVector &clustervalues, T background)
{
    for (typename volume<T>::nonsafe_fast_iterator it=this->nsfbegin(), itend=this->nsfend(); it != itend; ++it)
    {
	if ( (*it) != 0.0 )
	{
	    *it = clustervalues((*it));
	}
	else
	{
	    *it = background;
	}
    }
}


template <class T>
void volumeClust<T>::sumClusterRPVValues(volume<T> &RPV, ColumnVector &sums)
{
    typedef typename volume<T>::nonsafe_fast_iterator volumeIter;

    sums.ReSize( MISCMATHS::round(this->max()) );
    sums = 0;

    ColumnVector allVox, nonanVox;
    allVox.ReSize( MISCMATHS::round(this->max()) );
    nonanVox.ReSize( MISCMATHS::round(this->max()) );
    allVox = 0; nonanVox = 0;

    //TODO validate summation correctness

    for (std::pair<volumeIter, volumeIter> i(RPV.nsfbegin(), this->nsfbegin()); i.second != this->nsfend(); ++i.first, ++i.second)
    {
        if ( (*i.second) != 0 )
        {
            allVox((*i.second))++;
            if ( !isnan2((*i.first)) ) {
                nonanVox((*i.second))++;
                sums((*i.second)) += *i.first;
            }
        }
    }

    for( int i = 1; i <= MISCMATHS::round(this->max()); ++i )
    {
//cout << "*SUM RPV*************************************** " << sums(i) << "correction factor for nan" << allVox(i)/nonanVox(i) << endl;
        sums(i) *= allVox(i) / nonanVox(i);
    }

    this->projectClusterValues(sums, 0.0f);
}


}



#endif //CLUSTERSIZE_H
