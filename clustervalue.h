#ifndef CLUSTERSIZE_H
#define CLUSTERSIZE_H

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"

namespace NEWIMAGE{

template <class T>
class volumeClust: public volume<T>
{
  public:
    void projectClusterValues(ColumnVector &clustervalues, T background);
    void sumClusterValues(volume<T> &RPV, ColumnVector &sums);
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
void volumeClust<T>::sumClusterValues(volume<T> &RPV, ColumnVector &sums)
{
    typedef typename volume<T>::nonsafe_fast_iterator volumeIter;

    sums.ReSize( MISCMATHS::round(this->max()) );
    sums = 0;

    for (std::pair<volumeIter, volumeIter> i(RPV.nsfbegin(), this->nsfbegin()); i.second != this->nsfend(); ++i.first, ++i.second)
    {
        if ( (*i.second) != 0 )
        {
            sums((*i.second)) += *i.first;
        }
    }

    this->projectClusterValues(sums, 0.0f);
}


}



#endif //CLUSTERSIZE_H
