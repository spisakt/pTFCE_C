#ifndef CLUSTERSIZE_H
#define CLUSTERSIZE_H

#include "newimage/newimageall.h"

namespace NEWIMAGE{

template <class T>
class volumeClust: public volume<T>
{
  public:
    void projectClusterValues(ColumnVector &clustervalues, T background);
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


}



#endif //CLUSTERSIZE_H
