#ifndef PVALUTIL_H
#define PVALUTIL_H

double fwerp2z(double nresels, double p=0.05, bool twotailed=false, bool grf=true);

double fwerz2p(double nresels, double z, bool twotailed=false, bool grf=true);

#endif //PVALUTIL_H
