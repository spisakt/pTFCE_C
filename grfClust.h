#ifndef GRFCLUST_H
#define GRFCLUST_H


double dvox(double h);
double pvox(double h);
double Es(double h, unsigned int V, double Rd);
double dcl(double h, void * p);
double dclust(double h, void * p);
double dvox_dclust(double h, void * p);
double dvox_clust(double h, void * p);
double pvox_clust(double actH, void * p);


struct dcl_params{
    unsigned long V;
    double Rd;
    double c;
    double ZestThr;
};

#endif //GRFCLUST_H
