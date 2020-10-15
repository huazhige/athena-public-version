#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

int locate(double const *xx, double x, int n);
void interpn(double *val, double const *coor, double const *data, double const *axis, 
  int const *len, int ndim, int nval = 1);
double interp1(double x, double const *data, double const *axis, int len);

#endif
