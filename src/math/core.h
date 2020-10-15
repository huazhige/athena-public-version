#ifndef CORE_H_
#define CORE_H_
#include <math.h>

inline double sqr(double x) { 
  return x*x; 
}
inline double cub(double x) {
  return x*x*x;
}
inline double allmax(double *a, int n) {
  double v = a[0];
  for (int i = 1; i < n; ++i)
    if (v < a[i]) v = a[i];
  return v;
}
inline double allmin(double *a, int n) {
  double v = a[0];
  for (int i = 1; i < n; ++i)
    if (v > a[i]) v = a[i];
  return v;
}
inline int sign(double x) {
  return x < 0. ? -1 : 1;
}
inline double rad2deg(double phi) {
  return phi*180./M_PI;
}
inline double deg2rad(double phi) {
  return phi*M_PI/180.;
}

#endif
