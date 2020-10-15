#ifndef LINALG_H_
#define LINALG_H_

double vvdot(double const *a, double const *b, int n);
void mvdot(double *r, double **m, double const *v, int n1, int n2);
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double *b);
void leastsq(double **A, double *b, int n1, int n2);

#endif
