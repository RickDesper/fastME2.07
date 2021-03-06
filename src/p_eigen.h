#ifndef P_EIGEN_H
#define P_EIGEN_H

#include "p_utils.h"

#define BASE        2    // base of floating point arithmetic
#define DIGITS     40    // no. of digits to the base BASE in the fraction
#define MAXITER    30    // max2. no. of iterations to converge

#define pos(i,j,n)      ((i)*(n)+(j))


int Eigen(int job, double *A, int n, double *rr, double *ri,
          double *vr, double *vi, double *w);
void balance(double *mat, int n, int *low, int *hi, double *scale);
void unbalance(int n, double *vr, double *vi, int low, int hi,
               double *scale);
int realeig(int job, double *mat, int n,int low, int hi, double *valr,
            double *vali, double *vr, double *vi);
void elemhess(int job, double *mat, int n, int low, int hi, 
            double *vr, double *vi, int *work);
int ludcmp(double **a, int n, double *d);
//void det(double **a, int n, double *d);

/* complex functions */

typedef struct { double re, im; } complex;
#define csize(a) (fabs(a.re)+fabs(a.im))

complex compl (double re,double im);
complex _conj (complex a);
complex cplus (complex a, complex b);
complex cminus (complex a, complex b);
complex cby (complex a, complex b);
complex cdiv (complex a,complex b);
/* complex local_cexp (complex a); */
complex cfactor (complex x, double a);
int cxtoy (complex *x, complex *y, int n);
int cmatby (complex *a, complex *b, complex *c, int n,int m,int k);
int cmatout (FILE * fout, complex *x, int n, int m);
int cmatinv( complex *x, int n, int m, double *space);


#endif /*P_EIGEN_H_*/

