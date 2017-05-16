#ifndef MVR_H_
#define MVR_H_

#include "BIONJ.h"

tree *unj (double **D, set *species, int n);
int SymmetrizeMVR(double **delta, int n);
double Finish_branch_length_MVR(int i, int j, int k, double **delta, int n);
void FinishStrMVR (double **delta, int n, POINTERS *trees, char *StrTree);
void Branch_lengthMVR(int a, int b, double *la, double *lb, double **delta, int n);
double Reduction10MVR(int a, double la, int b, double lb, int i, double lamda, double **delta);
double Reduction11MVR(int x, int y, int i, double **delta);
double LamdaMVR(int x, int y, int i, double **delta);
double WeightMVR(int x, int y, int i, double **delta, double u);
double mu (int a, int b, double **delta, int n);

#endif /*MVR_H_*/

