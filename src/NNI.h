#ifndef NNI_H_
#define NNI_H_

#include "gme.h"
#include "heap.h"

double **buildAveragesTable(tree *T, double **D);
double wf2(double lambda, double D_AD, double D_BC, double D_AC, double D_BD,
	double D_AB, double D_CD);
int NNIEdgeTest(edge *e, tree *T, double **A, double *weight);
void NNIupdateAverages(double **A, edge *e, edge *par, edge *skew, 
	edge *swap, edge *fixed, tree *T);
void NNItopSwitch(tree *T, edge *e, int direction, double **A);
void NNIRetestEdge(int *p, int *q, edge *e,tree *T, double **avgDistArray, 
	double *weights, int *location, int *possibleSwaps);
void NNI(tree *T, double **avgDistArray, int *count, FILE *statfile);
/*
void NNIwithoutMatrix(tree *T, double **D, int *count);
void NNIWithPartialMatrix(tree *T,double **D,double **A,int *count);
*/

#endif /*NNI_H_*/

