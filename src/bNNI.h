#ifndef BNNI_H_
#define BNNI_H_

#include "bme.h"
#include "heap.h"

void bNNIRetestEdge(int *p, int *q, edge *e,tree *T, double **avgDistArray, 
	double *weights, int *location, int *possibleSwaps);
void bNNItopSwitch(tree *T, edge *e, int direction, double **A);
void bNNI(tree *T, double **avgDistArray, int *count, FILE *statfile);
void updateSubTreeAfterNNI(double **A, node *v, edge *rootEdge, node *closer, node *further,
	double dcoeff, int direction);
void bNNIupdateAverages(double **A, node *v, edge *par, edge *skew, 
	edge *swap, edge *fixed);
double wf5(double D_AD, double D_BC, double D_AC, double D_BD,
	double D_AB, double D_CD);
int bNNIEdgeTest(edge *e, tree *T, double **A, double *weight);
void limitedFillTableUp(edge *e, edge *f, double **A, edge *trigger);

#endif /*BNNI_H_*/

