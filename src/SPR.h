#ifndef SPR_H_
#define SPR_H_

#include "bme.h"
#include "traverse.h"
#include "newick.h"
#include "inputs.h"

void zero3DMatrix(double ***X, int h, int l, int w);
void findTableMin(int *imin, int *jmin, int *kmin, int n, double ***X, double *min);
void SPR(tree *T, double **D, double **A, int *count, FILE *statfile);
void assignSPRWeights(node *vtest, double **A, double ***swapWeights);
void assignDownWeightsUp(edge *etest, node *vtest, node *va, edge *back, node *cprev,
	double oldD_AB, double coeff, double **A, double ***swapWeights);
void assignDownWeightsSkew(edge *etest, node *vtest, node *va, edge *back, node *cprev,
	double oldD_AB, double coeff, double **A, double ***swapWeights);
void assignDownWeightsDown(edge *etest, node *vtest, node *va, edge *back, node *cprev,
	double oldD_AB, double coeff, double **A, double ***swapWeights);
void assignUpWeights(edge *etest, node *vtest, node *va, edge *back, node *cprev,
	double oldD_AB, double coeff, double **A, double ***swapWeights);
void pruneSubtree(edge *p, edge *u, edge *d);
void SPRsplitEdge(edge *e, edge *p, edge *d);
void SPRDownShift(tree *T, node *v, edge *e);
void SPRUpShift(tree *T, node *vmove, edge *esplit);
void SPRTopShift(tree *T, node *vmove, edge *esplit, int UpOrDown);

#endif /*SPR_H_*/

