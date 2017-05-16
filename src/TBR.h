#ifndef TBR_H_
#define TBR_H_

#include "SPR.h"
	
void assignTBRUpWeights(edge *ebottom, node *vtest, node *va, edge *back, node *cprev,
	double oldD_AB, double coeff, double **A, double **dXTop, double ***swapWeights,
	edge *etop, double *bestWeight, edge **bestSplit, edge **bestTop, edge **bestBottom);
void assignTBRDownWeightsUp(edge *etest, node *vtest, node *va, edge *back, node *cprev,
	double oldD_AB, double coeff, double **A, double ***swapWeights, double *bestWeight,
	edge **bestSplitEdge, edge **bestTop, edge **bestBottom);
void assignTBRDownWeightsSkew(edge *etest, node *vtest, node *va, edge *back, node *cprev,
	double oldD_AB, double coeff, double **A, double ***swapWeights, double *bestWeight,
	edge **bestSplitEdge, edge **bestTop, edge **bestBottom);
void assignTBRDownWeightsDown(edge *etest, node *vtest, node *va, edge *back, node *cprev,
	double oldD_AB, double coeff, double **A, double ***swapWeights, double *bestWeight,
	edge **bestSplitEdge, edge **bestTop, edge **bestBottom);
void calcTBRTopBottomAverage(node *vbottom, double **A, double **dXTop, double dXOut,
	edge *esplit, edge *etop, edge *eback, int UpOrDown);
void calcTBRaverages(tree *T, edge *esplit, double **A, double **dXTop);
void TBR(tree *T, double **A, double **D, int *count, FILE *statfile);
void TBRswitch(tree *T, edge *es, edge *et, edge *eb);

#endif /*TBR_H_*/

