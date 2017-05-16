#ifndef BME_H_
#define BME_H_

#include "traverse.h"

void BalWFext(edge *e, double **A);
void BalWFint(edge *e, double **A);
void assignBMEWeights(tree *T, double **A);
void BMEcalcDownAverage(tree *T, node *v, edge *e, double **D, double **A);
void BMEcalcUpAverage(tree *T, node *v, edge *e, double **D, double **A);
void BMEcalcNewvAverages(tree *T, node *v, double **D, double **A);
void updatePair(double **A, edge *nearEdge, edge *farEdge, node *v,
	node *root, double dcoeff, int direction);
void updateSubTree(double **A, edge *nearEdge, node *v, node *root,
	node *newNode, double dcoeff, int direction);
void BMEupdateAveragesMatrix(double **A, edge *e, node *v,node *newNode);
double wf3(double D_AB, double D_AC, double D_kB, double D_kC);
void BMEtestEdge(edge *e, node *v, double **A);
void BMEsplitEdge(tree *T, node *v, edge *e, double **A);
tree *BMEaddSpecies(tree *T,node *v, double **D, double **A);
void calcUpAverages(double **D, double **A, edge *e, edge *g);
void makeBMEAveragesTable(tree *T, double **D, double **A);

#endif /*BME_H_*/

