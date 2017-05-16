#ifndef GME_H_
#define GME_H_

#include "traverse.h"

void fillTableUp(edge *e, edge *f, double **A, double **D, tree *T);
void OLSext(edge *e, double **A);
void OLSint(edge *e, double **A);
void assignOLSWeights(tree *T, double **A);
void makeOLSAveragesTable(tree *T, double **D, double **A);
void GMEcalcDownAverage(node *v, edge *e, double **D, double **A);
void GMEcalcUpAverage(node *v, edge *e, double **D, double **A);
void GMEcalcNewvAverages(tree *T, node *v, double **D, double **A);
double wf4(double lambda, double lambda2, double D_AB, double D_AC, 
	double D_BC, double D_Av, double D_Bv, double D_Cv);
void testEdge(edge *e, node *v, double **A);
tree *GMEaddSpecies(tree *T,node *v, double **D, double **A);
void GMEupdateAveragesMatrix(double **A, edge *e, node *v, node *newNode);
void GMEsplitEdge(tree *T, node *v, edge *e, double **A);
void updateSubTreeAverages(double **A, edge *e, node *v, int direction);
void assignBottomsize(edge *e);
void assignTopsize(edge *e, int numLeaves);
void assignAllSizeFields(tree *T);

#endif /*GME_H_*/

