#ifndef INPUTS_H_
#define INPUTS_H_

#include "graph.h"
#include "heap.h"

void compareSets(tree *T, set *S);
void freeCharMatrix(char **D, int size);
void freeMatrix(double **D, int size);
double **loadMatrix(FILE *ifile, int *size, set *S);
void partitionSizes(tree *T);
char **loadSequences(FILE *ifile, int *numSeqs, set *taxa, int *inputSize);
double **loadScoreMatrix(int *d, FILE *ifile, char *alphabet);

#endif /*INPUTS_H_*/

