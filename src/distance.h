#ifndef DISTANCE_H_
#define DISTANCE_H_

#include "FastDist.h"
#include "random.h"
#include "inputs.h"

int countStateChanges(char *s, char *t, int length, char c1, char c2, int *filter);
void ijFilter(int *filter, char *s1, char *s2, int itype, int seqlength);
void calcDNATransitions(double **P,char *s1, char *s2, int length, int *filter, int numSelected);
int seqCharMatches(char *s, int length, char c, int *filter);
int matrixCharMatches(char **s, int numSeqs, int length, char c, int *filter);
double *calcDNAStationaryProbs(char **s, int numSeqs, int length, int *filter, int numSelected);
double *calcProteinStationaryProbs(char **s, int numSeqs, int length, int *filter, int numSelected);
void calcDNATransitionProbs(double **P,char *s1, char *s2, int length, int *filter, int numSelected);
void calcProteinTransitionProbs(double **P,char *s1, char *s2, int length, int *filter, int numSelected);
double calcTransitionRate(double **P);
double calcTransversionRate(double **P);
double calcK2P(double a, double b, double gamma);
void calcF84AuxProbs(double *Pi, double *Mloc, double *Nloc, double *Ploc);
double calcF84(double a, double b, double gamma, double Mloc, double Nloc, double Ploc);
void calcTNAuxProbs(double *Pi,double *Mloc, double *Nloc, double *PR, double *PY);
double calcTN93(double aR, double aY, double b, double PR, double PY, double PAPG, double PCPT, double gamma);
double XX(double PR, double PY, double b, double gamma);
int factorial(int n);
int *nextPerm(int *p, int index, int size, int length);
double permDiagProduct(double **P, int *p, int d);
double det(double **P, int d);
double logdet(double **P, double *Pi1, double *Pi2);
int support(int *v, int length);
double HammingDistance(char *v1, char *v2, int *filter, int length, int numSelected);
double protDiff(double *P);
double protFormula(double b, double gamma, double pdiff);
int aaIndex(char s, char *alphabet, int d);
double simScore(char *s1, char *s2, double **scoreMatrix, int seqlength, char *alphabet, int alphabetSize);
double expectedProtSimScore(double *P, double **scoreMatrix, int alphabetSize);
double scoreDistij(int i,int j,char *si, char *sj, int seqLength, double simExp,
	double *simDiags, double **scoreMatrix, char *alphabet, int alphabetSize);
void scoreDist(double *P, char **data, int numSpecies, int seqLength,
	double **scoreMatrix, double **D, char *alphabet, int alphabetSize);
void gapCheckFilter(int *filter, int itype, int seqlength, int numSeqs, char **data);
double **makeDistMatrix(char **data, int numSeqs, int numSites, double gamma,
	int model, int itype, int *filter, double **scoreMatrix, int scoreMatrixSize,
	char *alphabet, boolean gapCheck);
void symmetrizeDoubleMatrix(double **X, int n);

#endif /*DISTANCE_H_*/

