#ifndef __FastDist__
#define __FastDist__

#include "CompDist.h"

void ResolveAmbiguity( Sseq* seq1, Sseq* seq2, strIJ* Sij, int m, int model, int *filter);
void CompAmbStatistic(double* tp, double* AmbVec1, double* AmbVec2, double *puTs,
	double *pyTs, double *Tv);
void CorrectStatisticsIJ(strIJ* Sij, Sseq* seq1, Sseq* seq2, int m, int model, int *filter);
double** FastDist(char** data, int n, int m, int model, double gamma, int *filter);
void UpdateProbSij(strIJ* Sij, int m);
strIJ* StatisticIJ(char* s1, char* s2, int m, int *filter);
int IndDNA(char c);
Sseq*  Statistic(char* s, int m, int *filter);
double *GetFre(char c);
void Free(Sseq* S, int m);
int NearestNeighbor(int i, double** D, int n );

#endif

