#ifndef __CompDist__
#define __CompDist__

#include "utils.h"

/*
Input: alignment (char** data, 
                   int n, number of sequences
                   int m, number of sites
                   model M, 
                   double gamma);

output: distance matrix: double** D

- Statistics
- Inital Dist
- Ambiguity update
- Recompute Distance
 */

//structure of a sequence
// 0 --> A
// 1 --> C
// 2 --> G
// 3 --> T
// 4 --> unknown
// 5 --> # ambiguity

//statistics for sequence s[i]
typedef struct {
  char *p;
  int Num[10];
  //  int* mark_amb;   // mark ambiguous site
  double ** AmbVec; // ambigouus vectors of site: AmbVec[i]{0,1,2,3} = { pA, pC, pG, pT}
} Sseq;

// statistics for sequences s[i], s[j]
typedef struct {
  double puTs, pyTs, Tv; // purine Transitions, pyrimidine Transitions, Transversions
  double ts_purine_prob, ts_pyrimidine_prob, tv_prob;
  double idprob;
  int del; //number of deleted Positions
} strIJ;


double ComputeDistance(int model, strIJ* Sij, int* Num, int m, double gamma); // computing distance between i, j; Num[0,1,2,3] = {#A, #C, #G,
double Compute_TN93(strIJ* Sij, int* Num, int m, double gamma); // computing distance between i, j; Num[0,1,2,3] = {#A, #C, #G, #T}
double Compute_JC69(strIJ* Sij, int m, double gamma);
double Compute_K2P(strIJ* Sij, int m, double gamma);
double Compute_F84(strIJ* Sij, int* Num, int m, double gamma);
double* ComputeTP(int model, strIJ * Sij, int m);
double* Compute_TP_TN93(double ts_purine_prob, double ts_pyrimidine_prob, double tv_prob);
double* Compute_TP_K2P(double ts_purine_prob, double ts_pyrimidine_prob, double tv_prob, double idprob);
double *Compute_TP_F84(double ts_purine_prob, double ts_pyrimidine_prob, double tv,  double idprob);
double *Compute_TP_JC69(double ts_purine_prob, double ts_pyrimidine_prob, double tv_prob, double idprob);

#endif  /*__CompDist__*/

