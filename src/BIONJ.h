#ifndef BIONJ_H_
#define BIONJ_H_

#include "newick.h"

typedef struct word
{
  char name[MAX_LABEL_LENGTH];
  struct word *suiv;
}WORD;

typedef struct pointers
{
  WORD *head;
  WORD *tail;
}POINTERS;


tree *bionj (double **D, set *species, int n, boolean isNJ);
void Initialize(double **D, set *species, double **delta, POINTERS *trees, int n);
void Print_outputChar(int i, POINTERS *trees, char *output);
boolean Symmetrize(double **delta, int n);
void Concatenate(char chain1[MAX_LABEL_LENGTH], int ind, POINTERS *trees, int post);
double Distance(int i, int j, double **delta);
double Variance(int i, int j, double **delta);
int Emptied(int i, double **delta);
double Sum_S(int i, double **delta);
void Compute_sums_Sx(double **delta, int n);
void Best_pair(double **delta, int r, int *a, int *b, int n);
double Finish_branch_length(int i, int j, int k, double **delta);
void FinishStr (double **delta, int n, POINTERS *trees, char *StrTree);
double Agglomerative_criterion(int i, int j, double **delta, int r);
double Branch_length(int a, int b, double **delta, int r);
double Reduction4(int a, double la, int b, double lb, int i, double lamda, double **delta);
double Reduction10(int a, int b, int i, double lamda, double vab, double **delta);
double Lamda(int a, int b, double vab, double **delta, int n, int r);

#endif /*BIONJ_H_*/

