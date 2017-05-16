#ifndef FASTME_H_
#define FASTME_H_

#include "interface_options.h"
#include "NNI.h"
#include "bNNI.h"
#include "MVR.h"
#include "TBR.h"
#include "distance.h"
#include "p_lk.h"
#include "p_bootstrap.h"

int numSpecies;
int edgeCount;

allseq *p_bootstraps (allseq *alldata, int ns, int *site_num);
void OpenFiles (Options *options);
void printOptions (Options *options);
void printMatrix(double **D, int size, set *nodes, FILE *ofile);
void Print_Time_Info(time_t t_beg, time_t t_end);

#endif /*FASTME_H_*/

