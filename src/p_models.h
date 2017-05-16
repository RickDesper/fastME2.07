/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef P_MODELS_H
#define P_MODELS_H

#include "p_eigen.h"

void Init_Model(allseq *data, model *mod, FILE *in);
int Init_Qmat_Dayhoff(double *daa, phydbl *pi);
int Init_Qmat_JTT(double *daa, phydbl *pi);
int Init_Qmat_MtREV(double *daa, phydbl *pi);
int Init_Qmat_LG(double *daa, phydbl *pi);
int Init_Qmat_WAG(double *daa, phydbl *pi);
int Init_Qmat_DCMut(double *daa, phydbl *pi);
int Init_Qmat_RtREV(double *daa, phydbl *pi);
int Init_Qmat_CpREV(double *daa, phydbl *pi);
int Init_Qmat_VT(double *daa, phydbl *pi);
void PMat(phydbl l, model *mod, double ***Pij);
void PMat_Zero_Br_Len(model  *mod, double ***Pij);
void PMat_Empirical(phydbl l, model *mod, double ***Pij);

#endif /*P_MODELS_H_*/

