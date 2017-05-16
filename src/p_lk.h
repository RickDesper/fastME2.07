/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef P_LK_H
#define P_LK_H

#include "p_optimiz.h"

matrix *ML_Dist(allseq *data,model *mod);
phydbl Lk_Dist(phydbl *F, phydbl dist, model *mod);


#endif /*P_LK_H_*/






