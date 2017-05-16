/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef P_OPTIMIZ_H
#define P_OPTIMIZ_H

#include "p_models.h"
#include "p_lk.h"

void Opt_Dist_F(phydbl *dist, phydbl *F, model *mod);
phydbl Dist_F_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max, 
	phydbl *param, phydbl *F, model *mod);


#endif /*P_OPTIMIZ_H_*/
