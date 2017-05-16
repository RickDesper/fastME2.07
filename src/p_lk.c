/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "p_lk.h"

/*********************************************************/

matrix *ML_Dist(allseq *data, model *mod)
{
  int i,j,k,l;
  phydbl init;
//  int n_catg;
  phydbl d_max,sum;
  matrix *mat;
  allseq *twodata,*tmpdata;
  int state0, state1,len;
  phydbl *F;
  eigen *eigen_struct;

  tmpdata             = (allseq *)mCalloc(1,sizeof(allseq));
  tmpdata->c_seq      = (seq **)mCalloc(2,sizeof(seq *));
  tmpdata->b_frq      = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  tmpdata->ambigu     = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  F                   = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl ));
  eigen_struct        = (eigen *)Make_Eigen_Struct(mod);

  tmpdata->n_otu      = 2;

  tmpdata->crunch_len = 0;
  tmpdata->init_len   = 0;
  for (i=0; i<data->crunch_len; i++)
    {
      if(data->wght[i] > .0)
	{
	  tmpdata->crunch_len++;
	  tmpdata->init_len+=(int)data->wght[i];
	}
    }

  mat = JC69_Dist(data,mod);

  for (i=0; i<mod->n_catg; i++) // Don't use gamma distribution
    {
      mod->gamma_rr[i]      = 1.0;
      mod->gamma_r_proba[i] = 1.0;
    }
/*
  n_catg = mod->n_catg;
  mod->n_catg = 1;
*/
  for (j=0; j<data->n_otu-1; j++)
    {
      tmpdata->c_seq[0]       = data->c_seq[j];
      tmpdata->c_seq[0]->name = data->c_seq[j]->name;
      tmpdata->wght           = data->wght;
      
      for(k=j+1;k<data->n_otu;k++)
	{
	  tmpdata->c_seq[1]       = data->c_seq[k];
	  tmpdata->c_seq[1]->name = data->c_seq[k]->name;
	  
	  twodata = Compact_CSeq(tmpdata,mod);
	  for (l=0; l<mod->ns; l++)
	    twodata->b_frq[l] = data->b_frq[l];
	  Check_Ambiguities(twodata,mod->datatype,1);
	  Hide_Ambiguities(twodata);
	  
	  init = mat->dist[j][k];
	  if((init == DIST_MAX) || (init < .0)) init = 0.1;
	  	  
	  d_max = init;
	  
	  for (i=0; i<mod->ns*mod->ns; i++)
	    F[i]=.0;
	  len = 0;
	  for (l=0; l<twodata->c_seq[0]->len; l++)
	    {
	      state0 = Assign_State(twodata->c_seq[0]->state+l);
	      state1 = Assign_State(twodata->c_seq[1]->state+l);
	      if((state0 > -1) && (state1 > -1))
		{
		  F[mod->ns*state0+state1] += twodata->wght[l];
		  len += (int)twodata->wght[l];
		}
	    }
	  if(len > .0)
	    for (i=0; i<mod->ns*mod->ns; i++)
	      F[i] /= (phydbl)len;
	  	  
	  sum = 0.;
	  for (i=0; i<mod->ns*mod->ns; i++)
	    sum += F[i];
	  if(sum < .001) d_max = -1.;
	  else if((sum > 1. - .001) && (sum < 1. + .001)) Opt_Dist_F(&(d_max),F,mod);
	  else
	    {
	      printf("\n . sum = %f.\n",sum);
	      printf("\n Error: file %s at line %d.\n",__FILE__,__LINE__);
	      exit(EXIT_FAILURE);
	    }

	  // BRENT
	  
	  d_max = (d_max >= DIST_MAX ? DIST_MAX : d_max);
//	  if(d_max >= DIST_MAX)
//	    {
// 	      printf("\n. Large distance encountered between %s and %s sequences.",
// 		     tmpdata->c_seq[1]->name,
// 		     tmpdata->c_seq[0]->name);
//	      d_max = DIST_MAX;
//	    }
	  
	  // Do not correct for dist < BL_MIN, otherwise Fill_Missing_Dist 
	  // will not be called
	  
	  mat->dist[j][k] = d_max;
	  mat->dist[k][j] = mat->dist[j][k];
	  Free_Cseq(twodata);
	}
    }

//  mod->n_catg = n_catg;
  free(tmpdata->ambigu);
  free(tmpdata->b_frq);
  free(tmpdata->c_seq);
  free(tmpdata);
  Free_Eigen(eigen_struct);
  free(F);

  return mat;
}

/*********************************************************/

phydbl Lk_Dist(phydbl *F, phydbl dist, model *mod)
{
  int i,j,k;
  phydbl len,lnL,tmp;

  len = -1.;
  for (k=0; k<mod->n_catg; k++)
    {
      len = dist*mod->gamma_rr[k];
      if(len < BL_MIN)
        len = BL_MIN;
      else if(len > BL_MAX)
        len = BL_MAX;
      PMat(len,mod,&(mod->Pij_rr[k]));
    }

  lnL = .0;
  for (i=0; i<mod->ns; i++)
    {
      for (j=0; j<mod->ns; j++)
	{
	  tmp = .0;
	  for (k=0; k<mod->n_catg; k++)
	    {
	      tmp +=
		mod->gamma_r_proba[k] *
		mod->pi[i] *
		(phydbl)(mod->Pij_rr[k][i][j]);
	    }
	  lnL += F[mod->ns*i+j] * (phydbl)log(tmp);
	}
    }
  return lnL;
}

/*********************************************************/


