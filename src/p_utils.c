/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                         Utilities functions for ditance matrix            ;
;                         computation from amino acids sequences            ;
;                                                                           ;
;                         Vincent Lefort                                    ;
;                                                                           ;
;                         LIRMM - Montpellier- France                       ;
;                         vincent.lefort@lirmm.fr                           ;
;                                                                           ;
;                         code is heavily borrowed from PhyML :             ;
;                                                                           ;
;                         Stéphane Guindon                                  ;
;                                                                           ;
;                         LIRMM - Montpellier- France                       ;
;                         guindon@lirmm.fr                                  ;
;                                                                           ;
;                         Olivier Gascuel                                   ;
;                                                                           ;
;                         LIRMM - Montpellier- France                       ;
;                         gascuel@lirmm.fr                                  ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

#include "p_utils.h"


/*********************************************************/

double **Copy_PMat_to_DMat (matrix *mat)
{
  int i,j;
  double **D;
  D = initDoubleMatrix(mat->n_otu);
  for (i=0; i<mat->n_otu; i++)
    for (j=0; j<mat->n_otu; j++)
      D[i][j] = (double)mat->dist[i][j];

  return D;
}

/*********************************************************/

void Free_Seq(seq **d, int n_otu)
{
  int i;
  for (i=0; i<n_otu; i++)
    {
      free(d[i]->name);
      free(d[i]->state);
      if ((NULL != d[i]->is_ambigu) && (d[i]->is_ambigu))
        free(d[i]->is_ambigu);
      free(d[i]);
    }
  free(d);
  return;
}

/*********************************************************/

void Free_Cseq(allseq *data)
{
  int i;
  
  free(data->invar);
  free(data->wght);
  free(data->ambigu);
  free(data->b_frq);
  free(data->sitepatt);
  for (i=0; i<data->n_otu; i++)
    {
      free(data->c_seq[i]->name);
      if(data->c_seq[i]->state) 
	{
	  free(data->c_seq[i]->state);
	  if((NULL != data->c_seq[i]->is_ambigu) && (data->c_seq[i]->is_ambigu))
	    free(data->c_seq[i]->is_ambigu);
	}
      free(data->c_seq[i]);
    }
  free(data->c_seq);
  free(data);
  return;
}

/*********************************************************/

void Free_Prefix_Tree(pnode *n, int size)
{
  int i;
  
  for (i=0; i<size; i++)
    {
      if(n->next[i])
	{
	  Free_Prefix_Tree(n->next[i],size);
	}
    }
  Free_Pnode(n);
  return;
}

/*********************************************************/

void Free_Pnode(pnode *n)
{
  free(n->next);
  free(n);
  return;
}

/*********************************************************/

void Free_Eigen(eigen *eigen_struct)
{
  free(eigen_struct->space_int);
  free(eigen_struct->space);
  free(eigen_struct->e_val);
  free(eigen_struct->e_val_im);
  free(eigen_struct->r_e_vect);
  free(eigen_struct->r_e_vect_im);
  free(eigen_struct->l_e_vect);
  free(eigen_struct->q);
  free(eigen_struct);
  return;
}

/*********************************************************/

seq **Get_Seq (FILE *in, boolean interleaved, int *n_otu, int *len, int itype, set *taxa)
{
  seq **data;
  int i,j;
  char **buff;
  int n_unkn,n_removed,pos;
  int *remove;  

  if (interleaved)
    data = Read_Seq_Interleaved(in,n_otu,taxa);
  else
    data = Read_Seq_Sequential(in,n_otu,taxa);

  *len = data[0]->len;
  if (verbose)
    {
      printf("  . Number of sequences: %d.\n",*n_otu);
      printf("  . Length of sequences: %d.\n",*len);
    }

  if(data)
    {
      buff = (char **)mCalloc(*n_otu,sizeof(char *));
      for (i=0; i<*n_otu; i++)
        buff[i] = (char *)mCalloc(data[0]->len,sizeof(char));
      remove = (int *)mCalloc(data[0]->len,sizeof(int));

      n_removed = 0;

      for (i=0; i<data[0]->len; i++)
	{
	  for (j=0; j<*n_otu; j++)
	    {
	      if((data[j]->state[i] == '?') || (data[j]->state[i] == '-')) data[j]->state[i] = 'X';
	      if((itype == DNA) && (data[j]->state[i] == 'N')) data[j]->state[i] = 'X';
	      if(data[j]->state[i] == 'U') data[j]->state[i] = 'T';
	    }

	  n_unkn = 0;
	  for (j=0; j<*n_otu; j++)
	    if(data[j]->state[i] == 'X')
	      n_unkn++;

	  if(n_unkn == *n_otu)
	    {
	      remove[i] = 1;
	      n_removed++;
	    }

	  for (j=0; j<*n_otu; j++)
	    buff[j][i] = data[j]->state[i];
	}

      if(n_removed > 0)
	{
	  if(itype == DNA)
	    printf("\n . %d sites are made from completely undetermined states ('X', '-', '?' or 'N')...\n",n_removed);
	  else
	    printf("\n . %d sites are made from completely undetermined states ('X', '-', '?')...\n",n_removed);
	}
	
      pos = 0;
      for (i=0; i<data[0]->len; i++)
	{
// 	  if(!remove[i])
// 	    {
	      for (j=0; j<*n_otu; j++)
	        data[j]->state[pos] = buff[j][i];
	      pos++;
// 	    }
	}

      for (i=0; i<*n_otu; i++)
        {
         data[i]->len = pos;
         if (NULL != buff[i])
           free(buff[i]);
        }
      if (NULL != buff)
        free(buff);
      if (NULL != remove)
        free(remove);
    }
  return data;
}

/*********************************************************/

seq **Read_Seq_Sequential(FILE *in, int *n_otu, set *taxa)
{
  int i;
  node *v;
  char *line;
  int len,readok;
  seq **data;
  char c;
  char *format = (char *)mCalloc(MAX_NAME_LENGTH, sizeof(char));

  line = (char *)mCalloc(MAX_LINE_LENGTH,sizeof(char));

  readok = len = 0;
  do
    {
      if(fscanf(in,"%s",line) == EOF)
	{
	  free(line);
	  return NULL;
	}
      else
	{
	  if(strncmp(line,"\n",1) && strncmp(line,"\n",1) && strncmp(line,"\t",1))
	    {
	      *n_otu = atoi(line);
	      data = (seq **)mCalloc(*n_otu,sizeof(seq *));
	      if(*n_otu <= 0)
	        Exit("Problem with sequence format.");
	      fscanf(in,"%s",line);
	      len = atoi(line);
	      if(len <= 0)
	        Exit("Problem with sequence format.");
	      else readok = 1;
	    }
	}
    }while(!readok);

  while(((c=fgetc(in))!='\n') && (c != ' ') && (c != '\r') && (c != '\t'));

  for (i=0; i<*n_otu; i++)
    {
      data[i] = (seq *)mCalloc(1,sizeof(seq));
      data[i]->len = 0;
      data[i]->name = (char *)mCalloc(MAX_NAME_LENGTH,sizeof(char));
      data[i]->state = (char *)mCalloc(MAX_SEQ_LENGTH,sizeof(char));
      data[i]->is_ambigu = NULL;
      snprintf(format, MAX_NAME_LENGTH, "%%%ds", MAX_NAME_LENGTH);
      fscanf(in, format, data[i]->name);

      v = makeNode(data[i]->name,-1);
      v->index2 = i;
      taxa = addToSet(v,taxa);
      
      while(data[i]->len < len)
	Read_One_Line_Seq(&data,i,in);

      if(data[i]->len != len)
	{
	  printf("\n Err: Problem with species %s's sequence (check the format).\n",
		 data[i]->name);
	  exit(EXIT_FAILURE);
	}
    }

  //   fgets(line,MAX_LINE_LENGTH,in);
  // inter data sets

  free(format);
  free(line);
  return data;
}

/*********************************************************/

void Free_Mat(matrix *mat)
{
  int i;

  for (i=0; i<mat->n_otu; i++)
    {
      free(mat->P[i]);
      free(mat->Q[i]);
      free(mat->dist[i]);
      free(mat->name[i]);
    }

  free(mat->P);
  free(mat->Q);
  free(mat->dist);
  free(mat->name);
  free(mat->on_off);
  free(mat);
  
  return;
}

/*********************************************************/

seq **Read_Seq_Interleaved(FILE *in, int *n_otu, set *taxa)
{
  int i,end,num_block;
  node *v;
  char *line;
  int len,readok;
  seq **data;
  char c;
  char *format;

  line = (char *)mCalloc(MAX_LINE_LENGTH,sizeof(char));
  format = (char *)mCalloc(MAX_NAME_LENGTH, sizeof(char));

  readok = len = 0;
  do
    {
      if(fscanf(in,"%s",line) == EOF)
	{
	  free(format);
	  free(line);
	  return NULL;
	}
      else
	{
	  if(strncmp(line,"\n",1) && strncmp(line,"\r",1) && strncmp(line,"\t",1))
	    {
	      *n_otu = atoi(line);
	      data = (seq **)mCalloc(*n_otu,sizeof(seq *));
	      if(*n_otu <= 0)
	        Exit("Problem with sequence format.");
	      fscanf(in,"%s",line);
	      len = atoi(line);
	      if(len <= 0)
	        Exit("Problem with sequence format.");
	      else readok = 1;
	    }
	}
    }while(!readok);


  while(((c=fgetc(in))!='\n') && (c != ' ') && (c != '\r') && (c != '\t'));

  end = 0;
  for (i=0; i<*n_otu; i++)
    {
      data[i] = (seq *)mCalloc(1,sizeof(seq));
      data[i]->len = 0;
      data[i]->name = (char *)mCalloc(MAX_NAME_LENGTH,sizeof(char));
      data[i]->state = (char *)mCalloc(MAX_SEQ_LENGTH,sizeof(char));
      data[i]->is_ambigu = NULL;
      snprintf(format, MAX_NAME_LENGTH, "%%%ds", MAX_NAME_LENGTH);
      fscanf(in, format, data[i]->name);
      
      v = makeNode(data[i]->name,-1);
      v->index2 = i;
      taxa = addToSet(v,taxa);
      
      if(!Read_One_Line_Seq(&data,i,in))
	{
	  end = 1;
	  if((i != *n_otu) && (i != *n_otu-1))
	    {
	      printf("\n Err: Problem with species %s's sequence.\n",data[i]->name);
	      exit(EXIT_FAILURE);
	    }
	  break;
	}
    }

  if(data[0]->len == len)
    end = 1;

  if(!end)
    {
      end = 0;

      num_block = 1;
      do
	{
	  num_block++;

	  // interblock
	  if(!fgets(line,MAX_LINE_LENGTH,in)) break;

	  if(line[0] != 13 && line[0] != 10)
	    {
	      printf("\n Error: One or more missing sequences in block %d.\n",num_block-1);
	      exit(EXIT_FAILURE);
	    }
	  
	  for (i=0; i<*n_otu; i++)
	    if(data[i]->len != len)
	      break;
	  
	  if(i == *n_otu) break;
	  
	  for (i=0; i<*n_otu; i++)
	    {
	      if(data[i]->len > len)
		{
		  printf("\n Error: Problem with species %s's sequence.\n",data[i]->name);
		  printf("\n        Observed length=%d expected length=%d.\n",data[i]->len,len);
		  exit(EXIT_FAILURE);
		}
	      else if(!Read_One_Line_Seq(&data,i,in))
		{
		  end = 1;
		  if((i != *n_otu) && (i != *n_otu-1))
		    {
		      printf("\n Error: Problem with species %s's sequence.\n",data[i]->name);
		      exit(EXIT_FAILURE);
		    }
		  break;
		}
	    }
	}while(!end);
    }

  for (i=0; i<*n_otu; i++)
    {
      if(data[i]->len != len)
	{
	  printf("\n Error: Check sequence '%s' length.\n",data[i]->name);
	  exit(EXIT_FAILURE);
	}
    }

  free(format);
  free(line);
  
  return data;
}

/*********************************************************/

int Read_One_Line_Seq(seq ***data, int num_otu, FILE *in)
{
  char c;

  c=' ';
  while(1)
    {
//       if((c == EOF) || (c == '\n') || (c == '\r')) break;
      if((c == EOF) || (c == 13) || (c == 10))
        break;
      else if((c==' ') || (c=='\t'))
        {
          c=(char)fgetc(in);
          continue;
        }
      Uppercase(&c);

      if (strchr("ABCDEFGHIKLMNOPQRSTUVWXYZ?-.", c) == NULL)
	{
	  printf("\n Error: bad symbol: \"%c\" at position %d of species %s.\n",
		 c,(*data)[num_otu]->len,(*data)[num_otu]->name);
	  exit(EXIT_FAILURE);
	}

      if(c == '.')
	{
	  c = (*data)[0]->state[(*data)[num_otu]->len];
	  if(!num_otu)
	    Exit("Symbol \".\" should not appear in the first sequence.");
	}
      (*data)[num_otu]->state[(*data)[num_otu]->len]=c;
      (*data)[num_otu]->len++;
      c = (char)fgetc(in);
    }
  if(c == EOF)
    return 0;
  else
    return 1;
}

/*********************************************************/

model *Make_Model_Basic()
{
  model *mod;

  mod                     = (model *)mCalloc(1,sizeof(model));
//  mod->modelname          = (char *)mCalloc(MAX_NAME_LENGTH,sizeof(char));
//  mod->custom_mod_string  = (char *)mCalloc(MAX_NAME_LENGTH,sizeof(char));
  mod->user_b_freq        = (phydbl *)mCalloc(MAX_NAME_LENGTH,sizeof(phydbl));

  mod->rr                 = (phydbl *)mCalloc(6,sizeof(phydbl));
  mod->rr_val             = (phydbl *)mCalloc(6,sizeof(phydbl));
  mod->rr_num             = (int *)mCalloc(6,sizeof(int *));
  mod->n_rr_per_cat       = (int *)mCalloc(6,sizeof(int));
//  mod->s_opt              = (optimiz *)Alloc_Optimiz();

  return mod;
}

/*********************************************************/

void Make_Model_Complete(model *mod)
{
  int i,j;

  mod->pi             = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  
  mod->gamma_r_proba  = (phydbl *)mCalloc(mod->n_catg,sizeof(phydbl));
  mod->gamma_rr       = (phydbl *)mCalloc(mod->n_catg,sizeof(phydbl));
  mod->pi_unscaled    = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));

  mod->Pij_rr   = (double***)mCalloc(mod->n_catg,sizeof(double **));

  for (i=0; i<mod->n_catg; i++)
    {
      mod->Pij_rr[i] = (double **)mCalloc(mod->ns,sizeof(double *));
      for (j=0; j<mod->ns; j++)
        mod->Pij_rr[i][j] = (double *)mCalloc(mod->ns,sizeof(double));
    }

  mod->qmat      = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
  mod->qmat_buff = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
  mod->eigen     = (eigen *)Make_Eigen_Struct(mod);
  
  if(mod->n_rr_branch)
    {
      mod->rr_branch   = (phydbl *)mCalloc(mod->n_rr_branch,sizeof(phydbl));
      mod->p_rr_branch = (phydbl *)mCalloc(mod->n_rr_branch,sizeof(phydbl));
    }
  return;
}

/*********************************************************/
/*
optimiz *Alloc_Optimiz()
{
  optimiz *s_opt;
  s_opt = (optimiz *)mCalloc(1,sizeof(optimiz));
  return s_opt;
}
*/
/*********************************************************/

void Set_Defaults_Model(model *mod, double alpha)
{
  int i;

//  strncpy(mod->modelname,"LG",2);
//  strncpy(mod->custom_mod_string,"000000",6);
  mod->datatype                = PROTEIN;
  mod->whichmodel              = LG;
  mod->n_catg                  = 1;
  mod->kappa                   = 4.0;
//  mod->alpha                   = 1.0;
  mod->alpha                   = (phydbl)alpha;
  mod->lambda                  = 1.0;
  mod->bootstrap               = 0;
  mod->invar                   = 0;
  mod->pinvar                  = 0.0;
  mod->stepsize                = 1;
  mod->ns                      = 20;
  mod->n_diff_rr               = 0;
  for (i=0; i<6; i++)
    mod->rr_val[i]             = 1.0;
  for (i=0; i<4; i++)
    mod->user_b_freq[i]        = -1.;
//  mod->m4mod                   = NULL;
//  mod->use_m4mod               = 0;
  mod->n_rr_branch             = 0;
  mod->rr_branch_alpha         = 0.1;
  return;
}

/*********************************************************/
/*
void Set_Defaults_Optimiz(optimiz *s_opt)
{
  s_opt->print                = 1;
//  s_opt->last_opt             = 1;
  s_opt->opt_alpha            = 0;
  s_opt->opt_kappa            = 0;
  s_opt->opt_bl               = 1;
  s_opt->opt_lambda           = 0;
  s_opt->opt_pinvar           = 0;
  s_opt->opt_num_param        = 0;
  s_opt->opt_cov_delta        = 0;
//  s_opt->opt_cov_alpha        = 0;
//  s_opt->opt_cov_free_rates   = 0;
//  s_opt->opt_rr               = 0;
//  s_opt->init_lk              = UNLIKELY;
//  s_opt->n_it_max             = 1000;
//  s_opt->opt_topo             = 1;
//  s_opt->topo_search          = NNI_MOVE;
//  s_opt->random_input_tree    = 0;
//  s_opt->n_rand_starts        = 5;
//  s_opt->brent_it_max         = 500;
//  s_opt->steph_spr            = 1;
//  s_opt->user_state_freq      = 0;
//  s_opt->min_diff_lk_local    = 1.E-05;
//  s_opt->min_diff_lk_global   = 1.E-03;
//  s_opt->spr_step_after_nnis  = 0;
//  s_opt->p_moves_to_examine   = 0.1;
//  s_opt->fast_nni             = 0;
//  s_opt->greedy               = 0;
//  s_opt->general_pars         = 0;
//  s_opt->tree_size_mult       = 1;
//  s_opt->wim_n_rgrft          = -1;
//  s_opt->wim_n_globl          = -1;
//  s_opt->wim_max_dist         = -1;
//  s_opt->wim_n_optim          = -1;
//  s_opt->wim_n_best           = -1;
//  s_opt->wim_inside_opt       =  0;
  return;
}
*/
/*********************************************************/

eigen *Make_Eigen_Struct(model *mod)
{
  eigen *eig;

  eig              = (eigen *)mCalloc(1,sizeof(eigen));
  eig->size        = mod->ns;
  eig->space       = (double *)mCalloc(2*mod->ns,sizeof(double));
  eig->space_int   = (int *)mCalloc(2*mod->ns,sizeof(int));
  eig->e_val       = (double *)mCalloc(mod->ns,sizeof(double));
  eig->e_val_im    = (double *)mCalloc(mod->ns,sizeof(double));
  eig->r_e_vect    = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
  eig->r_e_vect_im = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
  eig->l_e_vect    = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
  eig->q           = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));

  return eig;
}

/*********************************************************/

allseq *Compact_Seq(seq **data, model *mod, boolean rm_ambigu)
{
  allseq *alldata_tmp,*alldata;
  int i,j,k,site;
  int n_patt,which_patt,n_invar;
  char **sp_names;
  int n_otu, n_sites;
  pnode *proot;
  int compress;
  int n_ambigu,is_ambigu;

  n_otu        = mod->n_otu;
  n_patt       = 0;
  which_patt   = 0;

  sp_names = (char **)mCalloc(n_otu,sizeof(char *));
  for (i=0; i<n_otu; i++)
    {
      sp_names[i] = (char *)mCalloc(MAX_NAME_LENGTH,sizeof(char));
      strncpy(sp_names[i],data[i]->name, MAX_NAME_LENGTH);
    }
  alldata_tmp = Make_Cseq(n_otu,data[0]->len,data[0]->len,sp_names);
  proot       = (pnode *)Create_Pnode(MAX_NAME_LENGTH);
 
  for (i=0; i<n_otu; i++)
    free(sp_names[i]);
  free(sp_names);

  if(data[0]->len%mod->stepsize)
    {
      printf("\n Error: Sequence length is not a multiple of %d.\n",mod->stepsize);
      exit(EXIT_FAILURE);
    }
  
  compress = 1;
/*   compress = 0; */
  n_ambigu = 0;
  is_ambigu = 0;

//  Fors(site,data[0]->len,io->mod->stepsize)
  for(site=0; site<data[0]->len; site+=mod->stepsize)
    {
      if(rm_ambigu)
	{
	  is_ambigu = 0;
	  for (j=0; j<n_otu; j++)
	    {
	      //if(Is_Ambigu(data[j]->state+site,mod->datatype,mod->stepsize))
	      if(Is_Ambigu(data[j]->state+site))
		{
		  break;
		}
	    }
	  if(j != n_otu)
	    {
	      is_ambigu = 1;
	      n_ambigu++;
	    }
	}

      if(!is_ambigu)
	{
	  if(compress)
	    {
	      which_patt = -1;
	      Traverse_Prefix_Tree(site,-1,&which_patt,&n_patt,data,mod,proot);
	      if(which_patt == n_patt-1) /* New pattern found */
		{
		  n_patt--;
		  k=n_patt;
		}
	      else
		{
		  k = n_patt-10;
		}
	    }
	  else
	    {
	      Warning ("Sequences are not compressed.");
	      k = n_patt;
	    }
	  
	  if(k == n_patt) /* add a new site pattern */
	    {
	      for (j=0; j<n_otu; j++)
		Copy_One_State(data[j]->state+site,
			       alldata_tmp->c_seq[j]->state+n_patt,
			       mod->stepsize);
	      
	      for (i=0; i<n_otu; i++)
		{
		  for (j=0; j<n_otu; j++)
		    {
		      if(!(Are_Compatible(alldata_tmp->c_seq[i]->state+n_patt,
					  alldata_tmp->c_seq[j]->state+n_patt))) break;
		    }
		  if(j != n_otu) break;
		}
	      
	      if((j == n_otu) && (i == n_otu)) /* all characters at that site are compatible -> the site is invariant */
		{
		  for (j=0; j<n_otu; j++)
		    {
		      alldata_tmp->invar[n_patt] = Assign_State(alldata_tmp->c_seq[j]->state+n_patt);
		      if(alldata_tmp->invar[n_patt] > -1.) break;
		    }
		}
	      else alldata_tmp->invar[n_patt] = -1;
	      
	      alldata_tmp->sitepatt[site] = n_patt;
	      alldata_tmp->wght[n_patt]  += 1;
	      n_patt                     += mod->stepsize;
	    }
	  else
	    {
	      alldata_tmp->sitepatt[site]    = which_patt;
	      alldata_tmp->wght[which_patt] += 1;
	    }
	}
    }
  
  data[0]->len -= n_ambigu;
  
  alldata_tmp->init_len                   = data[0]->len;
  alldata_tmp->crunch_len                 = n_patt;
  for (i=0; i<n_otu; i++)
    alldata_tmp->c_seq[i]->len = n_patt;
  
  if (verbose) {
    printf("\n  . %d patterns found (out of a total of %d sites).\n",n_patt,data[0]->len);
    if((rm_ambigu) && (n_ambigu))
      {
        printf("\n  . Removed %d columns of the alignment as they contain ambiguous characters (e.g., gaps). \n",n_ambigu);
      }
  }

  n_invar=0;
  for (i=0; i<alldata_tmp->crunch_len; i++)
    if(alldata_tmp->invar[i] > -1.)
      n_invar+=(int)alldata_tmp->wght[i];

  if (verbose)
    printf("\n  . %d sites without polymorphism (%.2f%c).\n\n",n_invar,100.*(phydbl)n_invar/data[0]->len,'%');

  alldata_tmp->obs_pinvar = (phydbl)n_invar/data[0]->len;

  n_sites = 0;
  for (i=0; i<alldata_tmp->crunch_len; i++)
    n_sites += alldata_tmp->wght[i];
  if(n_sites != data[0]->len)
    {
      printf("\n Error: file %s at line %d\n",__FILE__,__LINE__);
      exit(EXIT_FAILURE);
    }

//  if(mod->datatype == NT) Get_Base_Freqs(alldata_tmp);
//  else                        Get_AA_Freqs(alldata_tmp);
  Get_AA_Freqs(alldata_tmp);

  alldata = Copy_Cseq(alldata_tmp, alldata_tmp->crunch_len, mod->ns);

  Free_Cseq(alldata_tmp);
  Free_Prefix_Tree(proot,MAX_NAME_LENGTH);

  return alldata;
}

/*********************************************************/

allseq *Compact_CSeq(allseq *data, model *mod)
{
  allseq *alldata;
  int i,j,k,site;
  int n_patt,which_patt;
  int n_otu;

  n_otu = data->n_otu;

  alldata         = (allseq *)mCalloc(1,sizeof(allseq));
  alldata->n_otu  = n_otu;
  alldata->c_seq  = (seq **)mCalloc(n_otu,sizeof(seq *));
  alldata->wght   = (int *)mCalloc(data->crunch_len,sizeof(int));
  alldata->b_frq  = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  alldata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  alldata->invar  = (short int *)mCalloc(data->crunch_len,sizeof(short int));

  alldata->crunch_len = alldata->init_len = -1;
  for (j=0; j<n_otu; j++)
    {
      alldata->c_seq[j]            = (seq *)mCalloc(1,sizeof(seq));
      alldata->c_seq[j]->name      = (char *)mCalloc(MAX_NAME_LENGTH,sizeof(char));
      strcpy(alldata->c_seq[j]->name,data->c_seq[j]->name);
      alldata->c_seq[j]->state     = (char *)mCalloc(data->crunch_len,sizeof(char));
      alldata->c_seq[j]->is_ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
      alldata->c_seq[j]->state[0]  = data->c_seq[j]->state[0];
    }

  n_patt = which_patt =  0;

//  Fors(site,data->crunch_len,mod->stepsize)
  for (site=0; site<data->crunch_len; site+=mod->stepsize)
    {
      if(data->wght[site])
	{
//	  Fors(k,n_patt,mod->stepsize)
          for (k=0; k<n_patt; k+=mod->stepsize)
	    {
	      for (j=0; j<n_otu; j++)
		{
		  if(strncmp(alldata->c_seq[j]->state+k,
			     data->c_seq[j]->state+site,
			     mod->stepsize))
		    break;
		}
	      if(j == n_otu)
		{
		  which_patt = k;
		  break;
		}
	    }
	  
	  /*       /\* TO DO *\/ */
	  /*       k = n_patt; */
	  
	  if(k == n_patt)
	    {
	      for (j=0; j<n_otu; j++)
	        Copy_One_State(data->c_seq[j]->state+site,
				alldata->c_seq[j]->state+n_patt,
				mod->stepsize);

	      for (i=0; i<n_otu; i++)
		{
		  for (j=0; j<n_otu; j++)
		    {
		      if(!(Are_Compatible(alldata->c_seq[i]->state+n_patt,
					  alldata->c_seq[j]->state+n_patt))) break;
		    }
		  if(j != n_otu) break;
		}
	      
	      if((j == n_otu) && (i == n_otu)) 
		{
		  for (j=0; j<n_otu; j++)
		    {
		      alldata->invar[n_patt] = Assign_State(alldata->c_seq[j]->state+n_patt);
		      if(alldata->invar[n_patt] > -1.) break;
		    }
		}
	      else alldata->invar[n_patt] = -1;
	      
	      alldata->wght[n_patt] += data->wght[site];
	      n_patt+=mod->stepsize;
	    }
	  else
	    alldata->wght[which_patt] += data->wght[site];
	  /*       Print_Site(alldata,k,n_otu,"\n",mod->stepsize); */
	}
    }
  
  alldata->init_len   = data->crunch_len;
  alldata->crunch_len = n_patt;
  for (i=0; i<n_otu; i++)
    alldata->c_seq[i]->len = n_patt;

  Get_AA_Freqs(alldata);

  return alldata;
}

/*********************************************************/

allseq *Make_Cseq(int n_otu, int crunch_len, int init_len, char **sp_names)
{
  allseq *alldata;
  int j;

  alldata                        = (allseq *)mCalloc(1,sizeof(allseq));
  alldata->n_otu                 = n_otu;
  alldata->c_seq                 = (seq **)mCalloc(n_otu,sizeof(seq *));
  alldata->b_frq                 = (phydbl *)mCalloc(MAX_NAME_LENGTH,sizeof(phydbl));
  alldata->wght                  = (int *)mCalloc(crunch_len,sizeof(int));
  alldata->ambigu                = (short int *)mCalloc(crunch_len,sizeof(short int));
  alldata->invar                 = (short int *)mCalloc(crunch_len,sizeof(short int));
  alldata->sitepatt              = (int *)mCalloc(  init_len,sizeof(int ));

  alldata->crunch_len = crunch_len;
  alldata->init_len   = init_len;
  alldata->obs_pinvar = .0;

  for (j=0; j<n_otu; j++)
    {
      alldata->c_seq[j]            = (seq *)mCalloc(1,sizeof(seq));
      alldata->c_seq[j]->name      = (char *)mCalloc(MAX_NAME_LENGTH,sizeof(char));
      strncpy(alldata->c_seq[j]->name, sp_names[j], MAX_NAME_LENGTH);
      alldata->c_seq[j]->state     = (char *)mCalloc(crunch_len,sizeof(char));
      alldata->c_seq[j]->is_ambigu = (short int *)mCalloc(crunch_len,sizeof(short int));
    }

  return alldata;
}

/*********************************************************/

void Traverse_Prefix_Tree(int site, int seqnum, int *patt_num, int *n_patt,
	seq **data, model *mod, pnode *n)
{
  int ret_val;

  ret_val = -1;

  if(seqnum == mod->n_otu-1)
    {
      n->weight++;
      if(n->weight == 1)
	{
	  n->num = *n_patt;
	  (*n_patt) += 1;
	}
      (*patt_num) = n->num;
      return;
    }
  else
    {
      int next_state;

      next_state = -1;
      next_state = Assign_State_With_Ambiguity(data[seqnum+1]->state+site);

      if(!n->next[next_state])
	{
	  n->next[next_state] = Create_Pnode(MAX_NAME_LENGTH);
	}
      Traverse_Prefix_Tree(site,seqnum+1,patt_num,n_patt,data,mod,n->next[next_state]);
    }
  return;
}

/*********************************************************/

pnode *Create_Pnode(int size)
{
  pnode *n;
  int i;

  n = (pnode *)mCalloc(1,sizeof(pnode ));
  n->next = (pnode **)mCalloc(size,sizeof(pnode *));
  for (i=0; i<size; i++)
    n->next[i] = NULL;
  n->weight = 0;
  n->num = -1;
  return n;
}

/*********************************************************/

//int Is_Ambigu(char *state, int datatype, int stepsize)
int Is_Ambigu(char *state)
{
//  int i;

//  if(datatype == NT)
//    {
//      For(i,stepsize)
//	{
//	  if(strchr("MRWSYKBDHVNXO?-.",state[i]))
//	    return 1;
//	}
//    }
//  else
//    {
//      if(strchr("X?-.",state[0])) return 1;
//    }
  if (strchr("X?-.",state[0]))
    return 1;
  else
    return 0;
}

/*********************************************************/

void Copy_One_State(char *from, char *to, int state_size)
{
  int i;

  for (i=0; i<state_size; i++)
    to[i] = from[i];
  return;
}

/*********************************************************/

int Are_Compatible(char *statea, char *stateb)
{
  char a,b;
      a = statea[0]; b = stateb[0];
      switch(a)
	{
	case 'A' :
	  {
	    switch(b)
	      {
	      case 'A' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'R' :
	  {
	    switch(b)
	      {
	      case 'R' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'N' :
	  {
	    switch(b)
	      {
	      case 'N' :
	      case 'B' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'B' :
	  {
	    switch(b)
	      {
	      case 'N' :
	      case 'B' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'D' :
	  {
	    switch(b)
	      {
	      case 'D' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'C' :
	  {
	    switch(b)
	      {
	      case 'C' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'Q' :
	  {
	    switch(b)
	      {
	      case 'Q' :
	      case 'Z' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'Z' :
	  {
	    switch(b)
	      {
	      case 'Q' :
	      case 'Z' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'E' :
	  {
	    switch(b)
	      {
	      case 'E' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'G' :
	  {
	    switch(b)
	      {
	      case 'G' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'H' :
	  {
	    switch(b)
	      {
	      case 'H' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'I' :
	  {
	    switch(b)
	      {
	      case 'I' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'L' :
	  {
	    switch(b)
	      {
	      case 'L' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'K' :
	  {
	    switch(b)
	      {
	      case 'K' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'M' :
	  {
	    switch(b)
	      {
	      case 'M' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'F' :
	  {
	    switch(b)
	      {
	      case 'F' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'P' :
	  {
	    switch(b)
	      {
	      case 'P' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'S' :
	  {
	    switch(b)
	      {
	      case 'S' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'T' :
	  {
	    switch(b)
	      {
	      case 'T' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'W' :
	  {
	    switch(b)
	      {
	      case 'W' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'Y' :
	  {
	    switch(b)
	      {
	      case 'Y' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'V' :
	  {
	    switch(b)
	      {
	      case 'V' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'X' :
	  {
	    switch(b)
	      {
	      case 'A':case 'R':case 'N' :case 'B' :case 'D' :
	      case 'C':case 'Q':case 'Z' :case 'E' :case 'G' :
	      case 'H':case 'I':case 'L' :case 'K' :case 'M' :
	      case 'F':case 'P':case 'S' :case 'T' :case 'W' :
	      case 'Y':case 'V': case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	default :
	  {
            printf("\n Error: Please check that characters `%c` and `%c`\n",a,b);
            printf("          correspond to existing amino-acids.\n");
            exit(EXIT_FAILURE);
	    return 0;
	  }
	}
  return 1;
}

/*********************************************************/

int Assign_State(char *c)
{
  int state[3];
  state[0] = state[1] = state[2] = -1;
      switch(c[0])
	{
	case 'A' : {state[0]=0 ; break;}
	case 'R' : {state[0]=1 ; break;}
	case 'N' : {state[0]=2 ; break;}
	case 'D' : {state[0]=3 ; break;}
	case 'C' : {state[0]=4 ; break;}
	case 'Q' : {state[0]=5 ; break;}
	case 'E' : {state[0]=6 ; break;}
	case 'G' : {state[0]=7 ; break;}
	case 'H' : {state[0]=8 ; break;}
	case 'I' : {state[0]=9 ; break;}
	case 'L' : {state[0]=10; break;}
	case 'K' : {state[0]=11; break;}
	case 'M' : {state[0]=12; break;}
	case 'F' : {state[0]=13; break;}
	case 'P' : {state[0]=14; break;}
	case 'S' : {state[0]=15; break;}
	case 'T' : {state[0]=16; break;}
	case 'W' : {state[0]=17; break;}
	case 'Y' : {state[0]=18; break;}
	case 'V' : {state[0]=19; break;}

	case 'B' : {state[0] = 2; break;}
	case 'Z' : {state[0] = 5; break;}
	default  : {state[0]=-1;  break;}
	}
      return state[0];
}

/*********************************************************/

int Assign_State_With_Ambiguity(char *c)
{
  int state[3];
  state[0] = state[1] = state[2] = -1;
      switch(c[0])
	{
	case 'A' : {state[0]= 0; break;}
	case 'R' : {state[0]= 1; break;}
	case 'N' : {state[0]= 2; break;}
	case 'D' : {state[0]= 3; break;}
	case 'C' : {state[0]= 4; break;}
	case 'Q' : {state[0]= 5; break;}
	case 'E' : {state[0]= 6; break;}
	case 'G' : {state[0]= 7; break;}
	case 'H' : {state[0]= 8; break;}
	case 'I' : {state[0]= 9; break;}
	case 'L' : {state[0]=10; break;}
	case 'K' : {state[0]=11; break;}
	case 'M' : {state[0]=12; break;}
	case 'F' : {state[0]=13; break;}
	case 'P' : {state[0]=14; break;}
	case 'S' : {state[0]=15; break;}
	case 'T' : {state[0]=16; break;}
	case 'W' : {state[0]=17; break;}
	case 'Y' : {state[0]=18; break;}
	case 'V' : {state[0]=19; break;}
	case 'B' : {state[0]= 2; break;}
	case 'Z' : {state[0]= 5; break;}
	case 'X' : case '?' : case '-' : {state[0]=20; break;}
	default  : 
	  {
	    printf("\n Error: Unknown character state : %c.\n",state[0]);
	    Exit("Init failed (check the data type).");
	    break;
	  }
	}
      return state[0];
}

/*********************************************************/

void Get_AA_Freqs(allseq *data)
{
  int i,j,k;
  phydbl A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y;
  phydbl fA,fC,fD,fE,fF,fG,fH,fI,fK,fL,fM,fN,fP,fQ,fR,fS,fT,fV,fW,fY;
  int w;
  phydbl sum;
  
  fA = fC = fD = fE = fF = fG = fH = fI = fK = fL =
  fM = fN = fP = fQ = fR = fS = fT = fV = fW = fY = 1./20.;
  
  for (k=0; k<8; k++)
    {
      A = C = D = E = F = G = H = I = K = L =
      M = N = P = Q = R = S = T = V = W = Y = .0;

      for (i=0; i<data->n_otu; i++)
	{
	  for (j=0; j<data->crunch_len; j++)
	    {
	      w = data->wght[j];
	      if(w)
		{
		  switch(data->c_seq[i]->state[j])
		    {
		    case 'A' : A+=w;		break;
		    case 'C' : C+=w;		break;
		    case 'D' : D+=w;		break;
		    case 'E' : E+=w;		break;
		    case 'F' : F+=w;		break;
		    case 'G' : G+=w;		break;
		    case 'H' : H+=w;		break;
		    case 'I' : I+=w;		break;
		    case 'K' : K+=w;		break;
		    case 'L' : L+=w;		break;
		    case 'M' : M+=w;		break;
		    case 'N' : N+=w;		break;
		    case 'P' : P+=w;		break;
		    case 'Q' : Q+=w;		break;
		    case 'R' : R+=w;		break;
		    case 'S' : S+=w;		break;
		    case 'T' : T+=w;		break;
		    case 'V' : V+=w;		break;
		    case 'W' : W+=w;		break;
		    case 'Y' : Y+=w;		break;
		    case 'Z' : Q+=w;		break;
		    case 'X' : case '?' : case 'O' : case '-' :
		      A+=w*fA;
		      C+=w*fC;
		      D+=w*fD;
		      E+=w*fE;
		      F+=w*fF;
		      G+=w*fG;
		      H+=w*fH;
		      I+=w*fI;
		      K+=w*fK;
		      L+=w*fL;
		      M+=w*fM;
		      N+=w*fN;
		      P+=w*fP;
		      Q+=w*fQ;
		      R+=w*fR;
		      S+=w*fS;
		      T+=w*fT;
		      V+=w*fV;
		      W+=w*fW;
		      Y+=w*fY;
		      break;
		    default : break;
		    }
		}
	    }
	}
      sum = (A+C+D+E+F+G+H+I+K+L+M+N+P+Q+R+S+T+V+W+Y);
      fA = A/sum;      fC = C/sum;      fD = D/sum;      fE = E/sum;
      fF = F/sum;      fG = G/sum;      fH = H/sum;      fI = I/sum;
      fK = K/sum;      fL = L/sum;      fM = M/sum;      fN = N/sum;
      fP = P/sum;      fQ = Q/sum;      fR = R/sum;      fS = S/sum;
      fT = T/sum;      fV = V/sum;      fW = W/sum;      fY = Y/sum;
    }

  data->b_frq[0]  = fA;  data->b_frq[1]  = fR;  data->b_frq[2]  = fN;  data->b_frq[3]  = fD;
  data->b_frq[4]  = fC;  data->b_frq[5]  = fQ;  data->b_frq[6]  = fE;  data->b_frq[7]  = fG;
  data->b_frq[8]  = fH;  data->b_frq[9]  = fI;  data->b_frq[10] = fL;  data->b_frq[11] = fK;
  data->b_frq[12] = fM;  data->b_frq[13] = fF;  data->b_frq[14] = fP;  data->b_frq[15] = fS;
  data->b_frq[16] = fT;  data->b_frq[17] = fW;  data->b_frq[18] = fY;  data->b_frq[19] = fV;
  return;
}

/*********************************************************/

allseq *Copy_Cseq(allseq *ori, int len, int ns)
{
  allseq *new;
  int i,j,n_otu;
  char **sp_names;

  n_otu = ori->n_otu;

  sp_names = (char **)mCalloc(n_otu,sizeof(char *));
  for (i=0; i<n_otu; i++)
    {
      sp_names[i] = (char *)mCalloc(MAX_NAME_LENGTH,sizeof(char));
      strncpy(sp_names[i], ori->c_seq[i]->name, MAX_NAME_LENGTH);
    }

  new = Make_Cseq(n_otu, len, ori->init_len, sp_names);
  new->obs_pinvar = ori->obs_pinvar;

  for (i=0; i<ori->init_len; i++)
    new->sitepatt[i] = ori->sitepatt[i];

  for (j=0; j<ori->crunch_len; j++)
    {
      for (i=0; i<ori->n_otu; i++)
	{
	  new->c_seq[i]->state[j]     = ori->c_seq[i]->state[j];
	  new->c_seq[i]->is_ambigu[j] = ori->c_seq[i]->is_ambigu[j];
	}

      new->wght[j]   = ori->wght[j];
      new->ambigu[j] = ori->ambigu[j];
      new->invar[j]  = ori->invar[j];
    }

  for (i=0; i<ori->n_otu; i++)
    {
      new->c_seq[i]->len = ori->c_seq[i]->len;
      strncpy(new->c_seq[i]->name,ori->c_seq[i]->name, MAX_NAME_LENGTH);
    }

  new->init_len           = ori->init_len;
  new->clean_len          = ori->clean_len;
  new->crunch_len         = ori->crunch_len;
  for (i=0; i<ns; i++)
    new->b_frq[i] = ori->b_frq[i];
  new->n_otu              = ori->n_otu;

  for (i=0; i<n_otu; i++)
    free(sp_names[i]);
  free(sp_names);

  return new;
}

/*********************************************************/

void Check_Ambiguities(allseq *data, int datatype, int stepsize)
{
  int i,j;

//  Fors(j,data->crunch_len,stepsize)
  for (j=0; j<data->crunch_len; j+=stepsize)
    {
      for (i=0; i<data->n_otu; i++)
	{
	  data->ambigu[j]              = 0;
	  data->c_seq[i]->is_ambigu[j] = 0;
	}
      for (i=0; i<data->n_otu; i++)
	{
	  if(Is_Ambigu(data->c_seq[i]->state+j))
	    {
	      data->ambigu[j]              = 1;
	      data->c_seq[i]->is_ambigu[j] = 1;
	    }
	}
    }
  return;
}

/*********************************************************/

void Hide_Ambiguities(allseq *data)
{
  int i;
  for (i=0; i<data->crunch_len; i++)
    if(data->ambigu[i])
      data->wght[i] = 0;
  return;
}

/*********************************************************/

int Matinv(double *x, int n, int m, double *space)
{
/* x[n*m]  ... m>=n
*/
   int i,j,k;
   int *irow;
   double ee, t,t1,xmax;
   double det;
   
   ee = 1.0E-10;
   det = 1.0;

   irow = (int *)mCalloc(n,sizeof(int));

   for (i=0; i<n; i++)
     {
       xmax = 0.;
       for (j=i; j<n; j++)
         if (xmax < fabs(x[j*m+i])) 
	   { 
	     xmax = fabs(x[j*m+i]); 
	     irow[i]=j; 
	   }

      det *= xmax;
      if (xmax < ee)   
	{
	  free(irow);
	  printf("\n . Determinant becomes zero at %3d.\n", i+1);
	  printf("\n . Cannot invert the matrix of eigen vectors.\n");
	  return(0);
	}
      if (irow[i] != i) 
	{
	  for (j=0; j<m; j++)
	    {
	      t = x[i*m+j];
	      x[i*m+j] = x[irow[i]*m+j];
	      x[irow[i]*m+j] = t;
	    }
	}
      t = 1./x[i*m+i];
      for (j=0; j<n; j++)
	{
	  if (j == i) continue;
	  t1 = t*x[j*m+i];
	  for (k=0; k<m; k++)
	    x[j*m+k] -= t1*x[i*m+k];
	  x[j*m+i] = -t1;
	}
      for (j=0; j<m; j++)
        x[i*m+j] *= t;
      x[i*m+i] = t;
   }                            /* i  */
   for (i=n-1; i>=0; i--) 
     {
       if (irow[i] == i) continue;
       for (j=0; j<n; j++)
	 {
	   t = x[j*m+i];
	   x[j*m+i] = x[j*m + irow[i]];
	   x[j*m + irow[i]] = t;
	 }
     }
   free(irow);
   return (1);
}

/*********************************************************/

matrix *JC69_Dist(allseq *data, model *mod)
{
  int site,i,j,k;
  phydbl unc_len;
  matrix *mat;
  phydbl **len;

  len = (phydbl **)mCalloc(data->n_otu,sizeof(phydbl *));
  for (i=0; i<data->n_otu; i++)
    len[i] = (phydbl *)mCalloc(data->n_otu,sizeof(phydbl));

  unc_len = .0;

  mat = Make_Mat(data->n_otu);
  Init_Mat(mat,data);

//  Fors(site,data->c_seq[0]->len,mod->stepsize)
  for (site=0; site<data->c_seq[0]->len; site+=mod->stepsize)
    {
      for (j=0; j<data->n_otu-1; j++)
	{
	  for(k=j+1;k<data->n_otu;k++)
	    {
	      if((!Is_Ambigu(data->c_seq[j]->state+site)) &&
		 (!Is_Ambigu(data->c_seq[k]->state+site)))
		{
		  len[j][k]+=data->wght[site];
		  len[k][j]=len[j][k];
		  if(strncmp(data->c_seq[j]->state+site,
			     data->c_seq[k]->state+site,
			     mod->stepsize))
		    mat->P[j][k]+=data->wght[site];
		}
	    }
	}
    }

  for (i=0; i<data->n_otu-1; i++)
    for(j=i+1;j<data->n_otu;j++)
      {
	if(len[i][j])
	  {
	    mat->P[i][j] /= len[i][j];
	  }
	else
	  {
	    mat->P[i][j] = 1.;
	  }

	mat->P[j][i] = mat->P[i][j];

	if((1.-(mod->ns)/(mod->ns-1.)*mat->P[i][j]) < .0)
	  {
	    mat->dist[i][j] = DIST_MAX;
	  }
	else
	  mat->dist[i][j] = -(mod->ns-1.)/(mod->ns)*(phydbl)log(1.-(mod->ns)/(mod->ns-1.)*mat->P[i][j]);

	if(mat->dist[i][j] > DIST_MAX)
	  {
	    mat->dist[i][j] = DIST_MAX;
	  }
	mat->dist[j][i] = mat->dist[i][j];
      }

  for (i=0; i<data->n_otu; i++)
    free(len[i]);
  free(len);

  return mat;
}

/*********************************************************/

matrix *Make_Mat(int n_otu)
{
  matrix *mat;
  int i;

  mat = (matrix *)mCalloc(1,sizeof(matrix));

  mat->n_otu = n_otu;

  mat->P        = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->Q        = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->dist     = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->on_off   = (int *)mCalloc(n_otu,sizeof(int));
  mat->name     = (char **)mCalloc(n_otu,sizeof(char *));
//  mat->tip_node = (node **)mCalloc(n_otu,sizeof(node *));

  for (i=0; i<n_otu; i++)
    {
      mat->P[i]    = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->Q[i]    = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->dist[i] = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->name[i] = (char *)mCalloc(MAX_NAME_LENGTH,sizeof(char));
    }

  return mat;
}

/*********************************************************/

void Init_Mat(matrix *mat, allseq *data)
{
  int i;

  mat->n_otu = data->n_otu;
  mat->r = mat->n_otu;
  mat->curr_int = mat->n_otu;
  mat->method = 1;

  for (i=0; i<data->n_otu; i++)
    {
      strncpy(mat->name[i],data->c_seq[i]->name, MAX_NAME_LENGTH);
      mat->on_off[i] = 1;
    }
  return;
}

/*********************************************************/

void Fill_Missing_Dist(matrix *mat)
{
  int i,j;
  for (i=0; i<mat->n_otu; i++)
    {
      for(j=i+1;j<mat->n_otu;j++)
	{
	  if(i != j)
	    {
	      if(mat->dist[i][j] < .0) 
		{
		  Fill_Missing_Dist_XY(i,j,mat);
		  mat->dist[j][i] = mat->dist[i][j];
		}
	    }
	}
    }
  return;
}

/*********************************************************/

void Fill_Missing_Dist_XY(int x, int y, matrix *mat)
{
  int i,j;
  phydbl *local_mins,**S1S2;
  int cpt;
  int pos_best_estimate;
  phydbl min_crit, curr_crit;

  local_mins = (phydbl *)mCalloc(mat->n_otu*mat->n_otu,sizeof(phydbl ));
  S1S2       = (phydbl **)mCalloc(mat->n_otu*mat->n_otu,sizeof(phydbl *));
  for (i=0; i<mat->n_otu*mat->n_otu; i++)
    S1S2[i] = (phydbl *)mCalloc(2,sizeof(phydbl));

  cpt = 0;
  for (i=0; i<mat->n_otu; i++)
    {
      if((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
	{
	  for (j=0; j<mat->n_otu; j++)
	    {
	      if((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
		{
		  if((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
		    {
		      S1S2[cpt][0] = MIN(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
		      S1S2[cpt][1] = MAX(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
		      cpt++;
		    }
		}
	    }
	}
    }

  Qksort_Matrix(S1S2,0,0,cpt-1);

  local_mins[0] = S1S2[0][1];
  for(i=1;i<cpt;i++) local_mins[i] = (i*local_mins[i-1] + S1S2[i][1])/(phydbl)(i+1);
 
  pos_best_estimate = 0;
  min_crit = curr_crit = MDBL_MAX;

  for (i=0; i<cpt-1; i++)
    {
      if((local_mins[i] < S1S2[i+1][0]) && (local_mins[i] > S1S2[i][0]))
	{
	  curr_crit = Least_Square_Missing_Dist_XY(x,y,local_mins[i],mat);
	  if(curr_crit < min_crit)
	    {
	      min_crit = curr_crit;
	      pos_best_estimate = i;
	    }
	}
    }

  mat->dist[x][y] = local_mins[pos_best_estimate];
  mat->dist[y][x] = mat->dist[x][y];

  for (i=0; i<mat->n_otu*mat->n_otu; i++)
    free(S1S2[i]);
  free(S1S2);
  free(local_mins);

  return;
}

/********************************************************/

void Qksort_Matrix(phydbl **A, int col, int ilo, int ihi)
{
    phydbl pivot;	// pivot value for partitioning array
    int ulo, uhi;	// indices at ends of unpartitioned region
    int ieq;		// least index of array entry with value equal to pivot
    phydbl *tempEntry;	// temporary entry used for swapping

    tempEntry = NULL;

    if (ilo >= ihi) {
	return;
    }
    // Select a pivot value.
    pivot = A[(ilo + ihi)/2][col];
    // Initialize ends of unpartitioned region and least index of entry
    // with value equal to pivot.
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi) {
	if (A[uhi][col] > pivot) {
	    // Here, we can reduce the size of the unpartitioned region and
	    // try again.
	    uhi--;
	} else {
	    // Here, A[uhi] <= pivot, so swap entries at indices ulo and
	    // uhi.
	    tempEntry = A[ulo];
	    A[ulo] = A[uhi];
	    A[uhi] = tempEntry;
	    // After the swap, A[ulo] <= pivot.
	    if (A[ulo][col] < pivot) {
		// Swap entries at indices ieq and ulo.
		tempEntry = A[ieq];
		A[ieq] = A[ulo];
		A[ulo] = tempEntry;
		// After the swap, A[ieq] < pivot, so we need to change
		// ieq.
		ieq++;
		// We also need to change ulo, but we also need to do
		// that when A[ulo] = pivot, so we do it after this if
		// statement.
	    }
	    // Once again, we can reduce the size of the unpartitioned
	    // region and try again.
	    ulo++;
	}
    }
    // Now, all entries from index ilo to ieq - 1 are less than the pivot
    // and all entries from index uhi to ihi + 1 are greater than the
    // pivot.  So we have two regions of the array that can be sorted
    // recursively to put all of the entries in order.
    Qksort_Matrix(A, col, ilo, ieq - 1);
    Qksort_Matrix(A, col, uhi + 1, ihi);
  return;
}

/*********************************************************/

phydbl Least_Square_Missing_Dist_XY(int x, int y, phydbl dxy, matrix *mat)
{
  int i,j;
  phydbl fit;

  fit = .0;
  for (i=0; i<mat->n_otu; i++)
    {
      if((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
	{
	  for (j=0; j<mat->n_otu; j++)
	    {
	      if((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
		{
		  if((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
		    {
		      if(dxy < MIN(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]))
			{
			  fit += pow((mat->dist[i][x] + mat->dist[j][y]) - (mat->dist[i][y] + mat->dist[j][x]),2);
			}
		      else if((mat->dist[i][x] + mat->dist[j][y]) < (mat->dist[i][y] + mat->dist[j][x]))
			{
			  fit += pow(dxy - (mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]),2);
			}
		      else
			{
			  fit += pow(dxy - (mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j]),2);
			}
		    }
		}
	    }
	}
    }
  return fit;
}

/*********************************************************/

void Read_Qmat(double *daa, phydbl *pi, FILE *fp)
{
  int i,j;
  phydbl sum;

  for(i=1;i<20;i++)
    {
      for (j=0; j<19; j++)
	{
	  fscanf(fp,"%lf",&(daa[i*20+j]));
	  daa[j*20+i] = daa[i*20+j];
	  if(j == i-1) break; 
	}
    }

  for(i=1;i<20;i++)
    fscanf(fp,"%lf",pi+i);
  sum = .0;
  for(i=1;i<20;i++)
    sum += pi[i];
  if(fabs(sum - 1.) > 1.E-06)
    Exit("The rate matrix format is incorrect.");
  return;
}

/*********************************************************/
/*
int DiscreteGamma (phydbl *freqK, phydbl *rK,
	phydbl alfa, phydbl beta, int K, int median)
{
  // discretization of gamma distribution with equal proportions
  // in each category
  int i;
  phydbl gap05=1.0/(2.0*K), t, factor=alfa/beta*K, lnga1;

  if(K==1)
    {
      rK[0] = 1.0;
      return 0;
    }

  if (median) 
  {
    for (i=0; i<K; i++)
      rK[i]=PointGamma((i*2.0+1)*gap05, alfa, beta);
    for (i=0,t=0; i<K; i++)
      t+=rK[i];
    for (i=0; i<K; i++)
      rK[i]*=factor/t;
  }
  else {
    lnga1=LnGamma(alfa+1);
    for (i=0; i<K-1; i++)
      freqK[i]=PointGamma((i+1.0)/K, alfa, beta);
    for (i=0; i<K-1; i++)
      freqK[i]=IncompleteGamma(freqK[i]*beta, alfa+1, lnga1);
    rK[0] = freqK[0]*factor;
    rK[K-1] = (1-freqK[K-2])*factor;
    for (i=1; i<K-1; i++)
      rK[i] = (freqK[i]-freqK[i-1])*factor;
  }
  for (i=0; i<K; i++)
    freqK[i]=1.0/K;
  return (0);
}
*/
/*********************************************************/

phydbl PointChi2 (phydbl prob, phydbl v)
{
// returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
//   returns -1 if in error.   0.000002<prob<0.999998
//   RATNEST FORTRAN by
//       Best DJ & Roberts DE (1975) The percentage points of the
//       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
//   Converted into C by Ziheng Yang, Oct. 1993.

   phydbl e=.5e-6, aa=.6931471805, p=prob, g;
   phydbl xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

   if (p<.000002 || p>.999998 || v<=0) return (-1);

   g = LnGamma (v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*(phydbl)log(p)) goto l1;

   ch=pow((p*xx*(phydbl)exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=(phydbl)log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-(phydbl)exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;

l3:
   x=PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*((phydbl)log(1-p)-c*(phydbl)log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=IncompleteGamma (p1, xx, g))<0) {
      printf ("\n Error: IncompleteGamma.\n");
      return (-1);
   }
   p2=p-t;
   t=p2*(phydbl)exp(xx*aa+g+p1-c*(phydbl)log(ch));
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (fabs(q/ch-1) > e) goto l4;

   return (ch);
}

/*********************************************************/

phydbl IncompleteGamma(phydbl x, phydbl alpha, phydbl ln_gamma_alpha)
{
// returns the incomplete gamma ratio I(x,alpha) where x is the upper
//	   limit of the integration and alpha is the shape parameter.
// returns (-1) if in error
// ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
// (1) series expansion     if (alpha>x || x<=1)
// (2) continued fraction   otherwise
// RATNEST FORTRAN by
// Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
// 19: 285-287 (AS32)

   int i;
   phydbl p=alpha, g=ln_gamma_alpha;
   phydbl accurate=1e-8, overflow=1e30;
   phydbl factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);

   factor=(phydbl)exp(p*(phydbl)log(x)-x-g);
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (gin);
}

/*********************************************************/

phydbl PointNormal (phydbl prob)
{
// returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
// returns (-9999) if in error
// Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
// Applied Statistics 22: 96-97 (AS70)
//
// Newer methods:
//   Wichura MJ (1988) Algorithm AS 241: the percentage points of the
//     normal distribution.  37: 477-484.
//   Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
//     points of the normal distribution.  26: 118-121.

   phydbl a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   phydbl a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   phydbl b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   phydbl y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) return (-9999);

   y = sqrt ((phydbl)log(1/(p1*p1)));
   z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   return (p<0.5 ? -z : z);
}

/*********************************************************/

phydbl LnGamma (phydbl alpha)
{
// returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
// Stirling's formula is used for the central polynomial part of the procedure.
// Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
// Communications of the Association for Computing Machinery, 9:684

   phydbl x=alpha, f=0, z;

   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-(phydbl)log(f);
   }
   z = 1/(x*x);
   return  f + (x-0.5)*(phydbl)log(x) - x + .918938533204673
	  + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
	       +.083333333333333)/x;
}

/*********************************************************/

