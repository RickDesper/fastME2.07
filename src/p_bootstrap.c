/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                         Functions to update the bootstrap score of a tree ;
;                         comparing it with pseudotrees                     ;
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

#include "p_bootstrap.h"


void boot(char *bestTree, char **bootTree, int n_boot, FILE *treeFile)
{
  arbre *tree, *treeBoot;
  int i, n;
  i = n = 0;
  
  tree = Read_Tree(bestTree);  
  Alloc_Bip(tree);
  Get_Bip(tree->noeud[0],tree->noeud[0]->v[0],tree);

  while (n < n_boot) {
    i=0;
    treeBoot = Read_Tree(bootTree[n]);
    Alloc_Bip(treeBoot);
    Get_Bip(treeBoot->noeud[0],treeBoot->noeud[0]->v[0],treeBoot);
    Compare_Bip(tree,treeBoot);
    Free_Tree(treeBoot);
    n++;
  }
  
  fprintf (treeFile, Write_Tree(tree));
  Free_Tree(tree);

  return;
}

/*********************************************************/

arbre *Read_Tree(char *s_tree)
{
  char **subs;
  int i,n_ext,n_int,n_otu;
  arbre *tree;
  int degree;

  n_otu=0;
  for (i=0; i<(int)strlen(s_tree); i++)
    if(s_tree[i] == ',')
      n_otu++;
  n_otu+=1;

  tree = (arbre *)Make_Tree(n_otu);
//  Init_Tree(tree,tree->n_otu);
  Make_All_Tree_Nodes(tree);
  Make_All_Tree_Edges(tree);
//  Make_Tree_Path(tree);
//  Make_List_Of_Reachable_Tips(tree);

  tree->noeud[n_otu]->num = n_otu;
  tree->noeud[n_otu]->tax = 0;

  subs = Sub_Trees(s_tree,&degree);
  Clean_Multifurcation(subs,degree,3);
  if(degree == 2)
    Unroot_Tree(subs);
  degree = 3;

  tree->has_branch_lengths = 0;
  tree->num_curr_branch_available = 0;
  n_int = n_ext = 0;
  for (i=0; i<degree; i++)
    R_rtree(s_tree,subs[i],tree->noeud[n_otu],tree,&n_int,&n_ext);

  for (i=0; i<NODE_DEG_MAX; i++)
    free(subs[i]);
  free(subs);
  return tree;
}

/*********************************************************/

arbre *Make_Tree(int n_otu)
{
  arbre *tree;
//  int i;
  tree = (arbre *)mCalloc(1,sizeof(arbre ));
  Init_Tree(tree,n_otu);
/*
  tree->t_dir = (int **)mCalloc(2*n_otu-2,sizeof(int *));
  for (i=0; i<2*n_otu-2; i++)
    tree->t_dir[i] = (int *)mCalloc(2*n_otu-2,sizeof(int));
*/
  return tree;
}

/*********************************************************/

void Init_Tree(arbre *tree, int n_otu)
{
  tree->n_otu                     = n_otu;
//  tree->best_tree                 = NULL;
//  tree->old_tree                  = NULL;
//  tree->mat                       = NULL;
  tree->n_root                    = NULL;
  tree->e_root                    = NULL;
//  tree->ps_tree                   = NULL;

//  tree->depth_curr_path           = 0;
  tree->has_bip                   = 0;
//  tree->n_moves                   = 0;
//  tree->n_improvements            = 0;
//  tree->number_of_lk_calls        = 0;
//  tree->number_of_branch_lk_calls = 0;
//  tree->bl_from_node_stamps       = 0;
//  tree->lock_topo                 = 0;
//  tree->ps_page_number            = 0;
//  tree->init_lnL                  = UNLIKELY;
//  tree->best_lnL                  = UNLIKELY;
//  tree->c_lnL                     = UNLIKELY;
//  tree->n_swap                    = 0;

//  tree->n_pattern                 = -1;
//  tree->prop_of_sites_to_consider = 1.;
//  tree->n_root_pos                = -1.;
//  tree->print_labels              = 1;

//  tree->print_boot_val            = 0;
//  tree->print_alrt_val            = 0;
  tree->num_curr_branch_available = 0;
  return;
}

/*********************************************************/

void Make_All_Tree_Nodes(arbre *tree)
{
  int i;
  tree->noeud          = (bnode **)mCalloc(2*tree->n_otu-2,sizeof(bnode *));
//  tree->t_dead_nodes   = (bnode **)mCalloc(2*tree->n_otu-2,sizeof(bnode *));

  for (i=0; i<2*tree->n_otu-2; i++)
    {
      tree->noeud[i] = (bnode *)Make_Node_Light(i);
      if(i < tree->n_otu) tree->noeud[i]->tax = 1;
      else                tree->noeud[i]->tax = 0;
    }
  return;
}

/*********************************************************/

bnode *Make_Node_Light(int num)
{
  bnode *n;
  n        = (bnode *)mCalloc(1,sizeof(bnode));
  n->v     = (bnode **)mCalloc(3,sizeof(bnode *));
  n->l     = (phydbl *)mCalloc(3,sizeof(phydbl));
  n->b     = (bedge **)mCalloc(3,sizeof(bedge *));
  n->name  = (char *)mCalloc(MAX_LABEL_LENGTH,sizeof(char));
//  n->score = (phydbl *)mCalloc(3,sizeof(phydbl));
  Init_Node_Light(n,num);
  return n;
}

/*********************************************************/

void Init_Node_Light(bnode *n, int num)
{
//  n->list_of_reachable_tips = NULL;
  n->num                    = num;
  n->tax                    = -1;
//  n->dist_to_root           = .0;
  return;
}

/*********************************************************/

void Make_All_Tree_Edges(arbre *tree)
{
  int i;

  tree->t_edges      = (bedge **)mCalloc(2*tree->n_otu-3,sizeof(bedge *));
//  tree->t_dead_edges = (bedge **)mCalloc(2*tree->n_otu-3,sizeof(bedge *));

  for (i=0; i<2*tree->n_otu-3; i++)
    tree->t_edges[i] = (bedge *)Make_Edge_Light(NULL,NULL,i);
  return;
}

/*********************************************************/

bedge *Make_Edge_Light(bnode *a, bnode *d, int num)
{
  bedge *b;
  b = (bedge *)mCalloc(1,sizeof(bedge));

  Init_Edge_Light(b,num);

  if(a && b)
    {
      b->left = a;  b->rght = d;
      if(a->tax) {b->rght = a; b->left = d;} /* root */
      /* a tip is necessary on the right side of the edge */

      (b->left == a)?
	(Make_Edge_Dirs(b,a,d)):
	(Make_Edge_Dirs(b,d,a));

      b->l                    = a->l[b->l_r];
      if(a->tax) b->l         = a->l[b->r_l];
      if(b->l < BL_MIN)  b->l = BL_MIN;
      else if(b->l > BL_MAX) b->l = BL_MAX;
      b->l_old                = b->l;
    }
  else
    {
      b->left = NULL;
      b->rght = NULL;
    }

  return b;
}

/*********************************************************/

void Init_Edge_Light(bedge *b, int num)
{
//  b->num                  = num;
  b->bip_score            = 0;
//  b->dist_btw_edges       = .0;
//  b->topo_dist_btw_edges  = 0;
//  b->has_zero_br_len      = 0;
//  b->is_p_lk_l_u2d        = 0;
//  b->is_p_lk_r_u2d        = 0;

//  b->p_lk_left            = NULL;
//  b->p_lk_rght            = NULL;
//  b->Pij_rr               = NULL;
  return;
}

/*********************************************************/

void Make_Edge_Dirs(bedge *b, bnode *a, bnode *d)
{
  int i;

  if(a == b->rght)
    {
      printf("\n. a->num = %3d ; d->num = %3d",a->num,d->num);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      exit(EXIT_FAILURE);
    }
  if(d == b->left)
    {
      printf("\n. a->num = %3d ; d->num = %3d",a->num,d->num);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      exit(EXIT_FAILURE);
    }

  b->l_r = b->r_l = -1;
  for (i=0; i<3; i++)
    {
      if((a->v[i]) && (a->v[i] == d))
	{
	  b->l_r  = i; /* we consider here that 'a' is on the left handside of 'b'*/
	  a->b[i] = b;
	}
      if((d->v[i]) && (d->v[i] == a))
	{
	  b->r_l  = i; /* we consider here that 'd' is on the right handside of 'b'*/
	  d->b[i] = b;
	}
    }

  if(a->tax) {
    b->r_l = 0;
    for (i=0; i<3; i++)
      if(d->v[i]==a) {
        b->l_r = i;
        break;
      }
    }

  b->l_v1 = b->l_v2 = b->r_v1 = b->r_v2 = -1;
  for (i=0; i<3; i++)
    {
      if(b->left->v[i] != b->rght)
	{
	  if(b->l_v1 < 0) b->l_v1 = i;
	  else            b->l_v2 = i;
	}

      if(b->rght->v[i] != b->left)
	{
	  if(b->r_v1 < 0) b->r_v1 = i;
	  else            b->r_v2 = i;
	}
    }
  return;
}

/*********************************************************/
/*
void Make_Tree_Path(arbre *tree)
{
  tree->curr_path = (bnode **)mCalloc(tree->n_otu,sizeof(bnode *));
  return;
}
*/
/*********************************************************/
/*
void Make_List_Of_Reachable_Tips(arbre *tree)
{
  int i,j;

  for (i=0; i<2*tree->n_otu-2; i++)
    {
      tree->noeud[i]->list_of_reachable_tips = (bnode ***)mCalloc(3,sizeof(bnode **));
      tree->noeud[i]->n_of_reachable_tips    = (int *)mCalloc(3,sizeof(int));
      for (j=0; j<3; j++)
	tree->noeud[i]->list_of_reachable_tips[j] = (bnode **)mCalloc(tree->n_otu,sizeof(bnode *));
    }
  return;
}
*/
/*********************************************************/

char **Sub_Trees(char *tree, int *degree)
{
  char **subs;
  int posbeg,posend;
  int i;

  if(tree[0] != '(') {*degree = 1; return NULL;}

  subs=(char **)mCalloc(NODE_DEG_MAX,sizeof(char *));

  for (i=0; i<NODE_DEG_MAX; i++)
    subs[i]=(char *)mCalloc(strlen(tree)+1,sizeof(char));

  posbeg=posend=1;
  (*degree)=0;
  do
    {
      posbeg = posend;
      if(tree[posend] != '(')
	{
	  while((tree[posend] != ',' ) &&
		(tree[posend] != ':' ) &&
		(tree[posend] != '#' ) &&
		(tree[posend] != ')' ))
	    {
	      posend++ ;
	    }
	  posend -= 1;
	}
      else posend=Next_Par(tree,posend);

      while((tree[posend+1] != ',') &&
	    (tree[posend+1] != ':') &&
	    (tree[posend+1] != '#') &&
	    (tree[posend+1] != ')')) {posend++;}

      strncpy(subs[(*degree)],tree+posbeg,posend-posbeg+1);
      strncat(subs[(*degree)],"\0",NODE_DEG_MAX);

      posend += 1;
      while((tree[posend] != ',') &&
	    (tree[posend] != ')')) {posend++;}
      posend+=1;

      (*degree)++;
      if((*degree) == NODE_DEG_MAX)
	{
	  for (i=0; i<(*degree); i++)
	    printf("\n. Subtree %d : %s\n",i+1,subs[i]);
	  printf("\n. The degree of a node cannot be greater than %d\n",NODE_DEG_MAX);
	  exit(EXIT_FAILURE);
	}
    }
  while(tree[posend-1] != ')');

  return subs;
}

/*********************************************************/

int Next_Par(char *s, int pos)
{
  int curr;

  curr=pos+1;

  while(*(s+curr) != ')')
    {
      if(*(s+curr) == '(') curr=Next_Par(s,curr);
      curr++;
    }

  return curr;
}

/*********************************************************/

void Clean_Multifurcation(char **subtrees, int current_deg, int end_deg)
{
  if(current_deg <= end_deg) return;
  else
    {
      char *s_tmp;
      int i;

      s_tmp = (char *)mCalloc(MAX_LINE_LENGTH,sizeof(char));

      strncat(s_tmp,"(",1);
      strncat(s_tmp,subtrees[0],strlen (subtrees[0]));
      strncat(s_tmp,",",1);
      strncat(s_tmp,subtrees[1],strlen (subtrees[1]));
      strncat(s_tmp,")",1);
      free(subtrees[0]);
      subtrees[0] = s_tmp;

      for(i=1;i<current_deg-1;i++)
        strncpy(subtrees[i],subtrees[i+1],strlen (subtrees[i+1]));

      Clean_Multifurcation(subtrees,current_deg-1,end_deg);
    }
  return;
}

/*********************************************************/

void Unroot_Tree(char **subtrees)
{
  char **tmp_sub;
  int degree,i,j;

//  printf("\n. Removing the root...\n");

  tmp_sub = Sub_Trees(subtrees[0],&degree);
  if(degree >= 2)
    {
      strncpy(subtrees[2],subtrees[1],strlen(subtrees[1]));
      Clean_Multifurcation(tmp_sub,degree,2);
      for (j=0; j<2; j++)
        strncpy(subtrees[j],tmp_sub[j],strlen(tmp_sub[j]));
    }
  else
    {
      tmp_sub = Sub_Trees(subtrees[1],&degree);
      strncpy(subtrees[2],subtrees[0],strlen(subtrees[0]));
      Clean_Multifurcation(tmp_sub,degree,2);
      for (j=0; j<2; j++)
        strncpy(subtrees[j],tmp_sub[j],strlen(tmp_sub[j]));
    }

  for (i=0; i<degree; i++)
    free(tmp_sub[i]);
  free(tmp_sub);
  return;
}

/*********************************************************/
/* 'a' in node a stands for ancestor. 'd' stands for descendant */
void R_rtree(char *s_tree_a, char *s_tree_d, bnode *a, arbre *tree, int *n_int, int *n_ext)
{
  int i;
  bnode *d;
  int n_otu = tree->n_otu;

  if(strstr(s_tree_a," "))
    Exit("Tree must not contain a ' ' character.");

  if(s_tree_d[0] == '(')
    {
      char **subs;
      int degree;

      (*n_int)+=1;
      d      = tree->noeud[n_otu+*n_int];
      d->num = n_otu+*n_int;
      d->tax = 0;

      Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]);
      Read_Branch_Length(s_tree_d,s_tree_a,tree);

      for (i=0; i<3; i++)
       {
	 if(!a->v[i])
	   {
	     a->v[i]=d;
	     d->l[0]=a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l;
	     break;
	   }
       }
      d->v[0]=a;

      Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
      tree->num_curr_branch_available++;

      subs=Sub_Trees(s_tree_d,&degree);
      Clean_Multifurcation(subs,degree,2);
      R_rtree(s_tree_d,subs[0],d,tree,n_int,n_ext);
      R_rtree(s_tree_d,subs[1],d,tree,n_int,n_ext);
      for (i=0; i<NODE_DEG_MAX; i++)
        free(subs[i]);
      free(subs);
    }

  else
    {
      int i;

      d      = tree->noeud[*n_ext];
      d->tax = 1;

      Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]);
      Read_Branch_Length(s_tree_d,s_tree_a,tree);
      Read_Node_Name(d,s_tree_d,tree);

      for (i=0; i<3; i++)
	{
	 if(!a->v[i])
	   {
	     a->v[i]=d;
	     d->l[0]=a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l;
	     break;
	   }
	}
      d->v[0]=a;

      Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
      tree->num_curr_branch_available++;

      d->num=*n_ext;
      (*n_ext)+=1;
    }
  return;
}

/*********************************************************/

void Read_Branch_Label(char *s_d, char *s_a, bedge *b)
{
  char *sub_tp;
  char *p;
  int i,pos;

  sub_tp = (char *)mCalloc(MAX_LINE_LENGTH,sizeof(char));

  strncpy(sub_tp,s_d,strlen(s_d));
  strncat(sub_tp,"#",1);
  p = strstr(s_a,sub_tp);
  i = 0;
  b->n_labels = 0;
  if(p)
    {
      if(!(b->n_labels%BLOCK_LABELS))
        Make_New_Edge_Label(b);
      b->n_labels++;

      pos = 0;
      do
	{
	  b->labels[b->n_labels-1][pos] = p[i+strlen(s_d)+1];
	  i++;
	  pos++;
	  if(p[i+strlen(s_d)+1] == '#')
	    {
	      b->labels[b->n_labels-1][pos] = '\0';
	      b->n_labels++;
	      if(!(b->n_labels%BLOCK_LABELS))
	        Make_New_Edge_Label(b);
	      i++;
	      pos=0;
	    }
	}
      while((p[i+strlen(s_d)+1] != ':') &&
	    (p[i+strlen(s_d)+1] != ',') &&
	    (p[i+strlen(s_d)+1] != '('));

      b->labels[b->n_labels-1][pos] = '\0';
    }
/*
  if(p)
    {
      if(b->n_labels == 1)
	printf("\n. Found label '%s' on edge %3d.",b->labels[0],b->num);
      else
	{
	  printf("\n. Found labels ");
	  For(i,b->n_labels) printf("'%s' ",b->labels[i]);
	  printf("on edge %3d.",b->num);
	}
    }
*/
  free(sub_tp);
  return;
}

/*********************************************************/

void Make_New_Edge_Label(bedge *b)
{
  int i;

  b->labels = (char **)realloc(b->labels,(b->n_labels+BLOCK_LABELS)*sizeof(char *));

  if(!b->labels)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      exit (EXIT_FAILURE);
    }
  else
    {
      for(i=b->n_labels;i<b->n_labels+100;i++)
        b->labels[i] = (char *)mCalloc(MAX_LABEL_LENGTH,sizeof(char));
    }
  return;
}

/*********************************************************/

void Read_Branch_Length(char *s_d, char *s_a, arbre *tree)
{
  char *sub_tp;
  char *p;
  bedge *b;
  int i;

  b = tree->t_edges[tree->num_curr_branch_available];

  sub_tp = (char *)mCalloc(MAX_LINE_LENGTH,sizeof(char));

  for (i=0; i<b->n_labels; i++)
    {
      strncat(s_d,"#",1);
      strncat(s_d,b->labels[i],strlen(b->labels[i]));
    }

  strncpy(sub_tp,s_d,strlen(s_d));
  strncat(sub_tp,":",1);
  p = strstr(s_a,sub_tp);
  if(p)
    {
      b->l = atof((char *)p+(int)strlen(sub_tp)+1);
      tree->has_branch_lengths = 1;
    }
  free(sub_tp);
  return;
}

/*********************************************************/

void Connect_One_Edge_To_Two_Nodes(bnode *a, bnode *d, bedge *b, arbre *tree)
{
  int i,dir_a_d;

  dir_a_d = -1;
  for (i=0; i<3; i++)
    if(a->v[i] == d) {dir_a_d = i; break;}

  a->b[dir_a_d] = b;
//  b->num        = tree->num_curr_branch_available;
  b->left       = a;
  b->rght       = d;
  if(a->tax) {b->rght = a; b->left = d;} /* root */
  /* a tip is necessary on the right side of the edge */

  (b->left == a)?
    (Make_Edge_Dirs(b,a,d)):
    (Make_Edge_Dirs(b,d,a));

  b->l                    = a->l[b->l_r];
  if(a->tax) b->l         = a->l[b->r_l];
  if(b->l < BL_MIN)  b->l = BL_MIN;
  else if(b->l > BL_MAX) b->l = BL_MAX;
  b->l_old                = b->l;
  return;
}

/*********************************************************/

void Read_Node_Name(bnode *d, char *s_tree_d, arbre *tree)
{
  int i;

  if(!tree->t_edges[tree->num_curr_branch_available]->n_labels)
    {
      strncpy(d->name,s_tree_d,strlen(s_tree_d));
    }
  else
    {
      i = 0;
      do
	{
	  d->name[i] = s_tree_d[i];
	  i++;
	}
      while(s_tree_d[i] != '#');
      d->name[i] = '\0';
    }
  return;
}

/*********************************************************/

void Alloc_Bip(arbre *tree)
{
  int i,j,k;

  tree->has_bip = 1;

  for (i=0; i<2*tree->n_otu-2; i++)
    {
      tree->noeud[i]->bip_size = (int *)mCalloc(3,sizeof(int));
      tree->noeud[i]->bip_node = (bnode ***)mCalloc(3,sizeof(bnode **));
      tree->noeud[i]->bip_name = (char ***)mCalloc(3,sizeof(char **));
      for (j=0; j<3; j++)
	{
	  tree->noeud[i]->bip_node[j] =
	    (bnode **)mCalloc(tree->n_otu,sizeof(bnode *));

	  tree->noeud[i]->bip_name[j] =
	    (char **)mCalloc(tree->n_otu,sizeof(char *));

          for (k=0; k<tree->n_otu; k++)
	    tree->noeud[i]->bip_name[j][k] =
	    (char *)mCalloc(MAX_LABEL_LENGTH,sizeof(char ));
	}
    }
  return;
}

/*********************************************************/

void Get_Bip(bnode *a, bnode *d, arbre *tree)
{
  if(d->tax)
    {
      d->bip_node[0][0] = d;
      d->bip_size[0]    = 1;
      return;
    }
  else
    {
      int i,j,k;
      int d_a;

      d_a = -1;

      for (i=0; i<3; i++)
	{
	  if(d->v[i] != a)
	    Get_Bip(d,d->v[i],tree);
	  else d_a = i;
	}

      d->bip_size[d_a] = 0;
      for (i=0; i<3; i++)
	if(d->v[i] != a)
	  {
	    for (j=0; j<3; j++)
	      {
		if(d->v[i]->v[j] == d)
		  {
		    for (k=0; k<d->v[i]->bip_size[j]; k++)
		      {
			d->bip_node[d_a][d->bip_size[d_a]] = d->v[i]->bip_node[j][k];
			strncpy(d->bip_name[d_a][d->bip_size[d_a]],d->v[i]->bip_node[j][k]->name,strlen(d->v[i]->bip_node[j][k]->name));
			d->bip_size[d_a]++;
		      }
		    break;
		  }
	      }
	  }

      qsort(d->bip_name[d_a],d->bip_size[d_a],sizeof(char *),Sort_String);

      for (i=0; i<3; i++)
	if(a->v[i] == d)
	  {
	    a->bip_size[i] = 0;
	    for (j=0; j<tree->n_otu; j++)
	      {
	        for (k=0; k<d->bip_size[d_a]; k++)
		  {
		    if(d->bip_node[d_a][k] == tree->noeud[j])
		      break;
		  }

		if(k == d->bip_size[d_a])
		  {
		    a->bip_node[i][a->bip_size[i]] = tree->noeud[j];
		    strncpy(a->bip_name[i][a->bip_size[i]],tree->noeud[j]->name,strlen(tree->noeud[j]->name));
		    a->bip_size[i]++;
		  }
	      }

	    qsort(a->bip_name[i],a->bip_size[i],sizeof(char *),Sort_String);

	    if(a->bip_size[i] != tree->n_otu - d->bip_size[d_a])
	      {
		printf("%d %d \n",a->bip_size[i],tree->n_otu - d->bip_size[d_a]);
		Exit("Problem in counting bipartitions.");
	      }
	    break;
	  }
    }
  return;
}

/*********************************************************/

int Sort_String(const void *a, const void *b)
{
  return(strcmp((*(const char **)(a)), (*(const char **)(b))));
}

/*********************************************************/

void Compare_Bip(arbre *tree1, arbre *tree2)
{
  int i,j,k;
  bedge *b1,*b2;
  char **bip1,**bip2;
  int bip_size;

  for (i=0; i<2*tree1->n_otu-3; i++)
    {
      if((!tree1->t_edges[i]->left->tax) &&
	 (!tree1->t_edges[i]->rght->tax))
	{
	  b1 = tree1->t_edges[i];
          for (j=0; j<2*tree2->n_otu-3; j++)
	    {
	      if((!tree2->t_edges[j]->left->tax) &&
		 (!tree2->t_edges[j]->rght->tax))
		{
		  b2 = tree2->t_edges[j];
		  if(MIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]) ==
		     MIN(b2->left->bip_size[b2->l_r],b2->rght->bip_size[b2->r_l]))
		    {
		      bip_size = MIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]);
		      if(b1->left->bip_size[b1->l_r] == b1->rght->bip_size[b1->r_l])
			{
			  if(b1->left->bip_name[b1->l_r][0][0] < b1->rght->bip_name[b1->r_l][0][0])
			    {
			      bip1 = b1->left->bip_name[b1->l_r];
			    }
			  else
			    {
			      bip1 = b1->rght->bip_name[b1->r_l];
			    }
			}
		      else if(b1->left->bip_size[b1->l_r] < b1->rght->bip_size[b1->r_l])
			{
			  bip1 = b1->left->bip_name[b1->l_r];
			}
		      else
			{
			  bip1 = b1->rght->bip_name[b1->r_l];
			}
		      if(b2->left->bip_size[b2->l_r] == b2->rght->bip_size[b2->r_l])
			{
			  if(b2->left->bip_name[b2->l_r][0][0] < b2->rght->bip_name[b2->r_l][0][0])
			    {
			      bip2 = b2->left->bip_name[b2->l_r];
			    }
			  else
			    {
			      bip2 = b2->rght->bip_name[b2->r_l];
			    }
			}
		      else if(b2->left->bip_size[b2->l_r] < b2->rght->bip_size[b2->r_l])
			{
			  bip2 = b2->left->bip_name[b2->l_r];
			}
		      else
			{
			  bip2 = b2->rght->bip_name[b2->r_l];
			}
		      if(bip_size == 1)
		        Exit("Problem in Compare_Bip.");

                      for (k=0; k<bip_size; k++)
			{
			  if(strncmp(bip1[k],bip2[k],MAX_LABEL_LENGTH)) break;
			}
		      if(k == bip_size)
			{
			  b1->bip_score++;
			  b2->bip_score++;
			  break;
			}
		    }
		}
	    }
	}
    }
  return;
}

/*********************************************************/

void Free_Tree(arbre *tree)
{
  int i,j,k;
  bedge *b;
  bnode *n;

/*
  for (i=0; i<2*tree->n_otu-2; i++)
    free(tree->t_dir[i]);
  free(tree->t_dir);
*/
  if(tree->has_bip)
    {
      for (i=0; i<2*tree->n_otu-2; i++)
	{
	  free(tree->noeud[i]->bip_size);
	  for (j=0; j<3; j++)
	    {
	      free(tree->noeud[i]->bip_node[j]);
	      for (k=0; k<tree->n_otu; k++)
	        free(tree->noeud[i]->bip_name[j][k]);
	      free(tree->noeud[i]->bip_name[j]);
	    }
	  free(tree->noeud[i]->bip_node);
	  free(tree->noeud[i]->bip_name);
	}
    }
//  free(tree->curr_path);

  for (i=0; i<2*tree->n_otu-3; i++)
    {
      b = tree->t_edges[i];
      Free_Edge(b);
    }
  free(tree->t_edges);

  for (i=0; i<2*tree->n_otu-2; i++)
    {
      n = tree->noeud[i];
      Free_bNode(n);
    }
  free(tree->noeud);
/*
  for (i=0; i<tree->n_dead_edges; i++)
    Free_Edge(tree->t_dead_edges[i]);

  for (i=0; i<tree->n_dead_nodes; i++)
    Free_bNode(tree->t_dead_nodes[i]);

  free(tree->t_dead_edges);
  free(tree->t_dead_nodes);
*/
  free(tree);
  return;
}

/*********************************************************/

void Free_bNode(bnode *n)
{
//  int i;

  free(n->b);
  free(n->v);
  free(n->l);
//  free(n->score);
  free(n->name);
/*
  if(n->list_of_reachable_tips)
    {
      for (i=0; i<3; i++)
        free(n->list_of_reachable_tips[i]);
      free(n->list_of_reachable_tips);
      free(n->n_of_reachable_tips);
    }
*/
  free(n);
  return;
}

/*********************************************************/

void Free_Edge_Labels(bedge *b)
{
  int i;
  for (i=0; i<b->n_labels+b->n_labels%BLOCK_LABELS; i++)
    free(b->labels[i]);
  free(b->labels);
  b->labels = NULL;
  return;
}

/*********************************************************/

void Free_Edge(bedge *b)
{
  Free_Edge_Labels(b);
  free(b);
  return;
}

/*********************************************************/

char *Write_Tree(arbre *tree)
{
  char *s;
  int i;

  s=(char *)mCalloc(MAX_LINE_LENGTH,sizeof(char));

  s[0]='(';
/*
  #ifdef PHYML
  tree->n_root = NULL;
  tree->e_root = NULL;
  #endif
*/
/*   if(!tree->n_root) */
/*     { */
      i = 0;
      while((!tree->noeud[tree->n_otu+i]->v[0]) ||
	    (!tree->noeud[tree->n_otu+i]->v[1]) ||
	    (!tree->noeud[tree->n_otu+i]->v[2])) i++;

      R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[0],s,tree);
      R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[1],s,tree);
      R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[2],s,tree);
/*     } */
/*   else */
/*     { */
/*       R_wtree(tree->n_root,tree->n_root->v[0],s,tree); */
/*       R_wtree(tree->n_root,tree->n_root->v[1],s,tree); */
/*     } */

  s[(int)strlen(s)-1]=')';
  s[(int)strlen(s)]=';';

  return s;
}

/*********************************************************/

void R_wtree(bnode *pere, bnode *fils, char *s_tree, arbre *tree)
{
  int i,p;

  p = -1;
  if(fils->tax)
    {
      strncat(s_tree,fils->name,strlen(fils->name));
      if((fils->b[0]) && (fils->b[0]->l != -1))
	{
//	  if(tree->print_labels)
//	    {
	      if(fils->b[0]->n_labels < 10)
	        for (i=0; i<fils->b[0]->n_labels; i++)
	          sprintf(s_tree+(int)strlen(s_tree),"#%s",fils->b[0]->labels[i]);
	      else
		sprintf(s_tree+(int)strlen(s_tree),"#%d_labels",fils->b[0]->n_labels);
//	    }

	  strncat(s_tree,":",1);
 	  sprintf(s_tree+(int)strlen(s_tree),"%f",fils->b[0]->l);
	}
      sprintf(s_tree+(int)strlen(s_tree),",");
   }
  else
    {
      s_tree[(int)strlen(s_tree)]='(';
      for (i=0; i<3; i++)
	{
	  if(fils->v[i] != pere)
	    R_wtree(fils,fils->v[i],s_tree,tree);
	  else p=i;
	}
      s_tree[(int)strlen(s_tree)-1]=')';
      if(fils->b[0]->l != -1)
	{
	  sprintf(s_tree+(int)strlen(s_tree),"%d",fils->b[p]->bip_score);
//	  if(tree->print_labels)
//	    {
	      if(fils->b[p]->n_labels < 10)
	        for (i=0; i<fils->b[p]->n_labels; i++)
	          sprintf(s_tree+(int)strlen(s_tree),"#%s",fils->b[p]->labels[i]);
	      else
		sprintf(s_tree+(int)strlen(s_tree),"#%d_labels",fils->b[p]->n_labels);
//	    }

	  strncat(s_tree,":",1);
	  sprintf(s_tree+(int)strlen(s_tree),"%f",fils->b[p]->l);

	  strncat(s_tree,",",1);
	}
    }
  return;
}

/*********************************************************/


