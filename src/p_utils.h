#ifndef P_UTILS_H
#define P_UTILS_H

#include "graph.h"
#include "interface_utilities.h"


#ifndef MAX
#define MAX(a,b)                     ((a)>(b)?(a):(b))
#endif
#ifndef MIN
#define MIN(a,b)                     ((a)<(b)?(a):(b))
#endif
#ifndef SIGN
#define SIGN(a,b)                    ((b) > 0.0 ? fabs(a) : -fabs(a))
#endif
#ifndef SHFT
#define SHFT(a,b,c,d)                (a)=(b);(b)=(c);(c)=(d);
#endif
#ifndef PointGamma
#define PointGamma(prob,alpha,beta)  PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#endif

#ifndef UNLIKELY
#define UNLIKELY -1.e10
#endif
#ifndef DIST_MAX
#define DIST_MAX 5.00
#endif
#ifndef BL_MIN
#define BL_MIN 1.e-06
#endif
#ifndef BL_MAX
#define BL_MAX 100.0
#endif
#ifndef MDBL_MIN
#define MDBL_MIN 2.225074E-308
#endif
#ifndef MDBL_MAX
#define MDBL_MAX   1.797693E+308
#endif
#ifndef BRENT_ITMAX
#define BRENT_ITMAX 10000
#endif
#ifndef BRENT_ZEPS
#define BRENT_ZEPS 1.e-10
#endif
#ifndef BRENT_CGOLD
#define BRENT_CGOLD 0.3819660
#endif
#ifndef NODE_DEG_MAX
#define NODE_DEG_MAX 50
#endif
#ifndef BLOCK_LABELS
#define BLOCK_LABELS 100
#endif


typedef	double phydbl;

/*********************************************************/

typedef struct __Seq {
  char           *name; /* sequence name */
  int              len; /* sequence length */
  char          *state; /* sequence itself */
  short int *is_ambigu; /* is_ambigu[site] = 1 if state[site] is an ambiguous character. 0 otherwise */
}seq;

/*********************************************************/

typedef struct __AllSeq {
  struct __Seq **c_seq;             /* compressed sequences      */
  phydbl        *b_frq;             /* observed state frequencies */
  short int     *invar;             /* 1 -> states are identical, 0 states vary */
  int            *wght;             /* # of each site in c_seq */
  short int    *ambigu;             /* ambigu[i]=1 is one or more of the sequences at site
				       i display an ambiguous character */
  phydbl    obs_pinvar;
  int            n_otu;             /* number of taxa */
  int        clean_len;             /* uncrunched sequences lenghts without gaps */
  int       crunch_len;             /* crunched sequences lengths */
  int         init_len;             /* length of the uncompressed sequences */
  int        *sitepatt;             /* this array maps the position of the patterns in the
				       compressed alignment to the positions in the uncompressed
				       one */
}allseq;

/*********************************************************/

typedef struct __Model {
//  struct __Optimiz  *s_opt; /* pointer to parameters to optimize */
  struct __Eigen    *eigen;
//  struct __M4       *m4mod;

//  char          *modelname;
//  char  *custom_mod_string; /* string of characters used to define custom models of substitution */

  int              *rr_num; 
  int        *n_rr_per_cat; /* number of rate parameters in each category */
  int            n_diff_rr; /* number of different relative substitution rates in the custom model */
  int         update_eigen; /* update_eigen=1-> eigen values/vectors need to be updated */
  int           whichmodel;
  int             datatype; /* 0->DNA, 1->AA */
  int               n_catg; /* number of categories in the discrete gamma distribution */
  int                invar; /* =1 iff the substitution model takes into account invariable sites */
  int                   ns; /* number of states (4 for ADN, 20 for AA) */
  int              seq_len; /* sequence length */
  int             stepsize; /* stepsize=1 for nucleotide models, 3 for codon models */
  int                n_otu; /* number of taxa */
  int            bootstrap; /* Number of bootstrap replicates (0 : no bootstrap analysis is launched) */
//  int            use_m4mod; /* Use a Makrkov modulated Markov model ? */

  phydbl               *pi; /* states frequencies */
  phydbl      *pi_unscaled; /* states frequencies (unscaled) */

  phydbl    *gamma_r_proba; /* probabilities of the substitution rates defined by the discrete gamma distribution */
  phydbl         *gamma_rr; /* substitution rates defined by the discrete gamma distribution */
  phydbl             kappa; /* transition/transversion rate */
  phydbl            lambda; /* parameter used to define the ts/tv ratios in the F84 and TN93 models */
  phydbl             alpha; /* gamma shapa parameter */
  phydbl            pinvar; /* proportion of invariable sites */
  phydbl         alpha_old;
  phydbl         kappa_old;
  phydbl        lambda_old;
  phydbl        pinvar_old;

  phydbl               *rr; /* relative rate parameters of the GTR or custom model (given by rr_val[rr_num[i]]) */
  phydbl           *rr_val; /* relative rate parameters of the GTR or custom model */
  double         ***Pij_rr; /* matrix of change probabilities */
  phydbl                mr; /* mean rate = branch length/time interval  mr = -sum(i)(vct_pi[i].mat_Q[ii]) */
  phydbl      *user_b_freq; /* user-defined nucleotide frequencies */
  phydbl             *qmat;
  phydbl        *qmat_buff;

  phydbl        *rr_branch; /* relative substitution rates on each branch, for the whole set of sites */
  phydbl      *p_rr_branch; /* corresponding frequencies */
  int          n_rr_branch; /* number of classes */
  phydbl   rr_branch_alpha; /* Shape of the gamma distribution that defines the rr_branch and p_rr_branch values */
}model;

/*********************************************************/

typedef struct __Eigen{
  int              size;
  double             *q; /* matrix which eigen values and vectors are computed */
  double         *space;
  int        *space_int;
  double         *e_val; /* eigen values (vector), real part. */
  double      *e_val_im; /* eigen values (vector), imaginary part */
  double      *r_e_vect; /* right eigen vector (matrix), real part */
  double   *r_e_vect_im; /* right eigen vector (matrix), imaginary part */
  double      *l_e_vect; /* left eigen vector (matrix), real part */
}eigen;

/*********************************************************/

//typedef struct __Optimiz { /* parameters to be optimised (mostly used in 'optimiz.c') */
//  int                 print; /* =1 -> verbose mode  */

//  int             opt_alpha; /* =1 -> the gamma shape parameter is optimised */
//  int             opt_kappa; /* =1 -> the ts/tv ratio parameter is optimised */
//  int            opt_lambda; /* =1 -> the F84|TN93 model specific parameter is optimised */
//  int            opt_pinvar; /* =1 -> the proportion of invariants is optimised */
//  int        opt_state_freq; /* =1 -> the nucleotide frequencies are optimised */
//  int                opt_rr; /* =1 -> the relative rate parameters of the GTR or the customn model are optimised */
//  int         opt_num_param; /* if opt_topo=0 and opt_num_param=1 -> the numerical parameters of the
//				model are optimised. if opt_topo=0 and opt_free_param=0 -> no parameter is
//				optimised */
//  int         opt_cov_delta;
//  int         opt_cov_alpha;
//  int    opt_cov_free_rates;


//  int                opt_bl; /* =1 -> the branch lengths are optimised */
//  int              opt_topo; /* =1 -> the tree topology is optimised */
//  int           topo_search;
//  phydbl            init_lk; /* initial loglikelihood value */
//  int              n_it_max; /* maximum bnumber of iteration during an optimisation step */
//  int              last_opt; /* =1 -> the numerical parameters are optimised further while the
//			       tree topology remains fixed */
//  int     random_input_tree; /* boolean */
//  int         n_rand_starts; /* number of random starting points */
//  int          brent_it_max;
//  int             steph_spr;
//  int       user_state_freq;
//  int   spr_step_after_nnis;


//  phydbl     tree_size_mult; /* tree size multiplier */
//  phydbl  min_diff_lk_local;
//  phydbl min_diff_lk_global;
//  phydbl p_moves_to_examine;
//  int              fast_nni;
//  int                greedy;
//  int          general_pars;
//  int           wim_n_rgrft;
//  int           wim_n_globl;
//  int          wim_max_dist;
//  int           wim_n_optim;
//  int            wim_n_best;
//  int        wim_inside_opt;

//}optimiz;

/*********************************************************/

//typedef struct __M4 {
//  int                  n_h; /* number of hidden states */
//  int                  n_o; /* number of observable states  */
//  int        use_cov_alpha;
//  int         use_cov_free;

//  phydbl          **o_mats; /* set of matrices of substitution rates across observable states */
//  phydbl          *multipl; /* vector of values that multiply each o_mats matrix */
//  phydbl             *o_rr; /* relative rates (symmetric) of substitution between observable states */
//  phydbl             *h_rr; /* relative rates (symmetric) of substitution between hidden states */
//  phydbl            *h_mat; /* matrix that describes the substitutions between hidden states (aka switches) */
//  phydbl             *o_fq; /* equilibrium frequencies for the observable states */
//  phydbl             *h_fq; /* equilibrium frequencies for the hidden states */
//  phydbl    *h_fq_unscaled; /* unscaled equilibrium frequencies for the hidden states */
//  phydbl *multipl_unscaled; /* unscaled  vector of values that multiply each o_mats matrix */

//  phydbl             delta; /* switching rate */
//  phydbl             alpha; /* gamma shape parameter */
//}m4;

/*********************************************************/

typedef struct __Pnode{
  struct __Pnode **next;
  int weight;
  int num;
}pnode;

/*********************************************************/

typedef struct __Matrix { /* mostly used in BIONJ */
  phydbl    **P,**Q,**dist; /* observed proportions of transition, transverion and  distances
			       between pairs of  sequences */
//  arbre              *tree; /* tree... */
  int              *on_off; /* on_off[i]=1 if column/line i corresponds to a node that has not
			       been agglomerated yet */
  int                n_otu; /* number of taxa */
  char              **name; /* sequence names */
  int                    r; /* number of nodes that have not been agglomerated yet */
//  struct __Node **tip_node; /* array of pointer to the leaves of the tree */
  int             curr_int; /* used in the NJ/BIONJ algorithms */
  int               method; /* if method=1->NJ method is used, BIONJ otherwise */
}matrix;

/*********************************************************/

typedef struct __bNode {
  struct __bNode                       **v; /* table of pointers to neighbor nodes. Dimension = 2 x n_otu - 3 */
  struct __bNode               ***bip_node; /* three lists of pointer to tip nodes. One list for each direction */
  struct __bEdge                       **b; /* table of pointers to neighbor branches */

  int                           *bip_size; /* Size of each of the three lists from bip_node */
  int                                 num; /* node number */
  int                                 tax; /* tax = 1 -> external node, else -> internal node */
  char                        ***bip_name; /* three lists of tip node names. One list for each direction */
  char                              *name; /* taxon name (if exists) */
  phydbl                               *l; /* lengths of the (three or one) branche(s) connected to this node */
}bnode;

/*********************************************************/

typedef struct __bEdge {
  /*
    syntax :  (node) [edge]
(left_1) .                   .(right_1)
          \ (left)  (right) /
           \._____________./
           /    [b_fcus]   \
          /                 \
(left_2) .                   .(right_2)

  */

  struct __bNode               *left,*rght; /* node on the left/right side of the edge */
  int         l_r,r_l,l_v1,l_v2,r_v1,r_v2;
  /* these are directions (i.e., 0, 1 or 2): */
  /* l_r (left to right) -> left[b_fcus->l_r] = right */
  /* r_l (right to left) -> right[b_fcus->r_l] = left */
  /* l_v1 (left node to first node != from right) -> left[b_fcus->l_v1] = left_1 */
  /* l_v2 (left node to secnd node != from right) -> left[b_fcus->l_v2] = left_2 */
  /* r_v1 (right node to first node != from left) -> right[b_fcus->r_v1] = right_1 */
  /* r_v2 (right node to secnd node != from left) -> right[b_fcus->r_v2] = right_2 */

  phydbl                                l; /* branch length */
  phydbl                            l_old; /* old branch length */
  int                           bip_score; /* score of the bipartition generated by the corresponding edge
					      bip_score = 1 iif the branch is fond in both trees to be compared,
					      bip_score = 0 otherwise. */
  char                           **labels; /* string of characters that labels the corresponding edge */
  int                            n_labels; /* number of labels */

}bedge;

/*********************************************************/

typedef struct __Arbre {
  struct __bNode                       *n_root; /* root node */
  struct __bEdge                       *e_root; /* edge on which lies the root */
  struct __bNode                       **noeud; /* array of nodes that defines the tree topology */
  struct __bEdge                     **t_edges; /* array of edges */

  int                                 has_bip; /*if has_bip=1, then the structure to compare
						 tree topologies is allocated, has_bip=0 otherwise */
  int                                   n_otu; /* number of taxa */
  int                      has_branch_lengths; /* =1 iff input tree displays branch lengths */
  int               num_curr_branch_available; /*gives the number of the next cell in t_edges that is free to receive a pointer to a branch */

}arbre;

/*********************************************************/

double **Copy_PMat_to_DMat (matrix *mat);
void Free_Seq(seq **d, int n_otu);
void Free_Cseq(allseq *data);
void Free_Prefix_Tree(pnode *n, int size);
void Free_Pnode(pnode *n);
void Free_Eigen(eigen *eigen_struct);
void Free_Mat(matrix *mat);
seq **Get_Seq(FILE *in, boolean interleaved, int *n_otu, int *len, int itype, set *taxa);
seq **Read_Seq_Sequential(FILE *in,int *n_otu, set *taxa);
seq **Read_Seq_Interleaved(FILE *in,int *n_otu, set *taxa);
int Read_One_Line_Seq(seq ***data,int num_otu,FILE *in);
model *Make_Model_Basic();
void Make_Model_Complete(model *mod);
//optimiz *Alloc_Optimiz();
void Set_Defaults_Model(model *mod, double alpha);
//void Set_Defaults_Optimiz(optimiz *s_opt);
eigen *Make_Eigen_Struct(model *mod);
allseq *Compact_Seq(seq **data, model *mod, boolean rm_ambigu);
allseq *Compact_CSeq(allseq *data, model *mod);
allseq *Make_Cseq(int n_otu, int crunch_len, int init_len, char **sp_names);
void Traverse_Prefix_Tree(int site, int seqnum, int *patt_num, int *n_patt,
	seq **data, model *mod, pnode *n);
pnode *Create_Pnode(int size);
int Is_Ambigu(char *state);
void Copy_One_State(char *from, char *to, int state_size);
int Are_Compatible(char *statea, char *stateb);
int Assign_State(char *c);
int Assign_State_With_Ambiguity(char *c);
void Get_AA_Freqs(allseq *data);
allseq *Copy_Cseq(allseq *ori, int len, int ns);
void Check_Ambiguities(allseq *data, int datatype, int stepsize);
void Hide_Ambiguities(allseq *data);
int Matinv(double *x, int n, int m, double *space);
matrix *JC69_Dist(allseq *data, model *mod);
matrix *Make_Mat(int n_otu);
void Init_Mat(matrix *mat, allseq *data);
void Fill_Missing_Dist(matrix *mat);
void Fill_Missing_Dist_XY(int x, int y, matrix *mat);
void Qksort_Matrix(phydbl **A, int col, int ilo, int ihi);
phydbl Least_Square_Missing_Dist_XY(int x, int y, phydbl dxy, matrix *mat);
void Read_Qmat(double *daa, phydbl *pi, FILE *fp);
//int DiscreteGamma (phydbl *freqK, phydbl *rK,
//	phydbl alfa, phydbl beta, int K, int median);
phydbl PointChi2 (phydbl prob, phydbl v);
phydbl IncompleteGamma(phydbl x, phydbl alpha, phydbl ln_gamma_alpha);
phydbl PointNormal (phydbl prob);
phydbl LnGamma (phydbl alpha);

#endif /*P_UTILS_H_*/

