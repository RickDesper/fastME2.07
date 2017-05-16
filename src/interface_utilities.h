#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "utils.h"

#ifndef MAX_NAME_LENGTH
#define	MAX_NAME_LENGTH	64
#endif

typedef struct __Options {			/* mostly used in 'interface_options.c' */
	char		*I_data_file;		/* required input file containing user data (sequence alignment or distance matrix) */
	char		*I_tree_file;		/* optional input topology file */
	char		*I_seqmat_file;		/* optional input scoring matrix file (only used for SCOREDIST method) */
	char		*O_tree_file;		/* output file for the resulting tree */
	char		*O_mat_file;		/* optional output file for the distance matrix */
	char		*O_stat_file;		/* output file for execution informations */
	char		*O_boot_file;		/* output file for bootstrapped trees */
	FILE		*fpI_data_file;
	FILE		*fpI_tree_file;
	FILE		*fpI_seqmat_file;
	FILE		*fpO_tree_file;
	FILE		*fpO_mat_file;
	FILE		*fpO_stat_file;
	FILE		*fpO_boot_file;
	boolean		use_O_mat_file;
	char		*open_mode;
	boolean		is_interleaved;		/* input sequence file format (interleaved or sequential) */
	int		nb_datasets;		/* number of datasets (alignments or matrices, default is 1)*/
	int		nb_bootstraps;		/* number of replicates when bootstrapping */
	int		input_type;		/* input file data type (MATRIX, DNA, PROTEIN, SCOREDIST, default is MATRIX) */
	int		method;			/* method for building initial tree (BAL, OLS, NJ, UNJ, BIONJ or USER, default is BAL) */
	int		model;			/* evolutionary model (F84, TN93, K2P, JC69, Transversion Only, LogDet, or SCOREDIST, default is NONE which corresponds to distance matrix input) */
	double		gamma;			/* gamma variation across sites */
	int		NNI;			/* type of tree swapping (BAL, OLS, NONE, default is BAL) */
	int		branch;			/* branch lengths to assign to a topology (only when no NNI, BAL or OLS) */
	int		seed;			/* seed for randomization (if doing bootstrapping) */
	boolean		no_gap;			/* remove any site which has a gap in any sequence */
	boolean		use_SPR;		/* SPR postprocessing */
	boolean		use_TBR;		/* TBR postprocessing */
} Options;


void Free_Input (Options *input);
int Filexists (char *filename);
void Getstring_Stdin (char *file_name);
FILE *Openfile (char *filename, char *mode);
void Uppercase (char *ch);
void *mCalloc (int nb, size_t size);
void askOption (char *question, char *c);
boolean testM (char *c);
int getM (char *c);
boolean testN (char *c);
int getN (char *c);
boolean testW (char *c, boolean none);
int getW (char *c);
boolean testD (char *c);
boolean testP (char *c);
boolean testI (char *c);
int getI (char *c);
boolean testF (char *c);
boolean getF (char *c);
int getModel_DNA (char *c);
int getModel_PROTEIN (char *c);
void constantToStr (int c, char *str);

#endif /*UTILITIES_H_*/

