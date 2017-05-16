#include "fastme.h"

int main(int argc, char **argv)
{
  Options *options;
  set *species, *slooper, *species_bk;
  node *addNode;
  
  time_t t_beg,t_end;
  
  //if calling exhaustive search
  //T0 = BME
  //T1 = BME + BNNI
  //T2 = BME + BSPR
  //T3 = BME + BNNI + BSPR
  tree *T0, *T1, *T2, *T3;
  T0 = T1 = T2 = T3 = NULL;
  
  //specific variables for PhyML distance computation from protein alignement
  seq **sequences = NULL;
  model *mod = NULL;
  allseq *alldata = NULL;
  matrix *mat = NULL;
  //specific variables for PhyML bootstraps
  allseq *boot_data = NULL;
  int n_site = 0;
  int *site_num = NULL;
  
  int i, j, numSpecies, scoreMatrixSize;
  int setCounter = 0;
  int repCounter = 0;
  int nniCount = 0;
  int sprCount = 0;
  int tbrCount = 0;
  int seqLength = -1;
  int *filter;
  
  double **D, **A, **scoreMatrix;
  D = NULL;
  A = NULL;
  scoreMatrix = NULL;
  
  char *scoreMatrixAlphabet = NULL;
  
  // DNA sequences from alignment
  char **data = NULL;

  // Strings containing the best tree and the bootstraped trees
  char * bestTree;
  char ** bootTrees;
  
  options = chooseSettings (argc, argv);

  time(&t_beg);

  OpenFiles (options);
  if (SCOREDIST == options->model) {
    scoreMatrixAlphabet = (char *)calloc(PROTEIN_ALPHABET_SIZE,sizeof(char));
    scoreMatrix = loadScoreMatrix(&scoreMatrixSize,options->fpI_seqmat_file,scoreMatrixAlphabet);
  }


//PhyML model for protein distance computation
  if (PROTEIN == options->input_type && SCOREDIST != options->model) {
    //mod = (model *)Make_Model_Basic();
    mod = Make_Model_Basic();
    Set_Defaults_Model(mod, options->gamma);
//    Set_Defaults_Optimiz(mod->s_opt);
    mod->whichmodel = options->model;
    Make_Model_Complete(mod);
  }
//end of PhyML model

  printf ("\n"); fflush (stdout); fflush (stderr);
  sgenrand(options->seed);
  while ( setCounter < options->nb_datasets ) {
    repCounter=0;
    setCounter++;
    species = (set *) malloc(sizeof(set));
    species->firstNode = NULL;
    species->secondNode = NULL;
    printf ("\n#  Analysing dataset n° %d\n", setCounter);
    
    if (setCounter == 1)
      printOptions (options);
    fprintf (options->fpO_stat_file,"Dataset n° %d\n", setCounter);
    
    if (options->nb_datasets > 1) {
      if (setCounter > 1)
        fprintf (options->fpO_tree_file,"\n");
      fprintf (options->fpO_tree_file,"Dataset n° %d\n", setCounter);
      if (options->nb_bootstraps > 0) {
        if (setCounter > 1)
          fprintf (options->fpO_boot_file,"\n");
        fprintf (options->fpO_boot_file,"Dataset n° %d\n", setCounter);
      }
      if (options->use_O_mat_file)
        fprintf (options->fpO_mat_file,"Dataset n° %d\n", setCounter);
    }


  //Get data from Matrix, DNA alignment or Protein alignment
    if (MATRIX == options->input_type)
      D = loadMatrix(options->fpI_data_file,&numSpecies,species);
    else {
      // read sequences from input file
      sequences = Get_Seq (options->fpI_data_file, options->is_interleaved,
                  &numSpecies, &seqLength, options->input_type, species);
      // if DNA sequences then use FastDist data structure
      if (DNA == options->input_type || SCOREDIST == options->model) {
        data = (char **)calloc(numSpecies,sizeof(char *));
        for (i=0; i<numSpecies; i++) {
          data[i] = (char *)calloc(MAX_SEQ_LENGTH,sizeof(char));
          strncpy (data[i], sequences[i]->state, MAX_SEQ_LENGTH);
        }
      }
      // if protein sequences then use PhyML data structure
      else {
        mod->n_otu = numSpecies;
	alldata = Compact_Seq(sequences,mod,options->no_gap);
	Check_Ambiguities(alldata,mod->datatype,mod->stepsize);
        Init_Model(alldata,mod,options->fpI_seqmat_file);
        if (options->nb_bootstraps > 0) {
          site_num = (int *)mCalloc(alldata->init_len,sizeof(int));
          n_site = 0;
          for (i=0; i<alldata->crunch_len; i++)
            for (j=0; j<alldata->wght[i]; j++) {
              site_num[n_site] = i;
              n_site++;
            }
	}
      }
      Free_Seq (sequences, numSpecies);
    }
  //end of Get data

    // mem alloc for best tree and bootstraped trees strings
    bestTree = (char *)calloc(MAX_INPUT_SIZE,sizeof(char));
    bootTrees = (char **)calloc(options->nb_bootstraps,sizeof(char *));
    for (i=0; i<options->nb_bootstraps; i++) {
      bootTrees[i] = (char *)calloc(MAX_INPUT_SIZE,sizeof(char));
    }
    
    species_bk = copySet(species);
    while (repCounter <= options->nb_bootstraps)
      {
	if (repCounter > 0) {
	  species = copySet(species_bk);
	  fprintf(options->fpO_stat_file,"\t##### Boostraps n°%d #####\n", repCounter);
	}
	if ((!verbose) && (options->nb_bootstraps > 1)) {
	  if (repCounter == 1)
	    printf ("\n . Non parametric bootstrap analysis\n\n  [");
	  if (repCounter > 0) {
            printf (".");
            if (0 == (repCounter) % 20){
              printf ("] %d/%d\n", repCounter, options->nb_bootstraps);
              if (repCounter < options->nb_bootstraps)
	        printf ("  [");
	    }
	    else if (repCounter == options->nb_bootstraps)
	      printf ("] %d/%d\n", repCounter, options->nb_bootstraps);
	  }
	}
	fflush(stdout);


	if (MATRIX != options->input_type)
	  { //sequence input
	    if ((! verbose) && (repCounter == 0))
	      printf ("\n . Computing pairwise distances...\n");
	    //DNA sequences => FastDist distance computation
	    if (DNA == options->input_type || SCOREDIST == options->model) {
	      filter = (int *)malloc(seqLength*sizeof(int));
	      if (options->nb_bootstraps > 1)
	        {
		  filter = initZeroArray(seqLength);
		  bootstrapSelect(seqLength,filter,&options->seed);
	        }
	      else
	        filter = initOneArray(seqLength);
	      D = makeDistMatrix (data, numSpecies, seqLength, options->gamma,
	          options->model, options->input_type, filter, scoreMatrix,
	          scoreMatrixSize, scoreMatrixAlphabet, options->no_gap);
	      free (filter);
	    }
	    //Protein sequences => PhyML distance computation
	    else if (PROTEIN == options->input_type) {
	      if (repCounter > 0) {
	        boot_data = p_bootstraps (alldata, mod->ns, site_num);
	        Init_Model(boot_data,mod,options->fpI_seqmat_file);
	        mat = ML_Dist(boot_data,mod);
	      }
	      else
	        mat = ML_Dist(alldata,mod);
	      Fill_Missing_Dist(mat);
	      D = Copy_PMat_to_DMat (mat);
	      Free_Mat(mat);
	    }
	    if (options->use_O_mat_file)
    	      printMatrix(D,numSpecies,species,options->fpO_mat_file);
	  }

	//using the matrix to build a tree
        if ((! verbose) && (repCounter == 0))
	  printf ("\n . Computing tree...\n");
	A = initDoubleMatrix(2*numSpecies-2);  
	switch(options->method)
	  {
	  case USER:
	    T0 = loadNewickTree(options->fpI_tree_file, numSpecies);
	    T0 = detrifurcate(T0);
	    compareSets(T0,species);
	    partitionSizes(T0);
	    break;
	  case OLS:
	    for(slooper = species; NULL != slooper; slooper = slooper->secondNode)
	      {
		addNode = copyNode(slooper->firstNode);
		T0 = GMEaddSpecies(T0,addNode,D,A);
	      }
	    break;
	  case BAL:
	    for(slooper = species; NULL != slooper; slooper = slooper->secondNode)
	      {
		addNode = copyNode(slooper->firstNode);
		T0 = BMEaddSpecies(T0,addNode,D,A);
	      }
	    break;
	  case NJ:
            T0 = bionj (D, species, numSpecies, TRUE);
            compareSets(T0,species);
            partitionSizes(T0);
            break;
          case BIONJ:
            T0 = bionj (D, species, numSpecies, FALSE);
            compareSets(T0,species);
	    partitionSizes(T0);
	    break;
          case UNJ:
            T0 = unj (D, species, numSpecies);
            compareSets(T0,species);
	    partitionSizes(T0);
            break;
	  }

        if ((! verbose) && (repCounter == 0) && (options->NNI != NONE))
	  printf ("\n . Performing NNI...\n");
	switch(options->NNI)
	  {
	  case OLS:
	    if (OLS != options->method)
	      assignAllSizeFields(T0);
	    makeOLSAveragesTable(T0,D,A);
	    NNI(T0,A,&nniCount,options->fpO_stat_file);
	    assignOLSWeights(T0,A);
	    break;
	  case BAL:
	    if (BAL != options->method)
	      makeBMEAveragesTable(T0,D,A);
	    if ((options->use_SPR) || (options->use_TBR))
	      {
		T1 = copyTree(T0);
		bNNI(T1,A,&nniCount,options->fpO_stat_file);
		assignBMEWeights(T1,A);
	      }
	    else
	      {
		T1 = NULL;
		bNNI(T0,A,&nniCount,options->fpO_stat_file);
		assignBMEWeights(T0,A);
	      }
	    break;
	  case NONE:
	    switch(options->branch)
	      {
	      case OLS:
		if (OLS != options->method)
		  assignAllSizeFields(T0);
		makeOLSAveragesTable(T0,D,A);
		assignOLSWeights(T0,A);
		break;
	      case BAL:
		if (BAL != options->method)
		  makeBMEAveragesTable(T0,D,A);
		assignBMEWeights(T0,A);
		break;
	      default:
	        if ((options->method == OLS) || (options->method == BAL))
	          Exit ("Invalid value for branch length.");
	      }
	    break;
	  default:
	    Exit ("Invalid value for NNI type.");
	  }

	if (options->use_SPR) //perform SPRs where desired
	  {
            if ((! verbose) && (repCounter == 0))
	      printf ("\n . Performing SPR...\n");
	    if (options->use_TBR)
	      {
		T2 = copyTree(T0);
		makeBMEAveragesTable(T2,D,A);
		SPR(T2,D,A,&sprCount,options->fpO_stat_file);
		if (NULL != T1)
		  {
		    T3 = copyTree(T1);
		    makeBMEAveragesTable(T3,D,A);
		    SPR(T3,D,A,&sprCount,options->fpO_stat_file);
		  }
		else
		  T3 = NULL; //probably unnecessary
	      }
	    else //we ignore T2, T3 and will have T0 = BSPR, T1 = BNNI + BSPR
	      {
		T2 = T3 = NULL;
		makeBMEAveragesTable(T0,D,A);
		SPR(T0,D,A,&sprCount,options->fpO_stat_file);
		assignBMEWeights(T0,A);
		if (NULL != T1)
		  {
		    makeBMEAveragesTable(T1,D,A);
		    SPR(T1,D,A,&sprCount,options->fpO_stat_file);
		    assignBMEWeights(T1,A);
		  }
	      }
	  } //if options->use_SPR

	if (options->use_TBR)
	  {
            if ((! verbose) && (repCounter == 0))
	      printf ("\n . Performing TBR...\n");

	    makeBMEAveragesTable(T0,D,A);
	    TBR(T0,A,D,&tbrCount,options->fpO_stat_file); //original tree if no SPRs,
				                          //otherwise this is SPR tree
	    assignBMEWeights(T0,A);
	    if (NULL != T1) //needed if there is a BNNI tree or a BSPR tree
	      {
		makeBMEAveragesTable(T1,D,A); //if SPRs have been done to T, this is original tree
		TBR(T1,A,D,&tbrCount,options->fpO_stat_file);
		assignBMEWeights(T1,A);
	      }		    
	    if (NULL != T2) //needed if there are both BNNI tree and BSPR tree
	      {
		makeBMEAveragesTable(T2,D,A);
		TBR(T2,A,D,&tbrCount,options->fpO_stat_file); //original tree
		assignBMEWeights(T2,A);
	      }			   
	    if (NULL != T3) //needed if there are both BNNI tree and BSPR tree
	      {
		makeBMEAveragesTable(T3,D,A);
		TBR(T3,A,D,&tbrCount,options->fpO_stat_file); //original tree
		assignBMEWeights(T3,A);
	      }			   
	  } //if (options->use_TBR)

	if (NULL != T1)
	  weighTree(T1);
	if (NULL != T2)
	  weighTree(T2);
	if (NULL != T3)
	  weighTree(T3);
	if ((NULL != T1) || (NULL != T2) || (NULL != T3))
	  weighTree(T0);

	if ((NULL != T1) && (T0->weight > T1->weight))
	  {
	    freeTree(T0);
	    T0 = T1; //using T0 as the place to store the minimum evolution tree in all cases
	    T1 = NULL;
	  }
	else if (NULL != T1)
	    freeTree(T1);
	
	if ((NULL != T2) && (T0->weight > T2->weight))
	  {
	    freeTree(T0);
	    T0 = T2;
	    T2 = NULL;
	  }
	else if (NULL != T2)
	  freeTree(T2);
	if ((NULL != T3) && (T0->weight > T3->weight))
	  {
	    freeTree(T0);
	    T0 = T3;
	    T3 = NULL;
	  }
	else if (NULL != T3)
	  freeTree(T3);

        if (options->nb_bootstraps > 0) {
          if (repCounter == 0)
            NewickPrintTreeStr (T0, bestTree);
          else {
            NewickPrintTreeStr (T0, bootTrees[repCounter-1]);
            NewickPrintTree (T0, options->fpO_boot_file);
          }
        }
        else
	  NewickPrintTree (T0, options->fpO_tree_file);

	freeMatrix(A,2*numSpecies - 2);
	freeTree(T0);
	T3 = T2 = T1 = T0 = NULL;
	if (NULL != boot_data)
	  Free_Cseq(boot_data);
	
	if ((verbose) && ((options->NNI) || (options->use_SPR)))
	  printf("\n  . Performed %d NNI(s), %d SPR(s) and %d TBR(s) on dataset n°%d.\n", nniCount, sprCount, tbrCount, setCounter);
	fprintf (options->fpO_stat_file, "\tPerformed %d NNI(s), %d SPR(s) and %d TBR(s).\n\n", nniCount, sprCount, tbrCount);
	nniCount = sprCount = tbrCount = 0;
	repCounter++;
      } //end bootstrap loop
    
    //
    if (options->nb_bootstraps > 0)
      boot(bestTree, bootTrees, options->nb_bootstraps, options->fpO_tree_file);
    
    // free mem allocated for best tree and bootstrapped trees strings
    free (bestTree);
    for (i=0; i<options->nb_bootstraps; i++) {
      free (bootTrees[i]);
    }
    free (bootTrees);
      
    fprintf (options->fpO_tree_file, "\n");
    freeMatrix(D,numSpecies);
    freeSet(species);
    if (NULL != data)
      freeCharMatrix(data,numSpecies);
  } //end datasets loop
  printf("\n");fflush (stdout); fflush (stderr);

  if (NULL != site_num)
    free (site_num);
  if (NULL != alldata)
    Free_Cseq(alldata);
  Free_Input (options);
  
  time(&t_end);
  Print_Time_Info(t_beg,t_end);
  exit(EXIT_SUCCESS);
  return 0;
}

/*********************************************************/

allseq *p_bootstraps (allseq *alldata, int ns, int *site_num)
{
  int i, position;
  int init_len = 0;
  phydbl buff;
  
  allseq *boot_data = Copy_Cseq (alldata, alldata->crunch_len, ns);
  
  for (i=0; i<boot_data->crunch_len; i++)
    boot_data->wght[i] = 0;

  for (i=0; i<boot_data->init_len; i++) {
    buff  = rand();
    buff /= RAND_MAX;
    buff *= (phydbl)(alldata->init_len-1.0);
    if(buff-(int)(buff) > 0.5-MDBL_MAX)
      position = (int)(buff)+1;
    else
      position = (int)(buff);
    boot_data->wght[site_num[position]] += 1;
    init_len++;
  }

  if(init_len != alldata->init_len)
    Exit("Problem when copying sequences for bootstrap.");

  Get_AA_Freqs(boot_data);

  return (boot_data);
}

/*********************************************************/

void OpenFiles (Options *options)
{
  options->fpI_data_file = Openfile (options->I_data_file, "r");
  options->fpO_tree_file = Openfile (options->O_tree_file, options->open_mode);
  options->fpO_stat_file = Openfile (options->O_stat_file, options->open_mode);

  if (options->nb_bootstraps > 0)
    options->fpO_boot_file = Openfile (options->O_boot_file, options->open_mode);
  if (options->use_O_mat_file)
    options->fpO_mat_file = Openfile (options->O_mat_file, options->open_mode);
  if (USER == options->method)
    options->fpI_tree_file = Openfile (options->I_tree_file, "r");
  if (SCOREDIST == options->model)
    options->fpI_seqmat_file = Openfile (options->I_seqmat_file, "r");

  return;
}

/*********************************************************/

void printOptions (Options *options)
{
  char *tmp;
  tmp = (char *) mCalloc (MAX_NAME_LENGTH, sizeof(char));
  
  FILE *f = options->fpO_stat_file;
  
  if (options->open_mode[0] == 'w') {
    fprintf (f, "\n - FastME %s - \n\n", VERSION);
    fprintf (f, "\nPapers to be cited:\n");
    fprintf (f, "\nFastME algorithms (balanced and OLS GME, balanced and OLS NNI, balanced and OLS branch length optimization):");
    fprintf (f, "\n\tDesper R., Gascuel O. 2002. Fast and accurate phylogeny reconstruction algorithms based on the minimum-evolution principle.");
    fprintf (f, "\n\tJournal of Computational Biology. 9(5):687-705.");
    fprintf (f, "\nBIONJ algorithm:");
    fprintf (f, "\n\tGascuel O. 1997. BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data.");
    fprintf (f, "\n\tMolecular Biology and Evolution, 14(7):685-695");
    fprintf (f, "\nNJ algorithm:");
    fprintf (f, "\n\tSaitou N., Nei M. 1987. The neighbor-joining method: a new method for reconstructing phylogenetic trees.");
    fprintf (f, "\n\tMolecular Biology and Evolution, 4(4):406-25");
    fprintf (f, "\nUNJ algorithm:");
    fprintf (f, "\n\tGascuel O. 1997. Concerning the NJ algorithm and its unweighted version, UNJ.");
    fprintf (f, "\n\tMathematical Hierarchies and Biology,");
    fprintf (f, "\n\tB. Mirkin, F.R. McMorris, F.S. Roberts and A. Rzetsky (eds.),");
    fprintf (f, "\n\tAmerican Mathematical Society, Providence, 149-170\n\n");
  }

  fprintf (f, "-------------------------------------------------------------------------------\n");
  fprintf (f, "Settings for this run:\n\n");

//----------------------------------------------------------------------//

  constantToStr (options->method, tmp);
  fprintf (f, "  M "
    "                             Initial building method "
    " %-15s \n", tmp);

//----------------------------------------------------------------------//

  constantToStr (options->input_type, tmp);
  fprintf (f, "  I "
    "                                     Input data type "
    " %-15s \n", tmp);

//----------------------------------------------------------------------//

  if (options->input_type != MATRIX)
  {
    constantToStr (options->model, tmp);
    if (options->input_type == DNA)
      fprintf (f, "  E "
        "                              DNA evolutionary model "
        " %-15s \n", tmp);
    else if (options->input_type == PROTEIN)
      fprintf (f, "  E "
        "                          PROTEIN evolutionary model "
        " %-15s \n", tmp);
    
//----------------------------------------------------------------------//

    fprintf (f, "  B "
      "                     Bootstrap: number of replicates "
      " %-15d \n", options->nb_bootstraps);

//----------------------------------------------------------------------//

    fprintf (f, "  G "
      "                   Gamma rate variation across sites "
      " %-15f \n", options->gamma);

//----------------------------------------------------------------------//

    fprintf (f, "  R "
      "                             Remove sites whith gaps "
      " %-15s \n", (options->no_gap ? "yes" : "no"));

//----------------------------------------------------------------------//
			
    fprintf (f, "  O "
      "                   Output calculated distance matrix "
      " %-15s \n", (options->use_O_mat_file ? "yes" : "no"));
			
    fprintf (f, "\n");
  }

//----------------------------------------------------------------------//

  fprintf (f, "  D "
    "                                  Number of datasets "
    " %-15d \n", options->nb_datasets);

//----------------------------------------------------------------------//

  constantToStr (options->NNI, tmp);
  fprintf (f, "  N "
    " Type of tree swapping (NNI) (balanced, OLS or none) "
    " %-15s \n", tmp);

//----------------------------------------------------------------------//

  if (NONE == options->NNI)
  {
    constantToStr (options->branch, tmp);
    fprintf (f, "  W "
      "             Branch lengths assigned to the topology "
      " %-15s \n", tmp);
  }

//----------------------------------------------------------------------//

  fprintf (f, "  S "
    "                                  SPR postprocessing "
    " %-15s \n", (options->use_SPR ? "yes" : "no"));

//----------------------------------------------------------------------//

  fprintf (f, "  T "
    "                                  TBR postprocessing "
    " %-15s \n", (options->use_TBR ? "yes" : "no"));


  fprintf (f, "\n-------------------------------------------------------------------------------\n");
  
  return;
}

/*********************************************************/

void printMatrix(double **D, int size, set *nodes, FILE *ofile)
{
  node *v,*w;
  set *S,*T;
  fprintf(ofile,"%d\n",size);
  for(S=nodes;NULL!=S;S=S->secondNode)
    {
      v=S->firstNode;
      //fprintf(ofile,"%-9.9s",v->label);
      fprintf(ofile,"%s ",v->label);
        for(T=nodes;NULL!=T;T=T->secondNode)
	  {
	    w=T->firstNode;
	    fprintf(ofile,"%8.5lf",D[v->index2][w->index2]);
	  }
	fprintf(ofile,"\n");
    }
  fprintf(ofile,"\n");
  return;
}

/*********************************************************/

void Print_Time_Info(time_t t_beg, time_t t_end)
{
  div_t hour,min;

  hour = div(t_end-t_beg,3600);
  min  = div(t_end-t_beg,60  );
  min.quot -= hour.quot*60;

  printf ("\n . Time used %dh%dm%ds\n\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
  
  return;
}

