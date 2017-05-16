#include "interface_options.h"

Options *chooseSettings(int argc, char **argv)
{
	boolean ask;
	char *tmp;
	char choix;
	int n_trial = 0;
 	Options *input = (Options *) mCalloc (1, sizeof (Options));
//	putchar('\n');

	Set_Defaults_Input (input);

	if (argc == 1)
		Get_Input_Interactive (input);
	else
		Get_Input_CommandLine (input, argc, argv);

	if (! input->I_data_file)
	  Exit ("You must provide an input file.");

	/* ouput files default values */
	if ((NULL != input->O_tree_file) && (strlen (input->O_tree_file) == 0)) {
		strncpy (input->O_tree_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 17);
		strncat (input->O_tree_file, "_fastme_tree.txt", 16);
	}
	if ((NULL != input->O_stat_file) && (strlen (input->O_stat_file) == 0)) {
		strncpy (input->O_stat_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 17);
		strncat (input->O_stat_file, "_fastme_stat.txt", 16);
	}
	if ((input->use_O_mat_file) && (NULL != input->O_mat_file) && (strlen (input->O_mat_file) == 0)) {
		strncpy (input->O_mat_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 16);
		strncat (input->O_mat_file, "_fastme_mat.txt", 15);
	}	
	if ((input->input_type == MATRIX) && (input->nb_bootstraps > 0)) {
		Warning ("Bootstraps can only be used when inputting aligned sequences.");
		input->nb_bootstraps = 0;
	}
	if ((input->nb_bootstraps > 0) && (NULL != input->O_boot_file) && (strlen (input->O_boot_file) == 0)) {
		strncpy (input->O_boot_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 17);
		strncat (input->O_boot_file, "_fastme_boot.txt", 16);
	}

	ask = FALSE;
	tmp = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	if (Filexists (input->O_tree_file)) {
		strncat (tmp, basename (input->O_tree_file), strlen (basename (input->O_tree_file)));
		ask = TRUE;
	}
	if (Filexists (input->O_stat_file)) {
		if (ask) strncat (tmp, "'\n  '", 5);
		strncat (tmp, basename (input->O_stat_file), strlen (basename (input->O_stat_file)));
		ask = TRUE;
	}
	if (Filexists (input->O_boot_file)) {
		if (ask) strncat (tmp, "'\n  '", 5);
		strncat (tmp, basename (input->O_boot_file), strlen (basename (input->O_boot_file)));
		ask = TRUE;
	}
	if ((input->use_O_mat_file) && (Filexists (input->O_mat_file))) {
		if (ask) strncat (tmp, "'\n  '", 5);
		strncat (tmp, basename (input->O_mat_file), strlen (basename (input->O_mat_file)));
		ask = TRUE;
	}
	if (ask) {
		printf ("\n The file(s):\n  '%s'\n already exist(s).\n\n", tmp);
		printf (" Do you want to Replace or Append ? > ");
		n_trial = 0;
		do
		{
			if (n_trial > 0)
				printf ("Please type R or A > ");
			scanf ("%c", &choix);
			if (choix == '\n')
				choix = 'X'; 
			else
				getchar();
			Uppercase (&choix);
	  		if (++n_trial > 10)
	  			Exit("\n");
		}
		while ((choix != 'R') && (choix != 'A'));
		if (choix == 'R')
			strncpy (input->open_mode, "w", 3);
		else
			strncpy (input->open_mode, "a", 3);
	}
	free (tmp);
	return input;
}

/*********************************************************/

void Usage()
{
	printf (BOLD"NAME\n"
		FLAT"\tFastME - A distance based phylogeny reconstruction algorithm.\n\n"
		FLAT"\tFastME showed better topological accuracy than NJ,\n"
		FLAT"\tBIONJ, WEIGHBOR and FITCH, in all evolutionary\n"
		FLAT"\tconditions we tested, which include large range\n"
		FLAT"\tdeviations from molecular clock and substitution rates.\n\n"
		FLAT"\tRichard Desper and Olivier Gascuel,\n"
		FLAT"\tJournal of Computational Biology 19(5), 687-705, 2002.\n"
		FLAT"\tMolecular Biology and Evolution 21(3), 587-598, 2004.\n"
		FLAT"\tPlease cite these papers if you use this software in your publications.\n");

	printf(BOLD"\nSYNOPSIS\n"
		BOLD"\tfastme"
		FLAT"  ["BOLD"-m "LINE"method"FLAT"]"
		FLAT"  ["BOLD"-i "LINE"input data file"FLAT"]"
		FLAT"  ["BOLD"-T "LINE"input tree file"FLAT"]"
		FLAT"  ["BOLD"-f "LINE"input data file format"FLAT"]"
		FLAT"\n\t["BOLD"-D "LINE"model"FLAT"|"BOLD"-P "LINE"model"FLAT"|"BOLD"-S "LINE"socring matrix file"FLAT"]"
		FLAT"\n\t["BOLD"-o "LINE"output tree file"FLAT"]"
		FLAT"  ["BOLD"-I "LINE"output information file"FLAT"]"
		FLAT"\n\t["BOLD"-B "LINE"output bootstrap trees file"FLAT"]"
		FLAT"  ["BOLD"-O "LINE"output matrix file"FLAT"]"
		FLAT"\n\t["BOLD"-g "LINE"gamma"FLAT"]"
		FLAT"  ["BOLD"-r]"
		FLAT"  ["BOLD"-d "LINE"datasets"FLAT"]"
		FLAT"  ["BOLD"-b "LINE"replicates"FLAT"]"
		FLAT"  ["BOLD"-z "LINE"seed"FLAT"]"
		FLAT"\n\t["BOLD"-n "LINE"NNI"FLAT"]"
		FLAT"  ["BOLD"-w "LINE"branch"FLAT"]"
		FLAT"  ["BOLD"-s"FLAT"]"
		FLAT"  ["BOLD"-t"FLAT"]"
		FLAT"  ["BOLD"-v"FLAT"]"
		FLAT"  ["BOLD"-h"FLAT"] \n");

	printf(FLAT"\n\tYou can use fastme with no arguments, in this case change the value of\n"
		FLAT"\ta parameter by typing its corresponding character as shown on screen.\n");

	printf(BOLD"\nOPTIONS\n"
		BOLD"\n\t-b "LINE"replicates"BOLD", --bootstrap="LINE"replicates"
		FLAT"\n\t\tUse this option to indicate the number of "LINE"replicates"BOLD" fastme "FLAT"will"
		FLAT"\n\t\tdo for bootstrapping. By default, the number of "LINE"replicates"FLAT" is 0."
		FLAT"\n\t\tOnly helpful when the input data file contains sequences alignment(s).\n"
		
		BOLD"\n\t-B "LINE"output bootstrap trees file"BOLD", --output_boot="LINE"output bootstrap trees file"
		FLAT"\n\t\tUse this option if you want "BOLD"fastme "FLAT"to write bootstrap trees"
		FLAT"\n\t\tin the "LINE"bootstrap trees file"FLAT".\n"
		
		BOLD"\n\t-d "LINE"datasets"BOLD", --datasets="LINE"datasets"
		FLAT"\n\t\tUse this option to indicate the number of "LINE"datasets"FLAT" in your input"
		FLAT"\n\t\tdata file. By default, the number of "LINE"datasets"FLAT" is 1.\n"
		
		BOLD"\n\t-D "LINE"model"BOLD", --DNA="LINE"model"
		FLAT"\n\t\tUse this option if your input data file contains DNA sequences alignment(s)."
		FLAT"\n\t\tYou may also indicate the evolutionary "LINE"model"FLAT" which can be choosen from:"
		FLAT"\n\t\t"BOLD"(F)84"FLAT" (default), "BOLD"(J)C69"FLAT", "BOLD"(K)2P"FLAT", "BOLD"(L)ogDet"FLAT", "BOLD"(T)N93"FLAT" or "BOLD"T(r)ansversion only"FLAT".\n"
		
		BOLD"\n\t-f "LINE"format"BOLD", --format="LINE"format"
		FLAT"\n\t\tUse this option if your input data file contains sequences alignment(s)."
		FLAT"\n\t\tYou may indicate if the "LINE"format"FLAT" is "BOLD"(I)nterleaved"FLAT" (default) or "BOLD"(S)equential"FLAT".\n"
		
		BOLD"\n\t-g "LINE"gamma"BOLD", --gamma="LINE"gamma"
		FLAT"\n\t\tThis option sets the "LINE"gamma"FLAT" value for gamma variation across sites."
		FLAT"\n\t\tOnly helpful when the input data file contains sequences alignment(s).\n"
		
		BOLD"\n\t-h, --help"
		FLAT"\n\t\tDisplay this usage.\n"
		
		BOLD"\n\t-i "LINE"input data file"BOLD", --input_data="LINE"input data file"
		FLAT"\n\t\tThe "LINE"input data file"FLAT" contains sequence alignment(s)"
		FLAT"\n\t\tor a distance matrix(ces).\n"
		
		BOLD"\n\t-I "LINE"output information file"BOLD", --output_info="LINE"output information file"
		FLAT"\n\t\tUse this option if you want "BOLD"fastme "FLAT"to write information"
		FLAT"\n\t\tabout its execution in the "LINE"output information file"FLAT".\n"
		
		BOLD"\n\t-m "LINE"method"BOLD", --method="LINE"method"
		FLAT"\n\t\t"BOLD"fastme "FLAT"computes a tree using a distance algorithm."
		FLAT"\n\t\tYou may choose this "LINE"method"FLAT" from:"
		FLAT"\n\t\t"BOLD"(b)alanced_GME"FLAT" (default), "BOLD"(O)LS_GME"FLAT", "BOLD"B(I)ONJ"FLAT","
		FLAT"\n\t\t"BOLD"(N)J"FLAT" or "BOLD"(U)NJ"FLAT".\n"
		
		BOLD"\n\t-n "LINE"NNI"BOLD", --NNI="LINE"NNI"
		FLAT"\n\t\tThis option sets the type of tree swapping ("LINE"NNI"FLAT")"
		FLAT"\n\t\tYou may choose the "LINE"NNI"FLAT" type from:"
		FLAT"\n\t\t"BOLD"(b)alanced_NNI"FLAT" (default), "BOLD"(O)LS_NNI"FLAT" or "BOLD"(n)one"FLAT".\n"
		
		BOLD"\n\t-o "LINE"output tree file"BOLD", --output_tree="LINE"output tree file"
		FLAT"\n\t\t"BOLD"fastme "FLAT"will write the infered tree into the "LINE"output tree file"FLAT".\n"
		
		BOLD"\n\t-O "LINE"output matrix file"BOLD", --output_matrix="LINE"output matrix file"
		FLAT"\n\t\tUse this option if you want "BOLD"fastme "FLAT"to write the distance"
		FLAT"\n\t\tmatrix in the "LINE"output matrix file"FLAT".\n"
		
		BOLD"\n\t-P "LINE"model"BOLD", --protein="LINE"model"
		FLAT"\n\t\tUse this option if your input data file contains protein sequences alignment(s)."
		FLAT"\n\t\tYou may also indicate the evolutionary "LINE"model"FLAT" which can be choosen from:"
		FLAT"\n\t\t"BOLD"(C)pRev"FLAT", "BOLD"(D)CMut"FLAT", "BOLD"Day(h)off"FLAT", "BOLD"(J)TT"FLAT", "BOLD"(L)G"FLAT" (default),"
		FLAT"\n\t\t"BOLD"(M)tREV"FLAT", "BOLD"(R)tREV"FLAT", "BOLD"(V)T"FLAT" or "BOLD"(W)AG"FLAT".\n"

		BOLD"\n\t-r, --remove_gap"
		FLAT"\n\t\tUse this option to completely remove any site which has a gap in"
		FLAT"\n\t\tany sequence. By default, "BOLD"fastme "FLAT"is doing pairwise deletion of gaps.\n"
		
		BOLD"\n\t-s, --SPR"
		FLAT"\n\t\tUse this option to do "LINE"SPR"FLAT" postprocessing.\n"
		
		BOLD"\n\t-S "LINE"scoring matrix file"BOLD", --input_scoring="LINE"scoring matrix file"
		FLAT"\n\t\tUse this option if you want "BOLD"fastme "FLAT"to compute distance from alignment"
		FLAT"\n\t\tusing data in the "LINE"scoring matrix file"FLAT".\n"
		
		BOLD"\n\t-t, --TBR"
		FLAT"\n\t\tUse this option to do "LINE"TBR"FLAT" postprocessing.\n"
		
		BOLD"\n\t-T "LINE"input tree file"BOLD", --input_topology="LINE"input tree file"
		FLAT"\n\t\t"BOLD"fastme "FLAT"may use an existing topology available in the "LINE"input tree file"
		FLAT"\n\t\twhich corresponds to the input dataset.\n"
		
		BOLD"\n\t-v, --verbose"
		FLAT"\n\t\tUse this option to turn "BOLD"fastme"FLAT" in verbose mode.\n"
		
		BOLD"\n\t-w "LINE"branch"BOLD", --branch_length="LINE"branch"
		FLAT"\n\t\tUse this option to indicate the "LINE"branch"FLAT" length to assign to the tree."
		FLAT"\n\t\tYou may choose the "LINE"branch"FLAT" length from: "
		BOLD"(b)alanced"FLAT" (default), "BOLD"(O)LS"
		FLAT"\n\t\tor "BOLD"(n)one"FLAT". "BOLD"(n)one "FLAT"is only available with BIONJ, NJ or UNJ."
		FLAT"\n\t\tOnly helpful when not doing NNI.\n"
		
		BOLD"\n\t-z "LINE"seed"BOLD", --seed="LINE"seed"
		FLAT"\n\t\tUse this option to initialize randomization with "LINE"seed"FLAT" value."
		FLAT"\n\t\tOnly helpful when bootstrapping.\n");
		
	printf(FLAT"\n%s\n", VERSION);
	exit (EXIT_SUCCESS);
	return;
}

/*********************************************************/

void Set_Defaults_Input (Options *input)
{
	input->I_data_file	= (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->I_tree_file	= (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->I_seqmat_file	= (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->O_tree_file	= (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->O_mat_file	= (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->O_stat_file	= (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->O_boot_file	= (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->open_mode	= (char *) mCalloc (3, sizeof(char));
	
	input->fpI_data_file	= NULL;
	input->fpI_tree_file	= NULL;
	input->fpI_seqmat_file	= NULL;
	input->fpO_tree_file	= NULL;
	input->fpO_mat_file	= NULL;
	input->fpO_stat_file	= NULL;
	input->fpO_boot_file	= NULL;
	
	strncpy (input->open_mode, "w", 3);
	input->use_O_mat_file		= FALSE;
	input->is_interleaved		= TRUE;
	input->nb_datasets		= 1;
	input->nb_bootstraps		= 0;
	input->input_type		= MATRIX;
	input->method			= BAL;	
	input->model			= NONE;
	input->gamma			= 1.0;
	input->NNI			= BAL;
	input->branch			= BAL;
	input->seed			= time (NULL);
	input->no_gap			= FALSE;
	input->use_SPR			= FALSE;
	input->use_TBR			= FALSE;
	verbose				= FALSE;

	return;
}

/*********************************************************/

Options *Get_Input (int argc, char **argv)
{
 	Options *input = (Options *) mCalloc (1, sizeof (Options));
	putchar('\n');

	Set_Defaults_Input (input);

	if (argc == 1)
		Get_Input_Interactive (input);
	else
		Get_Input_CommandLine (input, argc, argv);

	return input;
}

/*********************************************************/

void Get_Input_CommandLine (Options *input, int argc, char **argv)
{
       int c;
       while (1)
         {
           static struct option long_options[] =
             {
               // These options don't set a flag.
               // We distinguish them by their indices.
               {"input_data",     required_argument, 0, 'i'},	// input file (sequence alignment or distance matrix)
               {"input_topology", required_argument, 0, 'T'},	// input topology file
               {"input_scoring",  required_argument, 0, 'S'},	// input scoring matrix file (sets input data type to SCOREDIST)
               {"output_tree",    required_argument, 0, 'o'},	// output file for the resulting tree
               {"output_matrix",  required_argument, 0, 'O'},	// output file for the distance matrix
               {"output_info",    required_argument, 0, 'I'},	// output file for execution informations
               {"output_boot",    required_argument, 0, 'B'},	// output file for bootstrap trees
               {"format",         required_argument, 0, 'f'},	// input sequence file format (interleaved or sequential)
               {"datasets",       required_argument, 0, 'd'},	// number of datasets (alignments or matrices)
               {"bootstrap",      required_argument, 0, 'b'},	// number of replicates when bootstrapping
               {"method",         required_argument, 0, 'm'},	// method for building initial tree
               {"NNI",            required_argument, 0, 'n'},	// type of NNI
               {"branch_length",  required_argument, 0, 'w'},	// branch lengths to assign to a topology
               {"seed",           required_argument, 0, 'z'},	// seed for randomization
               {"SPR",            no_argument,       0, 's'},	// SPR postprocessing
               {"TBR",            no_argument,       0, 't'},	// TBR postprocessing
               {"verbose",        no_argument,       0, 'v'},
               {"help",           no_argument,       0, 'h'},
               {"dna",            required_argument, 0, 'D'},	// evolutionary model for DNA sequence input (sets input data type to DNA)
               {"protein",        required_argument, 0, 'P'},	// sets input data type to PROTEIN
               {"gamma",          required_argument, 0, 'g'},	// gamma variation across sites
               {"remove_gap",     no_argument,       0, 'r'},	// remove any site which has a gap in any sequence
               {0, 0, 0, 0}
             };
           // getopt_long stores the option index here.
           int option_index = 0;
     
           c = getopt_long (argc, argv, "i:T:S:o:O:I:B:f:d:b:m:n:w:z:stvhD:P:g:r", long_options, &option_index);
     
           // Detect the end of the options.
           if (c == -1)
             break;
     
           switch (c)
             {
             case 0:
               printf ("option %s", long_options[option_index].name);
               if (optarg)
                 printf (" with arg %s", optarg);
               printf ("\n");
               break;
     
             case 'B':
               // output bootstrap trees file
               if (NULL != optarg)
                 strncpy (input->O_boot_file, optarg, MAX_FILE_NAME_LENGTH);
               break;

             case 'b':
               // number of replicates when bootstrapping
               input->nb_bootstraps = atoi (optarg);
	       if (input->nb_bootstraps < 1) {
                 fprintf (stderr, "\n. Error: -b option: '%d' invalid value for number of bootstrap replicates.\n", input->nb_bootstraps);
                 exit(EXIT_FAILURE);
               }
               break;
     
             case 'd':
               // number of datasets (alignments or matrices)
               input->nb_datasets = atoi (optarg);
	       if (input->nb_datasets < 1) {
                 fprintf (stderr, "\n. Error: -d option: '%d' invalid value for number of datasets.\n", input->nb_datasets);
                 exit(EXIT_FAILURE);
               }
               break;
     
             case 'D':
               // evolutionary model for DNA sequence input (sets input data type to DNA)
               input->input_type = DNA;
               if (NULL != optarg) {
                 if (testD (optarg))
                   input->model = getModel_DNA (optarg);
                 else {
                   fprintf (stderr, "\n. Error: -D option: '%s' invalid evolutionary model.\n", optarg);
                   exit(EXIT_FAILURE);
                 }
               }
               else
                 input->model = LOGDET;
               break;
     
             case 'f':
               // input sequence file format (interleaved or sequential)
               if (testF (optarg))
                 input->is_interleaved = getF(optarg);
               else {
                 fprintf (stderr, "\n. Error: -f option: '%s' invalid format.\n", optarg);
                 exit(EXIT_FAILURE);
               }
               break;
     
             case 'g':
               // gamma variation across sites
               input->gamma = atof(optarg);
               if (input->gamma <= 0) {
	         fprintf(stderr,"\n. Error: -g option: '%s' invalid value for gamma rate variation\n", optarg);
	         exit(EXIT_FAILURE);
	       }
               break;
     
             case 'h':
               Usage();
               break;
     
             case 'i':
               // input file (sequence alignment or distance matrix)
               if (Filexists (optarg))
                 strncpy (input->I_data_file, optarg, MAX_FILE_NAME_LENGTH);
               else {
                 fprintf (stderr, "\n. Error: -i option: '%s' file does not exist.\n", optarg);
                 exit(EXIT_FAILURE);
               }
               break;
     
             case 'I':
               // output file for execution informations
               if (NULL != optarg)
                 strncpy (input->O_stat_file, optarg, MAX_FILE_NAME_LENGTH);
               break;
     
             case 'm':
               // method for building initial tree
               if (testM (optarg))
                 input->method = getM (optarg);
               else {
                 fprintf (stderr, "\n. Error: -m option: '%s' invalid method.\n", optarg);
                 exit(EXIT_FAILURE);
               }
               break;
     
             case 'n':
               // type of NNI
               if (testN (optarg))
                 input->NNI = getN (optarg);
               else {
                 fprintf (stderr, "\n. Error: -n option: '%s' invalid NNI type.\n", optarg);
                 exit(EXIT_FAILURE);
               }
               break;
     
             case 'o':
               // output file for the resulting tree
               if (NULL != optarg)
                 strncpy (input->O_tree_file, optarg, MAX_FILE_NAME_LENGTH);
               break;
             
             case 'O':
               // output file for the distance matrix
               if (NULL != optarg)
                 strncpy (input->O_mat_file, optarg, MAX_FILE_NAME_LENGTH);
               input->use_O_mat_file = TRUE;
               break;
     
             case 'P':
               input->input_type = PROTEIN;
               // evolutionary model for PROTEIN sequence input (sets input data type to PROTEIN)
               if (NULL != optarg) {
                 if (testP (optarg))
                   input->model = getModel_PROTEIN (optarg);
                 else {
                   fprintf (stderr, "\n. Error: -P option: '%s' invalid evolutionary model.\n", optarg);
                   exit(EXIT_FAILURE);
                 }
               }
               else
                 input->model = LG;
               break;
     
             case 'r':
               // remove any site which has a gap in any sequence
               input->no_gap = TRUE;
               break;

             case 's':
               // SPR postprocessing
               input->use_SPR = TRUE;
               break;
     
             case 'S':
               strncpy (input->I_seqmat_file, optarg, MAX_FILE_NAME_LENGTH);
               input->model = SCOREDIST;
               input->input_type = PROTEIN;
               break;
     
             case 't':
               // TBR postprocessing
               input->use_TBR = TRUE;
               break;
     
             case 'T':
               // input topology file
               strncpy (input->I_tree_file, optarg, MAX_FILE_NAME_LENGTH);
               input->method = USER;
               break;
     
             case 'v':
               // verbose
               verbose = TRUE;
               break;

             case 'w':
               // branch lengths to assign to a topology
               if (testW (optarg, TRUE))
                 input->branch = getW (optarg);
               else {
                 fprintf (stderr, "\n. Error: -w option: '%s' invalid branch length type.\n", optarg);
                 exit(EXIT_FAILURE);
               }
               break;
     
             case 'z':
               // seed for randomization
               input->seed = atoi (optarg);
               break;
     
             case '?':
               // getopt_long already printed an error message.
               break;
     
             default:
               abort ();
             }
         }

     
       // Print any remaining command line arguments (not options).
/*
       if (optind < argc)
         {
           printf ("non-option ARGV-elements: ");
           while (optind < argc)
             printf ("%s ", argv[optind++]);
           putchar ('\n');
         }
*/

	return;
}

/*********************************************************/

void Get_Input_Interactive (Options *input)
{
	int n_trial = 0;
	char choix;
	char *tmp;
	
	printf ("Enter your input data file name > "); fflush (NULL);
	Getstring_Stdin (input->I_data_file);
	while (! Filexists (input->I_data_file))
	{
		tmp = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
		if (++n_trial > 10) {
			strncpy (tmp, "\nErr : the file '", 17);
			strncat (tmp, input->I_data_file, strlen (input->I_data_file));
			strncat (tmp, "' doesn't exist\n", 16);
			Exit (tmp);
		}
		strncpy (tmp, "The file '", 10);
		strncat (tmp, input->I_data_file, strlen (input->I_data_file));
		strncat (tmp, "' doesn't exist\n", 16);
		printf (tmp);
		printf ("Enter your input data file name > "); fflush (NULL);
		Getstring_Stdin (input->I_data_file);
		free (tmp);
	}

	strncpy (input->O_tree_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 17);
	strncat (input->O_tree_file, "_fastme_tree.txt", 16);
	strncpy (input->O_stat_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 17);
	strncat (input->O_stat_file, "_fastme_stat.txt", 16);
	strncpy (input->O_boot_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 17);
	strncat (input->O_boot_file, "_fastme_boot.txt", 16);

	choix = 0;
	do
	{
	#ifdef WIN32
		system ("cls");
	#elif UNIX
		printf ("\033[2J\033[H");
	#endif
	
		printf ("\n - FastME %s - \n\n\n", VERSION);
		printf ("Settings for this run:\n\n");
		
		tmp = (char *) mCalloc (MAX_NAME_LENGTH, sizeof(char));
		constantToStr (input->method, tmp);
		printf ("  M "
			"                                      Initial tree: build method\n"
			"              (balanced_GME, OLS_GME, BIONJ, NJ or UNJ) or user tree "
			" %-15s \n", tmp);
		free (tmp);
	
		printf ("  I "
			"         Input data type (distance matrix or sequence alignment) "
			" %-15s \n", (input->input_type == MATRIX ? "Distance matrix" : (input->input_type == DNA ? "DNA sequence alignment" : "Protein sequence alignment")));
	
		if (input->input_type != MATRIX) {
			tmp = (char *) mCalloc (MAX_NAME_LENGTH, sizeof(char));
			if (input->input_type == DNA) {
				if (input->model == NONE)
					input->model = F84;
				constantToStr (input->model, tmp);
				printf ("  E "
				"                        DNA evolutionary model or Scoring matrix\n"
				"                   (F84, TN93, K2P, JC69, Transversion Only, LogDet) "
				" %-15s \n", tmp);
			}
			else if (input->input_type == PROTEIN) {
				if (input->model == NONE)
					input->model = LG;
				constantToStr (input->model, tmp);
				printf ("  E "
				"                    Protein evolutionary model or Scoring matrix\n"
				"             (WAG, Dayhoff, JTT, MtREV, RtREV, CpREV, DCMut, VT, LG) "
				" %-15s \n", tmp);
			}
			free (tmp);
			
			printf ("  B "
				"                                 Bootstrap: number of replicates "
				" %-15d \n", input->nb_bootstraps);
	
			printf ("  F "
				"                                  Format of input sequences file "
				" %-15s \n", (input->is_interleaved ? "Interleaved" : "Sequential"));
				
			printf ("  G "
				"                               Gamma rate variation across sites "
				" %-15f \n", input->gamma);
			
			printf ("  R "
				"                                         Remove sites whith gaps "
				" %-15s \n", (input->no_gap ? "yes" : "no"));
				
			printf ("  O "
				"                               Output calculated distance matrix "
				" %-15s \n", (input->use_O_mat_file ? "yes" : "no"));
			
			printf ("\n");
		}
	
		printf ("  D "
			"                                              Number of datasets "
			" %-15d \n", input->nb_datasets);

		tmp = (char *) mCalloc (MAX_NAME_LENGTH, sizeof(char));
		constantToStr (input->NNI, tmp);
		printf ("  N "
			"           Type of tree swapping (balanced_NNI, OLS_NNI or none) "
			" %-15s \n", tmp);
		free (tmp);

		if (input->NNI == NONE)
		{
			tmp = (char *) mCalloc (MAX_NAME_LENGTH, sizeof(char));
			constantToStr (input->branch, tmp);
			if (input->method == BAL || input->method == OLS || input->method == USER)
				printf ("  W "
				"       Branch lengths assigned to the topology (balanced or OLS) "
				" %-15s \n", tmp);
			else
				printf ("  W "
				" Branch lengths assigned to the topology (balanced, OLS or none) "
				" %-15s \n", tmp);
			free (tmp);
		}
		
		printf ("  S "
			"                                              SPR postprocessing "
			" %-15s \n", (input->use_SPR ? "yes" : "no"));
	
		printf ("  T "
			"                                              TBR postprocessing "
			" %-15s \n", (input->use_TBR ? "yes" : "no"));
		

		printf ("\n");
		printf ("\nAre these settings correct? "
			"(type  Y  or letter for one to change)  ");

		scanf ("%c", &choix);
		if (choix == '\n') choix = 'X'; 
		else getchar(); // \n

		Uppercase (&choix);

		if ((choix == 'Y'))
			break;
	
		switch (choix)
		{
		
			case 'B':
			{
				char *c;
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
				askOption ("Number of replicates > ", c);
				n_trial = 0;
				while ((!atoi (c)) || (atoi (c) < 0))
				{
					if (++n_trial > 10)
						Exit ("\nErr : invalid number of replicates choosen\n");
					printf ("\nInvalid number of replicates choosen\n");
					askOption ("Enter a new value > ", c);
				}
				input->nb_bootstraps = atoi (c);
				free (c);
				break;
			}
			
			case 'D':
			{
				char *c;
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
				askOption ("Number of datatsets > ", c);
				n_trial = 0;
				while ((!atoi (c)) || (atoi (c) <= 0))
				{
					if (++n_trial > 10)
						Exit ("\nErr : the number of datasets must be a positive integer\n");
					printf ("\nThe number of datasets must be a positive integer\n");
					askOption ("Enter a new value > ", c);
				}
				input->nb_datasets = atoi (c);
				free (c);
				break;
			}
			
			case 'E':
			{
				char *c;
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
				if (input->input_type == DNA) {
					askOption ("Choose the DNA evolutionary model:\n  - (F)84, (T)N93, (K)2P, (J)C69, t(r)ansversion only, (L)ogDet\n  - (S)coring matrix\n> ", c);
					n_trial = 0;
					while ((! testD (c)) && (*c != 'S') && (*c != 's'))
					{
						if (++n_trial > 10)
							Exit ("\nErr : invalid evolutionary model choosen\n");
						printf ("\nInvalid evolutionary model choosen\n");
						askOption ("Enter a new value > ", c);
					}
					input->model = getModel_DNA(c);
				}
				else if (input->input_type == PROTEIN) {
					askOption ("Choose the PROTEIN evolutionary model:\n  - (C)pREV, (D)CMut, Day(h)off, (J)TT, (L)G, (M)tREV, (R)tREV, (V)T, (W)AG\n  - (S)coring matrix\n> ", c);
					n_trial = 0;
					while ((! testP (c)) && (*c != 'S') && (*c != 's'))
					{
						if (++n_trial > 10)
							Exit ("\nErr : invalid evolutionary model choosen\n");
						printf ("\nInvalid evolutionary model choosen\n");
						askOption ("Enter a new value > ", c);
					}
					input->model = getModel_PROTEIN(c);
				}
				else
					input->model = NONE;
				if (input->model == SCOREDIST) {
					printf ("Enter the scoring matrix file name > "); fflush (NULL);
					Getstring_Stdin (input->I_seqmat_file);
					n_trial = 0;
					while (! Filexists (input->I_seqmat_file))
					{
						tmp = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
						if (++n_trial > 10) {
							strncpy (tmp, "\nErr : the file '", 17);
							strncat (tmp, input->I_seqmat_file, strlen (input->I_seqmat_file));
							strncat (tmp, "' doesn't exist\n", 16);
							Exit (tmp);
						}
						strncpy (tmp, "The file '", 10);
						strncat (tmp, input->I_seqmat_file, strlen (input->I_seqmat_file));
						strncat (tmp, "' doesn't exist\n", 16);
						printf (tmp);
						printf ("Enter the scoring matrix file name > "); fflush (NULL);
						Getstring_Stdin (input->I_seqmat_file);
						free (tmp);
					}
				}
				free (c);
				break;
			}
			
			case 'F':
			{
				input->is_interleaved = abs (input->is_interleaved - 1);
				break;
			}
			
			case 'G':
			{
				char *c;
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
				askOption ("Gamma rate variation > ", c);
				n_trial = 0;
				while ((!atof (c)) || (atof (c) <= 0))
				{
					if (++n_trial > 10)
						Exit ("\nErr : invalid value for gamma rate variation\n");
					printf ("\nInvalid value for gamma rate variation\n");
					askOption ("Enter a new value > ", c);
				}
				input->gamma = atof (c);
				free (c);
				break;
			}
			
			case 'I':
			{
				char *c;
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
				askOption ("Choose your input data type: distance (M)atrix, (D)NA alignment, (P)rotein alignment > ", c);
				n_trial = 0;
				while (! testI (c))
				{
					if (++n_trial > 10)
						Exit ("\nErr : invalid datatype choosen\n");
					printf ("\nInvalid datatype choosen\n");
					askOption ("Enter a new value > ", c);
				}
                		input->input_type = getI (c);
                		if (input->input_type == DNA)
                			input->model = F84;
                		else if (input->input_type == PROTEIN)
                			input->model = LG;
                		else
                			input->model = NONE;
				free (c);
				break;
			}
			
			case 'M':
			{
                		boolean ask_I_tree_file = FALSE;
				char *c;
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
				askOption ("Choose your method: (b)alanced_GME, (O)LS_GME, B(I)ONJ, (N)J, (U)NJ or u(s)er > ", c);
				n_trial = 0;
				while (! testM (c))
				{
					if (++n_trial > 10)
						Exit ("\nErr : invalid method choosen\n");
					printf ("\nInvalid method choosen\n");
					askOption ("Enter a new value > ", c);
				}
                		input->method = getM (c);
				if (input->method == BAL || input->method == OLS)
					input->branch = input->method;
				else if (input->method == USER) {
					if (strlen (input->I_tree_file) > 0) {
						printf ("\nStarting tree topology file: '%s'\n", input->I_tree_file);
						askOption ("Do you wish to change ? (Y/n) > ", c);
						if ((*c == 'y') || (*c == 'Y') || (*c == '\n'))
							ask_I_tree_file = TRUE;
					}
					else
						ask_I_tree_file = TRUE;
					if (ask_I_tree_file) {
						printf ("\nEnter the starting tree topology file name > "); fflush (NULL);
						Getstring_Stdin (input->I_tree_file);
						n_trial = 0;
						while (! Filexists (input->I_tree_file))
						{
							tmp = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
							if (++n_trial > 10) {
								strncpy (tmp, "\nErr : the file '", 17);
								strncat (tmp, input->I_tree_file, strlen (input->I_tree_file));
								strncat (tmp, "' doesn't exist\n", 16);
								Exit (tmp);
							}
							strncpy (tmp, "The file '", 10);
							strncat (tmp, input->I_tree_file, strlen (input->I_tree_file));
							strncat (tmp, "' doesn't exist\n", 16);
							printf (tmp);
							printf ("Enter the starting tree topology file name > "); fflush (NULL);
							Getstring_Stdin (input->I_tree_file);
							free (tmp);
						}
					}
				}
				else
					input->branch = NONE;
				free (c);
				break;
			}

			case 'N':
			{
				char *c;
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
				askOption ("Choose your tree swapping: (b)alanced_NNI, (O)LS_NNI or (n)one > ", c);
				n_trial = 0;
				while (! testN (c))
				{
					if (++n_trial > 10)
						Exit ("\nErr : invalid NNI choosen\n");
					printf ("\nInvalid NNI choosen\n");
					askOption ("Enter a new value > ", c);
				}
                		input->NNI = getN (c);
				free (c);
				break;
			}
			
			case 'O':
			{
				input->use_O_mat_file = abs (input->use_O_mat_file - 1);
				break;
			}
			
			case 'R':
			{
				input->no_gap = abs (input->no_gap - 1);
				break;
			}
			
			case 'S':
			{
				input->use_SPR = abs (input->use_SPR - 1);
				break;
			}
			
			case 'T':
			{
				input->use_TBR = abs (input->use_TBR - 1);
				break;
			}

			case 'W':
			{
				boolean none;
				char *c;
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
				if (input->NNI == NONE)
				{
					if ((input->method == BAL) || (input->method == OLS) || (input->method == USER)) {
						askOption ("Choose your branch lengths to assign to the topology: (b)alanced or (O)LS) > ", c);
						none = FALSE;
					}
					else {
						askOption ("Choose your branch lengths to assign to the topology: (b)alanced, (O)LS) or (n)one > ", c);
						none = TRUE;
					}
					n_trial = 0;
					while (! testW (c, none))
					{
						if (++n_trial > 10)
							Exit ("\nErr : invalid branch length\n");
						printf ("\nInvalid branch length\n");
						askOption ("Enter a new value > ", c);
					}
					input->branch = getW (c);
				}
				free (c);
				break;
			}
			
			default:
				break;
		}

	} while (1);

	return;
}

