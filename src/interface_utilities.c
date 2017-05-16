#include "interface_utilities.h"

/*********************************************************/

void Free_Input (Options *input)
{
	if (NULL != input->I_data_file)
		free (input->I_data_file);
	if (NULL != input->I_tree_file)
		free (input->I_tree_file);
	if (NULL != input->I_seqmat_file)
		free (input->I_seqmat_file);
	if (NULL != input->O_tree_file)
		free (input->O_tree_file);
	if (NULL != input->O_mat_file)
		free (input->O_mat_file);
	if (NULL != input->O_stat_file)
		free (input->O_stat_file);
	if (NULL != input->open_mode)
		free (input->open_mode);
	if (NULL != input->fpI_data_file)
		fclose (input->fpI_data_file);
	if (NULL != input->fpI_tree_file)
		fclose (input->fpI_tree_file);
	if (NULL != input->fpI_seqmat_file)
		fclose (input->fpI_seqmat_file);
	if (NULL != input->fpO_tree_file)
		fclose (input->fpO_tree_file);
	if (NULL != input->fpO_mat_file)
		fclose (input->fpO_mat_file);
	if (NULL != input->fpO_stat_file)
		fclose (input->fpO_stat_file);
	if (NULL != input->fpO_boot_file)
		fclose (input->fpO_boot_file);
	if (NULL != input)
		free (input);
	return;
}

/*********************************************************/
	
int Filexists (char *filename)
{ 
	FILE *fp;
	fp = fopen (filename, "r");
	if (fp) {
		fclose(fp);
		return 1;
	}
	else
		return 0;
}

/*********************************************************/

void Getstring_Stdin (char *file_name)
{ 
	fgets (file_name, MAX_FILE_NAME_LENGTH, stdin);
	if (strchr (file_name, '\n') != NULL)
		*strchr (file_name, '\n') = '\0';
	return;
}

/*********************************************************/

FILE *Openfile (char *filename, char *mode)
{
	FILE *fp;
	fp = NULL;
	char *tmp;
	tmp = (char *) mCalloc (MAX_NAME_LENGTH, sizeof(char));
	snprintf (tmp, MAX_NAME_LENGTH, "Cannot open file '%s'.", filename);

	fp = fopen(filename, mode);
	if (!fp)
		Exit (tmp);
	free (tmp);
	return fp;
}

/*********************************************************/

void Uppercase (char *ch)
{
	/* convert ch to upper case -- either ASCII or EBCDIC */
	int i;
	for (i = 0; i < strlen (ch); ++i) {
		ch[i] = isupper((int)ch[i]) ? ch[i] : toupper((int)ch[i]);
	}
	return;
}

/*********************************************************/

void *mCalloc (int nb, size_t size)
{
	void *allocated = NULL;
	if ((allocated = calloc ((size_t) nb, (size_t) size)) == NULL)
		Exit ("\n. Err: low memory\n");
	return allocated;
}

/*********************************************************/

void askOption (char *question, char *c)
{
	printf (question);
	Getstring_Stdin (c);
//	Uppercase (c);
	return;
}

/*********************************************************/

boolean testM (char *c)
{
	boolean ret = FALSE;
	Uppercase (c);
	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "N", 1) == 0 || strncmp (c, "B", 1) == 0 || strncmp (c, "O", 1) == 0
		 || strncmp (c, "I", 1) == 0 || strncmp (c, "U", 1) == 0 || strncmp (c, "S", 1) == 0)
			ret = TRUE;
		break;
	case 2:
		if (strncmp (c, "NJ", 2) == 0)
			ret = TRUE;
		break;
	case 3:
		if (strncmp (c, "OLS", 3) == 0 || strncmp (c, "UNJ", 3) == 0)
			ret = TRUE;
		break;
	case 4:
		if (strncmp (c, "USER", 4) == 0)
			ret = TRUE;
		break;
	case 5:
		if (strncmp (c, "BIONJ", 5) == 0)
			ret = TRUE;
		break;
	case 7:
		if (strncmp (c, "OLS_GME", 7) == 0)
			ret = TRUE;
		break;
	case 8:
		if (strncmp (c, "BALANCED", 8) == 0)
			ret = TRUE;
		break;
	case 12:
		if (strncmp (c, "BALANCED_GME", 12) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}
	return ret;
}

/*********************************************************/

int getM (char *c)
{
	int ret = BAL;
	Uppercase (c);
	if (strncmp (c, "B", 1) == 0 || strncmp (c, "BALANCED", 8) == 0 || strncmp (c, "BALANCED_GME", 12) == 0) {
		ret = BAL;
	}
	else if (strncmp (c, "O", 1) == 0 || strncmp (c, "OLS", 3) == 0 || strncmp (c, "OLS_GME", 7) == 0) {
		ret = OLS;
	}
	else if (strncmp (c, "N", 1) == 0 || strncmp (c, "NJ", 2) == 0) {
		ret = NJ;
	}
	else if (strncmp (c, "I", 1) == 0 || strncmp (c, "BIONJ", 5) == 0) {
		ret = BIONJ;
	}
	else if (strncmp (c, "U", 1) == 0 || strncmp (c, "UNJ", 3) == 0) {
		ret = UNJ;
	}
	else if (strncmp (c, "S", 1) == 0 || strncmp (c, "USER", 4) == 0) {
		ret = USER;
	}
	return ret;
}

/*********************************************************/

boolean testN (char *c)
{
	boolean ret = FALSE;
	Uppercase (c);
	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "B", 1) == 0 || strncmp (c, "O", 1) == 0 || strncmp (c, "N", 1) == 0)
			ret = TRUE;
		break;
	case 3:
		if (strncmp (c, "OLS", 3) == 0)
			ret = TRUE;
		break;
	case 4:
		if (strncmp (c, "NONE", 4) == 0)
			ret = TRUE;
		break;
	case 7:
		if (strncmp (c, "OLS_NNI", 7) == 0)
			ret = TRUE;
		break;
	case 8:
		if (strncmp (c, "BALANCED", 8) == 0)
			ret = TRUE;
		break;
	case 12:
		if (strncmp (c, "BALANCED_NNI", 12) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}
	return ret;
}

/*********************************************************/

int getN (char *c)
{
	int ret = BAL;
	Uppercase (c);
	if (strncmp (c, "B", 1) == 0 || strncmp (c, "BALANCED", 8) == 0 || strncmp (c, "BALANCED_NNI", 12) == 0) {
		ret = BAL;
	}
	else if (strncmp (c, "O", 1) == 0 || strncmp (c, "OLS", 3) == 0 || strncmp (c, "OLS_NNI", 7) == 0) {
		ret = OLS;
	}
	else if (strncmp (c, "N", 1) == 0 || strncmp (c, "NONE", 4) == 0) {
		ret = NONE;
	}
	return ret;
}

/*********************************************************/

boolean testW (char *c, boolean none)
{
	boolean ret = FALSE;
	Uppercase (c);
	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "B", 1) == 0 || strncmp (c, "O", 1) == 0)
			ret = TRUE;
		else if  (strncmp (c, "N", 1) == 0 && none)
			ret = TRUE;
		break;
	case 3:
		if (strncmp (c, "OLS", 3) == 0)
			ret = TRUE;
		break;
	case 4:
		if (strncmp (c, "NONE", 4) == 0 && none)
			ret = TRUE;
		break;
	case 8:
		if (strncmp (c, "BALANCED", 8) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}
	return ret;
}

/*********************************************************/

int getW (char *c)
{
	int ret = BAL;
	Uppercase (c);
	if (strncmp (c, "B", 1) == 0 || strncmp (c, "BALANCED", 8) == 0) {
		ret = BAL;
	}
	else if (strncmp (c, "O", 1) == 0 || strncmp (c, "OLS", 3) == 0) {
		ret = OLS;
	}
	else if (strncmp (c, "N", 1) == 0 || strncmp (c, "NONE", 4) == 0) {
		ret = NONE;
	}
	return ret;
}

/*********************************************************/

boolean testD (char *c)
{
//(F)84, (T)N93, (K)2P, (J)C69, T(R)ANSVERSIONONLY, (L)OGDET
	boolean ret = FALSE;
	Uppercase (c);
	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "F", 1) == 0 || strncmp (c, "T", 1) == 0 || strncmp (c, "K", 1) == 0
		 || strncmp (c, "J", 1) == 0 || strncmp (c, "R", 1) == 0 || strncmp (c, "L", 1) == 0)
			ret = TRUE;
		break;
	case 3:
		if (strncmp (c, "K2P", 3) == 0 || strncmp (c, "F84", 3) == 0)
			ret = TRUE;
		break;
	case 4:
		if (strncmp (c, "TN93", 4) == 0 || strncmp (c, "JC69", 4) == 0)
			ret = TRUE;
		break;
	case 6:
		if (strncmp (c, "LOGDET", 6) == 0)
			ret = TRUE;
		break;
	case 12:
		if (strncmp (c, "TRANSVERSION", 12) == 0)
			ret = TRUE;
		break;
	case 16:
		if (strncmp (c, "TRANSVERSIONONLY", 16) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}
	return ret;
}

/*********************************************************/

boolean testP (char *c)
{
//(C)pREV, (D)CMut, Day(h)off, (J)TT, (L)G, (M)tREV, (R)tREV, (V)T, (W)AG, (S)coredsit
	boolean ret = FALSE;
	Uppercase (c);
	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "C", 1) == 0 || strncmp (c, "D", 1) == 0 || strncmp (c, "H", 1) == 0
		 || strncmp (c, "J", 1) == 0 || strncmp (c, "L", 1) == 0 || strncmp (c, "M", 1) == 0
		 || strncmp (c, "R", 1) == 0 || strncmp (c, "V", 1) == 0 || strncmp (c, "W", 1) == 0 || strncmp (c, "S", 1) == 0)
			ret = TRUE;
		break;
	case 2:
		if (strncmp (c, "LG", 2) == 0 || strncmp (c, "VT", 2) == 0)
			ret = TRUE;
		break;
	case 3:	
		if (strncmp (c, "JTT", 3) == 0 || strncmp (c, "WAG", 3) == 0)
			ret = TRUE;
		break;
	case 5:
		if (strncmp (c, "CPREV", 5) == 0 || strncmp (c, "DCMUT", 5) == 0 || strncmp (c, "MTREV", 5) == 0
		 || strncmp (c, "RTREV", 5) == 0)
			ret = TRUE;
		break;
	case 7:
		if (strncmp (c, "DAYHOFF", 7) == 0)
			ret = TRUE;
		break;
	case 9:
		if (strncmp (c, "SCOREDIST", 9) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}
	return ret;
}
/*********************************************************/

boolean testI (char *c)
{
//(M)ATRIX, (D)NA, (P)ROTEIN
	boolean ret = FALSE;
	Uppercase (c);
	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "M", 1) == 0 || strncmp (c, "D", 1) == 0 || strncmp (c, "P", 1) == 0)
			ret = TRUE;
		break;
	case 3:
		if (strncmp (c, "DNA", 3) == 0)
			ret = TRUE;
		break;
	case 6:
		if (strncmp (c, "MATRIX", 6) == 0)
			ret = TRUE;
		break;
	case 7:
		if (strncmp (c, "PROTEIN", 7) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}
	return ret;
}

/*********************************************************/

int getI (char *c)
{
	int ret = MATRIX;
	Uppercase (c);
	if (strncmp (c, "M", 1) == 0 || strncmp (c, "MATRIX", 6) == 0) {
		ret = MATRIX;
	}
	else if (strncmp (c, "D", 1) == 0 || strncmp (c, "DNA", 3) == 0) {
		ret = DNA;
	}
	else if (strncmp (c, "P", 1) == 0 || strncmp (c, "PROTEIN", 7) == 0) {
		ret = PROTEIN;
	}
	return ret;
}

/*********************************************************/

boolean testF (char *c)
{
//(I)nterleaved, (S)equential
	boolean ret = FALSE;
	Uppercase (c);
	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "I", 1) == 0 || strncmp (c, "S", 1) == 0 )
			ret = TRUE;
		break;
	case 10:
		if (strncmp (c, "SEQUENTIAL", 10) == 0)
			ret = TRUE;
		break;
	case 11:
		if (strncmp (c, "INTERLEAVED", 11) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}
	return ret;
}

/*********************************************************/

boolean getF (char *c)
{
	int ret = TRUE;
	Uppercase (c);
	if (strncmp (c, "I", 1) == 0 || strncmp (c, "INTERLEAVED", 11) == 0) {
		ret = TRUE;
	}
	else if (strncmp (c, "S", 1) == 0 || strncmp (c, "SEQUENTIAL", 10) == 0) {
		ret = FALSE;
	}
	return ret;
}

/*********************************************************/

int getModel_DNA (char *c)
{
//(F)ELSENSTEIN, (T)AMURA, (K)IMURA, (J)C, T(R)ANSVERSIONONLY, (L)OGDET, (P)ROTEIN, or (S)COREDIST
	int ret = NONE;
	Uppercase (c);
	if (strncmp (c, "F", 1) == 0 || strncmp (c, "F84", 3) == 0) {
		ret = F84;
	}
	else if (strncmp (c, "T", 1) == 0 || strncmp (c, "TN93", 4) == 0) {
		ret = TN93;
	}
	else if (strncmp (c, "K", 1) == 0 || strncmp (c, "K2P", 3) == 0) {
		ret = K2P;
	}
	else if (strncmp (c, "J", 1) == 0 || strncmp (c, "JC69", 4) == 0) {
		ret = JC69;
	}
	else if (strncmp (c, "R", 1) == 0 || strncmp (c, "TRANSVERSION", 12) == 0 || strncmp (c, "TRANSVERSIONONLY", 16) == 0) {
		ret = TRANSVERSIONONLY;
	}
	else if (strncmp (c, "L", 1) == 0 || strncmp (c, "LOGDET", 6) == 0) {
		ret = LOGDET;
	}
/*	else if (strncmp (c, "P", 1) == 0 || strncmp (c, "PROTEIN", 7) == 0) {
		ret = PROTEIN;
	}
*/	else if (strncmp (c, "S", 1) == 0 || strncmp (c, "SCOREDIST", 9) == 0) {
		ret = SCOREDIST;
	}
	return ret;
}

/*********************************************************/

int getModel_PROTEIN (char *c)
{
//(C)pREV, (D)CMut, Day(h)off, (J)TT, (L)G, (M)tREV, (R)tREV, (V)T, (W)AG, or (S)COREDIST
	int ret = NONE;
	Uppercase (c);
	if (strncmp (c, "C", 1) == 0 || strncmp (c, "CPREV", 5) == 0) {
		ret = CPREV;
	}
	else if (strncmp (c, "D", 1) == 0 || strncmp (c, "DCMUT", 5) == 0) {
		ret = DCMUT;
	}
	else if (strncmp (c, "H", 1) == 0 || strncmp (c, "DAYHOFF", 7) == 0) {
		ret = DAYHOFF;
	}
	else if (strncmp (c, "J", 1) == 0 || strncmp (c, "JTT", 3) == 0) {
		ret = JTT;
	}
	else if (strncmp (c, "L", 1) == 0 || strncmp (c, "LG", 2) == 0) {
		ret = LG;
	}
	else if (strncmp (c, "M", 1) == 0 || strncmp (c, "MTREV", 5) == 0) {
		ret = MTREV;
	}
	else if (strncmp (c, "R", 1) == 0 || strncmp (c, "RTREV", 5) == 0) {
		ret = RTREV;
	}
	else if (strncmp (c, "V", 1) == 0 || strncmp (c, "VT", 2) == 0) {
		ret = VT;
	}
	else if (strncmp (c, "W", 1) == 0 || strncmp (c, "WAG", 2) == 0) {
		ret = WAG;
	}
	else if (strncmp (c, "S", 1) == 0 || strncmp (c, "SCOREDIST", 9) == 0) {
		ret = SCOREDIST;
	}
	return ret;
}

/*********************************************************/

void constantToStr (int c, char *str)
{
	switch (c)
	{
		case BAL:
			strncpy (str, "Balanced", MAX_NAME_LENGTH);
			break;
		case OLS:
			strncpy (str, "OLS", MAX_NAME_LENGTH);
			break;
		case NJ:
			strncpy (str, "NJ", MAX_NAME_LENGTH);
			break;
		case UNJ:
			strncpy (str, "UNJ", MAX_NAME_LENGTH);
			break;
		case BIONJ:
			strncpy (str, "BIONJ", MAX_NAME_LENGTH);
			break;
		case NONE:
			strncpy (str, "None", MAX_NAME_LENGTH);
			break;
		case USER:
			strncpy (str, "User", MAX_NAME_LENGTH);
			break;
		case F84:
			strncpy (str, "F84", MAX_NAME_LENGTH);
			break;
		case TN93:
			strncpy (str, "TN93", MAX_NAME_LENGTH);
			break;
		case K2P:
			strncpy (str, "K2P", MAX_NAME_LENGTH);
			break;
		case JC69:
			strncpy (str, "JC69", MAX_NAME_LENGTH);
			break;
		case TRANSVERSIONONLY:
			strncpy (str, "TransversionOnly", MAX_NAME_LENGTH);
			break;
		case LOGDET:
			strncpy (str, "LogDet", MAX_NAME_LENGTH);
			break;
		case PROTEIN:
			strncpy (str, "Protein", MAX_NAME_LENGTH);
			break;
		case MATRIX:
			strncpy (str, "Matrix", MAX_NAME_LENGTH);
			break;
		case DNA:
			strncpy (str, "DNA", MAX_NAME_LENGTH);
			break;
		case SCOREDIST:
			strncpy (str, "Scoring matrix", MAX_NAME_LENGTH);
			break;
		case WAG:
			strncpy (str, "WAG", MAX_NAME_LENGTH);
			break;
		case DAYHOFF:
			strncpy (str, "Dayhoff", MAX_NAME_LENGTH);
			break;
		case JTT:
			strncpy (str, "JTT", MAX_NAME_LENGTH);
			break;
		case BLOSUM62:
			strncpy (str, "BLOSUM62", MAX_NAME_LENGTH);
			break;
		case MTREV:
			strncpy (str, "MtREV", MAX_NAME_LENGTH);
			break;
		case RTREV:
			strncpy (str, "RtREV", MAX_NAME_LENGTH);
			break;
		case CPREV:
			strncpy (str, "CpREV", MAX_NAME_LENGTH);
			break;
		case DCMUT:
			strncpy (str, "DCMut", MAX_NAME_LENGTH);
			break;
		case VT:
			strncpy (str, "VT", MAX_NAME_LENGTH);
			break;
		case LG:
			strncpy (str, "LG", MAX_NAME_LENGTH);
			break;
		default:
			break;
	}
	return;
}

