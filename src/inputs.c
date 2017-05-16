#include "inputs.h"

void compareSets(tree *T, set *S)
{
  edge *e;
  node *v,*w;
  set *X;
  e = depthFirstTraverse(T,NULL);
  while (NULL != e)
    {
      v = e->head;
      for(X = S; NULL != X; X = X->secondNode)
	{
	  w = X->firstNode;
	  if (0 == strcmp(v->label,w->label))
	    {
	      v->index2 = w->index2;
	      w->index2 = -1;
	      break;
	    }
	}
      e = depthFirstTraverse(T,e);
    }
  v = T->root;
  for(X = S; NULL != X; X = X->secondNode)
    {
      w = X->firstNode;
      if (0 == strcmp(v->label,w->label))
	{
	  v->index2 = w->index2;
	  w->index2 = -1;
	  break;
	}
    }
  if (-1 == v->index2)
    {
      fprintf(stderr,"\n Error: leaf '%s' in tree not in distance matrix.\n",v->label);
      exit(EXIT_FAILURE);
    }
  e = depthFirstTraverse(T,NULL);
  while (NULL != e)
    {
      v = e->head;
      if ((leaf(v)) && (-1 == v->index2))
	{
	  fprintf(stderr,"\n Error: leaf '%s' in tree not in distance matrix.\n",v->label);
	  exit(EXIT_FAILURE);
	}
      e = depthFirstTraverse(T,e);
      }
  for(X = S; NULL != X; X = X->secondNode)
    if (X->firstNode->index2 > -1)
      {
	fprintf(stderr,"\n Error: node '%s' in matrix but not a leaf in tree.\n",X->firstNode->label);
	exit(EXIT_FAILURE);
      }
  return;
}

void freeCharMatrix(char **D, int size)
{
  int i;
  for(i=0;i<size;i++) {
    if (NULL != D[i])
      free(D[i]);
  }
  if (NULL != D)
    free(D);
  return;
}

void freeMatrix(double **D, int size)
{
  int i;
  for(i=0;i<size;i++) {
    if (NULL != D[i])
      free(D[i]);
  }
  if (NULL != D)
    free(D);
  return;
}

double **loadMatrix(FILE *ifile, int *size, set *S)
{
  char nextString[MAX_EVENT_NAME];
  node *v;
  double **table;
  double val;
  int i,j, zsize;
  if (!(fscanf(ifile,"%s",nextString)))
    Exit ("Cannot load input matrix.");
  zsize = atoi(nextString);
  if ((zsize < 0) || (zsize > MAXSIZE))
    Exit ("Invalid input size.");
  table = (double **) calloc(zsize,sizeof(double *));

  // zsize is used to avoid stack smashing which may occur
  // with a buffer overrun on the *size allocated memory
  // in case of wrong matrix format
  for(i=0;i<zsize;i++)
    {
      j = 0;
      table[i] = (double *) calloc(zsize,sizeof(double));
      if (!(fscanf(ifile,"%s",nextString)))
	{
	  fprintf(stderr,"\n Error loading label %d.\n",i);
	  exit(EXIT_FAILURE);
	}
      v = makeNode(nextString,-1);
      v->index2 = i;
      S = addToSet(v,S);

      while (j < zsize)
	{
	  if (!(fscanf(ifile,"%s",nextString)))
	    {
	      fprintf(stderr,"\n Error loading (%d,%d)-entry.\n",i,j);
	      exit(EXIT_FAILURE);
	    }
	  if ((nextString[0] < '0' || nextString[0] > '9') && nextString[0] != '.')
	  //if (nextString[0] < '0' || nextString[0] > '9')
	    {
	      fprintf(stderr,"\n Error: invalid distance matrix.");
	      fprintf(stderr,"\n Numerical value expected for taxon '%s' instead of '%s'.\n", v->label, nextString);
	      exit(EXIT_FAILURE);
	    }
	  val = atof(nextString);
	  if ((i != j) && (0 > val))
	    {
	      fprintf(stderr,"\n Error: program is expecting distance matrix.");
	      fprintf(stderr,"\n Input of %s off diagonal is inappropriate.\n",nextString);
	      exit(EXIT_FAILURE);
	    }
	  table[i][j++] = val;	  
	}	  
    }
  *size = zsize;
  return(table);
}

void partitionSizes(tree *T)
{
  edge *e;
  e = depthFirstTraverse(T,NULL);
  while (NULL != e)
    {
      if (leaf(e->head))
	e->bottomsize = 1;
      else
	e->bottomsize = e->head->leftEdge->bottomsize 
	  + e->head->rightEdge->bottomsize;
      e->topsize = (T->size + 2)/2 - e->bottomsize;
      e = depthFirstTraverse(T,e);
    }
  return;
}

char **loadSequences(FILE *ifile, int *numSeqs, set *taxa, int *inputSize)
{
  int i,j;
  node *v;
  char *nextString;
  char **data;
  nextString = (char *)malloc(MAX_LABEL_LENGTH*sizeof(char));  
  if (!(fscanf(ifile,"%s",nextString)))
    Exit ("Cannot load input file.");
  *numSeqs = atoi(nextString);
  fscanf(ifile,"%s",nextString);
  *inputSize = atoi(nextString);
  data = (char **)calloc(*numSeqs,sizeof(char *));
  
  if (*inputSize < 0)
    Exit ("Invalid sequence length.");
  if (verbose) {
    printf("  . numSeqs is %d.\n",*numSeqs);
    printf("  . inputSize is %d.\n",*inputSize);
  }
  for(i=0;i<*numSeqs;i++)
    {
      if (!(fscanf(ifile,"%s",nextString)))
	{
	  fprintf(stderr,"\n Error reading label of sequence %d.\n",i+1);
	  exit(EXIT_FAILURE);
	}
      v = makeNode(nextString,-1);
      v->index2 = i;
      taxa = addToSet(v,taxa);
      data[i] = (char *)calloc(*inputSize+MAX_LABEL_LENGTH,sizeof(char));
      for(j=0;j<*inputSize;j++)
	{
	  data[i][j] = getc(ifile);
	  if (('\t' == data[i][j]) || (' ' == data[i][j]))
	    j--;
	  else if (('\n' == data[i][j]) || ('\r' == data[i][j]))
	    {
	      data[i][j] = '\0';
	      break;
	    }
	}		     
    }
  free(nextString);
  return(data);
}

double **loadScoreMatrix(int *d, FILE *ifile, char *alphabet)
{
  char nextString[MAX_EVENT_NAME];
  double **table;
  int i,j;
  if (!(fscanf(ifile,"%s",nextString)))
    Exit ("Cannot load size of score matrix.");
  *d = atoi(nextString);
  if ((*d < 1) || (*d > PROTEIN_ALPHABET_SIZE))
    Exit ("Invalid size of score matrix.");
  table = initDoubleMatrix(*d);
  /* *d is actual size of input matrix, which may be less than PROTEIN_ALPHABET_SIZE, 
     since the ALPHABET contains special characters*/
  for(i=0;i<*d;i++)
    {
      j = 0;
      if (!(fscanf(ifile,"%s",nextString)))
	{
	  fprintf(stderr,"\n Error loading label %d.\n",i);
	  exit(EXIT_FAILURE);
	}
      alphabet[i] = nextString[0];
      while (j < *d)
	{
	  if (!(fscanf(ifile,"%s",nextString)))
	    {
	      fprintf(stderr,"\n Error loading (%d,%d)-entry.\n",i,j);
	      exit(EXIT_FAILURE);
	    }
	  table[i][j++] = atof(nextString);
	} /*j loop*/	  
    }  /*i loop*/
  if (*d < PROTEIN_ALPHABET_SIZE)
    alphabet[*d]='\0';
  return(table);
}

