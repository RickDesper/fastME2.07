#include "utils.h"

int *initZeroArray(int l)
{
  int *x;
  int i;
  x = (int *)calloc(l,sizeof(int));
  for(i=0;i<l;i++)
    x[i] = 0;
  return(x);
}

int *initOneArray(int l)
{
  int *x;
  int i;
  x = (int *)calloc(l,sizeof(int));
  for(i=0;i<l;i++)
    x[i] = 1;
  return(x);
}

double **initDoubleMatrix(int d)
{
  int i,j;
  double **A;
  A = (double **) malloc(d*sizeof(double *));
  for(i=0;i<d;i++)
    {
      A[i] = (double *) malloc(d*sizeof(double));
      for(j=0;j<d;j++)
	A[i][j] = 0.0;
    }
  return(A);
}

boolean whiteSpace(char c)
{
  if ((' ' == c) || ('\t' == c) || ('\n' == c))
    return(TRUE);
  else
    return(FALSE);
}

void Exit (char *message)
{
	fprintf (stderr, "\n Error: %s\n", message);
	exit (EXIT_FAILURE);
	return;
}

void Warning (char *message)
{
	printf ("\n Warning: %s\n", message);
	return;
}

