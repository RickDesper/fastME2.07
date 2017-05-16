#include "distance.h"

/*edited 7/7/2007 to change default gamma from gamma=1.0 to gamma=0.0*/

int countStateChanges(char *s, char *t, int length, char c1, char c2, int *filter) /*returns the number of c1 -> c2 changes in s->t*/
{
  int i;
  int matches = 0;
  for(i=0;i<length;i++)
    if((c1 == s[i]) && (c2 == t[i]))
      matches+=filter[i];
  return(matches);
}

void ijFilter(int *filter, char *s1, char *s2, int itype, int seqlength)
{
	int i;
	int activeSiteCount;
	activeSiteCount = seqlength;
//	if ((PROTEIN == itype) || (SCOREDIST == itype))
	if (SCOREDIST == itype)
	  {
	    for(i=0;i<seqlength;i++)
	      {
		if ((NULL == strchr(PROTEIN_ALPHABET,s1[i])) || (NULL == strchr(PROTEIN_ALPHABET,s2[i])))
		  {
		    filter[i] = 0;
		    if (verbose) 
		      printf("  . Removing site %d.\n",i);
		    activeSiteCount--;
		    break;
		  }
		else if (('*' == s1[i]) || ('?' == s1[i]) || ('-' == s1[i]) ||
			 ('*' == s2[i]) || ('?' == s2[i]) || ('-' == s2[i]))
		  {
		    filter[i] = 0;
		    if (verbose)
		      printf("  . Removing site %d.\n",i);
		    activeSiteCount--;
		    break;
		  }
		else
		  filter[i] = 1;
	      }
	  }
	else
	  {
	    for(i=0;i<seqlength;i++)
	      {
		if ((NULL == strchr(DNA_ALPHABET,s1[i])) || (NULL == strchr(DNA_ALPHABET,s2[i])))
		  {
		    filter[i] = 0;
		    if (verbose)
		      printf("  . Removing site %d.\n",i);
		    activeSiteCount--;
		    break;
		  }
		else if (('*' == s1[i]) || ('?' == s1[i]) || ('-' == s1[i]) ||
			 ('*' == s2[i]) || ('?' == s2[i]) || ('-' == s2[i]))
		  {
		    filter[i] = 0;
		    if (verbose)
		      printf("  . Removing site %d.\n",i);
		    activeSiteCount--;
		    break;
		  }
	      }
	  }
	if (verbose)
	  printf("  . %d sites are active.\n",activeSiteCount);
  return;
}

void calcDNATransitions(double **P,char *s1, char *s2, int length, int *filter, int numSelected)
{
  int i,j;
  for(i=0;i<DNA_ALPHABET_SIZE;i++)
    {
      for(j=0;j<DNA_ALPHABET_SIZE;j++)
	P[i][j] =  countStateChanges(s1,s2,length,DNA_ALPHABET[i],
				     DNA_ALPHABET[j],filter);
    } /*note: P is NOT a probability matrix for this program
	sum_j P[i][j] = Pi[i]*/
  return;
}

int seqCharMatches(char *s, int length, char c, int *filter)
{
  int i;
  int matches = 0;
  for(i=0;i<length;i++)		
    if (c == s[i])
      matches+=filter[i];
  return(matches);
}

/*called when calculating stationary probabilities*/
int matrixCharMatches(char **s, int numSeqs, int length, char c, int *filter)
{
  int matches = 0;
  int i;
  for(i=0;i<numSeqs;i++)
    matches += seqCharMatches(s[i],length,c,filter);
  return(matches);
}

double *calcDNAStationaryProbs(char **s, int numSeqs, int length, int *filter, int numSelected)
{
  double *p;	
  int i;
  p = (double *)malloc(DNA_ALPHABET_SIZE*sizeof(double));
  for(i=0;i<DNA_ALPHABET_SIZE;i++)
    p[i] = (double) (matrixCharMatches(s,numSeqs,length,DNA_ALPHABET[i],filter)) / (numSeqs*numSelected);
  return(p);
}

double *calcProteinStationaryProbs(char **s, int numSeqs, int length, int *filter, int numSelected)
{
  double *p;	
  int i;
  p = (double *)malloc(PROTEIN_ALPHABET_SIZE*sizeof(double));
  for(i=0;i<PROTEIN_ALPHABET_SIZE;i++)
    p[i] = (double) (matrixCharMatches(s,numSeqs,length,PROTEIN_ALPHABET[i],filter)) /(numSeqs* numSelected);
  return(p);
}

void calcDNATransitionProbs(double **P,char *s1, char *s2, int length, int *filter, int numSelected)
{
  int i,j;
  for(i=0;i<DNA_ALPHABET_SIZE;i++)
    {
      for(j=0;j<DNA_ALPHABET_SIZE;j++)
	P[i][j] = (double) countStateChanges(s1,s2,length,DNA_ALPHABET[i],
				   DNA_ALPHABET[j],filter) / numSelected;
    } /*note: P is NOT a probability matrix for this program
	sum_j P[i][j] = Pi[i]*/
  return;
}

void calcProteinTransitionProbs(double **P,char *s1, char *s2, int length, int *filter, int numSelected)
{
  int i,j;
  for(i=0;i<PROTEIN_ALPHABET_SIZE;i++)
    {
      for(j=0;j<PROTEIN_ALPHABET_SIZE;j++)
	P[i][j] = (double) countStateChanges(s1,s2,length,PROTEIN_ALPHABET[i],
			   PROTEIN_ALPHABET[j],filter) / numSelected;
    } /*note: P is NOT a probability matrix for this program
	sum_j P[i][j] = Pi[i]*/
  return;
}

/*should check formulas for following two functions*/
/*transition rate is probability of witnessing a transition change from sequence i to 
  sequence j.  This rate should be equal to the sum of the four probabilities below*/
double calcTransitionRate(double **P)
{	
  double a;
  a = P[ADENINE][GUANINE] + P[GUANINE][ADENINE]  + P[CYTOSINE][THYMINE] 
    + P[THYMINE][CYTOSINE];
  return(a);
}

/*the transversion rate is the probability of seeing a non-transition change from sequence i
  to sequence j.  Essentially is the probability of a change minus the transitition probability*/
double calcTransversionRate(double **P)
{
  double b=0.0;
  b += P[ADENINE][CYTOSINE];
  b += P[ADENINE][THYMINE]; 
  b += P[GUANINE][CYTOSINE];
  b += P[GUANINE][THYMINE];  
  b += P[CYTOSINE][ADENINE];
  b += P[CYTOSINE][GUANINE]; 
  b += P[THYMINE][GUANINE];
  b += P[THYMINE][ADENINE];	
  return(b);
}

double calcK2P(double a, double b, double gamma) /*Kimura 2-parameter distance*/
{
  double returnValue, loc1, loc2;
  if ((0 == a) && (0 == b))
    return(0);
  loc1 = 1.0 - 2.0*a - b;
  loc2 = 1.0 - 2.0*b;
  if ((0 >= loc1) || (0 >= loc2))
    {
      Warning ("Distance goes to infinity, substituting large distance.");
      return(LARGE);
    }
  if (0.0 == gamma)
    returnValue = -0.5*(log(loc1)) - 0.25*(log(loc2));
  else
    returnValue = gamma*(0.5*(pow(loc1,-1.0/gamma))
			 +0.25*(pow(loc2,-1.0/gamma))
			 -0.75);
  return(returnValue);
}

void calcF84AuxProbs(double *Pi, double *Mloc, double *Nloc, double *Ploc)
{
  double PAG, PCT, PR, PY;
  
  PAG = Pi[ADENINE]*Pi[GUANINE]; 
  PCT = Pi[CYTOSINE]*Pi[THYMINE];
  PR = Pi[ADENINE] + Pi[GUANINE];
  PY = Pi[CYTOSINE] + Pi[THYMINE];
  
  *Mloc = PAG/PR + PCT/PY;
  *Nloc = PAG + PCT;
  *Ploc = PR*PY;
  return;
}

double calcF84(double a, double b, double gamma, double Mloc, double Nloc, 
	   double Ploc) 
     /*Pi stationary frequencies P transition probabilities*/
{
  double returnValue;
  double loc1;
  double loc2;
  if ((0 == a) && (0 == b))
    return(0);
  loc1 = 1.0 - a/(2.0*Mloc) - b*(Mloc - Nloc)/(2.0*Mloc*Ploc);
  loc2 = 1.0 - b/(2.0*Ploc);
  if ((0 >= loc1) || (0 >= loc2))
    {
      Warning ("Distance goes to infinity, substituting large distance.");
      return(LARGE);
    }   
  if (0.0 == gamma)
    returnValue = -2.0*Mloc*log(loc1) - 2.0*(Nloc + Ploc - Mloc)*log(loc2);
  else
    returnValue = 2*gamma*(Mloc*pow(loc1,-1.0/gamma)
			   + (Nloc + Ploc - Mloc)*pow(loc2,-1.0/gamma)
			   - Nloc - Ploc);
  return(returnValue);	
}

void calcTNAuxProbs(double *Pi,double *Mloc, double *Nloc, 
		    double *PR, double *PY)
{
  *PR = Pi[ADENINE] + Pi[GUANINE];
  *PY = Pi[CYTOSINE] + Pi[THYMINE];
  *Mloc = Pi[ADENINE]*Pi[GUANINE];
  *Nloc = Pi[CYTOSINE]*Pi[THYMINE];
  return;
}	

double calcTN93(double aR, double aY, double b, double PR, double PY, 
	  double PAPG, double PCPT, double gamma)
{
  double loc1, loc2;
  double logme;
  double returnValue;
  if ((0 == aR) && (0 == aY) && (0 == b))
    return(0);
  loc1 = PAPG/PR;
  loc2 = PCPT/PY;
  if (0.0 == gamma)
    {
      logme = (1-aR/(2*loc1) - b/(2*PR));
      if (logme <= 0)
	{
          Warning ("Distance goes to infinity, substituting large distance.");
	  return(LARGE);
	}
      returnValue = -2.0*loc1*log(logme);
      logme = (1-aY/(2*loc2) - b/(2*PY));
      if (logme <= 0)
	{
	  if (verbose)
	    Warning ("Logarithm of negative number.");
	  return(LARGE);
	}
      returnValue += -2.0*loc1*log(logme);
      returnValue -= 2.0*(PR*PY - PAPG*PY/PR - PCPT*PR/PY)*log(1 - b/(2*PR*PY));
    }
  else
    {
      returnValue = (2*gamma*loc1)*pow((1 - aR/(2*loc1) - b/(2*PR)),-1.0/gamma);
      returnValue += (2*gamma*loc2)*pow((1 - aY/(2*loc2) - b/(2*PY)),-1.0/gamma);	
      returnValue += 2*gamma*(PR*PY - loc1*PY - loc2*PR)*pow(1 - b/(PR*PY*2.0),-1.0/gamma);
      returnValue += -2*gamma*(PAPG + PCPT + PR*PY);
    }
  return(returnValue);
}

/*transversion only distance*/
double XX(double PR, double PY, double b, double gamma)
{
  double returnValue;
  double Z;
  Z = 1 - PR*PR - PY*PY;
  if (b/Z >= 1)
    {
      Warning ("Distance goes to infinity, substituting large distance.");
      return(LARGE);
    }
  if (0.0 == gamma)
    returnValue = -Z*log(1- b/Z);
  else
    returnValue = gamma*Z*(pow(1-b/Z,-1/gamma)-1);
  return(returnValue);
}

int factorial(int n) /*called only with positive integers*/
{
  if (1 == n)
    return(n);
  else
    return(n*factorial(n-1));
}

/*index refers to index of the permutation in the (I think) Gray code
  length is number of permtuations considered*/
/*index runs from 1 to factorial(length)*/
int *nextPerm(int *p, int index, int size, int length)
{
  int temp;
  if (0 ==  index % factorial(size))
    {
      temp = p[length-1];
      p[length-1] = p[length-1-size];
      p[length-1-size] = temp;
      return(p);
    }
  else
    return(nextPerm(p,index,size-1,length));
}

double permDiagProduct(double **P, int *p, int d)
{
  int i;
  double prod = 1.0;
  for(i=0;i<d;i++)
    prod = prod*P[i][p[i]];
  return(prod);
}

double det(double **P, int d)
{
  int *p;
  int signum = 1;
  int i;
  int numPerms;
  double det = 0;
  p = initPerm(d);
  numPerms = factorial(d);
  for(i=0;i<numPerms;i++)
    {
      det += signum*permDiagProduct(P,p,d);
      p = nextPerm(p,i+1,d-1,d);
      signum = -1*signum;
    }
  free(p);
  return(det);
}

double logdet(double **P, double *Pi1, double *Pi2)
{
  int i;
  double returnValue, detP;
  detP = det(P,4);
  if (0 >= detP)
    {
      Warning ("Distance goes to infinity, substituting large distance.");
      return(LARGE);
    }
  returnValue = -0.5* log(detP);
  for(i=0;i<4;i++)
    {
      if ((0 >= Pi1[i]) || (0 >= Pi2[i]))
	{
	  fprintf(stderr,"\n Error in logdet: value of Pi1[i] is %lf, of Pi2[i] is %lf, i is %d.\n",Pi1[i],Pi2[i],i);
	  exit(EXIT_FAILURE);
	}
      returnValue += (log(Pi1[i]) + log(Pi2[i]))/8;
    }
  return(returnValue);
}

int support(int *v, int length)
{
  int i;
  int count=0;
  for(i=0;i<length;i++)
    if (v[i])
      count++;
  return(count);
}

double HammingDistance(char *v1, char *v2, int *filter, int length, 
		       int numSelected)
{
  int d=0;
  int i;
  for(i=0;i<length;i++)
    if(v1[i] != v2[i])
      d+=filter[i];
  return((double) d/numSelected);
}

double protDiff(double *P)
{
  double sum = 0.0;
  int i;
  for(i=0;i<PROTEIN_ALPHABET_SIZE;i++)
    sum += P[i]*P[i];
  return(1-sum);
}

double protFormula(double b, double gamma, double pdiff)
{
  double d, y;
  y = 1 - b/pdiff;  
  if (y <= 0)
    return(LARGE);
  if (0.0 == gamma)
    d=-1.0*pdiff*log(y);
  else
    d = gamma*pdiff*(-1 + pow(y,-1/gamma));
  return(d);
}

int aaIndex(char s, char *alphabet, int d)
{
  int i;
  for(i=0;i<d;i++)
    if(s==alphabet[i])
      return(i);
  fprintf(stderr,"\n Error looking for character %c in protein string\n",s);
  exit(EXIT_FAILURE);
}

double simScore(char *s1, char *s2, double **scoreMatrix, int seqlength, char *alphabet, int alphabetSize)
{
  int i;
  double sum = 0.0;
  for(i=0;i<seqlength;i++)
    sum+=scoreMatrix[aaIndex(s1[i],alphabet,alphabetSize)][aaIndex(s2[i],alphabet,alphabetSize)];
  return(sum);
}

double expectedProtSimScore(double *P, double **scoreMatrix, int alphabetSize)
{
  int i,j;
  double sum=0;
  for(i=0;i<alphabetSize;i++)
    for(j=0;j<alphabetSize;j++)
      sum+=P[i]*P[j]*scoreMatrix[i][j];
  return(sum);
}

double scoreDistij(int i,int j,char *si, char *sj, int seqLength, double simExp,
	double *simDiags, double **scoreMatrix, char *alphabet, int alphabetSize)
{
  double simij, simUpper,ratio;
  simij = simScore(si,sj,scoreMatrix,seqLength,alphabet,alphabetSize);
  simUpper = 0.5*(simDiags[i]+simDiags[j]);
  ratio = (simij - simExp)/(simUpper - simExp);
  if (ratio <= 0 )
    return(LARGE);
  else
    return(-100.0*log(ratio));
}

void scoreDist(double *P, char **data, int numSpecies, int seqLength,
	double **scoreMatrix, double **D, char *alphabet, int alphabetSize)
{
  int i,j;
  double simExp;
  double *simDiags;
  simDiags = (double *)malloc(numSpecies*sizeof(double));
  simExp = seqLength*expectedProtSimScore(P,scoreMatrix,alphabetSize);
  for(i=0;i<numSpecies;i++)
    {
      D[i] = (double *)malloc(numSpecies*sizeof(double));
      D[i][i] = 0.0;
      simDiags[i] = simScore(data[i],data[i],scoreMatrix,seqLength,alphabet,alphabetSize);
    }
  for(i=0;i<numSpecies-1;i++)
    {
      for(j=i+1;j<numSpecies;j++)
	D[i][j] = D[j][i] = scoreDistij(i,j,data[i],data[j],seqLength,simExp,simDiags,scoreMatrix,alphabet,alphabetSize);
    }
  free(simDiags);
  return;
}

void gapCheckFilter(int *filter, int itype, int seqlength, int numSeqs, char **data)
{
  int i,j;
  int activeSiteCount;
  activeSiteCount = seqlength;
//  if ((PROTEIN == itype) || (SCOREDIST == itype))
  if (SCOREDIST == itype)
    {
      for(i=0;i<seqlength;i++)
	for(j=0;j<numSeqs;j++)
	  {
	    if (NULL == strchr(PROTEIN_ALPHABET,data[j][i]))
	      {
		filter[i] = 0;
		if (verbose) 
		  printf("  . Removing site %d.\n",i);
		activeSiteCount--;
		break;
	      }
	    else if (('*' == data[j][i]) || ('?' == data[j][i]) || ('-' == data[j][i]))
	      {
		filter[i] = 0;
		if (verbose)
		  printf("  . Removing site %d.\n",i);
		activeSiteCount--;
		break;
	      }
	  }
    }
  else
    {
      for(i=0;i<seqlength;i++)
	for(j=0;j<numSeqs;j++)
	  {
	    if (NULL == strchr(DNA_ALPHABET,data[j][i]))
	      {
		filter[i] = 0;
		if (verbose)
		  printf("  . Removing site %d.\n",i);
		activeSiteCount--;
		break;
	      }
	    else if (('*' == data[j][i]) || ('?' == data[j][i]) || ('-' == data[j][i]))
	      {
		filter[i] = 0;
		if (verbose)
		  printf("  . Removing site %d.\n",i);
		activeSiteCount--;
		break;
	      }
	  }
    }
  if (verbose)
    printf("  . %d sites are active.\n",activeSiteCount);
  return;
}

double **makeDistMatrix(char **data, int numSeqs, int numSites, double gamma,
	int model, int itype, int *filter, double **scoreMatrix, int scoreMatrixSize,
	char *alphabet, boolean gapCheck)
{
  double *PStationary; 
  double **PStateChanges, **Pi2;
  double a, b, a2;
  double Mloc, Nloc, Ploc;
  double PR, PY, pdiff;
  int i,j;
  int numSelected;
  double maxentry, minentry;
  double **D;
  PStationary = NULL;
  PStateChanges = Pi2 = NULL;
  Mloc = Nloc = Ploc = PR = PY = pdiff = 0.0;
  
  if (model == K2P || model == F84 || model == TN93 || model == JC69 ) 
    return  FastDist(data, numSeqs, numSites, model, gamma, filter);
      //double** FastDist(char** data, long n, long m, int model, double gamma)
  
  D = initDoubleMatrix(numSeqs);
  maxentry = 0.0;
  minentry = LARGE;
  if (gapCheck)
  {
	gapCheckFilter(filter,itype,numSites,numSeqs,data);
	numSelected = support(filter,numSites);
  }
  else
    numSelected = numSites;
  
  switch(model)
    /*calculations for various constants done outside i,j loop*/
    {
    case K2P:
      PStateChanges = initDoubleMatrix(DNA_ALPHABET_SIZE);
      break;
    case F84:
      PStateChanges  = initDoubleMatrix(DNA_ALPHABET_SIZE);
      PStationary  = calcDNAStationaryProbs(data,numSeqs,numSites,filter,numSelected);
      calcF84AuxProbs(PStationary,&Mloc,&Nloc,&Ploc);
      break;
    case TN93:
      PStateChanges  = initDoubleMatrix(DNA_ALPHABET_SIZE);
      PStationary  = calcDNAStationaryProbs(data,numSeqs,numSites,filter,numSelected);
      calcTNAuxProbs(PStationary,&Mloc,&Nloc,&PR,&PY);	  
      break;
    case TRANSVERSIONONLY:
      PStateChanges  = initDoubleMatrix(DNA_ALPHABET_SIZE);
      PStationary  = calcDNAStationaryProbs(data,numSeqs,numSites,filter,numSelected);
      PR = PStationary[ADENINE] + PStationary[GUANINE];
      PY = PStationary[CYTOSINE] + PStationary[THYMINE];
      break;
    case LOGDET:
      PStateChanges = initDoubleMatrix(DNA_ALPHABET_SIZE);
      PStationary  = calcDNAStationaryProbs(data,numSeqs,numSites,filter,numSelected);
      Pi2 = (double **) calloc(numSeqs,sizeof(double *));
      for(i=0;i<numSeqs;i++)
	Pi2[i] = calcDNAStationaryProbs(data+i,1,numSites,filter,numSelected);
      break;
/*
    case PROTEIN:
      PStateChanges = initDoubleMatrix(PROTEIN_ALPHABET_SIZE);    
      PStationary = calcProteinStationaryProbs(data,numSeqs,numSites,filter,numSelected);
      pdiff = protDiff(PStationary);
      break;
*/
    case SCOREDIST: /*SCOREDIST is shunted off to a separate function*/
      PStateChanges = initDoubleMatrix(PROTEIN_ALPHABET_SIZE);    
      PStationary = calcProteinStationaryProbs(data,numSeqs,numSites,filter,numSelected);
      scoreDist(PStationary,data,numSeqs,numSites,scoreMatrix,D,alphabet,scoreMatrixSize);
      free(PStationary);
      freeMatrix(PStateChanges,PROTEIN_ALPHABET_SIZE);
      return(D);
    default:
      fprintf(stderr,"\n Error: please specify model for sequence data\n");
      exit(EXIT_FAILURE);
    }
  for(i=0;i<numSeqs-1;i++){ /*i loop*/
    for(j=i;j<numSeqs;j++){ /*j loop*/
      if (i==j)
	D[i][j] = 0.0;
      else /*i ~=j */ {
	ijFilter(filter,data[i],data[j],itype,numSites);
	switch(model)
	  {
	  case K2P:
	    calcDNATransitionProbs(PStateChanges,data[i],data[j],numSites,filter,numSelected);
	    /*calculates PStateChanges matrix for data i,j*/
	    a = calcTransitionRate(PStateChanges);
	    b = calcTransversionRate(PStateChanges);
	    D[j][i] = D[i][j] = calcK2P(a,b,gamma);
	    break;
	  case F84:
	    calcDNATransitionProbs(PStateChanges,data[i],data[j],numSites,filter,numSelected);
	    a = calcTransitionRate(PStateChanges);
	    b = calcTransversionRate(PStateChanges);
	    D[i][j] = D[j][i] = calcF84(a,b,gamma,Mloc,Nloc,Ploc);
	    break; 
	  case TN93:
	    calcDNATransitionProbs(PStateChanges,data[i],data[j],numSites,filter,numSelected);
	    a = PStateChanges[ADENINE][GUANINE] + PStateChanges[GUANINE][ADENINE]; /*purine transition rate*/
	    a2 = PStateChanges[CYTOSINE][THYMINE] + PStateChanges[THYMINE][CYTOSINE]; /*pyrimidine transition rate*/
	    /*for testing purposes...*/
	    /*a2=a=calcTransitionRate(P);*/	    
	    b = calcTransversionRate(PStateChanges);
	    D[i][j] = D[j][i] = calcTN93(a,a2,b,PR,PY,Mloc,Nloc,gamma);
	    break;
	  case TRANSVERSIONONLY:
	    calcDNATransitionProbs(PStateChanges,data[i],data[j],numSites,filter,numSelected);
	    b = calcTransversionRate(PStateChanges);
	    D[i][j] = D[j][i] = XX(PR,PY,b,gamma);
	    break;
	  case LOGDET:
	    calcDNATransitionProbs(PStateChanges,data[i],data[j],numSites,filter,numSelected);
	    D[i][j] = D[j][i] = logdet(PStateChanges,Pi2[i],Pi2[j]);
	    break;
/*
	  case PROTEIN:
	    b = HammingDistance(data[i],data[j],filter,numSites,numSelected);
	    D[i][j] = D[j][i] = protFormula(b,gamma,pdiff);
	    break;	    
*/
	  } /*switch(model)*/
	if (D[i][j] > LARGE)
	  D[i][j] = D[j][i] = LARGE;
	if ( D[i][j] < minentry)
	  minentry = D[i][j];
	if ( D[i][j] > maxentry)
	  maxentry = D[i][j];
      }/*else*/
    }/*j loop*/
  }/*i loop*/
  free(PStationary);
/*
  if (PROTEIN == model)
    freeMatrix(PStateChanges,PROTEIN_ALPHABET_SIZE);
  else
*/
    freeMatrix(PStateChanges,DNA_ALPHABET_SIZE);
  if (LOGDET == model)
    freeMatrix(Pi2,numSeqs);  
  if (verbose)
    printf("  . Distance matrix: Min value is %lf, max value is %lf.\n",minentry,maxentry);
  return(D);
}

void symmetrizeDoubleMatrix(double **X, int n)
{
  int i,j;
  for(i=0;i<n-1;i++)
    for(j=i+1;j<n;j++)
      X[i][j]=X[j][i] = 0.5*(X[i][j] + X[j][i]);
  return;
}

