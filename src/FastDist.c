#include "FastDist.h"

void ResolveAmbiguity( Sseq* seq1, Sseq* seq2, strIJ* Sij, int m, int model, int *filter)
{
  int i, j;
  int filterRepeat;
  double* pt;

  pt = ComputeTP(model, Sij, m);  //makes transition matrix

  for (i = 0; i != m; ++i)
    for (filterRepeat=0 ; filterRepeat < filter[i]; filterRepeat++) /*kludge fix: loop through this bit filter[i] times*/
      if ((filter[i] && (IndDNA(seq1->p[i]) == 5))) // ambiguous
	{
	  double t[5], totalprob;
	  int ii, ok;
	  int id;
	  double *amb, *AA2;
	  
	  for (ii = 0; ii != 4; ++ii) t[ii] = seq1->AmbVec[i][ii];
	  
	  id  = IndDNA(seq2->p[i]);
	  if (id >= 4) // unknown or ambiguous --> bye bye
	    continue;
	  
	  //reach this point only if seq1 has ambiguous character, seq2 does not
	  amb = seq1->AmbVec[i];
	  AA2 = seq2->AmbVec[i];
	  
	  ok = 1; // check amb[j] > 0 becomes 0
	  for (j = 0; j != 4; ++j)
	    if (AA2[j] > 0 && amb[j] == 0) ok = 0;  //reject if amb[j] == 0 for j=AA2[i] base
	  if (!ok) continue;
	  totalprob = 0;
	  for (j = 0; j != 4; ++j)
	    totalprob += amb[j] * pt[j*4 + id]; // P[j] * pt[j->id]      
	  
	  for (j = 0; j != 4; ++j)
	    amb[j] = (double) amb[j] * pt[j*4 + id] / totalprob;     
	}
  free(pt);
  return;
}

// compute puTs, pyTs, and Tv between ambiguity and  ambiguity
void CompAmbStatistic(double* tp, double* AmbVec1, double* AmbVec2, double *puTs,
	double *pyTs, double *Tv)
{
  double match;
  double puTs1, pyTs1, Tv1, total;
  int i, j;

  match = 0;
  for (i = 0; i != 4; ++i)
    match += AmbVec1[i] * AmbVec2[i] * tp[i*4+i];
  
  puTs1 = 0; pyTs1 = 0, Tv1 = 0;
  for (i  = 0; i != 4; ++i)
    for (j = 0; j != 4; ++j)
      if (i != j)
	{
	double temp;
	temp = AmbVec1[i] *   AmbVec2[j] * tp[i*4+j];
	if ((i == 0 && j == 2) || (i == 2 && j == 0)) { puTs1 += temp ; continue; }
	if ((i == 1 && j == 3) || (i == 3 && j == 1)) { pyTs1 += temp; continue; }	
	Tv1 += temp;
      }

  total = match + puTs1 + pyTs1 + Tv1;

  if (total == 0)
    {
      puTs1 = 0, pyTs1 = 0, Tv1 = 0;
      match = 0;
      for (i  = 0; i != 4; ++i)
	for (j = 0; j != 4; ++j)
	  if (i != j)
	  {
	    double temp;
	    temp = AmbVec1[i] *   AmbVec2[j];
	    if ((i == 0 && j == 2) || (i == 2 && j == 0)) { puTs1 += temp ; continue; }
	    if ((i == 1 && j == 3) || (i == 3 && j == 1)) { pyTs1 += temp; continue; }	
	    Tv1 += temp;
	  }
    }
  
  total = match + puTs1 + pyTs1 + Tv1;
  if (total > 1E-20)
    {
      (*puTs) = puTs1/total;
      (*pyTs) = pyTs1/total;
      (*Tv)   = Tv1/total;
    }
  if (total <= 1E-20)
    {
      (*puTs) = 0;
      (*pyTs) = 0;
      (*Tv)   = 0;
    }
  return;
}

void CorrectStatisticsIJ(strIJ* Sij, Sseq* seq1, Sseq* seq2, int m, int model, int *filter)
{
  int i;
  int filterRepeat;
  for (i = 0; i != m; ++i)
    for(filterRepeat=0;filterRepeat<filter[i];filterRepeat++)
      {
	int id1, id2;
	double puTs, pyTs, Tv;
	double *tp;
	
	id1 = IndDNA(seq1->p[i]);
	id2 = IndDNA(seq2->p[i]);
	if (id1 <= 3 && id2 <= 3) continue; // not ambiguous --> byebye
	if (id1 == 4 || id2 == 4) continue; // if unknown --> byebye
	
	tp = ComputeTP(model, Sij, m);
	
	CompAmbStatistic(tp, seq1->AmbVec[i], seq2->AmbVec[i], &puTs, &pyTs, &Tv);
	Sij->puTs += puTs;
	Sij->pyTs += pyTs;
	Sij->Tv   += Tv;
	Sij->del  --;
	free(tp);
      }
  UpdateProbSij(Sij, m);
  return;
}

//n=numSeqs m =numSites
double** FastDist(char** data, int n, int m, int model, double gamma, int *filter)
{
  int i, j;
  double **D;
  Sseq** Seq; //the sequences
  strIJ **StrIJ; //pairs of sequences
//  int *mark, tempnum, numID; //number of identical sites on the whole dataset

  int Num[4]; // Num {0, 1, 2, 3} = total # {A, C, G, T}
  int Warnings;

  Warnings = 0;
  D = initDoubleMatrix(n);
  Seq = (Sseq**)calloc(n, sizeof(Sseq*)); 
  StrIJ = (strIJ**)calloc(n*n, sizeof(strIJ*));
  //mark = (int*)calloc(m, sizeof(int));

  if (verbose)
    printf("  . Getting statistics of input data.\n");
  
  for (i = 0; i != n; ++i)
    Seq[i] = Statistic(data[i],  m, filter); // make statistic for sequence s
  
  for (i = 0; i != 4; ++i) Num[i] = 0;
  for (i = 0; i != n; ++i)
    for (j = 0; j != 4; ++j) Num[j] += Seq[i]->Num[j];
  //base counts

  for (i = 0; i != n; ++i)
    for (j = i+1; j != n; ++j)
      {
	StrIJ[i*n+j] = StatisticIJ(data[i], data[j], m, filter);
      }
  // i,j transition, transversion rates calculated ^ i< j v i>=j
  for (i = 0; i != n; ++i)
    for (j = 0; j != i+1; ++j)
      {
	StrIJ[i*n+j] = StatisticIJ(data[i], data[j], m, filter);
	if (StrIJ[i*n +j]->del * 2 > m) Warnings = 1;
      }
  
  if (verbose)
    printf("  . Start computing distance matrix of model TN93, F84, K2P or JC69.\n");
  
  //calculating distance matrix using analytic formula
  for (i = 0; i != n; ++i)
    for (j = i+1; j != n; ++j) 
      {
	double JCdist;
	
	D[i][j] = ComputeDistance(model, StrIJ[i*n+j], Num, m, gamma);
	D[j][i] = D[i][j];
	
	JCdist = ComputeDistance(JC69, StrIJ[i*n+j], Num, m, gamma);
	if (JCdist >= 5) Warnings = 2;
      	
	if (verbose)
	  printf("  . Distance %d %d = %f.\n", i, j, D[i][j]);
      }

  //for (i = 0;i != n; ++i) D[i][i] = 0;
  //zeros along diagonal (is this necessary?) InitDoubleMatrix sets initial values to zero already

  // compute ambiguity
  if (verbose)
    printf("  . Computing Ambiguities probabilities using nearest neighbor.\n");
  for (i = 0; i != n; ++i)
    {
      int k;
      int hasAmbiguity;
      hasAmbiguity = 0;
      for (k = 0; k!= m; ++k)
	{
	  int ind;
	  ind = IndDNA(Seq[i]->p[k]);
	  if ((ind == 5) && filter[k]) 
	    {
	      hasAmbiguity = 1;
	      break; 
	    }
	}
      if(hasAmbiguity == 1)
	{
	  j = NearestNeighbor(i, D, n);
	  if (verbose)
	    printf("  . Nearest Neightbor of %d is %d distance = %15.10f.\n", i, j, D[i][j]);
	  ResolveAmbiguity(Seq[i], Seq[j], StrIJ[i*n+j], m, model, filter);
	}
    }

  if(verbose)
    printf("  . Correcting statistics using ambiguity probabilities.\n");
  for (i = 0; i != n; ++i)
    for (j = 0; j != n; ++j)
      if (i != j)
	{
	  int k;
	  int hasAmbiguity;
	  hasAmbiguity = 0;
	  for (k = 0; k!= m; ++k)
	    {
	      int ind1 = IndDNA(Seq[i]->p[k]);
	      int ind2 = IndDNA(Seq[j]->p[k]);
	      if (ind1 == 5 || ind2 == 5) 
		{ 
		  hasAmbiguity = 1;
		  break;
		}
	    }
	  
	  if(hasAmbiguity == 1)
	    {
	      CorrectStatisticsIJ(StrIJ[i*n+j], Seq[i], Seq[j], m, model, filter);
	    }
	}

  if(verbose)
    printf("  . Finalizing distance computation.\n");
  for (i = 0; i != n; ++i)
    for (j = i+1; j != n; ++j)
	{
	  int k;
	  int hasAmbiguity = 0;
	  for (k = 0; k!= m; ++k)
	    if (filter[k])
	      {
		int ind1 = IndDNA(Seq[i]->p[k]);
		int ind2 = IndDNA(Seq[j]->p[k]);
		if (ind1 == 5 || ind2 == 5) 
		  { 
		    hasAmbiguity = 1;
		    break;
		  }
	      }
	  if(hasAmbiguity == 1)
	    {
	      
	      if (StrIJ[i*n +j]->del == m)
		{
		  D[i][j] = 1001;
		}
	      else 
		D[i][j] = ComputeDistance(model, StrIJ[i*n+j], Num, m, gamma);
	      
	      if (verbose)
	        printf("  . i = %d j = %d  strlen = %d transition = %f transversion = %f deleted = %d dist = %f.\n", 
	      			  i, j,  m, StrIJ[i*n + j]->puTs+ StrIJ[i*n + j]->pyTs, StrIJ[i*n + j]->Tv, StrIJ[i*n + j]->del, D[i][j]);
	      
	      D[j][i] = D[i][j];
	    }
	}

  //Free pointers 
  for (i = 0; i != n; ++i)
    Free(Seq[i], m);
  free(Seq);
  for (i = 0;i != n*n; ++i)
    free(StrIJ[i]);
  free(StrIJ);
  
  for (i = 0; i != n; ++i)
    D[i][i] = 0.0;
  
  for (i = 0; i != n; ++i)
    for (j = 0; j != n; ++j)
      {
	if (D[i][j] < 0)
	  {
	    printf("\n Warning: Distance between %d %d is too large, putting 5 in place.\n", i, j);
	    D[i][j] = 5;
	  }
	if (D[i][j] > 5)
	  {
	    printf("\n Warning: Distance between %d %d is too large, putting 5 in place.\n", i, j);
	    D[i][j] = 5;
	  }
	
      }
  if(Warnings == 1)
    printf("\n Warning: Very long gaps in this data set.\n");
  if(Warnings == 2)
    printf("\n Warning: Give up this dataset because JC69 dist >5.\n");
  return D;
}

void UpdateProbSij(strIJ* Sij, int m)
{
  double strlen;
  strlen = m - Sij -> del;
  Sij->ts_purine_prob     = (double) Sij->puTs/strlen;
  Sij->ts_pyrimidine_prob = (double) Sij->pyTs/strlen;
  Sij->tv_prob            = (double) Sij->Tv/strlen;
  Sij->idprob             = 1 - (Sij->ts_purine_prob + Sij->ts_pyrimidine_prob +  Sij->tv_prob);
  return;
}

//StrIJ[i*n+j] = StatisticIJ(data[i], data[j], m);
strIJ* StatisticIJ(char* s1, char* s2, int m, int *filter)
{
  int i;
  strIJ *Sij;
  Sij = (strIJ*)malloc(sizeof(strIJ)); //RD: calloc -> malloc

  Sij->puTs = 0;
  Sij->pyTs = 0;
  Sij->Tv   = 0;
  Sij->del  = 0;

  for (i = 0; i != m; ++i)
    {
      int id1, id2;
      id1 = IndDNA(s1[i]); id2 = IndDNA(s2[i]);
      
      if (id1 >=4 || id2 >= 4) { Sij->del += 1; continue; } // ambiguous and unknown
      if (id1 == id2) continue;

      if ( (id1 == 0 && id2 == 2) || (id1 == 2 && id2 == 0)) 
	{ 
	  Sij->puTs += filter[i]; 
	  continue;
	} // A->G, G->A

      if ( (id1 == 1 && id2 == 3) || (id1 == 3 && id2 == 1)) 
	{ 
	  Sij->pyTs += filter[i]; 
	  continue;
	} // C->T, T->C
      Sij->Tv += filter[i];     // others
    }

  UpdateProbSij(Sij, m);

  return Sij;
}

/*

0 <--> A
1 <--> C
2 <--> G
3 <--> T
4 <--> Unknown
5 <--> Amb

 */

int IndDNA(char c) /*returns numeric value for character 'c' 
		     a=A=0,c=C=1,g=G=2,t=T=3,.=_=4
		     other characters return 5*/
{
  char DNA_A[10];
  int i;

  strncpy(DNA_A,"ACGT-\0",6);

  if (c >= 'a' && c <= 'z') c -= 32;
  if (c == '.') c= '-';
  for (i = 0;i != 5; ++i)
    if (DNA_A[i] == c) return i;

  return 5;
}

Sseq*  Statistic(char* s, int m, int *filter) // make statistic for sequence s
     //RD: by "statistic" Hang means the count for each possible base
{
  Sseq *S;
  int i;

  S = (Sseq*)calloc(1, sizeof(Sseq));

  S->p = s;

  for (i = 0; i != 10; ++i)    S->Num[i] = 0;
  for (i = 0; i != m; ++i)     S->Num[IndDNA(s[i])] += filter[i];

  // init AmbVec[m][4]
  S->AmbVec = (double**)calloc(m, sizeof(double*));
  for (i = 0; i != m; ++i) S->AmbVec[i] = NULL;
		    
  for (i = 0; i != m; ++i)
    if (IndDNA(s[i]) != 4)
      S->AmbVec[i] = GetFre(s[i]);

  return S;
}

double *GetFre(char c) /*returns frequency vector for ambiguity character*/
{
  char DNA_Amb[20]; 
  double const v13= (double)1/3;
  double Fre[16][4] = {
    {1.0, 0.0, 0.0, 0.0}, // A
    {0.0, 1.0, 0.0, 0.0}, // C
    {0.0, 0.0, 1.0, 0.0}, // G
    {0.0, 0.0, 0.0, 1.0}, // T
    {0.5, 0.5, 0.0, 0.0}, // M
    {0.5, 0.0, 0.5, 0.0}, // R
    {0.5, 0.0, 0.0, 0.5}, // W
    {0.0, 0.5, 0.5, 0.0}, // S
    {0.0, 0.5, 0.0, 0.5}, // Y
    {0.0, 0.0, 0.5, 0.5}, // K
    {v13, v13, v13, 0.0}, // V
    {v13, 0.0, v13, v13}, // D
    {0.0, v13, v13, v13}, // B
    {.25, .25, .25, .25}, // N
    {v13, v13, 0.0, v13},  // H
    {1,1,1,1}
  };

  int i, li;
  double *t;

  li = -1;
  strncpy(DNA_Amb, "ACGTMRWSYKVDBNH\0", 16);
  if (c >= 'a' && c <= 'z') c -= 32;
  if (c == 'X' || c == '?') c = 'N';

/*  for (i = 0; i < (strlen(DNA_Amb)); ++i)*/
  for(i=0; i < 15; ++i)
    {
      if (DNA_Amb[i] == c)  { li = i; break;}
      //      if (c == 'H') printf("%c\n", DNA_Amb[i]);
    }

  if (li == -1)
    {
      printf("\n What is [%c, %d]???\n", c, c);
      exit(EXIT_FAILURE);
    }

  t = (double*)calloc(4, sizeof(double));
  for (i = 0; i != 4; ++i)
    t[i] = Fre[li][i];

  return t;
}

void Free(Sseq* S, int m)
{
  int i;
  for (i = 0; i != m; ++i)
    free(S->AmbVec[i]);
  free(S->AmbVec);

  free(S);
  return;
}

int NearestNeighbor(int i, double** D, int n )
{
  int j , k;
  double Min;
  j = 0;

  if (i == 0) j = 1;
  Min = 10000000.0;

  for (k = 0; k != n; ++k)
    if (k != i)
      if (D[i][k] < Min && D[i][k] >=0) {
	j = k;
	Min = D[i][k];
      } 
  return j;
}

