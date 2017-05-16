#include "CompDist.h"

double ComputeDistance(int model, strIJ* Sij, int* Num, int m, double gamma) // computing distance between i, j; Num[0,1,2,3] = {#A, #C, #G, #T}
     //returns analytic formula distance using one of the four models
{
  //  if (rate >= 0.75) printf("WARNING: Distance between two sequence is too large\n");
  //printf("puTs = %f pyTs = %f del = %d\n", Sij->puTs, Sij->pyTs, Sij->del);
  if(Sij->del == m) return 5;
  if(Sij->puTs == 0 && Sij->pyTs ==0 && Sij->Tv ==0 && Sij->del < m) return 0;

  if (model == TN93)  return  Compute_TN93(Sij, Num, m, gamma);
  if (model == K2P)   return  Compute_K2P(Sij, m, gamma);
  if (model == F84)   return  Compute_F84(Sij, Num, m, gamma);
  if (model == JC69)  return  Compute_JC69(Sij, m, gamma);

  Exit ("Model must be TN93, K2P, F84 or JC69");
  return -1;
}

//-----------------------------------------------
//TAMURA NEI - regular
//ML_string_distance compute_TN93(int strlen, double gamma, TN_string_distance sd,
//			      int numAs, int numCs, int numGs, int numTs){
//RD:formula verified
double Compute_TN93(strIJ* Sij, int* Num, int m, double gamma) // computing distance between i, j; Num[0,1,2,3] = {#A, #C, #G, #T}
{
  int strlen;
  int i;
  double norm, piA, piC, piG, piT, piR, piY;
  double ts_purine_prob, ts_pyrimidine_prob, tv_prob ;
  double loc1, loc2, loc3, loc4, distance;
 
  strlen = m - Sij->del;
  norm = 0;

  for (i = 0; i != 4; ++i) norm += Num[i];

  piA  = 0.000001+((double)Num[0])/norm;
  piC = 0.000001+((double)Num[1])/norm;
  piG = 0.000001+((double)Num[2])/norm;
  piT = 0.000001+((double)Num[3])/norm;

  piR = piA+piG;
  piY = piC+piT;

  ts_purine_prob     = Sij->ts_purine_prob;
  ts_pyrimidine_prob = Sij->ts_pyrimidine_prob;
  tv_prob            = Sij->tv_prob;

  if (verbose)
    printf("  . TN parameters: pR = %f pY = %f ts_purine_prob= %f ts_pyrimidine_prob= %f tv_prob= %f.\n", piR, piY, ts_purine_prob, ts_pyrimidine_prob, tv_prob);

  loc1 = 1.0 -  piR*ts_purine_prob/(2.0*piA*piG) - tv_prob/(2.0*piR);
  loc2 = 1.0  -  piY*ts_pyrimidine_prob/(2.0*piC*piT) - tv_prob/(2.0*piY);
  loc3 = piR*piY - piA*piG*piY/piR - piT*piC*piR/piY;
  loc4 = 1.0 - tv_prob/(2.0*piR*piY);

  distance = 0;

  if (gamma == 0)
    distance = - 2.0 * (piA*piG/piR)   *log(loc1)
               - 2.0 * (piT* piC/piY) *log(loc2)
               - 2.0 * loc3 * log(loc4);

  if (gamma > 0)
    distance = 2.0 * gamma * (piA*piG/piR) * pow(loc1, -1/gamma)
      + 2.0 * gamma * (piT * piC/piY) *pow(loc2, -1/gamma)
      + 2.0 * loc3 * pow(loc4, -1/gamma)
      - 2.0 * gamma * (piA * piG + piC * piT + piR * piY);
  
  return distance;
}

double Compute_JC69(strIJ* Sij, int m, double gamma) 
{
  double P, loc, distance; 
  P =  Sij->ts_purine_prob + Sij->ts_pyrimidine_prob + Sij->tv_prob; 
  loc = 1 - 4.0 / 3.0 * P;
  if (verbose)
    printf("  . JC69 Parameter : P = %f.\n", P);
  
  distance = -(3.0/4.0) * log(loc);
  if (gamma > 0)
    distance = 0.75 * gamma * (pow(loc, -1/gamma) - 1);
  
  return distance;
}

//verified
double Compute_K2P(strIJ* Sij, int m, double gamma)
{
  double tsprob, tvprob, distance;

  tsprob = Sij->ts_purine_prob + Sij->ts_pyrimidine_prob;  //sd.transitions/strlen;
  tvprob = Sij->tv_prob; //sd.transversions/strlen;

  if (verbose)
    printf("  . K2P parameters a = %f b = %f.\n", tsprob, tvprob);

  distance = - 0.5 * log( (1.0 - 2.0*tsprob-tvprob)*sqrt(1.0-2.0*tvprob) );
  
  if (gamma !=0)    
    distance = gamma * (0.5 * pow(1.0 - 2.0 * tsprob - tvprob, -1/gamma) + 0.25 * pow(1.0 - 2.0 * tvprob, -1/gamma) - 0.75);

  return distance;
}		   

double Compute_F84(strIJ* Sij, int* Num, int m, double gamma)
{
  int strlen; 
  int i;
  double norm, piA, piC, piG, piT, piR, piY;
  double M, N, P, A, B, loc1, loc2, distance;

  strlen = m - Sij->del;
  norm = 0;

  for (i = 0; i != 4; ++i) norm += Num[i];

  piA = 0.000001+((double)Num[0])/norm;
  piC = 0.000001+((double)Num[1])/norm;
  piG = 0.000001+((double)Num[2])/norm;
  piT = 0.000001+((double)Num[3])/norm;

  piR = piA+piG;
  piY = piC+piT;

  M   = piC*piT/piY + piA*piG/piR; //M
  N   = piC*piT + piA*piG;         //N
  P   = piR*piY;                   //P
  
  A = Sij->ts_purine_prob + Sij->ts_pyrimidine_prob;
  B = Sij->tv_prob;

  if (verbose)
    printf("  . F84 parameter: A,B,M, N, P  m %f %f %f %f %f %d transition %f transversion %f.\n", A, B, M ,N, P, m, Sij->puTs + Sij->pyTs, Sij->Tv);

  loc1 = 1.0  - A / (2.0 * M) - (M - N) * B /( 2.0 * M * P);
  loc2 = 1.0 - B /(2.0 * P);

  distance = -2.0 * M * log(loc1)- 2.0 * (N + P - M) * log(loc2);
 
  if (gamma > 0.0)
    distance = 2.0* gamma * ( M * pow(loc1, -1/gamma) + (N + P - M) * pow(loc2, -1/gamma) - N - P);
  return distance;
}

double* ComputeTP(int model, strIJ * Sij, int m)
{
  //  tp.pAA = 0
  //  tp.pAC = 1
  //  tp.pAG = 2
  //  tp.pAT = 3
  //  tp.pCC = 5
  //  tp.pCG = 6
  //  tp.pCT = 7
  //  tp.pGG = 10
  //  tp.pGT = 11
  //  tp.pTT = 15

  double ts_purine_prob, ts_pyrimidine_prob, tv_prob, idprob;
  
  ts_purine_prob     =  Sij->ts_purine_prob;
  ts_pyrimidine_prob =  Sij->ts_pyrimidine_prob;
  tv_prob            =  Sij->tv_prob;
  idprob             =  Sij->idprob;
 
  if (model == TN93)   return Compute_TP_TN93(ts_purine_prob, ts_pyrimidine_prob, tv_prob);
  if (model == K2P)    return Compute_TP_K2P(ts_purine_prob, ts_pyrimidine_prob, tv_prob, idprob);
  if (model == F84)    return Compute_TP_F84(ts_purine_prob, ts_pyrimidine_prob, tv_prob, idprob);
  if (model == JC69)   return Compute_TP_JC69(ts_purine_prob, ts_pyrimidine_prob, tv_prob, idprob);
  return NULL;
}

double* Compute_TP_TN93(double ts_purine_prob, double ts_pyrimidine_prob, double tv_prob) 
{
  double* out;
  out = (double*)calloc(16, sizeof(double));
  
  out[0] = 1.0-ts_purine_prob-tv_prob; 
  out[1] = tv_prob*0.5;  
  out[2] = ts_purine_prob; 
  out[3] = tv_prob*0.5;
  
  out[5] = 1.0-ts_pyrimidine_prob-tv_prob;  
  out[6] = tv_prob*0.5; 
  out[7] = ts_pyrimidine_prob;

  out[10]= 1.0-ts_purine_prob-tv_prob;  
  out[11] = tv_prob*0.5;
  out[15] = 1.0-ts_pyrimidine_prob-tv_prob;
  
  out[4] = out[1];
  out[8] = out[2];
  out[9] = out[6];
  out[12] = out[3];
  out[13] = out[7];
  out[14] = out[11];

  return out;
}

double* Compute_TP_K2P(double ts_purine_prob, double ts_pyrimidine_prob, double tv_prob, double idprob) 
{
  double* out; 
  double tsprob, tvprob;
  
  out = (double*)calloc(16, sizeof(double));
  tsprob = ts_purine_prob + ts_pyrimidine_prob;
  tvprob = tv_prob;

  out[0]  = idprob; 
  out[1]  = tvprob;
  out[2]  = tsprob;
  out[3]  = tvprob;
  
  out[5]  = idprob;
  out[6]  = tvprob;
  out[7]  = tsprob; 

  out[10] = idprob;
  out[11] = tvprob;
  out[15] = idprob;
  
  out[4] = out[1];
  out[8] = out[2];
  out[9] = out[6];
  out[12] = out[3];
  out[13] = out[7];
  out[14] = out[11];

  /*
  tp.pAA = idprob; 
  tp.pAC = tvprob;  
  tp.pAG = tsprob; 
  tp.pAT = tvprob;
  tp.pCC = idprob;  
  tp.pCG = tvprob; 
  tp.pCT = tsprob;
  tp.pGG = idprob;  
  tp.pGT = tvprob;
  tp.pTT = idprob;
  */
  return out;
}

double *Compute_TP_F84(double ts_purine_prob, double ts_pyrimidine_prob, double tv,  double idprob) 
{
  double* out; 
  double tsprob, tvprob;

  out = (double*)calloc(16, sizeof(double));

  tsprob = ts_purine_prob + ts_pyrimidine_prob;
  tvprob = tv/2;//(1.0 - idprob -tsprob)/2.0;

  out[0]  = idprob; //AA 
  out[1]  = tvprob; //AC
  out[2]  = tsprob; //AG
  out[3]  = tvprob; //AT
   
  out[5]  = idprob; //CC
  out[6]  = tvprob; //CG 
  out[7]  = tsprob; //CT

  out[10] = idprob; //GG
  out[11] = tvprob; //GT
  
  out[15] = idprob; //TT
  
  out[4] = out[1];
  out[8] = out[2];
  out[9] = out[6];
  out[12] = out[3];
  out[13] = out[7];
  out[14] = out[11];

  return out;
}
 
double *Compute_TP_JC69(double ts_purine_prob, double ts_pyrimidine_prob, double tv_prob, double idprob)
{
  double* out;
  double changeprob;

  out = (double*)calloc(16, sizeof(double));
  changeprob = (ts_purine_prob + ts_pyrimidine_prob +tv_prob)/3.0;

  out[0]  = idprob; 
  out[1]  = changeprob;
  out[2]  = changeprob;
  out[3]  = changeprob;
  
  out[5]  = idprob;
  out[6]  = changeprob;
  out[7]  = changeprob; 

  out[10] = idprob;
  out[11] = changeprob;
  
  out[15] = idprob;
  
  out[4] = out[1];
  out[8] = out[2];
  out[9] = out[6];
  out[12] = out[3];
  out[13] = out[7];
  out[14] = out[11];

  return out;
}

