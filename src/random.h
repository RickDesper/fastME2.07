#ifndef RANDOM_H_
#define RANDOM_H_

#include "utils.h"

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1 -1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >18)

void sgenrand(unsigned long seed);
double uniformGenerator(int *junk);
int getIntRandom(int *idum, int maxValue);
int getWeightedIntRandom(int *idum, int maxValue, double *weights);
double getMatrixMean(double **mat, double **P, int height, int width);
double StandardExponential(int *idum);
double transformedExponential(int *idum, double scale, double lag);
double StandardGaussian(int *idum);
double transformedGaussian(double mean, double stddev, int *idum);
boolean coinToss(double p, int *idum);
void subsetSelect(double p, int seqlength, int *selected, int *seed);
void bootstrapSelect(int seqlength, int *selected, int *seed);

#endif /*RANDOM_H_*/

