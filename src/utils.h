#ifndef UTILS_H_
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fcntl.h>
//#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
//#include <libgen.h>
//#include <sys/types.h>
//#include <sys/stat.h>

typedef char boolean;
boolean verbose;

#ifndef VERSION
#define VERSION		"v2.07"
#endif

#ifndef true
#define true 1
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef false
#define false 0
#endif
#ifndef FALSE
#define FALSE 0
#endif

//constants for input sequences
#ifndef MAX_LINE_LENGTH
#define MAX_LINE_LENGTH 2000000
#endif
#ifndef MAX_SEQ_LENGTH
#define MAX_SEQ_LENGTH 2000000
#endif

//constants for tree nodes and edges
#ifndef MAX_LABEL_LENGTH
#define MAX_LABEL_LENGTH 20
#endif
#ifndef MAX_FILE_NAME_LENGTH
#define MAX_FILE_NAME_LENGTH 256
#endif
#ifndef MAX_EVENT_NAME
#define MAX_EVENT_NAME 24
#endif

#ifndef INFTY
#define INFTY 10000000
#endif
#ifndef NEGINFTY
#define NEGINFTY -10000000
#endif
#ifndef NONE
#define NONE 0
#endif
#ifndef UP
#define UP 1
#endif
#ifndef DOWN
#define DOWN 2
#endif
#ifndef LEFT
#define LEFT 3
#endif
#ifndef RIGHT
#define RIGHT 4
#endif
#ifndef SKEW
#define SKEW 5
#endif

#ifndef ReadOpenParenthesis
#define ReadOpenParenthesis 0
#endif
#ifndef ReadSubTree
#define ReadSubTree 1
#endif
#ifndef ReadLabel
#define ReadLabel 2
#endif
#ifndef ReadWeight
#define ReadWeight 3
#endif
#ifndef AddEdge
#define AddEdge 4
#endif
#ifndef ReadSize
#define ReadSize 5
#endif
#ifndef ReadEntries
#define ReadEntries 6
#endif
#ifndef Done
#define Done 7
#endif

//constant for Newick strings
#ifndef INPUT_SIZE
#define INPUT_SIZE 100
#endif
#ifndef MAX_INPUT_SIZE
#define MAX_INPUT_SIZE 10000000
#endif
#ifndef MAX_DIGITS
#define MAX_DIGITS 20
#endif

//maximum distance matrix size
#ifndef MAXSIZE
#define MAXSIZE 70000
#endif

#ifndef OLS
#define OLS 1
#endif
#ifndef BAL
#define BAL 2
#endif
#ifndef BALSPR
#define BALSPR 3
#endif
#ifndef NJ
#define NJ 5
#endif
#ifndef UNJ
#define UNJ 6
#endif
#ifndef BIONJ
#define BIONJ 7
#endif

#ifndef NONE
#define NONE 9
#endif
#ifndef USER
#define USER 10
#endif

//constant for minimum improvement...to avoid looping problem caused
//  by round-off error
#ifndef EPSILON
#define EPSILON .0000001
#endif
#ifndef LARGE
#define LARGE 3.0
#endif


#ifndef DNA_ALPHABET_SIZE
#define DNA_ALPHABET_SIZE 4
#endif
#ifndef DNA_ALPHABET
#define DNA_ALPHABET "ACGT"
#endif

#ifndef PROTEIN_ALPHABET_SIZE
#define PROTEIN_ALPHABET_SIZE 26
#endif
#ifndef PROTEIN_ALPHABET
#define PROTEIN_ALPHABET "ABCDEFGHIKLMNPQRSTVWYZX*?-"
#endif

//constants for DNA alphabet
#ifndef ADENINE
#define ADENINE 0
#endif
#ifndef CYTOSINE
#define CYTOSINE 1
#endif
#ifndef GUANINE
#define GUANINE 2
#endif
#ifndef THYMINE
#define THYMINE 3
#endif

//nucleotides model constants
#ifndef F84
#define F84 11
#endif
#ifndef TN93
#define TN93 12
#endif
#ifndef K2P
#define K2P 13
#endif
#ifndef JC69
#define JC69 14
#endif
#ifndef TRANSVERSIONONLY
#define TRANSVERSIONONLY 15
#endif
#ifndef LOGDET
#define LOGDET 16
#endif
//input data type constants
#ifndef MATRIX
#define MATRIX 21
#endif
#ifndef DNA
#define DNA 22
#endif
#ifndef PROTEIN
#define PROTEIN 23
#endif
#ifndef SCOREDIST
#define SCOREDIST 24
#endif
//protein models constants
#ifndef WAG
#define WAG 31
#endif
#ifndef DAYHOFF
#define DAYHOFF 32
#endif
#ifndef JTT
#define JTT 33
#endif
#ifndef BLOSUM62
#define BLOSUM62 34
#endif
#ifndef MTREV
#define MTREV 35
#endif
#ifndef RTREV
#define RTREV 36
#endif
#ifndef CPREV
#define CPREV 37
#endif
#ifndef DCMUT
#define DCMUT 38
#endif
#ifndef VT
#define VT 39
#endif
#ifndef LG
#define LG 40
#endif

int *initZeroArray(int l);
int *initOneArray(int l);
double **initDoubleMatrix(int d);
boolean whiteSpace(char c);
void Exit (char *message);
void Warning (char *message);

#endif /*UTILS_H_*/

