#ifndef HEAP_H_
#define HEAP_H_

#include "utils.h"

int *initPerm(int size);
void permInverse(int *p, int *q, int length);
void swap(int *p, int *q, int i, int j);
void heapify(int *p, int *q, double *HeapArray, int i, int n);
void reHeapElement(int *p, int *q, double *v, int length, int i);
void popHeap(int *p, int *q, double *v, int length, int i);
void pushHeap(int *p, int *q, double *v, int length, int i);
int makeThreshHeap(int *p, int *q, double *v, int arraySize, double thresh);

#endif /*HEAP_H_*/

