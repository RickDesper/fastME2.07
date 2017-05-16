#ifndef NEWICK_H_
#define NEWICK_H_

#include "graph.h"

node *decodeNewickSubtree(char *treeString, tree *T, int *uCount);
tree *readNewickString (char *str, int numLeaves);
tree *loadNewickTree(FILE *ifile, int numLeaves);
void NewickPrintSubtree(tree *T, edge *e, FILE *ofile);
void NewickPrintBinaryTree(tree *T, FILE *ofile);
void NewickPrintTrinaryTree(tree *T, FILE *ofile);
void NewickPrintTree(tree *T, FILE *ofile);
void NewickPrintSubtreeStr(tree *T, edge *e, char *str);
void NewickPrintBinaryTreeStr(tree *T, char *str);
void NewickPrintTrinaryTreeStr(tree *T, char *str);
void NewickPrintTreeStr(tree *T, char *str);
#endif /*NEWICK_H_*/
