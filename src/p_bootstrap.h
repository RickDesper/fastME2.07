
#ifndef P_BOOTSTRAPS_H
#define P_BOOTSTRAPS_H

#include "p_utils.h"


void boot(char *bestTree, char **bootTree, int n_boot, FILE *treeFile);
arbre *Read_Tree(char *s_tree);
arbre *Make_Tree(int n_otu);
void Init_Tree(arbre *tree, int n_otu);
void Make_All_Tree_Nodes(arbre *tree);
bnode *Make_Node_Light(int num);
void Init_Node_Light(bnode *n, int num);
void Make_All_Tree_Edges(arbre *tree);
bedge *Make_Edge_Light(bnode *a, bnode *d, int num);
void Init_Edge_Light(bedge *b, int num);
void Make_Edge_Dirs(bedge *b, bnode *a, bnode *d);
//void Make_Tree_Path(arbre *tree);
//void Make_List_Of_Reachable_Tips(arbre *tree);
char **Sub_Trees(char *tree, int *degree);
int Next_Par(char *s, int pos);
void Clean_Multifurcation(char **subtrees, int current_deg, int end_deg);
void Unroot_Tree(char **subtrees);
void R_rtree(char *s_tree_a, char *s_tree_d, bnode *a, arbre *tree, int *n_int, int *n_ext);
void Read_Branch_Label(char *s_d, char *s_a, bedge *b);
void Make_New_Edge_Label(bedge *b);
void Read_Branch_Length(char *s_d, char *s_a, arbre *tree);
void Connect_One_Edge_To_Two_Nodes(bnode *a, bnode *d, bedge *b, arbre *tree);
void Read_Node_Name(bnode *d, char *s_tree_d, arbre *tree);
void Alloc_Bip(arbre *tree);
void Get_Bip(bnode *a, bnode *d, arbre *tree);
int Sort_String(const void *a, const void *b);
void Compare_Bip(arbre *tree1, arbre *tree2);
void Free_Tree(arbre *tree);
void Free_bNode(bnode *n);
void Free_Edge_Labels(bedge *b);
void Free_Edge(bedge *b);
char *Write_Tree(arbre *tree);
void R_wtree(bnode *pere, bnode *fils, char *s_tree, arbre *tree);


#endif /*P_BOOTSTRAPS_H_*/

