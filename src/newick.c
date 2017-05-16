#include "newick.h"

int nodeCount;
int edgeCount;

/*decodeNewickSubtree is used to turn a string of the form
  "(v1:d1,v2:d2,(subtree) v3:d3....vk:dk) subroot:d," into a subtree
  rooted at subrooted, with corresponding subtrees and leaves at v1
  through vk.  It is called recursively on subtrees*/
node *decodeNewickSubtree(char *treeString, tree *T, int *uCount)
{
  node *thisNode = NULL;
  node *centerNode;
  double thisWeight = 0;
  edge *thisEdge;
//  char label[MAX_LABEL_LENGTH];
  char stringWeight[MAX_LABEL_LENGTH];
  int state;
  int i = 0;
  int j;
  int left,right;
  int parcount;
//  snprintf(label, 13, "Default_Label");
  left = right = 0;
  parcount = 0;
  state = ReadOpenParenthesis;
  if('(' == treeString[0])
    parcount++;
//  centerNode = makeNode(label,0);
//  centerNode = makeNode("",0);
  centerNode = makeNode(" ",nodeCount++);
  T->size++;
  while(parcount > 0)
    {
      while(whiteSpace(treeString[i]))
	i++;
      switch(state) 
	{
	case(ReadOpenParenthesis):
	  if('(' != treeString[0])
	    Exit ("Invalid Newick string.");
	  i++;
	  state = ReadSubTree;
	  break;
	case(ReadSubTree):
	  if('(' == treeString[i])  /*if treeString[i] is a left parenthesis,
				      we scan down the string until we find its partner.
				      the relative number of '('s vs. ')'s is counted
				      by the variable parcount*/
	    {
	      left = i++;
	      parcount++;
	      while(parcount > 1)
		{
		  while (('(' != treeString[i]) && (')' != treeString[i]))
		    i++;  /*skip over characters which are not parentheses*/
		  if('(' == treeString[i])
		    parcount++;
		  else if (')' == treeString[i])
		    parcount--;
		  i++;
		}  /*end while */
	      right = i;  /*at this point, the subtree string goes from 
			    treeString[left] to treeString[right - 1]*/
	      thisNode = decodeNewickSubtree(treeString + left,T,uCount);  /*note that this
								      step will put 
								      thisNode in T*/
	      i = right;  /*having created the node for the subtree, we move
			    to assigning the label for the new node.
			    treeString[right] will be the start of this label */
	    } /* end if ('(' == treeString[i]) condition */
	  else
	    {	  
//	      thisNode = makeNode(label,nodeCount++);
	      thisNode = makeNode(" ",nodeCount++);
	      T->size++;
	    }
	  state = ReadLabel;
	  break;
	case(ReadLabel):
	  left = i;  /*recall "left" is the left marker for the substring, "right" the right*/
	  if (':' == treeString[i])   /*this means an internal node?*/
	    {
	      //sprintf(thisNode->label,"");
	      snprintf(thisNode->label,1," ");
	      right = i;
	    }
	  else
	    {
	      while((':' != treeString[i]) && (',' != treeString[i]) && (')' != treeString[i]))
		i++;
	      right = i;
	      j = 0;
	      for(i = left; i < right;i++)
		if(!(whiteSpace(treeString[i])))
		  thisNode->label[j++] = treeString[i];
	      thisNode->label[j] = '\0';
	    }	      
	  if(':' == treeString[right])
	    state = ReadWeight;
	  else
	    {
	      state = AddEdge;
	      thisWeight = 0.0;
	    }
	  i = right + 1;
	  break;
	case(ReadWeight):
	  left = i;
	  while
	    ((',' != treeString[i]) && (')' != treeString[i]))
	    i++;
	  right = i;
	  j = 0;
	  for(i = left; i < right; i++)
	    stringWeight[j++] = treeString[i];
	  stringWeight[j] = '\0';
	  thisWeight = atof(stringWeight);
	  state=AddEdge;
	  break;
	case(AddEdge):
//	  thisEdge = makeEdge(label,centerNode,thisNode,thisWeight);
	  thisEdge = makeEdge(" ",centerNode,thisNode,thisWeight);
	  thisNode->parentEdge = thisEdge;
	  if (NULL == centerNode->leftEdge)
	    centerNode->leftEdge = thisEdge;
	  else if (NULL == centerNode->rightEdge)
	    centerNode->rightEdge = thisEdge;
	  else if (NULL == centerNode->middleEdge)
	    centerNode->middleEdge = thisEdge;
	  else
	    {
	      fprintf(stderr,"\n Error: node %s has too many (>3) children.\n",centerNode->label);
	      exit(EXIT_FAILURE);
	    }
	  snprintf(thisEdge->label,MAX_LABEL_LENGTH,"E%d",edgeCount++);
//	  edgeCount++;
	  i = right + 1;
	  if (',' == treeString[right])
	    state = ReadSubTree;
	  else
	    parcount--;
	  break;
	}
    }
  if (NULL != centerNode->parentEdge && verbose)
    Warning ("Decode Newick sub-tree routine exception.");
  return(centerNode);
}

tree *readNewickString (char *str, int numLeaves)
{
  tree *T;
  node *centerNode;
  int i = 0;
  int j = 0;
  int inputLength;
  int uCount = 0;
  int parCount = 0;
  char rootLabel[MAX_LABEL_LENGTH];
  nodeCount = edgeCount = 0;
  
  T = newTree();

  if ('(' != str[0])
    Exit ("Error reading generated tree - does not start with '('.");
  inputLength = strlen (str)+1;
  for(i = 0; i < inputLength; i++)
    {
      if ('(' == str[i])
	parCount++;
      else if (')' == str[i])
	parCount--;
      if (parCount > 0)
	;
      else if (0 == parCount)
	{
	  i++;
/*	  if(';' == str[i])
	    snprintf(rootLabel,MAX_LABEL_LENGTH,"URoot");
	  else
	    {*/
//	      while(';' != str[i])
	      while((';' != str[i]) && (! (whiteSpace (str[i]))) && (j < MAX_LABEL_LENGTH)) {
//	        if (!(whiteSpace (str[i++])))
//	          rootLabel[j++] = str[i-1];  /*be careful here */
	          rootLabel[j++] = str[i++];
	      }
	      rootLabel[j] = '\0';
//	    }
	  i = inputLength;
	}
      else if (parCount < 0)
        Exit ("Error reading generated tree. Too many right parentheses.");
    }
  centerNode = decodeNewickSubtree (str, T, &uCount);
  snprintf (centerNode->label, MAX_LABEL_LENGTH, rootLabel);
  T->root = centerNode;
  return (T);
}

tree *loadNewickTree(FILE *ifile, int numLeaves)
{
//  char label[] = "EmptyEdge";
  tree *T;
  node *centerNode;
  int i = 0;
  int j = 0;
  int inputLength;
  int uCount = 0;
  int parCount = 0;
  char c;
  char rootLabel[MAX_LABEL_LENGTH];
  boolean Comment;
  char *nextString;
  nodeCount = edgeCount = 0;
  T = newTree();
  nextString = (char *) malloc(numLeaves*INPUT_SIZE*sizeof(char));
  if (NULL == nextString)
    nextString = (char *) malloc(MAX_INPUT_SIZE*sizeof(char));
  Comment = FALSE;
  while(1 == fscanf(ifile,"%c",&c))
    {
      if('[' == c)
	Comment = TRUE;
      else if (']' == c)
	Comment = FALSE;
      else if (!(Comment))
	{
	  if(whiteSpace(c)) 
	    {
	      if (i > 0)
		nextString[i++] = ' ';
	    }
	  else  /*note that this else goes with if(whiteSpace(c))*/
	    nextString[i++] = c;
	  if (';' == c)
	    break;
	}
    }
  if ('(' != nextString[0])
    Exit ("Invalid input tree file format. Does not start with '('.");
  inputLength = i;
  for(i = 0; i < inputLength;i++)
    {
      if ('(' == nextString[i])
	parCount++;
      else if (')' == nextString[i])
	parCount--;
      if (parCount > 0)
	;
      else if (0 == parCount)
	{
	  i++;
//	  if(';' == nextString[i])
//	   ; //snprintf(rootLabel,MAX_LABEL_LENGTH,"URoot")
//	  else
//	    {
	      while(';' != nextString[i]) 
		if(!(whiteSpace(nextString[i++])))
		  rootLabel[j++] = nextString[i-1];  /*be careful here */
	      rootLabel[j] = '\0';
//	    }
	  i = inputLength;
	}
      else if (parCount < 0)
        Exit ("Invalid input tree file format. Too many right parentheses.");
    }
  centerNode = decodeNewickSubtree(nextString,T,&uCount);
  snprintf(centerNode->label,MAX_LABEL_LENGTH,"%s",rootLabel);
  T->root = centerNode;
  free(nextString);
  if (NULL != T->root->parentEdge)
    Warning ("Tree poorly rooted.");
  return(T);
}

void NewickPrintSubtree(tree *T, edge *e, FILE *ofile)
{
  if (NULL == e)
    Exit ("Newick Printing routine error.");
  if(!(leaf(e->head)))
    {
      fprintf(ofile,"(");
      NewickPrintSubtree(T,e->head->leftEdge,ofile);
      fprintf(ofile,",");
      NewickPrintSubtree(T,e->head->rightEdge,ofile);
      fprintf(ofile,")");
    }
  fprintf(ofile,"%s",e->head->label);
  fprintf(ofile,":%lf",e->distance);
  return;
}

void NewickPrintBinaryTree(tree *T, FILE *ofile)
{
  edge *e, *f;
  node *rootchild;
  e = T->root->leftEdge;
  rootchild = e->head;
  fprintf(ofile,"(");
  f = rootchild->leftEdge;
  if (NULL != f)
    {
      NewickPrintSubtree(T,f,ofile);
      fprintf(ofile,",");
    }
  f = rootchild->rightEdge;
  if (NULL != f)
    {
      NewickPrintSubtree(T,f,ofile);
      fprintf(ofile,",");
    }
  fprintf(ofile,"%s:%lf",T->root->label,e->distance);  
  fprintf(ofile,")");
  if (NULL != rootchild->label)
       fprintf(ofile,"%s",rootchild->label);
  fprintf(ofile,";\n");
  return;
}

void NewickPrintTrinaryTree(tree *T, FILE *ofile)
{
  edge *f;
  f = T->root->leftEdge;
  fprintf(ofile,"(");
  if (NULL != f)
    {
      NewickPrintSubtree(T,f,ofile);
      fprintf(ofile,",");
    }
  f = T->root->rightEdge;
  if (NULL != f)
    {
      NewickPrintSubtree(T,f,ofile);
      fprintf(ofile,",");
    }
  f = T->root->middleEdge;
  if (NULL != f)
    {
      NewickPrintSubtree(T,f,ofile);
      fprintf(ofile,")");
    }
  if (NULL != T->root->label)
       fprintf(ofile,"%s",T->root->label);
  fprintf(ofile,";\n");
  return;
}

void NewickPrintTree(tree *T, FILE *ofile)
{
  if (leaf(T->root))
    NewickPrintBinaryTree(T,ofile);
  else
    NewickPrintTrinaryTree(T,ofile);
  return;
}

void NewickPrintSubtreeStr(tree *T, edge *e, char *str)
{
  char *tmp;
  if (NULL == e)
    Exit ("Newick Printing routine error.");
  if(!(leaf(e->head)))
    {
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, "(", 1);
      NewickPrintSubtreeStr(T,e->head->leftEdge,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ",", 1);
      NewickPrintSubtreeStr(T,e->head->rightEdge,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ")", 1);
    }
  if (strlen (str) < MAX_INPUT_SIZE - strlen (e->head->label) -1)
    strncat (str, e->head->label, strlen (e->head->label));
  
  if (strlen (str) < MAX_INPUT_SIZE - 2)
    strncat (str, ":", 1);
  
  tmp = (char *)calloc(INPUT_SIZE, sizeof(char));
  if (strlen(tmp))
    strncpy(tmp, "", strlen(tmp));
  snprintf (tmp, INPUT_SIZE, "%lf", e->distance);
  if (strlen (str) < MAX_INPUT_SIZE - strlen (tmp) -1)
    strncat (str, tmp, strlen (tmp));
  free (tmp);
  return;
}

void NewickPrintBinaryTreeStr(tree *T, char *str)
{
  edge *e, *f;
  node *rootchild;
  char *tmp;
  e = T->root->leftEdge;
  rootchild = e->head;
  if (strlen (str) < MAX_INPUT_SIZE -2)
    strncat (str, "(", 1);
  f = rootchild->leftEdge;
  if (NULL != f)
    {
      NewickPrintSubtreeStr(T,f,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ",", 1);
    }
  f = rootchild->rightEdge;
  if (NULL != f)
    {
      NewickPrintSubtreeStr(T,f,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ",", 1);
    }
  if (strlen (str) < MAX_INPUT_SIZE - strlen (T->root->label) -1)
    strncat (str, T->root->label, strlen (T->root->label));
  
  if (strlen (str) < MAX_INPUT_SIZE - 2)
    strncat (str, ":", 1);
    
  tmp = (char *)calloc(INPUT_SIZE, sizeof(char));
  if (strlen(tmp))
    strncpy(tmp, "", strlen(tmp));
  snprintf (tmp, INPUT_SIZE, "%lf", e->distance);
  if (strlen (str) < MAX_INPUT_SIZE - strlen (tmp) -1)
    strncat (str, tmp, strlen (tmp));
  free (tmp);
  
  if (strlen (str) < MAX_INPUT_SIZE - 2)
    strncat (str, ")", 1);
  
  if (NULL != rootchild->label)
    if (strlen (str) < MAX_INPUT_SIZE - strlen (rootchild->label) -1)
      strncat (str, T->root->label, strlen (rootchild->label));

  if (strlen (str) < MAX_INPUT_SIZE - 3)
    strncat (str, ";\n", 2);
  return;
}

void NewickPrintTrinaryTreeStr(tree *T, char *str)
{
  edge *f;
  f = T->root->leftEdge;
  if (strlen (str) < MAX_INPUT_SIZE -2)
    strncat (str, "(", 1);
  if (NULL != f)
    {
      NewickPrintSubtreeStr(T,f,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ",", 1);
    }
  f = T->root->rightEdge;
  if (NULL != f)
    {
      NewickPrintSubtreeStr(T,f,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ",", 1);
    }
  f = T->root->middleEdge;
  if (NULL != f)
    {
      NewickPrintSubtreeStr(T,f,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ")", 1);
    }
  if (NULL != T->root->label)
    if (strlen (str) < MAX_INPUT_SIZE - strlen (T->root->label) -1)
      strncat (str, T->root->label, strlen (T->root->label));
  if (strlen (str) < MAX_INPUT_SIZE - 3)
    strncat (str, ";\n", 2);
  return;
}

void NewickPrintTreeStr(tree *T, char *str)
{
  if (leaf(T->root))
    NewickPrintBinaryTreeStr(T,str);
  else
    NewickPrintTrinaryTreeStr(T,str);
}

