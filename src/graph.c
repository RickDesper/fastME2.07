#include "graph.h"

boolean leaf(node *v)
{
  int count = 0;
  if (NULL != v->parentEdge)
    count++;
  if (NULL != v->leftEdge)
    count++;
  if (NULL != v->rightEdge)
    count++;
  if (NULL != v->middleEdge)
    count++;
  if (count > 1)
    return(FALSE);
  return(TRUE);
}

set *addToSet(node *v, set *X)
{
  if (NULL == X)
    {
      X = (set *) malloc(sizeof(set));
      X->firstNode = v;
      X->secondNode = NULL;
    }
  else if (NULL == X->firstNode) 
    X->firstNode = v;
  else
    X->secondNode = addToSet(v,X->secondNode);
  return(X);
}

set *copySet(set *X)
{
  set *ret = NULL;
  if (NULL != X) {
    ret = (set *) malloc(sizeof(set));
    ret->firstNode = copyNode(X->firstNode);
    ret->secondNode = copySet (X->secondNode);
  }
  return ret;
}

node *makeNode(char *label, int index)
{
  node *newNode;  /*points to new node added to the graph*/
  newNode = (node *) malloc(sizeof(node));
  strncpy(newNode->label, label, MAX_LABEL_LENGTH);
  newNode->index = index;
  newNode->index2 = -1;
  newNode->parentEdge = NULL;
  newNode->leftEdge = NULL;
  newNode->middleEdge = NULL;
  newNode->rightEdge = NULL;
  /*all fields have been initialized*/
  return(newNode);
}

edge *makeEdge(char *label, node *tail, node *head, double weight)
{
  edge *newEdge;
  newEdge = (edge *) malloc(sizeof(edge));
  strncpy(newEdge->label, label, MAX_LABEL_LENGTH);
  newEdge->tail = tail;
  newEdge->head = head;
  newEdge->distance = weight;
  newEdge->totalweight = 0.0;
  return(newEdge);
}

tree *newTree()
{
  tree *T;
  T = (tree *) malloc(sizeof(tree));
  T->root = NULL;
  T->size = 0;
  T->weight = -1;
  return(T);
}

void freeSubTree(edge *e)
     /*frees subtree below e->head, recursively, then frees e->head and e*/
{
  node *v;
  edge *e1, *e2;
  v = e->head;
  e1 = v->leftEdge;
  if (NULL != e1)
    freeSubTree(e1);
  e2 = v->rightEdge;
  if (NULL != e2)
    freeSubTree(e2);
  if (NULL != v)
    free(v);
  /*e->tail = NULL;
  e->head = NULL;*/
  if (NULL != e)
    free(e);
  return;
}

void freeTree(tree *T)
{
  node *v;
  v = T->root;
  if (NULL != v->leftEdge)
    freeSubTree(v->leftEdge);
  if (NULL != T->root)
    free(T->root);
  if (NULL != T)
    free(T);
  return;
}

void freeSet(set *S)
{
  if (NULL != S)
    {
      freeSet(S->secondNode);
      if (NULL != S->firstNode)
        free(S->firstNode);
    }
  if (NULL != S)
    free(S);
  return;
}

/*copyNode returns a copy of v which has all of the fields identical to those
of v, except the node pointer fields*/
node *copyNode(node *v)
{
  node *w;
  w = makeNode(v->label,v->index);
  w->index2 = v->index2;
  return(w);
}

/*copyEdge calls makeEdge to make a copy of a given edge */
/*does not copy all fields*/
edge *copyEdge (edge *e)
{
  edge *newEdge;
  newEdge = makeEdge(e->label,e->tail,e->head,e->distance);
  newEdge->topsize = e->topsize;
  newEdge->bottomsize = e->bottomsize;
  e->totalweight = 1.0;
  newEdge->totalweight = -1.0;
  return(newEdge);
}

/*detrifurcate takes the (possibly trifurcated) input tree
  and reroots the tree to a leaf*/
/*assumes tree is only trifurcated at root*/
tree *detrifurcate(tree *T)
{
  node *v, *w;
  edge *e, *f;
  v = T->root;
  w = NULL;
  if(leaf(v))
    return(T);
  if (NULL != v->parentEdge)
    {
      fprintf(stderr,"\n Error: root %s is poorly rooted.\n",v->label);
      exit(EXIT_FAILURE);
    }
  for(e = v->middleEdge, v->middleEdge = NULL; NULL != e; e = f )
    {
      w = e->head;
      v = e->tail;
      e->tail = w;
      e->head = v;
      f = w->leftEdge;
      v->parentEdge = e;
      w->leftEdge = e;
      w->parentEdge = NULL;      
    }
  T->root = w;
  return(T);
}

edge *siblingEdge(edge *e)
{
  if(e == e->tail->leftEdge)
    return(e->tail->rightEdge);
  else
    return(e->tail->leftEdge);
}

void updateSizes(edge *e, int direction)
{
  edge *f;
  switch(direction)
    {
    case UP:
      f = e->head->leftEdge;
      if (NULL != f)
	updateSizes(f,UP);
      f = e->head->rightEdge;
      if (NULL != f)
	updateSizes(f,UP);
      e->topsize++;
      break;
    case DOWN:
      f = siblingEdge(e);
      if (NULL != f)
	updateSizes(f,UP);
      f = e->tail->parentEdge;
      if (NULL != f)
	updateSizes(f,DOWN);
      e->bottomsize++;
      break;
    }
  return;
}      

node *copySubtree(node *v)
{
	node *newNode;
	newNode = copyNode(v);
	if (NULL != v->leftEdge)
	{
		newNode->leftEdge = copyEdge(v->leftEdge);
		newNode->leftEdge->tail = newNode;
		snprintf(newNode->leftEdge->label,MAX_LABEL_LENGTH,"%s",v->leftEdge->label);
		newNode->leftEdge->head = copySubtree(v->leftEdge->head);
		newNode->leftEdge->head->parentEdge = newNode->leftEdge;
		if (newNode->leftEdge->totalweight > 0.0)
			printf("oops\n");
	}
	if (NULL != v->rightEdge)
	{
		newNode->rightEdge = copyEdge(v->rightEdge);
		newNode->rightEdge->tail = newNode;
		snprintf(newNode->rightEdge->label,MAX_LABEL_LENGTH,"%s",v->rightEdge->label);
		newNode->rightEdge->head = copySubtree(v->rightEdge->head);
		newNode->rightEdge->head->parentEdge = newNode->rightEdge;
		if (newNode->rightEdge->totalweight > 0.0)
			printf("oops\n");
	}

	return(newNode);
}

tree *copyTree(tree *T)
{
	tree *Tnew;
	node *n1, *n2, *n3;
	edge *e1, *e2;
	n1 = copyNode(T->root);
	Tnew = newTree();
	Tnew->root = n1;
	if (NULL != T->root->leftEdge)
	{
			e1 = copyEdge(T->root->leftEdge);
			snprintf(e1->label,MAX_LABEL_LENGTH,"%s",T->root->leftEdge->label);
			n1->leftEdge = e1;
			n2 = copySubtree(e1->head);
			e1->head = n2;
			e1->tail = n1;
			n2->parentEdge = e1;
			if (n2->parentEdge->totalweight > 0.0)
				printf("oops\n");
	}
	if (NULL != T->root->rightEdge)
	{
		e2 = copyEdge(T->root->rightEdge);
		snprintf(e2->label,MAX_LABEL_LENGTH,"%s",T->root->rightEdge->label);
		n1->rightEdge = e2;
		n3 = copySubtree(e2->head);
		e2->tail = n1;
		e2->head = n3;
		n3->parentEdge = e2;
	}
	Tnew->size = T->size;
	Tnew->weight = T->weight;
	return(Tnew);
}

void weighTree(tree *T)
{
  edge *e;
  T->weight = 0;
  for(e = depthFirstTraverse(T,NULL);NULL!=e;e=depthFirstTraverse(T,e))
    T->weight += e->distance;
  return;
}

edge *findEdge(tree *T, edge *e)
{
	edge *f;
	for(f = depthFirstTraverse(T,NULL);NULL!=f;f=depthFirstTraverse(T,f))
	{
		if (0 == strcmp(e->label,f->label))
			return(f);
	}
	fprintf(stderr,"\n Error: Cannot find edge %s with tail %s and head %s\n",e->label,e->tail->label,e->head->label);
	exit(EXIT_FAILURE);
}

node *indexedNode(tree *T, int i)
{
  edge *e;
  for(e = depthFirstTraverse(T,NULL);NULL!=e;e=depthFirstTraverse(T,e))
    if (i == e->head->index)
      return(e->head);
  if (i == T->root->index)
    return(T->root);
  return(NULL);
}

edge *indexedEdge(tree *T, int i)
{
  edge *e;
  for(e = depthFirstTraverse(T,NULL);NULL!=e;e=depthFirstTraverse(T,e))
    if (i == e->head->index)
      return(e);
  return(NULL);
}

