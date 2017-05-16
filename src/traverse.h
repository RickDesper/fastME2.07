#ifndef TRAVERSE_H_
#define TRAVERSE_H_

#include "graph.h"

edge *findBottomLeft(edge *e);
edge *findBottomRight(edge *e);
edge *moveRight(edge *e);
edge *moveMiddle(edge *e);
edge *moveRightFirstMiddle(edge *e);
edge *moveLeft(edge *e);
edge *depthFirstTraverse(tree *T, edge *e);
edge *moveUpRight(edge *e);
edge *topFirstTraverse(tree *T, edge *e);
edge *depthRightFirstTraverse(tree *T, edge *e);

#endif /*TRAVERSE_H_*/

