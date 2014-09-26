#ifndef CODE4CHIEF_H#define CODE4CHIEF_H#include "node.h"#include "data_field.h"#include "tree_fns.h"using namespace std;void matrixToData(long **matrix, DATA *data);void typesToData(int *types, DATA *data);void treeToMatrix(NODE *tree, long **tree_matrix, int ntaxa);void treeToMatrixMachine(NODE *this_node, long  **tree_matrix, int ntaxa);void matrixToTree(long **tree_matrix, NODE *tree, int ntaxa);#endif