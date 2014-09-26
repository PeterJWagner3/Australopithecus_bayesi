#ifndef STEPWISE_H
#define STEPWISE_H

	#include <stdlib.h>
	#include "data_field.h"
	#include "node.h"
	#include "parsimony_class_field.h"
	#include "tree_fns.h"

	int *GetSimpleAdditionSequence(int *);
	void PutSimpleAdditionSequence(int *, int *);
	int *GetStratAdditionSequence (DATA * );
	void PutStratAdditionSequence (int * , DATA * );
	int *GetRandomAdditionSequence(int *);
	void PutRandomAdditionSequence(int *, int *);
	NODE *GetStepwiseAdditionTree(int *, DATA *);
	void GetStepwiseAdditionTree(NODE *, int *, DATA *);
	void PutDefaultLadder(NODE *tree, int );

#endif

