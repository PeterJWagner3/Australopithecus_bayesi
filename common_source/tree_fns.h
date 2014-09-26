#ifndef TREE_FNS_H
#define TREE_FNS_H

	#include <iostream>
	#include <fstream>
	#include "data_field.h"
	#include "node.h"
	#include "memory_koz.h"
	
	using namespace std;

	NODE *GetRoot(NODE *const , int *);
	NODE *GetRootPartialTree(NODE *const );
	void ReorderNodes(NODE *const, int *);
	void ReorderMachine(NODE *const , NODE *const , int *, int *, int &);
	void RerootWithoutReorder(NODE * const , NODE * const , int *);
	void RerootClade(NODE * const , NODE * const , NODE * const , int * );
	void Reroot(NODE * const , NODE * const , int *);
	void RerootMachine(NODE *, NODE *);
	void GetAncs(NODE *);
	void MakeStratBL(NODE * , DATA * );
	void CutCladeFromTree(NODE * , int * , int * );
	
	NODE *GetInitializedTree(int );
	NODE ** AllocateTreeList(int , int );
//	void ** AllocateTreeList(NODE ** const , int * , int * );
	void FreeTreeList(NODE **, int *);
	void PrintTree(NODE *, int *);
	void CopyClade(NODE *, NODE *);
	void CopyTree(NODE *const , NODE *const , int *);

	int CountTopologies(NODE **, int *, int *);

	void ReadCladeFromFile (ifstream * const , NODE * const , int * const );
	void ReadCladeMachine (ifstream * const , NODE * const , int , int &, int * const );

	void UnfixTree(NODE *, int *);
	void UnfixMachine(NODE *);

	void AddTreeToCube(NODE *, bool ** const, int *);
	void AddToCubeMachine(NODE *, bool ** const, int *);
	bool CompareCubes(bool ** const , bool ** const, int *);
	bool CompareTopologies(NODE * , NODE * , int * );
	void PrintCube(bool **, int * );
	bool CompareAncs(NODE *, NODE *, int * );
	void MakeRichness(NODE * );
	void LadderizeRight(NODE * , int * );
	void LadderizeRightMachine(NODE * );
	
	void MakeReconstructedRichnessInt(int ** , DATA *data, NODE *original_tree);
	void ReconstructedRichnessMachine(DATA *data, NODE * a, int **);
	void MakeReconstructedDisparityInt(float ** , DATA *data, NODE *tree);
	void MakeObservedMembership(DATA *data, int *members, int intv);
	void MakeReconstructedMembership(DATA *data, NODE *tree, int *members, int intv);
	float MakeDisparity(DATA *data, int *members);
	void MakeReconstructedChrEvInt(float **chr_ev, DATA *data, NODE *tree);
	void MakeReconstructedChrEvIntLA(float **chr_ev, DATA *data, NODE *tree);

#endif