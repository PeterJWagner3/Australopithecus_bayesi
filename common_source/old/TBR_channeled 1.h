#ifndef TBR_CHANNELED_H
#define TBR_CHANNELED_H

#include <iostream>
#include "parsimony_class_field.h"
#include "memory_koz.h"
#include "tree_fns.h"
#include "nexus_tree_write.h"
#include "stepwise.h"

class HeuristicSearchClass {
	public:
		HeuristicSearchClass(int * const ); 
		~HeuristicSearchClass();
		void HeuristicSearch(DATA * const, NODE ** const, int &);
//		NODE *Bootstrap(DATA * const, int *, int *);
		bool OutgroupSpecified;
		void SetOutgroup(int , bool );
		bool GetOutgroup(int ) const;
		void SetMaxTrees(int );
		int *GetMaxTrees();
		void SetChannel(int );
		int *GetChannel();
		void SetSearchReps(int );
		void SetAncSearchReps(int );
		void SetAncSearchType(bool );
		void SetMaxSavingsMethod(bool );
		void SetStartTrees(NODE ** , int * );
		
	private:
		bool AncSearchType;			// 0=exhaustive, 1=heuristic
		int AncSearchReps; 
		bool MaxSavingsMethod;
		NODE **tree_list;			// single rep optimal trees
		NODE *tree;					// dummy tree for Regraft
		NODE *tree2;				// dummy tree for ReRootPrunedNode
		NODE *tree3;				// dummy tree for GetPrunedNode
		NODE *tree4;				// dummy tree for Rooting within Regraft
		NODE **starttrees;			// user  defined  trees to use instead of random addition sequence
		 
		int *add_seq;				// addition sequence for heuristic anc search
		bool ***tree_cube;			// int matrix representation of topologies including all globally and this_rep optimal trees
		bool **this_cube;			// int matrix representation of this topology
		int *ntaxa;
		int root_taxon;
		int reps;
		int current_rep;
		int swapped_trees;			// number of trees swapped in the current rep
		int trees_this_rep;			// number of optimal trees found in the current rep
		float best_debt;			// best overall debt from all reps
		float best_debt_cladogram;	// best overall debt from all reps for the optimal cladogram
		float best_debt_this_rep;	// best debt during a single rep
		float best_debt_this_rep_cladogram;	// best debt during a single rep
//		int rearrangements;			// total rearrangements tried
		int MaxTrees;
		int ChannelSize;
		int nstarttrees;

		int nislands;
		class IslandData *island;
		ParsimonyClass parsimony;
		
		bool abort_swap;			// toggle to abort swapping on the "current" tree
		bool abort_rep;				// toggle to abort the entire rep
		bool added_topology;
		bool random; 				// toggle either runs with random addition sequences, or static (for debugging)

		bool *outgroup;	// 2*ntaxa long array of outgroup membership (0=ingroup, 1=outgroup)
		bool SingleOutgroup;
		void AllocateOutgroup();
		void MakeOutgroupFromTree(NODE * );
		void initializeOutgroup(DATA *);

		void SingleHeuristicSearch(DATA * const, int &);										// performs a single replicate heuristic search
		void GetPrunedNode (NODE * const, DATA * const , int &);									// recursively finds all possible subtrees to break off
		void ReRootPrunedNode(NODE * const , NODE * const , DATA * const data, int &);
		void GetTargetNode (NODE * const, NODE * const, NODE * const , DATA * const, int &);	// recursively finds all possible nodes to which to reattach subtree
		void RegraftSubtree (NODE * const, int * const, int * const, DATA * const , int &);		// regrafts subtree to target node and determines the optimality of the new tree

		void beginAncSearch (DATA * , int &);
		void FullAncSearch(NODE * const, NODE * const, int , DATA * const, int &);	// finds all optimal ancestorizations
		void HeuristicAncSearch(NODE * const, NODE * const, int , DATA * const, int &);	// finds all optimal ancestorizations
		void FinishAncSearch(NODE * const, NODE * const , DATA * const, int &);	// finds all optimal ancestorizations

		float MaxStratSavings(NODE *, DATA *, int *, bool);
		void DoMinStratDebt(NODE *, NODE *, int , DATA * const , float &);
		void DoMinTotalDebt(NODE *, NODE *, int , DATA * const , float &, int * const);
		void QuickMorphDebt(NODE * const , DATA * const , float *);
		void QuickStratSavings(NODE * const , DATA * const , float *);

		void PrintIslandData();
};
	
class IslandData {
	public:
		int ntrees;
		int ntopologies;
		int first_tree;
		int last_tree;
		int first_rep;
		int times_hit;
		float debt;
};

#endif

