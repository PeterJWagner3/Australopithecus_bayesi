#include "stepwise.h"

int *GetSimpleAdditionSequence (int *ntaxa) {
	int i;
	int *addseq;
	
	addseq=new int [*ntaxa+1];
	for (i=1; i<=*ntaxa; ++i) addseq[i]=i;	// sequential order
	return addseq;
}

int *GetStratAdditionSequence (DATA *data) {
	int i, j, k;
	int *addseq;
	
	k=1;
	addseq=new int [*data->GetNTaxa()+1];
	for (j=1;j<=*data->GetNStates(data->GetStratChr());++j) {
		for (i=1; i<=*data->GetNTaxa(); ++i) {
			if (*data->GetFInt(i)==j) {
				addseq[k]=i;
				++k;
			}
		}
	}
	return addseq;
}

void PutStratAdditionSequence (int *addseq, DATA *data) {
	int i, j, k;
	
	k=1;
	for (j=0;j<*data->GetNStates(data->GetStratChr());++j) {
		for (i=1; i<=*data->GetNTaxa(); ++i) {
			if (*data->GetFInt(i)==j) {
				addseq[k]=i;
				++k;
			}
		}
	}
}

void PutSimpleAdditionSequence(int *addseq, int *ntaxa) {
	for (int i=1; i<=*ntaxa; ++i) addseq[i]=i;	// sequential order
}

int *GetRandomAdditionSequence(int *ntaxa) {
	int i;
	bool *used=0;
	int *addseq=0;
	int selected_sp;
	float x;
	
	addseq=new int [*ntaxa+1];
	used=new bool [*ntaxa+1];
	for (i=0;i<=*ntaxa;++i) used[i]=0;

	used[0]=true;								// leave this in here - necessary for proper function
	for (i=1;i<=*ntaxa;++i) {
		addseq[i]=selected_sp=0;
		while (used[selected_sp] || selected_sp>*ntaxa) {	// select random taxa until we find one we haven't used
			x=0;
			while (x==0 || x==1) x=((float)rand()/(float)RAND_MAX)*(*ntaxa);
			selected_sp=(int)x+1;
		}
		addseq[i]=selected_sp;		// add it to the addition sequence
		used[selected_sp]=true;					// change the status to "used"
	}
	delete [] used;
	return addseq;
}

void PutRandomAdditionSequence(int *addseq, int *ntaxa) {
	int i;
	bool *used=0;
	int selected_sp;
	float x;
	
	used=new bool [*ntaxa+1];
	for (i=0;i<=*ntaxa;++i) used[i]=0;

	used[0]=true;								// leave this in here - necessary for proper function
	for (i=1;i<=*ntaxa;++i) {
		addseq[i]=selected_sp=0;
		while (used[selected_sp] || selected_sp>*ntaxa) {	// select random taxa until we find one we haven't used
			x=0;
			while (x==0 || x==1) x=((float)rand()/(float)RAND_MAX)*(*ntaxa);
			selected_sp=(int)x+1;
		}
		addseq[i]=selected_sp;		// add it to the addition sequence
		used[selected_sp]=true;					// change the status to "used"
	}
	delete [] used;
}

void GetStepwiseAdditionTree(NODE *tree, int *addseq, DATA *data)
{
	NODE *hold_node=0, *root=0, *old_root=0;
	ParsimonyClass parsimony;
	float best_length, current_length;
	int i, ntaxa, branch, best_branch=-1, nnodes, node, best_node=-1;
	bool taxon_added;
	
	ntaxa=*data->GetNTaxa();
	tree[ntaxa+1].left=&tree[addseq[1]];			// first node's left descendant is first taxon in addition sequence
	tree[ntaxa+1].right=&tree[addseq[2]];			// first node's right descendant is next taxon in addition sequence
	nnodes=1;														// start with just one node
	root=&tree[ntaxa+1];
		
	for (i=3;i<=ntaxa;++i) {									// for all remaining taxa in addition sequence (tree begins with two taxa)
		best_length=10e10;
		for (node=0;node<=nnodes;++node) {						// for all nodes
			for (branch=0;branch<=1;++branch) {					// for each branch from that node
				taxon_added=false;
				tree[nnodes+ntaxa+1].left=&tree[addseq[i]];		// new node's left is the newly added taxon
				if (node==0) {									// try new taxon at the base of this tree
					tree[nnodes+ntaxa+1].right=root;			// new node's right is the previous root node
					tree[nnodes+ntaxa+1].anc=0;					// new node's ancestor is the ROOT (0)
					root->anc=&tree[nnodes+ntaxa+1];			// previous node's new anc is now the new node
					old_root=root;								// keep track of the old root node
					root=&tree[nnodes+ntaxa+1];					// make the new node the root_node
					taxon_added=1;
				}
				else if (node>0 && !branch) {					// try new node on left branch
					hold_node=tree[node+ntaxa].left;
					if (hold_node->isAnc) hold_node->isAnc=0;			// if branch is a fixed ancestor, unfix it
					tree[nnodes+ntaxa+1].right=tree[node+ntaxa].left;	// new node's right is the previous node's left descendant
					tree[node+ntaxa].left=&tree[nnodes+ntaxa+1];		// previous node's new left is now the new node
					taxon_added=1;
				}
				else {		 									// try new node on right branch
					hold_node=tree[node+ntaxa].right;
					if (hold_node->isAnc) hold_node->isAnc=0;			// if branch is a fixed ancestor, unfix it
					tree[nnodes+ntaxa+1].right=tree[node+ntaxa].right;	// new node's right is the previous node's right descendant
					tree[node+ntaxa].right=&tree[nnodes+ntaxa+1];		// previous node's new right is now the new node
					taxon_added=1;
				}
				GetAncs(root);

			// Check Treelength of Current Position
//				if (data->hasStratChr) parsimony.DateTree(root, data);
				current_length=parsimony.TotalDebt(root, data);
				if (current_length<best_length)	{		// if the current placement is optimal
					best_node=node;						// record the optimal node
					best_branch=branch;					// record the optimal branch on that node
					best_length=current_length;			// record best current length
				}

			// Remove the taxon
				if (taxon_added) {
					tree[addseq[i]].anc=0;		// removal of taxon removes it's anc
					tree[addseq[i]].sister=0;	// removal of taxon removes it's sister
					if (node==0) {
						root=old_root;						// return the root_node to it's previous owner
						root->anc=0;						// remove new taxon from the base of this tree
						root->sister=0;						// remove old root's sister
					}
					else if (node>0 && branch==0) tree[node+ntaxa].left=tree[nnodes+ntaxa+1].right;	// restore old node membership
					else tree[node+ntaxa].right=tree[nnodes+ntaxa+1].right;							// restore old node membership
				}
				if (node==0) break;										// prevent for loop from doing this twice
			} // end for possible branch
		} // end for each node
		
	  // Once found, place the taxon in at the optimal node and branch
		++nnodes;
		tree[nnodes+ntaxa].left=&tree[addseq[i]];		// add new taxon to new node
		if (best_node==0) {
			tree[nnodes+ntaxa].right=root;				// new node's new right descendant is the rest of the tree
			tree[nnodes+ntaxa].anc=0;					// new node's ancestor is the ROOT (0)
			tree[nnodes+ntaxa].sister=0;				// new node has no sister
			root->anc=&tree[nnodes+ntaxa];				// old root node's ancestor is the  new node
			root=&tree[nnodes+ntaxa];					// new node is the root node
		}
		else {
			if (best_branch==0) {
				tree[nnodes+ntaxa].right=tree[best_node+ntaxa].left;// transfer optimal node's old left descendant to new node
				tree[best_node+ntaxa].left=&tree[nnodes+ntaxa];		// optimal node's new left descendant is the new node
			}
			else {
				tree[nnodes+ntaxa].right=tree[best_node+ntaxa].right;// transfer optimal node's old left descendant to new node
				tree[best_node+ntaxa].right=&tree[nnodes+ntaxa];	// optimal node's new left descendant is the new node
			}
		}
		GetAncs(root);
	}
	LadderizeRight(tree, data->GetNTaxa());
}

void PutDefaultLadder(NODE *tree, int ntaxa) {
	int current_node=ntaxa+1;
	tree[current_node].left=&tree[ntaxa-1];
	tree[current_node].right=&tree[ntaxa];
	++current_node;
	for (int i=ntaxa-2;i>=1;--i) {
		tree[current_node].left=&tree[i];
		tree[current_node].right=&tree[current_node-1];
		++current_node;
	}
	NODE *root=&tree[(2*ntaxa)-1];
	GetAncs(root);
	ReorderNodes(tree, &ntaxa);
}
