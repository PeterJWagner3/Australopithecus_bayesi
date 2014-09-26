#include "tree_fns.h"

using namespace std;

NODE *GetRoot(NODE *const tree, int *ntaxa)
{
	NODE *root=0;
	
	root=&tree[*ntaxa+1];
	while (root->anc==&tree[0]) root=&tree[root->label+1];
	while (root->anc!=0) root=root->anc;
	return root;
}

NODE *GetRootPartialTree(NODE *const tree)
{
	NODE *root=0;
	
	root=&tree[1];
	while (root->anc==&tree[0]) root=&tree[root->label+1];
	while (root->anc!=0) root=root->anc;
	return root;
}

void ReorderNodes(NODE *const start_tree, int *ntaxa) {
	NODE *ordered_tree=0, *root=0;
	int o_nodes, root_node;
	
	ordered_tree=GetInitializedTree(*ntaxa);
	o_nodes=*ntaxa+1;
	root=GetRoot(start_tree, ntaxa);
	root_node=root->label;
	ReorderMachine(start_tree, ordered_tree, ntaxa, &root_node, o_nodes);		// copy it to a new tree with ordered nodes
	root=GetRoot(ordered_tree, ntaxa);
	GetAncs(root);
	CopyTree(ordered_tree, start_tree, ntaxa);
	delete [] ordered_tree;
}

void ReorderMachine(NODE *const orig_tree, NODE * const new_tree, int *ntaxa, int *orig_node, int &s_nodes)
{
	NODE *hold_node=0, *hold_node2=0;
	int current_node;
	int gobots;

	current_node=s_nodes;
	if (current_node>=2**ntaxa) {
		cout << "\n\n***\tThis node is messed up\n";
		cin >>gobots;
	}
	if (s_nodes==*ntaxa+1) new_tree[current_node].anc=0;		// if this is the root node of the new tree, the ancestor is the ROOT (0)

	hold_node=orig_tree[*orig_node].left;
	if (!hold_node->tip) {
		new_tree[current_node].left=&new_tree[s_nodes+1];		// new tree's left descendant is a new node
		new_tree[s_nodes+1].anc=&new_tree[current_node];		// that new node's ancestor is the current node
		++s_nodes;												// increment the number of nodes
		ReorderMachine(orig_tree, new_tree, ntaxa, &hold_node->label, s_nodes);
	}
	else {														// the left descendant is an OTU
		
		new_tree[current_node].left=&new_tree[hold_node->label];
		new_tree[hold_node->label].anc=&new_tree[current_node];
		hold_node2=new_tree[current_node].left;
		hold_node2->tip=1;
		hold_node2->isAnc=hold_node->isAnc;
		}

	hold_node=orig_tree[*orig_node].right;
	if (!hold_node->tip) {
		new_tree[current_node].right=&new_tree[s_nodes+1];	// source tree's right descendant is a new node
		new_tree[s_nodes+1].anc=&new_tree[current_node];		// that new node's ancestor is the current node
		++s_nodes;										// increment the number of nodes
		ReorderMachine(orig_tree, new_tree, ntaxa, &hold_node->label, s_nodes);
	}
	else {
		new_tree[current_node].right=&new_tree[hold_node->label];	// if the right descendant is an OTU
		new_tree[hold_node->label].anc=&new_tree[current_node];
		hold_node2=new_tree[current_node].right;
		hold_node2->tip=1;
		hold_node2->isAnc=hold_node->isAnc;
	}
	
	hold_node=new_tree[current_node].left;
	hold_node2=new_tree[current_node].right;
	hold_node->sister=hold_node2;
	hold_node2->sister=hold_node;
	if (hold_node->label > hold_node2->label) {
		new_tree[current_node].left=hold_node2;
		new_tree[current_node].right=hold_node;
	}
}

void RerootWithoutReorder(NODE *const tree, NODE *const reroot_node, int *ntaxa)
{
	NODE *hold_node=0, *hold_node2=0;
	NODE *left_desc=0, *right_desc=0;
	NODE *root;
	int i, j;
	bool go;
	int root_node;
//	int gobots;	

	root=GetRoot(tree, ntaxa);
	GetAncs(root);
	left_desc=root->left;
	right_desc=root->right;
	if (reroot_node==root || reroot_node==left_desc || reroot_node==right_desc) return;

//	cout << ("Tree Enters (ReRoot on %d)\n", reroot_node->label);
//	PrintTree(tree, ntaxa);
	go=1;
	for (i=0;i<=1;++i)
		{
		root=GetRoot(tree, ntaxa);
		if (!i) hold_node=root->left;
		else hold_node=root->right;
		if (hold_node->tip && hold_node->isAnc)	//	if one taxon is fixed
			{
			for (j=0;j<=1;++j)
				{
				if (!i) hold_node2=root->right;			// hold_node2 is hold_node's sister
				else hold_node2=root->left;
				if (!j) hold_node2=hold_node2->left;
				else hold_node2=hold_node2->right;
				if (hold_node2->tip && hold_node2->isAnc)
					{
//					cout << ("\n\n***\tRerooting this tree would have resulted in two fixed ancs on a single node, so rerooting was not performed\n");
//					go=0;

//					cout << ("***\tRerooting this tree would have resulted in two fixed ancs on a single node, so unanc'ing both\n");
					hold_node2->Unfix();
					hold_node->Unfix();
					go=1;
					
					}
				}
			}
		}

	if (go)
		{		
		root_node=root->label;
//		while (tree[root_node].anc!=0 && root<2*(*ntaxa)) ++root;	// find the current root of the input tree
//		hold_node=reroot_node->anc;
//		if (reroot_node->label==root || hold_node->label==root) return;	// abort if the requested root is the current root
//		else 
		
		RerootMachine(reroot_node->anc, reroot_node);

		tree[root_node].left=reroot_node;							// new root's left member reroot_node
		tree[root_node].right=reroot_node->anc;						// new root's right member reroot_node's anc
		tree[root_node].anc=0;										// new root's anc is ROOT
		GetAncs(&tree[root_node]);									// copy in new ancs and sisters
		LadderizeRight(tree, ntaxa);
//		cout << ("Tree Before Reorder (ReRoot)\n");
//		PrintTree(tree, ntaxa);
//		ReorderNodes(tree, ntaxa);							// reorder nodes in newly rerooted treee (THIS DOESN'T "STICK" UPON RETURN TO MAIN)
//		GetAncs(&tree[*ntaxa+1]);							// GetAncs should be performed in ReorderNodes
//		cout << ("Tree After Reorder (ReRoot)\n");
//		PrintTree(tree, ntaxa);
//		cin >>gobots;
		}
}

void RerootClade(NODE *const tree, NODE * const clade, NODE * const reroot_node, int *ntaxa)
{
	NODE *hold_node=0, *hold_node2=0;
	NODE *left_desc=0, *right_desc=0;
	NODE *root=0, *last_anc=0;
	int i, j;
	bool go;
//	int root_node;
//	int gobots;	

	if (clade->tip) return;

//	cout << ("ReRootClade:\tTree Enters (ReRoot clade %d on %d)\n", clade->label, reroot_node->label);
//	PrintTree(tree, ntaxa);

	left_desc=clade->left;
	right_desc=clade->right;
	if (reroot_node==left_desc || reroot_node==right_desc) return;
	last_anc=clade->anc;
	clade->anc=0;

	go=1;
	for (i=0;i<=1;++i)
		{
		if (!i) hold_node=clade->left;
		else hold_node=clade->right;
		if (hold_node->tip && hold_node->isAnc)	//	if one taxon is fixed
			{
			for (j=0;j<=1;++j)
				{
				if (!i) hold_node2=clade->right;			// hold_node2 is hold_node's sister
				else hold_node2=clade->left;
				if (!j) hold_node2=hold_node2->left;
				else hold_node2=hold_node2->right;
				if (hold_node2->tip && hold_node2->isAnc)
					{
//					cout << ("\n\n***\tRerooting this tree would have resulted in two fixed ancs on a single node, so rerooting was not performed\n");
//					go=0;

//					cout << ("***\tRerooting this tree would have resulted in two fixed ancs on a single node, so unanc'ing both\n");
					hold_node2->Unfix();
					hold_node->Unfix();
					go=1;
					}
				}
			}
		}

	if (go)
		{		
		RerootMachine(reroot_node->anc, reroot_node);

		clade->left=reroot_node;							// new root's left member reroot_node
		clade->right=reroot_node->anc;						// new root's right member reroot_node's anc
		clade->anc=last_anc;								// new root's anc is ROOT

//		cout << ("Tree Before GetAncs (ReRootClade):\n");
//		PrintTree(tree, ntaxa);

		for (j=*ntaxa+1;j<2**ntaxa;++j)
			{
			if (!tree[j].anc) root=&tree[j];
			}
		GetAncs(root);										// copy in new ancs and sisters
//		cout << ("Tree Leaves (ReRootClade):\n");
//		PrintTree(tree, ntaxa);
//		ReorderNodes(tree, ntaxa);							// reorder nodes in newly rerooted treee (THIS DOESN'T "STICK" UPON RETURN TO MAIN)
//		GetAncs(&tree[*ntaxa+1]);							// GetAncs should be performed in ReorderNodes
//		cout << ("Tree After Reorder (ReRoot)\n");
//		PrintTree(tree, ntaxa);
//		cin >>gobots;
		}
}

void Reroot(NODE *const tree, NODE * const reroot_node, int *ntaxa) {
	NODE *hold_node=0, *hold_node2=0, *left_desc=0, *right_desc=0, *root;
	int i, j, root_node;
	bool go;
//	int gobots;	

	root=GetRoot(tree, ntaxa);
	GetAncs(root);
	left_desc=root->left;
	right_desc=root->right;
	if (reroot_node==root || reroot_node==left_desc || reroot_node==right_desc) {
		LadderizeRight(tree, ntaxa);
		ReorderNodes(tree, ntaxa);
		return;
	}

//	cout << "Tree Enters (ReRoot on "<< reroot_node->label << endl;
//	PrintTree(tree, ntaxa);

	go=1;
	for (i=0;i<=1;++i) {
		root=GetRoot(tree, ntaxa);
		if (!i) hold_node=root->left;
		else hold_node=root->right;
		if (hold_node->tip && hold_node->isAnc)	//	if one taxon is fixed
			{
			for (j=0;j<=1;++j)
				{
				if (!i) hold_node2=root->right;
				else hold_node2=root->left;
				if (!j) hold_node2=hold_node2->left;
				else hold_node2=hold_node2->right;
				if (hold_node2->tip && hold_node2->isAnc)
					{
//					cout << ("\n\n***\tRerooting this tree would have resulted in two fixed ancs on a single node, so rerooting was not performed\n");
//					go=0;

//					cout << ("***\tRerooting this tree would have resulted in two fixed ancs on a single node, so unanc'ing both\n");
					hold_node2->Unfix();
					hold_node->Unfix();
					go=1;
					}
				}
			}
		}

	if (go) {		
		root_node=root->label;
//		while (tree[root_node].anc!=0 && root<2*(*ntaxa)) ++root;	// find the current root of the input tree
//		hold_node=reroot_node->anc;
//		if (reroot_node->label==root || hold_node->label==root) return;	// abort if the requested root is the current root
//		else 
		
		RerootMachine(reroot_node->anc, reroot_node);
		tree[root_node].left=reroot_node;							// new root's left member reroot_node
		tree[root_node].right=reroot_node->anc;						// new root's right member reroot_node's anc
		tree[root_node].anc=0;										// new root's anc is ROOT
		GetAncs(&tree[root_node]);									// copy in new ancs and sisters
//		cout << ("Tree Before Ladderize (ReRoot)\n");
//		PrintTree(tree, ntaxa);
		LadderizeRight(tree, ntaxa);
//	cin>>gobots;
//		cout << ("Tree Before Reorder (ReRoot)\n");
//		PrintTree(tree, ntaxa);
		ReorderNodes(tree, ntaxa);							// reorder nodes in newly rerooted treee (THIS DOESN'T "STICK" UPON RETURN TO MAIN)
//	cin>>gobots;
		GetAncs(&tree[*ntaxa+1]);							// GetAncs should be performed in ReorderNodes
//		cout << ("Tree After Reorder (ReRoot)\n");
//		PrintTree(tree, ntaxa);
//	cin >>gobots;
	}
}

void RerootMachine(NODE *current_node, NODE *previous_node) {
	NODE *anc_node=0, *hold_node=0;
	int gobots;

	anc_node=current_node->anc;								// find ancestor of ancestor (checking to see if current node's anc is the root of the input tree)
	if (current_node->label==anc_node->label) {
		cout << "\n\n***\tThis node is messed up\n";
		cin >>gobots;
	}
	if (anc_node->anc!=0) {									// if the current node's anc isn't the root of the input tree
		RerootMachine (current_node->anc, current_node);	// continue traversing the tree
		if (current_node->left==previous_node) hold_node=current_node->right;	// if we came from the left, keep track of the right
		else hold_node=current_node->left;					// if from right, keep track of left
		current_node->left=hold_node;			// "new" current node's left is the node from the previous two lines
		current_node->right=current_node->anc;	// "new" current node's right is its former ancestor
	}
	else {													// if current node's anc IS the root of the input tree, we need to cut out the former root
		if (current_node->left==previous_node) {				// if we came from the left
			if (anc_node->left==current_node) current_node->left=anc_node->right;
			else current_node->left=anc_node->left;
		}
		else {												// if we came from the right
			if (anc_node->left==current_node) current_node->right=anc_node->right;
			else current_node->right=anc_node->left;
		}
	}
	current_node->anc=previous_node;			// "new" current node's anc is it's former "upstream" node (from whence we came)
}

void GetAncs(NODE *node) {
	NODE *left_desc=0, *right_desc=0;
//	NODE *left_desc2=0, *right_desc2=0;
//	int gobots;
	
	left_desc=node->left;
	right_desc=node->right;
/*	left_desc2=left_desc->left;
	right_desc2=left_desc->right;
	if (left_desc2==node || right_desc2==node) {
		cout << "\n\n\t***\tThe descendant of the left node is also the anc";
		cin >>gobots;
	}
	left_desc2=right_desc->left;
	right_desc2=right_desc->right;
	if (left_desc2==node || right_desc2==node) {
		cout << "\n\n\t***\tThe descendant of the right node is also the anc";
		cin >>gobots;
		}
*/
	if (!left_desc->tip) GetAncs(left_desc);
	if (!right_desc->tip) GetAncs(right_desc);
	left_desc->anc=node;
	right_desc->anc=node;
	left_desc->sister=right_desc;
	right_desc->sister=left_desc;
}

void CutCladeFromTree(NODE *tree, int *clade, int *ntaxa)
{
	NODE *a, *aa, *d, *d2;
	int gobots;
	
	a=tree[*clade].anc;			// a is the selected clade's ancestor and will be removed from tree
	aa=a->anc;					// aa is a's ancestor
	d=tree[*clade].sister;		// d is the selected clade's sister
	
	d->anc=aa;					// d's anc removed, and new anc is aa
	if (aa && a==aa->left) 		// if a is not the root and is the left desc of aa
		{
		aa->left=d;
		d->sister=aa->right;
		d2=aa->right;
		d2->sister=d;
		}
	else if (aa && a==aa->right) // if a is not the root and is the right desc of aa
		{
		aa->right=d;
		d->sister=aa->left;
		d2=aa->left;
		d2->sister=d;
		}
	else if (!aa)				// if a IS the root
		{
		if (&tree[*clade]==a->left)
			{
			d=a->right;
			d->anc=0;
			d->sister=0;
			}
		else if (&tree[*clade]==a->right) {
			d=a->right;
			d->anc=0;
			d->sister=0;
		}
		else  {
			PrintTree(tree, ntaxa);
			cin>>gobots;
		}
	}
	else 
		{
		PrintTree(tree, ntaxa);
		cin>>gobots;
		}

	tree[*clade].left=tree[*clade].right=tree[*clade].anc=tree[*clade].sister=&tree[0];
	tree[*clade].tip=false;
	a->left=a->right=a->anc=a->sister=&tree[0];	
	
}

NODE * GetInitializedTree(int ntaxa) {
	NODE *tree=0;
	int i;
	
	tree=new NODE[2*ntaxa];
	for (i=0;i<2*ntaxa;++i) {
		tree[i].left=tree[i].right=tree[i].anc=tree[i].sister=0;	// set left, right, anc and sister to 0 to start
		tree[i].label=i;							// set the labels
		if (i>0 && i<=ntaxa) tree[i].tip=true;		// tip is 1 for extant taxa
		else tree[i].tip=false;						// zero for nodes
		tree[i].isAnc=false;						// none are ancestorized nor part of a lineage
	}
	return tree;
}

NODE **AllocateTreeList(int ntrees, int ntaxa) {
	NODE **tree_list=0;
	int i;

	// allocate pointers to rows 
	tree_list = new NODE * [ntrees+1];
	if (!tree_list) cout << "allocation error in tree_list - primary";
	
	// allocate rows and set pointers to them 
	for (i=0 ; i<=ntrees ; i++)  {
		tree_list[i] = GetInitializedTree(ntaxa);
		if (!tree_list[i]) cout << "allocation error in tree_list - row " << i;
	}

	return tree_list;
}

/*
void ** AllocateTreeList(NODE ** tree_list, int *ntrees, int *ntaxa)
{
	int i;

	// allocate pointers to rows 
	tree_list = new NODE * [*ntrees+1];
	if (!tree_list) cout << "allocation error in tree_list - primary";
	
	// allocate rows and set pointers to them 
	for (i=0 ; i<=*ntrees ; i++) 
		{
		tree_list[i] = GetInitializedTree(ntaxa);
		if (!tree_list[i]) cout << "allocation error in tree_list - row " << i;
		}
}
*/
void FreeTreeList(NODE **tree_list, int *ntrees)
{
	int i;
	
	for (i=0 ; i<=*ntrees ; i++) delete [] tree_list[i];
	delete [] tree_list;
}

void PrintTree(NODE *tree, int *ntaxa)
{
	int i;
	NODE *hold_node=0;
	
	cout << "\n\tNode\tDesc\tAnc\t\tSister\n";
	cout << "\t----\t----\t---\t\t-----\n";
	for (i=1;i<2*(*ntaxa);++i)
		{
		if (tree[i].anc!=&tree[0])
			{
			cout << "\t" << tree[i].label;
			if (tree[i].isAnc) cout << "*\t\t";
			else cout << "\t\t";
			if (!tree[i].tip) 
				{
				if (!tree[i].left) cout << "(-,";
				else
					{
					hold_node=tree[i].left;		
					cout <<  "(" << hold_node->label << ",";
					}
				if (!tree[i].right) cout << "-)\t";
				else
					{
					hold_node=tree[i].right;		
					cout << hold_node->label << ")\t";
					}
				}
			else cout << "-\t\t";
			if (!tree[i].anc) cout << "ROOT\t\t";
			else
				{
				hold_node=tree[i].anc;		
				cout << hold_node->label << "\t\t";
				if (!tree[i].sister) cout << "-)\t";
				else
					{
					hold_node=tree[i].sister;		
					cout << hold_node->label;
					}
				}
			cout << "\n";
			}
		}
	cout << "\n";
}

void CopyTree(NODE *const source, NODE *const target, int *ntaxa) {
	NODE *root=0;	
	root=GetRoot(source, ntaxa);
	CopyClade(root, target);
}

void CopyClade(NODE *const source_node, NODE *const target) {
	NODE *hold_node=0;

	if (!source_node->tip) {					// source node is not a tip
		hold_node=source_node->left;		
		target[source_node->label].left=&target[hold_node->label];
		CopyClade(hold_node, target);
		
		hold_node=source_node->right;
		target[source_node->label].right=&target[hold_node->label];
		CopyClade(hold_node, target);
	}
	else {									// source_node is a tip
		target[source_node->label].left=0;
		target[source_node->label].right=0;
	}
	if (!source_node->anc) {
		target[source_node->label].anc=0;
		target[source_node->label].sister=0;
	}
	else {
		hold_node=source_node->anc;		
		target[source_node->label].anc=&target[hold_node->label];
		
		hold_node=source_node->sister;		
		target[source_node->label].sister=&target[hold_node->label];
	}
	target[source_node->label].tip=source_node->tip;
	target[source_node->label].isAnc=source_node->isAnc;
	target[source_node->label].label=source_node->label;
}

void UnfixTree(NODE *tree, int *ntaxa)
{
	NODE *root=0;
	
	root=GetRoot(tree, ntaxa);
	UnfixMachine(root);
}

void UnfixMachine(NODE *node)
{
	NODE *left_desc=0, *right_desc=0;
	
	left_desc=node->left;
	right_desc=node->right;
	
	left_desc->Unfix();
	right_desc->Unfix();
	if (!left_desc->tip) UnfixMachine(left_desc);
	if (!right_desc->tip) UnfixMachine(right_desc);
}

int CountTopologies(NODE **trees, int *ntrees, int *ntaxa)
{
	int ntopologies;
	int i, j;
	bool identical;
	bool **cube1, **cube2;
	NODE *root1=0, *root2=0;

	cube1=bmatrix(2**ntaxa, *ntaxa+1);
	cube2=bmatrix(2**ntaxa, *ntaxa+1);
	ntopologies=1;
	for (i=2;i<=*ntrees;++i) {
		root1=GetRoot(trees[i], ntaxa);
		AddTreeToCube(trees[i], cube1, ntaxa);
		identical=false;
		for (j=1;j<i;++j) 
			{
			AddTreeToCube(trees[j], cube2, ntaxa);
			identical=CompareCubes(cube1, cube2, ntaxa);
			if (identical) break;
			}
		if (!identical) ++ntopologies;
		}
	free_matrix(cube1, 2**ntaxa);
	free_matrix(cube2, 2**ntaxa);

	return ntopologies;
}

bool CompareTopologies(NODE *tree1, NODE *tree2, int *ntaxa) {
	bool identical=false;
	bool **cube1, **cube2;

	cube1=bmatrix(2**ntaxa, *ntaxa+1);
	cube2=bmatrix(2**ntaxa, *ntaxa+1);
	AddTreeToCube(tree1, cube1, ntaxa);
	AddTreeToCube(tree2, cube2, ntaxa);
	identical=CompareCubes(cube1, cube2, ntaxa);
	free_matrix(cube1, 2**ntaxa);
	free_matrix(cube2, 2**ntaxa);

	return identical;
}

void MakeStratBL(NODE *tree, DATA *data) {
	NODE *a, *s;
	int i;
	
	for (i=1;i<2**data->GetNTaxa();++i) {
		a=tree[i].anc;
		s=tree[i].sister;
		if (i<=*data->GetNTaxa()) {
			if (!tree[i].isAnc && !s->isAnc) tree[i].bl=(float *)*data->GetLInt(i)-*data->GetFInt(a->label)+1; //tree[i].bl=*data->GetLInt(i)-*data->GetFInt(a->label)+1;
			else {
				if (*data->GetFInt(s->label)>=*data->GetLInt(i)) tree[i].bl=0;
				else tree[i].bl=(float *)*data->GetLInt(i)-*data->GetFInt(s->label);
			}
		}
		else {
			if (!a)	tree[i].bl=0;
			else tree[i].bl=(float *)*data->GetFInt(i)-*data->GetFInt(a->label)+1;
		}
	}
}

void ReadCladeFromFile (ifstream * const finNexus, NODE * const tree, int * const ntaxa)
{
	int nnodes=1;
	char *inchar;
	NODE *tree_holder=0;

	inchar=new char;
	
//	tree_holder=GetInitializedTree(ntaxa);
	while (*inchar!='(') *inchar=finNexus->get();
	ReadCladeMachine (finNexus, tree_holder, 1, nnodes, ntaxa);
	GetAncs(&tree_holder[*ntaxa+1]);
	CopyTree(tree_holder, tree, ntaxa);
	delete [] tree_holder;
	delete inchar;

}

void ReadCladeMachine (ifstream * const finNexus, NODE * const tree, int this_node, int &nnodes, int * const ntaxa)
{
	char *inchar;
	char *holder;
	int counter;
	int desc;
	int i;
	bool fix=false;
	
	holder=new char [4];
	inchar=new char;
	
	for (i=0;i<=1;++i)
		{
		*inchar=finNexus->get();
		if (*inchar=='(')												// if we're starting a new clade
			{
			++nnodes;													// increment the number of nodes
			if (this_node>0)											// if this is not the first node	
				{
				if (!i) tree[this_node+*ntaxa].left=&tree[nnodes+*ntaxa];// next node is the left desc
				else tree[this_node+*ntaxa].right=&tree[nnodes+*ntaxa];	// next node is the right desc
				}
			ReadCladeMachine(finNexus, tree, nnodes, nnodes, ntaxa);	// follow next node
			*inchar=finNexus->get();									// grabbing the ) that closes this node
			}
		else if (*inchar>=48 && *inchar<=57) 							// if the descendant is an OTU
			{
			counter=0;
			while (*inchar>=48 && *inchar<=57) 							// while we're still getting the numerical label
				{
				holder[counter]=*inchar;
				*inchar=finNexus->get();
				++counter;
				}
			finNexus->unget();											// put back the char that wasn't an integer
			holder[counter]='\0';										// NULL terminate the string
			desc=atoi(holder);											// convert string to int
			if (!i) tree[this_node+*ntaxa].left=&tree[desc];	// assign this OTU to left descendant of this node
			else tree[this_node+*ntaxa].right=&tree[desc];		// assign this OTU to right descendant of this node
			if (fix) 
				{
				tree[desc].Fix();								// fix it
				finNexus->unget();	// if fixed, this is the right node, and there will be no closing ), so we need to put the integer back for when this node is closed and a closing ) is expected
				}
			}
//		if (!i && *inchar==')') fix=true;						// if we've just entered the left descendant and we're closing the clade early, its sister is fixed
		if (!i) 
			{
			*inchar=finNexus->get();							// get the first chr after the left desc
			if (*inchar==',') fix=false;						// otherwise if it's a comma, then unfixed
			else if (*inchar==')') fix=true;					// ) if following an integer, it's fixed
			else if (*inchar>=48 && *inchar<=57) 				// the integer if following a ), it's fixed
				{
				fix=true;
				finNexus->unget();								// put that integer back so we can grab it next
				}
//			if (*inchar==')') *inchar=finNexus->get();
			}		
		}
//	if (!i && tree[this_node+*ntaxa].left) *inchar=finNexus->get();
	delete inchar;
	delete [] holder;
}

void AddTreeToCube(NODE *in_tree, bool ** const cube, int *ntaxa)
{
	NODE *root=0;
	int i;
	
	for (i=1;i<=*ntaxa;++i) cube[i][i]=true;
	root=GetRoot(in_tree, ntaxa);
	AddToCubeMachine(root, cube, ntaxa);
}

void AddToCubeMachine(NODE *a, bool ** const cube, int *ntaxa)
{
	NODE *l=0, *r=0;
	int i;
	
	l=a->left;
	if (!l->tip) AddToCubeMachine(l, cube, ntaxa);
	r=a->right;
	if (!r->tip) AddToCubeMachine(r, cube, ntaxa);
	for (i=1;i<=*ntaxa;++i) cube[a->label][i]=cube[l->label][i] or cube[r->label][i];
}

bool CompareCubes(bool ** const cube1, bool ** const cube2, int *ntaxa) {
	int i, j, k;
	bool identical=true;
	
	for (i=(2**ntaxa)-1;i>*ntaxa;--i) {		// for all nodes in cube1
		for (j=(2**ntaxa)-1;j>*ntaxa;--j) {	// for all nodes in cube2
			identical=true;
			k=1;
			while (k<=*ntaxa && cube1[i][k]==cube2[j][k]) ++k;
			if (k<=*ntaxa) identical=false;
			if (identical) break;
		}
		if (!identical) break;
	}	
	return identical;
}

void PrintCube(bool **cube, int *ntaxa)
{
	int i, j;

	for (i=*ntaxa+1;i<2**ntaxa;++i)
		{
		cout << "Node " << i << "\t";
		for (j=1;j<=*ntaxa;++j)
			{
			cout << cube[i][j];
			}
		cout << "\n";
		}
	cout << "\n";
}

bool CompareAncs(NODE *tree1, NODE *tree2, int *ntaxa)
{
	int i;
	bool identical;
	
	identical=true;
	for (i=1;i<=*ntaxa;++i)
		{
		if (tree1[i].isAnc!=tree2[i].isAnc) identical=false;
		if (!identical) break;
		}		
	return identical;
}

void MakeRichness(NODE *node)
{
	NODE *l, *r;
///	int gobots;
	
	l=node->left;
	r=node->right;
	if (!l->richness) l->richness=new int;
	if (!l->tip) MakeRichness(l);
	else *l->richness=1;
	if (!r->richness) r->richness=new int;
	if (!r->tip) MakeRichness(r);
	else *r->richness=1;
	
//	cout << "Node " << node->label << endl;
	if (!node->richness) node->richness=new int;
	*node->richness=*l->richness+*r->richness;
//	cin>>gobots;
}

void LadderizeRight(NODE *tree, int *ntaxa) {
	NODE *root;
//	int gobots;
	
	root=GetRoot(tree, ntaxa);
//	cin>>gobots;
	MakeRichness(root);
//	cin>>gobots;
	LadderizeRightMachine(root);
}

void LadderizeRightMachine(NODE *node) {
	NODE *l, *r;
	
	l=node->left;
	r=node->right;
	if (l->richness>r->richness) {
		node->left=r;
		node->right=l;
		l=node->left;
		r=node->right;
	}
	if (!l->tip) LadderizeRightMachine(l);
	if (!r->tip) LadderizeRightMachine(r);
}
	
void MakeReconstructedRichnessInt(int **richness, DATA *data, NODE *tree) {
	NODE *root;
//	int i;
		
	root=GetRoot(tree, data->GetNTaxa());
	ReconstructedRichnessMachine(data, root, richness);
//	for (i=0;i<*data->GetNStates(data->GetStratChr());++i) cout << "interval " << i << "\t" << richness[i][0] << "\t" << richness[i][1] << endl;
}

void ReconstructedRichnessMachine(DATA *data, NODE * a, int **richness) {
	NODE *l, *r, *f=0, *d=0;
	int i;
	
	l=a->left;
	r=a->right;
	if (l->isAnc) { f=l; d=r; }
	else if (r->isAnc) { f=r; d=l; }
		
	if (!a->anc) {																			// if it's the root
		if (f) {																			// if it's fixed as an anc
			++richness[*data->GetFInt(f)][4];												// record the FA of the anc
			if (d->tip) for (i=*data->GetFInt(f);i<=*data->GetLInt(d);++i) ++richness[i][1];// if the desc is a tip, add the anc from it's FA to the LA of the desc (this can only happen for a two taxon tree
			else for (i=*data->GetFInt(f);i<*data->GetFInt(d);++i) ++richness[i][1];		// else,add the anc from it's FA to just before the FA of the desc
			if (*data->GetLInt(f)>=*data->GetFInt(d)) {										// if the LA of the anc is ³ the FA of the desc (the two taxa overlap in their stratigraphic distributions)
				++richness[*data->GetFInt(d)][4];											// if the anc coexists with its desc, then add the OD of the desc (there are two FAs here, otherwise, there's only one)
				for (i=*data->GetFInt(d);i<=*data->GetLInt(f);++i) ++richness[i][1];		// add richness for the overlapped portion
			}
		}	
		else {																				// if it's the root, but unfixed 
			richness[*data->GetFInt(a)][4]+=2;												// record the two OD's of the two desc, because this node has not yet contributed to richess (neither lineage has been previously recorded)
			if (l->tip) for (i=*data->GetFInt(a);i<=*data->GetLInt(l);++i) ++richness[i][1];// if left is a tip, add richness from the node age to the LA of the left
			else for (i=*data->GetFInt(a);i<*data->GetFInt(l);++i) ++richness[i][1];		// if left is a node, add richness from the FA of the node to just before the FA of the left
			if (r->tip) for (i=*data->GetFInt(a);i<=*data->GetLInt(r);++i) ++richness[i][1];// if right is a tip, add richness from the node age to the LA of the right
			else for (i=*data->GetFInt(a);i<*data->GetFInt(r);++i) ++richness[i][1];		// if right is a node, add richness from the FA of the node to just before the FA of the right
		}
	}
	else {																					// if it's not the root
		if (f) {																			// if it is fixed as an anc
			if (d->tip) for (i=*data->GetFInt(f);i<=*data->GetLInt(d);++i) ++richness[i][1];// if desc is a tip, add richness from the FA of the fixed anc to the LA of the desc
			else for (i=*data->GetFInt(f);i<*data->GetFInt(d);++i) ++richness[i][1];		// else, add richness from the FA of the fixed anc to just before the FA of the desc
			if (*data->GetLInt(f)>=*data->GetFInt(d)) {										// if the LA of the anc is ³ the FA of the desc (the two taxa overlap in their stratigraphic distributions)
				++richness[*data->GetFInt(d)][4];											// add and extra OD for the desc
				for (i=*data->GetFInt(d);i<=*data->GetLInt(f);++i) ++richness[i][1];		// add richness for the overlapped portion
			}
		}	
		else {
			++richness[*data->GetFInt(a)][4];												// if it's not fixed, add the the OD's of the ONE desc - the other is part of the lineage that led to this node
			if (l->tip) for (i=*data->GetFInt(a);i<=*data->GetLInt(l);++i) ++richness[i][1];// if left is a tip, add richness from the node age to the LA of the left
			else for (i=*data->GetFInt(a);i<*data->GetFInt(l);++i) ++richness[i][1];		// if left is a node, add richness from the FA of the node to just before the FA of the left
			if (r->tip) for (i=*data->GetFInt(a);i<=*data->GetLInt(r);++i) ++richness[i][1];// if right is a tip, add richness from the node age to the LA of the right
			else for (i=*data->GetFInt(a);i<*data->GetFInt(r);++i) ++richness[i][1];		// if right is a node, add richness from the FA of the node to just before the FA of the right
		}
	}

	if (l->tip) for (i=*data->GetFInt(l);i<=*data->GetLInt(l);++i) {
		++richness[i][0];
		if (i==*data->GetFInt(l)) ++richness[i][3];
		if (i==*data->GetLInt(l)) ++richness[i][5];
	}
	else ReconstructedRichnessMachine(data, l, richness);	
	if (r->tip) for (i=*data->GetFInt(r);i<=*data->GetLInt(r);++i) {
		++richness[i][0];
		if (i==*data->GetFInt(r)) ++richness[i][3];
		if (i==*data->GetLInt(r)) ++richness[i][5];
	}
	else ReconstructedRichnessMachine(data, r, richness);
}

void MakeReconstructedDisparityInt(float **disparity, DATA *data, NODE *tree) {
	int i, intv;
	int *members;
	
	members=new int [2**data->GetNTaxa()];
	for (intv=0;intv<*data->GetNStates(data->GetStratChr());++intv) {
//		cout<<"Interval "<<intv<<"\t";

		for (i=0;i<2**data->GetNTaxa();++i) members[i]=0;			// initialize members
		MakeObservedMembership(data, members, intv);
		disparity[intv][0] = MakeDisparity(data, members);
//		cout<<"Members:\t";
//		for (i=1;i<=members[0];++i) cout<<members[i]<<" ";
//		cout<<endl;

		for (i=0;i<2**data->GetNTaxa();++i) members[i]=0;			// reinitialize members
		MakeReconstructedMembership(data, tree, members, intv);
		disparity[intv][1] = MakeDisparity(data, members);
//		cout<<"Members:\t";
//		for (i=1;i<=members[0];++i) cout<<members[i]<<" ";
//		cout<<endl;
	}
	delete [] members;
}

void MakeObservedMembership(DATA *data, int *members, int intv) {
	int i;

	for (i=1;i<=*data->GetNTaxa();++i) {
		if (*data->GetFInt(i)<=intv && *data->GetLInt(i)>=intv) {
			++members[0];
			members[members[0]]=i;
		}
	}
}

void MakeReconstructedMembership(DATA *data, NODE *tree, int *members, int intv) {
	int i;
	NODE *a, *s;

	for (i=1;i<2**data->GetNTaxa();++i) {	// for all taxa
		a=tree[i].anc;
		if (a) {							// if this desc is not part of the root
			s=tree[i].sister;
			if (!s->isAnc) {				// if the sister is not fixed
				if ((i<=*data->GetNTaxa() && (*data->GetFInt(a->label)<=intv && *data->GetLInt(i)>=intv)) ||// if this taxon is a tip, and intv is between its node's age and its LA
					 (i>*data->GetNTaxa() && (*data->GetFInt(a->label)<=intv && intv<*data->GetFInt(i)))) { // i has two descendants that form the instant its FInt begins, so its LInt is the interval before
					++members[0];
					members[members[0]]=i;
				}
			}
			else { 				// if the sister IS fixed
				if (*data->GetLInt(s->label)>=*data->GetFInt(i)) { 	// if the sister's LInt is equal or after the FInt of the desc
					if ((i<=*data->GetNTaxa() && (*data->GetFInt(i)<=intv && *data->GetLInt(i)>=intv)) // ||
						/* (i>*data->GetNTaxa() && (*data->GetFInt(i)<=intv && intv<*data->GetFInt(i)))*/) { // i has two descendants that form the instant its FInt begins, so its LInt is the interval before 
						++members[0];
						members[members[0]]=i;
					}
				}
				else {			// the sister's LInt is prior to the FInt of the desc
					if ((i<=*data->GetNTaxa() && (*data->GetLInt(s->label)<intv && *data->GetLInt(i)>=intv))  ||	// ith lineage goes between the interval after the anc's LA and i's LA
						 (i>*data->GetNTaxa() && (*data->GetLInt(s->label)<intv && intv<*data->GetFInt(i)))) { // i has two descendants that form the instant its FInt begins, so its LInt is the interval before
						++members[0];
						members[members[0]]=i;
					}
				}
			}
		}
	}
}

float MakeDisparity(DATA *data, int *members) {
	int i, j, k;
	float comps_c=0.0, comps_t=0.0, diffs=0.0, mpwd=0.0;
//	int gobots;

	for (i=1;i<members[0];++i) {		// for all taxa (except the last)
		for (j=i+1;j<=members[0];++j) {		// for all subsequent taxa
			++comps_t;							// assume this will be a valid comparison between two taxa in the same interval
			comps_c=diffs=0;					// set number of comparisons to zero
			for (k=1;k<=*data->GetNChrs();++k) {	// for all chrs
				if (data->GetInc(k) && *data->GetChrType(k)!=4) { // use only included chrs that are not strat chrs.
//					cout << k << endl;
					if (*data->GetMatrix(members[i], k) && *data->GetMatrix(members[j], k)) { // if both taxa do not have missing data
						++comps_c;									// this is a valid comparison
						if (!(*data->GetMatrix(members[i], k) & *data->GetMatrix(members[j], k))) {
							 ++diffs;	// if they don't share any states, then they're different
//							 cout << members[i] << "\t" << members[j] << "\t" << k << endl;
						}
					}
				}
			}
//			cin >> gobots;
			if (!comps_c) --comps_t;	// if there were no comparable chrs then this wasn't a valid comparison after all
			else mpwd+=(diffs/comps_c);	// average distance between i and j
		}
	}
	if (comps_t) mpwd=mpwd/comps_t;	// average distance over all taxon comparisons
	else mpwd=0;
	
	return mpwd;
}

void MakeReconstructedChrEvInt(float **chr_ev, DATA *data, NODE *tree) {
	int i, j, intv;
	int comps_c, changes, comps_t;
	int gobots;
	NODE *f, *d;
	
//	data->PrintStrat();
//	PrintTree(tree, data->GetNTaxa());
	
	for (intv=0;intv<*data->GetNStates(data->GetStratChr());intv++) { // for all intervals;
		chr_ev[intv][0]=0;
		comps_t=0;
		for (i=*data->GetNTaxa()+1;i<2**data->GetNTaxa();i++) {				// for all nodes on the tree
			f=d=0;
			if (tree[i].left->isAnc) { f=tree[i].left; d=tree[i].right; }
			else if (tree[i].right->isAnc) { f=tree[i].right; d=tree[i].left; }
			if (!f && *data->GetFInt(i)==intv) {	// if the node has a FInt in this interval and is not the root
//				cout << "Comparing Taxon " << i << " to Taxon " << tree[i].left->label << " in interval " << intv << endl;
				++comps_t;									// add this as a valid taxon comp
				comps_c=changes=0;							// reset changes and comparisons
				for (j=1;j<=*data->GetNChrs();j++) {		// for all possible characters
					if (*data->GetChrType(j)!=4 && !data->IsMissing(tree[i].left->label,j)) {			// if this node is not missing this character
						++comps_c;							// add this as a valid character comparison
						if (!data->ShareAnyState(i, tree[i].left->label, j)) ++changes;	// if there has been an unambiguous character change, increase the number of changes
					}
				}
				if (comps_c) chr_ev[intv][0]+=((float)changes/(float)comps_c); 	// number of changes per comparison
//				cout << "Comparing Taxon " << i << " to Taxon " << tree[i].right->label << " in interval " << intv << endl;
				++comps_t;									// add this as a valid taxon comp
				comps_c=changes=0;							// reset changes and comparisons
				for (j=1;j<=*data->GetNChrs();j++) {		// for all possible characters
					if (*data->GetChrType(j)!=4 && !data->IsMissing(tree[i].right->label,j)) {			// if this node is not missing this character
						++comps_c;							// add this as a valid character comparison
						if (!data->ShareAnyState(i, tree[i].right->label, j)) ++changes;	// if there has been an unambiguous character change, increase the number of changes
					}
				}
				if (comps_c) chr_ev[intv][0]+=((float)changes/(float)comps_c); 	// number of changes per comparison
			}
			else if (f && *data->GetFInt(d->label)==intv) {
//				cout << "Comparing Taxon " << f->label << "* to Taxon " << d->label << " in interval " << intv << endl;
				++comps_t;									// add this as a valid taxon comp
				comps_c=changes=0;							// reset changes and comparisons
				for (j=1;j<=*data->GetNChrs();j++) {		// for all possible characters
					if (*data->GetChrType(j)!=4 && !data->IsMissing(d->label,j)) {			// if this node is not missing this character
						++comps_c;							// add this as a valid character comparison
						if (!data->ShareAnyState(i, d->label, j)) ++changes;	// if there has been an unambiguous character change, increase the number of changes
					}
				}
				if (comps_c) chr_ev[intv][0]+=((float)changes/(float)comps_c); 	// number of changes per comparison
			}
		}
		if (comps_t) chr_ev[intv][0]/=(float)comps_t;						// average proportion of changes per lineage
		else chr_ev[intv][0]=0;
	}
//	cin >> gobots;
}

void MakeReconstructedChrEvIntLA(float **chr_ev, DATA *data, NODE *tree) {
	int i, j, intv;
	int comps_c, changes, comps_t;
	int gobots;
	NODE *f, *d;
	
//	data->PrintStrat();
//	PrintTree(tree, data->GetNTaxa());
	
	for (intv=0;intv<*data->GetNStates(data->GetStratChr());intv++) { // for all intervals;
		chr_ev[intv][1]=0;
		comps_t=0;
		for (i=*data->GetNTaxa()+1;i<2**data->GetNTaxa();i++) {				// for all nodes on the tree
			f=d=0;
			if (tree[i].left->isAnc) { f=tree[i].left; d=tree[i].right; }
			else if (tree[i].right->isAnc) { f=tree[i].right; d=tree[i].left; }
			if (!f) {
				if (*data->GetFInt(tree[i].label)==*data->GetFInt(tree[i].left->label)==intv || 
					(*data->GetFInt(tree[i].label)!=*data->GetFInt(tree[i].left->label) && intv==*data->GetFInt(tree[i].left->label)-1)) {	// if the node has a FInt in this interval and is not the root
	//				cout << "Comparing Taxon " << i << " to Taxon " << tree[i].left->label << " in interval " << intv << endl;
					++comps_t;									// add this as a valid taxon comp
					comps_c=changes=0;							// reset changes and comparisons
					for (j=1;j<=*data->GetNChrs();j++) {		// for all possible characters
						if (*data->GetChrType(j)!=4 && !data->IsMissing(tree[i].left->label,j)) {			// if this node is not missing this character
							++comps_c;							// add this as a valid character comparison
							if (!data->ShareAnyState(i, tree[i].left->label, j)) ++changes;	// if there has been an unambiguous character change, increase the number of changes
						}
					}
					if (comps_c) chr_ev[intv][1]+=((float)changes/(float)comps_c); 	// number of changes per comparison
				}
				if (*data->GetFInt(tree[i].label)==*data->GetFInt(tree[i].right->label)==intv || 
					(*data->GetFInt(tree[i].label)!=*data->GetFInt(tree[i].right->label) && intv==*data->GetFInt(tree[i].left->label)-1)) {	// if the node has a FInt in this interval and is not the root
	//				cout << "Comparing Taxon " << i << " to Taxon " << tree[i].right->label << " in interval " << intv << endl;
					++comps_t;									// add this as a valid taxon comp
					comps_c=changes=0;							// reset changes and comparisons
					for (j=1;j<=*data->GetNChrs();j++) {		// for all possible characters
						if (*data->GetChrType(j)!=4 && !data->IsMissing(tree[i].right->label,j)) {			// if this node is not missing this character
							++comps_c;							// add this as a valid character comparison
							if (!data->ShareAnyState(i, tree[i].right->label, j)) ++changes;	// if there has been an unambiguous character change, increase the number of changes
						}
					}
					if (comps_c) chr_ev[intv][1]+=((float)changes/(float)comps_c); 	// number of changes per comparison
				}
			}
			else if (f && *data->GetFInt(d->label)==intv) {
//				cout << "Comparing Taxon " << f->label << "* to Taxon " << d->label << " in interval " << intv << endl;
				++comps_t;									// add this as a valid taxon comp
				comps_c=changes=0;							// reset changes and comparisons
				for (j=1;j<=*data->GetNChrs();j++) {		// for all possible characters
					if (*data->GetChrType(j)!=4 && !data->IsMissing(d->label,j)) {			// if this node is not missing this character
						++comps_c;							// add this as a valid character comparison
						if (!data->ShareAnyState(i, d->label, j)) ++changes;	// if there has been an unambiguous character change, increase the number of changes
					}
				}
				if (comps_c) chr_ev[intv][1]+=((float)changes/(float)comps_c); 	// number of changes per comparison
			}
		}
		if (comps_t) chr_ev[intv][1]/=(float)comps_t;						// average proportion of changes per lineage
		else chr_ev[intv][1]=0;
	}
//	cin >> gobots;
}

