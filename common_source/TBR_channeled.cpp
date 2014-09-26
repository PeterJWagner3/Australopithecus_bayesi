#include "TBR_channeled.h"

using namespace std;
HeuristicSearchClass::HeuristicSearchClass(int * const no_taxa) {
	ntaxa=no_taxa;			// ntaxa become a fixed constant
	MaxTrees=1001;
	ChannelSize=0;
	reps=10;
	AncSearchType=1;		// 0=Exhaustive, 1=Heuristic
	AncSearchReps=1;		// Number of ancestral search reps
	MaxSavingsMethod=1;		// 0=Maximum Stratigraphic Savings, 1=Heuristically Estimated Total Savings
	OutgroupSpecified=false;
	SingleOutgroup=false;
	AllocateOutgroup();
	nstarttrees=0;
	starttrees=0;
	}
HeuristicSearchClass::~HeuristicSearchClass() {}
void HeuristicSearchClass::SetMaxTrees(int i) { MaxTrees=i; }
int *HeuristicSearchClass::GetMaxTrees() { return &MaxTrees; }
void HeuristicSearchClass::SetChannel(int i) { ChannelSize=i; }
int *HeuristicSearchClass::GetChannel() { return &ChannelSize; }
void HeuristicSearchClass::SetSearchReps(int i) { reps=i; }
void HeuristicSearchClass::SetAncSearchReps(int i) { AncSearchReps=i; }
void HeuristicSearchClass::SetAncSearchType(bool i) { AncSearchType=i; }
void HeuristicSearchClass::SetMaxSavingsMethod(bool i) { MaxSavingsMethod=i; }
void HeuristicSearchClass::SetStartTrees(NODE **start_trees, int *ntrees) {
	starttrees=AllocateTreeList(ntrees, ntaxa);							// allocate tree list
	for (int i=1;i<=*ntrees;++i) CopyTree(start_trees[i], starttrees[i], ntaxa);
	nstarttrees=*ntrees;
	reps=*ntrees;
	}

void HeuristicSearchClass::HeuristicSearch(DATA * const data, NODE ** const OptTrees, int &nopttrees){
	NODE *hold_node=0, *hold_node2=0;
	NODE *root=0;
	WriteNexus TempTreeOutput;
	char *tree_filename=0;
	bool island_found;
	int i, j;
		
	random=true;
	if (!random) {	
		srand(1);
		cout  << "\n\t*** WARNING: Random is off ***\n";
	}
	else srand((unsigned int) time((time_t*)NULL)); 
	
	cout<<"\n\n***********************************************\nBeginning Heuristic Search ("<<reps<<" RandSeqAdd reps";	// SUPPRESS ME!
	if (data->hasStratChr && AncSearchType) cout << ", "<<AncSearchReps<<" AncSearch reps)\n";	// SUPPRESS ME!
	else if (data->hasStratChr && !AncSearchType) cout << ", Exhaustive AncSearch)\n";	// SUPPRESS ME!
	else cout << ")\n";
	nislands=0;
	island=new IslandData [reps+1];
	tree_filename=new char [20];										// allocate name of temp tree file
	tree_filename="trees_temp.nex\0";									// set name of temp tree file
	MaxTrees+=1;
	tree_list=AllocateTreeList(&MaxTrees, ntaxa);							// allocate tree list
	MaxTrees-=1;
	tree=GetInitializedTree(ntaxa);
	tree2=GetInitializedTree(ntaxa);
	tree3=GetInitializedTree(ntaxa);
	if (data->hasStratChr) tree4=GetInitializedTree(ntaxa);
	tree_cube=bcube(2*MaxTrees, 2**ntaxa, *ntaxa+1);						// allocate matrix representation of tree list
	this_cube=bmatrix(2**ntaxa, *ntaxa+1);								// allocate matrix representation of this tree
	add_seq=new int [*ntaxa+1];										// allocate addition sequence array
	nopttrees=0;
	if (!ChannelSize) ChannelSize=MaxTrees;
	initializeOutgroup(data);
	
	for (current_rep=1;current_rep<=reps;++current_rep) {
		abort_rep=false;
		trees_this_rep=1;
		swapped_trees=0;
		if (!starttrees) {
			if (!random) PutSimpleAdditionSequence(add_seq, ntaxa);	// sequential order
			else PutRandomAdditionSequence(add_seq, ntaxa);			// random order
			GetStepwiseAdditionTree(tree_list[1], add_seq, data);	// starting tree for rep is next tree in list past currently stored optimal trees
		}
		else CopyTree(starttrees[current_rep], tree_list[1], ntaxa);
		Reroot(tree_list[1], &tree_list[1][root_taxon], ntaxa);		// reroot stepwise addition tree
		root=GetRoot(tree_list[1], ntaxa);
//		if (data->hasStratChr) parsimony.DateTree(root, data);		// nodes must have FA's for strat chrs to work
		best_debt_this_rep=parsimony.TotalDebt(root, data);			// get debt of stepwise tree
		AddTreeToCube(tree_list[1], tree_cube[nopttrees+1], ntaxa);	// add stepwise tree to "cube"
		tree_cube[nopttrees+1][0][0]=true;							// the topology of this tree is a new one
		if (current_rep==1) best_debt=best_debt_this_rep;			// if this is the first rep, the stepwise tree debt is the best overall debt
		cout << "\nAddSeqRep " << current_rep << " of " << reps << "\t" << nopttrees << " optimal trees\t" << best_debt << " total debt\n";	// SUPPRESS ME!

		if (data->hasStratChr) {						// if stratocladistic search, do initial ancestor search
			CopyTree(tree_list[1], tree, ntaxa);		// copy AddSeq tree to dummy tree to be rerooted
			if (SingleOutgroup) {
				AddTreeToCube(tree, this_cube, ntaxa);
				for (i=1;i<=nopttrees;++i) {			// compare stepwise tree with optimal trees from previous reps
					if (tree_cube[i][0][0] && CompareCubes(this_cube, tree_cube[i], ntaxa)) {	// tree was identical with one from a previous rep, dump this replicate
						abort_rep=true;
						for (j=1;j<=nislands;++j) if (island[j].first_tree<=i && island[j].last_tree>=i) ++island[j].times_hit;
						break;						
					}
				}
				if (!abort_rep) {
					root=GetRoot(tree, ntaxa);
					added_topology=true;	// we've already added an unanc'ed tree of this topology
					if (!AncSearchType) FullAncSearch(tree, root, 1, data, nopttrees);
					else {
						for (i=1;i<=AncSearchReps;++i) {
							if (!random) PutSimpleAdditionSequence(add_seq, ntaxa);			// sequential order
							else PutRandomAdditionSequence(add_seq, ntaxa);					// random order
							HeuristicAncSearch(tree, root, 1, data, nopttrees);
						}
					}
				}
			}
			else for (j=1;j<2**ntaxa;++j) {					// try rerooting on all possible nodes
				root=GetRoot(tree, ntaxa);
				if (tree[j].anc!=root) {
					CopyTree(tree, tree4, ntaxa);			// copy AddSeq tree to dummy tree to be rerooted
					Reroot(tree4, &tree4[j], ntaxa);		// reroot the dummy tree
					root=GetRoot(tree4, ntaxa);
//					if (data->hasStratChr) parsimony.DateTree(root, data);		// nodes must have FA's for strat chrs to work
					if (OutgroupSpecified) {				// if there is an outgroup
						MakeOutgroupFromTree(root);			// find all outgroup nodes on this tree
						hold_node=root->left;
						hold_node2=root->right;
					}
					if (!OutgroupSpecified || (OutgroupSpecified && (GetOutgroup(hold_node->label) || GetOutgroup(hold_node2->label)))) {
						AddTreeToCube(tree4, this_cube, ntaxa);
						for (i=1;i<=nopttrees;++i) {					// compare stepwise tree with optimal trees from previous reps
							abort_rep=CompareCubes(this_cube, tree_cube[i], ntaxa);	// tree was identical with one from a previous rep, dump this replicate
							if (abort_rep) for (j=1;j<=nislands;++j) if (island[j].first_tree<=i && island[j].last_tree>=i) ++island[j].times_hit;
							if (abort_rep) break;
						}
						if (!abort_rep) {
							added_topology=false;		// could allow two trees of same topology - not a big deal (?)
							if (!AncSearchType) FullAncSearch(tree4, root, 1, data, nopttrees);
							else {
								for (i=1;i<=AncSearchReps;++i) {
									if (!random) PutSimpleAdditionSequence(add_seq, ntaxa);			// sequential order
									else PutRandomAdditionSequence(add_seq, ntaxa);					// random order
									HeuristicAncSearch(tree4, root, 1, data, nopttrees);
								}
							}
						}
					}
				}
			if (abort_rep) break;
			} // for j
		}
		else {		// if there is no strat chr
			for (i=1;i<=nopttrees;++i) {	// compare stepwise tree with optimal trees from previous reps
				abort_rep=CompareCubes(tree_cube[nopttrees+1], tree_cube[i], ntaxa);	// tree was identical with one from a previous rep, dump this replicate
				if (abort_rep) {
					for (j=1;j<=nislands;++j) {
						if (island[j].first_tree<=i && island[j].last_tree>=i) ++island[j].times_hit;
					}
				}
				if (abort_rep) break;
			}
		}

		if (!abort_rep) {
			cout << "\tBeginning search with "<< trees_this_rep <<" tree(s) of length "<< best_debt_this_rep << "\n";	// SUPPRESS ME!
			SingleHeuristicSearch(data, nopttrees);
		}
		else cout << "\t*** Random Addition Topology was identical to one previously found\n";	// SUPPRESS ME!

		if (!abort_rep) {
			if (best_debt_this_rep<best_debt) {
				nopttrees=0;						// reset number of optimal trees
				best_debt=best_debt_this_rep;		// reset the best debt
				for (i=1;i<=nislands;++i) {
					island[i].first_tree=-1;
					island[i].last_tree=-1;
				}
			}
			if (best_debt_this_rep==best_debt) {
				++nislands;
				island[nislands].first_tree=nopttrees+1;
				island[nislands].last_tree=nopttrees+trees_this_rep;
				island[nislands].ntrees=trees_this_rep;
				island[nislands].ntopologies=CountTopologies(tree_list, &trees_this_rep, ntaxa);
				island[nislands].times_hit=1;
				island[nislands].debt=best_debt_this_rep;
				for (i=1;i<=trees_this_rep;++i) {
					if (nopttrees<MaxTrees-1) {
						added_topology=false;
						++nopttrees;
						for (j=1;j<*ntaxa;++j) tree_cube[nopttrees][j][0]=0;			// reset cube of nopttrees'th cube
						AddTreeToCube(tree_list[i], tree_cube[nopttrees], ntaxa);
						for (j=nopttrees-1;j>0;--j) {
							if (tree_cube[j][0][0] && CompareCubes(tree_cube[nopttrees], tree_cube[j], ntaxa)) {
								tree_cube[nopttrees][0][0]=false;
								added_topology=true;
							}
							if (added_topology) break;
						}
						if (!added_topology) tree_cube[nopttrees][0][0]=true;
						CopyTree(tree_list[i], OptTrees[nopttrees], ntaxa);
					}
//					else cout << "\n\t***\tMAXTREES met\n";
				}
				TempTreeOutput.WriteNexusTreeFile (OptTrees, &nopttrees, ntaxa, data->taxon_label, tree_filename);

 /*//added for debugging
				int gobots;
				PrintTree(OptTrees[1], ntaxa);
				root=GetRoot(OptTrees[1], ntaxa);
				data->PrintPolymorph();
				data->PrintNStatesPerTaxon();
//				parsimony.DateTree(root, data);
				cout << parsimony.TotalDebt(root, data) << endl;
				cout << best_debt_this_rep << endl;
				for (i=1;i<=*data->GetNChrs();++i) {
					if (*data->GetChrType(i)!=4) cout<<"\t"<<i<<"\t"<< parsimony.SingleChrDebt(OptTrees[1], data, &i)<<"\n";
					else cout<<"\t"<<i<<"\t"<<parsimony.GetStratDebt(root, data)<<"\n";
				}
				cin>>gobots;
				OptTrees[1][1].Fix();
				parsimony.DateTree(root, data);
				cout << parsimony.TotalDebt(root, data) << endl;
				for (i=1;i<=*data->GetNChrs();++i) {
					if (*data->GetChrType(i)!=4) cout<<"\t"<<i<<"\t"<< parsimony.SingleChrDebt(OptTrees[1], data, &i)<<"\n";
					else cout<<"\t"<<i<<"\t"<<parsimony.GetStratDebt(root, data)<<"\n";
				}
				cin>>gobots;
// end added for debugging
*/				
			}
			else {
				cout << "\t*** Unoptimal Island found with " << trees_this_rep << " trees of length " << best_debt_this_rep << "\n";	// SUPPRESS ME!
				island_found=false; 
				for (i=1;i<=nislands;++i) {
					if (best_debt_this_rep==island[i].debt && trees_this_rep==island[i].ntrees) island_found=true;
					if (island_found) ++island[i].times_hit;
					if (island_found) break;
				}
				if (!island_found) {
					++nislands;
					island[nislands].first_tree=-1;
					island[nislands].last_tree=-1;
					island[nislands].times_hit=1;
					island[nislands].debt=best_debt_this_rep;
					island[nislands].ntrees=trees_this_rep;
					island[nislands].ntopologies=CountTopologies(tree_list, &trees_this_rep, ntaxa);
				}
			}
		}
	}

	PrintIslandData();
	FreeTreeList(tree_list, &MaxTrees);
	free_cube(tree_cube, 2*MaxTrees, 2**ntaxa);
	free_matrix(this_cube, 2**ntaxa);
	delete [] add_seq;
	delete [] tree;
	delete [] tree2;
	delete [] tree3;
	if (data->hasStratChr) delete [] tree4;
	delete [] outgroup;
	delete [] island;
}

void HeuristicSearchClass::SingleHeuristicSearch(DATA * const data, int &nopttrees)
{
	NODE *root=0;
	int i;
	
	while (swapped_trees<trees_this_rep) {
		abort_swap=0;
		++swapped_trees;
		root=GetRoot(tree_list[swapped_trees], ntaxa);
		cout << "\tSwapping Tree " << swapped_trees <<" of "<< trees_this_rep;	// SUPPRESS ME!
		if (swapped_trees>1) {
			AddTreeToCube(tree_list[swapped_trees], this_cube, ntaxa);
			for (i=1;i<swapped_trees;++i) {	// do not swap on trees identical to ones we've already swapped in the same rep
//				if (tree_cube[nopttrees+i][0][0] && CompareCubes(this_cube, tree_cube[nopttrees+i], ntaxa)) abort_swap=true;	// tree has same cladistic topology as another one in the set (needed for strato trees with same topology, but different anc'ings - don't need to swap those again)
				if (!tree_cube[nopttrees+swapped_trees][0][0]) abort_swap=true;	// tree has same cladistic topology as another one in the set (needed for strato trees with same topology, but different anc'ings - don't need to swap those again)
				if (abort_swap) cout << "\t Identical topology\n";	// SUPPRESS ME!
				if (abort_swap) break;
			}
		}
		if (!abort_swap) {
			root=GetRoot(tree_list[swapped_trees], ntaxa);
			GetPrunedNode(root, data, nopttrees);
		}
		if (!abort_swap) cout << "\n";	// SUPPRESS ME!
		if (abort_rep) return;
	}
}

// Recursive function to find all possible subtrees to prune on a tree
void HeuristicSearchClass::GetPrunedNode(NODE * const current_node, DATA * const data, int &nopttrees)
{
	NODE *pruned_node=0, *tree_pruned_node=0, *target_node=0, *hold_node=0, *root=0;
	int i;
	int gobots;

	for (i=0;i<=1;++i) {
		if (!abort_swap) {
			if (!i) pruned_node=current_node->left;						// Swap subtree that is left branch of this node
			else pruned_node=current_node->right;						// Swap subtree that is right branch of this node

			if (!pruned_node->tip && !(pruned_node->left->tip && pruned_node->right->tip)) ReRootPrunedNode(pruned_node, pruned_node, data, nopttrees);
			else {
				CopyTree(tree_list[swapped_trees], tree3, ntaxa);			// copy original to dummy
				UnfixTree(tree3, ntaxa);										// unfix the entire tree before swapping branches
				tree_pruned_node=&tree3[pruned_node->label];					// copy the pruned node to the dummy tree
				RerootWithoutReorder(tree3, tree_pruned_node, ntaxa);		// reroot tree on the pruned node
				root=GetRoot(tree3, ntaxa);
//				if (data->hasStratChr) parsimony.DateNodes(root, data);		// nodes must have FA's for strat chrs to work
				hold_node=root->left;										// find node to start searching for target nodes (opposite side of the root from the pruned node)
				if (hold_node==tree_pruned_node) target_node=root->right;
				else {
					hold_node=root->right;
					if (hold_node==tree_pruned_node) target_node=root->left;
					else {
						cout << "\t\tUnable to obtain target node\n";
						PrintTree(tree3, ntaxa);
						cin >> gobots;
					}
				}
				if (!target_node->tip) GetTargetNode(tree3, tree_pruned_node, target_node, data, nopttrees);
				// do I need to add a step that attaches pruned to tips?
			}
			if (abort_swap) return;
			else if (!pruned_node->tip) GetPrunedNode(pruned_node, data, nopttrees);	// If this pruned node is not a tip, then follow the branch for more subtrees to swap
		}
		if (abort_swap) return;
	}
}

// Recursive function to find all possible rootings of the pruned node
void HeuristicSearchClass::ReRootPrunedNode(NODE * const current_node, NODE * const pruned_node, DATA * const data, int &nopttrees)
{
	NODE *new_root=0, *tree_new_root=0, *hold_node=0, *root=0, *tree_pruned_node=0, *target_node=0;
	int i;
	int gobots;

	for (i=0;i<=1;++i) {
		if (!abort_swap) {
			if (!i) new_root=current_node->left;					// Reroot on left branch of this current_node
			else new_root=current_node->right;						// Swap subtree that is right branch of this node

			if (current_node!=pruned_node || (current_node==pruned_node && !i)) {	// begin with the SPR swap (i.e., just this node without rerooting)
				CopyTree(tree_list[swapped_trees], tree2, ntaxa);		// copy original to dummy
				UnfixTree(tree2, ntaxa);								// unfix the entire tree before swapping branches
				tree_pruned_node=&tree2[pruned_node->label];			// copy the pruned node to the dummy tree
				RerootWithoutReorder(tree2, tree_pruned_node, ntaxa);	// reroot tree on the pruned node
				root=GetRoot(tree2, ntaxa);
//				if (data->hasStratChr) parsimony.DateNodes(root, data);		// nodes must have FA's for strat chrs to work

				tree_new_root=&tree2[new_root->label];
				
				hold_node=root->left;										// find node to start searching for target nodes (opposite side of the root from the pruned node)
				if (hold_node==tree_pruned_node) target_node=root->right;
				else {
					hold_node=root->right;
					if (hold_node==tree_pruned_node) target_node=root->left;
					else {
						cout << "\t\tUnable to obtain target node\n";
						PrintTree(tree2, ntaxa);
						cin >> gobots;
					}
				}
				RerootClade(tree2, tree_pruned_node, tree_new_root, ntaxa);
				if (!target_node->tip) GetTargetNode(tree2, tree_pruned_node, target_node, data, nopttrees);
			}

			if (abort_swap) return;
			else if (!new_root->tip) ReRootPrunedNode(new_root, tree_pruned_node, data, nopttrees);	// If this pruned node is not a tip, then follow the branch for more subtrees to swap
		}
		if (abort_swap) return;
	}
}

// Recursive function to find all possible places on the remaining tree to regraft the subtree
void HeuristicSearchClass::GetTargetNode(NODE * const in_tree, NODE * const pruned_node, NODE * const current_node, DATA * const data, int &nopttrees)
{
	int i;
	NODE *root=0, *target_node=0;
	
	for (i=0;i<=1;++i) {
		if (!abort_swap) {
		  	if (!i) target_node=current_node->left;		// first check left branch
			else target_node=current_node->right;		// then check the right branch
			root=GetRoot(in_tree, ntaxa);
			if (current_node->anc!=root) RegraftSubtree(in_tree, &pruned_node->label, &target_node->label, data, nopttrees);
			if (abort_swap) return;
			else if (!target_node->tip) GetTargetNode(in_tree, pruned_node, target_node, data, nopttrees);
		}
		if (abort_swap) return;
	}
}

void HeuristicSearchClass::RegraftSubtree(NODE * const current_tree, int * const pruned_node, int * const target_node, DATA * const data, int &nopttrees)
{
	NODE *root=0,*hold_node=0;		//
	NODE *cut_node=0;		// node that is removed when subtree is pruned
	NODE *drop_node=0;		// node that replaces cut node (pruned_node's sister)
	float current_debt=0;
	int cut_node_label=0;
	int i, j;
	bool lights;
	int gobots;
	
	lights=false;
	CopyTree(current_tree, tree, ntaxa);		// copy tree to be swapped - original tree is untouched
	
  // Get Cut_Node
	cut_node=tree[*pruned_node].anc;	// node that is cut out is the ancestor of the pruned_node
	cut_node_label=cut_node->label;		// keep track of it's label because we're going to reuse it for the new node
	if (lights) cout << "Cut_Node = " << cut_node->label << "\n";

  // Get Drop_Node
  	drop_node=tree[*pruned_node].sister;
	if (lights) cout << "Drop_Node = " << drop_node->label << "\n";

  // Get Hold_Node
	if (cut_node->anc) {								// as long as the pruned_node's anc is not ROOT
		hold_node=cut_node->anc;						// 	cut_node's ancestor will receive cut_node's other descendant (drop_node) 
		if (hold_node->left==cut_node) hold_node->left=drop_node;
		else if (hold_node->right==cut_node) hold_node->right=drop_node;
		else {
			cout << "\n\n***\tCan't Obtain Hold_Node\n";
			cin >> gobots;
		}
		if (lights) cout << "Hold_Node = " << hold_node->label << "\n";
		drop_node->anc=hold_node;					// move the cut node's ancestor to the drop node
	}
	else {												// however, if it is ROOT, 
		drop_node->anc=0;								// just cut it out, and make the drop_node the root
		root=drop_node;
		if (lights) cout << "Hold_Node = ROOT\n";
	}

  // Place Subtree on Requested Branch
	tree[cut_node_label].right=&tree[*target_node];	// new node's right branch is the target node
	tree[cut_node_label].left=&tree[*pruned_node];	// new node's left branch is the clipped subtree
	tree[cut_node_label].anc=tree[*target_node].anc;// new node's ancestor is the target node's OLD anc
	tree[*pruned_node].anc=&tree[cut_node_label];	// subtree's new anc is the new node
	hold_node=tree[*target_node].anc;				// hold node is the target_node's OLD anc for now
	if (hold_node->left==&tree[*target_node]) hold_node->left=&tree[cut_node_label];			// target node's old anc's left is the new node
	else if (hold_node->right==&tree[*target_node]) hold_node->right=&tree[cut_node_label];	// target node's old anc's right is the new node
	else {
		cout << "\n\n***\tRegraft Failed\n";
		cin >> gobots;
	}
	tree[*target_node].anc=&tree[cut_node_label];	// target node's new anc is also the new node

//	++rearrangements;
	Reroot(tree, &tree[root_taxon], ntaxa);							// reroot new topology on the outgroup
	if (data->hasStratChr) beginAncSearch (data, nopttrees);
	else {							// if there is no strat character
		root=GetRoot(tree, ntaxa);									// need root from this new topology
		current_debt=parsimony.TotalDebt(root, data);
		if (current_debt<=best_debt_this_rep) {
			AddTreeToCube(tree, this_cube, ntaxa);
			if (current_debt==best_debt) {
				for (i=1;i<=nopttrees;++i) {							// compare with previously discovered optimal trees
					if (CompareCubes(this_cube, tree_cube[i], ntaxa)) {	// if tree was identical with a previously discovered one
						abort_rep=true; 	// abort this rep
						abort_swap=true;	// abort swapping on this tree
						for (j=1;j<=nislands;++j) if (island[j].first_tree<=i && island[j].last_tree>=i) ++island[j].times_hit;
						cout << "\t*** Found a previous optimal cladogram\n";	// SUPPRESS ME!
						return;
					}
				}
			}
			else if (current_debt<best_debt) {
				nopttrees=0;
				best_debt=current_debt;
			}
			if (current_debt==best_debt_this_rep) {
				for (i=1;i<=trees_this_rep;++i) 
					if (CompareCubes(this_cube, tree_cube[nopttrees+i], ntaxa)) return;// if tree was identical with a previously discovered one from this rep
			}
			else if (current_debt<best_debt_this_rep) {
				swapped_trees=trees_this_rep=0;
				best_debt_this_rep=current_debt;
				cout << "\t Better cladogram of length " << best_debt_this_rep << "\n";	// SUPPRESS ME!
				abort_swap=1;						// do not continue swapping on the suboptimal tree
			}
			if (current_debt==best_debt_this_rep) {
				if (trees_this_rep<ChannelSize) {
					++trees_this_rep;
					AddTreeToCube(tree, tree_cube[nopttrees+trees_this_rep], ntaxa);
					tree_cube[nopttrees+trees_this_rep][0][0]=true;
					CopyTree(tree, tree_list[trees_this_rep], ntaxa);
				}
			}
		}
	}
}

void HeuristicSearchClass::beginAncSearch (DATA *data, int &nopttrees) {
	int i, j;
	float current_debt=0;
	float a;
	NODE *root;
	NODE *hold_node=0;		//
	NODE *hold_node2=0;		//

	if (SingleOutgroup) {
		AddTreeToCube(tree, this_cube, ntaxa);
		for (i=1;i<=nopttrees;++i) {				// compare with previously discovered optimal trees
			if (tree_cube[i][0][0] && CompareCubes(this_cube, tree_cube[i], ntaxa)) { // if tree was identical with a previously discovered topology
				abort_rep=true; 	// abort this rep
				abort_swap=true;	// abort swapping on this tree
				for (j=1;j<=nislands;++j) if (island[j].first_tree<=i && island[j].last_tree>=i) ++island[j].times_hit;
				cout << "\t*** Topology identical to a previous optimal tree\n";
				return;
			}
		}			
		for (i=1;i<=trees_this_rep;++i) if (tree_cube[nopttrees+i][0][0] && CompareCubes(this_cube, tree_cube[nopttrees+i], ntaxa)) 
			return; // if tree was identical with a previously discovered one from this rep
		root=GetRoot(tree, ntaxa);						// need root from this new topology
//		parsimony.DateTree(root, data);					// need to date the nodes on this new topology
		current_debt=parsimony.TotalDebt(root, data);	// this calls DateTree which is needed for MaxStratSavings
		a=MaxStratSavings(tree, data, &AncSearchReps, MaxSavingsMethod);
//			cout << ("\nCurrent Debt: %.1f\tMSS: %.1f\tBestDebt: %.1f\n", current_debt, a, best_debt_this_rep);
//			cin >> gobots;
		if ((current_debt-a)<=best_debt_this_rep) { // if the current debt minus the max strat debt savings is less than the best overall debt (i.e. this topology COULD be shorter)
			added_topology=false;
			if (!AncSearchType) FullAncSearch(tree, root, 1, data, nopttrees);
			else {
				for (i=1;i<=AncSearchReps;++i) {
					
					if (!random) PutSimpleAdditionSequence(add_seq, ntaxa);					// sequential order
					else PutRandomAdditionSequence(add_seq, ntaxa);							// random order
					HeuristicAncSearch(tree, root, 1, data, nopttrees);
				}
			}
		}
	}
	else for (j=1;j<2**ntaxa;++j) {
		root=GetRoot(tree, ntaxa);
		if (tree[j].anc!=root) {								// these two are redundant with the single rerooting on the current root 
			CopyTree(tree, tree4, ntaxa);						// copy tree to the dummy tree
			Reroot(tree4, &tree4[j], ntaxa);					// reroot dummy tree on taxon j
//			PrintTree(tree4, ntaxa);
			root=GetRoot(tree4, ntaxa);
			if (OutgroupSpecified) MakeOutgroupFromTree(root);	// if the outgroup has been specified, check for outgroup nodes
			hold_node=root->left;
			hold_node2=root->right;
			if (!OutgroupSpecified ||
			    (OutgroupSpecified && (GetOutgroup(hold_node->label) || 
			    					   GetOutgroup(hold_node2->label)))) {
				AddTreeToCube(tree4, this_cube, ntaxa);
				for (i=1;i<=nopttrees;++i) {			// compare with previously discovered optimal trees
					if (tree_cube[i][0][0] && CompareCubes(this_cube, tree_cube[i], ntaxa)) { // if tree was identical with a previously discovered one
						abort_rep=true; 	// abort this rep
						abort_swap=true;	// abort swapping on this tree
						for (j=1;j<=nislands;++j) if (island[j].first_tree<=i && island[j].last_tree>=i) ++island[j].times_hit;
						cout << "\t*** Topology identical to a previous optimal tree\n";
						return;
					}
				}
				for (i=1;i<=trees_this_rep;++i) 
					if (tree_cube[nopttrees+i][0][0] && CompareCubes(this_cube, tree_cube[nopttrees+i], ntaxa)) return; // if tree was identical with a previously discovered one from this rep
//				parsimony.DateTree(root, data);						// get dates for new rooting
				current_debt=parsimony.TotalDebt(root, data);
				a=MaxStratSavings(tree4, data, &AncSearchReps, MaxSavingsMethod);
				if (current_debt<=(best_debt_this_rep+a)) {
					added_topology=false;
					if (!AncSearchType) FullAncSearch(tree4, root, 1, data, nopttrees);
					else {
						for (i=1;i<=AncSearchReps;++i)
							{
							if (!random) PutSimpleAdditionSequence(add_seq, ntaxa);					// sequential order
							else PutRandomAdditionSequence(add_seq, ntaxa);							// random order
							if (i>1) UnfixTree(tree4, ntaxa);	// if there are multiple anc searches on the same tree, it must be "cleared" before performing again
							HeuristicAncSearch(tree4, root, 1, data, nopttrees);
						}
					}
				}
			}
		}
	}
}

void HeuristicSearchClass::FullAncSearch(NODE * const in_tree, NODE * const root, int taxon, DATA * const data, int &nopttrees) {
	NODE *s=0;
	bool go;
	int i;
	
	s=in_tree[taxon].sister;
	for (i=0;i<=1;++i) {
		go=0;
		if (!i && !s->isAnc) {
			in_tree[taxon].Fix();
			go=1;
		}
		else if (i) {
			in_tree[taxon].Unfix();
			go=1;
		}
		if (go) {
			if (taxon<*ntaxa) FullAncSearch(in_tree, root, taxon+1, data, nopttrees);
			else FinishAncSearch(in_tree, root, data, nopttrees);
		}
	}
}

void HeuristicSearchClass::HeuristicAncSearch(NODE * const in_tree, NODE * const root, int taxon, DATA * const data, int &nopttrees) {
	NODE *s=0;
	float savings=-1.0, sister_savings=-1.0;
	int i;
	bool go;
//	int gobots;
	
	parsimony.MakeFullFS(in_tree, data); 						// begin with a "full" matrix with reconstructed ancestral states
//	PrintTree(in_tree, data->GetNTaxa());
//	data->PrintFullMatrix();
//	cin>>gobots;

	s=in_tree[add_seq[taxon]].sister;

	if (!s->isAnc) {
		float morph_debt=-1.0;
		QuickMorphDebt(&in_tree[add_seq[taxon]], data, &morph_debt);// get the morph debt that will be incurred if this taxon is fixed
		if (*data->GetFInt(add_seq[taxon])>=0) QuickStratSavings(&in_tree[add_seq[taxon]], data, &savings);	// get the strat savings if this taxon is fixed and not missing strat data
		else savings=0;
		if (savings<morph_debt) savings=-1;						// if the savings isn't better than the debt
		else {
			savings-=morph_debt;								// subtract incurred morph debt from stratigraphic savings
			if (s->tip)	{		 								// if this node's sister is a tip
				QuickMorphDebt(s, data, &morph_debt);			// get the morph debt if its sister is fixed
				if (*data->GetFInt(s)>=0) QuickStratSavings(s, data, &sister_savings);	// get the strat savings of its sister if fixed and not missing strat data
				else sister_savings=0;
				sister_savings-=morph_debt;						// sister's savings is the difference
			}
		}
	}

	for (i=0;i<=1;++i)	{	// first time through (i=0), try the taxon fixed, then (i=1) unfixed
		go=false;
		if (!i && savings>=0 && savings>=sister_savings) {	// if savings is positive AND greater than or equal to its sister
			in_tree[add_seq[taxon]].Fix();					// fix the taxon
			go=true;										// and continue checking other taxa
		}
		else if (i && (savings<=0 || (savings>0 && savings<=sister_savings))) {	// if savings is less than or equal to zero, or it's greater than zero, but less than or equal to it's sister's savings
			in_tree[add_seq[taxon]].Unfix();			// unfix the taxon
			go=true;										// and continue checking taxa
		}
 		if (go) {
			if (taxon<*ntaxa) HeuristicAncSearch(in_tree, root, taxon+1, data, nopttrees);
			else FinishAncSearch(in_tree, root, data, nopttrees);
		}
	}
	in_tree[add_seq[taxon]].Unfix();						// unfix the taxon
}

void HeuristicSearchClass::FinishAncSearch(NODE * const in_tree, NODE * const root, DATA * const data, int &nopttrees) {
	float debt;
	int i;
//	int gobots;
	
	debt=parsimony.TotalDebt(root, data);			// get the debt
	if (debt<=best_debt_this_rep) {
		if (debt<best_debt) {
			nopttrees=0;
			best_debt=debt;
		}
		if (debt<best_debt_this_rep) {
			swapped_trees=trees_this_rep=0;
			best_debt_this_rep=debt;
			cout << "\t Better strat tree of length " << best_debt_this_rep << "\n";
			abort_swap=true;						// do not continue swapping on the suboptimal tree
		}
		else if (AncSearchReps>1 && debt==best_debt_this_rep) {
			for (i=trees_this_rep;i>0;--i) { // new tree is most likely to match recently added trees 
				if (tree_cube[nopttrees+i][0][0] && CompareAncs(in_tree, tree_list[i], ntaxa) && CompareCubes(this_cube, tree_cube[nopttrees+i], ntaxa)) return;// if tree was identical with a previously discovered one from this rep
			}
		}
		if (debt==best_debt_this_rep && trees_this_rep<ChannelSize) {
			++trees_this_rep;
			LadderizeRight(in_tree, ntaxa);
			AddTreeToCube(in_tree, tree_cube[nopttrees+trees_this_rep], ntaxa);
			if (!added_topology) {	// this is to prevent branch swapping on trees with same topology, but different ancestorizations - only one topology needs ot be swapped
				tree_cube[nopttrees+trees_this_rep][0][0]=true;	// sets the first tree to be the topology to be swapped
				added_topology=true;	// prevents subsequent trees from being designated as such
			}
			else tree_cube[nopttrees+trees_this_rep][0][0]=false;	// tree will not be branch swapped
			CopyTree(in_tree, tree_list[trees_this_rep], ntaxa);
		}
	}
}

float HeuristicSearchClass::MaxStratSavings(NODE *start_tree, DATA *data, int *reps, bool method)
{
	NODE *root=0;
	int i;
	float overall_max_savings;
	float max_savings;
	
	root=GetRoot(start_tree, data->GetNTaxa());
	overall_max_savings=0.0;

	for (i=1;i<=*reps;++i) {
		max_savings=0.0;
		if (!random) PutSimpleAdditionSequence(add_seq, ntaxa);
		else PutRandomAdditionSequence(add_seq, data->GetNTaxa());							// random order
		if (!method) DoMinStratDebt(start_tree, root, 1, data, max_savings);
		else {
			parsimony.MakeFullFS(start_tree, data);
			DoMinTotalDebt(start_tree, root, 1, data, max_savings, add_seq);
		}
		if (max_savings>overall_max_savings) overall_max_savings=max_savings;
		UnfixTree(start_tree, data->GetNTaxa());
	}
	return overall_max_savings;
}

void HeuristicSearchClass::DoMinStratDebt(NODE *in_tree, NODE *root, int taxon, DATA * const data, float &max_savings) {
	NODE *s=0;
	float savings;
	
	s=in_tree[add_seq[taxon]].sister;

	if (*data->GetFInt(add_seq[taxon])<*data->GetFInt(s) && 	// this taxon's FInt is before its sister's (if they're the same, no strat debt can be saved)
		!s->isAnc && 											// and its sister isn't fixed
		*data->GetFInt(add_seq[taxon])>=0)	 					// and strat data is not missing
			QuickStratSavings(&in_tree[add_seq[taxon]], data, &savings);
	else savings=0;												// otherwise, it can't save any debt
	if (savings>0) {											// this taxon saved some strat debt
		max_savings+=savings;									// add savings to the max savings
		in_tree[add_seq[taxon]].Fix();							// and fix this taxon
	}
	if (taxon<*data->GetNTaxa()) DoMinStratDebt(in_tree, root, taxon+1, data, max_savings);
}

void HeuristicSearchClass::DoMinTotalDebt(NODE *in_tree, NODE *root, int taxon, DATA * const data, float &max_savings, int * const add_seq) {
	NODE *s=0;
	float savings=-1, sister_savings=-1;

	s=in_tree[add_seq[taxon]].sister;

	if (*data->GetFInt(add_seq[taxon])<=*data->GetFInt(s->label) && !s->isAnc) {
		float morph_debt;
		QuickMorphDebt(&in_tree[add_seq[taxon]], data, &morph_debt);
		if (*data->GetFInt(&in_tree[add_seq[taxon]])>=0) QuickStratSavings(&in_tree[add_seq[taxon]], data, &savings);
		else savings=0;
		if (savings<morph_debt) savings=-1;
		else {
			savings-=morph_debt;
			if (s->tip) { 					// if this node's sister is a tip, 
				QuickMorphDebt(s, data, &morph_debt);
				if (*data->GetFInt(s)>=0) QuickStratSavings(s, data, &sister_savings);
				else sister_savings=0;
				sister_savings-=morph_debt;
			}
		}
	}
	if (savings>0 && (savings>sister_savings || (savings==sister_savings && add_seq[taxon]<s->label))) {	
		max_savings+=savings;
		in_tree[add_seq[taxon]].Fix();
		parsimony.MakeFullFS(in_tree, data);
	}
	if (taxon<*data->GetNTaxa()) DoMinTotalDebt(in_tree, root, taxon+1, data, max_savings, add_seq);
}

void HeuristicSearchClass::QuickMorphDebt(NODE * const node, DATA * const data, float *morph_debt) {
	NODE *a;
	int j;
	
	*morph_debt=0;
	a=node->anc;
	
	for (j=1;j<=*data->GetNChrs();++j) {
		if (*data->GetMatrix(node->label, j) &&	// if the data is not missing
		    data->GetInc(j) && 					// if the data is included
		    *data->GetNStates(j)>1) {			// if the data has more than 1 state

			if( !(*data->GetMatrix(node->label, j) & *data->GetMatrix(a->label, j)) ) {
				if ( *data->GetChrType(j)==0) *morph_debt+=*data->GetChrWeight(j);	// if the taxon and its anc share no states in common, increase the debt
				else if (*data->GetChrType(j)==1) { 
					int max_state, min_state;
					if (*data->GetMatrix(node->label, j) > *data->GetMatrix(a->label, j)) {
						max_state=data->MaxState(a->label, j);
						min_state=data->MinState(node->label, j);
					}
					else {
						max_state=data->MaxState(node->label, j);
						min_state=data->MinState(a->label, j);
					}
					*morph_debt += (((float)min_state-(float)max_state) * *data->GetChrWeight(j));
				}
			}
		}
	}
}

void HeuristicSearchClass::QuickStratSavings(NODE * const node, DATA * const data, float *savings) {
	NODE *s=0;
	
	s=node->sister;
	*savings=0;
	if (*data->GetFInt(node)<*data->GetFInt(s)) {		// if node's FA is less than its sister's it CAN save strat debt  
		if (*data->GetFInt(s)>*data->GetLInt(node)) 	// if the FInt of the desc is after the LInt of the fixed anc
			 *savings=*data->GetChrWeight(data->GetStratChr())*(((float)*data->GetLInt(node)-(float)*data->GetFInt(node))+1);	// debt saved with be range of f (rest of ghost lineage will persist
		else *savings=*data->GetChrWeight(data->GetStratChr())*((float)*data->GetFInt(s)-(float)*data->GetFInt(node));		// this is the debt caused by the ghost lineage to d
	}
}

void HeuristicSearchClass::AllocateOutgroup() {
	outgroup = new bool [2**ntaxa];
	for (int i=0;i<2**ntaxa;++i) outgroup[i]=false;
}
void HeuristicSearchClass::SetOutgroup (int taxon, bool i) { outgroup[taxon]=i; }
bool HeuristicSearchClass::GetOutgroup(int taxon) const { return outgroup[taxon]; }

void HeuristicSearchClass::MakeOutgroupFromTree(NODE *node) {
	NODE *l, *r;

	l=node->left;	
	if (!l->tip) MakeOutgroupFromTree(l);
	r=node->right;	
	if (!r->tip) MakeOutgroupFromTree(r);
	
	if (outgroup[l->label] && outgroup[r->label]) outgroup[node->label]=true;
	else outgroup[node->label]=false;
}

void HeuristicSearchClass::initializeOutgroup(DATA *data) {
	int i, chr_number;

	root_taxon=1;
	for (i=1;i<=*ntaxa;++i) {
		if (GetOutgroup(i)) {
			if (!SingleOutgroup) {
				SingleOutgroup=true;
				root_taxon=i;
			}
			else {
				SingleOutgroup=false;
				root_taxon=1;
				break;
			}
		}
	}

	if (data->hasStratChr && !SingleOutgroup) {
		chr_number=data->GetStratChr();
		if (OutgroupSpecified) {
			for (i=1;i<=*ntaxa;++i) {
				if (GetOutgroup(i)) {
					if (data->MinState(i, chr_number)<data->MinState(root_taxon, chr_number)) root_taxon=i; // find oldest outgroup taxon to be arbitrary root
				}
			}
		}
		else {																// if no outgroup taxa
			for (i=1;i<=*ntaxa;++i) if (data->MinState(i, chr_number)<data->MinState(root_taxon, chr_number)) root_taxon=i;	// find oldest ingroup taxon to be arbitrary root
		}
	}
	else if (!data->hasStratChr && !SingleOutgroup && OutgroupSpecified) {
		for (i=1;i<=*ntaxa;++i) {
			if (GetOutgroup(i)) {
				root_taxon=i;
				break;
			}
		}
	}
	
	if (OutgroupSpecified) {
		cout << "\tOutgroup is:\t";	// SUPPRESS ME!
		for (i=1;i<=*ntaxa;++i) if (GetOutgroup(i)) cout << data->taxon_label[i] <<" ";
		cout << "\n";
	}
	else cout << "\tNo outgroup specified\n";	// SUPPRESS ME!

}

void HeuristicSearchClass::PrintIslandData() {
	int i;
	
/*	cout << "\n";
	cout<<"\tIsland\tDebt\tNTrees(NTops)\tFirst-Last\tTimes Hit\n";
	for (i=1;i<=nislands;++i) {
		cout.width(10);
		cout<<i;
		cout.width(6);
		cout<<island[i].debt;
		cout.width(14);
		cout<<island[i].ntrees<<"("<<island[i].ntopologies<<")";
		cout.width(11);
		if (island[i].debt==best_debt && island[i].ntrees==1) cout<<island[i].first_tree;
		else if (island[i].debt==best_debt && island[i].ntrees>1) cout<<island[i].first_tree<<"-"<<island[i].last_tree;
		else cout << "---";
		cout.width(10);
		cout <<island[i].times_hit<<"\n";
	}
*/
}

/*
NODE *HeuristicSearchClass::Bootstrap(DATA * const data, int *boot_reps, int *boot_MaxTrees)
{
	NODE ***boot_list;
	NODE *consensus;
	int i;
	
	boot_list= new NODE ** [*boot_reps+1];
	for (i=1;i<=*boot_reps;++i) boot_list[i]=AllocateTreeList(boot_MaxTrees, data->GetNTaxa());	// each row is a list of trees from each rep
	for (i=1;i<=*boot_reps;++i) FreeTreeList(boot_list[i], data->GetNTaxa());
	delete [] boot_list;
	
	return consensus;
}

*/
/*
void HeuristicSearchClass::RegraftSubtree(NODE * const current_tree, int * const pruned_node, int * const target_node, DATA * const data, int &nopttrees)
{
//	NODE *tree=0;
	NODE *root=0;
	NODE *hold_node=0;		//
	NODE *cut_node=0;		// node that is removed when subtree is pruned
	NODE *drop_node=0;		// node that replaces cut node (pruned_node's sister)
	float current_debt=0;
	int cut_node_label=0;
	int i, j;
	bool lights;
	int gobots;
	
	lights=false;
		
//	tree=GetInitializedTree(ntaxa);				// dummy tree whose branches will be swapped
	CopyTree(current_tree, tree, ntaxa);		// copy tree to be swapped - original tree is untouched

	if (lights)
		{
//		cout << "Original Tree\n";
//		PrintTree(current_tree, ntaxa);
//		cout << "Copied Tree\n";
//		PrintTree(tree, ntaxa);
		cout << "\tSwapping Node " << tree[*pruned_node].label << " below Node " << tree[*target_node].label << "\t";
		}

  // Get Cut_Node
	cut_node=tree[*pruned_node].anc;	// node that is cut out is the ancestor of the pruned_node
	cut_node_label=cut_node->label;						// keep track of it's label because we're going to reuse it for the new node
	if (lights) cout << "Cut_Node = " << cut_node->label << "\n";

  // Get Drop_Node
  	drop_node=tree[*pruned_node].sister;
	if (lights) cout << "Drop_Node = " << drop_node->label << "\n";

  // Get Hold_Node
	if (cut_node->anc!=0)								// as long as the pruned_node's anc is not ROOT
		{
		hold_node=cut_node->anc;						// 	cut_node's ancestor will receive cut_node's other descendant (drop_node) 
		if (hold_node->left==cut_node) hold_node->left=drop_node;
		else if (hold_node->right==cut_node) hold_node->right=drop_node;
		else
			{
			cout << "\n\n***\tCan't Obtain Hold_Node\n";
			cin >> gobots;
			}
		if (lights) cout << "Hold_Node = " << hold_node->label << "\n";
		drop_node->anc=hold_node;					// move the cut node's ancestor to the drop node
		}
	else 												// however, if it is ROOT, 
		{
		drop_node->anc=0;								// just cut it out, and make the drop_node the root
		root=drop_node;
		if (lights) cout << "Hold_Node = ROOT\n";
		}

  // Place Subtree on Requested Branch
	tree[cut_node_label].right=&tree[*target_node];	// new node's right branch is the target node
	tree[cut_node_label].left=&tree[*pruned_node];	// new node's left branch is the clipped subtree
	tree[cut_node_label].anc=tree[*target_node].anc;// new node's ancestor is the target node's OLD anc
	tree[*pruned_node].anc=&tree[cut_node_label];	// subtree's new anc is the new node
	hold_node=tree[*target_node].anc;				// hold node is the target_node's OLD anc for now
	if (hold_node->left==&tree[*target_node]) hold_node->left=&tree[cut_node_label];			// target node's old anc's left is the new node
	else if (hold_node->right==&tree[*target_node]) hold_node->right=&tree[cut_node_label];	// target node's old anc's right is the new node
	else {
		cout << "\n\n***\tRegraft Failed\n";
		cin >> gobots;
		}
	tree[*target_node].anc=&tree[cut_node_label];	// target node's new anc is also the new node

	if (lights)
		{
//		cout << "Tree was changed to:\n";
//		PrintTree(tree, ntaxa);
//		cin >> gobots;
		cout << "Tree " << swapped_trees << " is still:\n";
		PrintTree(current_tree, ntaxa);
//		cin >> gobots;
		}

//	++rearrangements;
	Reroot(tree, &tree[root_taxon], ntaxa);						// reroot new topology on the outgroup
	if (lights) {
		cout << "Tree after rerooting on " << root_taxon << ":\n";
		PrintTree(tree, ntaxa);
		cin >> gobots;
	}

	if (data->hasStratChr) beginAncSearch (data, nopttrees);
	else {							// if there is no strat character
		root=GetRoot(tree, ntaxa);									// need root from this new topology
		current_debt=parsimony.TotalDebt(root, data);
		if (current_debt<=best_debt_this_rep) {
			AddTreeToCube(tree, this_cube, ntaxa);
			if (current_debt==best_debt) {
				for (i=1;i<=nopttrees;++i) {							// compare with previously discovered optimal trees
					if (CompareCubes(this_cube, tree_cube[i], ntaxa)) {	// if tree was identical with a previously discovered one
						abort_rep=true; 	// abort this rep
						abort_swap=true;	// abort swapping on this tree
						for (j=1;j<=nislands;++j) if (island[j].first_tree<=i && island[j].last_tree>=i) ++island[j].times_hit;
						cout << "\t*** Found a previous optimal cladogram\n";
						return;
					}
				}
			}
			else if (current_debt<best_debt) {
				nopttrees=0;
				best_debt=current_debt;
			}
			if (current_debt==best_debt_this_rep) {
				for (i=1;i<=trees_this_rep;++i) 
					if (CompareCubes(this_cube, tree_cube[nopttrees+i], ntaxa)) return;// if tree was identical with a previously discovered one from this rep
			}
			else if (current_debt<best_debt_this_rep) {
				swapped_trees=trees_this_rep=0;
				best_debt_this_rep=current_debt;
				cout << "\t Better cladogram of length " << best_debt_this_rep << "\n";
				abort_swap=1;						// do not continue swapping on the suboptimal tree
			}
			if (current_debt==best_debt_this_rep) {
				if (trees_this_rep<MaxTrees-1) {
					++trees_this_rep;
					AddTreeToCube(tree, tree_cube[nopttrees+trees_this_rep], ntaxa);
					tree_cube[nopttrees+trees_this_rep][0][0]=true;
					CopyTree(tree, tree_list[trees_this_rep], ntaxa);
				}
	//			else cout << "\n\t***\tMAXTREES met\n";
			}
		}
	}
}
*/

/*
void HeuristicSearchClass::HeuristicAncSearch(NODE * const in_tree, NODE * const root, int taxon, DATA * const data, int &nopttrees)
{
	NODE *s=0;
	float savings=-1.0, sister_savings=-1.0, morph_debt=-1.0;
	int i;
	bool go;
	int gobots;
	

	parsimony.MakeFullFS(in_tree, data); 						// begin with a "full" matrix with reconstructed ancestral states
//	PrintTree(in_tree, data->GetNTaxa());
//	data->PrintFullMatrix();
//	cin>>gobots;

	s=in_tree[add_seq[taxon]].sister;

	if (!s->isAnc) {
		QuickMorphDebt(&in_tree[add_seq[taxon]], data, &morph_debt);// get the morph debt that will be incurred if this taxon is fixed
		if (*data->GetFInt(add_seq[taxon])>=0) QuickStratSavings(&in_tree[add_seq[taxon]], data, &savings);	// get the strat savings if this taxon is fixed and not missing strat data
		else savings=0;
		if (savings<morph_debt) savings=-1;						// if the savings isn't better than the debt
		else {
			savings-=morph_debt;								// subtract incurred morph debt from stratigraphic savings
			if (s->tip)	{		 								// if this node's sister is a tip
				QuickMorphDebt(s, data, &morph_debt);			// get the morph debt if its sister is fixed
				if (*data->GetFInt(s)>=0) QuickStratSavings(s, data, &sister_savings);	// get the strat savings of its sister if fixed and not missing strat data
				else sister_savings=0;
				sister_savings-=morph_debt;						// sister's savings is the difference
			}
		}
	}

	for (i=0;i<=1;++i)	{	// first time through (i=0), try the taxon fixed, then (i=1) unfixed
		go=false;
		if (!i && savings>=0 && savings>=sister_savings) {	// if savings is positive AND greater than or equal to its sister
			in_tree[add_seq[taxon]].Fix();					// fix the taxon
			go=true;										// and continue checking other taxa
		}
		else if (i && (savings<=0 || (savings>0 && savings<=sister_savings))) {	// if savings is less than or equal to zero, or it's greater than zero, but less than or equal to it's sister's savings
			in_tree[add_seq[taxon]].Unfix();			// unfix the taxon
			go=true;										// and continue checking taxa
		}
 		if (go) {
			if (taxon<*ntaxa) HeuristicAncSearch(in_tree, root, taxon+1, data, nopttrees);
			else FinishAncSearch(in_tree, root, data, nopttrees);
		}
	}
	in_tree[add_seq[taxon]].Unfix();						// unfix the taxon
//	if (tried_fixed) parsimony.MakeFullFS(in_tree, data);	// if it tried taxon fixed, but not unfixed, it won't FS again so we need this here
}
*/
