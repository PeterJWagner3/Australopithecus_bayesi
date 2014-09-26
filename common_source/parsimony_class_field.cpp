#include "parsimony_class_field.h"

using namespace std;

ParsimonyClass::ParsimonyClass() { matrix=0; }
ParsimonyClass::~ParsimonyClass() {};

float ParsimonyClass::TotalDebtPolymorph(NODE *root, DATA *data) {
	if (!matrix) matrix=data->PassMatrix();
	debt=0;
	for (chr_number=1;chr_number<=*data->GetNChrs();++chr_number) {
		if (*data->GetNStates(chr_number)>1 && data->GetInc(chr_number)) { 	/*&& *data->GetPattern(chr_number)>0*/ 
			if (*data->GetChrType(chr_number)==0) DownPassUnord(root, data);
			else if (*data->GetChrType(chr_number)==1) DownPassOrd(root, data);				
			else if (*data->GetChrType(chr_number)==4) {
				DateTree(root, data);
				DoStratDebtInt(root, data);
			}
			else DownPassUnord(root, data);
		}
	}
	return debt;
}

float ParsimonyClass::TotalDebt(NODE *root, DATA *data) {
	if (!matrix) matrix=data->PassMatrix();
	debt=0;
	for (chr_number=1;chr_number<=*data->GetNChrs();++chr_number) {
		if (*data->GetNStates(chr_number)>1 && data->GetInc(chr_number)) { 	/*&& *data->GetPattern(chr_number)>0*/ 
			if (*data->GetChrType(chr_number)==0) {							// unordered chrs
				if (!data->GetPolymorph(0, chr_number)) DownPassUnord(root, data);
				else {
					DownPassUnordPoly(root, data);
					for (int i=*data->GetNTaxa()+1;i<2**data->GetNTaxa();++i) data->SetPolymorph(i, chr_number, false);	// reset polymorphic status of all internal nodes
				}
			}
			else if (*data->GetChrType(chr_number)==1) {					// ordered chrs
				if (!data->GetPolymorph(0, chr_number)) DownPassOrd(root, data);
				else {
					DownPassOrdPoly(root, data);
					for (int i=*data->GetNTaxa()+1;i<2**data->GetNTaxa();++i) data->SetPolymorph(i, chr_number, false);	// reset polymorphic status of all internal nodes
				}
			}
			else if (*data->GetChrType(chr_number)==4) {
				DateTree(root, data);
				DoStratDebtInt(root, data);
			}
			else {												// unsupported chr types treated as unordered chrs.
				if (!data->GetPolymorph(0, chr_number)) DownPassUnord(root, data);
				else {
					DownPassUnordPoly(root, data);
					for (int i=*data->GetNTaxa()+1;i<2**data->GetNTaxa();++i) data->SetPolymorph(i, chr_number, false);	// reset polymorphic status of all internal nodes
				}
			}
		}
	}
	return debt;
}

float ParsimonyClass::SingleChrDebt(NODE * const tree, DATA * const data, int *chr) {
	NODE *root=0;
	
	if (!matrix) matrix=data->PassMatrix();
	debt=0;
	root=GetRoot(tree, data->GetNTaxa());
	chr_number=*chr;
	if (*data->GetChrType(chr_number)==0) {
		if (!data->GetPolymorph(0, chr_number)) DownPassUnord(root, data);
		else DownPassUnordPoly(root, data);
		for (int i=*data->GetNTaxa()+1;i<2**data->GetNTaxa();++i) data->SetPolymorph(i, chr_number, false);	// reset polymorphic status of all internal nodes
	}
	else if (*data->GetChrType(chr_number)==1) {
		if (!data->GetPolymorph(0, chr_number)) DownPassOrd(root, data);
		else {
			DownPassOrdPoly(root, data);
			for (int i=*data->GetNTaxa()+1;i<2**data->GetNTaxa();++i) data->SetPolymorph(i, chr_number, false);	// reset polymorphic status of all internal nodes		
		}
	}
	else if (*data->GetChrType(chr_number)==4) GetStratDebt (root, data);
	return debt;
}

void ParsimonyClass::GetDownPassSet(NODE * const tree, DATA * const data, int *chr) {
	NODE *root=0;
	chr_number=*chr;
	root=GetRoot(tree, data->GetNTaxa());
	DownPassUnord(root, data);
}

void ParsimonyClass::GetFinalSet(NODE * const root, DATA * const data, int *chr) {
//	int gobots;
	
	chr_number=*chr;
	if (*data->GetChrType(chr_number)==0) {
		if (!data->GetPolymorph(0,chr_number)) DownPassUnord(root, data);
		else {	
			DownPassUnordPoly(root, data);
//			for (int i=*data->GetNTaxa()+1;i<2**data->GetNTaxa();++i) data->SetPolymorph(i, chr_number, false);	// reset polymorphic status of all internal nodes
		}
		FinalPassUnord(root, data);
	}
	else if (*data->GetChrType(chr_number)==1) {
		DownPassOrd(root, data);
//		data->PrintFullMatrix();
		FinalPassOrd(root, data);
//		data->PrintFullMatrix();
//		cin>>gobots;
	}
}

void ParsimonyClass::MakeFullFS(NODE * const in_tree, DATA * const data) {
	NODE *root=0;

	if (!matrix) matrix=data->PassMatrix();
	root=GetRoot(in_tree, data->GetNTaxa());
	for (int j=1;j<=*data->GetNChrs();++j) if (*data->GetChrType(j)<4) GetFinalSet(root, data, &j);
}

void ParsimonyClass::DownPassUnord(NODE * const a, DATA * const data) {
	NODE *l=0, *r=0, *f=0, *d=0;

	l=a->left;
	r=a->right;
	if (!l->tip) DownPassUnord(l, data);
	if (!r->tip) DownPassUnord(r, data);
	if (l->isAnc) { f=l; d=r; }
	else if (r->isAnc) { f=r; d=l; }

	if (!f) {		// if this node is not fixed as an ancestor
		if (matrix[l->label][chr_number] && matrix[r->label][chr_number]) {
			matrix[a->label][chr_number] = matrix[l->label][chr_number] & matrix[r->label][chr_number];	// anc is intersection of desc's
			if (!matrix[a->label][chr_number]) {	// if there were no identical states
				debt+=*data->GetChrWeight(chr_number);	// increase the debt
				matrix[a->label][chr_number] = matrix[l->label][chr_number] | matrix[r->label][chr_number];
			}
		}
		else if (matrix[l->label][chr_number] && !matrix[r->label][chr_number]) matrix[a->label][chr_number] = matrix[l->label][chr_number];
		else matrix[a->label][chr_number] = matrix[r->label][chr_number];
	}
	else {			// if this node IS fixed as an ancestor
		if (matrix[f->label][chr_number]) {
			matrix[a->label][chr_number] = matrix[f->label][chr_number];
			if (matrix[d->label][chr_number] && !(matrix[f->label][chr_number] & matrix[d->label][chr_number])) debt+=*data->GetChrWeight(chr_number);
		}
		else matrix[a->label][chr_number] = matrix[d->label][chr_number];
	}
}

void ParsimonyClass::DownPassOrd(NODE * const a, DATA * const data) {
	NODE *l=0, *r=0, *f=0, *d=0;

	l=a->left;
	r=a->right;
	if (!l->tip) DownPassOrd(l, data);
	if (!r->tip) DownPassOrd(r, data);
	if (l->isAnc) { f=l; d=r; }
	else if (r->isAnc) { f=r; d=l; }

	if (!f) {		// if this node is not fixed as an ancestor
		if (matrix[l->label][chr_number] && matrix[r->label][chr_number]) matrix[a->label][chr_number] = OrdCombination(l, r, true, data);
		else if (matrix[l->label][chr_number] && !matrix[r->label][chr_number]) matrix[a->label][chr_number] = matrix[l->label][chr_number];	// only right is missing
		else matrix[a->label][chr_number] = matrix[r->label][chr_number];																		// only left is missing, or both are
	}
	else {			// if this node IS fixed as an ancestor
		if (matrix[f->label][chr_number]) {	// the fixed taxon is not missing data
			if (!data->isUncertain(f->label, chr_number) || 
				(data->isUncertain(f->label, chr_number) && !(matrix[f->label][chr_number] & matrix[d->label][chr_number])))
				 	matrix[a->label][chr_number] = matrix[f->label][chr_number];	// the ancestral states are that of the fixed taxon
			else matrix[a->label][chr_number] = matrix[f->label][chr_number] & matrix[d->label][chr_number];
			if (matrix[d->label][chr_number]) { // if the descendant is not missing data (i.e., both tips have data)
				OrdCombination(f, d, true, data);	// this will count up the debt, but not assign the return value (set of states) to anything
			}
		}
		else matrix[a->label][chr_number] = matrix[d->label][chr_number];	// the descendant's DPS, which can be empty (missing data)
	}
}

void ParsimonyClass::DownPassUnordPoly(NODE * const a, DATA * const data) {
	NODE *l=0, *r=0, *f=0, *d=0;
	int k;

	l=a->left;
	r=a->right;
	if (!l->tip) DownPassUnordPoly(l, data);
	if (!r->tip) DownPassUnordPoly(r, data);
	if (l->isAnc) { f=l; d=r; }
	else if (r->isAnc) { f=r; d=l; }

	data->SetPolymorph(a->label, chr_number,false);
	if (!f) { 		// if this node is not fixed as an ancestor
		if (matrix[l->label][chr_number] && matrix[r->label][chr_number]) {
			matrix[a->label][chr_number] = matrix[l->label][chr_number] & matrix[r->label][chr_number];	// anc has intersection of desc's states
			if (!data->GetPolymorph(l->label, chr_number) &&	// if left is not polymorphic and
				!data->GetPolymorph(r->label, chr_number)) {	// right is not polymorphic
					if (!matrix[a->label][chr_number]) debt+=*data->GetChrWeight(chr_number);	// with non-overlapping states, there can only be one step
			}
			else {	// at least one descendant is polymorphic
				if (data->GetPolymorph(l->label, chr_number)) debt+=((float)data->getNumStates(l->label, chr_number))**data->GetChrWeight(chr_number);
				if (data->GetPolymorph(r->label, chr_number)) debt+=((float)data->getNumStates(r->label, chr_number))**data->GetChrWeight(chr_number);
				if (data->GetPolymorph(l->label, chr_number) && data->GetPolymorph(r->label, chr_number)) {
					if (matrix[a->label][chr_number]) {
						int shared_states=0;
						for (k=0;k<*data->GetNStates(chr_number);++k) {
							if ((matrix[l->label][chr_number] & matrix[r->label][chr_number]) & (unsigned long long) pow(2,k)) ++shared_states;
						}
						if (shared_states) {	
							debt-=((2*(float)shared_states)**data->GetChrWeight(chr_number));	// it has already counted the shared states twice
							if (shared_states>1) data->SetPolymorph(a->label, chr_number, true);// if both desc are polymorphic with the same states, then the anc must be as well
						}
					}
					else debt-=*data->GetChrWeight(chr_number);	// at least one of the observed states is going to be the ancestral state
				}
				else if (matrix[a->label][chr_number]) debt-=*data->GetChrWeight(chr_number);// a polymorphic taxon shares a single state with a nonpolymorphic descendant if their states overlap
			}
			if (!matrix[a->label][chr_number]) matrix[a->label][chr_number] = matrix[l->label][chr_number] | matrix[r->label][chr_number];
		}
		else if (matrix[l->label][chr_number] && !matrix[r->label][chr_number]) {
			matrix[a->label][chr_number] = matrix[l->label][chr_number];	// only right is missing
			if  (data->GetPolymorph(l->label,chr_number)) data->SetPolymorph(a->label,chr_number, true);
		}
		else {
			matrix[a->label][chr_number] = matrix[r->label][chr_number];																		// only left is missing, or both are
			if (data->GetPolymorph(r->label,chr_number)) data->SetPolymorph(a->label,chr_number, true);
		}
	}
	else {			// if this node IS fixed as an ancestor
		if (matrix[f->label][chr_number]) {							// if the anc's data is not missing
			if (!data->isUncertain(f->label, chr_number) || 	// anc is either certain or polymorphic (not uncertain), or 
			   (data->isUncertain(f->label, chr_number) && !(matrix[f->label][chr_number] & matrix[d->label][chr_number]))) // anc is uncertain, but has no states overlapping the desc
			   		matrix[a->label][chr_number] = matrix[f->label][chr_number];
			else matrix[a->label][chr_number] = matrix[f->label][chr_number] & matrix[d->label][chr_number];
			if (data->GetPolymorph(f->label, chr_number)) data->SetPolymorph(a->label, chr_number, true);

			if (matrix[d->label][chr_number]) {						// if the desc's data is not missing
				if (!data->GetPolymorph(f->label, chr_number) && !data->GetPolymorph(d->label, chr_number)) {					// if both are certain or uncertain (not  polymorphic) 
					if (!(matrix[f->label][chr_number] & matrix[d->label][chr_number])) debt+=(*data->GetChrWeight(chr_number));	// with non-overlapping states, there can only be one step
				}
				else {	// at least one taxon is polymorphic
					if (data->GetPolymorph(f->label, chr_number)) debt+=((float)data->getNumStates(f->label, chr_number))**data->GetChrWeight(chr_number);
					if (data->GetPolymorph(d->label, chr_number)) debt+=((float)data->getNumStates(d->label, chr_number))**data->GetChrWeight(chr_number);
					if (data->GetPolymorph(f->label, chr_number) && data->GetPolymorph(d->label, chr_number)) {
						if  (matrix[f->label][chr_number] & matrix[d->label][chr_number]) {
							int shared_states=0;
							for (k=0;k<*data->GetNStates(chr_number);++k) {
								if ((matrix[f->label][chr_number] & matrix[d->label][chr_number]) & (unsigned long long) pow(2,k)) ++shared_states;
							}
							if (shared_states) debt-=(2*(float)shared_states)**data->GetChrWeight(chr_number);
						}
						else debt-=*data->GetChrWeight(chr_number);	// at least one of the observed states in f transforms to one of the states in d (don't  need to lose  that state and gain another)
					}
					else if (matrix[f->label][chr_number] & matrix[d->label][chr_number]) debt-=*data->GetChrWeight(chr_number);
				}
			}
		}		
		else {
			matrix[a->label][chr_number] = matrix[d->label][chr_number];																		// only left is missing, or both are
			if  (data->GetPolymorph(d->label,chr_number)) data->SetPolymorph(a->label,chr_number, true);
		}
	}
}

void ParsimonyClass::DownPassOrdPoly(NODE * const a, DATA * const data) {
	NODE *l=0, *r=0, *f=0, *d=0;

	l=a->left;
	r=a->right;
	if (!l->tip) DownPassOrdPoly(l, data);
	if (!r->tip) DownPassOrdPoly(r, data);
	if (l->isAnc) { f=l; d=r; }
	else if (r->isAnc) { f=r; d=l; }

//	if (l->label==9 || r->label==9 || l->label==10 || r->label==10) {
//		cout << "hi";
//	}

	data->SetPolymorph(a->label, chr_number, false);
	if (!f) {		// if this node is not fixed as an ancestor
		if (matrix[l->label][chr_number] && matrix[r->label][chr_number]) {
			if (!data->GetPolymorph(l->label, chr_number) && !data->GetPolymorph(r->label, chr_number)) matrix[a->label][chr_number] = OrdCombination(l, r, true, data);
			else {
				matrix[a->label][chr_number] = OrdCombinationPoly(l, r, true, data);
				if (data->GetPolymorph(l->label, chr_number) && data->GetPolymorph(r->label, chr_number)) {
					bool flag=false;
					for (int i=0;i<=*data->GetNStates(chr_number);++i) {
						if (!flag && ((matrix[l->label][chr_number] & matrix[r->label][chr_number]) & (unsigned long long)pow(2, i))) flag=true;
						else if (flag && ((matrix[l->label][chr_number] & matrix[r->label][chr_number]) & (unsigned long long)pow(2, i))) data->SetPolymorph(a->label, chr_number, true);
					}
				}
			}
		}
		else if (matrix[l->label][chr_number] && !matrix[r->label][chr_number]) {
			matrix[a->label][chr_number] = matrix[l->label][chr_number];	// only right is missing
			if  (data->GetPolymorph(l->label,chr_number)) data->SetPolymorph(a->label,chr_number, true);
		}
		else {
			matrix[a->label][chr_number] = matrix[r->label][chr_number];																		// only left is missing, or both are
			if  (data->GetPolymorph(r->label,chr_number)) data->SetPolymorph(a->label,chr_number, true);
		}
	}
	else {			// if this node IS fixed as an ancestor
		if (matrix[f->label][chr_number]) {	// the fixed taxon is not missing data
			if (!data->isUncertain(f->label, chr_number) || 
				(data->isUncertain(f->label, chr_number) && !(matrix[f->label][chr_number] & matrix[d->label][chr_number])))
				 	matrix[a->label][chr_number] = matrix[f->label][chr_number];	// the ancestral states are those of the fixed taxon
			else matrix[a->label][chr_number] = matrix[f->label][chr_number] & matrix[d->label][chr_number];
			if (data->GetPolymorph(f->label, chr_number)) data->SetPolymorph(a->label, chr_number, true);
				 
			if (matrix[d->label][chr_number]) { // if the descendant is not missing data (i.e., both tips have data)
				if (!data->GetPolymorph(l->label, chr_number) && !data->GetPolymorph(r->label, chr_number)) OrdCombination(f, d, true, data);	// this will count up the debt, but not assign the return value (set of states) to anything
				else 
					OrdCombinationPoly(f, d, true, data);	// this will count up the debt, but not assign the return value (set of states) to anything
			}
		}
		else {
			matrix[a->label][chr_number] = matrix[d->label][chr_number];																		// only left is missing, or both are
			if  (data->GetPolymorph(d->label,chr_number)) data->SetPolymorph(a->label,chr_number, true);
		}
	}
}

void ParsimonyClass::FinalPassUnord(NODE * const node, DATA * const data) {
	NODE *l=0, *r=0, *a=0, *f=0, *d=0;
	unsigned long long PA;
		
	l=node->left;
	r=node->right;

	if (node->anc) {	// this node is not the root
		PA=matrix[node->label][chr_number];
		a=node->anc;
		if (l->isAnc) { f=l; d=r; }		// if left desc is fixed
		else if (r->isAnc) { f=r; d=l; }	// if right desc is fixed
		if (!f)	{							// if this node is not fixed
			matrix[node->label][chr_number] = PA & matrix[a->label][chr_number];
			if (matrix[node->label][chr_number]!=matrix[a->label][chr_number]) {
			//  if both desc are nodes or tips without missing data
				if (matrix[l->label][chr_number] && matrix[r->label][chr_number]) {
					if (matrix[l->label][chr_number] & matrix[r->label][chr_number]) 
						matrix[node->label][chr_number]=(((matrix[l->label][chr_number] | matrix[r->label][chr_number]) & matrix[a->label][chr_number]) | PA);
					else matrix[node->label][chr_number]=(PA | matrix[a->label][chr_number]);
				}
			//  if both are missing data
				else if (!matrix[l->label][chr_number] && !matrix[r->label][chr_number])
					matrix[node->label][chr_number]=matrix[a->label][chr_number];
			// one branch is missing data, the other without missing data
				else matrix[node->label][chr_number]=(matrix[a->label][chr_number] | PA);
			}
		}
		else {			// this node is fixed
			if (!matrix[f->label][chr_number]) {
				matrix[node->label][chr_number] = matrix[a->label][chr_number] & matrix[node->label][chr_number];
				if (!matrix[node->label][chr_number]) matrix[node->label][chr_number]=(matrix[a->label][chr_number] | PA);
			}
		}			
	}

	if (!l->tip) FinalPassUnord(l, data);
	if (!r->tip) FinalPassUnord(r, data);
}

void ParsimonyClass::FinalPassOrd(NODE * const node, DATA * const data) {
	NODE *l=0, *r=0, *f=0, *d=0;
		
	l=node->left;
	r=node->right;
	if (node->anc)	{ // this node is not the root
		NODE *a=node->anc;
//		matrix[node->label][chr_number]=matrix[node->label][chr_number] & matrix[a->label][chr_number];
		if (matrix[node->label][chr_number]!=matrix[a->label][chr_number]) {	// if the final of the anc is the same as the current downpass set for this node, then it won't change and we don't have to check it
			if (l->isAnc) { f=l; d=r; }		// if left desc is fixed
			else if (r->isAnc) { f=r; d=l; }	// if right desc is fixed
			if (!f) {							// if this node is NOT fixed
				NODE *s=node->sister;
				matrix[node->label][chr_number] = OrdCombination(s, a, false, data);// ***not exactly the uppass set of this node because it includes the final (not uppass) of the anc
				matrix[node->label][chr_number] = OrdCombination(l, r, node, data);	// use the downpasses from the desc and the ***uppass of this node
			}
			else {			// this node is fixed
				if (!matrix[f->label][chr_number]) {
					matrix[node->label][chr_number] = matrix[a->label][chr_number] & matrix[d->label][chr_number];
					if (!matrix[node->label][chr_number]) matrix[node->label][chr_number]=OrdCombination(a, d, false, data);
				}
			}
		}
	}
	if (!l->tip) FinalPassOrd(l, data);
	if (!r->tip) FinalPassOrd(r, data);
}

void ParsimonyClass::DateTree(NODE *root, DATA *data) {
	if (!data->hasStratChr) {
		cout<<"\n\nNo Stratigraphic Character Present (DateNodes)\n";
		return;
	}
	else {
		chr_number=data->GetStratChr();
		if (!matrix) matrix=data->PassMatrix();
		uppass_required=false;	
		StratDownPassInt(root, data);
		if (uppass_required) StratUpPassInt(root, data);
	}
}

void ParsimonyClass::StratDownPassInt(NODE *a, DATA *data) {
	NODE *l=0, *r=0;

	l=a->left;
	r=a->right;
	if (!l->tip) StratDownPassInt(l, data);			// if left desc is a node, check it first
	if (!r->tip) StratDownPassInt(r, data);			// if right desc is a node, check it first
	if (*data->GetFInt(l)>=0 && *data->GetFInt(r)>=0) {
		if (*data->GetFInt(l)<=*data->GetFInt(r)) data->SetFInt(a, *data->GetFInt(l));	// if left is older than right set FInt anc to FInt left
		else data->SetFInt(a, *data->GetFInt(r));	// set FInt anc to FInt right
	}
	else {
		if (*data->GetFInt(l)<0) data->SetFInt(a, *data->GetFInt(r));
		else data->SetFInt(a, *data->GetFInt(l));
		uppass_required=true;											// uppass is required when there is missing stratigraphic data
	}
}

void ParsimonyClass::StratUpPassInt(NODE *node, DATA *data) {
	NODE *l, *r, *d=0, *a=0;
		
	l=node->left;
	r=node->right;
	if (*data->GetFInt(l)<0 && *data->GetFInt(r)<0) data->SetFInt(node, -1);
	else if (*data->GetFInt(l)<0) 
		d=node->right;
	else if (*data->GetFInt(r)<0) 
		d=node->left;
	
	if (d && node->anc) {
		a=node->anc;
		if (*data->GetFInt(d)>=*data->GetFInt(a) && *data->GetFInt(d)<=*data->GetFInt(a)) data->SetFInt(node, *data->GetFInt(d));
		else {
			data->SetFInt(node, *data->GetFInt(a));
			data->SetLInt(node, *data->GetFInt(d));
		}							
		}		
	if (!l->tip) StratUpPassInt(l, data);
	if (!r->tip) StratUpPassInt(r, data);
}

void ParsimonyClass::StratDownPassChr(NODE *a, DATA *data) {
	NODE *l=0, *r=0;

	l=a->left;
	r=a->right;
	if (!l->tip) StratDownPassChr(l, data);			// if left desc is a node, check it first
	if (!r->tip) StratDownPassChr(r, data);			// if right desc is a node, check it first
	if (matrix[l->label][chr_number] && matrix[r->label][chr_number]) {
		if (data->MinState(l->label, chr_number)<=data->MinState(r->label, chr_number)) matrix[a->label][chr_number]=pow(2, data->MinState(l->label, chr_number));	// if left is older than right set FInt anc to FInt left
		else matrix[a->label][chr_number] = pow (2, data->MinState(r->label, chr_number));	// set FInt anc to FInt right
	}
	else {
		if (!matrix[l->label][chr_number]) matrix[a->label][chr_number] = pow (2, data->MinState(r->label, chr_number));
		else matrix[a->label][chr_number] = pow (2, data->MinState(l->label, chr_number));
		uppass_required=true;
	}
}

void ParsimonyClass::StratUpPassChr(NODE *node, DATA *data)
{
	NODE *l, *r, *d=0, *a=0;
	int i;
		
	l=node->left;
	r=node->right;
	if (!matrix[l->label][chr_number] && !matrix[r->label][chr_number]) matrix[node->label][chr_number]=0;
	else if (!matrix[l->label][chr_number]) d=node->right;
	else if (!matrix[r->label][chr_number]) d=node->left;
	
	if (d) {
		a=node->anc;
		matrix[node->label][chr_number] = pow (2, data->MinSharedState(a->label, d->label, chr_number));
		if (!matrix[node->label][chr_number]) {
			if (data->MinState(a->label, chr_number) <= data->MinState(d->label, chr_number)) {
				i=0;
				while (i>=data->MinState(a->label, chr_number) && i<=data->MinState(d->label, chr_number)) {
					matrix[node->label][chr_number]+=pow(2, i);
					++i;
				}
			}
			else {
				i=0;
				while (i>=data->MinState(d->label, chr_number) && i<=data->MinState(a->label, chr_number)) {
					matrix[node->label][chr_number]+=pow(2, i);
					++i;
				}
			}
		}
	}		
	if (!l->tip) StratUpPassChr(l, data);
	if (!r->tip) StratUpPassChr(r, data);
}

float ParsimonyClass::GetStratDebt(NODE *root, DATA *data) {
	debt=0;
	DateTree(root, data);
	DoStratDebtInt(root, data);
	return debt;
}

void ParsimonyClass::DoStratDebtInt(NODE *a, DATA *data) {
	NODE *l=0, *r=0, *d=0, *f=0;
//	int gobots;
	
	l=a->left;
	r=a->right;
	if (!l->tip) DoStratDebtInt(l, data);
	if (!r->tip) DoStratDebtInt(r, data);

	if (l->isAnc) {		// if left branch is fixed as an ancestor
		f=a->left;		// keep track of label of fixed ancestor
		d=a->right;		// keep track of unfixed descendant
	}
	else if (r->isAnc) {// if right branch is fixed
		f=a->right;		// keep track of label of fixed ancestor
		d=a->left;		// keep track of unfixed descendant
	}

	if (!f || *data->GetFInt(f)<0) {	// if this node is not fixed or the fixed taxon is missing strat data
		if (*data->GetFInt(l)>=0) debt+=*data->GetChrWeight(data->GetStratChr())*((float)*data->GetFInt(l)-(float)*data->GetFInt(a));		// debt is FInt anc minus FInt desc
		if (*data->GetFInt(r)>=0) debt+=*data->GetChrWeight(data->GetStratChr())*((float)*data->GetFInt(r)-(float)*data->GetFInt(a));		// debt is FInt anc minus FInt desc
	}
	else {				// one taxon is fixed as anc
		if (*data->GetFInt(f)>*data->GetFInt(d)) 			// if FInt of desc is older than that of anc
			debt+=*data->GetChrWeight(data->GetStratChr())*((float)*data->GetFInt(f)-(float)*data->GetFInt(d));		// debt is FInt anc minus FInt desc

		else if	(*data->GetFInt(d)>*data->GetLInt(f))		// of if FInt desc is younger than LInt anc
			debt+=*data->GetChrWeight(data->GetStratChr())*(((float)*data->GetFInt(d)-(float)*data->GetLInt(f))-1);	// debt is FInt desc minus LInt anc minus one
	}
}
/*
void ParsimonyClass::DoStratDebtInt(NODE *a, DATA *data)
{
	int branch;
	bool fixed;
	NODE *descendant=0, *fixed_taxon=0;
//	int gobots;
	
	descendant=a->left;
	if (!descendant->tip) DoStratDebtInt(descendant, data);
	descendant=a->right;
	if (!descendant->tip) DoStratDebtInt(descendant, data);

	fixed=0;
	for (branch=0;branch<=1;++branch)				// check current node for fixed ancestors
		{
		if (!branch) descendant=a->left;		// check left branch first 
		else descendant=a->right;			// then right
		if (descendant->isAnc)						// if this taxon is fixed as an ancestor
			{
			fixed=1;								// label this node as fixed
			if (branch==0) 							// left branch							
				{
				fixed_taxon=a->left;			// keep track of label of fixed ancestor
				descendant=a->right;			// keep track of unfixed descendant
				}
			else									// right branch
				{
				fixed_taxon=a->right;		// keep track of label of fixed ancestor
				descendant=a->left;			// keep track of unfixed descendant
				}
			}
		if (fixed) break;
		}

	if (!fixed)	// if this node is not fixed
		{
		for (branch=0;branch<=1;++branch)
			{
			if (!branch) descendant=a->left;
			else descendant=a->right;
			debt+=*data->GetChrWeight(*data->GetStratChr())*(*data->GetFInt(descendant)-*data->GetFInt(a));		// debt is FInt anc minus FInt desc
			}
		}
	else		// one taxon is fixed as anc
		{
//		cout << "\nFA desc:"<<*data->GetFInt(descendant)<<"\tLA anc:"<<*data->GetLInt(fixed_taxon)<<"\n";
		if (*data->GetFInt(fixed_taxon)>*data->GetFInt(descendant)) 			// if FInt of desc is older than that of anc
			debt+=*data->GetChrWeight(*data->GetStratChr())*(*data->GetFInt(fixed_taxon)-*data->GetFInt(descendant));		// debt is FInt anc minus FInt desc

		else if	(*data->GetFInt(descendant)>*data->GetLInt(fixed_taxon))		// of if FInt desc is younger than LInt anc
			debt+=*data->GetChrWeight(*data->GetStratChr())*((*data->GetFInt(descendant)-*data->GetLInt(fixed_taxon))-1);	// debt is FInt desc minus LInt anc minus one
		}
}
*/
void ParsimonyClass::DoStratDebtChr(NODE *a, DATA *data) {
	bool is_fixed;
	NODE *l=0, *r=0, *d=0, *f=0;
//	int gobots;
	
	l=a->left;
	r=a->right;
	if (!l->tip) DoStratDebtChr(l, data);
	if (!r->tip) DoStratDebtChr(r, data);

	is_fixed=false;
	if (l->isAnc) {				// if this taxon is fixed as an ancestor
		is_fixed=true;			// label this node as fixed
		f=a->left;		// keep track of label of fixed ancestor
		d=a->right;		// keep track of unfixed descendant
	}
	else if (r->isAnc) {		// right branch
		is_fixed=true;			// label this node as fixed
		f=a->right;		// keep track of label of fixed ancestor
		d=a->left;		// keep track of unfixed descendant
	}

	if (!is_fixed){ // if this node is not fixed
		debt+=*data->GetChrWeight(chr_number)*((float)data->MinState(l->label, chr_number)-(float)data->MinState(a->label, chr_number));		// debt is FInt anc minus FInt desc
		debt+=*data->GetChrWeight(chr_number)*((float)data->MinState(r->label, chr_number)-(float)data->MinState(a->label, chr_number));		// debt is FInt anc minus FInt desc
	}
	else {		// one taxon is fixed as anc
		if (data->MinState(f->label, chr_number)>data->MinState(d->label, chr_number)) 			// if FInt of desc is older than that of anc
			debt+=*data->GetChrWeight(chr_number)*((float)data->MinState(f->label, chr_number)-(float)data->MinState(d->label, chr_number));		// debt is FInt anc minus FInt desc

		else if	(data->MinState(d->label, chr_number)>data->MaxState(f->label, chr_number))		// of if FInt desc is younger than LInt anc
			debt+=*data->GetChrWeight(chr_number)*(((float)data->MinState(d->label, chr_number)-(float)data->MaxState(f->label, chr_number))-1);	// debt is FInt desc minus LInt anc minus one
	}
}

float ParsimonyClass::SingleChrConsistencyIndex(NODE *tree, int *chr, DATA *data) {
	float CI, min_debt;
	NODE *root=0;
	
	chr_number=*chr;
	root=GetRoot(tree, data->GetNTaxa());
	min_debt=*data->GetNStates(chr_number)-1;
	debt=0;
	if (*data->GetChrType(*chr)==0) DownPassUnord(root, data);
	else if (*data->GetChrType(*chr)==1) DownPassOrd(root, data);
	else DownPassUnord(root, data);
	CI=min_debt/debt;
	
	return CI;
}

float ParsimonyClass::FullConsistencyIndex(NODE *tree, DATA *data) {
	float CI, min_debt;
	NODE *root=0;
	
	min_debt=debt=0.0;
	root=GetRoot(tree, data->GetNTaxa());
	for (int j=1;j<=*data->GetNChrs();++j) {
		if (*data->GetChrType(j)<4) {
			min_debt+=*data->GetNStates(j)-1;
			chr_number=j;
			if (*data->GetChrType(j)==0) DownPassUnord(root, data);
			else if (*data->GetChrType(j)==1) DownPassOrd(root, data);
		}
	}
	CI=min_debt/debt;
	return CI;
}

float ParsimonyClass::SingleChrRetentionIndex(NODE *tree, int *chr, DATA *data) {
	float RI, min_debt, max_debt;
	int i, k, max_state;
	NODE *root=0;
	int *state_count;
	
	chr_number=*chr;
	root=GetRoot(tree, data->GetNTaxa());
	state_count=new int[*data->GetNStates(chr_number)];
	for (k=0;k<*data->GetNStates(chr_number);++k) state_count[k]=0;
	for (i=1;i<=*data->GetNTaxa();++i) {
		for (k=0;k<*data->GetNStates(chr_number);++k)  if (data->HasState(i, chr_number, k)) ++state_count[k];
	}
	max_state=0;
	for (k=0;k<*data->GetNStates(chr_number);++k) if (state_count[k]>state_count[max_state]) max_state=k; 
	max_debt=*data->GetNTaxa()-state_count[max_state];	// this is only true for unordered chrs
	min_debt=*data->GetNStates(chr_number)-1;
	debt=0;
	DownPassUnord(root, data);
	delete [] state_count;
	if (min_debt==max_debt) RI=0;
	else RI=(max_debt-debt)/(max_debt-min_debt);
	
	return RI;
}

float ParsimonyClass::FullRetentionIndex(NODE *tree, DATA *data) {
	float RI, min_debt, max_debt;
	int i, j, k, max_state;
	NODE *root=0;
	int *state_count;
	
	root=GetRoot(tree, data->GetNTaxa());
	min_debt=max_debt=debt=0;
	for (j=1;j<=*data->GetNChrs();++j) {
		if (*data->GetChrType(j)==0) {
			state_count=new int[*data->GetNStates(j)];
			for (k=0;k<*data->GetNStates(j);++k) state_count[k]=0;
			for (i=1;i<=*data->GetNTaxa();++i) {
				for (k=0;k<*data->GetNStates(j);++k)  if (data->HasState(i, j, k)) ++state_count[k];
//				for (k=0;k<*data->GetNStates(j);++k)  if (*data->GetMatrix(i, j, k)==2) ++state_count[k];
			}
		//	for (k=0;k<*data->GetNStates(j);++k) printf ("\tState %d observed in %d taxa\n", k, state_count[k]);

			max_state=0;
			for (k=0;k<*data->GetNStates(j);++k) if (state_count[k]>state_count[max_state]) max_state=k; 
			max_debt+=*data->GetNTaxa()-state_count[max_state];		// this is only true for unordered chrs
			min_debt+=*data->GetNStates(j)-1;						// this is only true for unordered chrs
			chr_number=j;
//			debt+=SingleChrParsDebt(tree, root, data);
			DownPassUnord(root, data);			
			delete [] state_count;
		}
	}
	if (min_debt==max_debt) RI=0;
	else RI=(max_debt-debt)/(max_debt-min_debt);
	
	return RI;
}

float ParsimonyClass::SingleChrStratRetentionIndex(NODE *tree, int chr, DATA *data) {
	float SRI, min_debt, max_debt;
	int i, k, min_int;
	bool flag;
	NODE *root=0;
	
	chr_number=chr;
	
	root=GetRoot(tree, data->GetNTaxa());
	max_debt=0;
	for (i=1;i<=*data->GetNTaxa();++i) max_debt+=*data->GetFInt(i);
//	min_int=*data->GetNStates(data->GetStratChr());
	flag=false;
	k=0;
	while (!flag) {
		for (i=1;i<=*data->GetNTaxa();++i) {
			if (data->HasState(i, data->GetStratChr(), k)) {
				min_int=k;
				flag=true;
			}
			if (flag) break;
			else ++k;
		}
	}
	min_debt=0;
	for (k=min_int;k<*data->GetNStates(data->GetStratChr());++k) {
		flag=false;
		for (i=1;i<=*data->GetNTaxa();++i) {
			if (*data->GetFInt(i)<=k && *data->GetLInt(i)>=k) flag=true;
			if (flag) break;
		}
		if (!flag) ++min_debt;
	}
	debt=0;
	DoStratDebtInt(root, data);
	
	if (min_debt==max_debt) SRI=0;
	else SRI=(max_debt-debt)/(max_debt-min_debt);
	
	return SRI;
}

float ParsimonyClass::SingleChrRescaledCI(NODE *tree, int *chr, DATA *data) {
	float RCI;
	
	RCI=SingleChrConsistencyIndex(tree, chr, data)*SingleChrRetentionIndex(tree, chr, data);
	return RCI;
}

float ParsimonyClass::FullRescaledCI(NODE *tree, DATA *data) {
	float RCI;
	
	RCI=FullConsistencyIndex(tree, data)*FullRetentionIndex(tree, data);
	return RCI;
}

unsigned long long ParsimonyClass::OrdCombination(NODE *l, NODE *r, bool doDebt, DATA *data) {
	unsigned long long combination;
	
	combination = matrix[l->label][chr_number] & matrix[r->label][chr_number];
	if (!combination) {
		int max_state, min_state;
		if (matrix[l->label][chr_number]>matrix[r->label][chr_number]) {
			max_state=data->MaxState(r->label, chr_number);
			min_state=data->MinState(l->label, chr_number);
		}
		else {
			max_state=data->MaxState(l->label, chr_number);
			min_state=data->MinState(r->label, chr_number);
		}
		for (int i=max_state;i<=min_state;++i) combination += pow(2, i);
		if (doDebt) debt += (((float)min_state-(float)max_state) * *data->GetChrWeight(chr_number));
	}
	return combination;
}

unsigned long long ParsimonyClass::OrdCombinationPoly(NODE *l, NODE *r, bool doDebt, DATA *data) {
	unsigned long long combination;
	int i;
//	int gobots;
	
	combination = matrix[l->label][chr_number] & matrix[r->label][chr_number];
	if (!combination) {
		int max_state, min_state;
		if (matrix[l->label][chr_number]>matrix[r->label][chr_number]) {
			max_state=data->MaxState(r->label, chr_number);
			min_state=data->MinState(l->label, chr_number);
		}
		else {
			max_state=data->MaxState(l->label, chr_number);
			min_state=data->MinState(r->label, chr_number);
		}
		for (i=max_state;i<=min_state;++i) combination += pow(2, i);
		if (doDebt) debt += (((float)min_state-(float)max_state) * *data->GetChrWeight(chr_number));
	}

	if (doDebt) {
		if (data->GetPolymorph(r->label, chr_number))	// if r is polymorphic, add in debt for each state lower than the highest state
			for (i=0;i<*data->GetNStates(chr_number);++i) if ((matrix[r->label][chr_number] & combination) & (unsigned long long)pow(2, i)) debt+=*data->GetChrWeight(chr_number);
		if (data->GetPolymorph(l->label, chr_number))	// if l is polymorphic, add in debt for each state higher than the lowest state
			for (i=0;i<*data->GetNStates(chr_number);++i) if ((matrix[l->label][chr_number] & combination) & (unsigned long long)pow(2, i)) debt+=*data->GetChrWeight(chr_number);
	}
	return combination;
}

unsigned long long ParsimonyClass::OrdCombination(unsigned long long l, unsigned long long r, DATA *data) {
	unsigned long long combination;
	
	combination = l & r;
	if (!combination) {
		int max_state=0, min_state=0, i;
		if (l > r) {
			for (i=*data->GetMaxObsState();i>=0;--i) {
				if (r & (unsigned long long) (pow(2, i))) max_state=i;
				if (max_state) break;
			}
			for (i=0;i<=*data->GetMaxObsState();++i) {
				if (l & (unsigned long long) (pow(2, i))) min_state=i;
				if (min_state) break;
			}
		}
		else {
			for (i=*data->GetMaxObsState();i>=0;--i) {
				if (l & (unsigned long long) (pow(2, i))) max_state=i;
				if (max_state) break;
			}
			for (i=0;i<=*data->GetMaxObsState();++i) {
				if (r & (unsigned long long) (pow(2, i))) min_state=i;
				if (min_state) break;
			}
		}
		for (i=max_state;i<=min_state;++i) combination += pow(2, i);
	}
	return combination;
}

unsigned long long ParsimonyClass::OrdCombination(NODE *l, NODE *r, NODE *a, DATA *data) {
	unsigned long long A, B, C, answer;
//	int gobots;
	
	A=matrix[l->label][chr_number];
	B=matrix[r->label][chr_number];
	C=matrix[a->label][chr_number];
	
	answer=(A & B & C);
	if (answer) return answer;
	else return  ((OrdCombination(OrdCombination(A, B, data), C, data)) & 
				  (OrdCombination(OrdCombination(A, C, data), B, data)) &
				  (OrdCombination(OrdCombination(B, C, data), A, data)));

}

