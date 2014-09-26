#include "chr_likelihood.h"

using namespace std;

LikelihoodClass::LikelihoodClass() { data=0; }
LikelihoodClass::LikelihoodClass(DATA * const in_data) { data=in_data; }
LikelihoodClass::~LikelihoodClass() {};

double LikelihoodClass::getFullML(NODE *this_node) {
	int k;	
	pÆ=new double [*data->GetNChrs()+1];    	/* vector of change probabilities per character */
//	jL=new double [2**data->GetNTaxa()];
	for (int i=1;i<=*data->GetNChrs();i++) pÆ[i]=0.05;		/* sets probs at 0.05: change this to read file	*/
	/* conditional likelihood for each possible state of each character at each node */
	cL =dcube(2**data->GetNTaxa(), *data->GetNChrs()+1, *data->GetMaxObsState());

	/* conditional likelihoods of terminal taxa, states within taxa within characters loop	*/
	for (chr_number=1;chr_number<=*data->GetNChrs();chr_number++) {
		if (*data->GetChrType(chr_number)!=4) {				/* 4 is a stratigraphic character	*/
			for (int i=1;i<=*data->GetNTaxa();i++) {		/* taxon loop						*/
				for (int k=0;k<*data->GetNStates(chr_number);k++) {		/* state loop within character	*/
					if (data->HasState(i, chr_number, k)) cL[i][chr_number][k]=1;
					else cL[i][chr_number][k]=0;				
				}	/* end state loop	*/
			}	/* end taxon loop	*/
			getFullMLMachine(this_node);					/* recursive alogirthm to get conditional likelihood of node + descendant nodes */
		}	/* make sure that this is a morphological character!	*/
	}	/* end character loop */

	double holdL;			/* holder for joint likelihood values	*/
	lnL=0;					/* log-likelihood of analyzed portion (usually whole tree	*/
	/* joint likelihood = sum of the conditional likelihoods for all possible states	*/
	for (chr_number=1;chr_number<=*data->GetNChrs();chr_number++) {
		if (*data->GetChrType(chr_number)!=4) {
			holdL=0.0;								/* hold log-likelihood sums	*/
			// sum likelihoods over all possible ancestral state reconstructions */
			for (k=0;k<*data->GetNStates(chr_number);k++) holdL+=cL[this_node->label][chr_number][k];
			lnL+=log(holdL); 						// multiply (add log) likelihoods over all characters
		}	/* analyze only morphological characters	*/
	}	/* end character loop	*/

	delete [] pÆ;												/* frees memory of a vector	*/
	free_cube(cL, 2**data->GetNTaxa(), *data->GetNChrs()+1);
	
	return -lnL;
}

void LikelihoodClass::getFullMLMachine(NODE *this_node) {
	NODE *l, *r;
	double leftL, rightL;
	int k, ks;
	
	l=this_node->left;	
	r=this_node->right;	
	if (!l->tip) getFullMLMachine(l);
	if (!r->tip) getFullMLMachine(r);

	for (k=0;k<*data->GetNStates(chr_number);k++) {	// ancestral state of this node
		leftL=rightL=0.0;
		for (ks=0;ks<*data->GetNStates(chr_number);++ks) {	// sum likelihoods over all states of descendant node
			if (ks==k) leftL+=cL[l->label][chr_number][ks] * PoissonProbablity(false, l);	/* MODIFY THIS LATER FOR HIDDEN REVERSALS: Poisson(pÆ[chr_number],l->bl,0)+((paths_back(1,*data->GetNStates(chr_number))/(pow(GetNStates(chr_number)-1,1)))*Poisson(pÆ[chr_number],l->bl,1)+((paths_back(2,*data->GetNStates(chr_number))/(pow(GetNStates(chr_number)-1,2)))*Poisson(pÆ[chr_number],l->bl,2)+((paths_back(3,*data->GetNStates(chr_number))/(pow(GetNStates(chr_number)-1,3)))*Poisson(pÆ[chr_number],l->bl,3)+((paths_back(4,*data->GetNStates(chr_number))/(pow(GetNStates(chr_number)-1,4)))*Poisson(pÆ[chr_number],l->bl,4)+((paths_back(5,*data->GetNStates(chr_number))/(pow(GetNStates(chr_number)-1,5)))*Poisson(pÆ[chr_number],l->bl,5);*/
			else leftL+=cL[l->label][chr_number][ks] * PoissonProbablity(true, l);			/* MODIFY THIS LATER FOR MULTIPLE CHANGES: =Poisson(pÆ[chr_number],l->bl,0)+((pow(GetNStates(chr_number)-1,1))+(paths_back(1,*data->GetNStates(chr_number)))/(pow(GetNStates(chr_number)-1,1)))*Poisson(pÆ[chr_number],l->bl,1)+((pow(GetNStates(chr_number)-1,2))+(paths_back(2,*data->GetNStates(chr_number)))/(pow(GetNStates(chr_number)-1,2)))*Poisson(pÆ[chr_number],l->bl,2)+((pow(GetNStates(chr_number)-1,3))+(paths_back(3,*data->GetNStates(chr_number)))/(pow(GetNStates(chr_number)-1,3)))*Poisson(pÆ[chr_number],l->bl,3)+((pow(GetNStates(chr_number)-1,4))+(paths_back(4,*data->GetNStates(chr_number)))/(pow(GetNStates(chr_number)-1,4)))*Poisson(pÆ[chr_number],l->bl,4)+((pow(GetNStates(chr_number)-1,5))+(paths_back(5,*data->GetNStates(chr_number)))/(pow(GetNStates(chr_number)-1,5)))*Poisson(pÆ[chr_number],l->bl,5)
	*/
		}
		for (ks=0;ks<*data->GetNStates(chr_number);++ks) {	// sum likelihoods over all states of descendant node
			if (ks==k) rightL+=cL[r->label][chr_number][ks] * PoissonProbablity(false, r);	/* MODIFY THIS LATER FOR HIDDEN REVERSALS	*/
			else rightL+=cL[r->label][chr_number][ks] * PoissonProbablity(true, r);		/* MODIFY THIS LATER FOR MULTIPLE CHANGES	*/
		}
	cL[this_node->label][chr_number][k]+=leftL*rightL;	// conditional likelihood for state K - multiply likelihoods from each branch
	}
}

double LikelihoodClass::PoissonProbablity(bool change, NODE *l)	{
	if (!change)	{
		return Poisson(pÆ[chr_number],l->bl,0)
			+((paths_back(1,*data->GetNStates(chr_number))/(pow(*data->GetNStates(chr_number)-1,1)))*Poisson(pÆ[chr_number],l->bl,1))
			+((paths_back(2,*data->GetNStates(chr_number))/(pow(*data->GetNStates(chr_number)-1,2)))*Poisson(pÆ[chr_number],l->bl,2))
			+((paths_back(3,*data->GetNStates(chr_number))/(pow(*data->GetNStates(chr_number)-1,3)))*Poisson(pÆ[chr_number],l->bl,3))
			+((paths_back(4,*data->GetNStates(chr_number))/(pow(*data->GetNStates(chr_number)-1,4)))*Poisson(pÆ[chr_number],l->bl,4))
			+((paths_back(5,*data->GetNStates(chr_number))/(pow(*data->GetNStates(chr_number)-1,5)))*Poisson(pÆ[chr_number],l->bl,5));
		}
	else	{
		return Poisson(pÆ[chr_number],l->bl,0)
			+((pow(*data->GetNStates(chr_number)-1,1))+(paths_back(1,*data->GetNStates(chr_number)))/(pow(*data->GetNStates(chr_number)-1,1)))*Poisson(pÆ[chr_number],l->bl,1)
			+((pow(*data->GetNStates(chr_number)-1,2))+(paths_back(2,*data->GetNStates(chr_number)))/(pow(*data->GetNStates(chr_number)-1,2)))*Poisson(pÆ[chr_number],l->bl,2)
			+((pow(*data->GetNStates(chr_number)-1,3))+(paths_back(3,*data->GetNStates(chr_number)))/(pow(*data->GetNStates(chr_number)-1,3)))*Poisson(pÆ[chr_number],l->bl,3)
			+((pow(*data->GetNStates(chr_number)-1,4))+(paths_back(4,*data->GetNStates(chr_number)))/(pow(*data->GetNStates(chr_number)-1,4)))*Poisson(pÆ[chr_number],l->bl,4)
			+((pow(*data->GetNStates(chr_number)-1,5))+(paths_back(5,*data->GetNStates(chr_number)))/(pow(*data->GetNStates(chr_number)-1,5)))*Poisson(pÆ[chr_number],l->bl,5);
		}
	}

void LikelihoodClass::makeJointLikelihoods() {
	int i,k;
	double holdL;
	for (i=*data->GetNTaxa()+1;i<2**data->GetNTaxa();i++) {
		jL[i]=0.0;
		for (chr_number=1;chr_number<=*data->GetNChrs();chr_number++) {
			if (*data->GetChrType(chr_number)!=4) {
				holdL=0.0;
				for (k=0;k<*data->GetNStates(chr_number);k++) holdL+=cL[i][chr_number][k];	// sum likelihoods over all possible ASRs
				jL[i]+=log(holdL); 						// multiply likelihoods over all characters
			}
		}
	}
}
