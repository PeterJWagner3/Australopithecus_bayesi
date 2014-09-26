#ifndef CHR_LIKELIHOOD_H
#define CHR_LIKELIHOOD_H

//	#include <iostream>
	#include "math.h"
	#include "data_field.h"
	#include "tree_fns.h"
	#include "probability.h"
	#include "memory_koz.h"

class LikelihoodClass
	{
	public:
		LikelihoodClass();
		LikelihoodClass(DATA * const );
		~LikelihoodClass();
		
		double getFullML(NODE * );
		double PoissonProbablity(bool change, NODE *l);

	private:
		double lnL;						// overall likelihood of the tree
		int chr_number;
		unsigned long long ** matrix;
		double ***cL;				// node x chr x state cube of conditional likelihoods
		double *jL;					// node-long array of joint likelihoods for each node
		double *p∆;					// nchr-long array of rates of change
		DATA *data;
				
		void getFullMLMachine(NODE * );
		void makeJointLikelihoods();

	} ;

#endif

