#ifndef PARISMONY_CLASS_FIELD_H
#define PARISMONY_CLASS_FIELD_H

	#include <iostream>
	#include "data_field.h"
	#include "tree_fns.h"

class ParsimonyClass
	{
	public:
		ParsimonyClass();
		ParsimonyClass(DATA * const );
		~ParsimonyClass();
		float TotalDebt (NODE *, DATA *);
		float TotalDebtPolymorph(NODE *, DATA *);
//		float TotalDebt (NODE *, NODE *, DATA *, unsigned long long *);
		void DateTree(NODE *, DATA *);
		float GetStratDebt (NODE *, DATA *);
		float SingleChrDebt (NODE * const, DATA * const, int * );
		
		void GetDownPassSet(NODE * , DATA * , int *);
//		void GetUpPassSet(NODE * , DATA * , bool **, int *);
		void GetFinalSet(NODE * const , DATA * const , int *);
		void MakeFullFS(NODE * const , DATA * const );

	  // Chr Stats
		float SingleChrConsistencyIndex(NODE *tree, int *chr, DATA *data);
		float FullConsistencyIndex(NODE *tree, DATA *data);
		float SingleChrRetentionIndex(NODE *tree, int *chr, DATA *data);
		float FullRetentionIndex(NODE *tree, DATA *data);
		float SingleChrRescaledCI(NODE *tree, int *chr_number, DATA *data);
		float FullRescaledCI(NODE *tree, DATA *data);
		float SingleChrStratRetentionIndex(NODE *tree, int chr_number, DATA *data);

	private:
		float debt;
		int chr_number;
		bool uppass_required;
		unsigned long long ** matrix;
//		DATA * const data;
				
		void DownPassUnord(NODE * const , DATA * const );
		void DownPassOrd(NODE * const , DATA * const );
		void DownPassUnordPoly(NODE * const , DATA * const );
		void DownPassOrdPoly(NODE * const , DATA * const );

//		void UpPassUnord (NODE * const , DATA * const );

		void FinalPassUnord (NODE * const , DATA * const );
		void FinalPassOrd (NODE * const , DATA * const );

		void DoStratDebtInt(NODE *, DATA *);
		void StratDownPassInt(NODE *, DATA *);
		void StratUpPassInt(NODE *, DATA *);
		void DoStratDebtChr(NODE *, DATA *);
		void StratDownPassChr(NODE *, DATA *);
		void StratUpPassChr(NODE *, DATA *);

//		float SingleChrParsDebt (NODE *, NODE *, DATA *, float **);
//		float SingleChrParsDebt (NODE *, NODE *, DATA *);
//		void GetLengthSet (NODE * const, float **, float **, DATA *);

		float ***MasterMatrix(int *);
		float **unordered_stepmatrix(int );
		float **ordered_stepmatrix(int );
		float **irrevup_stepmatrix(int );
		float **irrevdn_stepmatrix(int );
		
		unsigned long long OrdCombination(NODE *l, NODE *r, bool , DATA * );
		unsigned long long OrdCombinationPoly(NODE *l, NODE *r, bool , DATA * );
		unsigned long long OrdCombination(unsigned long long l, unsigned long long r, DATA * );
		unsigned long long OrdCombination(NODE *l, NODE *r, NODE *a, DATA * );


	} ;

#endif

