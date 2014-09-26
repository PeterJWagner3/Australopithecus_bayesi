#ifndef MAXSTATES
#define MAXSTATES 34
#endif

#ifndef MAXLABELSIZE
#define MAXLABELSIZE 30
#endif

#ifndef DATA_FIELD_H
#define DATA_FIELD_H

#include "node.h"
#include "math.h"
#include "memory.h"
#include "memory_koz.h"
#include <iostream>
#include <sstream>
#include <stdio.h>

class DATA
	{
	public:
		DATA();
		DATA(int , int );
		~DATA();

		void MakeMatrixFromFile(char *);
		void MakeMatrixFromFile(char *, long);
		void MakeStrat();
		void MakeStratSVP();
		void MakePattern();
		void MakeNStates();
		void DefaultTaxonLabels();
		void AddChrsToMatrix(int );
		void BootstrapMatrixByTaxa();
		void BootstrapMatrixByChrs();

		void SetDNA(bool );
		void SetNTaxa(int );
		void SetNChrs(int );
		void SetStratChr(int );
		void SetMatrix(int , int, unsigned long long );
		void SetNStates(int , int);
		void SetChrType(int , int);
		void SetFInt(int , int);
		void SetLInt(int, int);
		void SetFInt(NODE *, int);
		void SetLInt(NODE *, int);
		void IncludeChr(int);
		void ExcludeChr(int);
		void SetInfChr(int, bool);
		void SetChrWeight(int, float);
		void SetPolymorph(int , int , bool );

		bool IsDNA() const;
		int *GetNTaxa() const;
		int *GetNChrs() const;
		int GetStratChr() const;
		unsigned long long *GetMatrix(int , int ) const;
		int *GetNStates(int ) const;
		int *GetMaxObsState();
		int *GetChrType(int ) const;
		int *GetFInt(int) const;
		int *GetLInt(int) const;
		int *GetFInt(NODE *) const;
		int *GetLInt(NODE *) const;
		int *GetPattern(int) const;
		bool GetInc(int ) const;
		bool GetInfChr(int ) const;
		bool GetPolymorph(int, int ) const;
		char *GetTaxonName() const; 
		float *GetChrWeight(int ) const; 

		bool HasState(int , int , int );
		bool ShareAnyState(int , int , int );
		bool IsMissing(int , int ); 
		bool isUncertain(int taxon1, int chr);
		int getNumStates(int , int );
		int MinState(int , int );
		int MaxState(int , int );
		int GetMaxInt() const;
		int MinSharedState(int , int , int );
		int **outputStratAsMatrix();

		unsigned long long **PassMatrix() const;
		unsigned long long **PassBootMatrix() const;
		
		void PrintMatrix();
		void PrintFullMatrix();
		void PrintSingleChr (int );
		void PrintStrat() const;
		void PrintInclusion() const;
		void PrintInfChr() const;
		void PrintPolymorph () const;
		void PrintNStatesPerTaxon ();
		void PrintNStates() const;		
		bool hasPolymorph;
		bool hasStratChr;
		bool FAChrPresent;
		bool LAChrPresent;
		char **taxon_label;// ntaxa x MAXSTATELABEL array of taxon names
	
	private:
		bool DNA;		// 1 if data is DNA sequence data, 0 if categorical
		int *ntaxa;		// number of taxa in matrix
		int *nchrs;		// number of characters in matrix
		int stratchr;	// label (char number) of stratigraphic character
		int fachr;		// label (char number) of stratigraphic character
		int lachr;		// label (char number) of stratigraphic character
		int MaxInt;		// label of the maximum stratigraphic interval
		unsigned long long **matrix;	// ntaxa x nchrs matrix of binary representations of data matrix
		unsigned long long **boot_matrix;	// bootstrapped ntaxa x nchrs matrix of binary representations of data matrix
		int *nstates;	// nchrs-long array of observed numbers of states per character
		int *chrtype;	// nchrs-long array of character types
		bool *included;	// nchrs-long array of character inclusion status
		bool *infchr;	// nchrs-long array of informative chr status
		float *weight;	// nchrs-long array of character weights
		bool **polymorph;// ntaxa x nchrs matrix of polymophism
		int maxobsstate;// the maximum observed number of states in the matrix
		int *pattern;	// number of times this state pattern is observed		
		int *fint;		// 2*ntaxa-long array of first appearance intervals of taxa
		int *lint;		// 2*ntaxa-long array of last appearance intervals of taxa
		
		void AllocateChrType();
		void AllocateMatrix();
		void AllocateTaxonLabels();
		void AllocateChrWeights();
		void AllocateIncludedChrs();

		void MakeInfChrs();
		void MakePolymorph();

		void GetTaxonNameFromFile(FILE *, int * const );
		bool GetStateFromFile(FILE *, int *, int *);
		bool ChrToState(char *, int *);
		void StateToChr(int *, char *);
		
 	};
	
#endif

