/* Functions to evolve phylogenies and character states on phylogenies
/*		Written by The Wagner
/*	v 0.0	12/1995
/*  v 0.1   05/1996
/*  v 0.2	01/1997
/*  v 0.3	07/1997
/*  v 0.4	03/1998
/*	v 0.5	05/1999
/*	v 0.6	12/1999
/*	v 0.7	02/2000
/*	v 0.8	05/2002
/*	v 0.9	02/2003
/*	v 1.0	06/2003
/*	v 1.1	09/2003 - added evolveorderedinclade & evolveunorderedinclade
/*	v 1.2	05/2005 - added ensurecharstates & evolveunorderedinclade
/*****************************************************************************************************************/ 
#ifdef MontCarloPhylogenyFunctions

	#include <stdlib.h>
	#include <stdio.h>
	#include <time.h>
	#include <math.h>
	#include <string.h>
	#define matrixchange
	#include "matrixchange.h"
	#define matrixreading
	#include "matrixreading.h"
	#define memory
	#include "memory.h"
	#define minmax
	#include "minmax.h"
	#define Probability
	#include "probability.h"
	#define TreeRead
	#include "tree_read.h"

	long **evolvetree(int notu, double *mbl, int fossils);
	long **evolvetreeVenn(int notu, double *mbl, int fossils);
	long **evolvecladogram(int notu, long **tree);
	long **evolvematrix(long **tree, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int ttlstp, int UNKNOWN, int INAP);
	long *evolvecompat(long **tree, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int minstp, int maxstp, int comptype, int UNKNOWN, int INAP);
	long **evolvetocompat(long **tree, int tcomp, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int *depend, int comptype, int UNKNOWN, int INAP);
	long **evolvetocompatdelt(long **tree, int tcomp, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int *depend, int comptype, int UNKNOWN, int INAP, int /*&*/deltas);
	long **evolvetocompatsavessteps(long **tree, int tcomp, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int *depend, int *chsteps, int comptype, int UNKNOWN, int INAP);
	long **EvolveCharacterRate(int charn, double pi, long **tree, int notu, long **matrix, int *nstates, int *ctype, int *bias, int UNKNOWN, int INAP);
	long *evolvecharacterNsteps(long **tree, int N, int notu, long *chvector, int nstates, int ctype, int bias, int UNKNOWN, int INAP);
	long *evolvecharacterNstepsovertime(long **tree, int N, int notu, long *chvector, int nstates, int ctype, int bias, int UNKNOWN, int INAP);
	void evolveorderedinclade(long *clade, long **matrix, int ch, int delta, int INAP, int UNKNOWN);
	void multistatevett(long **matrix, int ch, int st, int notu, int INAP, int UNKNOWN);
	void evolveunorderedinclade(long *clade, long **matrix, int ch, int state, int INAP, int UNKNOWN);
	void evolvebinaryinclade(long *clade, long **matrix, int ch);
	long swap(long n);
	long mswap(long n, int X);
	long multi(long c, int bias);
	double MinRate(int notu);
	double MaxRate(int notu, int MXST);
	int branchnumber(long **tree, int notu);
	long **evolveadditivedependent(int charn, int steps, long **tree, int notu, long **matrix, int *nstates, int *ctype, int *bias, int UNKNOWN, int INAP, int APS);
	int	**descendantnodes (long ** tree, int notu);
	int number_of_nodes(long ** tree, int notu);
	double *getfossilparams();
	void derive_each_char(long **matrix, long **invmatrix, int **taxachange, int *nstates, int *ctype, int *maxch, int *steps, int *bias, int *sychos, int *aptaxa, int *mpd, int *trpd, int notu, int nchars, int maxstp, int mxdel, int UNKNOWN, int INAP, long **tree, int *branches, int *available, int nodes, int ttlbr, int *deltas);
	void ensure_charstates(long **matrix, long **invmatrix, int **taxachange, int *nstates, int *ctype, int *maxch, int *steps, int *bias, int *mpd, int *trpd, int notu, int nchars, int UNKNOWN, int INAP, long **tree, int nodes/*, int *deltas*/);
	void ensure_apomorphies(long **matrix, long **invmatrix, int **taxachange, int *nstates, int *ctype, int *maxch, int *steps, int *bias, int *mpd, int *trpd, int notu, int nchars, int UNKNOWN, int INAP, long **tree, int nodes/*, int *deltas*/);
	void evolvedependent(long **matrix,long **tree,int *trpd,int ttlbr,int notu,int *steps,int ind,int dep,int nstates,int ctype,int bias,int **taxachange,int *deltas,int UNKNOWN,int INAP);
	void evolvedependentNsteps(int N, int iN, long **matrix,long **tree, int *trpd, int **taxachange, int *steps, int notu, int ttlbr, int ind,int dep,int nstates,int ctype,int bias,int UNKNOWN,int INAP, int *deltas);
	void recursivechange(long **tree, long **matrix, int node, int chr, int newst, int notu);
	
#else

	extern long **evolvetree(int notu, double *mbl, int fossils);
	extern long **evolvetreeVenn(int notu, double *mbl, int fossils);
	extern long **evolvecladogram(int notu, long **tree);
	extern long **evolvematrix(long **tree, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int ttlstp, int UNKNOWN, int INAP);
	extern long *evolvecompat(long **tree, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int minstp, int maxstp, int comptype, int UNKNOWN, int INAP);
	extern long **evolvetocompat(long **tree, int tcomp, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int *depend, int comptype, int UNKNOWN, int INAP);
	extern long **evolvetocompatsavessteps(long **tree, int tcomp, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int *depend, int *chsteps, int comptype, int UNKNOWN, int INAP);
	extern long **evolvetocompatdelt(long **tree, int tcomp, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int *depend, int comptype, int UNKNOWN, int INAP, int /*&*/deltas);
	extern long **EvolveCharacterRate(int charn, double pi, long **tree, int notu, long **matrix, int *nstates, int *ctype, int *bias, int UNKNOWN, int INAP);
	extern long *evolvecharacterNsteps(long **tree, int N, int notu, long *chvector, int nstates, int ctype, int bias, int UNKNOWN, int INAP);
	extern long *evolvecharacterNstepsovertime(long **tree, int N, int notu, long *chvector, int nstates, int ctype, int bias, int UNKNOWN, int INAP);
	extern void evolveorderedinclade(long *clade, long **matrix, int ch, int delta, int INAP, int UNKNOWN);
	extern void multistatevett (long **matrix, int ch, int st, int notu, int INAP, int UNKNOWN);
	extern void evolveunorderedinclade(long *clade, long **matrix, int ch, int state, int INAP, int UNKNOWN);
	extern void evolvebinaryinclade(long *clade, long **matrix, int ch);
	extern long swap(long n);
	extern long mswap(long n, int X);
	extern long multi(long c, int bias);
	extern double MinRate(int notu);
	extern double MaxRate(int notu, int MXST);
	extern int branchnumber(long **tree, int notu);
	extern long **evolveadditivedependent(int charn, int steps, long **tree, int notu, long **matrix, int *nstates, int *ctype, int *bias, int UNKNOWN, int INAP, int APS);
	extern int	**descendantnodes (long ** tree, int notu);
	extern int number_of_nodes(long ** tree, int notu);
	extern double *getfossilparams();
	extern void derive_each_char(long **matrix, long **invmatrix, int **taxachange, int *nstates, int *ctype, int *maxch, int *steps, int *bias, int *sychos, int *aptaxa, int *mpd, int *trpd, int notu, int nchars, int maxstp, int mxdel, int UNKNOWN, int INAP, long **tree, int *branches, int *available, int nodes, int ttlbr, int *deltas);
	extern void ensure_charstates(long **matrix, long **invmatrix, int **taxachange, int *nstates, int *ctype, int *maxch, int *steps, int *bias, int *mpd, int *trpd, int notu, int nchars, int UNKNOWN, int INAP, long **tree, int nodes/*, int *deltas*/);
	extern void ensure_apomorphies(long **matrix, long **invmatrix, int **taxachange, int *nstates, int *ctype, int *maxch, int *steps, int *bias, int *mpd, int *trpd, int notu, int nchars, int UNKNOWN, int INAP, long **tree, int nodes/*, int *deltas*/);
	extern void evolvedependent(long **matrix,long **tree,int *trpd,int ttlbr,int notu,int *steps,int ind,int dep,int nstates,int ctype,int bias,int **taxachange,int *deltas,int UNKNOWN,int INAP);
	extern void evolvedependentNsteps(int N, int iN, long **matrix,long **tree, int *trpd, int **taxachange, int *steps, int notu, int ttlbr, int ind,int dep,int nstates,int ctype,int bias,int UNKNOWN,int INAP, int *deltas);
	extern void recursivechange(long **tree, long **matrix, int node, int chr, int newst, int notu);

#endif

