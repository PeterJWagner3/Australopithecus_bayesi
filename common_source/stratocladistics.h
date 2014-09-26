/*************************************************************v0.0 - 2001.05.03 BY WRITTEN BY P.J.WAGNER III v0.1 - 2003.01.28 - SCI, GER, RCI, MSM, etc. added.  v0.2 - 2005.02.19 - stratigraphic compatibility added.  *************************************************************/#ifdef stratocladistics	#define BUD			0	#define BIF			1	#define TREE_END	-1	#include <stdlib.h>	#include <stdio.h>	#include <math.h>//	#include "treeread.h"	#include "matrixreading.h"	#include "memory.h"	#include "sort.h"	long *dateclade(long *fa, long **tree, int clades, int notu);	long **datecladefull(long **range, long **tree, int clades, int notu);	long **datecladesim(long **tree, int clades, int notu);	double *datecladereal(double *fa, long **tree, int clades, int notu, int sign);	void datecladerealape(long *ape, double **ranges, double **divergences, int nodes, int notu, unsigned long *ancestral);	void datecladerealape_pulley(int HTU, long *ape, double **divergences, int nodes, int notu, unsigned long *ancestral);	void datecladerealape_one_node(int HTU, int desc[], double **divergences, int f1, int notu, unsigned long *ancestral);	double sum_divergence_times(double **divergences, int notu, int nodes);	void datecladerealape_addexp(long *ape, double **ranges, double **divergences, int nodes, int notu, double lambda, unsigned long *ancestral);	void dateclade_stratlikelihood_ape_addexp(long *ape, double **divergences, int nodes, int notu, double lambda, unsigned long *ancestral, unsigned long *minbr, double **ranges, double ***thetas, double **theta_mods, int freqs, double *stratlikes, double **timescale, long **initsampling, long **stagesampling, double **geoglikes, int areas, double iota);	void dateclade_stratolikelihood_one_node(int HTU, double ***thetas, double **theta_mods, double **divergences, double **timescale, long **initsampling, long **stagesampling, double **geoglikes, int desc[], int notu, int freqs, double *stratlikes, int areas, int f1, double iota, int foote);	void dateclade_stratlikelihood_ape_pulley(int HTU, long *ape, double **divergences, int notu, int nodes, unsigned long *ancestral, double ***thetas, double **theta_mods, double **timescale, long **initsampling, long **stagesampling, double **geoglikes, int freqs, double *stratlikes, int stages, int areas, double iota);	long *datecladefile(char *origins, long **tree, int clades, int notu);	double *datecladerealfile(char *Origins, long **tree, int clades, int notu);	void datecladeadd(long *fa, long **tree, int clades, int notu);	long *datecladeextra(int **Diverg, long *fa, long **tree, int clades, int notu);	long *datecladediverge(long **DT, long *fa, long **tree, int clades, int notu);	long *blbinnaive(long **tree, long *fas, int clades, int notu);	double *blmanaive(long **tree, double *fas, int clades, int notu);	double *branchlngma(long **tree, double **ranges, int *bl, int notu, int sign);	double **branchlength_ma_to_draw(long **tree, double **ranges, int *bl, int notu, int sign, double min);	void CladeSurviveAdd(long *la, long **tree, int clades, int notu);	long *CladeSurvive(long *la, long **tree, int clades, int notu);	double *temporalbranchlength(long **tree, double *fa, int clades, int notu);	double **temporalbranchlengthequal(long **tree, long *fba, double *byr, int clades, int notu);	double *temporalbranchlengthApo(long **tree, long *fa, int *bl, int clades, int notu);	double *AdjustTemporalBlap(long **tree, long *fa, int *bl, double *BT,double rate, int clades, int notu);	int naivestratdebt(long **tree, long *fa, int notu);	int stratdebt(long **tree, long**ranges, int *bl, int notu);	long *stratdebtnodes(long **tree, long**ranges, int *bl, int notu);	double stratdebtreal(long **tree, double **ranges, int *bl, int notu);	double calcSCI(long **tree, long *fa, int notu);	double calcSCIm(long **tree, long *fa, int notu);	double calcGER(long **tree, long *fa, int notu);	double calcRCI(long **tree, long**ranges, int notu);	double calcMSM(long **tree, long *fa, int notu);	double calcGERcor(long **tree, long**ranges, int *bl, int notu);	double calcRCIcor(long **tree, long**ranges, int *bl, int notu);	double calcMSMcor(long **tree, long**ranges, int *bl, int notu);	unsigned long *stratcompatibility (long **ranges, long **matrix, int *types, int nchars, int notu, int UNKNOWN, int INAP);	unsigned long *stratcompatfull(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	unsigned long *stratcompatfullplus(long **ranges, long **matrix, int *states, unsigned long *charcomps, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	double *stratcompatfullplusplus(long **ranges, long **matrix, int *states, unsigned long *charcomps, int ch1, int ch2, int notu, int tiebreaker, int UNKNOWN, int INAP);	double *stratcompatfullplusplusplus(long **ranges, long **matrix, int *states, unsigned long *charcomps, int ch1, int ch2, int notu, int onset, int term, int tiebreaker, int UNKNOWN, int INAP);	double *stratcompatfullplusdisparity(long **ranges, long **matrix, int *states, unsigned long *charcomps, int ch1, int ch2, int notu, int onset, int term, int tiebreaker, int UNKNOWN, int INAP);	long **statepairranges (long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	long **statepairfinds (long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	long **stateranges (long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP);	long **statefinds (long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP);	unsigned long stratconsiststate(long **staterng, long **chfnd, int s);	unsigned long stratconsistchar(long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP);	unsigned long stratconsiststpair(long **pairrng, long **pairfnd, int pr);	unsigned long stratconsistchpair(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	double stratfrconsistchpair(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	double stratfrconsistchar(long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP);	unsigned long *nu_stratcompat(unsigned long **charmat, unsigned long **compmat, long **ranges, int nstates, int notu, int strong, int hier, int UNKNOWN, int INAP);	void cleancladerangedata(long **ranges, int notu);	double *stratcompatparambstest(char taxonname[90], char citation[90], double *summary, double *mbl, long **omatrix, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, char excl, int tiebreaker, int debug, int UNKNOWN, int INAP);	double *stratcompatparambstest2(char taxonname[90], char citation[90], double *summary, double *mbl, long **omatrix, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, char excl, int tiebreaker, int debug, int UNKNOWN, int INAP);	double *stratcompatparambstest3(char taxonname[90], char citation[90], double *summary, double *mbl, long **omatrix, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, char excl, int tiebreaker, int debug, int UNKNOWN, int INAP);	double stratlikelihoodbranch(double fa, double origin, double *r, double *vr);	double gaplikelihoodgeog (double top, double bottom, double **preservation, double *timescale, double *geogpr, int stages, int areas, int rates);	void gaplikelihoodgeogsep (double **geogstrlnl, double top, double bottom, double **preservation, double *timescale, int stages, int areas, int rates, int taxon);	double **datebranchestodraw(long **tree, double **ranges, int *bl, int notu, int clades);	long **dateallbranchesfull(long **ranges, long **tree, long *bl, int clades, int notu);	int	*genstratcompat(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	double	propgenstratcompat(long **ranges, long **matrix, int *states, unsigned long **compmatrix, int notu, int nchars, int UNKNOWN, int INAP);//	double *stratcompatparambstest(char taxonname[90], char citation[90], double *summary, double *mbl, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, int debug, int UNKNOWN, int INAP);	#else	extern long *dateclade(long *fa, long **tree, int clades, int notu);	extern long **datecladefull(long **range, long **tree, int clades, int notu);	extern long **datecladesim(long **tree, int clades, int notu);	extern double *datecladereal(double *fa, long **tree, int clades, int notu, int sign);	extern void datecladerealape(long *ape, double **ranges, double **divergences, int nodes, int notu, unsigned long *ancestral);	extern void datecladerealape_pulley(int HTU, long *ape, double **divergences, int nodes, int notu, unsigned long *ancestral);	extern void datecladerealape_one_node(int HTU, int desc[], double **divergences, int f1, int notu, unsigned long *ancestral);	extern double sum_divergence_times(double **divergences, int notu, int nodes);	extern void datecladerealape_addexp(long *ape, double **ranges, double **divergences, int nodes, int notu, double lambda, unsigned long *ancestral);	extern long *datecladefile(char *origins, long **tree, int clades, int notu);	extern double *datecladerealfile(char *Origins, long **tree, int clades, int notu);	extern void datecladeadd(long *fa, long **tree, int clades, int notu);	extern long *datecladeextra(int **Diverg, long *fa, long **tree, int clades, int notu);	extern long *datecladediverge(long **DT, long *fa, long **tree, int clades, int notu);	extern long *blbinnaive(long **tree, long *fas, int clades, int notu);	extern double *blmanaive(long **tree, double *fas, int clades, int notu);	extern double *branchlngma(long **tree, double **ranges, int *bl, int notu, int sign);	extern double **branchlength_ma_to_draw(long **tree, double **ranges, int *bl, int notu, int sign, double min);	extern void CladeSurviveAdd(long *la, long **tree, int clades, int notu);	extern long *CladeSurvive(long *la, long **tree, int clades, int notu);	extern double *temporalbranchlength(long **tree, double *fa, int clades, int notu);	extern double **temporalbranchlengthequal(long **tree, long *fba, double *byr, int clades, int notu);	extern double *temporalbranchlengthApo(long **tree, long *fa, int *bl, int clades, int notu);	extern double *AdjustTemporalBlap(long **tree, long *fa, int *bl, double *BT,double rate, int clades, int notu);	extern int naivestratdebt(long **tree, long *fa, int notu);	extern int stratdebt(long **tree, long**ranges, int *bl, int notu);	extern long *stratdebtnodes(long **tree, long**ranges, int *bl, int notu);	extern double stratdebtreal(long **tree, double **ranges, int *bl, int notu);	extern double calcSCI(long **tree, long *fa, int notu);	extern double calcSCIm(long **tree, long *fa, int notu);	extern double calcGER(long **tree, long *fa, int notu);	extern double calcRCI(long **tree, long**ranges, int notu);	extern double calcMSM(long **tree, long *fa, int notu);	extern double calcGERcor(long **tree, long**ranges, int *bl, int notu);	extern double calcRCIcor(long **tree, long**ranges, int *bl, int notu);	extern double calcMSMcor(long **tree, long**ranges, int *bl, int notu);	extern unsigned long *stratcompatibility (long **ranges, long **matrix, int *types, int nchars, int notu, int UNKNOWN, int INAP);	extern unsigned long *stratcompatfull(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	extern unsigned long *stratcompatfullplus(long **ranges, long **matrix, int *states, unsigned long *charcomps, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	extern double *stratcompatfullplusplus(long **ranges, long **matrix, int *states, unsigned long *charcomps, int ch1, int ch2, int notu, int tiebreaker, int UNKNOWN, int INAP);	extern double *stratcompatfullplusplusplus(long **ranges, long **matrix, int *states, unsigned long *charcomps, int ch1, int ch2, int notu, int onset, int term, int tiebreaker, int UNKNOWN, int INAP);	extern double *stratcompatfullplusdisparity(long **ranges, long **matrix, int *states, unsigned long *charcomps, int ch1, int ch2, int notu, int onset, int term, int tiebreaker, int UNKNOWN, int INAP);	extern long **statepairranges (long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	extern long **statepairfinds (long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	extern long **stateranges (long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP);	extern long **statefinds (long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP);	extern unsigned long stratconsiststate(long **staterng, long **chfnd, int s);	extern unsigned long stratconsistchar(long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP);	extern unsigned long stratconsiststpair(long **pairrng, long **pairfnd, int pr);	extern unsigned long stratconsistchpair(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	extern double stratfrconsistchpair(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	extern double stratfrconsistchar(long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP);	extern unsigned long *nu_stratcompat(unsigned long **charmat, unsigned long **compmat, long **ranges, int nstates, int notu, int strong, int hier, int UNKNOWN, int INAP);	extern void cleancladerangedata(long **ranges, int notu);	extern double *stratcompatparambstest(char taxonname[90], char citation[90], double *summary, double *mbl, long **omatrix, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, char excl, int tiebreaker, int debug, int UNKNOWN, int INAP);	extern double *stratcompatparambstest2(char taxonname[90], char citation[90], double *summary, double *mbl, long **omatrix, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, char excl, int tiebreaker, int debug, int UNKNOWN, int INAP);	extern double *stratcompatparambstest3(char taxonname[90], char citation[90], double *summary, double *mbl, long **omatrix, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, char excl, int tiebreaker, int debug, int UNKNOWN, int INAP);	extern double stratlikelihoodbranch(double fa, double origin, double *r, double *vr);	extern double gaplikelihoodgeog (double top, double bottom, double **preservation, double *timescale, double *geogpr, int stages, int areas, int rates);	extern void gaplikelihoodgeogsep (double **geogstrlnl, double top, double bottom, double **preservation, double *timescale, int stages, int areas, int rates, int taxon);	extern double **datebranchestodraw(long **tree, double **ranges, int *bl, int notu, int clades);	extern long **dateallbranchesfull(long **ranges, long **tree, long *bl, int clades, int notu);	extern int	*genstratcompat(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);	extern double	propgenstratcompat(long **ranges, long **matrix, int *states, unsigned long **compmatrix, int notu, int nchars, int UNKNOWN, int INAP);//	extern double *stratcompatparambstest(char taxonname[90], char citation[90], double *summary, double *mbl, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, int debug, int UNKNOWN, int INAP);	extern void dateclade_stratlikelihood_ape_addexp(long *ape, double **divergences, int nodes, int notu, double lambda, unsigned long *ancestral, unsigned long *minbr, double **ranges, double ***thetas, double **theta_mods, int freqs, double *stratlikes, double **timescale, long **initsampling, long **stagesampling, double **geoglikes, int areas, double iota);	extern void dateclade_stratolikelihood_one_node(int HTU, double ***thetas, double **theta_mods, double **divergences, double **timescale, long **initsampling, long **stagesampling, double **geoglikes, int desc[], int notu, int freqs, double *stratlikes, int areas, int f1, double iota, int foote);	extern 	void dateclade_stratlikelihood_ape_pulley(int HTU, long *ape, double **divergences, int notu, int nodes, unsigned long *ancestral, double ***thetas, double **theta_mods, double **timescale, long **initsampling, long **stagesampling, double **geoglikes, int freqs, double *stratlikes, int stages, int areas, double iota);	#endif