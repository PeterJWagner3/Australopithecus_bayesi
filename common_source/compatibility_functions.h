#ifdef compatibility_functions#include <stdlib.h>#include <stdio.h>#include <time.h>#include <math.h>#include <string.h>#define distribution_fits#include "distribution_fits.h"#define filereading#include "filereading.h"#define matrixanalysis#include "matrixanalysis.h"#define matrixchange#include "matrixchange.h"#define matrixreading#include "matrixreading.h"#define memory#include "memory.h"//#define minmax#include "minmax.h"#define MonteCarloPhylogenyFunctions#include "MonteCarloPhylogenyFunctions.h"#define TreeRead#include "tree_read.h"int charcombs(long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);int *charcombos(long *t1, long *t2, int st, int notu, int UNKNOWN, int INAP);int findswingpair(long **chpairs, int pairs, int maxstates);int unordcompatible(int *t1, int *t2, int st1, int st2, int notu, int UNKNOWN, int INAP);int ordcompatible(int *t1, int *t2, int st1, int st2, int notu, int UNKNOWN, int INAP);int binarycompatible(int *t1, int *t2, int notu, int UNKNOWN, int INAP);unsigned long **compatible(int *states, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);void **makecompatibilitymatrix(unsigned long ** comatrix, int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);void **remakecompatibilitymatrixforonechar(int onechar, unsigned long ** comatrix, int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);int nu_comp(int *states, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);int nu_comp_rem(int ch1, int ch2, int *states, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);unsigned long *char_comp(int *states, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);int countcomp(unsigned long **comatrix, int nchars);void *countcharcomps(unsigned long **comatrix, int nchars, unsigned long *charcomps);unsigned long **mutualcomp(unsigned long **comatrix, int nchars);int possiblemutualcompatibilities(unsigned long **comatrix, unsigned long *charcomps, int nchars);double *propmutualcomp(unsigned long **comatrix, unsigned long *charcomps, int *autaps, int *nstates, int nchars, char excl);int char_nu_comp(int ch, int *states, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);int pair_comp(int ch1, int ch2, int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);long *PossibleCharCombinations(int *states, int nchars);long **charcombinations(long **combomatrix, int *nstates, int notu, long **chmatrix, int nchars, int UNKNOWN, int INAP);double *likeallstepsgivencomp(long **matrix, int *ctype, int *nstates, int nchars, int notu, int empcomp, int comptype, int fossils, double *mbl, int *bias, int *maxd, int pars, int debug, int UNKNOWN, int INAP);double **probcompgivenallsteps(long **matrix, int *ctype, int *nstates, int nchars, int notu, int empcomp, int comptype, int fossils, double *mbl, int *bias, int *maxd, int pars, int max, int debug, int UNKNOWN, int INAP);double **likelystepsperchar(long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);double **likelystepsperch(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);/*double **likelystepspercharacter(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int *ststeps, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);	*/int statepairs(long **mat, int notu, int nchar, int maxstate/*, int UNKNOWN, int INAP*/);int **obsstatepairs(int *t1, int *t2, int *combos, int st, int notu);void likelystepspercharacter(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);void **likelypoisratessperch(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);void **likelygammadist(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);void **likelyvariablerates(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);void **likelyvarratesremainder(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);void **likelyvarratesremainderinap(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);double **likelyvarratesremnooutput(/*char *taxonname, */long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP, int truns, int mruns, int cruns);double *mutualcompatibilitytest(double *summary, double *mbl, long **omatrix, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, char excl, int basis, int debug, int UNKNOWN, int INAP);#elseextern int charcombs(long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP);extern int *charcombos(long *t1, long *t2, int st, int notu, int UNKNOWN, int INAP);extern int findswingpair(long **chpairs, int pairs, int maxstates);extern int unordcompatible(int *t1, int *t2, int st1, int st2, int notu, int UNKNOWN, int INAP);extern int ordcompatible(int *t1, int *t2, int st1, int st2, int notu, int UNKNOWN, int INAP);extern int binarycompatible(int *t1, int *t2, int notu, int UNKNOWN, int INAP);extern unsigned long **compatible(int *states, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);extern void **makecompatibilitymatrix(unsigned long ** comatrix, int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);extern void **remakecompatibilitymatrixforonechar(int onechar, unsigned long ** comatrix, int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);extern int nu_comp(int *states, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);extern int nu_comp_rem(int ch1, int ch2, int *states, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);extern unsigned long *char_comp(int *states, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);extern int countcomp(unsigned long **comatrix, int nchars);extern void *countcharcomps(unsigned long **comatrix, int nchars, unsigned long *charcomps);extern unsigned long **mutualcomp(unsigned long **comatrix, int nchars);extern int possiblemutualcompatibilities(unsigned long **comatrix, unsigned long *charcomps, int nchars);extern double *propmutualcomp(unsigned long **comatrix, unsigned long *charcomps, int *autaps, int *nstates, int nchars, char excl);extern int char_nu_comp(int ch, int *states, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);extern int pair_comp(int ch1, int ch2, int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);extern long *PossibleCharCombinations(int *states, int nchars);extern long **charcombinations(long **combomatrix, int *nstates, int notu, long **chmatrix, int nchars, int UNKNOWN, int INAP);extern double *likeallstepsgivencomp(long **matrix, int *ctype, int *nstates, int nchars, int notu, int empcomp, int comptype, int fossils, double *mbl, int *bias, int *maxd, int pars, int debug, int UNKNOWN, int INAP);extern double **probcompgivenallsteps(long **matrix, int *ctype, int *nstates, int nchars, int notu, int empcomp, int comptype, int fossils, double *mbl, int *bias, int *maxd, int pars, int max, int debug, int UNKNOWN, int INAP);extern double **likelystepsperchar(long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);extern double **likelystepsperch(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);/*extern double **likelystepspercharacter(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int *ststeps, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);	*/extern int statepairs(long **mat, int notu, int nchar, int maxstate/*, int UNKNOWN, int INAP*/);extern int **obsstatepairs(int *t1, int *t2, int *combos, int st, int notu);/*extern double **likelystepsperch(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);	*/extern void likelystepspercharacter(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);extern void **likelypoisratessperch(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);extern void **likelygammadist(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);extern void **likelyvariablerates(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);extern void **likelyvarratesremainder(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);extern void **likelyvarratesremainderinap(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP);extern double **likelyvarratesremnooutput(/*char *taxonname, */long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP, int truns, int mruns, int cruns);extern double *mutualcompatibilitytest(double *summary, double *mbl, long **omatrix, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, char excl, int basis, int debug, int UNKNOWN, int INAP);#endif