#ifdef mcmc#include <stdlib.h>#include <stdio.h>#include <time.h>#include <math.h>#include <string.h>#include "probability.h"//#include <gsl/gsl_rng.h>//#include <gsl/gsl_randist.h>//int GibbsSample (int argc, char *argv);double **GibbsSample (int n, double rho);void swapbr (long *tree, int ttlbr, int notu);void swapbrX(long *ape, int ttlbr, int notu, int br);void swapbr2 (long tree[], int ttlbr, int notu);void changebd(long **tree, double *bd, int nodes);double exponentialbranchdurations (double lambda);void swap_and_recalculate(long *ape, int notu, int nodes, long **statechars, long *charswstates, int maxst, double *alphas, double *betas, int rates, double ***statelikes, double lambda, double **ranges, double **divergences, unsigned long *ancestral, unsigned long *minbr);void change_branch_length(long *ape, int notu, int nodes, long **statechars, long *charswstates, int nchars, int maxst, double *alphas, double *betas, int rates, double ***statelikes, double lambda, double **ranges, double **divergences, unsigned long *ancestral, unsigned long *minbr);void collapse_node_and_recalculate(long *ape, int notu, int *nodes, long **statechars, long *charswstates, int nchars, int maxst, double *alphas, double *betas, int rates, double ***statelikes, double lambda, double **ranges, double **divergences, unsigned long *ancestral, unsigned long *minbr);//void collapse_node_and_recalculate_strat(long *ape, int notu, int *nodes, long **statechars, long *charswstates, int nchars, int maxst, double *alphas, double *betas, int rates, double ***statelikes, double lambda, double **ranges, double **divergences, unsigned long *ancestral, unsigned long *minbr, double **geoglikes, int areas);void collapse_node_and_recalculate_strat(long *ape, int notu, int *nodes, long **statechars, long *charswstates, int nchars, int maxst, double *alphas, double *betas, int rates, double ***statelikes, double lambda, double **ranges, double **divergences, unsigned long *ancestral, unsigned long *minbr);void change_char_rate(double *alphas, int rates);#else//extern int GibbsSample (int argc, char *argv);extern double **GibbsSample (int n, double rho);extern void swapbr (long * tree, int ttlbr, int notu);extern void swapbrX(long *ape, int ttlbr, int notu, int br);extern void swapbr2 (long [], int ttlbr, int notu);extern void changebd(long **tree, double *bd, int nodes);extern double exponentialbranchdurations (double lambda);extern void swap_and_recalculate(long *ape, int notu, int nodes, long **statechars, long *charswstates, int maxst, double *alphas, double *betas, int rates, double ***statelikes, double lambda, double **ranges, double **divergences, unsigned long *ancestral, unsigned long *minbr);extern void change_branch_length(long *ape, int notu, int nodes, long **statechars, long *charswstates, int nchars, int maxst, double *alphas, double *betas, int rates, double ***statelikes, double lambda, double **ranges, double **divergences, unsigned long *ancestral, unsigned long *minbr);extern void collapse_node_and_recalculate(long *ape, int notu, int *nodes, long **statechars, long *charswstates, int nchars, int maxst, double *alphas, double *betas, int rates, double ***statelikes, double lambda, double **ranges, double **divergences, unsigned long *ancestral, unsigned long *minbr);//extern void collapse_node_and_recalculate_strat(long *ape, int notu, int *nodes, long **statechars, long *charswstates, int nchars, int maxst, double *alphas, double *betas, int rates, double ***statelikes, double lambda, double **ranges, double **divergences, unsigned long *ancestral, unsigned long *minbr, double **geoglikes, int areas);extern void collapse_node_and_recalculate_strat(long *ape, int notu, int *nodes, long **statechars, long *charswstates, int nchars, int maxst, double *alphas, double *betas, int rates, double ***statelikes, double lambda, double **ranges, double **divergences, unsigned long *ancestral, unsigned long *minbr);extern void change_char_rate(double *alphas, int rates);#endif