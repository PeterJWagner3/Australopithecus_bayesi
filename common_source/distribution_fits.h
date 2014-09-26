/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
 *	Peter Wagner: pwagner@fmnh.org
 *	Matthew Kosnik: mkosnik@uchicago.edu
 *
 *	This file is copyright (C) 2002-2009 Peter Wagner, Matthew Kosnik
 *
 *	This program is free software; you can redistribute it and/or modify it 
 *	under the terms of version 2 the GNU General Public License as published 
 *	by the Free Software Foundation.
 *
 *	This program is distributed in the hope that it will be useful, but WITHOUT
 *	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 *	FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 *	more details.
 *
 *	To view a copy of the license go to:
 *	http://www.fsf.org/copyleft/gpl.html
 *	To receive a copy of the GNU General Public License write the Free Software
 * 	Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *	Copies of this source code are available without cost from:
 *	http://geosci.uchicago.edu/paleo/csource/
 *	
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/***********************************************************************
 v0.0 	2001.05.04 BY P.J.WAGNER III
 CODE WRITTEN - CODEWARRIOR 6.1 / MacOS 9.04
 v0.1	2001.05.04 BY M. KOSNIK
 FUNCTIONS MADE USING WAGNER'S CODE
 VARIABLE NAMES CHANGED TO MAKE THEM MORE UNDERSTANDABLE
 COMMENTS ADDED
 v0.2	2001.07.30 BY M. KOSNIK
 ADDITIONAL COMMENTS...
 REMOVAL OF UNUSED CODE....
 FILE REORGANIZATION....
 CLEANED SYNTAX TO COMPILE WITH PROJECT BUILDER IN OS X
 Apple Computer, Inc. version gcc-926, based on gcc version 2.95.2 19991024 (release)
 v0.3	2001.09.03 BY M. KOSNIK
 ADDITIONAL COMMENTS...
 REMOVAL OF UNUSED CODE....
 CREATED double *ranks_w_ties(int *abundance, int num_taxa);
 FIXED:
 - array overflow errors linked to while() loops
 
 v. 0.5 2001.12.25 BY M. KOSNIK
 massive revissions to code.
 renamed all functions.
 searches to user defined number of decimal points.
 - steps through until fit is worse
 - then steps through at next decimal until fit is worse and repeats.
 - does check to make sure peak is not passed
 
 v. 0.6 2001.12.27 BY M. KOSNIK
 major revissions to code.
 - move contants to header file
 - revise optimization incrementing
 v. 0.7 2002.01.09 BY M. KOSNIK
 major revissions to log-normal code.
 - rewrite log normal to use revised wagner log-normal generation code 
 v. 0.8 2002.02.07 BY M. KOSNIK
 major revisions to gs / zm code.
 - rewrite to go backward when passed peaks
 - search to minimum change in support
 v. 0.9 2002.07.10 BY P. WAGNER
 create new likelihood function, move old one to "foote" function
 v. 0.10 2002.07.14 BY M. KOSNIK
 create new aic,aic_c,bic functions.
 v. 0.11 2002.07.15 BY M. KOSNIK
 bug fix... at certain richnesses (50-105) slope was exploring values less than 1.
 these values are clearly out of the range we should be exploring.
 limited the subtract ei option in the gs / zm loops to for be used in cases
 where e-ei>emin. Seems to solve problem. Left in code to print error message if
 we explore values less than emin.
 v. 0.11 2002.07.29 BY M. KOSNIK
 cleaned up formatting.
 v. 0.12 2002.07.29 BY M. KOSNIK
 added fit log-series.
 added fit log-power.
 v. 0.13 2002.07.31 BY M. KOSNIK
 log-power produces good results at excessive sample size (10000).
 fit coef & exp are both highly correlated with sim coef & exp.
 log-series produces reasonable results at excessive sample size.
 v. 0.14 2003.11.14 BY P. WAGNER
 created new likelihood functions based on octave logic.  These determine the expected
 proportions of species with 0�max finds given hypothesized distribution and 
 richness.  Log-likelihood based on #with x finds x log expected proportion.  This
 includes the hypothesized unsampled taxa, too.  
 v. 0.15 2004.01.07 BY P. WAGNER
 implemented Poisson test based on octave logic, looking at observed and expected numbers 
 of taxa with 1, 2, 3, etc. finds.
 v. 0.16 2004.07.30 BY P. WAGNER
 removed "faux" log-normal routines and replaced these with routines that calculate
 the log-normal by dividing the area under the normal curve into X+1 equal units.
 The position of partition J on the X-axis become the relative abundance of taxon J
 on the logged-Y axis of a Whittaker plot.  
 v. 0.17 2005.01.02 BY P. WAGNER
 added Hubbell's zero sum multinomial, using Olzewski's routine.
 v. 0.18 2005.01.28 BY P. WAGNER
 added routines to put support �error� bars on richness.
 v. 0.19 2005.03.12 BY P. WAGNER
 added separate routines for using Poisson and multinomial distributions.
 *********************************************************************************************/

#define SUPINC 0.05f				/* minimum change in support 	*/
#define FITINC 1					/* increment to fit by			*/

#ifdef distribution_fits
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "matrixchange.h"
#include "memory.h"
#include "sort.h"
#include "distribution_calc.h"
#include "historicaldiversity.h"

#define ZERO 0.00001f
#define	e	2.718281828

/* log normal parameters	*/
#define NSTDV 4				/* number of octaves per standard deviation */
#define OMIN 4				/* min number of octaves = num octaves		*/

/*	double calc_aic_c(int parameters, double fit, int measurements); 
 double calc_aic(int parameters, double fit); 
 double calc_bic(int parameters, double fit, int measurements); 	*/

double calc_likelihood(double *A, int *S, int sizeofS, int Samplesize);
double calc_likelihood_Foote(double *A, int *S, int sizeofS);
double calc_likelihood_exp(double *X, long *nufinds, long *histogram, int unqfnds);

double *count_fit_mul_uni(int *empdist, int ntaxa, int nspec);
double *count_fit_mul_geo(int *empdist, int ntaxa, int nspec);
double *count_fit_mul_geos(int *empdist, int ntaxa, int nspec);
double *count_fit_mul_zipf(int *empdist, int ntaxa, int nspec);
double *count_fit_mul_lgn(int *empdist, int ntaxa, int nspec);
double *count_fit_mul_lnt(int *empdist, int ntaxa, int nspec);
double *count_fit_mul_zsm(int *empdist, int ntaxa, int nspec);
double *count_fit_mul_zsm_test(int *empdist, int ntaxa, int nspec);
double count_fit_mul_best(int *empdist, int ntaxa);

double **support_bars_mul_unif(double bestS, int ntaxa, int bestR, int *empdist, double bar);
double **support_bars_mul_geo(double bestS, double bestM, int ntaxa, int *empdist, double bar);
double **support_bars_mul_geos(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
double **support_bars_mul_zipf(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
double **support_bars_mul_lgn(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
double **support_bars_mul_zsm(double besttheta, double bestm, double bestS, int ntaxa, int *empdist, double bar);
double **support_bars_mul_zipf_full(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
double **support_bars_mul_lgn_full(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);

double *count_fit_poi_uni(int *empdist, int ntaxa, int nspec);
double *count_fit_poi_geo(int *empdist, int ntaxa, int nspec);
double *count_fit_poi_geos(int *empdist, int ntaxa, int nspec);
double *count_fit_poi_zipf(int *empdist, int ntaxa, int nspec);
double *count_fit_poi_lgn(int *empdist, int ntaxa, int nspec);
double *count_fit_poi_lnt(int *empdist, int ntaxa, int nspec);
double *count_fit_poi_zsm(int *empdist, int ntaxa, int nspec);
double *count_fit_poi_zsm_test(int *empdist, int ntaxa, int nspec);
double count_fit_mul_zsm_exact(int *empdist, int ntaxa, int nspec, double theta, double m);
double count_fit_poi_best(int *empdist, int ntaxa);

double **support_bars_poi_unif(double bestS, int ntaxa, int bestR, int *empdist, double bar);
double **support_bars_poi_geo (double bestS, double bestM, int ntaxa, int *empdist, double bar);
double **support_bars_poi_geos(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
double **support_bars_poi_zipf(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
double **support_bars_poi_lgn(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
double **support_bars_poi_zsm(double besttheta, double bestm, double bestS, int ntaxa, int *empdist, double bar, double omega);

double *wht_fit_gs(int *empdist, int ntaxa, int nspec);
double *wht_fit_zm(int *empdist, int ntaxa, int nspec);
double *wht_fit_lnt(int *empdist, int ntaxa, int nspec);

double *expfinds(double *dist, int S, int mxfds, int ttfds);
double *expoccurrences(double *dist, int S, int posfds);
double *expfindspart(double *dist, int S, long *nofds, int unfds, int ttfds);
double *expfindsprop(double *dist, int S, int mxfds, int t);
double *expfindsfromexpind(double *oct, int mxfds, int ttl);
double *expfindsfromexpindprop (double *oct, int mxfds, int ttl);
double *expfindsln(double mag, int S, int mxfds, int t);

double **richness_support_bars_mul_geo(double ubS, double lbS, double inc, double minf, int ntaxa, int *empdist);

double *fuzz_fit_mul_uniform(double *dist, int obsr, int mxstp);
double *fuzz_fit_mul_onepgamma(double *dist, int obsr, int mxstp, double base, int pts);
double *fuzz_fit_mul_twopgamma(double *dist, int obsr, int mxstp, double base, int pts);
double *fuzz_fit_mul_lognormal(double *dist, int obsr, int mxstp, double base, int pts);
double *fuzz_fit_mul_twoplgnorm(double *dist, int obsr, int mxstp, double base, int pts);
double *fuzz_fit_mul_unif2param(double *dist, int obsr, int mxstp);

double occurrence_fit_mul_best(int *empdist, int ntaxa);
double *occurrence_fit_mul_uni(int *empdist, int ntaxa, int ncoll);
double *occurrence_fit_mul_onepgamma(int *empdist, int ntaxa, int ncoll);
double *occurrence_fit_mul_geo(int *empdist, int ntaxa, int ncoll);
double *occurrence_fit_mul_lgn(int *empdist, int ntaxa, int ncoll);
double *occurrence_fit_mul_zipf(int *empdist, int ntaxa, int ncoll);

double *rate_fit_mul_uni(int *empdist, int ntaxa, int ncoll);
double *rate_fit_mul_geo(int *empdist, int ntaxa, int ncoll);
double *rate_fit_mul_lgn(int *empdist, int ntaxa, int ncoll);
double *rate_fit_mul_Zipf(int *empdist, int ntaxa, int ncoll);
double *rate_fit_mul_gammaone(int *empdist, int ntaxa, int ncoll);

double rate_fit_Pois_best(int *empdist, int ntaxa, int ncoll);
double *rate_fit_Pois_uni(int *empdist, int ntaxa, int ncoll);
double *rate_fit_Pois_geo(int *empdist, int ntaxa, int ncoll);
double *rate_fit_Pois_lgn(int *empdist, int ntaxa, int ncoll);
double *rate_fit_Pois_Zipf(int *empdist, int ntaxa, int ncoll);
double *rate_fit_Pois_gammaone(int *empdist, int ntaxa, int ncoll);

double **modelshared(double *tad, int S1, int S2, int S1S2, int N);
double *sharedgivenTADs(unsigned long *shad1, unsigned long *shad2, int shared, int model1, int model2, int maxtest, int ncoll1, int ncoll2);
double *sharedgivenTADs_supportbars(unsigned long *shad1, unsigned long *shad2, int shared, int model1, int model2, int maxS, double supportbar, int ncoll1, int ncoll2);
double *mul_lgn_at_set_S_occur(int *empdist, int ntaxa, int Sh, int ncoll);
double *mul_geo_at_set_S_occur(int *empdist, int Sobs, int Shyp, int ncoll);

double *sharedrichnesssupport(int sSobs, int sSmax, int binA, int binB, int possFindsA, int possFindsB, int modelA, int modelB, int richA, int richB, double **params);

#else

/*	extern double calc_aic_c(int parameters, double fit, int measurements); 
 extern double calc_aic(int parameters, double fit); 
 extern double calc_bic(int parameters, double fit, int measurements); */

extern double calc_likelihood(double *A, int *S, int sizeofS, int Samplesize);
extern double calc_likelihood_Foote(double *A, int *S, int sizeofS);
extern double calc_likelihood_exp(double *X, long *nufinds, long *histogram, int unqfnds);

extern double *count_fit_mul_uni(int *empdist, int ntaxa, int nspec);
extern double *count_fit_mul_geo(int *empdist, int ntaxa, int nspec);
extern double *count_fit_mul_geos(int *empdist, int ntaxa, int nspec);
extern double *count_fit_mul_zipf(int *empdist, int ntaxa, int nspec);
extern double *count_fit_mul_lgn(int *empdist, int ntaxa, int nspec);
extern double *count_fit_mul_lnt(int *empdist, int ntaxa, int nspec);
extern double *count_fit_mul_zsm(int *empdist, int ntaxa, int nspec);
extern double *count_fit_mul_zsm_test(int *empdist, int ntaxa, int nspec);
extern double count_fit_mul_zsm_exact(int *empdist, int ntaxa, int nspec, double theta, double m);
extern double count_fit_mul_best(int *empdist, int ntaxa);

extern double **support_bars_mul_unif(double bestS, int ntaxa, int bestR, int *empdist, double bar);
extern double **support_bars_mul_geo (double bestS, double bestM, int ntaxa, int *empdist, double bar);
extern double **support_bars_mul_geos(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
extern double **support_bars_mul_zipf(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
extern double **support_bars_mul_lgn(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
extern double **support_bars_mul_zsm(double besttheta, double bestm, double bestS, int ntaxa, int *empdist, double bar);
extern double **support_bars_mul_zipf_full(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
extern double **support_bars_mul_lgn_full(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);

extern double *count_fit_poi_uni(int *empdist, int ntaxa, int nspec);
extern double *count_fit_poi_geo(int *empdist, int ntaxa, int nspec);
extern double *count_fit_poi_geos(int *empdist, int ntaxa, int nspec);
extern double *count_fit_poi_zipf(int *empdist, int ntaxa, int nspec);
extern double *count_fit_poi_lgn(int *empdist, int ntaxa, int nspec);
extern double *count_fit_poi_lnt(int *empdist, int ntaxa, int nspec);
extern double *count_fit_poi_zsm(int *empdist, int ntaxa, int nspec);
extern double *count_fit_poi_zsm_test(int *empdist, int ntaxa, int nspec);
extern double *count_fit_poi_zsm_test_exact(int *empdist, int ntaxa, int nspec, double theta, double m);
extern double count_fit_poi_best(int *empdist, int ntaxa);

extern double **support_bars_poi_unif(double bestS, int ntaxa, int bestR, int *empdist, double bar);
extern double **support_bars_poi_geo (double bestS, double bestM, int ntaxa, int *empdist, double bar);
extern double **support_bars_poi_geos(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
extern double **support_bars_poi_zipf(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
extern double **support_bars_poi_lgn(double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
extern double **support_bars_poi_zsm(double besttheta, double bestm, double bestS, int ntaxa, int *empdist, double bar, double omega);

extern double *wht_fit_gs(int *empdist, int ntaxa, int nspec);
extern double *wht_fit_zm(int *empdist, int ntaxa, int nspec);
extern double *wht_fit_lnt(int *empdist, int ntaxa, int nspec);

extern double *expfinds(double *dist, int S, int mxfds, int ttfds);
extern double *expoccurrences(double *dist, int S, int posfds);
extern double *expfindspart(double *dist, int S, long *nofds, int unfds, int ttfds);
extern double *expfindsprop(double *dist, int S, int mxfds, int t);
extern double *expfindsfromexpind(double *oct, int mxfds, int ttl);
extern double *expfindsfromexpindprop (double *oct, int mxfds, int ttl);
extern double *expfindsln(double mag, int S, int mxfds, int t);

extern double **richness_support_bars_mul_geo(double ubS, double lbS, double inc, double minf, int ntaxa, int *empdist);

extern double *fuzz_fit_mul_uniform(double *dist, int obsr, int mxstp);
extern double *fuzz_fit_mul_onepgamma(double *dist, int obsr, int mxstp, double base, int pts);
extern double *fuzz_fit_mul_twopgamma(double *dist, int obsr, int mxstp, double base, int pts);
extern double *fuzz_fit_mul_lognormal(double *dist, int obsr, int mxstp, double base, int pts);
extern double *fuzz_fit_mul_twoplgnorm(double *dist, int obsr, int mxstp, double base, int pts);
extern double *fuzz_fit_mul_unif2param(double *dist, int obsr, int mxstp);

extern double occurrence_fit_mul_best(int *empdist, int ntaxa);
extern double *occurrence_fit_mul_onepgamma(int *empdist, int ntaxa, int ncoll);
extern double *occurrence_fit_mul_uni(int *empdist, int ntaxa, int ncoll);
extern double *occurrence_fit_mul_geo(int *empdist, int ntaxa, int ncoll);
extern double *occurrence_fit_mul_lgn(int *empdist, int ntaxa, int ncoll);
extern double *occurrence_fit_mul_zipf(int *empdist, int ntaxa, int ncoll);

extern double *rate_fit_mul_uni(int *empdist, int ntaxa, int ncoll);
extern double *rate_fit_mul_geo(int *empdist, int ntaxa, int ncoll);
extern double *rate_fit_mul_lgn(int *empdist, int ntaxa, int ncoll);
extern double *rate_fit_mul_Zipf(int *empdist, int ntaxa, int ncoll);

extern double rate_fit_Pois_best(int *empdist, int ntaxa, int ncoll);
extern double *rate_fit_Pois_uni(int *empdist, int ntaxa, int ncoll);
extern double *rate_fit_Pois_geo(int *empdist, int ntaxa, int ncoll);
extern double *rate_fit_Pois_lgn(int *empdist, int ntaxa, int ncoll);
extern double *rate_fit_Pois_Zipf(int *empdist, int ntaxa, int ncoll);

extern double **modelshared(double *tad, int S1, int S2, int S1S2, int N);
extern double *sharedgivenTADs(unsigned long *shad1, unsigned long *shad2, int shared, int model1, int model2, int maxtest, int ncoll1, int ncoll2);
extern double *sharedgivenTADs_supportbars(unsigned long *shad1, unsigned long *shad2, int shared, int model1, int model2, int maxS, double supportbar, int ncoll1, int ncoll2);
extern double *mul_lgn_at_set_S_occur(int *empdist, int ntaxa, int Sh, int ncoll);
extern double *mul_geo_at_set_S_occur(int *empdist, int Sobs, int Shyp, int ncoll);
extern double *sharedrichnesssupport(int sSobs, int sSmax, int binA, int binB, int possFindsA, int possFindsB, int modelA, int modelB, int richA, int richB, double **params);

#endif
