/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
/*	Peter Wagner: pwagner@fmnh.org
/*	Matthew Kosnik: mkosnik@uchicago.edu
/*
/*	This file is copyright (C) 2002-2005 Matthew Kosnik, Peter Wagner
/*
/*	This program is free software; you can redistribute it and/or modify it 
/*	under the terms of version 2 the GNU General Public License as published 
/*	by the Free Software Foundation.
/*
/*	This program is distributed in the hope that it will be useful, but WITHOUT
/*	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
/*	FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
/*	more details.
/*
/*	To view a copy of the license go to:
/*	http://www.fsf.org/copyleft/gpl.html
/*	To receive a copy of the GNU General Public License write the Free Software
/* 	Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
/*
/*	Copies of this source code are available without cost from:
/*	http://geosci.uchicago.edu/paleo/csource/
/*	
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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
			proportions of species with 0Émax finds given hypothesized distribution and 
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
*********************************************************************************************/

#define SUPINC 0.01f				/* minimum change in support 	*/
#define FITINC 1					/* increment to fit by			*/

#ifdef distribution_fits
	#include <math.h>
	#include <float.h>
	#include <stdlib.h>
	#include <stdio.h>
	#include "memory.h"
	#include "sort.h"
	#include "distribution_calc.h"
	#include "historicaldiversity.h"
	

	#define ZERO 0.00001f
	#define	e	2.718281828

	/* log normal parameters	*/
	#define NSTDV 4				/* number of octaves per standard deviation */
	#define OMIN 4				/* min number of octaves = num octaves		*/

	double calc_aic_c(int parameters, double fit, int measurements); 
	double calc_aic(int parameters, double fit); 
	double calc_bic(int parameters, double fit, int measurements); 

	double calc_likelihood(double *A, int *S, int sizeofS, int Samplesize);
	double calc_likelihood_Foote(double *A, int *S, int sizeofS);
	double calc_likelihood_exp(double *X, long *nufinds, long *histogram, int unqfnds);

	double *fit_gs(int *abundance, int ntaxa, int Samplesize);
	double *fit_zm(int *abundance, int ntaxa, int Samplesize);
	double *fit_ln(int *abundance, int ntaxa, int Samplesize);
	double *fit_fln(int *abundance, int ntaxa, int Samplesize);
	double *fit_ls(int *abundance, int ntaxa, int Samplesize);
	double *fit_lp(int *abundance, int ntaxa, int Samplesize);

	double *chi_fit_gs(int *empdist, int ntaxa, int nspec);
	
	double *poi_fit_uni(int *empdist, int ntaxa, int nspec);
	double *poi_fit_geo(int *empdist, int ntaxa, int nspec);
	double *poi_fit_geos(int *empdist, int ntaxa, int nspec);
	double *poi_fit_zipf(int *empdist, int ntaxa, int nspec);
	double *poi_fit_lnu(int *empdist, int ntaxa, int nspec);
	double *poi_fit_lnt(int *empdist, int ntaxa, int nspec);
	double *poi_fit_zsm(int *empdist, int ntaxa, int nspec);
	double *poi_fit_zsm_test(int *empdist, int ntaxa, int nspec);
	
	double *wht_fit_gs(int *empdist, int ntaxa, int nspec);
	double *wht_fit_zm(int *empdist, int ntaxa, int nspec);
	double *wht_fit_lnt(int *empdist, int ntaxa, int nspec);
	
	double *expfinds (double *dist, int S, int mxfds, int ttfds);
	double *expfindspart (double *dist, int S, long *nofds, int unfds, int ttfds);
	double *expfindsln (double mag, int S, int mxfds, int t);
	
	
	double **support_bars_unif (double bestS, int ntaxa, int bestR, int *empdist, double bar);
	double **support_bars_geos (double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
	double **support_bars_zipf (double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
	double **support_bars_lnu (double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);

#else

	extern double calc_aic_c(int parameters, double fit, int measurements); 
	extern double calc_aic(int parameters, double fit); 
	extern double calc_bic(int parameters, double fit, int measurements); 

	extern double calc_likelihood(double *A, int *S, int sizeofS, int Samplesize);
	extern double calc_likelihood_Foote(double *A, int *S, int sizeofS);
	extern double calc_likelihood_exp(double *X, long *nufinds, long *histogram, int unqfnds);

	extern double *fit_gs(int *abundance, int num_taxa, int Samplesize);
	extern double *fit_zm(int *abundance, int num_taxa, int Samplesize);
	extern double *fit_ln(int *abundance, int num_taxa, int Samplesize);
	extern double *fit_fln(int *abundance, int ntaxa, int Samplesize);
	extern double *fit_ls(int *abundance, int num_taxa, int Samplesize);
	extern double *fit_lp(int *abundance, int num_taxa, int Samplesize);

	extern double *chi_fit_gs(int *empdist, int ntaxa, int nspec);
	
	extern double *poi_fit_uni(int *empdist, int ntaxa, int nspec);
	extern double *poi_fit_geo(int *empdist, int ntaxa, int nspec);
	extern double *poi_fit_geos(int *empdist, int ntaxa, int nspec);
	extern double *poi_fit_zipf(int *empdist, int ntaxa, int nspec);
	extern double *poi_fit_lnu(int *empdist, int ntaxa, int nspec);
	extern double *poi_fit_lnt(int *empdist, int ntaxa, int nspec);
	extern double *poi_fit_zsm(int *empdist, int ntaxa, int nspec);
	extern double *poi_fit_zsm_test(int *empdist, int ntaxa, int nspec);
	
	extern double *wht_fit_gs(int *empdist, int ntaxa, int nspec);
	extern double *wht_fit_zm(int *empdist, int ntaxa, int nspec);
	extern double *wht_fit_lnt(int *empdist, int ntaxa, int nspec);
	
	extern double *expfinds (double *dist, int S, int mxfds, int ttfds);
	extern double *expfindspart (double *dist, int S, long *nofds, int unfds, int ttfds);
	extern double *expfindsln (double mag, int S, int mxfds, int t);

	extern double **support_bars_unif (double bestS, int ntaxa, int bestR, int *empdist, double bar);
	extern double **support_bars_geos (double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
	extern double **support_bars_zipf (double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);
	extern double **support_bars_lnu (double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar);

#endif
