/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
/*	Peter J. Wagner: pwagner@fmnh.org
/*	Matthew Kosnik: mkosnik@uchicago.edu
/*
/*	This file is copyright (C) 2001 Peter J. Wagner & Matthew Kosnik
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

#define distribution_fits
#include "distribution_fits.h"
#include "matrixanalysis.h"
#include "matrixreading.h"
#include "minmax.h"
#include "Probability.h"


/* calculates Akaike's Information Criterion (aic)
	L = model's log likelihood
	k = number of parameters.
	2 * L + 2 * k
***********************************************************************/
double calc_aic( int k, double L) {

double aic;

aic = (-2*L) + (2*k);

return aic;
}

/* calculates modified Akaike's Information Criterion (aic-c) Sugiura 1978 
	aic adjusted to accout for sample size (n)
	L = model's log likelihood
 	k = number of parameters.
	n = sample size.
	
	Burnham & Anderson 1998 p. 51 suggest using this when (n / K)<40
***********************************************************************/
double calc_aic_c( int k, double L, int n) {

double aic_c;

aic_c = (-2*L) + (2*k)*((float) n / ((float) n - k - 1));

return aic_c;
}

/* calculates Bayesian Information Criterion (bic) for sent data.
	L = model's log likelihood
	k = number of parameters.
	n = sample size.

	Burnham & Anderson 1998 p. 68 suggest not using this one...
***********************************************************************/
double calc_bic( int k, double L, int n) {

double bic;

bic = -2*L + log(n)*k;

return bic;
}


/* CALCULATES THE LOG LIKELIHOOD OF A GIVEN DISTRIBUTION
Multinomial

NEEDS:
 -	A = theoretical distribution
 -	R = size of A 
 -	S = actual distribution
 -	N = abundance of S
 -	Z = size of S
RETURNS:
 -	L = Summed log likelihood
***********************************************************************/
double calc_likelihood(double *A, int *S, int N, int Z) {

double L = 0.0f;							// Summed Log Likelihood
double a = 0.0f, b = 0.0f;
int i = 0;								// loop variable

// CALCULATE LOG LIKELIHOOD FOR EACH TAXON IN BOTH LISTS
for (i=0; i<Z; i++) {
 
	a = log(A[i]) * S[i];					// a = pi ^ni
	b = log(1-A[i]) * (N-S[i]);				// b = (1 - pi) ^ (N-ni)
	L+=(a+b);							// Sum L for all taxa
	}
if (L>0) L=-1*DBL_MAX;						// if c blows out = very bad fit

return L;
}

/* CALCULATES THE LOG LIKELIHOOD OF A GIVEN DISTRIBUTION
Probability of A finds from N specimens given hypothesized frequency based on Foote (in prep.)
 
 NEEDS:
 -	A[i] = proportion for taxon rank i from the hypothetical distribution
 -	n[i] = number of specimens/occurrences for specimen i
 -	N = total number of specimens
 RETURNS:
 -	L = Summed log likelihood
***********************************************************************/
double calc_likelihood_Foote(double *A, int *n, int N) {

double L = 0.0f;						// Summed Log Likelihood
double lnp = 0.0f;						// taxon support
double a = 0.0f, b = 0.0f, c = 0.0f;	// temp variables
int i;									// loop variable

/* CALCULATE LOG LIKELIHOOD FOR EACH TAXON WITH n > 1	*/
for (i=0; n[i]>1; i++) {
	
	a = log(A[i]) * n[i];				// a = pi ^ni
	b = log(1-A[i]) * (N-n[i]);			// b = (1 - pi) ^ (N-ni)
	c = log(1 - pow((1-A[i]),N));		// c = 1 - ( 1 - pi ) ^ N - the probability of being sampled
	if (c>0)
		c=0;							// sometimes this gets wonky - blows out
	lnp = (a + b) - c;					// L = (a * b) / c
	L+=lnp;								// Sum L for all taxa
	}
if (L>0) L=-1*DBL_MAX;					// if c blows out = very bad fit

return L;
}

/* CALCULATES THE LOG LIKELIHOOD OF A DISTRIBUTION GIVEN OBSERVED AND EXPECTED NUMBERS OF TAXA WITH X FINDS
 P. Wagner 11/14/2003
 NEEDS:
 -	X[i] = expected proportion for taxa with nufinds[i] finds
 -	nufinds[i] = number of finds (of the pool of observed unique numbers of finds)
 -	histogram[n] = observed number of taxa with n finds
 -  unqfinds = the number of number of finds (= richness+1 if all taxa have a unique number, = 2 if all taxa have the same number of finds [+1 for unsampled taxa])
 RETURNS:
 -	L = Summed log likelihood
***********************************************************************/
double calc_likelihood_exp(double *X, long *nufinds, long *histogram, int unqfnds)
{
double L = 0.0f;						// Summed Log Likelihood
double a = 0.0f;	// temp variables
int i, n;								// loop variable

for (i=1; i<unqfnds; ++i)	{
	n=nufinds[i];
	a=histogram[n];
	L = L+(a*log(X[i]));
	}

if (L>0) L=-1*DBL_MAX;					// if c blows out = very bad fit

return L;
}

/* CALCULATE LIKELIHOOD THAT A GIVEN DISTRIBUTION FITS THE GEOMETRIC SERIES. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: slope of geometric series
	- result[2]: optimal richness			
COMMENTS:
	- SUPINC, FITINC defined in header file
***********************************************************************/
double *fit_gs(int *empdist, int ntaxa, int nspec) {
int i = 0;					// LOOP VARIABLE
int r = 0;					// LOOP RICHNESS
double ev = 0.000f;			// LOOP SLOPE
double emin = 1.0f;			// min slope
double ei = 0.000f;			// how much to increment ev in each loop
double es[3];				// previous log likelihoods (cell number = num previous).
double rs[3];				// previous log likelihoods (cell number = num previous).
double bep[3];				// BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format
double *brp;				// BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array
double *fitdist;				// fit distribution
double pbes;				// PREVIOUS BEST SUPPORT (E VALUE)

brp=dvector(3);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;
// increment true richness until that fails to improve likelihood
for (r = ntaxa; ((rs[2] == 0.0f) || (rs[1] > (rs[2] + SUPINC))); r++) {	

	pbes = 0.0f;
	ei = (double) FITINC / 10;
	// increment starting at .1 since we are almost never going to find slopes of 2+
	
	for (i=0; i<3; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = 0.0f;
	// increment slope until that fails to improve likelihood or resolution limit reached
	for (ev = emin; ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))); ev += ei) {	
		
		fitdist = proportional_gs_distribution(ev,r);				// MAKE DISTRIBUTION
		es[0] = calc_likelihood_Foote(fitdist, empdist, nspec);	// CALCULATE SUPPORT
		free_dvector(fitdist);
		
		//Debugging line
		if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {							// IF BETTER THAN BEST FIT
			pbes = bep[0];
			bep[0] = es[0];								// STORE FIT
			bep[1] = ev;									// STORE SLOPE
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
			ev -= ei;									// STEP BACK ONE UNIT
			ei *= -1;									// STEP BACKWARD TO FIND PEAK
			es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP
			} 
		else {									// NOT IMPROVING TRY SMALLER INCREMENT
			ev -= ei;									// STEP BACK ONE UNIT
			ei /= 10;									// SET SMALLER UNIT
		}
		es[2] = es[1];									// Store last 2 attempts to identify
		es[1] = es[0];									// when the peak is past
		}
	if (bep[0] >= brp[0]) {								// IF BETTER THAN BEST FIT
		brp[2] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		}
	rs[2] = rs[1];
	rs[1] = bep[0];
//	printf("\nGS R=%f, ev=%1.10f, S=%4.14f\n",brp[2],brp[1],brp[0]);
	}
return brp;
}

/* CALCULATE LIKELIHOOD THAT A GIVEN DISTRIBUTION FITS THE ZIPF-MANDELBROT. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: slope in log-log space
	- result[2]: optimal richness			
COMMENTS:
	- SUPINC, FITINC defined in header file
***********************************************************************/
double *fit_zm(int *empdist, int ntaxa, int nspec) {
int i = 0;					// LOOP VARIABLE
int r = 0;					// LOOP RICHNESS
double ev = 0.000f;			// LOOP SLOPE
double ein = 1.0f;			// initival slope values
double emin = 1.0f;			// cannot have an evalue less than 1.0
double ei = 0.000f;			// how much to increment ev in each loop
double es[3];				// previous log likelihoods (cell number = num previous).
double rs[3];				// previous log likelihoods (cell number = num previous).
double bep[3];				// BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format
double *brp;				// BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array
double *fitdist;				// fit distribution
double pbes;				// PREVIOUS BEST SUPPORT (E VALUE)

brp=dvector(3);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;

ein=log(ntaxa-1)/(log(empdist[0])-log(empdist[ntaxa-1]));

// increment true richness until that fails to improve likelihood
for (r = ntaxa; ((rs[2] == 0.0f) || (rs[1] > (rs[2] + SUPINC))); r++) {	

	pbes = 0.0f;
	/* mak 2002.04.28 - Sampler indicates that we spend 5 of 7 samples in this function
	So increased increment value to try and improve convergence speed.
	*/
	ei = (double) FITINC;
	for (i=0; i<3; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = 0.0f;
	// increment slope until that fails to improve likelihood or resolution limit reached
	for (ev = ein; ev>emin && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))); ev += ei) {	
		
		fitdist = proportional_zm_distribution(ev,r);				// MAKE DISTRIBUTION
		es[0] = calc_likelihood_Foote(fitdist, empdist, nspec);	// CALCULATE SUPPORT
		free_dvector(fitdist);
		
		//Debugging line
		if (ev<=emin) printf("\nDANGER: ZM R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {							// IF BETTER THAN BEST FIT
			pbes = bep[0];
			bep[0] = es[0];								// STORE FIT
			bep[1] = ev;									// STORE SLOPE
			} 
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
			ev -= ei;									// STEP BACK ONE UNIT
			ei *= -1;									// STEP BACKWARD TO FIND PEAK
			es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP
			} 
		else {									// NOT IMPROVING TRY SMALLER INCREMENT
			ev -= ei;									// STEP BACK ONE UNIT
			ei /= 10;									// SET SMALLER UNIT
			}
		es[2] = es[1];
		es[1] = es[0];
		}
	if (bep[0] > brp[0]) {								// IF BETTER THAN BEST FIT
		brp[2] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		}
	rs[2] = rs[1];
	rs[1] = bep[0];
//	printf("\nZM R=%f, ev=%f, S=%f\n",brp[2],brp[1],brp[0]);
	}
return brp;
}

/* CALCULATE LIKELIHOOD THAT A GIVEN DISTRIBUTION FITS THE LOG SERIES. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: slope of geometric series
	- result[2]: optimal richness			
COMMENTS:
	- SUPINC, FITINC defined in header file
***********************************************************************/
double *fit_ls(int *empdist, int ntaxa, int nspec) {
int i = 0;					// LOOP VARIABLE
int r = 0;					// LOOP RICHNESS
double ev = 0.000f;			// LOOP SLOPE
double emin = ZERO;			// min slope
double ei = 0.000f;			// how much to increment ev in each loop
double es[3];				// previous log likelihoods (cell number = num previous).
double rs[3];				// previous log likelihoods (cell number = num previous).
double bep[3];				// BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format
double *brp;				// BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array
double *fitdist;			// fit distribution
double pbes;				// PREVIOUS BEST SUPPORT (E VALUE)

//printf("\nENTER fit_ls"); fflush(stdout);

brp=dvector(3);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;
// increment true richness until that fails to improve likelihood
for (r = ntaxa; ((rs[2] == 0.0f) || (rs[1] > (rs[2] + SUPINC))); r++) {	

//	printf("\nr\t%d",r); fflush(stdout);

	pbes = 0.0f;
	ei = (double) FITINC / 10;
	// increment starting at .1 since ev need to be between 0 - 1
	
	for (i=0; i<3; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = 0.0f;
	// increment slope until that fails to improve likelihood or resolution limit reached
	for (ev = emin; ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))); ev += ei) {	

//		printf("\t%2.1f",ev); fflush(stdout);

// function removed by P. Wagner 2002.09
		fitdist = proportional_gs_distribution(ev,r);				// MAKE DISTRIBUTION
		es[0] = calc_likelihood_Foote(fitdist, empdist, nspec);	// CALCULATE SUPPORT
		free_dvector(fitdist);
		
//		printf("\nLS R=%d, ev=%1.10f, S=%4.14f\n",r,ev,es[0]);  fflush(stdout);

		//Debugging line
		if (ev<=emin) printf("\nDANGER: LS R=%d, ev=%f, S=%f ",r,ev,es[0]);

//		printf("\nes0 %4.3f\tes1 %4.3f\tes2 %4.3f",es[0],es[1],es[2]); fflush(stdout);
//		printf("\ne\t%4.3f\tes0\t%4.3f",ev,es[0]); fflush(stdout);
		if (es[0] >= bep[0]) {							// IF BETTER THAN BEST FIT
			pbes = bep[0];
			bep[0] = es[0];								// STORE FIT
			bep[1] = ev;									// STORE SLOPE
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
			ev -= ei;									// STEP BACK ONE UNIT
			ei *= -1;									// STEP BACKWARD TO FIND PEAK
			es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP
			} 
		else {									// NOT IMPROVING TRY SMALLER INCREMENT
			ev -= ei;									// STEP BACK ONE UNIT
			ei /= 10;									// SET SMALLER UNIT
			}
		es[2] = es[1];									// Store last 2 attempts to identify
		es[1] = es[0];									// when the peak is past
		}
	if (bep[0] >= brp[0]) {								// IF BETTER THAN BEST FIT
		brp[2] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		}
	rs[2] = rs[1];
	rs[1] = bep[0];
//	printf("\nLS R=%f, ev=%1.10f, S=%4.14f\n",brp[2],brp[1],brp[0]); fflush(stdout);
	}
return brp;
}

/* CALCULATE LIKELIHOOD THAT A GIVEN DISTRIBUTION FITS THE LOG NORMAL ERIES. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- brp[0]: log likelihood
	- brp[1]: magnitude of increase per octave
	- brp[2]: "modal" octave (given in SD's away from Normal Mean)
	- brp[3]: number of octaves
	- brp[4]: optimal richness
***********************************************************************/
double *fit_ln(int *empdist, int ntaxa, int nspec)	{
int i = 0;					// simple loops
int r = 0;					// current richness
int o = OMIN;				// current number of octaves
int ops = 4;				// Octaves per standard deviation
int s = 0;					// loop value of offset
int j, dud;
double ev = 0.00f;			// loop value of mag
double ei = 0.00f;			// how much in increment mag by in loop
double rs[3];				// r previous log likelihoods (cell number = num previous).
double es[3];				// m previous log likelihoods (cell number = num previous).
//double os[3];				// o previous log likelihoods (cell number = num previous).
double ss[3];				// s octave previous log likelihoods (cell number = num previous).
double *brp;				// BEST r PARAMETERS (RICHNESS) - return array format
double bep[5];				// BEST ev PARAMETERS (MAGNITUDE OF CHANGE) - return array format
//double bop[5];			// BEST o PARAMETERS (NUMBER OF OCTAVES) - return array format
double bsp[5];				// BEST S PARAMETERS (STARTING OCATVE) - return array format
double pbes;				// PREVIOUS BEST SUPPORT (E VALUE)
double *fitdist;				// distribution being evaluated
double *ddistribution;			// Normal distribution for Log Normal
int *idistribution;
int *lastidist;				// when richness is low, different modes produce the same distribution
							// this is used to identify and skip these cases

brp = dvector(5);

lastidist=ivector(NSTDV*ops);

for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<5; i++) brp[i] = -1.0*DBL_MAX;
// LOOP - NUMBER OF TAXA
for (r = ntaxa; ((rs[2] == 0.0f) || (rs[1] > (rs[2] + SUPINC))); r++) {	

	for (i=0; i<3; i++) ss[i] = 0.0f;
	for (i=0; i<5; i++) bsp[i] = -1.0*DBL_MAX;
// LOOP - STARTING OCTAVE
	/* changed to ops*NSTDV from o*NSTDV */
/*	for (s = 0; (((ss[2] == 0.0f) || (ss[1] > (ss[2] + SUPINC))) && (s<ops*NSTDV)); s++) {	*/
/* The rounding vagaries mean that mO = 2 is better than mO=3 when mO actually is 4 with low richness */
	dud=0;
	for (s = 0; (dud<2 && (s<ops*NSTDV)); s++) {
		/* wrong numbers were being entered here */
//		ddistribution = NormalDistribution(NSTDV, o, s);		// make basic distribution to speed things up abit
		ddistribution = NormalDistribution(NSTDV, ops, s);		// make basic distribution to speed things up abit
		idistribution = OctaveRichness(r, ddistribution, NSTDV, o);
		free_dvector(ddistribution);
		
		/* Make sure that idistribution is different - it is not always when R is low */
		for (j=0; (j==0 && s<(ops*NSTDV)); s=s)	{
			for (i=0; i<(NSTDV*ops); ++i)	{
				j=abs(idistribution[i]-lastidist[i]);
				if (j>0)	i=(NSTDV*ops);
				}
			if (j==0)	{
				++s;
				free_ivector(idistribution);
				ddistribution = NormalDistribution(NSTDV, ops, s);		// make basic distribution to speed things up abit
				idistribution = OctaveRichness(r, ddistribution, NSTDV, o);
				free_dvector(ddistribution);
				}
			}
		
		ei = (double) FITINC;
		pbes = 0.0f;
		
		/* INSERT ROUTINE TO MAKE SURE THAT THE SAME DISTRIBUTION IS NOT GENERATED */
		/*	THIS WILL HAPPEN WHEN RICHNESS IS LOW */
		
		for (i=0; i<3; i++) es[i] = 0.0f;
		for (i=0; i<5; i++) bep[i] = -1.0*DBL_MAX;
// LOOP - MAGNITUDE OF CHANGE BETWEEN OCTAVES
		for (ev = 1.0f; ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))); ev += ei) {	
			/* wrong numbers were being entered here */
//			fitdist = LNAbundCalc(r, NSTDV, o, ev, idistribution);	// MAKE DISTRIBUTION
			fitdist = LNAbundCalc(r, NSTDV, ops, ev, idistribution);	// MAKE DISTRIBUTION
			es[0] = calc_likelihood_Foote(fitdist, empdist, nspec);	// CALCULATE SUPPORT
			free_dvector(fitdist);
			
			if (es[0] >= bep[0]) {							// IF BETTER THAN BEST FIT
				pbes = bep[0];
				bep[0] = es[0];								// STORE FIT
				bep[1] = ev;									// STORE SLOPE
			} else if ((es[2] > es[1]) && (es[0] > es[1])) {		// TOO FAR: THE PEAK HAS BEEN PASSED
				ev -= ei;									// STEP BACK ONE UNIT
				ei *= -1;									// STEP BACKWARD TO FIND PEAK
				es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP
			} else {									// NOT IMPROVING TRY SMALLER INCREMENT
				ev -= ei;									// STEP BACK ONE UNIT
				ei /= 10;									// SET SMALLER UNIT
			}
			es[2] = es[1];									// Store last 2 attempts to identify
			es[1] = es[0];									// when the peak is past
		}	/* End examination of magnitudes of increase at this modal octave and richness */
		for (i=0; i<NSTDV*ops; ++i)	lastidist[i]=idistribution[i];
		free_ivector(idistribution);
		if (bep[0] >= bsp[0]) {								// IF BETTER THAN BEST FIT
			bsp[2] = s;
			for (i=0; i<2; i++)
				bsp[i] = bep[i];
			dud=0;
			}
		else	++dud;		/* for going up and down hills - allow one dud before trashing things */
		ss[2] = ss[1];
		ss[1] = bep[0];
		
	}	/* End examination of modal octaves at this richness */
	/* was bop instead of bsp! */
	if (bsp[0] >= brp[0]) {								// IF BETTER THAN BEST FIT
		brp[3] = o;
		brp[4] = r;
		for (i=0; i<3; i++)
			brp[i] = bsp[i];
	}
	rs[2] = rs[1];
	rs[1] = bsp[0];
//	printf("\nLNr r=%f, o=%f, s=%f, ev=%f, S=%f\n",brp[4],brp[3],brp[2],brp[1],brp[0]);
}
free_ivector(lastidist);
return brp;
}

/* CALCULATE LIKELIHOOD THAT A GIVEN DISTRIBUTION FITS THE LOG POWER. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- brp[0]: log likelihood
	- brp[1]: ev / c
	- brp[2]: x 
	- brp[3]: optimal richness
	
MAK 2002.07.31 - found to give excellent results (correlations of sim vs. fit)
***********************************************************************/
double *fit_lp(int *empdist, int ntaxa, int nspec) {

int i = 0;					// LOOP VARIABLE
double *fitdist;				// fit distribution

int r = 0;					// LOOP RICHNESS
double rs[3];				// previous log likelihoods (cell number = num previous).
int rc = 0;
double *brp;				// BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array

double x = 0.000f;			// LOOP SLOPE
double xmin = ZERO;			// min slope
double xi = 0.000f;			// how much to increment ev in each loop
int xc = 0;
double xs[3];				// previous log likelihoods (cell number = num previous).
double bxp[4];				// BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format
double pbxs = 0.0f;			// PREVIOUS BEST SUPPORT (E VALUE)

double ev = 0.000f;			// LOOP SLOPE
double emin = ZERO;		// min slope
double ei = 0.000f;			// how much to increment ev in each loop
double es[3];				// previous log likelihoods (cell number = num previous).
double bep[4];				// BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format
double pbes = 0.0f;			// PREVIOUS BEST SUPPORT (E VALUE)

bep[0] = -1.0*DBL_MAX;
bep[1] = emin;
bep[2] = xmin;

brp=dvector(4);
for (i=0; i<2; i++) rs[i] = 0.0f;
for (i=0; i<4; i++) brp[i]= -1.0*DBL_MAX;
// increment true richness until that fails to improve likelihood
for (r = ntaxa; ((rc<1) || (rs[1] > (rs[2] + SUPINC))); r++) {	
	
	rc++;
	xc=0;
	pbxs = 0.0f;
	xi = (double) FITINC;
	for (i=0; i<2; i++) xs[i] = 0.0f;
	for (i=0; i<4; i++) bxp[i] = -1.0*DBL_MAX;
	// increment slope until that fails to improve likelihood or resolution limit reached
	for (x = xmin; ((xc<1) || (bxp[0] > (pbxs + SUPINC))); x += xi) {	

		if (x<=xmin) printf("\nDANGER: LP R=%d, ev=%f, x=%f ",r,ev,x);

		xc++;
		pbes = 0.0f;
		ei = (double) FITINC;
		for (i=0; i<2; i++) es[i] = 0.0f;
		for (i=0; i<4; i++) bep[i] = -1.0*DBL_MAX;
		// increment slope until that fails to improve likelihood or resolution limit reached
		for (ev = emin; ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))); ev += ei) {	
	
			// DO COMPARISON WITH THESE PARAMETERS
			fitdist = proportional_lp_distribution (ev, x, r);
			es[0] = calc_likelihood_Foote(fitdist, empdist, nspec);	// CALCULATE SUPPORT
			free_dvector(fitdist);
						
			//Debugging line
			//printf("\nLP 0 r=%d, ev=%1.10f, x=%1.10f, S=%4.14f",r,ev,x,es[0]); fflush(stdout);
			if (ev<=emin) printf("\nDANGER: LP R=%d, ev=%f, x=%f, S=%f ",r,ev,x,es[0]);

			if (es[0] >= bep[0]) {							// IF BETTER THAN BEST FIT
				pbes = bep[0];
				bep[0] = es[0];								// STORE FIT
				bep[1] = ev;									// STORE SLOPE
				bep[2] = x;
				bep[3] = r;
				}
			else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
				ev -= ei;									// STEP BACK ONE UNIT
				ei *= -1;									// STEP BACKWARD TO FIND PEAK
				es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP
				}
			else {									// NOT IMPROVING TRY SMALLER INCREMENT
				ev -= ei;									// STEP BACK ONE UNIT
				ei /= 10;									// SET SMALLER UNIT
				}
			es[2] = es[1];									// Store last 2 attempts to identify
			es[1] = es[0];									// when the peak is past

			//Debugging line
			//printf("\nLP E r= %d, ev= %1.10f, x= %1.10f, S= %4.14f",(int) bep[3],bep[1],bep[2],bep[0]); fflush(stdout);
			}

		if (bep[0] >= bxp[0]) {							// IF BETTER THAN BEST FIT
			pbxs = bxp[0];
			for (i=0; i<4; i++) bxp[i] = bep[i];
			}
		else if ((xs[2] > xs[1]) && (bep[0] > xs[1]) && ((x-xi)>xmin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
			x -= xi;									// STEP BACK ONE UNIT
			xi *= -1;									// STEP BACKWARD TO FIND PEAK
			xs[1] = xs[2];								// SET PREVIOUS S TO IGNORE OVER STEP
			}
		else {									// NOT IMPROVING TRY SMALLER INCREMENT
			x -= xi;									// STEP BACK ONE UNIT
			xi /= 10;									// SET SMALLER UNIT
			}
		xs[2] = xs[1];									// Store last 2 attempts to identify
		xs[1] = bep[0];									// when the peak is past

		//Debugging line
		//printf("\nLP X r=%d, ev=%1.10f, x=%1.10f, S=%4.14f",(int) bxp[3],bxp[1],bxp[2],bxp[0]); fflush(stdout);
		}
	
	rs[2] = rs[1];
	rs[1] = xs[1];
	// IF BETTER THAN BEST FIT - MAKE THIS THE BEST FIT
	if (bxp[0] >= brp[0])
		for (i=0; i<4; i++) brp[i] = bxp[i];

	//Debugging line
	//printf("\nLP 3 bxp: %1.4f brp: %1.4f",bxp[0],brp[0]); fflush(stdout);
	}

//Debugging line
//printf("\nLP F r=%d, ev=%1.10f, x=%1.10f, S=%4.14f",(int) brp[3],brp[1],brp[2],brp[0]); fflush(stdout);
return brp;
}


/* CALCULATE LIKELIHOOD THAT A GIVEN DISTRIBUTION FITS THE FAUX LOG-NORMAL SERIES. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- brp[0]: log likelihood
	- brp[1]: slope between taxon 1 and 2
	- brp[2]: slope at mode
	- brp[3]: optimal richness			
COMMENTS:
	- SUPINC, FITINC defined in header file
***********************************************************************/
double *fit_fln(int *empdist, int ntaxa, int nspec) {

int i = 0;					// LOOP VARIABLE
int r = 0;					// LOOP RICHNESS
double me = 1.0001f;		// GEOMETRIC DECAY (=SLOPE AT MODE)
double ie = 0.000f;			// INITIAL SLOPE
double iemin = 1.0002f;		// min slope
double iei = 0.000f;		// how much to increment ie in each loop
double ies[3];				// previous log likelihoods (cell number = num previous).
double rs[3];				// previous log likelihoods (cell number = num previous).
double bdep[3];				// BEST ie PARAMETERS (DISTRIBUTION SLOPE) - return array format
double *brp;				// BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array
double *fitdist;			// fit distribution
double pbes;				// PREVIOUS BEST SUPPORT (E VALUE)
int	m=0;

brp=dvector(4);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;
/* increment true richness until that fails to improve likelihood	*/
for (r = ntaxa; ((rs[2] == 0.0f) || (rs[1] > (rs[2] + SUPINC))); r++) {	

	iemin=me+(me-1.0f);
	pbes = 0.0f;
	iei = (double) FITINC / 10;
	/* increment starting at .1 since we are almost never going to find slopes of 2+	*/
	
		for (i=0; i<3; i++) bdep[i] = -1.0*DBL_MAX;
		for (i=0; i<3; i++) ies[i] = 0.0f;
		/* increment slope until that fails to improve likelihood or resolution limit reached	*/
		for (ie = iemin; ((pbes == 0.0f) || (bdep[0] > (pbes + SUPINC))); ie += iei) {	
			
			fitdist = proportional_fln_distribution(me,ie,r,m);				// MAKE DISTRIBUTION
			ies[0] = calc_likelihood_Foote(fitdist, empdist, nspec);	// CALCULATE SUPPORT
			free_dvector(fitdist);
			
			//Debugging line
			if (ie<=iemin) printf("\nDANGER: FLN R=%d, ie=%f, me=%f, S=%f ",r,ie,me,ies[0]);

			if (ies[0] >= bdep[0]) {							// IF BETTER THAN BEST FIT
				pbes = bdep[0];
				bdep[0] = ies[0];							// STORE FIT
				bdep[1] = ie;								// STORE INITIAL SLOPE
				bdep[2] = me;								// STORE SLOPE AT MODE
				}
			else if ((ies[2] > ies[1]) && (ies[0] > ies[1]) && ((ie-iei)>iemin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
				ie -= iei;									// STEP BACK ONE UNIT
				iei *= -1;									// STEP BACKWARD TO FIND PEAK
				ies[1] = ies[2];							// SET PREVIOUS S TO IGNORE OVER STEP
				}
			else {											// NOT IMPROVING TRY SMALLER INCREMENT
				ie -= iei;									// STEP BACK ONE UNIT
				iei /= 10;									// SET SMALLER UNIT
			}
			ies[2] = ies[1];								// Store last 2 attempts to identify
			ies[1] = ies[0];								// when the peak is past
			}
	if (bdep[0] >= brp[0]) {							/* IF BETTER THAN BEST FIT */
		brp[3] = r;										/* New best richness */
		for (i=0; i<3; i++)						/* new best fit, best init slope, best modal slope */
			brp[i] = bdep[i];
		}
	rs[2] = rs[1];
	rs[1] = bdep[0];
//	printf("\nGS R=%f, ie=%1.10f, S=%4.14f\n",brp[2],brp[1],brp[0]);
	}
/**************************
brp[0]: Log likelihood
brp[1]: Initial Slope
brp[2]: Modal Slope
brp[3]: Best richness
**************************/
return brp;
}

/* CALCULATE LIKELIHOOD THAT A GIVEN DISTRIBUTION FITS THE GEOMETRIC SERIES ASSUMING A POISSON DISTRIBUTION FOR EXPECTED SPECIES WITH X FINDS. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: slope of geometric series
	- result[2]: optimal richness			
COMMENTS:
	- SUPINC, FITINC defined in header file
***********************************************************************/
double *poi_fit_gs(int *empdist, int ntaxa, int nspec) 
{
int i = 0;				/* LOOP VARIABLE	*/
int r = 0;					/* LOOP RICHNESS	*/
int ri;						/* richness increment	*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
double ev = 0.000f;			/* LOOP SLOPE 		*/
double emin = 1.0f;			/* min slope		*/
double ein = 0.000f;		/* initial slope	*/
double ei = 0.000f;			/* how much to increment ev in each loop							*/
double es[3];				/* previous log likelihoods (cell number = num previous).		*/
double rs[3];				/* previous log likelihoods (cell number = num previous).		*/
double bep[3];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/
double *brp;				/* BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array	*/
double *fitdist;			/* fit distribution												*/
double *expect;				/* expected number of species with 0Émax finds					*/
double *obsrvd;				/* observed number of species with 0Émax finds					*/
double pbes;				/* previous best support for decay rate							*/

brp=dvector(3);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
ein=pow(empdist[0],(((double) empdist[ntaxa-1])/((double) ntaxa)));

ri=ntaxa/10;
if (ri<2)	ri=2;
// increment true richness until that fails to improve likelihood
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (r=ntaxa; abs(ri)>0 && r>=ntaxa; r+=ri)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;
	ei = (double) FITINC / 10;
	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	
	for (i=0; i<3; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = 0.0f;
	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; (ev>=emin && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC)))); ev += ei) {	
		
		/* generate geometric distribution with richness r and decay of ev */
		fitdist = proportional_gs_distribution(ev,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,2*empdist[0],nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector(expect,obsrvd,2*empdist[0]);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {					/* IF BETTER THAN BEST FIT */
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								// STORE SLOPE
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
			ev -= ei;									// STEP BACK ONE UNIT
			ei *= -1;									// STEP BACKWARD TO FIND PEAK
			es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP
			}
		else {											// NOT IMPROVING TRY SMALLER INCREMENT
			ev -= ei;									// STEP BACK ONE UNIT
			ei /= 10;									// SET SMALLER UNIT
			}
		
		while (ev+ei<= emin)	ei/=2;
		
		es[2] = es[1];									// Store last 2 attempts to identify
		es[1] = es[0];									// when the peak is past
		}
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								// IF BETTER THAN BEST FIT
		brp[2] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		ein=bep[1];										/* set initial decay to best decay found so far */ 
		}
	/* optimal richness is overshot */
	else if ((rs[2]>rs[1] && (rs[0]>rs[1])) && (abs(ri)>1))	{
		r-=ri;				/* step back one unit	 */
		ri*=-1;				/* step backwards to peak */
		rs[1]=rs[2];		/* set to prior ln L to ignore over step */
		}
	else	{
		r-=ri;
		if (abs(ri)>1)	ri/=2;
		else			ri=0;
		}
	
	/* geometric will go on and on forever with miniscule increases when it is a very poor fit */
	if ((rs[0]>=rs[1] && rs[0]<(rs[1]+SUPINC)) && (rs[1]>=rs[2] && rs[1]<(rs[2]+SUPINC)))	ri=0;
	
	rs[2] = rs[1];
	rs[1] = bep[0];
	}
return brp;
}

/* CALCULATE LIKELIHOOD THAT A GIVEN DISTRIBUTION FITS THE GEOMETRIC SERIES ASSUMING A POISSON DISTRIBUTION FOR EXPECTED SPECIES WITH X FINDS. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: slope of geometric series
	- result[2]: optimal richness			
COMMENTS:
	- SUPINC, FITINC defined in header file
***********************************************************************/
double *poi_fit_zm(int *empdist, int ntaxa, int nspec)
{
int i = 0;				/* LOOP VARIABLE	*/
int r = 0;					/* LOOP RICHNESS	*/
int ri;						/* richness increment	*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
/*int mxspc;					/* maximum number of possible finds to consider when calculating P */
double ev = 0.000f;			/* LOOP SLOPE 		*/
double emin = 1.0f;			/* min slope		*/
double ein = 1.000f;		/* initial slope	*/
double ei = 0.000f;			/* how much to increment ev in each loop							*/
double es[3];				/* previous log likelihoods (cell number = num previous).		*/
double rs[3];				/* previous log likelihoods (cell number = num previous).		*/
double bep[3];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/
double *brp;				/* BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array	*/
double *fitdist;			/* fit distribution												*/
double *expect;				/* expected number of species with 0Émax finds					*/
double *obsrvd;				/* observed number of species with 0Émax finds					*/
double pbes;				/* previous best support for decay rate							*/

brp=dvector(3);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
ein=log(ntaxa-1)/(log(empdist[0])-log(empdist[ntaxa-1]));

ri=ntaxa/10;
if (ri<2)	ri=2;
// increment true richness until that fails to improve likelihood
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (r=ntaxa; abs(ri)>0 && r>=ntaxa; r+=ri)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;
	ei = (double) FITINC / 10;
	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	
	for (i=0; i<3; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = 0.0f;
	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; (ev>=emin && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC)))); ev += ei) {	
		
		/*Debugging line */
		if (ev<=emin)	printf("\nDANGER: ZM R=%d, ev=%f, S=%f ",r,ev,es[0]);

		/* generate Zipf-Mandelbrot distribution with richnes r and decay of ev */
		fitdist = proportional_zm_distribution(ev,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,2*empdist[0],nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector(expect,obsrvd,2*empdist[0]);
		free_dvector(expect);

		if (es[0] >= bep[0]) {					// IF BETTER THAN BEST FIT
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								// STORE SLOPE
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
			ev -= ei;									// STEP BACK ONE UNIT
			ei *= -1;									// STEP BACKWARD TO FIND PEAK
			es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP
			}
		else {											// NOT IMPROVING TRY SMALLER INCREMENT
			ev -= ei;									// STEP BACK ONE UNIT
			ei /= 10;									// SET SMALLER UNIT
			}
		es[2] = es[1];									// Store last 2 attempts to identify
		es[1] = es[0];									// when the peak is past
		}
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								// IF BETTER THAN BEST FIT
		brp[2] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		ein=bep[1];										/* set initial decay to best decay found so far */ 
		}
	/* optimal richness is overshot */
	else if ((rs[2]>rs[1] && (rs[0]>rs[1])) && (abs(ri)>1))	{
		r-=ri;				/* step back one unit	 */
		ri*=-1;				/* step backwards to peak */
		rs[1]=rs[2];		/* set to prior ln L to ignore over step */
		}
	else	{
		r-=ri;
		if (abs(ri)>1)	ri/=2;
		else			ri=0;
		}
	rs[2] = rs[1];
	rs[1] = bep[0];
//	printf("\nGS R=%f, ev=%1.10f, S=%4.14f\n",brp[2],brp[1],brp[0]);
	}
return brp;
}

/* FIND THE MOST-LIKELY FAUX-LOG-NORMAL SERIES, VARYING INTIAL SLOPE BUT ASSUMING MODAL SLOPE TO BE 0. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: slope of geometric series
	- result[2]: optimal richness			
COMMENTS:
	- SUPINC, FITINC defined in header file
***********************************************************************/
double *poi_fit_fln_1(int *empdist, int ntaxa, int nspec) 
{
int i = 0;					/* LOOP VARIABLE	*/
int r = 0;					/* LOOP RICHNESS	*/
int ri;						/* richness increment	*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
double ev = 0.000f;			/* LOOP SLOPE 		*/
double emin = 1.0f;			/* min slope		*/
double ein = 0.000f;		/* initial slope	*/
double ei = 0.000f;			/* how much to increment ev in each loop							*/
double es[3];				/* previous log likelihoods (cell number = num previous).			*/
double rs[3];				/* previous log likelihoods (cell number = num previous).			*/
double bep[3];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/
double *brp;				/* BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array		*/
double *fitdist;			/* fit distribution													*/
double *expect;				/* expected number of species with 0Émax finds						*/
double *obsrvd;				/* observed number of species with 0Émax finds						*/
double pbes;				/* previous best support for decay rate								*/
double mode;				/* modal of log-normal, given as a taxon number						*/

brp=dvector(3);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the initial slope, based on the slope between the first and second species */
ein=((double) empdist[0])/((double) empdist[1]);

ri=ntaxa/10;
if (ri<2)	ri=2;
// increment true richness until that fails to improve likelihood
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (r=ntaxa; abs(ri)>0 && r>=ntaxa; r+=ri)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;
	ei = (double) FITINC / 10;
	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	
	for (i=0; i<3; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = 0.0f;
	
	if (r%2==1)	mode = r/2;
	else		mode = (r/2)-0.5;

	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; (ev>=emin && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC)))); ev += ei) {	
		
		/* generate geometric distribution with richnes r and decay of ev */
		fitdist = proportional_fln_distribution(1.0f, ev,r,mode);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,2*empdist[0],nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector(expect,obsrvd,2*empdist[0]);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: fLN R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {							// IF BETTER THAN BEST FIT
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								// STORE SLOPE
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
			ev -= ei;									// STEP BACK ONE UNIT
			ei *= -1;									// STEP BACKWARD TO FIND PEAK
			es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP
			}
		else {											// NOT IMPROVING TRY SMALLER INCREMENT
			ev -= ei;									// STEP BACK ONE UNIT
			ei /= 10;									// SET SMALLER UNIT
			}
		es[2] = es[1];									// Store last 2 attempts to identify
		es[1] = es[0];									// when the peak is past
		}
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								// IF BETTER THAN BEST FIT
		brp[2] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		ein=bep[1];										/* set initial decay to best decay found so far */ 
		}
	/* optimal richness is overshot */
	else if ((rs[2]>rs[1] && (rs[0]>rs[1])) && (abs(ri)>1))	{
		r-=ri;				/* step back one unit	 */
		ri*=-1;				/* step backwards to peak */
		rs[1]=rs[2];		/* set to prior ln L to ignore over step */
		}
	else	{
		r-=ri;
		if (abs(ri)>1)	ri/=2;
		else			ri=0;
		}
	rs[2] = rs[1];
	rs[1] = bep[0];
//	printf("\nGS R=%f, ev=%1.10f, S=%4.14f\n",brp[2],brp[1],brp[0]);
	}
return brp;
}

/* FIND THE MOST-LIKELY FAUX-LOG-NORMAL SERIES, VARYING INTIAL AND MODAL SLOPE, BUT HOLDING THE MODE IN THE MIDDLE. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: slope of geometric series
	- result[2]: optimal richness			
COMMENTS:
	- SUPINC, FITINC defined in header file
***********************************************************************/
double *poi_fit_fln_2(int *empdist, int ntaxa, int nspec) 
{
int i = 0;				/* LOOP VARIABLE	*/
int r = 0;					/* LOOP RICHNESS	*/
int ri;						/* richness increment	*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
double dev = 0.000f;		/* LOOP SHIFT IN SLOPE 		*/
double dmin = 1.0f;			/* min change in ev slope		*/
double din = 0.000f;		/* initial change ev slope	*/
double di = 0.000f;			/* how much to increment change in ev in each loop					*/
double ev = 0.000f;			/* LOOP SLOPE 														*/
double emin = 1.0f;			/* min slope														*/
double ein = 0.000f;		/* initial slope													*/
double ei = 0.000f;			/* how much to increment ev in each loop							*/
double ds[3];				/* previous log likelihoods (cell number = num previous).			*/
double es[3];				/* previous log likelihoods (cell number = num previous).			*/
double rs[3];				/* previous log likelihoods (cell number = num previous).			*/
double bdp[3];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format 	*/
double bep[4];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format 	*/
double *brp;				/* BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array		*/
double *fitdist;			/* fit distribution													*/
double *expect;				/* expected number of species with 0Émax finds						*/
double *obsrvd;				/* observed number of species with 0Émax finds						*/
double pbds;				/* previous best support for decay shift rate						*/
double pbes;				/* previous best support for decay rate								*/
double mode;				/* modal of log-normal, given as a taxon number						*/

brp=dvector(4);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<4; i++) brp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
//ein=pow(empdist[0],(((double) empdist[ntaxa-1])/((double) ntaxa)));
ein = emin;		/* begin searching for just the shift parameter */
din=((double) empdist[0])/((double) empdist[1]);
if (din<ein)	din=ein+0.01;

ri=ntaxa/10;
if (ri<2)	ri=2;
// increment true richness until that fails to improve likelihood
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (r=ntaxa; abs(ri)>0 && r>=ntaxa; r+=ri)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;
	ei = (double) FITINC / 10;
	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	
	for (i=0; i<4; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i]=0.0f;

	if (r%2==1)	mode = r/2;
	else		mode = (r/2)-0.5;

	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev=ein; ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))); ev += ei) {	
		
		for (i=0; i<3; i++) bdp[i]=-1.0*DBL_MAX;
		for (i=0; i<3; i++) ds[i]=0.0f;

		di = (double) FITINC / 10;
		/*Debugging line */
/*		if (ev<=emin) printf("\nDANGER: fln R=%d, ev=%f, S=%f ",r,ev,es[0]);	*/

		pbds = 0.0f;
		for (dev=din; ((pbds==0.0f) || (bdp[0]>(pbds+SUPINC))); dev+=di)	{
			/* generate geometric distribution with richnes r and decay of ev */
			fitdist = proportional_fln_distribution(ev,dev,r,mode);				/* MAKE DISTRIBUTION */
			/* find the expected proportions of taxa with 0Éx finds */
			expect=expfinds(fitdist,r,2*empdist[0],nspec);
			free_dvector(fitdist);
			ds[0] = lnPoisson_vector(expect,obsrvd,2*empdist[0]);
			free_dvector(expect);

			if (ds[0] >= bdp[0]) {						/* IF BETTER THAN BEST FIT */
				pbds = bdp[0];									/* save last best ssq for evenness */
				es[0]=bdp[0]=ds[0];								/* STORE FIT */
				bdp[1] = dev;									/* STORE shift in slope	*/
				}
			else if ((ds[2] > ds[1]) && (ds[0] > ds[1]) && ((dev-di)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
				dev -= di;									// STEP BACK ONE UNIT
				di *= -1;									// STEP BACKWARD TO FIND PEAK
				ds[1] = ds[2];								// SET PREVIOUS S TO IGNORE OVER STEP
				}
			else {											// NOT IMPROVING TRY SMALLER INCREMENT
				dev -= di;									// STEP BACK ONE UNIT
				di /= 10;									// SET SMALLER UNIT
				}
			/* it might go on and on forever with miniscule increases when it is a very poor fit */
			if ((ds[0]>=ds[1] && ds[0]<(ds[1]+SUPINC)) && (ds[1]>=ds[2] && ds[1]<(ds[2]+SUPINC)))	break;

			ds[2] = ds[1];									// Store last 2 attempts to identify
			ds[1] = ds[0];									// when the peak is past

			}

		if (es[0] >= bep[0]) {							// IF BETTER THAN BEST FIT
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0] = bep[0] = es[0];						/* STORE FIT */
			din = bep[1] = bdp[1];
			bep[2] = ev;								// STORE SLOPE
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
			ev -= ei;									// STEP BACK ONE UNIT
			ei *= -1;									// STEP BACKWARD TO FIND PEAK
			es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP
			}
		else {											// NOT IMPROVING TRY SMALLER INCREMENT
			ev -= ei;									// STEP BACK ONE UNIT
			ei /= 10;									// SET SMALLER UNIT
			}

		/* it might go on and on forever with miniscule increases when it is a very poor fit */
		if ((es[0]>=es[1] && es[0]<(es[1]+SUPINC)) && (es[1]>=es[2] && es[1]<(es[2]+SUPINC)))	break;
		es[2] = es[1];									// Store last 2 attempts to identify
		es[1] = es[0];									// when the peak is past

		}
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								// IF BETTER THAN BEST FIT
		brp[3] = r;
		for (i=0; i<3; i++)
			brp[i] = bep[i];
		din=bep[1];										/* set initial decay to best decay found so far */ 
		ein=bep[2];										/* set initial decay to best decay found so far */ 
		}
	/* optimal richness is overshot */
	else if ((rs[2]>rs[1] && (rs[0]>rs[1])) && (abs(ri)>1))	{
		r-=ri;				/* step back one unit	 */
		ri*=-1;				/* step backwards to peak */
		rs[1]=rs[2];		/* set to prior ln L to ignore over step */
		}
	else	{
		r-=ri;
		if (abs(ri)>1)	ri/=2;
		else			ri=0;
		}

	/* it might go on and on forever with miniscule increases when it is a very poor fit */
	if ((rs[0]>=rs[1] && rs[0]<(rs[1]+SUPINC)) && (rs[1]>=rs[2] && rs[1]<(rs[2]+SUPINC)))	ri=0;

	rs[2] = rs[1];
	rs[1] = bep[0];
//	printf("\nGS R=%f, ev=%1.10f, S=%4.14f\n",brp[2],brp[1],brp[0]);
	}
return brp;
}

/* FIND THE MOST-LIKELY FAUX-LOG-NORMAL SERIES, VARYING INTIAL AND MODAL SLOPE AND THE MODE. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: slope of geometric series
	- result[2]: optimal richness			
COMMENTS:
	- SUPINC, FITINC defined in header file
***********************************************************************/
double *poi_fit_fln_3(int *empdist, int ntaxa, int nspec) 
{
int i = 0;					/* LOOP VARIABLE	*/
int r = 0;					/* LOOP RICHNESS	*/
int ri;						/* richness increment													*/
int rin;					/* initial richness to use in each search (begins as ntaxa)				*/
double mode;				/* mode of log-normal, given as a taxon number							*/
int mi=((double) ntaxa)/4;	/* how much to increment the mode										*/
double imode;				/* initial mode (set to the median of taxa ranks)						*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
double dev = 0.000f;		/* LOOP SHIFT IN SLOPE 													*/
double dmin = 1.0f;			/* min change in ev slope												*/
double din = 0.000f;		/* initial change ev slope												*/
double di = 0.000f;			/* how much to increment change in ev in each loop						*/
double ev = 0.000f;			/* LOOP SLOPE 															*/
double emin = 1.00f;		/* min slope															*/
double ein = 0.000f;		/* initial slope														*/
double ei = 0.000f;			/* how much to increment ev in each loop								*/
double ds[3];				/* previous initial decay log likelihoods (cell number = num previous).	*/
double es[3];				/* previous modal decay log likelihoods (cell number = num previous).	*/
double rs[3];				/* previous richness log likelihoods (cell number = num previous).		*/
double ms[3];				/* previous mode log likelihoods (cell number = num previous).			*/
double bdp[3];				/* BEST initial decay parameters - return array format 					*/
double bep[4];				/* BEST modal decay parameters - return array format 					*/
double brp[5];				/* BEST richness parameters - return array								*/
double *bmp;				/* Best mode parameters - return array format							*/
double *fitdist;			/* fit distribution														*/
double *expect;				/* expected number of species with 0Émax finds							*/
double *obsrvd;				/* observed number of species with 0Émax finds							*/
double pbds;				/* previous best support for initial decay								*/
double pbes;				/* previous best support for modal decay								*/
double pbrs;				/* previous best support for richness									*/
double pbms;				/* previous best support for mode location								*/

bmp=dvector(5);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<5; i++) bmp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
//ein=pow(empdist[0],(((double) empdist[ntaxa-1])/((double) ntaxa)));
ein = emin;		/* begin searching for just the shift parameter */
din=((double) empdist[0])/((double) empdist[1]);
if (din<ein)	din=ein+0.01;

ri=ntaxa/10;
if (ri<2)	ri=2;
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
if (r%2==1)	imode = ntaxa/2;
else		imode = (ntaxa/2)-0.5;

pbms=0.0f;
rin=ntaxa;
/* adjust mode until that fails to improve likelihood	*/
for (mode=imode; abs(mi)>0 && ((pbms == 0.0f) || (bmp[0] > (pbms + SUPINC))); mode+=mi)	{
	pbrs=0.0f;
	for (i=0; i<4; i++) brp[i]= -1.0*DBL_MAX;
	/* adjust true richness until that fails to improve likelihood	*/
	ri=ntaxa/10;
	if (ri<2)	ri=2;
	for (r=ntaxa; abs(ri)>0 && r>=ntaxa; r+=ri)	{
		obsrvd[0]=r-ntaxa;
		pbes = 0.0f;
		ei = (double) FITINC / 10;
		/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
		
		for (i=0; i<4; i++) bep[i] = -1.0*DBL_MAX;
		for (i=0; i<3; i++) es[i]=0.0f;

		/* increment slope until that fails to improve likelihood or resolution limit reached */
		for (ev=ein; ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))); ev += ei) {	
			
			for (i=0; i<3; i++) bdp[i]=-1.0*DBL_MAX;
			for (i=0; i<3; i++) ds[i]=0.0f;

			di = (double) FITINC / 10;
			/*Debugging line */
	/*		if (ev<=emin) printf("\nDANGER: fln R=%d, ev=%f, S=%f ",r,ev,es[0]);	*/

			pbds = 0.0f;
			for (dev=din; ((pbds==0.0f) || (bdp[0]>(pbds+SUPINC))); dev+=di)	{
				/* generate geometric distribution with richnes r and decay of ev */
				fitdist = proportional_fln_distribution(ev,dev,r,mode);				/* MAKE DISTRIBUTION */
				/* find the expected proportions of taxa with 0Éx finds */
				expect=expfinds(fitdist,r,2*empdist[0],nspec);
				free_dvector(fitdist);
				ds[0] = lnPoisson_vector(expect,obsrvd,2*empdist[0]);
				free_dvector(expect);

				if (ds[0] >= bdp[0]) {						/* IF BETTER THAN BEST FIT */
					pbds = bdp[0];									/* save last best ssq for evenness */
					es[0]=bdp[0]=ds[0];								/* STORE FIT */
					bdp[1] = dev;									/* STORE shift in slope	*/
					}
				else if ((ds[2] > ds[1]) && (ds[0] > ds[1]) && ((dev-di)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
					dev -= di;									// STEP BACK ONE UNIT
					di *= -1;									// STEP BACKWARD TO FIND PEAK
					ds[1] = ds[2];								// SET PREVIOUS S TO IGNORE OVER STEP
					}
				else {											// NOT IMPROVING TRY SMALLER INCREMENT
					dev -= di;									// STEP BACK ONE UNIT
					di /= 10;									// SET SMALLER UNIT
					}

				/* it might go on and on forever with miniscule increases when it is a very poor fit */
				if ((ds[0]>=ds[1] && ds[0]<(ds[1]+SUPINC)) && (ds[1]>=ds[2] && ds[1]<(ds[2]+SUPINC)))	break;
				ds[2] = ds[1];									// Store last 2 attempts to identify
				ds[1] = ds[0];									// when the peak is past
				}

			if (es[0] >= bep[0]) {							// IF BETTER THAN BEST FIT
				pbes = bep[0];								/* save last best ssq for evenness */
				rs[0] = bep[0] = es[0];						/* STORE FIT */
				din = bep[1] = bdp[1];
				bep[2] = ev;								// STORE SLOPE
				}
			else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
				ev -= ei;									// STEP BACK ONE UNIT
				ei *= -1;									// STEP BACKWARD TO FIND PEAK
				es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP
				}
			else {											// NOT IMPROVING TRY SMALLER INCREMENT
				ev -= ei;									// STEP BACK ONE UNIT
				ei /= 10;									// SET SMALLER UNIT
				}
			/* it might go on and on forever with miniscule increases when it is a very poor fit */
			if ((es[0]>=es[1] && es[0]<(es[1]+SUPINC)) && (es[1]>=es[2] && es[1]<(es[2]+SUPINC)))	break;
			es[2] = es[1];									// Store last 2 attempts to identify
			es[1] = es[0];									// when the peak is past
			}
		/* if this richness is better than the last */
		if (bep[0] >= brp[0]) {								// IF BETTER THAN BEST FIT
			brp[3] = r;
			for (i=0; i<3; i++)
				brp[i] = bep[i];
			ms[0]=brp[0];
			din=bep[1];										/* set initial decay to best decay found so far */ 
			ein=bep[2];										/* set initial decay to best decay found so far */ 
			}
		/* optimal richness is overshot */
		else if ((rs[2]>rs[1] && (rs[0]>rs[1])) && (abs(ri)>1))	{
			r-=ri;				/* step back one unit	 */
			ri*=-1;				/* step backwards to peak */
			rs[1]=rs[2];		/* set to prior ln L to ignore over step */
			}
		else	{
			r-=ri;
			if (abs(ri)>1)	ri/=2;
			else			ri=0;
			}
		rs[2] = rs[1];
		rs[1] = bep[0];
	//	printf("\nGS R=%f, ev=%1.10f, S=%4.14f\n",brp[2],brp[1],brp[0]);
		}

	/* if this mode is better than the last */
	if (brp[0] >= bmp[0]) {								// IF BETTER THAN BEST FIT
		bmp[4] = mode;
		for (i=0; i<4; i++)
			bmp[i] = brp[i];
		din=bmp[1];										/* set initial initial slope to best initial slope found so far */ 
		ein=bmp[2];										/* set initial modal decay to best modal decay found so far 	*/ 
		rin=bmp[3];										/* set initial richness to best richness found so far			*/
		}
	/* optimal richness is overshot */
	else if ((ms[2]>ms[1] && (ms[0]>ms[1])) && (abs(mi)>1))	{
		mode-=mi;				/* step back one unit	 */
		mi*=-1;					/* step backwards to peak */
		ms[1]=ms[2];		/* set to prior ln L to ignore over step */
		}
	else	{
		mode-=mi;
		if (abs(mi)>1)	mi/=2;
		else			mi=0;
		}
	ms[2] = ms[1];
	ms[1] = brp[0];
	}

return bmp;
}

/* CALCULATE LIKELIHOOD THAT A GIVEN DISTRIBUTION FITS THE FAUX LOG-NORMAL SERIES. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- brp[0]: log likelihood
	- brp[1]: slope between taxon 1 and 2
	- brp[2]: slope at mode
	- brp[3]: optimal richness			
COMMENTS:
	- SUPINC, FITINC defined in header file
***********************************************************************/
double *poi_fit_fln(int *empdist, int ntaxa, int nspec)
{
int i = 0;					// LOOP VARIABLE
int r = 0;					// LOOP RICHNESS
double me = 1.0001f;		// GEOMETRIC DECAY (=SLOPE AT MODE)
double ie = 0.000f;			// INITIAL SLOPE
double iemin = 1.0002f;		// min slope
double iei = 0.000f;		// how much to increment ie in each loop
double ies[3];				// previous log likelihoods (cell number = num previous).
double rs[3];				// previous log likelihoods (cell number = num previous).
double bdep[3];				// BEST ie PARAMETERS (DISTRIBUTION SLOPE) - return array format
double *brp;				// BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array
double *fitdist;			// fit distribution
double pbes;				// PREVIOUS BEST SUPPORT (E VALUE)
double *expect;				/* expected number of species with 0Émax finds					*/
double *obsrvd;				/* observed number of species with 0Émax finds					*/
int	m=0;

brp=dvector(4);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;
/* increment true richness until that fails to improve likelihood	*/
for (r = ntaxa; ((rs[2] == 0.0f) || (rs[1] > (rs[2] + SUPINC))); r++) {	

	iemin=me+(me-1.0f);
	pbes = 0.0f;
	iei = (double) FITINC / 10;
	/* increment starting at .1 since we are almost never going to find slopes of 2+	*/
	
		for (i=0; i<3; i++) bdep[i] = -1.0*DBL_MAX;
		for (i=0; i<3; i++) ies[i] = 0.0f;
		/* increment slope until that fails to improve likelihood or resolution limit reached	*/
		for (ie = iemin; ((pbes == 0.0f) || (bdep[0] > (pbes + SUPINC))); ie += iei) {	
			
			fitdist = proportional_fln_distribution(me,ie,r,m);				// MAKE DISTRIBUTION
			/* find the expected proportions of taxa with 0Éx finds */
			expect=expfinds(fitdist,r,2*empdist[0],nspec);
			free_dvector(fitdist);
			ies[0] = lnPoisson_vector(expect,obsrvd,2*empdist[0]);
			free_dvector(expect);
			
			//Debugging line
			if (ie<=(iemin-SUPINC)) printf("\nDANGER: FLN R=%d, ie=%f, me=%f, S=%f ",r,ie,me,ies[0]);

			if (ies[0] >= bdep[0]) {							// IF BETTER THAN BEST FIT
				pbes = bdep[0];
				bdep[0] = ies[0];							// STORE FIT
				bdep[1] = ie;								// STORE INITIAL SLOPE
				bdep[2] = me;								// STORE SLOPE AT MODE
				}
			else if ((ies[2] > ies[1]) && (ies[0] > ies[1]) && ((ie-iei)>iemin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
				ie -= iei;									// STEP BACK ONE UNIT
				iei *= -1;									// STEP BACKWARD TO FIND PEAK
				ies[1] = ies[2];							// SET PREVIOUS S TO IGNORE OVER STEP
				}
			else {											// NOT IMPROVING TRY SMALLER INCREMENT
				ie -= iei;									// STEP BACK ONE UNIT
				iei /= 10;									// SET SMALLER UNIT
				}
			ies[2] = ies[1];								// Store last 2 attempts to identify
			ies[1] = ies[0];								// when the peak is past
			}
	if (bdep[0] >= brp[0]) {							/* IF BETTER THAN BEST FIT */
		brp[3] = r;										/* New best richness */
		for (i=0; i<3; i++)						/* new best fit, best init slope, best modal slope */
			brp[i] = bdep[i];
		}
	rs[2] = rs[1];
	rs[1] = bdep[0];
//	printf("\nGS R=%f, ie=%1.10f, S=%4.14f\n",brp[2],brp[1],brp[0]);
	}
/**************************
brp[0]: Log likelihood
brp[1]: Initial Slope
brp[2]: Modal Slope
brp[3]: Best richness
**************************/
return brp;
}


/* CALCULATE LIKELIHOOD THAT A GIVEN DISTRIBUTION FITS THE GEOMETRIC SERIES. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: slope of geometric series
	- result[2]: optimal richness			
COMMENTS:
	- SUPINC, FITINC defined in header file
***********************************************************************/
double *chi_fit_gs(int *empdist, int ntaxa, int nspec) 
{
int i = 0;				/* LOOP VARIABLE	*/
int r = 0;					/* LOOP RICHNESS	*/
int ri;						/* richness increment	*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
/*int mxspc;					/* maximum number of possible finds to consider when calculating P */
double ev = 0.000f;			/* LOOP SLOPE 		*/
double emin = 1.0f;			/* min slope		*/
double ein = 0.000f;		/* initial slope	*/
double ei = 0.000f;			/* how much to increment ev in each loop							*/
double es[3];				/* previous log likelihoods (cell number = num previous).		*/
double rs[3];				/* previous log likelihoods (cell number = num previous).		*/
double bep[3];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/
double *brp;				/* BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array	*/
double *fitdist;			/* fit distribution												*/
double *expect;				/* expected number of species with 0Émax finds					*/
//double *exp;				/* expected number of species with 0Émax finds					*/
double *obsrvd;				/* observed number of species with 0Émax finds					*/
//long *obsfnd;				/* array of unique values for numbers of finds					*/
double pbes;				/* previous best support for decay rate							*/

brp=dvector(3);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
ein=pow(empdist[0],(((double) empdist[ntaxa-1])/((double) ntaxa)));

ri=ntaxa/10;
if (ri<2)	ri=2;
// increment true richness until that fails to improve likelihood
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (r=ntaxa; abs(ri)>0; r+=ri)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;
	ei = (double) FITINC / 10;
	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	
	for (i=0; i<3; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = 0.0f;
	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))); ev += ei) {	
		
		/* generate geometric distribution with richnes r and decay of ev */
		fitdist = proportional_gs_distribution(ev,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,empdist[0],nspec);
		free_dvector(fitdist);
		es[0] = -1*dsumsqdiffs(expect,obsrvd,1+empdist[0]);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=(emin-SUPINC)) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {							// IF BETTER THAN BEST FIT
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								// STORE SLOPE
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED
			ev -= ei;									// STEP BACK ONE UNIT
			ei *= -1;									// STEP BACKWARD TO FIND PEAK
			es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP
			}
		else {											// NOT IMPROVING TRY SMALLER INCREMENT
			ev -= ei;									// STEP BACK ONE UNIT
			ei /= 10;									// SET SMALLER UNIT
			}
		es[2] = es[1];									// Store last 2 attempts to identify
		es[1] = es[0];									// when the peak is past
		}
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								// IF BETTER THAN BEST FIT
		brp[2] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		ein=bep[1];										/* set initial decay to best decay found so far */ 
		}
	/* optimal richness is overshot */
	else if ((rs[2]>rs[1] && (rs[0]>rs[1])) && (abs(ri)>1))	{
		r-=ri;				/* step back one unit	 */
		ri*=-1;				/* step backwards to peak */
		rs[1]=rs[2];		/* set to prior ln L to ignore over step */
		}
	else	{
		r-=ri;
		if (abs(ri)>1)	ri/=2;
		else			ri=0;
		}
	rs[2] = rs[1];
	rs[1] = bep[0];
//	printf("\nGS R=%f, ev=%1.10f, S=%4.14f\n",brp[2],brp[1],brp[0]);
	}
return brp;
}


/* expectedfinds - finds the expected number of taxa with 0...max finds given binomial probabilities, total finds and distribution
/* Requires:
		dist - an array giving relative abundances of S taxa (must sum to 1.0);
		S - the length of dist (i.e., richness);
		mxfds - the number of unique sample numbers - this equals S+1 if all species have a unique number of samples, 2 if all are sampled 
				X times (I am considering 0 to be one of the number of times sampled for likely unsampled species
		t - the total number sampled
/* Returns:
		expected - an array in which e[x] gives the expected number of species sampled x times
	NOTE: sometimes this is really slow - if the maximum number of finds is really high, then use expectedfindspart and modify
		routines accordingly.  
******************************************************************************************************************************************/
double *expfinds (double *dist, int S, int mxfds, int t)
{
int j, n, m, sp;
double *expected;
double	lp=0.0000000f, lnc=0.0000000f, y=0.0000000f;

expected=dvector(mxfds+1);

/* For each possible number of finds (i.e., 0 to the maximum observed) calculate the expected number of species with n finds */
for (n=0; n<mxfds; ++n)	{
	expected[n]=0;

	/* calculate combinations in logarithms */
	lnc=0.0000000f;
	m=n;
	if (m>(t-n))	m=t-n;
	for (j=(t-m)+1; j<=t; ++j)		lnc=lnc+log(j);
	for (j=2; j<=m; ++j)			lnc=lnc-log(j);

	for (sp=0; sp<S; ++sp)	{
		y=1-dist[sp];
		lp=lnc+((double)n)*log(dist[sp])+((double)(t-n)*log(y));
		expected[n]=expected[n]+pow(e,lp);
//		expected[n]=expected[n]+binomexact(n,t,dist[sp]);
		}
	}
return expected;
}

/* expectedfindspart - finds the expected number of taxa with observed numbers finds given binomial probabilities, total finds and distribution
/* Requires:
		dist - an array giving relative abundances of S taxa (must sum to 1.0);
		S - the length of dist (i.e., richness);
		nofds - an array giving the sample numbers
		unfds - the number of unique sample numbers - this equals S+1 if all species have a unique number of samples, 2 if all are sampled 
				X times (I am considering 0 to be one of the number of times sampled for likely unsampled species
		t - the total number sampled
/* Returns:
		expected - an array in which e[x] gives the expected number of species sampled x times
	NOTE: this is done only for observed numbers of samples because it was too slow to do them all for some reason.....
******************************************************************************************************************************************/
double *expfindspart (double *dist, int S, long *nofds, int unfds, int t)
{
int i, j, n, m, sp;
double *expected;
double	lp=0.0000000f, lnc=0.0000000f, y=0.0000000f;

expected=dvector(unfds);

/* For each possible number of finds (i.e., 0 to the maximum observed) calculate the expected number of species with n finds */
for (i=0; i<unfds; ++i)	{
	n=nofds[i];
	expected[n]=0;

	/* calculate combinations in logarithms */
	lnc=0.0000000f;
	m=n;
	if (m>(t-n))	m=t-n;
	for (j=(t-m)+1; j<=t; ++j)		lnc=lnc+log(j);
	for (j=2; j<=m; ++j)			lnc=lnc-log(j);

	for (sp=0; sp<S; ++sp)	{
		y=1-dist[sp];
		lp=lnc+((double)n)*log(dist[sp])+((double)(t-n)*log(y));
		expected[i]=expected[i]+pow(e,lp);
//		expected[i]=expected[i]+binomexact(n,t,dist[sp]);
		}
	}
return expected;
}