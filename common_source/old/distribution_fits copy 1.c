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

double L = 0.0f;							/* Summed Log Likelihood */
double a = 0.0f, b = 0.0f;
int i = 0;									/* loop variable */

/* CALCULATE LOG LIKELIHOOD FOR EACH TAXON IN BOTH LISTS	*/
for (i=0; i<Z; i++) {
 
	a = log(A[i]) * S[i];					/* a = pi ^ni */
	b = log(1-A[i]) * (N-S[i]);				/* b = (1 - pi) ^ (N-ni) */
	L+=(a+b);								/* Sum L for all taxa */
	}
if (L>0) L=-1*DBL_MAX;						/* if c blows out = very bad fit */

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

double L = 0.0f;						/* Summed Log Likelihood */
double lnp = 0.0f;						/* taxon support */
double a = 0.0f, b = 0.0f, c = 0.0f;	/* temp variables */
int i;									/* loop variable */

/* CALCULATE LOG LIKELIHOOD FOR EACH TAXON WITH n > 1	*/
for (i=0; n[i]>1; i++) {
	
	a = log(A[i]) * n[i];				/* a = pi ^ni */
	b = log(1-A[i]) * (N-n[i]);			/* b = (1 - pi) ^ (N-ni) */
	c = log(1 - pow((1-A[i]),N));		/* c = 1 - ( 1 - pi ) ^ N - the probability of being sampled */
	if (c>0)
		c=0;							/* sometimes this gets wonky - blows out */
	lnp = (a + b) - c;					/* L = (a * b) / c */
	L+=lnp;								/* Sum L for all taxa */
	}
if (L>0) L=-1*DBL_MAX;					/* if c blows out = very bad fit */

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
*************************************************************************************************/
double calc_likelihood_exp(double *X, long *nufinds, long *histogram, int unqfnds)
{
double L = 0.0f;						/* Summed Log Likelihood */
double a = 0.0f;						/* temp variables */
int i, n;								/* loop variable */

for (i=1; i<unqfnds; ++i)	{
	n=nufinds[i];
	a=histogram[n];
	L = L+(a*log(X[i]));
	}

if (L>0) L=-1*DBL_MAX;					/* if c blows out = very bad fit */

return L;
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
*************************************************************************************************/
double *poi_fit_gs(int *empdist, int ntaxa, int nspec) 
{
int i = 0;					/* LOOP VARIABLE													*/
int r = 0;					/* LOOP RICHNESS													*/
int ri;						/* richness increment												*/
int rin;					/* initial richness													*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
double ev = 0.000f;			/* LOOP SLOPE 														*/
double emin = 1.0f;			/* min slope														*/
double ein = 0.000f;		/* initial slope													*/
double ei = 0.000f;			/* how much to increment ev in each loop							*/
double es[3];				/* previous log likelihoods (cell number = num previous).			*/
double rs[3];				/* previous log likelihoods (cell number = num previous).			*/
double bep[3];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/
double *brp;				/* BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array		*/
double *fitdist;			/* fit distribution													*/
double *expect;				/* expected number of species with 0Émax finds						*/
double *obsrvd;				/* observed number of species with 0Émax finds						*/
double pbes;				/* previous best support for decay rate								*/

brp=dvector(3);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
ein=pow(empdist[0],(((double) empdist[ntaxa-1])/((double) ntaxa)));
if (ein<1.0)	ein=1.0f;

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);
ri=ntaxa/2;
if (ri<2)	ri=2;
/* increment true richness until that fails to improve likelihood */
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{
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
		es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {					/* IF BETTER THAN BEST FIT */
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								/* STORE SLOPE */
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei *= -1;									/* STEP BACKWARD TO FIND PEAK */
			es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
			}
		else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei /= 10;									/* SET SMALLER UNIT */
			}
		
		while (ev+ei<= emin)	ei/=2;
		
		es[2] = es[1];									/* Store last 2 attempts to identify */
		es[1] = es[0];									/* when the peak is past */
		}
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
		brp[2] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		ein=bep[1];										/* set initial decay to best decay found so far */ 
		}
	/* optimal richness is overshot */
	else if ((rs[2]>rs[1] && (rs[0]>rs[1])) && (r-ri)>=ntaxa)	{
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
free_dvector(obsrvd);
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
*************************************************************************************************/
double *poi_fit_zm(int *empdist, int ntaxa, int nspec)
{
int i = 0;				/* LOOP VARIABLE	*/
int r = 0;					/* LOOP RICHNESS	*/
int ri;						/* richness increment	*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
int rin;					/* initial richness */
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

/* find a good seed value for the Zipf-Mandelbrot, based on the log-log slopes between the initial taxa */
ein=(log(empdist[0])/log(empdist[2]));
if (ein<emin)
	ein=pow((log(empdist[0])/log(empdist[6])),0.5);	
else if (empdist[2]==1)
	ein=pow((log(empdist[0])/log(empdist[1])),2);
if (ein<emin)	ein=emin;

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);
ri=ntaxa/2;
if (ri<2)	ri=2;

/* increment true richness until that fails to improve likelihood */
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;
	ei = (double) FITINC / 10;
	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	
	for (i=0; i<3; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = 0.0f;
	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; (ev>=emin && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC)))); ev += ei) {	
		
		/*Debugging line */
		if (ev<=emin)	{
			ev=emin+(ei/10);
			if (ev<emin)	{
				ei/=-10;
				ev=emin+ei;
				}
			}

		/* generate Zipf-Mandelbrot distribution with richnes r and decay of ev */
		fitdist = proportional_zm_distribution(ev,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,2*empdist[0],nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);

		if (es[0] >= bep[0]) {					/* IF BETTER THAN BEST FIT */
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								/* STORE SLOPE */
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei *= -1;									/* STEP BACKWARD TO FIND PEAK */
			es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
			}
		else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei /= 10;									/* SET SMALLER UNIT */
			}
		es[2] = es[1];									/* Store last 2 attempts to identify */
		es[1] = es[0];									/* when the peak is past */
		}
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
		brp[2] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		ein=bep[1];										/* set initial decay to best decay found so far */ 
		}
	/* optimal richness is overshot */
/*	else if ((rs[2]>rs[1] && (rs[0]>rs[1])) && (abs(ri)>1))	{	*/
	else if ((rs[2]>rs[1] && (rs[0]>rs[1])) && (r-ri)>=ntaxa)	{
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
	}
free_dvector(obsrvd);
return brp;
}

/* FIND THE MOST-LIKELY FAUX-LOG-NORMAL SERIES, VARYING INTIAL AND MODAL SLOPE AND THE MODE. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: initial slope of shifting geometric series
	- result[2]: modal slope of shifting geometric series
	- result[3]: optimal richness			
	- result[4]: mode (in #taxa from median taxon rank)			
COMMENTS:
	- SUPINC, FITINC defined in header file
*****************************************************************************************************/
double *poi_fit_ln(int *empdist, int ntaxa, int nspec) 
{
int j=0,i=0;				/* LOOP VARIABLE														*/
int r = 0;					/* LOOP RICHNESS														*/
int ri;						/* richness increment													*/
int lri;					/* previous richness increment											*/
int rin;					/* initial richness to use in each search (begins as ntaxa)				*/
double mode;				/* mode of log-normal, given as a taxon number							*/
double x=1.0f;				/* x is the absolute value of the mode - cannot do abs(double)			*/
double mi=-1.0f;			/* how much to increment the mode										*/
double mimin=0.125;			/* minimum mode change													*/
double mdmx=1.5f;
double mdmn=-2.5f;
double ocs=3.0f;
/*double imode;				/* initial mode (set to the median of taxa ranks)						*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
double ev = 0.000f;			/* LOOP SLOPE 															*/
double emin = 1.00f;		/* min slope															*/
double mxev=75.0f;			/* maximum magnitude change to consider									*/
double ein = 0.000f;		/* initial slope														*/
double ei = 0.000f;			/* how much to increment ev in each loop								*/
double es[3];				/* previous modal decay log likelihoods (cell number = num previous).	*/
double rs[3];				/* previous richness log likelihoods (cell number = num previous).		*/
double ms[3];				/* previous mode log likelihoods (cell number = num previous).			*/
double bep[2];				/* BEST modal decay parameters - return array format 					*/
double brp[3];				/* BEST richness parameters - return array								*/
double pbes;				/* previous best support for modal decay								*/
double pbrs;				/* previous best support for richness									*/
double pbms;				/* previous best support for mode location								*/
double *expect;				/* expected number of species with 0Émax finds							*/
double *obsrvd;				/* observed number of species with 0Émax finds							*/
double *fitdist;			/* fit distribution												*/
double *bmp;				/* Best mode parameters - return array format							*/

bmp=dvector(8);
for (i=0; i<8; i++) bmp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);

if (rin<100)	j=4;
else			j=5;
ein=pow(empdist[0],1/((double) j));

if (2*rin<75)			mdmx=2.5;
else if (2*rin<350)		mdmx=3.0;
else if (2*rin<2000)	mdmx=3.5;
else if (2*rin<15000)	mdmx=4.0;
else 					mdmx=4.5;

mdmn=-1*mdmx;

/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (i=0; i<3; ++i)	ms[i]=-1.0*DBL_MAX;
pbms=0.0f;
/* adjust mode until that fails to improve likelihood	*/
	/* for some reason the program is disobeying the second part of the conditional */
for (mode=0; x>=mimin && (mode>=mdmn && mode<=mdmx)/* && ((pbms == 0.0f) || (bmp[0] > (pbms + SUPINC)))*/; mode+=mi)	{
	pbrs=0.0f;
	for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;
	for (i=0; i<3; ++i)	rs[i]=-1.0*DBL_MAX;

	/* adjust true richness until that fails to improve likelihood	*/
	lri=ri=ntaxa/2;
	if (ri<2)	ri=2;
		
	for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{

		if (r<75)			ocs=2.5;
		else if (r<350)		ocs=3.0;
		else if (r<2000)	ocs=3.5;
		else if (r<15000)	ocs=4.0;
		else 				ocs=4.5;

		if (mode>ocs)					ocs=mode;
		else if (mode<0 && ocs<-1*mode)	ocs=-1*mode;
		
		/* make sure that mode starts between beginning and end */
		while (abs(mode/2)>r)	++r;
		
		obsrvd[0]=r-ntaxa;
		pbes = 0.0f;
		ei = ein / 10;
		/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
		
		for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
		for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;

		if (ein+ei>=mxev)	ei*=-1;
		/* increment slope until that fails to improve likelihood or resolution limit reached */
		for (ev=ein; ev<=mxev && ((pbes==0.0f) || (bep[0]>(pbes+SUPINC))); ev+=ei)	{
			
			/* generate Zipf-Mandelbrot distribution with richnes r and decay of ev */
			fitdist = proportional_lgn_distribution(ev,mode,ocs,r);		/* MAKE DISTRIBUTION */
			/* find the expected proportions of taxa with 0Éx finds */
			expect=expfinds(fitdist,r,2*empdist[0],nspec);
			es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
			free_dvector(expect);			/* memory for expect freed! */
			free_dvector(fitdist);			/* memory for expect freed! */

			if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
				pbes = bep[0];								/* save last best ssq for evenness */
				rs[0] = bep[0] = es[0];						/* STORE FIT */
				bep[1] = ev;								/* STORE SLOPE */
				}
			/* check to make sure that this is correct elsewhere */
			else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>=emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
				ev -= ei;									/* STEP BACK ONE UNIT */
				ei *= -1;									/* STEP BACKWARD TO FIND PEAK */
				es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
				}
			else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
				ev -= ei;									/* STEP BACK ONE UNIT */
				ei /= 3;									/* SET SMALLER UNIT */
				}
			/* it might go on and on forever with miniscule increases when it is a very poor fit */
			if ((es[0]>=es[1] && es[0]<(es[1]+SUPINC)) && (es[1]>=es[2] && es[1]<(es[2]+SUPINC)))	break;
			es[2] = es[1];									/* Store last 2 attempts to identify */
			es[1] = es[0];									/* when the peak is past */
			}
		/* if this richness is better than the last */
		if (bep[0]>=brp[0]) {								/* IF BETTER THAN BEST FIT */
			pbrs=brp[0];
			brp[2] = r;
			for (i=0; i<2; i++)
				brp[i] = bep[i];
			ms[0]=brp[0];
			ein=bep[1];										/* set initial decay to best decay found so far */ 
			
			lri=ri;
			}
		/* optimal richness is overshot - back up */
		else if (abs(ri)==abs(lri))	{
			r-=ri;
			lri=ri;
			ri/=2;
			}
		/* two tries in this direction produced poorer results - try the other direction */
		else if (abs(lri)>abs(ri))	{
			r-=ri;
			lri=ri;
			/* do not look at values lower than observed richness */
			if ((r-ri)>=ntaxa)	ri*=-1;
			else				ri/=2;
			}
		/* optimal richness is overshot */
		else if (rs[2]>rs[1] && rs[0]>rs[1])	{	
			r-=ri;				/* step back one unit	 */
			ri*=-1;				/* step backwards to peak */
			while ((r+ri)<ntaxa && abs(ri)>0)	ri/=2;
			rs[1]=rs[2];		/* set to prior ln L to ignore over step */
			}
		else	{
			r-=ri;						/* put r back to where it was */
			if (abs(ri)>1)	ri/=2;
			else			ri=0;
			while ((r+ri)<=ntaxa && abs(ri)>0)	ri/=2;
			}
		/* it might go on and on forever with miniscule increases when it is a very poor fit */
		if ((rs[0]>=rs[1] && rs[0]<(rs[1]+SUPINC)) && (rs[1]>=rs[2] && rs[1]<(rs[2]+SUPINC)))	ri=0;

		rs[2] = rs[1];
		rs[1] = bep[0];
		}

	/* if this mode is better than the last */
	if (brp[0] >= bmp[0]) {								/* IF BETTER THAN BEST FIT */
		pbms=bmp[0];
		for (i=0; i<3; i++)
			bmp[i] = brp[i];
		bmp[3] = mode;
		/* save the best untruncated log-normal separately */
		if (mode==0)	{
			for (i=0; i<3; ++i)	bmp[i+4]=brp[i];
			bmp[7]=0;
			}
		}
	/* optimal richness is overshot */
	else if (ms[2]>ms[1] && ms[0]>ms[1])	{
		mode-=mi;				/* step back one unit	 */
		mi*=-1;					/* step backwards to peak */
		ms[1]=ms[2];		/* set to prior ln L to ignore over step */
		}
	else	{
		mode-=mi;
		mi/=2;
		}

	/* it might go on and on forever with miniscule increases when it is a very poor fit */
	if ((ms[0]>=ms[1] && ms[0]<(ms[1]+SUPINC)) && (ms[1]>=ms[2] && ms[1]<(ms[2]+SUPINC)))	ri=0;
	ms[2] = ms[1];
	ms[1] = brp[0];

	/* do not overshoot mode limits! */
	while ((mode+mi<mdmn || mode+mi>mdmx) && x>mimin)	{
		mi/=2;
		if (mi<1)	x=-1*mi;
		}
	
	ein=bmp[1];										/* set initial modal decay to best modal decay found so far 	*/ 
	rin=bmp[2];										/* set initial richness to best richness found so far			*/
	x=mi;
	if (x<0)	x*=-1;
	}

free_dvector(obsrvd);
return bmp;
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
*************************************************************************************************/
double *poi_fit_fln_1(int *empdist, int ntaxa, int nspec) 
{
int i = 0;					/* LOOP VARIABLE	*/
int r = 0;					/* LOOP RICHNESS	*/
int ri;						/* richness increment	*/
int rin;					/* initial richness	*/
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

brp=dvector(5);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<4; i++) brp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the initial slope, based on the slope between the first and second species */
ein=((double) empdist[0])/((double) empdist[1]);

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);
ri=ntaxa/2;
if (ri<2)	ri=2;

/* increment true richness until that fails to improve likelihood */
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{
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
		es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: fLN R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								/* STORE SLOPE */
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei *= -1;									/* STEP BACKWARD TO FIND PEAK */
			es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
			}
		else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei /= 10;									/* SET SMALLER UNIT */
			}
		es[2] = es[1];									/* Store last 2 attempts to identify */
		es[1] = es[0];									/* when the peak is past */
		}
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
		brp[4] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		if (r%2==1)	brp[3] = r/2;						/* mode is always the median in this case		*/
		else		brp[3] = (r/2)-0.5;
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

	/* it might go on and on forever with miniscule increases when it is a very poor fit */
	if ((rs[0]>=rs[1] && rs[0]<(rs[1]+SUPINC)) && (rs[1]>=rs[2] && rs[1]<(rs[2]+SUPINC)))	ri=0;
	rs[2] = rs[1];
	rs[1] = bep[0];
	}
brp[2]=1.0f;
free_dvector(obsrvd);
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
*************************************************************************************************/
double *poi_fit_fln_2(int *empdist, int ntaxa, int nspec) 
{
int i = 0;				/* LOOP VARIABLE	*/
int r = 0;					/* LOOP RICHNESS	*/
int ri;						/* richness increment	*/
int rin;					/* initial richness increment	*/
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

brp=dvector(5);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<5; i++) brp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
ein = emin;		/* begin searching for just the shift parameter */
din=((double) empdist[0])/((double) empdist[1]);
if (din<ein)	din=ein+0.01;

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);
ri=ntaxa/2;
if (ri<2)	ri=2;
/* increment true richness until that fails to improve likelihood */
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{
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
			ds[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
			free_dvector(expect);

			if (ds[0] >= bdp[0]) {						/* IF BETTER THAN BEST FIT */
				pbds = bdp[0];									/* save last best ssq for evenness */
				es[0]=bdp[0]=ds[0];								/* STORE FIT */
				bdp[1] = dev;									/* STORE shift in slope	*/
				}
			else if ((ds[2] > ds[1]) && (ds[0] > ds[1]) && ((dev-di)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
				dev -= di;									/* STEP BACK ONE UNIT */
				di *= -1;									/* STEP BACKWARD TO FIND PEAK */
				ds[1] = ds[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
				}
			else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
				dev -= di;									/* STEP BACK ONE UNIT */
				di /= 10;									/* SET SMALLER UNIT */
				}
			/* it might go on and on forever with miniscule increases when it is a very poor fit */
			if ((ds[0]>=ds[1] && ds[0]<(ds[1]+SUPINC)) && (ds[1]>=ds[2] && ds[1]<(ds[2]+SUPINC)))	break;

			ds[2] = ds[1];									/* Store last 2 attempts to identify */
			ds[1] = ds[0];									/* when the peak is past */

			}

		if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0] = bep[0] = es[0];						/* STORE FIT */
			din = bep[1] = bdp[1];
			bep[2] = ev;								/* STORE SLOPE */
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei *= -1;									/* STEP BACKWARD TO FIND PEAK */
			es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
			}
		else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei /= 10;									/* SET SMALLER UNIT */
			}

		/* it might go on and on forever with miniscule increases when it is a very poor fit */
		if ((es[0]>=es[1] && es[0]<(es[1]+SUPINC)) && (es[1]>=es[2] && es[1]<(es[2]+SUPINC)))	break;
		es[2] = es[1];									/* Store last 2 attempts to identify */
		es[1] = es[0];									/* when the peak is past */

		}
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
		brp[4] = r;
		for (i=0; i<3; i++)
			brp[i] = bep[i];
		din=bep[1];										/* set initial decay to best decay found so far */ 
		ein=bep[2];										/* set initial decay to best decay found so far */ 
		if (r%2==1)	brp[3] = r/2;						/* mode is always the median in this case		*/
		else		brp[3] = (r/2)-0.5;
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
	- result[1]: initial slope of shifting geometric series
	- result[2]: modal slope of shifting geometric series
	- result[3]: optimal richness			
	- result[4]: mode (in #taxa from median taxon rank)			
COMMENTS:
	- SUPINC, FITINC defined in header file
*****************************************************************************************************/
double *poi_fit_fln_3(int *empdist, int ntaxa, int nspec) 
{
int i = 0;					/* LOOP VARIABLE														*/
int r = 0;					/* LOOP RICHNESS														*/
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
for (i=0; i<5; i++) bmp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
ein = emin;		/* begin searching for just the shift parameter */
din=((double) empdist[0])/((double) empdist[1]);
if (din<ein)	din=ein+0.01;

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);

/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
pbms=0.0f;
/* adjust mode until that fails to improve likelihood	*/
for (mode=0; abs(mi)>0 && ((pbms == 0.0f) || (bmp[0] > (pbms + SUPINC))); mode+=mi)	{
	pbrs=0.0f;
	for (i=0; i<4; i++) brp[i]= -1.0*DBL_MAX;
	for (i=0; i<3; ++i)	rs[i]=-1.0*DBL_MAX;

	/* adjust true richness until that fails to improve likelihood	*/
	ri=ntaxa/2;
	if (ri<2)	ri=2;
	for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{

		/* make sure that mode starts between beginning and end */
		while (abs(mode/2)>r)	++r;
		
		if (r%2==1)	imode = r/2;
		else		imode = (r/2)-0.5;

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
				fitdist = proportional_fln_distribution(ev,dev,r,imode+mode);				/* MAKE DISTRIBUTION */
				/* find the expected proportions of taxa with 0Éx finds */
				expect=expfinds(fitdist,r,2*empdist[0],nspec);
				free_dvector(fitdist);
				ds[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
				free_dvector(expect);

				if (ds[0] >= bdp[0]) {						/* IF BETTER THAN BEST FIT */
					pbds = bdp[0];									/* save last best ssq for evenness */
					es[0]=bdp[0]=ds[0];								/* STORE FIT */
					bdp[1] = dev;									/* STORE shift in slope	*/
					}
				else if ((ds[2] > ds[1]) && (ds[0] > ds[1]) && ((dev-di)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
					dev -= di;									/* STEP BACK ONE UNIT */
					di *= -1;									/* STEP BACKWARD TO FIND PEAK */
					ds[1] = ds[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
					}
				else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
					dev -= di;									/* STEP BACK ONE UNIT */
					di /= 10;									/* SET SMALLER UNIT */
					}

				/* it might go on and on forever with miniscule increases when it is a very poor fit */
				if ((ds[0]>=ds[1] && ds[0]<(ds[1]+SUPINC)) && (ds[1]>=ds[2] && ds[1]<(ds[2]+SUPINC)))	break;
				ds[2] = ds[1];									/* Store last 2 attempts to identify */
				ds[1] = ds[0];									/* when the peak is past */
				}

			if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
				pbes = bep[0];								/* save last best ssq for evenness */
				rs[0] = bep[0] = es[0];						/* STORE FIT */
				din = bep[1] = bdp[1];
				bep[2] = ev;								/* STORE SLOPE */
				}
			else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
				ev -= ei;									/* STEP BACK ONE UNIT */
				ei *= -1;									/* STEP BACKWARD TO FIND PEAK */
				es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
				}
			else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
				ev -= ei;									/* STEP BACK ONE UNIT */
				ei /= 10;									/* SET SMALLER UNIT */
				}
			/* it might go on and on forever with miniscule increases when it is a very poor fit */
			if ((es[0]>=es[1] && es[0]<(es[1]+SUPINC)) && (es[1]>=es[2] && es[1]<(es[2]+SUPINC)))	break;
			es[2] = es[1];									/* Store last 2 attempts to identify */
			es[1] = es[0];									/* when the peak is past */
			}
		/* if this richness is better than the last */
		if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
			brp[4] = r;
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
		/* it might go on and on forever with miniscule increases when it is a very poor fit */
		if ((rs[0]>=rs[1] && rs[0]<(rs[1]+SUPINC)) && (rs[1]>=rs[2] && rs[1]<(rs[2]+SUPINC)))	ri=0;

		rs[2] = rs[1];
		rs[1] = bep[0];
		}

	/* if this mode is better than the last */
	if (brp[0] >= bmp[0]) {								/* IF BETTER THAN BEST FIT */
		for (i=0; i<5; i++)
			bmp[i] = brp[i];
		bmp[3] = mode;
		din=bmp[1];										/* set initial initial slope to best initial slope found so far */ 
		ein=bmp[2];										/* set initial modal decay to best modal decay found so far 	*/ 
		rin=bmp[4];										/* set initial richness to best richness found so far			*/
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

	/* it might go on and on forever with miniscule increases when it is a very poor fit */
	if ((ms[0]>=ms[1] && ms[0]<(ms[1]+SUPINC)) && (ms[1]>=ms[2] && ms[1]<(ms[2]+SUPINC)))	ri=0;
	ms[2] = ms[1];
	ms[1] = brp[0];
	
	}

free_dvector(obsrvd);
return bmp;
}

/* FIND THE MOST-LIKELY FAUX-LOG-NORMAL SERIES, VARYING INTIAL AND MODAL SLOPE AND THE MODE. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: initial slope of shifting geometric series
	- result[2]: modal slope of shifting geometric series
	- result[3]: optimal richness			
	- result[4]: mode (in #taxa from median taxon rank)			
COMMENTS:
	- SUPINC, FITINC defined in header file
*****************************************************************************************************/
double *poi_fit_fln_3b(int *empdist, int ntaxa, int nspec, double ge, int gr)
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
double emin = 0.75f;		/* min modal slope														*/
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
for (i=0; i<5; i++) bmp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* use previously calculated geometric for seed values */
rin=chao2(empdist,ntaxa);
din=ein = ge;


/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
pbms=0.0f;
/* adjust mode until that fails to improve likelihood	*/
for (mode=0; abs(mi)>0 && ((pbms == 0.0f) || (bmp[0] > (pbms + SUPINC))); mode+=mi)	{
	pbrs=0.0f;
	for (i=0; i<4; i++) brp[i]= -1.0*DBL_MAX;
	for (i=0; i<3; ++i)	rs[i]=0;

	/* adjust true richness until that fails to improve likelihood	*/
	ri=ntaxa/2;
	if (ri<2)	ri=2;
	for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{

		/* make sure that mode starts between beginning and end */
		while (abs(mode/2)>r)	++r;
		
		if (r%2==1)	imode = r/2;
		else		imode = (r/2)-0.5;

		obsrvd[0]=r-ntaxa;
		pbes = 0.0f;
		ei = -1.0 * (double) FITINC / 10;
		/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
		
		for (i=0; i<4; i++) bep[i] = -1.0*DBL_MAX;
		for (i=0; i<3; i++) es[i]=0.0f;

		/* increment slope until that fails to improve likelihood or resolution limit reached */
		for (ev=ein; ev>=emin && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))); ev += ei) {	
			
			for (i=0; i<3; i++) bdp[i]=-1.0*DBL_MAX;
			for (i=0; i<3; i++) ds[i]=0.0f;

			di = (double) FITINC / 10;
			/*Debugging line */
	/*		if (ev<=emin) printf("\nDANGER: fln R=%d, ev=%f, S=%f ",r,ev,es[0]);	*/

			pbds = 0.0f;
			for (dev=din; ((pbds==0.0f) || (bdp[0]>(pbds+SUPINC))); dev+=di)	{
				/* generate geometric distribution with richnes r and decay of ev */
				fitdist = proportional_fln_distribution(ev,dev,r,imode+mode);				/* MAKE DISTRIBUTION */
				/* find the expected proportions of taxa with 0Éx finds */
				expect=expfinds(fitdist,r,2*empdist[0],nspec);
				free_dvector(fitdist);
				ds[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
				free_dvector(expect);

				if (ds[0] >= bdp[0]) {						/* IF BETTER THAN BEST FIT */
					pbds = bdp[0];									/* save last best ssq for evenness */
					es[0]=bdp[0]=ds[0];								/* STORE FIT */
					bdp[1] = dev;									/* STORE shift in slope	*/
					}
				else if ((ds[2] > ds[1]) && (ds[0] > ds[1]) && ((dev-di)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
					dev -= di;									/* STEP BACK ONE UNIT */
					di *= -1;									/* STEP BACKWARD TO FIND PEAK */
					ds[1] = ds[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
					}
				else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
					dev -= di;									/* STEP BACK ONE UNIT */
					di /= 10;									/* SET SMALLER UNIT */
					}

				/* it might go on and on forever with miniscule increases when it is a very poor fit */
				if ((ds[0]>=ds[1] && ds[0]<(ds[1]+SUPINC)) && (ds[1]>=ds[2] && ds[1]<(ds[2]+SUPINC)))	break;
				ds[2] = ds[1];									/* Store last 2 attempts to identify */
				ds[1] = ds[0];									/* when the peak is past */
				}

			if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
				pbes = bep[0];								/* save last best ssq for evenness */
				rs[0] = bep[0] = es[0];						/* STORE FIT */
				din = bep[1] = bdp[1];
				bep[2] = ev;								/* STORE SLOPE */
				}
			else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
				ev -= ei;									/* STEP BACK ONE UNIT */
				ei *= -1;									/* STEP BACKWARD TO FIND PEAK */
				es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
				}
			else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
				ev -= ei;									/* STEP BACK ONE UNIT */
				ei /= 10;									/* SET SMALLER UNIT */
				}
			while ((ev+ei)<emin && ei>0.00001)	ei /= 10;
			
			/* it might go on and on forever with miniscule increases when it is a very poor fit */
			if ((es[0]>=es[1] && es[0]<(es[1]+SUPINC)) && (es[1]>=es[2] && es[1]<(es[2]+SUPINC)))	break;
			es[2] = es[1];									/* Store last 2 attempts to identify */
			es[1] = es[0];									/* when the peak is past */
			}
		/* if this richness is better than the last */
		if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
			brp[4] = r;
			for (i=0; i<3; i++)
				brp[i] = bep[i];
			ms[0]=brp[0];
			din=bep[1];										/* set initial decay to best decay found so far */ 
			ein=bep[2];										/* set initial decay to best decay found so far */ 
			if (ein<1.0)	ein+=0.1;
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
		}

	/* if this mode is better than the last */
	if (brp[0] >= bmp[0]) {								/* IF BETTER THAN BEST FIT */
		for (i=0; i<5; i++)
			bmp[i] = brp[i];
		bmp[3] = mode;
		din=bmp[1];										/* set initial initial slope to best initial slope found so far */ 
		ein=bmp[2];										/* set initial modal decay to best modal decay found so far 	*/ 
		rin=bmp[4];										/* set initial richness to best richness found so far			*/
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

	/* it might go on and on forever with miniscule increases when it is a very poor fit */
	if ((ms[0]>=ms[1] && ms[0]<(ms[1]+SUPINC)) && (ms[1]>=ms[2] && ms[1]<(ms[2]+SUPINC)))	mi=0;
	ms[2] = ms[1];
	ms[1] = brp[0];
	
	}

free_dvector(obsrvd);
return bmp;
}


/* FIND THE MOST-LIKELY FAUX-LOG-NORMAL SERIES, VARYING INTIAL AND MODAL SLOPE AND THE MODE. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: initial slope of shifting geometric series
	- result[2]: modal slope of shifting geometric series
	- result[3]: optimal richness			
	- result[4]: mode (in #taxa from median taxon rank)			
COMMENTS:
	- SUPINC, FITINC defined in header file
*****************************************************************************************************/
double *poi_fit_fln_3c(int *empdist, int ntaxa, int nspec)
{
int i = 0;					/* LOOP VARIABLE														*/
int r = 0;					/* LOOP RICHNESS														*/
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
double emin = 0.75f;		/* min modal slope														*/
double ein = 1.000f;		/* initial slope														*/
double ei = 0.000f;			/* how much to increment ev in each loop								*/
double ds[3];				/* previous initial decay log likelihoods (cell number = num previous).	*/
double es[3];				/* previous modal decay log likelihoods (cell number = num previous).	*/
double rs[3];				/* previous richness log likelihoods (cell number = num previous).		*/
double ms[3];				/* previous mode log likelihoods (cell number = num previous).			*/
double bdp[3];				/* BEST initial decay parameters - return array format 					*/
double bep[4];				/* BEST modal decay parameters - return array format 					*/
double brp[5];				/* BEST richness parameters - return array								*/
double *bmp;				/* Best mode parameters - return array format							*/
double *b2p;				/* Best two parameter (mode = median and mid-slope=1.0					*/
double *fitdist;			/* fit distribution														*/
double *expect;				/* expected number of species with 0Émax finds							*/
double *obsrvd;				/* observed number of species with 0Émax finds							*/
double pbds;				/* previous best support for initial decay								*/
double pbes;				/* previous best support for modal decay								*/
double pbrs;				/* previous best support for richness									*/
double pbms;				/* previous best support for mode location								*/

bmp=dvector(5);
for (i=0; i<5; i++) bmp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* use previously calculated geometric for seed values */
b2p=poi_fit_fln_1(empdist,ntaxa,nspec);

/* begin searching for just the shift parameter */
din=b2p[1];
rin=b2p[4];

/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
pbms=0.0f;
/* adjust mode until that fails to improve likelihood	*/
for (mode=0; abs(mi)>0 && ((pbms == 0.0f) || (bmp[0] > (pbms + SUPINC))); mode+=mi)	{
	pbrs=0.0f;
	for (i=0; i<4; i++) brp[i]= -1.0*DBL_MAX;
	for (i=0; i<3; ++i)	rs[i]=0;

	/* adjust true richness until that fails to improve likelihood	*/
	ri=ntaxa/2;
	if (ri<2)	ri=2;
	for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{

		/* make sure that mode starts between beginning and end */
		while (abs(mode/2)>r)	++r;
		
		if (r%2==1)	imode = r/2;
		else		imode = (r/2)-0.5;

		obsrvd[0]=r-ntaxa;
		pbes = 0.0f;
		ei = -1.0 * (double) FITINC / 10;
		/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
		
		for (i=0; i<4; i++) bep[i] = -1.0*DBL_MAX;
		for (i=0; i<3; i++) es[i]=0.0f;

		/* increment slope until that fails to improve likelihood or resolution limit reached */
		for (ev=ein; ev>=emin && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))); ev += ei) {	
			
			for (i=0; i<3; i++) bdp[i]=-1.0*DBL_MAX;
			for (i=0; i<3; i++) ds[i]=0.0f;

			di = (double) FITINC / 10;
			/*Debugging line */
	/*		if (ev<=emin) printf("\nDANGER: fln R=%d, ev=%f, S=%f ",r,ev,es[0]);	*/

			pbds = 0.0f;
			for (dev=din; ((pbds==0.0f) || (bdp[0]>(pbds+SUPINC))); dev+=di)	{
				/* generate geometric distribution with richnes r and decay of ev */
				fitdist = proportional_fln_distribution(ev,dev,r,imode+mode);				/* MAKE DISTRIBUTION */
				/* find the expected proportions of taxa with 0Éx finds */
				expect=expfinds(fitdist,r,2*empdist[0],nspec);
				free_dvector(fitdist);
				ds[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
				free_dvector(expect);

				if (ds[0] >= bdp[0]) {						/* IF BETTER THAN BEST FIT */
					pbds = bdp[0];									/* save last best ssq for evenness */
					es[0]=bdp[0]=ds[0];								/* STORE FIT */
					bdp[1] = dev;									/* STORE shift in slope	*/
					}
				else if ((ds[2] > ds[1]) && (ds[0] > ds[1]) && ((dev-di)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
					dev -= di;									/* STEP BACK ONE UNIT */
					di *= -1;									/* STEP BACKWARD TO FIND PEAK */
					ds[1] = ds[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
					}
				else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
					dev -= di;									/* STEP BACK ONE UNIT */
					di /= 10;									/* SET SMALLER UNIT */
					}

				/* it might go on and on forever with miniscule increases when it is a very poor fit */
				if ((ds[0]>=ds[1] && ds[0]<(ds[1]+SUPINC)) && (ds[1]>=ds[2] && ds[1]<(ds[2]+SUPINC)))	break;
				ds[2] = ds[1];									/* Store last 2 attempts to identify */
				ds[1] = ds[0];									/* when the peak is past */
				}

			if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
				pbes = bep[0];								/* save last best ssq for evenness */
				rs[0] = bep[0] = es[0];						/* STORE FIT */
				din = bep[1] = bdp[1];
				bep[2] = ev;								/* STORE SLOPE */
				}
			else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>=emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
				ev -= ei;									/* STEP BACK ONE UNIT */
				ei *= -1;									/* STEP BACKWARD TO FIND PEAK */
				es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
				}
			else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
				ev -= ei;									/* STEP BACK ONE UNIT */
				ei /= 10;									/* SET SMALLER UNIT */
				}
			while ((ev+ei)<emin && ei>0.00001)	ei /= 10;
			
			/* it might go on and on forever with miniscule increases when it is a very poor fit */
			if ((es[0]>=es[1] && es[0]<(es[1]+SUPINC)) && (es[1]>=es[2] && es[1]<(es[2]+SUPINC)))	break;
			es[2] = es[1];									/* Store last 2 attempts to identify */
			es[1] = es[0];									/* when the peak is past */
			}
		/* if this richness is better than the last */
		if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
			brp[4] = r;
			for (i=0; i<3; i++)
				brp[i] = bep[i];
			ms[0]=brp[0];
			din=bep[1];										/* set initial decay to best decay found so far */ 
			ein=bep[2];										/* set initial decay to best decay found so far */ 
			if (ein<1.0)	ein+=0.1;
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
		}

	/* if this mode is better than the last */
	if (brp[0] >= bmp[0]) {								/* IF BETTER THAN BEST FIT */
		for (i=0; i<5; i++)
			bmp[i] = brp[i];
		bmp[3] = mode;
		din=bmp[1];										/* set initial initial slope to best initial slope found so far */ 
		ein=bmp[2];										/* set initial modal decay to best modal decay found so far 	*/ 
		rin=bmp[4];										/* set initial richness to best richness found so far			*/
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

	/* it might go on and on forever with miniscule increases when it is a very poor fit */
	if ((ms[0]>=ms[1] && ms[0]<(ms[1]+SUPINC)) && (ms[1]>=ms[2] && ms[1]<(ms[2]+SUPINC)))	mi=0;
	ms[2] = ms[1];
	ms[1] = brp[0];
	
	}

free_dvector(obsrvd);
return bmp;
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
*****************************************************************************************************/
double *wht_fit_gs(int *empdist, int ntaxa, int nspec) 
{
int i = 0;					/* LOOP VARIABLE														*/
int r = 0;					/* LOOP RICHNESS														*/
int ri;						/* richness increment													*/
int	rin=ntaxa;				/* initial seed richness 												*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
double ev = 0.000f;			/* LOOP SLOPE 															*/
double emin = 1.0f;			/* min slope															*/
double ein = 0.000f;		/* initial slope														*/
double ei = 0.000f;			/* how much to increment ev in each loop								*/
double es[3];				/* previous log likelihoods (cell number = num previous).				*/
double rs[3];				/* previous log likelihoods (cell number = num previous).				*/
double bep[3];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format		*/
double *brp;				/* BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array			*/
double *fitdist;			/* fit distribution														*/
double pbes;				/* previous best support for decay rate									*/

brp=dvector(3);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;

/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
ein=pow(empdist[0],(((double) empdist[ntaxa-1])/((double) ntaxa)));

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);
ri=ntaxa/2;
if (ri<2)	ri=2;
/* increment true richness until that fails to improve likelihood */
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{
	pbes = 0.0f;
	ei = (double) FITINC / 10;
	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	
	for (i=0; i<3; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = 0.0f;
	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; (ev>=emin && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC)))); ev += ei) {	
		
		/* generate geometric distribution with richness r and decay of ev */
		fitdist = proportional_gs_distribution(ev,r);				/* MAKE DISTRIBUTION */
		es[0] = calc_likelihood_Foote(fitdist, empdist, nspec);		/* CALCULATE SUPPORTt */
		free_dvector(fitdist);

		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {					/* IF BETTER THAN BEST FIT */
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								/* STORE SLOPE */
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei *= -1;									/* STEP BACKWARD TO FIND PEAK */
			es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
			}
		else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei /= 10;									/* SET SMALLER UNIT */
			}
		
		while (ev+ei<= emin)	ei/=2;
		
		es[2] = es[1];									/* Store last 2 attempts to identify */
		es[1] = es[0];									/* when the peak is past */
		}
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
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
double *wht_fit_zm(int *empdist, int ntaxa, int nspec)
{
int i = 0;				/* LOOP VARIABLE	*/
int r = 0;					/* LOOP RICHNESS	*/
int ri;						/* richness increment	*/
int	rin=ntaxa;				/* initial seed richness */
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
double ev = 0.000f;			/* LOOP SLOPE 		*/
double emin = 1.0f;			/* min slope		*/
double ein = 1.000f;		/* initial slope	*/
double ei = 0.000f;			/* how much to increment ev in each loop							*/
double es[3];				/* previous log likelihoods (cell number = num previous).		*/
double rs[3];				/* previous log likelihoods (cell number = num previous).		*/
double bep[3];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/
double *brp;				/* BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array	*/
double *fitdist;			/* fit distribution												*/
double pbes;				/* previous best support for decay rate							*/

brp=dvector(3);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;

/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
ein=log(ntaxa-1)/(log(empdist[0])-log(empdist[ntaxa-1]));

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);
ri=ntaxa/2;
if (ri<2)	ri=2;
/* increment true richness until that fails to improve likelihood */
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{
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
		es[0] = calc_likelihood_Foote(fitdist, empdist, nspec);	/* CALCULATE SUPPORTt */
		free_dvector(fitdist);

		if (es[0] >= bep[0]) {					/* IF BETTER THAN BEST FIT */
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								/* STORE SLOPE */
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei *= -1;									/* STEP BACKWARD TO FIND PEAK */
			es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
			}
		else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei /= 10;									/* SET SMALLER UNIT */
			}
		es[2] = es[1];									/* Store last 2 attempts to identify */
		es[1] = es[0];									/* when the peak is past */
		}
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
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
	- result[1]: initial slope of shifting geometric series
	- result[2]: modal slope of shifting geometric series
	- result[3]: optimal richness			
	- result[4]: mode (in #taxa from median taxon rank)			
COMMENTS:
	- SUPINC, FITINC defined in header file
***********************************************************************/
double *wht_fit_fln_3(int *empdist, int ntaxa, int nspec) 
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
double pbds;				/* previous best support for initial decay								*/
double pbes;				/* previous best support for modal decay								*/
double pbrs;				/* previous best support for richness									*/
double pbms;				/* previous best support for mode location								*/

bmp=dvector(5);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<5; i++) bmp[i]= -1.0*DBL_MAX;


/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
ein = emin;		/* begin searching for just the shift parameter */
din=((double) empdist[0])/((double) empdist[1]);
if (din<ein)	din=ein+0.01;

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);

pbms=0.0f;
/* adjust mode until that fails to improve likelihood	*/
for (mode=0; abs(mi)>0 && ((pbms == 0.0f) || (bmp[0] > (pbms + SUPINC))); mode+=mi)	{
	pbrs=0.0f;
	for (i=0; i<4; i++) brp[i]= -1.0*DBL_MAX;
	/* adjust true richness until that fails to improve likelihood	*/
	ri=ntaxa/2;
	if (ri<2)	ri=2;
	for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{

		/* make sure that mode starts between beginning and end */
		while (abs(mode/2)>r)	++r;
		
		if (r%2==1)	imode = r/2;
		else		imode = (r/2)-0.5;

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
				fitdist = proportional_fln_distribution(ev,dev,r,imode+mode);		/* MAKE DISTRIBUTION */
				ds[0] = calc_likelihood_Foote(fitdist, empdist, nspec);				/* CALCULATE SUPPORTt */

				if (ds[0] >= bdp[0]) {							/* IF BETTER THAN BEST FIT */
					pbds = bdp[0];								/* save last best ssq for evenness */
					es[0]=bdp[0]=ds[0];							/* STORE FIT */
					bdp[1] = dev;								/* STORE shift in slope	*/
					}
				else if ((ds[2] > ds[1]) && (ds[0] > ds[1]) && ((dev-di)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
					dev -= di;									/* STEP BACK ONE UNIT */
					di *= -1;									/* STEP BACKWARD TO FIND PEAK */
					ds[1] = ds[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
					}
				else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
					dev -= di;									/* STEP BACK ONE UNIT */
					di /= 10;									/* SET SMALLER UNIT */
					}

				/* it might go on and on forever with miniscule increases when it is a very poor fit */
				if ((ds[0]>=ds[1] && ds[0]<(ds[1]+SUPINC)) && (ds[1]>=ds[2] && ds[1]<(ds[2]+SUPINC)))	break;
				ds[2] = ds[1];									/* Store last 2 attempts to identify */
				ds[1] = ds[0];									/* when the peak is past */
				}

			if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
				pbes = bep[0];								/* save last best ssq for evenness */
				rs[0] = bep[0] = es[0];						/* STORE FIT */
				din = bep[1] = bdp[1];
				bep[2] = ev;								/* STORE SLOPE */
				}
			else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
				ev -= ei;									/* STEP BACK ONE UNIT */
				ei *= -1;									/* STEP BACKWARD TO FIND PEAK */
				es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
				}
			else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
				ev -= ei;									/* STEP BACK ONE UNIT */
				ei /= 10;									/* SET SMALLER UNIT */
				}
			/* it might go on and on forever with miniscule increases when it is a very poor fit */
			if ((es[0]>=es[1] && es[0]<(es[1]+SUPINC)) && (es[1]>=es[2] && es[1]<(es[2]+SUPINC)))	break;
			es[2] = es[1];									/* Store last 2 attempts to identify */
			es[1] = es[0];									/* when the peak is past */
			}
		/* if this richness is better than the last */
		if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
			pbrs=brp[0];									/* save last best fit */
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
		/* it might go on and on forever with miniscule increases when it is a very poor fit */
		if ((rs[0]>=rs[1] && rs[0]<(rs[1]+SUPINC)) && (rs[1]>=rs[2] && rs[1]<(rs[2]+SUPINC)))	ri=0;
		rs[2] = rs[1];
		rs[1] = bep[0];
		}

	/* if this mode is better than the last */
	if (brp[0] >= bmp[0]) {								/* IF BETTER THAN BEST FIT */
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
		ms[1]=ms[2];			/* set to prior ln L to ignore over step */
		}
	else	{
		mode-=mi;
		if (abs(mi)>1)	mi/=2;
		else			mi=0;
		}
	/* it might go on and on forever with miniscule increases when it is a very poor fit */
	if ((ms[0]>=ms[1] && ms[0]<(ms[1]+SUPINC)) && (ms[1]>=ms[2] && ms[1]<(ms[2]+SUPINC)))	mi=0;
	ms[2] = ms[1];
	ms[1] = brp[0];
	}

return bmp;
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
int	rin=ntaxa;				/* initial seed richness */
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

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);
ri=ntaxa/2;
if (ri<2)	ri=2;
/* increment true richness until that fails to improve likelihood */
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (r=rin; abs(ri)>0; r+=ri)	{
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
		expect=expfinds(fitdist,r,2*empdist[0],nspec);
		free_dvector(fitdist);
		es[0] = -1*dsumsqdiffs(expect,obsrvd,1+empdist[0]);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=(emin-SUPINC)) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								/* STORE SLOPE */
			}
		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei *= -1;									/* STEP BACKWARD TO FIND PEAK */
			es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */
			}
		else {											/* NOT IMPROVING TRY SMALLER INCREMENT */
			ev -= ei;									/* STEP BACK ONE UNIT */
			ei /= 10;									/* SET SMALLER UNIT */
			}
		es[2] = es[1];									/* Store last 2 attempts to identify */
		es[1] = es[0];									/* when the peak is past */
		}
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
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
		}
	}
return expected;
}


/* expectedfindsln - finds the expected number of taxa with 0...max finds given binomial probabilities for a lognormal distribution
/* Requires:
		mag - an magnitude of increase for each octave;
		mode - where the mode is, given as a taxon number (median of 1ÉS if untruncated);
		S - where the mode is, given as a taxon number (median of 1ÉS if untruncated);
		mxfds - maximum finds with which to bother
/* Returns:
		expected - an array in which e[x] gives the expected number of species sampled x times
	NOTE: sometimes this is really slow - if the maximum number of finds is really high, then use expectedfindspart and modify
		routines accordingly.  
******************************************************************************************************************************************/
double *expfindsln (double mag, double mode, int S, int mxfds, int t)
{
int		j, m, n, sp, mxsp=0;
double	oct, moct, x, rnd=0.0f, ttl=0.0f;
double	lp=0.0000000f, lnc=0.0000000f, y=0.0000000f;
double	*fn;
double	*expected;

expected=dvector(mxfds+1);

x=S;
for (n=1; x>=normheight(0,0,1); ++n)	{
	x=S*normheight(n,0,1);
	}
	
/* make sure that mxsp>=mxfds */
/* if S<100, then there are 5 octaves that should have species */
if (S<100)	{
	mxsp=pow(mag,5);
/*	moct=mode + (((double) n) - 0.5);	*/
	moct=2.5;
	}
/* if S>100, then there are 6 octaves that should have species */
/* if S>1000, then there are 7 octaves that should have species - however, that 
		blows out memory, so we have to get by on 6 */
else	{
	mxsp=pow(mag,6);
	moct=3.0f;
	}
if (mxsp<mxfds)	mxsp=mxfds;

fn=dvector(mxsp+1);

/* find the probability of species beginning with n specimens*/
for (n=1; n<=mxsp; ++n)	{
	/* calculate how far along the x-axis we are given that the axis really is in a log-scale */
	oct=log(n)/log(mag);
	x=0.5+n;
	/* find the area within the histogram centred on n - that is the expected frequency of species with n specimens */
	if (n>1)	{
		fn[n]=((double) S)*normheight(oct,moct,1)*((log(x)/log(mag))-(log(x-1)/log(mag)));
		}
	else	{
		fn[n]=((double) S)*normheight(oct,moct,1)*(log(x)/log(mag));
		}
	/* estimate the total number of expected individuals */
	ttl+=fn[n]*((double) n);
	}

/* now calculate the expected number of species with n specimens 	*/
/* this is P[n | fn] x E[fn]										*/
for (n=0; n<=mxfds; ++n)	{
	expected[n]=0;

	/* calculate combinations in logarithms */
	lnc=0.0000000f;
	m=n;
	if (m>(t-n))	m=t-n;
	for (j=(t-m)+1; j<=t; ++j)		lnc=lnc+log(j);
	for (j=2; j<=m; ++j)			lnc=lnc-log(j);

	for (sp=1; sp<=mxsp; ++sp)	{
		y=1-(((double) sp)/ttl);	/* y = probability of not sampling the taxon */
		lp=lnc+((double)n)*log(((double) sp)/ttl)+((double)(t-n)*log(y));
		/* add a conditional probability because we do not expect simply one species */
		lp+=log(fn[sp]);
		expected[n]=expected[n]+pow(e,lp);
		}
	rnd+=expected[n];
	}

for (n=0; n<=mxfds; ++n)	expected[n]*=((double) S)/rnd;

free_dvector(fn);
return expected;
}