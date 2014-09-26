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
#include <time.h>
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

double *poi_fit_uni(int *empdist, int ntaxa, int nspec)
{
int i = 0;					/* LOOP VARIABLE													*/
int r = 0;					/* LOOP RICHNESS													*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
int rin;					/* initial slope													*/
int ri;						/* how much to increment ev in each loop							*/
int lri[2];					/* last two slope increments										*/
int iri;					/* initial slope increment at each richness							*/
double rs[3];				/* previous log likelihoods (cell number = num previous).			*/
double *brp;				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/
double *fitdist;			/* fit distribution													*/
double *expect;				/* expected number of species with 0Émax finds						*/
double *obsrvd;				/* observed number of species with 0Émax finds						*/
double pbes;				/* previous best support for decay rate								*/

brp=dvector(2);
for (i=0; i<2; i++) brp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
rin=chao2(empdist,ntaxa);

/* increment true richness until that fails to improve likelihood */
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/

pbes = 0.0f;
/*	ei = (double) FITINC / 10;	*/
/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/

for (i=0; i<2; i++) brp[i] = -1.0*DBL_MAX;
for (i=0; i<3; i++) rs[i] = 0.0f;
lri[0]=lri[1]=0.0f;
ri=iri=rin/2;
while (ri+rin<= ntaxa)	ri/=2;
for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{
	obsrvd[0]=r-ntaxa;
	
	fitdist=dvector(r);
	for (i=0; i<r; ++i)	fitdist[i]=1/((double) r);
	expect=expfinds(fitdist,r,nspec,nspec);
	free_dvector(fitdist);

	rs[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
	free_dvector(expect);

	if (rs[0] >= brp[0]) {					/* IF BETTER THAN BEST FIT */
		pbes = brp[0];						/* save last best ssq for evenness */
		brp[0] = rs[0];						/* STORE FIT */
		brp[1] = r;							/* STORE SLOPE */
		
		/* while we are getting better on the initial increment, just ride with it			*/
		if (ri==iri)	{
			lri[1]=lri[0];
			lri[0]=ri;
			}
		else if ((r-ri)>=ntaxa)	{
			lri[1]=lri[0];
			lri[0]=ri;
			ri*=-1;
			}
		/* if we have a later improvement, wander halfway back to the last improvement 	*/
		/* (remember, we always start at the most likely slope up to that point			*/
		else	{
			lri[1]=lri[0];
			lri[0]=ri;
			ri=(ri/abs(ri))*(abs(ri)+1)/2;
			}
		}

	/* if this richness is better than the last		*/
	/* optimal richness is overshot 				*/
	else {
		r=brp[1];				/* step back to best richness	 */		
		/* if we start off going downhill, then cut the increment in half */
		if (ri==iri && r==rin)	{
			lri[1]=lri[0];
			lri[0]=ri;
			if (ri%2==0 || abs(ri)==1)	ri/=2;
			else						ri=(ri/abs(ri))*(abs(ri)+1)/2;
			}
		/* if this increment is opposite of the last, then cut it in  half */
		else if (abs(ri)==abs(lri[0]))	{
			lri[1]=lri[0];
			lri[0]=ri;
			/* go towards the one with the higher likelihood */
			if (abs(ri)==1)				ri=0;
			else						ri=(ri/abs(ri))*(abs(ri)+1)/2;
/*			if (rs[1]>rs[0])	{	
			if (
				if (ri%2==0 || abs(ri)==1)	ri/=-2;
				else						ri=(ri/abs(ri))*(abs(ri)+1)/-2;
				}
			else	{
				if (ri%2==0 || abs(ri)==1)	ri/=2;
				else						ri=(ri/abs(ri))*(abs(ri)+1)/2;
				}	*/
			}
		/* if this increment is half of the last, then reverse it - unless that goes too low! */
		else if (abs(ri)<abs(lri[0]))	{
			lri[1]=lri[0];
			lri[0]=ri;
			if ((r-ri)>=ntaxa)			ri*=-1;
			else						ri=(ri/abs(ri))*(abs(ri)+1)/2;	
			}
		else	{
			if (ri%2==0 || abs(ri)==1)	ri/=2;
			else						ri=(ri/abs(ri))*(abs(ri)+1)/2;
			}
		}
	
	/* do not let r exceed 5000 - that is the maximum richness that we can evaluate */
	/* do not bother going for r < observed! */
	while ((r+ri)< ntaxa)	ri=(ri/abs(ri))*(abs(ri)+1)/2;
	
	rs[2] = rs[1];
	rs[1] = brp[0];
	
	if (abs(lri[0])==1 && abs(lri[1])==1)	ri=0;
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
double *poi_fit_geo(int *empdist, int ntaxa, int nspec) 
{
int i = 0;					/* LOOP VARIABLE													*/
int r = 0;					/* LOOP RICHNESS													*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
double ev = 0.000f;			/* LOOP SLOPE 														*/
double emin = 1.000000001f;	/* min slope														*/
double ein = 0.000f;		/* initial slope													*/
double ei = 0.000f;			/* how much to increment ev in each loop							*/
double lei[2];				/* last two slope increments										*/
double iei;					/* initial slope increment at each richness							*/
double mei=0.000001f;		/* minimum evenness increment										*/
double es[3];				/* previous log likelihoods (cell number = num previous).			*/
double *bep;				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/
double *fitdist;			/* fit distribution													*/
double *expect;				/* expected number of species with 0Émax finds						*/
double *obsrvd;				/* observed number of species with 0Émax finds						*/
double pbes;				/* previous best support for decay rate								*/

bep=dvector(2);
for (i=0; i<2; i++) bep[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */
if (empdist[ntaxa-1]>=1)
	ein=pow(empdist[0],(((double) empdist[ntaxa-1])/((double) ntaxa)));
else
	ein=pow(empdist[0],1/((double) ntaxa));

if (ein<1.0)	ein=1.0f;
if (empdist[0]==empdist[ntaxa-1])	ein=1.0f;

/* increment true richness until that fails to improve likelihood */
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
iei=ein-1;
if (ei==0)	ei=0.00001f;

pbes = 0.0f;
/*	ei = (double) FITINC / 10;	*/
/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/

for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
lei[0]=lei[1]=0.0f;
ei=iei;
while (ei+ein<= emin)	ei/=2;
for (ev = ein; ((ev>=emin && fabs(ei)>mei) && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC)))); ev += ei) {	
	
	/* generate geometric distribution with richness r and decay of ev */
	fitdist = proportional_geo_distribution(1/ev,pow(10,-9));				/* MAKE DISTRIBUTION */

	for (r=1; ((1-(1/ev))*pow((1/ev),(r-1)))>pow(10,-9); r=r)	++r;

	/* find the expected proportions of taxa with 0Éx finds */
	expect=expfinds(fitdist,r,nspec,nspec);
	free_dvector(fitdist);
	es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
	free_dvector(expect);
	/*Debugging line */
	if (ev<=emin) printf("\nDANGER: Geo R=%d, ev=%f, S=%f ",r,ev,es[0]);

	if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
		pbes = bep[0];								/* save last best ssq for evenness */
		bep[0] = es[0];								/* STORE FIT */
		bep[1] = ev;								/* STORE SLOPE */
		
		/* while we are getting better on the initial increment, just ride with it			*/
		if (ei==iei)	{
			lei[1]=lei[0];
			lei[0]=ei;
			}
		/* if we have a later improvement, wander halfway back to the last improvement 	*/
		/* (remember, we always start at the most likely slope up to that point			*/
		else	{
			lei[1]=lei[0];
			lei[0]=ei;
			ei/=-2;
			if ((ev+ei)<=emin)	ei*=-1;
			}
		}
				
	/* If likelihood has not increased, then reset and change the increment  value	*/
	else	{
		/* if we went from x -> -x, then we want to cut the increment in half */
		if (ei==-1*lei[0] || ei==iei)	{
			lei[1]=lei[0];
			lei[0]=ei;
			if (ei==iei)	{
				/* determine whether es[0] or es[2] is the second best - move towards that	*/
				if (es[0]>es[2])	ei/=2;
				else				ei/=-2;
				}
			else	{
				/* determine whether es[0] or es[1] is the second best - move towards that	*/
				if (es[0]>es[1])	ei/=2;
				else				ei/=-2;
				}
			}
		/* if we just divided in half, then we want to reverse the increment */
		else if (ei==-1*iei || (2*fabs(ei)==fabs(lei[0])))	{
			lei[1]=lei[0];
			lei[0]=ei;
			/* if ev+x/2 got closer than ev-x, back up and go to ev+x/4	*/
			if (es[0]>es[2])	ei/=2;
			/* otherwise go to ev-x/2									*/
			else				ei*=-1;
			}
		ev=bep[1];
		}

	/* make sure that ei does not take ev below 1.0 */
	while (ev+ei<= emin)	{
		if (fabs(lei[0])==fabs(lei[1]))	ei*=-1;
		else							ei=-1*lei[0];
		}
	es[2] = es[1];									/* Store last 2 attempts to identify	*/
	es[1] = es[0];									/* when the peak is past 				*/

	}

free_dvector(obsrvd);
return bep;
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
double *poi_fit_geos(int *empdist, int ntaxa, int nspec) 
{
int i = 0;					/* LOOP VARIABLE													*/
int r = 0;					/* LOOP RICHNESS													*/
int ri;						/* richness increment												*/
int iri;					/* initial richness increment each loop								*/
int lri[2];					/* last two richness increments										*/
int rin;					/* initial richness													*/
int mxr=5000;				/* this is as long as we can make an array							*/
/*int lri[2];					/* last two richness increments									*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
double ev = 0.000f;			/* LOOP SLOPE 														*/
double emin = 1.000000001f;	/* min slope														*/
double ein = 0.000f;		/* initial slope													*/
double ei = 0.000f;			/* how much to increment ev in each loop							*/
double lei[2];				/* last two slope increments										*/
double iei;					/* initial slope increment at each richness							*/
double mei=0.000001f;		/* minimum evenness increment										*/
double es[3];				/* previous log likelihoods (cell number = num previous).			*/
double rs[3];				/* previous log likelihoods (cell number = num previous).			*/
double bep[2];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/
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
if (empdist[ntaxa-1]>=1)
	ein=pow(empdist[0],(((double) empdist[ntaxa-1])/((double) ntaxa)));
else
	ein=pow(empdist[0],1/((double) ntaxa));
if (ein<1.0)	ein=1.0f;
if (empdist[0]==empdist[ntaxa-1])	ein=1.0f;

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);
if (rin<ntaxa)	rin=ntaxa;
/*if ((rin-ntaxa)%2==1)	++rin;	*/
iri=ri=ntaxa/4;
if (ri<2)	iri=ri=2;

/* increment true richness until that fails to improve likelihood */
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
iei=ein-1;
if (ei==0)	ei=0.00001f;

for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;
/*	ei = (double) FITINC / 10;	*/
	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	
	for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
	lei[0]=lei[1]=0.0f;
	ei=iei;
	while (ei+ein<= emin)	ei/=2;
	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; ((ev>=emin && fabs(ei)>mei) && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC)))); ev += ei) {	
		
		/* generate geometric distribution with richness r and decay of ev */
		fitdist = proportional_geos_distribution(ev,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,nspec,nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, just ride with it			*/
			if (ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				if ((ev+ei)<=emin)	ei/=-2;
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lei[1]=lei[0];
				lei[0]=ei;
				ei/=-2;
				if ((ev+ei)<=emin)	ei/=-1;
				}
			}
					
		/* If likelihood has not increased, then reset and change the increment  value	*/
		else	{
			/* if we went from x -> -x, then we want to cut the increment in half */
			if (ei==-1*lei[0] || ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				/* if we are using the original incrementer */
				if (fabs(ei)==fabs(iei))	{
					/* determine whether es[0] or es[2] is the second best - move towards that	
					if (es[0]>es[2])	ei/=2;
					else				ei/=-2;*/
					if ((ev-ei)>emin)		ei/=-2;
					else					ei/=2;
					}
				else	{
					/* determine whether es[0] or es[1] is the second best - move towards that	*/
				/*	if (es[0]>es[1])	ei/=2;	*/
					if (fabs(lei[0])==fabs(lei[1]))	{
						ei/=-2;
						if ((ev+ei)<emin)	ei/=-1;
						}
					else	{
						ei/=-1;
						if ((ev+ei)<emin)	ei/=-2;
						}
					}	/* end case where ei is no the original */
				}
			/* if we just divided in half, then we want to reverse the increment */
			else if (ei==-1*iei || 2*fabs(ei)==fabs(lei[0]))	{
				lei[1]=lei[0];
				lei[0]=ei;
				/* if ev+x/2 got closer than ev-x, back up and go to ev+x/4	*
				if (es[0]>es[2])	ei/=2;								*/
				if (fabs(lei[0])==fabs(lei[1]))	ei/=2;
				/* otherwise go to ev-x/2									*/
				else							ei*=-1;
				}
			ev=bep[1];
			}

		/* make sure that ei does not take ev below 1.0 */
		while (ev+ei<= emin)	ei/=2;
		es[2] = es[1];									/* Store last 2 attempts to identify	*/
		es[1] = es[0];									/* when the peak is past 				*/

		}

	/* reset evenness incrementer */
	if (r!=rin)
		iei=bep[1]-brp[1];			/* set ei to the difference bn. current & former best */
	else if (r==rin || ei==0)
		iei=bep[1]-ein;
	
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
		/* reset evenness incrementer */
		if (r!=rin)	{
			iei=bep[1]-brp[1];			/* set ei to the difference bn. current & former best */
			if (fabs(iei)<0.01)			iei=0.01*(iei/fabs(iei));
			}
		else if (r==rin || ei==0)	{
			iei=bep[1]-ein;
			if (fabs(iei)<0.01)			iei=0.01*(iei/fabs(iei));
			}

		brp[2] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		ein=bep[1];										/* set initial decay to best decay found so far */ 

		lri[1]=lri[0];
		lri[0]=ri;

		/* if we are past the initial ri, then we want to use it only once	*/
		/* otherwise, we repeat r's											*/
		if (abs(ri)!=abs(iri))	{
			/* do not get stuck on 1!	*/
			if (iri!=1 && (abs(lri[0])==1 && abs(ri)==1))	ri=0;
			/* if proceeding ri/2 does not take you below ntaxa, then cut ri in half	*/
			else if ((r+((ri/abs(ri))*(abs(ri)+1)/2))>=ntaxa)
															ri=(ri/abs(ri))*(abs(ri)+1)/2;
			/* if proceeding ri/2 does takes you below ntaxa, then reverse ri	*/
			else											ri/=-1*(ri/abs(ri))*(abs(ri)+1)/2;
			}
		/* if on the absolute original ri, then go with it unless we are in danger of dropping below the minimum	*/
		else if ((r+ri)<ntaxa)								ri/=-1*(ri/abs(ri))*(abs(ri)+1)/2;
		}
	/* optimal richness is overshot */
	else {
		r=brp[2];				/* step back to best richness	 */		
		/* if we start off going downhill, then go the other way */
		if (ri==iri && r==rin)	{
			lri[1]=lri[0];
			lri[0]=ri;
			if (abs(lri[1])==1 && abs(lri[0])==1)	ri=0;
			else									ri=(ri/abs(ri))*(abs(ri)+1)/2;
			}
		/* if this increment is opposite of the last, then cut it in  half */
		else if (abs(ri)==abs(lri[0]))	{
			lri[1]=lri[0];
			lri[0]=ri;
			/* go towards the one with the higher likelihood */
			if (abs(lri[0])==1 && abs(lri[1])==1)	ri=0;
/*			else if (rs[1]>rs[0])						ri=(ri/abs(ri))*(abs(ri)+1)/-2;	*/
			else									ri=(ri/abs(ri))*(abs(ri)+1)/2;
			}
		/* if this increment is less then the last, then reverse it - unless that takes us too low	*/
		else if (abs(ri)<abs(lri[0]))	{
			lri[1]=lri[0];
			lri[0]=ri;
			if ((r-ri)>=ntaxa)						ri*=-1;
			else									ri=(ri/abs(ri))*(abs(ri)+1)/2;	
			}
		}
	
	/* do not let r exceed 5000 - that is the maximum richness that we can evaluate */
	/* do not bother going for r < observed! */
	while ((r+ri)>mxr || (r+ri)< ntaxa)	{
		if (r==mxr)	ri*=-1;
		else	{	
			if (abs(ri)==1)	ri=0;
			else			ri=(ri/abs(ri))*(abs(ri)+1)/2;
			}
		}
	
	/* geometric might go on and on forever with miniscule increases when it is a very poor fit */
	if (ri==iri && ((rs[0]>=rs[1] && rs[0]<(rs[1]+SUPINC)) && (rs[1]>=rs[2] && rs[1]<(rs[2]+SUPINC))))	ri=0;
		
	if (ei<0)	iei*=-1;
	if (ei<mei)	iei=100*mei;
	
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
*****************************************************************************************************/
double *poi_fit_zipf(int *empdist, int ntaxa, int nspec)
{
int i = 0;					/* LOOP VARIABLE														*/
int r = 0;					/* LOOP RICHNESS														*/
int ri;						/* richness increment													*/
int iri;					/* initial richness increment each loop									*/
int	lri[2];					/* recent changes in richness											*/
int mxr=10000;				/* this is as long as we can make an array								*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
int rin;					/* initial richness 													*/
double ev = 0.000f;			/* LOOP SLOPE 															*/
double emin = 1.0000001f;	/* min slope															*/
double ein = 1.00001f;		/* initial slope														*/
double lei[2];				/* last two evenness incrementers										*/
double iei;					/* initial slope increment at each richness								*/
double ei = 0.000f;			/* how much to increment ev in each loop								*/
double mei=0.000001f;		/* minimum evenness increment											*/
double es[3];				/* previous log likelihoods (cell number = num previous).				*/
double rs[3];				/* previous log likelihoods (cell number = num previous).				*/
double bep[2];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format		*/
double *brp;				/* BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array			*/
double *fitdist;			/* fit distribution														*/
double *expect;				/* expected number of species with 0Émax finds							*/
double *obsrvd;				/* observed number of species with 0Émax finds							*/
double pbes;				/* previous best support for decay rate									*/

brp=dvector(3);
for (i=0; i<3; i++) rs[i] = 0.0f;
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* find a good seed value for the Zipf-Mandelbrot, based on the log-log slopes between the initial taxa */
if (empdist[0]>empdist[2])	ein=(log(empdist[0])/log(empdist[2]));
if (empdist[0]==empdist[ntaxa-1])	ein=1.0f;
if (ntaxa>=10)	{
	ev=((double) empdist[9])/((double) empdist[0]);
	ein=-1*log(ev)/log(10);
	}
if (ein<emin)
/*	ein=pow((log(empdist[0])/log(empdist[6])),0.5);	*/	
	if (empdist[2]==1 && empdist[0]>empdist[1])
		ein=pow((log(empdist[0])/log(empdist[1])),2);

if (ein<emin)	ein=1.01f;

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);
iri=ri=ntaxa/4;
if (ri<2)	iri=ri=2;

/* increment true richness until that fails to improve likelihood */
/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
iei=ein-1;
if (iei==0)	iei=0.00001f;

for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;

	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
	lei[0]=lei[1]=0.0f;
	ei=iei;
	while (ei+ein<= emin)	ei/=2;

	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; ((ev>=emin && fabs(ei)>mei) && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC)))); ev += ei) {	
		
		/* generate zipf distribution with richness r and log-decay of ev */
		fitdist = proportional_zipf_distribution(ev,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,nspec,nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, just ride with it			*/
			if (ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				if ((ev+ei)<=emin)	ei/=-2;
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lei[1]=lei[0];
				lei[0]=ei;
				ei/=-2;
				if ((ev+ei)<=emin)	ei/=-1;
				}
			}	/* end case where likelihood is improved	*/
					
		/* If likelihood has not increased, then reset and change the increment  value	*/
		else	{
			/* if we went from x -> -x, then we want to cut the increment in half */
			if (ei==-1*lei[0] || ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				/* if we are using the original incrementer */
				if (fabs(ei)==fabs(iei))	{
					/* determine whether es[0] or es[2] is the second best - move towards that	
					if (es[0]>es[2])	ei/=2;
					else				ei/=-2;*/
					if ((ev-ei)>emin)		ei/=-2;
					else					ei/=2;
					}
				else	{
					/* determine whether es[0] or es[1] is the second best - move towards that	*/
				/*	if (es[0]>es[1])	ei/=2;	*/
					if (fabs(lei[0])==fabs(lei[1]))	{
						ei/=-2;
						if ((ev+ei)<emin)	ei/=-1;
						}
					else	{
						ei/=-1;
						if ((ev+ei)<emin)	ei/=-2;
						}
					}	/* end case where ei is no the original */
				}
			/* if we just divided in half, then we want to reverse the increment */
			else if (ei==-1*iei || 2*fabs(ei)==fabs(lei[0]))	{
				lei[1]=lei[0];
				lei[0]=ei;
				/* if ev+x/2 got closer than ev-x, back up and go to ev+x/4	*
				if (es[0]>es[2])	ei/=2;								*/
				if (fabs(lei[0])==fabs(lei[1]))	ei/=2;
				/* otherwise go to ev-x/2									*/
				else							ei*=-1;
				
				if (ev+ei<= emin)	{
					if (fabs(lei[0])==fabs(lei[1]))	ei/=-2;
					else							ei/=-1;
					}
				}
			ev=bep[1];
			}	/* end case where likelihood is not improved	*/

		/* make sure that ei does not take ev below 1.0 */
/*		if (ev+ei<= emin)	ei/=2;	*/
		es[2] = es[1];									/* Store last 2 attempts to identify	*/
		es[1] = es[0];									/* when the peak is past 				*/
		}

	/* reset evenness incrementer */
	if (r!=rin)
		iei=bep[1]-brp[1];			/* set ei to the difference bn. current & former best */
	else if (r==rin || ei==0)
		iei=bep[1]-ein;
	
	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
		/* reset evenness incrementer */
		if (r!=rin)	{
			iei=bep[1]-brp[1];			/* set ei to the difference bn. current & former best */
			if (fabs(iei)<0.01)			iei=0.01*(iei/fabs(iei));
			}
		else if (r==rin || ei==0)	{
			iei=bep[1]-ein;
			if (fabs(iei)<0.01)			iei=0.01*(iei/fabs(iei));
			}

		brp[2] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		ein=bep[1];										/* set initial decay to best decay found so far */ 

		lri[1]=lri[0];
		lri[0]=ri;

		/* if we are past the initial ri, then we want to use it only once	*/
		/* otherwise, we repeat r's											*/
		if (abs(ri)!=abs(iri))	{
			/* do not get stuck on 1!	*/
			if (iri!=1 && (abs(lri[0])==1 && abs(ri)==1))	ri=0;
			/* if proceeding ri/2 does not take you below ntaxa, then cut ri in half	*/
			else if ((r+((ri/abs(ri))*(abs(ri)+1)/2))>=ntaxa)
															ri=(ri/abs(ri))*(abs(ri)+1)/2;
			/* if proceeding ri/2 does takes you below ntaxa, then reverse ri	*/
			else											ri/=-1*(ri/abs(ri))*(abs(ri)+1)/2;
			}
		/* if on the absolute original ri, then go with it unless we are in danger of dropping below the minimum	*/
		else if ((r+ri)<ntaxa)								ri/=-1*(ri/abs(ri))*(abs(ri)+1)/2;
		}

	/* optimal richness is overshot */	
	else {
		r=brp[2];				/* step back to best richness	 */		
		/* if we start off going downhill, then go the other way */
		if (ri==iri && r==rin)	{
			lri[1]=lri[0];
			lri[0]=ri;
			if (abs(lri[1])==1 && abs(lri[0])==1)	ri=0;
			else									ri=(ri/abs(ri))*(abs(ri)+1)/2;
			}
		/* if this increment is opposite of the last, then cut it in  half */
		else if (abs(ri)==abs(lri[0]))	{
			lri[1]=lri[0];
			lri[0]=ri;
			/* go towards the one with the higher likelihood */
			if (abs(lri[0])==1 && abs(lri[1])==1)	ri=0;
/*			else if (rs[1]>rs[0])						ri=(ri/abs(ri))*(abs(ri)+1)/-2;	*/
			else									ri=(ri/abs(ri))*(abs(ri)+1)/2;
			}
		/* if this increment is less then the last, then reverse it - unless that takes us too low	*/
		else if (abs(ri)<abs(lri[0]))	{
			lri[1]=lri[0];
			lri[0]=ri;
			if ((r-ri)>=ntaxa)						ri*=-1;
			else									ri=(ri/abs(ri))*(abs(ri)+1)/2;	
			}
		}
	/* do not let r exceed 5000 - that is the maximum richness that we can evaluate */
	/* do not bother going for r < observed! 
	while ((r+ri)>mxr || (r+ri)< ntaxa)	{
		if (r==mxr)	ri*=-1;
		else	ri/=2;
		}				*/
	
	/* geometric will go on and on forever with miniscule increases when it is a very poor fit */
	if ((rs[0]>=rs[1] && rs[0]<(rs[1]+SUPINC)) && (rs[1]>=rs[2] && rs[1]<(rs[2]+SUPINC)))	ri=0;
	
/*	if (abs(ri)%2==1 && abs(ri)>1)	{
		if (ri>0)	++ri;
		else		--ri;
		}	*/

	if (ei<0)	iei*=-1;
	if (ei<mei)	iei=100*mei;
	
	rs[2] = rs[1];
	rs[1] = bep[0];
	}

free_dvector(obsrvd);
return brp;
}

/* FIND THE MOST-LIKELY UNTRUNCATED LOG-NORMAL SERIES, VARYING MAGNITUDE OF INCREASE AMONG OCTAVES. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: magnitude of increase per octave
	- result[2]: optimal richness
COMMENTS:
	- SUPINC, FITINC defined in header file
*****************************************************************************************************/
double *poi_fit_lnu(int *empdist, int ntaxa, int nspec) 
{
int i=0;					/* LOOP VARIABLE														*/
int r = 0;					/* LOOP RICHNESS														*/
int ri;						/* richness increment													*/
int iri;					/* initial richness increment each loop									*/
int lri[2];					/* previous richness increment											*/
int rin=0;					/* initial richness to use in each search (begins as ntaxa)				*/
int mxr=5000;				/* this is as long as we can make an array								*/
double lei[2];				/* last two evenness incrementers										*/
double iei;					/* initial slope increment at each richness								*/
double ev = 0.000f;			/* LOOP SLOPE 															*/
double emin = 1.00f;		/* min slope															*/
double mxev=75.0f;			/* maximum magnitude change to consider									*/
double ein = 0.000f;		/* initial slope														*/
double ei = 0.000f;			/* how much to increment ev in each loop								*/
double mei=0.000001f;		/* minimum evenness increment											*/
double es[3];				/* previous modal decay log likelihoods (cell number = num previous).	*/
double rs[3];				/* previous richness log likelihoods (cell number = num previous).		*/
double bep[2];				/* BEST decay for hypothesized evenness at given S - return array format */
double *brp;				/* BEST richness parameters - return array								*/
double pbes;				/* previous best support for modal decay								*/
double pbrs;				/* previous best support for richness									*/
double *expect;				/* expected number of species with 0Émax finds							*/
double *obsrvd;				/* observed number of species with 0Émax finds							*/
double *fitdist;			/* fit distribution												*/
double mdmx;

brp=dvector(3);
for (i=0; i<3; i++) brp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);
iri=ri=ntaxa/4;
if (ri<2)	ri=2;
lri[0]=lri[1]=0;
	

if (2*rin<10)			mdmx=1.5;
else if (2*rin<20)		mdmx=2.0;
else if (2*rin<75)		mdmx=2.5;
else if (2*rin<350)		mdmx=3.0;
else if (2*rin<2000)	mdmx=3.5;
else if (2*rin<15000)	mdmx=4.0;
else 					mdmx=4.5;

if (rin%2==1)				r=empdist[rin/2];
else						r=(((double) empdist[rin/2])+((double) empdist[(rin/2)-1]))/2;

if (rin/2 < ntaxa)			ein = pow(e,(log(((double) empdist[0])/((double) empdist[rin/2]))/mdmx));
else	{
	if (empdist[ntaxa-1]>0)	ein = pow(e,(log(((double) empdist[0])/((double) empdist[ntaxa-1]))/mdmx));
	else					ein = pow(e,log(((double) empdist[0])/mdmx));
	}
if (ein==1)	ein=1.20f;

iei=ein-1;

pbrs=0.0f;
for (i=0; i<3; i++) brp[i]=-1.0*DBL_MAX;
for (i=0; i<3; ++i)	rs[i]= -1.0*DBL_MAX;

/* adjust true richness until that fails to improve likelihood	*/
for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;

	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
	lei[0]=lei[1]=0.0f;
	ei=iei;
	/* generally speaking, if you decrease richness, then you increase the probability of observation	*/
	/* 		by increasing evenness, which means looking at lower values of ev.							*/
	if (r!=rin && r<brp[2])			if (ei>0)	ei*=-1;
	else if (r!=rin && r>brp[2])	if (ei<0)	ei*=-1;	

	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; ((ev>=emin && fabs(ei)>mei) && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC)))); ev += ei) {	
		
		/* generate log-normal distribution with richness r and decay of magnitude of increase ev */
		fitdist=proportional_lgn_distribution(ev,0,0,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,nspec,nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: ln R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
			pbes = bep[0];								/* save last best ssq for evenness */
			rs[0]=bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;								/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, or its negative, just ride with it			*/
			if (fabs(ei)==fabs(iei))	{
				lei[1]=lei[0];
				lei[0]=ei;
				if ((ev+ei)<=emin)	ei/=-2;				/* 6/2/05: if we are dropping towards the min, back up just to make sure we search space */ 
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lei[1]=lei[0];
				lei[0]=ei;
				ei/=-2;
				if ((ev+ei)<=emin)	ei*=-1;
				}
			}
					
		/* If likelihood has not increased, then reset and change the increment  value	*/
		else	{
			/* if we went from x -> -x, then we want to cut the increment in half */
			if (ei==-1*lei[0] || ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				/* if we are using the original incrementer */
				if (ei==iei)	{
					/* determine whether es[0] or es[2] is the second best - move towards that	
					if (es[0]>es[2])	ei/=2;
					else				ei/=-2;*/
					if ((ev-ei)>emin)		ei*=-1;
					else					ei/=2;
					}
				else	{
					/* determine whether es[0] or es[1] is the second best - move towards that	*/
				/*	if (es[0]>es[1])	ei/=2;	*/
					if (fabs(lei[0])==fabs(lei[1]))	{
						ei/=-2;
						if ((ev+ei)<emin)	ei/=-1;
						}
					else	{
						ei/=-1;
						if ((ev+ei)<emin)	ei/=-2;
						}
					}
				}
			/* if we just divided in half, then we want to reverse the increment */
			else if (ei==-1*iei || 2*fabs(ei)==fabs(lei[0]))	{
				lei[1]=lei[0];
				lei[0]=ei;
				/* if ev+x/2 got closer than ev-x, back up and go to ev+x/4	*
				if (es[0]>es[2])	ei/=2;								*/
				if (fabs(lei[0])==fabs(lei[1]))	ei/=2;
				/* otherwise go to ev-x/2									*/
				else							ei*=-1;
				}
			ev=bep[1];
			}

		/* make sure that ei does not take ev below 1.0 */
/*		while (ev+ei<= emin)	ei/=2;	*/
		if ((ev+ei)<=emin)	{
			if (fabs(ei)==fabs(lei[0]))			ei=lei[0]/2;	/* this will use half of the last positive ei 	*/
			else if (2*fabs(ei)==fabs(lei[0]))	ei*=-1;			/* this will happen if there is improvment		*/
			}
		es[2] = es[1];									/* Store last 2 attempts to identify	*/
		es[1] = es[0];									/* when the peak is past 				*/
		}

	/* if this richness is better than the last */
	if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
		/* reset evenness incrementer */
		if (r!=rin)	{
			iei=bep[1]-brp[1];			/* set ei to the difference bn. current & former best */
			if (fabs(iei)<0.01)			iei=0.01*(iei/fabs(iei));
			}
		else if (r==rin || ei==0)	{
			iei=bep[1]-ein;
			if (fabs(iei)<0.01)			iei=0.01*(iei/fabs(iei));
			}

		brp[2] = r;
		for (i=0; i<2; i++)
			brp[i] = bep[i];
		ein=bep[1];										/* set initial decay to best decay found so far */ 

		lri[1]=lri[0];
		lri[0]=ri;

		/* if we are past the initial ri, then we want to use it only once	*/
		/* otherwise, we repeat r's											*/
		if (abs(ri)!=abs(iri))	{
			/* do not get stuck on 1!	*/
			if (iri!=1 && (abs(lri[0])==1 && abs(ri)==1))	ri=0;
			/* if proceeding ri/2 does not take you below ntaxa, then cut ri in half	*/
			else if ((r+((ri/abs(ri))*(abs(ri)+1)/2))>=ntaxa)
															ri=(ri/abs(ri))*(abs(ri)+1)/2;
			/* if proceeding ri/2 does takes you below ntaxa, then reverse ri	*/
			else											ri/=-1*(ri/abs(ri))*(abs(ri)+1)/2;
			}
		/* if on the absolute original ri, then go with it unless we are in danger of dropping below the minimum	*/
		else if ((r+ri)<ntaxa)								ri/=-1*(ri/abs(ri))*(abs(ri)+1)/2;
		}	/* end cases of improved richness parameter */

	/* optimal richness is overshot */
	else {
		r=brp[2];				/* step back to best richness	 */		
		/* if we start off going downhill, then go the other way */
		if (ri==iri && r==rin)	{
			lri[1]=lri[0];
			lri[0]=ri;
			if (abs(lri[1])==1 && abs(lri[0])==1)	ri=0;
			else									ri=(ri/abs(ri))*(abs(ri)+1)/2;
			}
		/* if this increment is opposite of the last, then cut it in  half */
		else if (abs(ri)==abs(lri[0]))	{
			lri[1]=lri[0];
			lri[0]=ri;
			/* go towards the one with the higher likelihood */
			if (abs(lri[0])==1 && abs(lri[1])==1)	ri=0;
/*			else if (rs[1]>rs[0])						ri=(ri/abs(ri))*(abs(ri)+1)/-2;	*/
			else									ri=(ri/abs(ri))*(abs(ri)+1)/2;
			}
		/* if this increment is less then the last, then reverse it - unless that takes us too low	*/
		else if (abs(ri)<abs(lri[0]))	{
			lri[1]=lri[0];
			lri[0]=ri;
			if ((r-ri)>=ntaxa)						ri*=-1;
			else									ri=(ri/abs(ri))*(abs(ri)+1)/2;	
			}
		}
	
	/* do not let r exceed 5000 - that is the maximum richness that we can evaluate */
	/* do not bother going for r < observed! 
	while ((r+ri)>mxr || (r+ri)< ntaxa)	{
		if (r==mxr)	ri*=-1;
		else	ri/=2;
		}	*/
	
	/* geometric will go on and on forever with miniscule increases when it is a very poor fit */
	if ((rs[0]>=rs[1] && rs[0]<(rs[1]+SUPINC)) && (rs[1]>=rs[2] && rs[1]<(rs[2]+SUPINC)))	ri=0;
	
/*	if (abs(ri)%2==1 && abs(ri)>1)	{
		if (ri>0)	++ri;
		else		--ri;
		}	*/

	rs[2] = rs[1];
	rs[1] = bep[0];
	}
free_dvector(obsrvd);
return brp;
}


/* FIND THE MOST-LIKELY LOG-NORMAL SERIES, VARYING INTIAL AND MODAL SLOPE AND THE MODE. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood (truncated)
	- result[1]: magnitude of increase per octave (truncated)
	- result[2]: optimal richness (truncated)
	- result[3]: mode (truncated)
	- result[4]: log likelihood (untruncated)
	- result[5]: magnitude of increase per octave (untruncated)
	- result[6]: optimal richness (untruncated)	
	- result[7]: mode of untruncated is 0
COMMENTS:
	- SUPINC, FITINC defined in header file
*****************************************************************************************************/
double *poi_fit_lnt(int *empdist, int ntaxa, int nspec) 
{
int j=0,i=0;			/* LOOP VARIABLE															*/
int r = 0;				/* LOOP RICHNESS															*/
int ri;					/* richness increment														*/
int iri;				/* initial richness increment each loop										*/
int lri[2];				/* previous richness increment												*/
int rin=0;				/* initial richness to use in each search (begins as ntaxa)					*/
int mxr=5000;			/* this is as long as we can make an array									*/
int unqfnd=0;			/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
double mode;			/* mode of log-normal, given as a taxon number								*/
double x=1.0f;			/* x is the absolute value of the mode - cannot do abs(double)				*/
double y=1.0f;
double z=1.0f;
double mi=-1.0f;		/* how much to increment the mode											*/
double imi=-1.0f;
double lmi[2];			/* prior increment															*/
double mimin=0.125;		/* minimum mode change													*/
double mdmx=1.5f;
double mdmn=-2.5f;
double tr;
double ev = 0.000f;		/* LOOP SLOPE 																*/
double emin = 1.00f;	/* min slope																*/
double mxev=75.0f;		/* maximum magnitude change to consider										*/
double ein = 0.000f;	/* initial slope															*/
double ei = 0.000f;		/* how much to increment ev in each loop									*/
double iei;				/* initial slope increment at each richness								*/
double aei=0.00f;		/* absolute evenness increment											*/
double mei=0.000001f;	/* minimum evenness increment											*/
double es[3];			/* previous modal decay log likelihoods (cell number = num previous).		*/
double rs[3];			/* previous richness log likelihoods (cell number = num previous).			*/
double ms[3];			/* previous mode log likelihoods (cell number = num previous).				*/
double bep[2];			/* BEST modal decay parameters - return array format 						*/
double brp[3];			/* BEST richness parameters - return array									*/
double pbes;			/* previous best support for modal decay									*/
double pbrs;			/* previous best support for richness										*/
double pbms;			/* previous best support for mode location									*/
double *expect;			/* expected number of species with 0Émax finds								*/
double *obsrvd;			/* observed number of species with 0Émax finds								*/
double *fitdist;		/* fit distribution															*/
double *bmp;			/* Best mode parameters - return array format								*/
double lei[2];			/* last evenness increment													*/

bmp=dvector(8);
for (i=0; i<8; i++) bmp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);

if (rin<100)			j=4;
else					j=5;

if (2*rin<75)			mdmx=2.5;
else if (2*rin<350)		mdmx=3.0;
else if (2*rin<2000)	mdmx=3.5;
else if (2*rin<15000)	mdmx=4.0;
else 					mdmx=4.5;

if (rin/2 < ntaxa)	ein = pow(e,log(empdist[0]/empdist[rin/2])/mdmx);
else				ein = pow(e,log(empdist[0]/empdist[ntaxa-1])/mdmx);

aei=iei=ein-1;

mdmn=-1*mdmx;

/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/
for (i=0; i<3; ++i)	ms[i]=-1.0*DBL_MAX;
pbms=0.0f;
/* adjust mode until that fails to improve likelihood	*/
	/* for some reason the program is disobeying the second part of the conditional */
imi=mi;
for (mode=0; x>=mimin && (mode>=mdmn && mode<=mdmx)/* && ((pbms == 0.0f) || (bmp[0] > (pbms + SUPINC)))*/; mode+=mi)	{

	if (mode==0)	{
/*		ei = ein-0.95;	*/
		tr=0;
		}
	else	{
/*		ei = ein/10;	*/
		tr=1;
		}
	
	iei=-1*(ein-1.05);
	
	pbrs=0.0f;
	for (i=0; i<3; i++) brp[i]=-1.0*DBL_MAX;
	for (i=0; i<3; ++i)	rs[i]= -1.0*DBL_MAX;

	/* adjust true richness until that fails to improve likelihood	*/
	lri[0]=ri=rin;
	if (ri<2)	ri=2;
		
	iri=ri;
	for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{
		obsrvd[0]=r-ntaxa;
		pbes = 0.0f;

		/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
		for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
		for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
		lei[0]=lei[1]=0.0f;
		ei=iei;
		/* increment slope until that fails to improve likelihood or resolution limit reached */
		for (ev = ein; ((ev>=emin && aei>mei) && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC)))); ev += ei) {	
			
			/* generate log-normal distribution with richness r and decay of magnitude of increase ev */
			fitdist=proportional_lgn_distribution(ev,0,0,r);				/* MAKE DISTRIBUTION */
			/* find the expected proportions of taxa with 0Éx finds */
			expect=expfinds(fitdist,r,nspec,nspec);
			free_dvector(fitdist);
			es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
			free_dvector(expect);
			/*Debugging line */
			if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

			if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */
				pbes = bep[0];								/* save last best ssq for evenness */
				rs[0]=bep[0] = es[0];						/* STORE FIT */
				bep[1] = ev;								/* STORE SLOPE */
				
				/* while we are getting better on the initial increment, just ride with it			*/
				if (ei==iei)	{
					lei[1]=lei[0];
					lei[0]=ei;
					}
				/* if we have a later improvement, wander halfway back to the last improvement 	*/
				/* (remember, we always start at the most likely slope up to that point			*/
				else	{
					lei[1]=lei[0];
					lei[0]=ei;
					ei/=-2;
					}
				}
						
			/* If likelihood has not increased, then reset and change the increment  value	*/
			else	{
				/* if we went from x -> -x, then we want to cut the increment in half */
				if (ei==-1*lei[0] || ei==iei)	{
					lei[1]=lei[0];
					lei[0]=ei;
					if (ei==iei)	{
						/* determine whether es[0] or es[2] is the second best - move towards that	*/
						if (es[0]>es[2])	ei/=2;
						else				ei/=-2;
						}
					else	{
						/* determine whether es[0] or es[1] is the second best - move towards that	*/
						if (es[0]>es[1])	ei/=2;
						else				ei/=-2;
						}
					}
				/* if we just divided in half, then we want to reverse the increment */
				else if (ei==-1*iei || (2*ei==lei[0] || -2*ei==lei[0]))	{
					lei[1]=lei[0];
					lei[0]=ei;
					/* if ev+x/2 got closer than ev-x, back up and go to ev+x/4	*/
					if (es[0]>es[2])	ei/=2;
					/* otherwise go to ev-x/2									*/
					else				ei*=-1;
					}
				ev=bep[1];
				}

			/* make sure that ei does not take ev below 1.0 */
			while (ev+ei<= emin)	ei/=2;
			es[2] = es[1];									/* Store last 2 attempts to identify	*/
			es[1] = es[0];									/* when the peak is past 				*/

			aei=ei;											/* tally absolute evenness increment	*/
			if (aei<0)	aei*=-1;
			}

		/* reset evenness incrementer */
		if (r!=rin)
			aei=iei=bep[1]-brp[1];			/* set ei to the difference bn. current & former best */
		else if (r==rin || ei==0)
			aei=iei=bep[1]-ein;
		
		/* if this richness is better than the last */
		if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */
			brp[2] = r;
			for (i=0; i<2; i++)
				brp[i] = bep[i];
			ein=bep[1];										/* set initial decay to best decay found so far */ 
			/* if we are past the initial ri, then we want to use it only once	*/
			/* otherwise, we repeat r's											*/
			if (r!=rin)	{
				lri[1]=lri[0];
				lri[0]=ri;
				}
			if (ri!=iri && ri!=(-1*iri))		ri/=2;
			}
		/* optimal richness is overshot */
		else {
			r=brp[2];				/* step back to best richness	 */		
			/* if we start off going downhill, then go the other way */
			if (ri==iri && r==rin)	{
				lri[1]=lri[0];
				lri[0]=ri;
				ri*=-1;
				while ((r+ri)<ntaxa)		ri/=-2;
				}
			/* if this increment is opposite of the last, then cut it in half */
			else if (abs(ri)==abs(lri[0]))	{
				lri[1]=lri[0];
				lri[0]=ri;
				/* go towards the one with the higher likelihood */
				if (rs[1]>rs[0])			ri/=-2;
				else						ri/=2;
				}
			/* if this increment is less then the last, then */
			else if (abs(ri)<abs(lri[0]))	{
				lri[1]=lri[0];
				lri[0]=ri;
				ri*=-1;
				}
			}
		
		/* do not let r exceed 5000 - that is the maximum richness that we can evaluate */
		/* do not bother going for r < observed! */
		while ((r+ri)>mxr || (r+ri)< ntaxa)	{
			if (r==mxr)	ri*=-1;
			else	ri/=2;
			}
		
		if ((rs[0]>=rs[1] && rs[0]<(rs[1]+SUPINC)) && (rs[1]>=rs[2] && rs[1]<(rs[2]+SUPINC)))	ri=0;
		
		if (abs(ri)%2==1 && abs(ri)>1)	{
			if (ri>0)	++ri;
			else		--ri;
			}

		if (ei<0)	iei*=-1;
		if (ei<mei)	iei=100*mei;
		aei=iei;
		if (aei<0)	aei*=-1;
		
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
		/* if we are truncating the log-normal, then we want to use the mode shift only once */
		else {
			lmi[1]=lmi[0];
			lmi[0]=mi;
			}
		if (mi!=imi && mi!=(-1*imi))	mi/=2;
		}
	/* optimal mode is overshot */
	else	{
		mode=bmp[3];
		if (mi==imi && mode==0)	{
			lmi[1]=lmi[0];
			lmi[0]=mi;
			mi*=-1;
			}
		else if (mi==(-1*lmi[0]) || mi==lmi[0])	{
			lmi[1]=lmi[0];
			lmi[0]=mi;
			if (ms[1]>ms[0])	mi/=-2;
			else				mi/=2;
			} 
		else if (abs(mi)<abs(lmi[0]))	{
			lmi[1]=lmi[0];
			lmi[0]=mi;
			mi*=-1;
			}
		}
	/* it might go on and on forever with miniscule increases when it is a very poor fit */
	if ((ms[0]>=ms[1] && ms[0]<(ms[1]+SUPINC)) && (ms[1]>=ms[2] && ms[1]<(ms[2]+SUPINC)))	ri=0;
	ms[2] = ms[1];
	ms[1] = brp[0];

/*	if (msw>=2 && mrn==1)	{
		mi/=2;
		if (mi<1)	x=-1*mi;
		mrn=0;
		}	*/
	
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


/* FIND THE MOST-LIKELY ZERO-SUM MULTINOMIAL, VARYING m and theta. 
NEEDS:
	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest
	- ntaxa (size of array)
	- nspec (sum of array)
RETURNS:
	- result[0]: log likelihood
	- result[1]: most likely m
	- result[2]: most likely theta
COMMENTS:
	- SUPINC, FITINC defined in header file
*****************************************************************************************************/
double *poi_fit_zsm(int *empdist, int ntaxa, int nspec)
{
int i=0;					/* LOOP VARIABLE															*/
/*int	r=0;					/* richness (calculated from distribution									*/
double m=0.0f;				/* m for loop																*/
double mi=0.25f;			/* m increment																*/
double imi=0.25f;			/* initial m increment each loop											*/
double lmi[2];				/* previous m increment														*/
double min=0.75f;			/* initial m to use in each search (begins as ntaxa)						*/
long double mmx=0.999999f;	/* maximum m to consider													*/
double mmn=0.85f;			/* minimum m																*/
double lti[2];				/* last two theta incrementers												*/
double iti;					/* initial slope increment at each m										*/
double theta=0.000f;		/* LOOP theta 																*/
double tmin=0.01f;			/* min theta slope															*/
double tin=10.000f;			/* initial theta															*/
double ti=10.000f;			/* how much to increment theta in each loop									*/
double mti=0.0001f;			/* minimum theta increment													*/
double ts[3];				/* previous modal decay log likelihoods (cell number = num previous).		*/
double ms[3];				/* previous m log likelihoods (cell number = num previous).					*/
double *btp;				/* BEST decay for hypothesized theta at given S - return array format 		*/
double bmp[2];				/* BEST m parameters - return array											*/
double pbms;				/* previous best support for m												*/
double pbts=-1*DBL_MAX;		/* previous best support for theta											*/
double *expect;				/* expected number of species with 0Émax finds								*/
double *obsrvd;				/* observed number of species with 0Émax finds								*/
/*double *fitdist;			/* fit distribution															*/

btp=dvector(3);
for (i=0; i<3; i++) btp[i]=-1.0*DBL_MAX;
for (i=0; i<3; ++i)	ts[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

tin=chao2(empdist,ntaxa);
tin*=0.3;
ti=iti=tin/4;
mti=0.0001*iti;

lti[0]=tin;
lti[1]=0.0f;
	
min=(mmx+mmn)/2;		/* it seems that m's on either side of 0.85 repeat each other....	*/

/*imi=(1-min)/2;	*/
imi=mmx-min;
/* adjust true m until that fails to improve likelihood	*/
for (theta = tin; ((theta>=tmin && fabs(ti)>mti) || pbts == -1*DBL_MAX); theta += ti) {	

	/*Debugging line */
	if (theta<=tmin) printf("\nDANGER: theta=%f, m=%f ",theta, m);

	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	for (i=0; i<3; i++) ms[i] = -1.0*DBL_MAX;
	bmp[0]=bmp[1]= -1.0*DBL_MAX;
	lmi[0]=lmi[1]=0.0f;
	mi=imi;
	
	pbms=-1*DBL_MAX;

	/* find the best theta for this m */
	for (m=min; ((mi!=0 && (m>=mmn && m<=mmx))  && (pbms==-1*DBL_MAX || bmp[0]>pbms +SUPINC)); m+=mi)	{
		obsrvd[0]=m-ntaxa;

		expect=zerosum(theta,m,nspec,empdist[0]);		
		
		ms[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);

		/* if this m is better than the last */
		if (ms[0] >= bmp[0]) {							/* IF BETTER THAN BEST FIT */
			pbms = bmp[0];								/* save last best ssq for m */
			ts[0]=bmp[0]=ms[0];							/* STORE FIT */
			bmp[1]=m;									/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, just ride with it			*/
			if (fabs(mi)==fabs(imi))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				if (m+mi<mmn || m+mi>mmx)				mi/=-2;
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi/=2;
				if (m+mi<mmn || m+mi>mmx)				mi/=-1;
				}
			}

		/* optimal m is overshot */
		else {
			m=bmp[1];				/* step back to best m	 */		
			/* if we start off going downhill, then go the other way */
			
			if (mi==imi)						mi*=-1;
			else if (fabs(mi)==fabs(lmi[0]))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi/=2;
				if (m+mi>mmx || m+mi<mmn)	{
					if (fabs(mi)==fabs(lmi[0]))	mi/=-2;
					else						mi/=-1;
					}
				}
			else if (mi==lmi[0]/2)	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi*=-1;
				if (m+mi>mmx || m+mi<mmn)		mi/=-2;
				}
			else if (mi==imi && m==min)	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi*=-1;
				if (m+mi>mmx || m+mi<mmn)		mi/=-2;
				}
			/* if this increment is opposite of the last, then cut it in half */
			else if (fabs(mi)==fabs(lmi[0]))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				/* go towards the one with the higher likelihood */
				if (ms[1]>ms[0])				mi/=-2;
				else							mi/=2;
				if (m+mi>mmx || m+mi<mmn)		mi/=-1;
				}
			/* if this increment is less then the last, then */
			else if (fabs(mi)<fabs(lmi[0]))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi*=-1;
				if (m+mi>mmx || m+mi<mmn)		mi/=-2;
				}
			}
		
		if ((ms[0]>=ms[1] && ms[0]<(ms[1]+SUPINC)) && (ms[1]>=ms[2] && ms[1]<(ms[2]+SUPINC)))	
			mi=0;
		
		/* if we basically are running into improvement up to the min/max, then just take the min/max	*/
		else if ((bmp[1]==mmx || bmp[1]==mmn) && (bmp[0]<(ms[0]+0.1) && (bmp[0]<(ms[1]+0.1) && bmp[0]<(ms[2]+0.1))))
			mi=0;
			
		/* if it is trying to get past the maximum, then kill it	*/
		else if ((bmp[1]==mmx && (bmp[0]<btp[0] && btp[1]>mmx)) && (bmp[0]>ms[0] && (bmp[0]>ms[1] && bmp[0]>ms[2])))
			mi=0;

		/* if log likelihood is hopelessly low for this theta, then skip to the next	*/
		else if (bmp[0]<1.25*btp[0])
			mi=0;
		
		ms[2] = ms[1];
		ms[1] = ms[0];
		}	/* end search for m */


	if (ts[0] >= btp[0]) {							/* IF BETTER THAN BEST FIT */
		pbts = btp[0];								/* save last best ssq for theta */
		btp[0] = ts[0];								/* STORE FIT */
		btp[1] = bmp[1];							/* STORE m */
		btp[2] = theta;								/* STORE m */
		
		/* while we are getting better on the initial increment, just ride with it			*/
		if (fabs(ti)==fabs(iti))	{
			lti[1]=lti[0];
			lti[0]=ti;
			if (theta+ti<tmin)				ti/=-2;
			
			/* to avoid local optima, at least on intial search */
			else if (btp[0] <= (pbts + SUPINC))	{
				/* return to prior theta	*/
				if (theta==(tin+iti))	{
					theta=tin;
					ti*=-1;
					pbts-=1.0f;
					}
				else if (theta==(tin-iti))	{
					theta=tin;
					ti/=2;
					pbts-=1.0f;
					}
				}	/* end case where we had trivial improvement within one step of initial theta */
			}
		/* if we have a later improvement, wander halfway back to the last improvement 	*/
		/* (remember, we always start at the most likely slope up to that point			*/
		else	{
			lti[1]=lti[0];
			lti[0]=ti;
			ti/=2;
			if (theta+ti<tmin)				ti/=-1;
			}
		}	/* end case where likelihood has improved	*/
				
	/* If likelihood has not increased, then reset and change the increment  value	*/
	else	{
		/* if we went from x -> -x, then we want to cut the increment in half */
		lti[1]=lti[0];
		lti[0]=ti;
		if (fabs(ti)==fabs(iti))	{
			/* if the first increment does not help, then reverse direction	*/
			if (theta==tin+iti)				ti*=-1;
			/* if an increment after the first stops helping, then go half way	*/
			else							ti/=2;
			if (theta+ti<tmin)				ti/=-2;
			}
		else	{
			if (fabs(lti[0])==fabs(lti[1]))	ti/=2;
			else							ti*=-1;
			if (theta+ti<tmin)				ti/=-2;
			}
/*			if (ti==-1*lti[0] || ti==iti)	{
			if (ti==iti)	{
				/* determine whether ts[0] or ts[2] is the second best - move towards that	/
				if (ts[0]>ts[2])	ti/=2;
				else				ti/=-2;
				}
			else	{
				/* determine whether ts[0] or ts[1] is the second best - move towards that	/
				if (ts[0]>ts[1])	ti/=2;
				else				ti/=-2;
				}
			}
		/* if we just divided in half, then we want to reverse the increment /
		else if (ti==-1*iti || (2*ti==lti[0] || -2*ti==lti[0]))	{
			/* if theta+x/2 got closer than theta-x, back up and go to theta+x/4	/
			if (fabs(ti)==fabs(lti[0]))	ti/=2;
			else						ti*=-1;

			} */
		theta=btp[2];
		}
		
	/* start the next m at the best m so far */
/*	if (btp[1]>0.5)	imi=(mmx-btp[1])/2;	*/
/*	else			imi=(btp[1]-mmn)/2;	*/
/*	min=btp[1];		*/
	/* make sure that imi is not too small */
	if (fabs(imi)<0.0001)	{
		if (imi<0)						imi=-0.0001;
		else							imi=0.0001;
		}
	/* make sure that imi does not take m out of bounds */
	if (min+imi>mmx || min+imi<mmn)		imi*=-1;

	/* make sure that ti does not take theta below 1.0 */
	if (theta+ti<= tmin)	{
		if (fabs(lti[0])==fabs(ti))			ti/=-2;
		else								ti/=-1;
		}

	ts[2] = ts[1];									/* Store last 2 attempts to identify	*/
	ts[1] = ts[0];									/* when the peak is past 				*/

	
	if (btp[0] < (pbts + SUPINC) && fabs(ti)<((fabs(iti)/20)))	ti=0.0f;
	}	/* end search for theta */


tin=btp[2];
ti=iti=tin/4;
mti=0.0001*iti;

mmn=0.0001f;
mmx=0.85;

lti[0]=tin;
lti[1]=0.0f;
	
min=(mmx+mmn)/2;		/* it seems that m's on either side of 0.85 repeat each other....	*/
imi=mmx-min;

for (theta = tin; ((theta>=tmin && fabs(ti)>mti) || pbts == -1*DBL_MAX); theta += ti) {	

	/*Debugging line */
	if (theta<=tmin) printf("\nDANGER: theta=%f, m=%f ",theta, m);

	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	for (i=0; i<3; i++) ms[i] = -1.0*DBL_MAX;
	bmp[0]=bmp[1]= -1.0*DBL_MAX;
	lmi[0]=lmi[1]=0.0f;
	mi=imi;
	
	pbms=-1*DBL_MAX;

	/* find the best theta for this m */
	for (m=min; ((mi!=0 && (m>=mmn && m<=mmx))  && (pbms==-1*DBL_MAX || bmp[0]>pbms +SUPINC)); m+=mi)	{
		obsrvd[0]=m-ntaxa;

		expect=zerosum(theta,m,nspec,empdist[0]);		
		
		ms[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);

		/* if this m is better than the last */
		if (ms[0] >= bmp[0]) {							/* IF BETTER THAN BEST FIT */
			pbms = bmp[0];								/* save last best ssq for m */
			ts[0]=bmp[0]=ms[0];							/* STORE FIT */
			bmp[1]=m;									/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, just ride with it			*/
			if (fabs(mi)==fabs(imi))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				if (m+mi<mmn || m+mi>mmx)				mi/=-2;
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi/=2;
				if (m+mi<mmn || m+mi>mmx)				mi/=-1;
				}
			}

		/* optimal m is overshot */
		else {
			m=bmp[1];				/* step back to best m	 */		
			/* if we start off going downhill, then go the other way */
			
			if (mi==imi)						mi*=-1;
			else if (fabs(mi)==fabs(lmi[0]))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi/=2;
				if (m+mi>mmx || m+mi<mmn)	{
					if (fabs(mi)==fabs(lmi[0]))	mi/=-2;
					else						mi/=-1;
					}
				}
			else if (mi==lmi[0]/2)	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi*=-1;
				if (m+mi>mmx || m+mi<mmn)		mi/=-2;
				}
			else if (mi==imi && m==min)	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi*=-1;
				if (m+mi>mmx || m+mi<mmn)		mi/=-2;
				}
			/* if this increment is opposite of the last, then cut it in half */
			else if (fabs(mi)==fabs(lmi[0]))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				/* go towards the one with the higher likelihood */
				if (ms[1]>ms[0])				mi/=-2;
				else							mi/=2;
				if (m+mi>mmx || m+mi<mmn)		mi/=-1;
				}
			/* if this increment is less then the last, then */
			else if (fabs(mi)<fabs(lmi[0]))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi*=-1;
				if (m+mi>mmx || m+mi<mmn)		mi/=-2;
				}
			}
		
		if ((ms[0]>=ms[1] && ms[0]<(ms[1]+SUPINC)) && (ms[1]>=ms[2] && ms[1]<(ms[2]+SUPINC)))	
			mi=0;
		
		/* if we basically are running into improvement up to the min/max, then just take the min/max	*/
		else if ((bmp[1]==mmx || bmp[1]==mmn) && (bmp[0]<(ms[0]+0.1) && (bmp[0]<(ms[1]+0.1) && bmp[0]<(ms[2]+0.1))))
			mi=0;

		/* if log likelihood is hopelessly low for this theta, then skip to the next	*/
		else if (bmp[0]<1.75*btp[0])
			mi=0;
		else if ((fabs(mi)<fabs(imi)/16) && bmp[0]<1.25*btp[0])
			mi=0;
		
		ms[2] = ms[1];
		ms[1] = ms[0];
		}	/* end search for m */


	if (ts[0] >= btp[0]) {							/* IF BETTER THAN BEST FIT */
		pbts = btp[0];								/* save last best ssq for theta */
		btp[0] = ts[0];								/* STORE FIT */
		btp[1] = bmp[1];							/* STORE m */
		btp[2] = theta;								/* STORE m */
		
		/* while we are getting better on the initial increment, just ride with it			*/
		if (fabs(ti)==fabs(iti))	{
			lti[1]=lti[0];
			lti[0]=ti;
			if (theta+ti<tmin)				ti/=-2;
			
			/* to avoid local optima, at least on intial search */
			else if (btp[0] <= (pbts + SUPINC))	{
				/* return to prior theta	*/
				if (theta==(tin+iti))	{
					theta=tin;
					ti*=-1;
					pbts-=1.0f;
					}
				else if (theta==(tin-iti))	{
					theta=tin;
					ti/=2;
					pbts-=1.0f;
					}
				}	/* end case where we had trivial improvement within one step of initial theta */
			}
		/* if we have a later improvement, wander halfway back to the last improvement 	*/
		/* (remember, we always start at the most likely slope up to that point			*/
		else	{
			lti[1]=lti[0];
			lti[0]=ti;
			ti/=2;
			if (theta+ti<tmin)				ti/=-1;
			}
		}	/* end case where likelihood has improved	*/
				
	/* If likelihood has not increased, then reset and change the increment  value	*/
	else	{
		/* if we went from x -> -x, then we want to cut the increment in half */
		lti[1]=lti[0];
		lti[0]=ti;
		if (fabs(ti)==fabs(iti))	{
			/* if the first increment does not help, then reverse direction	*/
			if (theta==tin+iti)				ti*=-1;
			/* if an increment after the first stops helping, then go half way	*/
			else							ti/=2;
			if (theta+ti<tmin)				ti/=-2;
			}
		else	{
			if (fabs(lti[0])==fabs(lti[1]))	ti/=2;
			else							ti*=-1;
			if (theta+ti<tmin)				ti/=-2;
			}
/*			if (ti==-1*lti[0] || ti==iti)	{
			if (ti==iti)	{
				/* determine whether ts[0] or ts[2] is the second best - move towards that	/
				if (ts[0]>ts[2])	ti/=2;
				else				ti/=-2;
				}
			else	{
				/* determine whether ts[0] or ts[1] is the second best - move towards that	/
				if (ts[0]>ts[1])	ti/=2;
				else				ti/=-2;
				}
			}
		/* if we just divided in half, then we want to reverse the increment /
		else if (ti==-1*iti || (2*ti==lti[0] || -2*ti==lti[0]))	{
			/* if theta+x/2 got closer than theta-x, back up and go to theta+x/4	/
			if (fabs(ti)==fabs(lti[0]))	ti/=2;
			else						ti*=-1;

			} */
		theta=btp[2];
		}
		
	/* start the next m at the best m so far */
/*	if (btp[1]>0.5)	imi=(mmx-btp[1])/2;	*/
/*	else			imi=(btp[1]-mmn)/2;	*/
/*	min=btp[1];		*/
	/* make sure that imi is not too small */
	if (fabs(imi)<0.0001)	{
		if (imi<0)						imi=-0.0001;
		else							imi=0.0001;
		}
	/* make sure that imi does not take m out of bounds */
	if (min+imi>mmx || min+imi<mmn)		imi*=-1;

	/* make sure that ti does not take theta below 1.0 */
	if (theta+ti<= tmin)	{
		if (fabs(lti[0])==fabs(ti))			ti/=-2;
		else								ti/=-1;
		}

	ts[2] = ts[1];									/* Store last 2 attempts to identify	*/
	ts[1] = ts[0];									/* when the peak is past 				*/

	
	if (btp[0] < (pbts + SUPINC) && fabs(ti)<((fabs(iti)/20)))	ti=0.0f;
	}	/* end search for theta */

free_dvector(obsrvd);
return btp;
}

double *poi_fit_zsm_test(int *empdist, int ntaxa, int nspec)
{
int i=0;					/* LOOP VARIABLE															*/
/*int	r=0;					/* richness (calculated from distribution									*/
double m=0.0f;				/* m for loop																*/
double mi=0.25f;			/* m increment																*/
double imi=0.25f;			/* initial m increment each loop											*/
double lmi[2];				/* previous m increment														*/
double min=0.75f;			/* initial m to use in each search (begins as ntaxa)						*/
long double mmx=0.999999f;	/* maximum m to consider													*/
double mmn=0.85f;			/* minimum m																*/
double lti[2];				/* last two theta incrementers												*/
double iti;					/* initial slope increment at each m										*/
double theta=0.000f;		/* LOOP theta 																*/
double tmin=0.01f;			/* min theta slope															*/
double tin=10.000f;			/* initial theta															*/
double ti=10.000f;			/* how much to increment theta in each loop									*/
double mti=0.0001f;			/* minimum theta increment													*/
double ts[3];				/* previous modal decay log likelihoods (cell number = num previous).		*/
double ms[3];				/* previous m log likelihoods (cell number = num previous).					*/
double *btp;				/* BEST decay for hypothesized theta at given S - return array format 		*/
double bmp[2];				/* BEST m parameters - return array											*/
double pbms;				/* previous best support for m												*/
double pbts=-1*DBL_MAX;		/* previous best support for theta											*/
double *expect;				/* expected number of species with 0Émax finds								*/
double *obsrvd;				/* observed number of species with 0Émax finds								*/
/*double *fitdist;			/* fit distribution															*/
FILE	*fopen();	
FILE 	*outfile;

btp=dvector(3);
for (i=0; i<3; i++) btp[i]=-1.0*DBL_MAX;
for (i=0; i<3; ++i)	ts[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

tin=chao2(empdist,ntaxa);
tin*=0.3;
ti=iti=tin/4;
mti=0.0001*iti;

lti[0]=tin;
lti[1]=0.0f;
	
min=(mmx+mmn)/2;		/* it seems that m's on either side of 0.85 repeat each other....	*/

outfile=fopen("zerosumtest.txt","w");
fprintf(outfile,"theta\tÆtheta\tm\tÆm\tlnL\ttime\n");

/*imi=(1-min)/2;	*/
imi=mmx-min;
/* adjust true m until that fails to improve likelihood	*/
for (theta = tin; ((theta>=tmin && fabs(ti)>mti) || pbts == -1*DBL_MAX); theta += ti) {	

	/*Debugging line */
	if (theta<=tmin) printf("\nDANGER: theta=%f, m=%f ",theta, m);

	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	for (i=0; i<3; i++) ms[i] = -1.0*DBL_MAX;
	bmp[0]=bmp[1]= -1.0*DBL_MAX;
	lmi[0]=lmi[1]=0.0f;
	mi=imi;
	
	pbms=-1*DBL_MAX;

	/* find the best theta for this m */
	for (m=min; ((mi!=0 && (m>=mmn && m<=mmx))  && (pbms==-1*DBL_MAX || bmp[0]>pbms +SUPINC)); m+=mi)	{
		obsrvd[0]=m-ntaxa;

		expect=zerosum(theta,m,nspec,empdist[0]);		
		
		ms[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);

		i=time((time_t *)NULL);
		fprintf(outfile,"%7.6f\t%7.6f\t%7.6f\t%7.6f\t%7.6f\t%d\n",theta,ti,m,mi,ms[0],i);
/*		fflush(stdout);	*/
		fclose(outfile);
		outfile=fopen("zerosumtest.txt","a");

		/* if this m is better than the last */
		if (ms[0] >= bmp[0]) {							/* IF BETTER THAN BEST FIT */
			pbms = bmp[0];								/* save last best ssq for m */
			ts[0]=bmp[0]=ms[0];							/* STORE FIT */
			bmp[1]=m;									/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, just ride with it			*/
			if (fabs(mi)==fabs(imi))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				if (m+mi<mmn || m+mi>mmx)				mi/=-2;
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi/=2;
				if (m+mi<mmn || m+mi>mmx)				mi/=-1;
				}
			}

		/* optimal m is overshot */
		else {
			m=bmp[1];				/* step back to best m	 */		
			/* if we start off going downhill, then go the other way */
			
			if (mi==imi)						mi*=-1;
			else if (fabs(mi)==fabs(lmi[0]))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi/=2;
				if (m+mi>mmx || m+mi<mmn)	{
					if (fabs(mi)==fabs(lmi[0]))	mi/=-2;
					else						mi/=-1;
					}
				}
			else if (mi==lmi[0]/2)	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi*=-1;
				if (m+mi>mmx || m+mi<mmn)		mi/=-2;
				}
			else if (mi==imi && m==min)	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi*=-1;
				if (m+mi>mmx || m+mi<mmn)		mi/=-2;
				}
			/* if this increment is opposite of the last, then cut it in half */
			else if (fabs(mi)==fabs(lmi[0]))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				/* go towards the one with the higher likelihood */
				if (ms[1]>ms[0])				mi/=-2;
				else							mi/=2;
				if (m+mi>mmx || m+mi<mmn)		mi/=-1;
				}
			/* if this increment is less then the last, then */
			else if (fabs(mi)<fabs(lmi[0]))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi*=-1;
				if (m+mi>mmx || m+mi<mmn)		mi/=-2;
				}
			}
		
		if ((ms[0]>=ms[1] && ms[0]<(ms[1]+SUPINC)) && (ms[1]>=ms[2] && ms[1]<(ms[2]+SUPINC)))	
			mi=0;
		
		/* if we basically are running into improvement up to the min/max, then just take the min/max	*/
		else if ((bmp[1]==mmx || bmp[1]==mmn) && (bmp[0]<(ms[0]+0.1) && (bmp[0]<(ms[1]+0.1) && bmp[0]<(ms[2]+0.1))))
			mi=0;
			
		/* if it is trying to get past the maximum, then kill it	*/
		else if ((bmp[1]==mmx && (bmp[0]<btp[0] && btp[1]>mmx)) && (bmp[0]>ms[0] && (bmp[0]>ms[1] && bmp[0]>ms[2])))
			mi=0;

		/* if log likelihood is hopelessly low for this theta, then skip to the next	*/
		else if (bmp[0]<1.25*btp[0])
			mi=0;
		
		ms[2] = ms[1];
		ms[1] = ms[0];
		}	/* end search for m */


	if (ts[0] >= btp[0]) {							/* IF BETTER THAN BEST FIT */
		pbts = btp[0];								/* save last best ssq for theta */
		btp[0] = ts[0];								/* STORE FIT */
		btp[1] = bmp[1];							/* STORE m */
		btp[2] = theta;								/* STORE m */
		
		/* while we are getting better on the initial increment, just ride with it			*/
		if (fabs(ti)==fabs(iti))	{
			lti[1]=lti[0];
			lti[0]=ti;
			if (theta+ti<tmin)				ti/=-2;
			
			/* to avoid local optima, at least on intial search */
			else if (btp[0] <= (pbts + SUPINC))	{
				/* return to prior theta	*/
				if (theta==(tin+iti))	{
					theta=tin;
					ti*=-1;
					pbts-=1.0f;
					}
				else if (theta==(tin-iti))	{
					theta=tin;
					ti/=2;
					pbts-=1.0f;
					}
				}	/* end case where we had trivial improvement within one step of initial theta */
			}
		/* if we have a later improvement, wander halfway back to the last improvement 	*/
		/* (remember, we always start at the most likely slope up to that point			*/
		else	{
			lti[1]=lti[0];
			lti[0]=ti;
			ti/=2;
			if (theta+ti<tmin)				ti/=-1;
			}
		}	/* end case where likelihood has improved	*/
				
	/* If likelihood has not increased, then reset and change the increment  value	*/
	else	{
		/* if we went from x -> -x, then we want to cut the increment in half */
		lti[1]=lti[0];
		lti[0]=ti;
		if (fabs(ti)==fabs(iti))	{
			/* if the first increment does not help, then reverse direction	*/
			if (theta==tin+iti)				ti*=-1;
			/* if an increment after the first stops helping, then go half way	*/
			else							ti/=2;
			if (theta+ti<tmin)				ti/=-2;
			}
		else	{
			if (fabs(lti[0])==fabs(lti[1]))	ti/=2;
			else							ti*=-1;
			if (theta+ti<tmin)				ti/=-2;
			}
/*			if (ti==-1*lti[0] || ti==iti)	{
			if (ti==iti)	{
				/* determine whether ts[0] or ts[2] is the second best - move towards that	/
				if (ts[0]>ts[2])	ti/=2;
				else				ti/=-2;
				}
			else	{
				/* determine whether ts[0] or ts[1] is the second best - move towards that	/
				if (ts[0]>ts[1])	ti/=2;
				else				ti/=-2;
				}
			}
		/* if we just divided in half, then we want to reverse the increment /
		else if (ti==-1*iti || (2*ti==lti[0] || -2*ti==lti[0]))	{
			/* if theta+x/2 got closer than theta-x, back up and go to theta+x/4	/
			if (fabs(ti)==fabs(lti[0]))	ti/=2;
			else						ti*=-1;

			} */
		theta=btp[2];
		}
		
	/* start the next m at the best m so far */
/*	if (btp[1]>0.5)	imi=(mmx-btp[1])/2;	*/
/*	else			imi=(btp[1]-mmn)/2;	*/
/*	min=btp[1];							*/
	/* make sure that imi is not too small */
	if (fabs(imi)<0.0001)	{
		if (imi<0)						imi=-0.0001;
		else							imi=0.0001;
		}
	/* make sure that imi does not take m out of bounds */
	if (min+imi>mmx || min+imi<mmn)		imi*=-1;

	/* make sure that ti does not take theta below 1.0 */
	if (theta+ti<= tmin)	{
		if (fabs(lti[0])==fabs(ti))			ti/=-2;
		else								ti/=-1;
		}

	ts[2] = ts[1];									/* Store last 2 attempts to identify	*/
	ts[1] = ts[0];									/* when the peak is past 				*/

	
	if (btp[0] < (pbts + SUPINC) && fabs(ti)<((fabs(iti)/20)))	ti=0.0f;
	}	/* end search for theta */


tin=btp[2];
ti=iti=tin/4;
mti=0.0001*iti;

mmn=0.0001f;
mmx=0.85;

lti[0]=tin;
lti[1]=0.0f;
	
min=(mmx+mmn)/2;		/* it seems that m's on either side of 0.85 repeat each other....	*/
imi=mmx-min;

for (theta = tin; ((theta>=tmin && fabs(ti)>mti) || pbts == -1*DBL_MAX); theta += ti) {	

	/*Debugging line */
	if (theta<=tmin) printf("\nDANGER: theta=%f, m=%f ",theta, m);

	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	for (i=0; i<3; i++) ms[i] = -1.0*DBL_MAX;
	bmp[0]=bmp[1]= -1.0*DBL_MAX;
	lmi[0]=lmi[1]=0.0f;
	mi=imi;
	
	pbms=-1*DBL_MAX;

	/* find the best theta for this m */
	for (m=min; ((mi!=0 && (m>=mmn && m<=mmx))  && (pbms==-1*DBL_MAX || bmp[0]>pbms +SUPINC)); m+=mi)	{
		obsrvd[0]=m-ntaxa;

		expect=zerosum(theta,m,nspec,empdist[0]);		
		
		ms[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);

		i=time((time_t *)NULL);
		fprintf(outfile,"%7.6f\t%7.6f\t%7.6f\t%7.6f\t%7.6f\t%d\n",theta,ti,m,mi,ms[0],i);
/*		fflush(stdout);	*/
		fclose(outfile);
		outfile=fopen("zerosumtest.txt","a");

		/* if this m is better than the last */
		if (ms[0] >= bmp[0]) {							/* IF BETTER THAN BEST FIT */
			pbms = bmp[0];								/* save last best ssq for m */
			ts[0]=bmp[0]=ms[0];							/* STORE FIT */
			bmp[1]=m;									/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, just ride with it			*/
			if (fabs(mi)==fabs(imi))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				if (m+mi<mmn || m+mi>mmx)				mi/=-2;
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi/=2;
				if (m+mi<mmn || m+mi>mmx)				mi/=-1;
				}
			}

		/* optimal m is overshot */
		else {
			m=bmp[1];				/* step back to best m	 */		
			/* if we start off going downhill, then go the other way */
			
			if (mi==imi)						mi*=-1;
			else if (fabs(mi)==fabs(lmi[0]))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi/=2;
				if (m+mi>mmx || m+mi<mmn)	{
					if (fabs(mi)==fabs(lmi[0]))	mi/=-2;
					else						mi/=-1;
					}
				}
			else if (mi==lmi[0]/2)	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi*=-1;
				if (m+mi>mmx || m+mi<mmn)		mi/=-2;
				}
			else if (mi==imi && m==min)	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi*=-1;
				if (m+mi>mmx || m+mi<mmn)		mi/=-2;
				}
			/* if this increment is opposite of the last, then cut it in half */
			else if (fabs(mi)==fabs(lmi[0]))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				/* go towards the one with the higher likelihood */
				if (ms[1]>ms[0])				mi/=-2;
				else							mi/=2;
				if (m+mi>mmx || m+mi<mmn)		mi/=-1;
				}
			/* if this increment is less then the last, then */
			else if (fabs(mi)<fabs(lmi[0]))	{
				lmi[1]=lmi[0];
				lmi[0]=mi;
				mi*=-1;
				if (m+mi>mmx || m+mi<mmn)		mi/=-2;
				}
			}
		
		if ((ms[0]>=ms[1] && ms[0]<(ms[1]+SUPINC)) && (ms[1]>=ms[2] && ms[1]<(ms[2]+SUPINC)))	
			mi=0;
		
		/* if we basically are running into improvement up to the min/max, then just take the min/max	*/
		else if ((bmp[1]==mmx || bmp[1]==mmn) && (bmp[0]<(ms[0]+0.1) && (bmp[0]<(ms[1]+0.1) && bmp[0]<(ms[2]+0.1))))
			mi=0;

		/* if log likelihood is hopelessly low for this theta, then skip to the next	*/
		else if (bmp[0]<1.75*btp[0])
			mi=0;
		else if ((fabs(mi)<fabs(imi)/16) && bmp[0]<1.25*btp[0])
			mi=0;
		
		ms[2] = ms[1];
		ms[1] = ms[0];
		}	/* end search for m */


	if (ts[0] >= btp[0]) {							/* IF BETTER THAN BEST FIT */
		pbts = btp[0];								/* save last best ssq for theta */
		btp[0] = ts[0];								/* STORE FIT */
		btp[1] = bmp[1];							/* STORE m */
		btp[2] = theta;								/* STORE m */
		
		/* while we are getting better on the initial increment, just ride with it			*/
		if (fabs(ti)==fabs(iti))	{
			lti[1]=lti[0];
			lti[0]=ti;
			if (theta+ti<tmin)				ti/=-2;
			
			/* to avoid local optima, at least on intial search */
			else if (btp[0] <= (pbts + SUPINC))	{
				/* return to prior theta	*/
				if (theta==(tin+iti))	{
					theta=tin;
					ti*=-1;
					pbts-=1.0f;
					}
				else if (theta==(tin-iti))	{
					theta=tin;
					ti/=2;
					pbts-=1.0f;
					}
				}	/* end case where we had trivial improvement within one step of initial theta */
			}
		/* if we have a later improvement, wander halfway back to the last improvement 	*/
		/* (remember, we always start at the most likely slope up to that point			*/
		else	{
			lti[1]=lti[0];
			lti[0]=ti;
			ti/=2;
			if (theta+ti<tmin)				ti/=-1;
			}
		}	/* end case where likelihood has improved	*/
				
	/* If likelihood has not increased, then reset and change the increment  value	*/
	else	{
		/* if we went from x -> -x, then we want to cut the increment in half */
		lti[1]=lti[0];
		lti[0]=ti;
		if (fabs(ti)==fabs(iti))	{
			/* if the first increment does not help, then reverse direction	*/
			if (theta==tin+iti)				ti*=-1;
			/* if an increment after the first stops helping, then go half way	*/
			else							ti/=2;
			if (theta+ti<tmin)				ti/=-2;
			}
		else	{
			if (fabs(lti[0])==fabs(lti[1]))	ti/=2;
			else							ti*=-1;
			if (theta+ti<tmin)				ti/=-2;
			}
/*			if (ti==-1*lti[0] || ti==iti)	{
			if (ti==iti)	{
				/* determine whether ts[0] or ts[2] is the second best - move towards that	/
				if (ts[0]>ts[2])	ti/=2;
				else				ti/=-2;
				}
			else	{
				/* determine whether ts[0] or ts[1] is the second best - move towards that	/
				if (ts[0]>ts[1])	ti/=2;
				else				ti/=-2;
				}
			}
		/* if we just divided in half, then we want to reverse the increment /
		else if (ti==-1*iti || (2*ti==lti[0] || -2*ti==lti[0]))	{
			/* if theta+x/2 got closer than theta-x, back up and go to theta+x/4	/
			if (fabs(ti)==fabs(lti[0]))	ti/=2;
			else						ti*=-1;

			} */
		theta=btp[2];
		}
		
	/* start the next m at the best m so far */
/*	if (btp[1]>0.5)	imi=(mmx-btp[1])/2;	*/
/*	else			imi=(btp[1]-mmn)/2;	*/
/*	min=btp[1];		*/
	/* make sure that imi is not too small */
	if (fabs(imi)<0.0001)	{
		if (imi<0)						imi=-0.0001;
		else							imi=0.0001;
		}
	/* make sure that imi does not take m out of bounds */
	if (min+imi>mmx || min+imi<mmn)		imi*=-1;

	/* make sure that ti does not take theta below 1.0 */
	if (theta+ti<= tmin)	{
		if (fabs(lti[0])==fabs(ti))			ti/=-2;
		else								ti/=-1;
		}

	ts[2] = ts[1];									/* Store last 2 attempts to identify	*/
	ts[1] = ts[0];									/* when the peak is past 				*/

	
	if (btp[0] < (pbts + SUPINC) && fabs(ti)<((fabs(iti)/20)))	ti=0.0f;
	}	/* end search for theta */

free_dvector(obsrvd);
fclose(outfile);
return btp;
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
double bep[2];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format		*/
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
	
	for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; (ev>=emin && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC)))); ev += ei) {	
		
		/* generate geometric distribution with richness r and decay of ev */
		fitdist = proportional_geos_distribution(ev,r);				/* MAKE DISTRIBUTION */
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
double bep[2];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/
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
	
	for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; (ev>=emin && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC)))); ev += ei) {	
		
		/*Debugging line */
		if (ev<=emin)	printf("\nDANGER: ZM R=%d, ev=%f, S=%f ",r,ev,es[0]);

		/* generate Zipf-Mandelbrot distribution with richnes r and decay of ev */
		fitdist = proportional_zipf_distribution(ev,r);				/* MAKE DISTRIBUTION */
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
	- result[0]: log likelihood (truncated)
	- result[1]: magnitude of increase per octave (truncated)
	- result[2]: optimal richness (truncated)
	- result[3]: mode (truncated)
	- result[4]: log likelihood (untruncated)
	- result[5]: magnitude of increase per octave (untruncated)
	- result[6]: optimal richness (untruncated)	
	- result[7]: mode of untruncated is 0
COMMENTS:
	- SUPINC, FITINC defined in header file
*****************************************************************************************************/
double *wht_fit_lnt(int *empdist, int ntaxa, int nspec) 
{
int j=0,i=0;				/* LOOP VARIABLE														*/
int r = 0;					/* LOOP RICHNESS														*/
int ri;						/* richness increment													*/
int iri;						/* initial richness increment each loop								*/
int lri;					/* previous richness increment											*/
int rin=0;					/* initial richness to use in each search (begins as ntaxa)				*/
int mxr=5000;				/* this is as long as we can make an array								*/
double mode;				/* mode of log-normal, given as a taxon number							*/
double x=1.0f;				/* x is the absolute value of the mode - cannot do abs(double)			*/
double y=1.0f;
double z=1.0f;
double mi=-1.0f;			/* how much to increment the mode										*/
double lmi=2.0f;			/* prior increment														*/
double mimin=0.125;			/* minimum mode change													*/
double mdmx=1.5f;
double mdmn=-2.5f;
double tr;
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
double *obsrvd;				/* observed number of species with 0Émax finds							*/
double *fitdist;			/* fit distribution														*/
double *bmp;				/* Best mode parameters - return array format							*/

bmp=dvector(8);
for (i=0; i<8; i++) bmp[i]= -1.0*DBL_MAX;

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);

/* use Chao 2 estimator to get seed richness */
rin=chao2(empdist,ntaxa);

if (rin<100)			j=4;
else					j=5;
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

	if (mode==0)	{
		ei = ein-1;
		tr=0;
		}
	else	{
		ei = ein/10;
		tr=1;
		}
	
	pbrs=0.0f;
	for (i=0; i<3; i++) brp[i]=-1.0*DBL_MAX;
	for (i=0; i<3; ++i)	rs[i]= -1.0*DBL_MAX;

	/* adjust true richness until that fails to improve likelihood	*/
	lri=ri=rin;
	if (ri<2)	ri=2;
		
	iri=ri;
	for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{

		/* make sure that mode starts between beginning and end */
		while (abs(mode/2)>r)	++r;
		
		obsrvd[0]=r-ntaxa;
		pbes = 0.0f;
		/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
		
		for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
		for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;

		if (ein+ei>=mxev)	ei*=-1;
		/* increment slope until that fails to improve likelihood or resolution limit reached */
		for (ev=ein; ev<=mxev && ((pbes==0.0f) || (bep[0]>(pbes+SUPINC))); ev+=ei)	{
			
			/* generate Log-Normal distribution with richness r with each octave increasing ev times */
/*			fitdist = proportional_lgn_distributionN(ev,mode,ocs,r);		/* MAKE DISTRIBUTION */
			fitdist=proportional_lgn_distribution(ev,tr,mode,r);
			es[0] = calc_likelihood_Foote(fitdist, empdist, nspec);	/* CALCULATE SUPPORTt */
			free_dvector(fitdist);

			if (es[0] >= bep[0]) {					/* IF BETTER THAN BEST FIT */
				pbes = bep[0];						/* save last best likelihood for evenness */
				rs[0] = bep[0] = es[0];				/* STORE FIT */
				bep[1] = ev;						/* STORE SLOPE */
				}
			/* check to make sure that this is correct elsewhere */
			else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>=emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */
				ev = bep[1];						/* step back to best magnitude */
				ei *= -1;							/* STEP BACKWARD TO FIND PEAK */
				es[1] = es[2];						/* SET PREVIOUS S TO IGNORE OVER STEP */
				}
			else {									/* NOT IMPROVING TRY SMALLER INCREMENT */
				ev = bep[1];						/* step back to best magnitude */
				ei /= 2;							/* SET SMALLER UNIT */
				}
			/* it might go on and on forever with miniscule increases when it is a very poor fit */
			if ((es[0]>=es[1] && es[0]<(es[1]+SUPINC)) && (es[1]>=es[2] && es[1]<(es[2]+SUPINC)))	break;
			es[2] = es[1];							/* Store last 2 attempts to identify */
			es[1] = es[0];							/* when the peak is past */
			}

		/* reset evenness incrementer */
		if (r!=rin)
			ei=bep[1]-brp[1];			/* set ei to the difference bn. current & former best */
		else if (r==rin)
			ei=bep[1]-ein;
		
		if (ei==0)	ei=bep[1]/10;

		/* if this richness is better than the last */
		if (bep[0]>=brp[0]) {						/* IF BETTER THAN BEST FIT */
			pbrs=brp[0];
			brp[2] = r;
			for (i=0; i<2; i++)
				brp[i] = bep[i];
			ms[0]=brp[0];
			ein=bep[1];								/* set initial decay to best decay found so far */ 
			
			/* if we are past the initial ri, then we want to use it only once	*/
			/* otherwise, we repeat r's											*/
			if (ri!=iri && ri!=(-1*iri))	ri/=2;
			lri=ri;
			}
		/* optimal richness is overshot - back up */
		else if (abs(ri)==abs(lri))	{
			r=brp[2];				/* step back to best richness	 */
			lri=ri;
			ri/=2;
			}
		/* two tries in this direction produced poorer results - try the other direction */
		else if (abs(lri)>abs(ri))	{
			r=brp[2];				/* step back to best richness	 */
			lri=ri;
			/* do not look at values lower than observed richness */
			if ((r-ri)>=ntaxa)	ri*=-1;
			else				ri/=2;
			}
		/* optimal richness is overshot */
		else if (rs[2]>rs[1] && rs[0]>rs[1])	{	
			r=brp[2];				/* step back to best richness	 */
			ri*=-1;				/* step backwards to peak */
			while ((r+ri)<ntaxa && abs(ri)>0)	ri/=2;
			rs[1]=rs[2];		/* set to prior ln L to ignore over step */
			}
		else	{
			r=brp[2];				/* step back to best richness	 */
			if (abs(ri)>1)	ri/=2;
			else			ri=0;
			while ((r+ri)<=ntaxa && abs(ri)>0)	ri/=2;
			}
		/* it might go on and on forever with miniscule increases when it is a very poor fit */
		if ((ri==iri || ri==(-1*iri)) && ((rs[0]>=rs[1] && rs[0]<(rs[1]+SUPINC)) && (rs[1]>=rs[2] && rs[1]<(rs[2]+SUPINC))))	ri=0;

		/* do not let r exceed 5000 - that is the maximum richness that we can evaluate */
		/* do not bother going for r < observed! */
		while ((r+ri)>mxr)	{
			if (r==mxr)	{
				ri*=-1;
				ri/=2;
				}
			else	ri=mxr-r;
			}
		
		while ((r+ri)< ntaxa)	{
			if (r==ntaxa)	{
				ri*=-1;
				ri/=2;
				}
			else	ri=ntaxa-r;
			}

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
		/* once we drop to half a mode, then we do not want to go two steps in the same direction */
		if (lmi>-1 && lmi<1)	{
			lmi=mi;
			mi/=2;
			}
		else lmi=mi;
		}
	else	{
		y=mi;
		if (y<0)	y*=-1;
		z=lmi;
		if (z<0)	z*=-1;
		
		if (mode==mdmx)			mdmx-=mimin;
		else if (mode==mdmn)	mdmn+=mimin;
		
		mode-=mi;
		lmi=mi;

		if (y==z)	/* this means we'll go twice in the same direction - that is bad... */
			mi/=2;
			
		else if (z>y)		/* if we've stepped down, then we want to step opposite */
			mi*=-1;

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
		}

	/* it might go on and on forever with miniscule increases when it is a very poor fit */
	if ((ms[0]>=ms[1] && ms[0]<(ms[1]+SUPINC)) && (ms[1]>=ms[2] && ms[1]<(ms[2]+SUPINC)))	ri=0;
	ms[2] = ms[1];
	ms[1] = brp[0];

/*	if (msw>=2 && mrn==1)	{
		mi/=2;
		if (mi<1)	x=-1*mi;
		mrn=0;
		}	*/
	
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
double bep[2];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/
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
	
	for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))); ev += ei) {	
		
		/* generate geometric distribution with richnes r and decay of ev */
		fitdist = proportional_geos_distribution(ev,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,nspec,nspec);
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
	if (expected[n]==0)	expected[n]=pow(10,-323);	/* this is about as low as you can go */
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
double *expfindsln (double mag, int S, int mxfds, int t)
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


/*support_bars_geos - routine to put error bars on richness given a Geometric distribution
*****************************************************************************************/
double **support_bars_unif (double bestS, int ntaxa, int bestR, int *empdist, double bar)
{
int i = 0;					/* LOOP VARIABLE														*/
int r = 0;					/* LOOP RICHNESS														*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
int rlb;					/* lower bound for richness												*/
int rub;					/* upper bound for richness												*/
int mtaxa=1000;				/* maximum number of taxa to consider									*/
int nspec;					/* number of specimens													*/
double diff;				/* difference in log-likelihoods										*/
double *fitdist;			/* fit distribution														*/
double *expect;				/* expected number of species with 0Émax finds							*/
double *obsrvd;				/* observed number of species with 0Émax finds							*/
double pbes;				/* previous best support for decay rate									*/
double supportlo[1000];		/* support for richness values below the ML value						*/
double supportup[1000];		/* support for richness values above the ML value						*/
double **sbars;				/* support and slopes for richness values within bar units of support	*/

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);
nspec=sumivector(empdist,ntaxa);

supportlo[0]=bestS;
diff=0;

rub=rlb=bestR;
for (r=bestR-1; diff<=bar && r>=ntaxa; --r)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;

	fitdist=dvector(r);
	for (i=0; i<r; ++i)	fitdist[i]=1/((double) r);

	expect=expfinds(fitdist,r,nspec,nspec);
	free_dvector(fitdist);
	supportlo[abs(bestR-r)-1] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
	free_dvector(expect);
	
	diff=fabs(bestS-supportlo[abs(bestR-r)-1]);
	rlb=r;
	}

diff=0;
for (r=bestR+1; diff<=bar && r<=mtaxa; ++r)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;

	fitdist=dvector(r);
	for (i=0; i<r; ++i)	fitdist[i]=1/((double) r);

	expect=expfinds(fitdist,r,nspec,nspec);
	free_dvector(fitdist);
	supportup[abs(bestR-r)-1] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
	free_dvector(expect);

	diff=fabs(bestS-supportup[abs(bestR-r)-1]);
	rub=r;
	}

free_dvector(obsrvd);

sbars=dmatrix(2+(rub-rlb),3);
for (r=rlb; r<=rub; ++r)	{
	sbars[r-rlb][0]=r;
	if (r<bestR)	{
		sbars[r-rlb][1]=1.0;
		sbars[r-rlb][2]=supportlo[abs(bestR-r)-1];
		}
	else if (r>bestR)	{
		sbars[r-rlb][1]=1.0;
		sbars[r-rlb][2]=supportup[abs(bestR-r)-1];
		}
	else	{
		sbars[r-rlb][1]=1.0;
		sbars[r-rlb][2]=bestS;
		}
	}

sbars[(rub-rlb)+1][0]=sbars[(rub-rlb)+1][1]=sbars[(rub-rlb)+1][2]=100;
return sbars;
}

/*support_bars_geos - routine to put error bars on richness given a Geometric distribution
*****************************************************************************************/
double **support_bars_geos (double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar)
{
int i = 0;					/* LOOP VARIABLE														*/
int r = 0;					/* LOOP RICHNESS														*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
int rlb;					/* lower bound for richness												*/
int rub;					/* upper bound for richness												*/
int mtaxa=1000;				/* maximum number of taxa to consider									*/
int nspec;					/* number of specimens													*/
double ev = 0.000f;			/* LOOP SLOPE 															*/
double emin = 1.0000001f;	/* min slope															*/
double ein = 1.00001f;		/* initial slope														*/
double lei[2];				/* last two evenness incrementers										*/
double iei;					/* initial slope increment at each richness								*/
double ei = 0.000f;			/* how much to increment ev in each loop								*/
double mei=0.000001f;		/* minimum evenness increment											*/
double es[3];				/* previous log likelihoods (cell number = num previous).				*/
double bep[2];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format		*/
double *fitdist;			/* fit distribution														*/
double *expect;				/* expected number of species with 0Émax finds							*/
double *obsrvd;				/* observed number of species with 0Émax finds							*/
double pbes;				/* previous best support for decay rate									*/
double supportlo[1000];		/* support for richness values below the ML value						*/
double supportup[1000];		/* support for richness values above the ML value						*/
double evlo[1000];			/* evenness maximizing likelihood of richness values below ML richness	*/
double evup[1000];			/* evenness maximizing likelihood of richness values above ML richness	*/
double **sbars;				/* support and slopes for richness values within bar units of support	*/

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);
nspec=sumivector(empdist,ntaxa);

ein=bestM;
iei=-1*(1/((double) nspec));
if (iei==0)	iei=-0.00001f;

bep[0]=bestS;
for (r=bestR-1; fabs(bestS-bep[0])<=bar && r>=ntaxa; --r)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;

	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
	lei[0]=lei[1]=0.0f;
	ei=iei;
	while (ei+ein<= emin)	ei/=2;

	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; ((ev>=emin && fabs(ei)>mei) && ((pbes == 0.0f) || (bep[0] > (pbes + (SUPINC/100))))); ev += ei) {	
		
		/* generate Geometric distribution with richness r and log-decay of ev */
		fitdist = proportional_geos_distribution(ev,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,nspec,nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {					/* IF BETTER THAN BEST FIT */
			pbes = bep[0];						/* save last best ssq for evenness */
			bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;						/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, just ride with it			*/
			if (ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lei[1]=lei[0];
				lei[0]=ei;
				ei/=-2;
				}
			}
					
		/* If likelihood has not increased, then reset and change the increment  value	*/
		else	{
			/* if we went from x -> -x, then we want to cut the increment in half */
			if (ei==-1*lei[0] || ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				if (ei==iei)	{
					/* determine whether es[0] or es[2] is the second best - move towards that	*/
					if (es[0]>es[2])	ei/=2;
					else				ei/=-2;
					}
				else	{
					/* determine whether es[0] or es[1] is the second best - move towards that	*/
					if (es[0]>es[1])	ei/=2;
					else				ei/=-2;
					}
				}
			/* if we just divided in half, then we want to reverse the increment */
			else if (ei==-1*iei || (2*ei==lei[0] || -2*ei==lei[0]))	{
				lei[1]=lei[0];
				lei[0]=ei;
				/* if ev+x/2 got closer than ev-x, back up and go to ev+x/4	*/
				if (es[0]>es[2])	ei/=2;
				/* otherwise go to ev-x/2									*/
				else				ei*=-1;
				}
			ev=bep[1];
			}

		/* make sure that ei does not take ev below 1.0 */
		while (ev+ei<= emin)	ei/=2;
		es[2] = es[1];									/* Store last 2 attempts to identify	*/
		es[1] = es[0];									/* when the peak is past 				*/
		}

	supportlo[abs(bestR-r)-1]=bep[0];
	evlo[abs(bestR-r)-1]=bep[1];

	rlb=r;
	}

ein=bestM;
iei=(2/((double) nspec));
if (iei==0)	iei=0.00001f;

bep[0]=bestS;
for (r=bestR+1; fabs(bestS-bep[0])<=bar && r<=mtaxa; ++r)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;

	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
	lei[0]=lei[1]=0.0f;
	ei=iei;
	while (ei+ein<= emin)	ei/=2;

	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; ((ev>=emin && fabs(ei)>mei) && ((pbes == 0.0f) || (bep[0] > (pbes + (SUPINC/100))))); ev += ei) {	
		
		/* generate zipf distribution with richness r and log-decay of ev */
		fitdist = proportional_geos_distribution(ev,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,nspec,nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {					/* IF BETTER THAN BEST FIT */
			pbes = bep[0];						/* save last best ssq for evenness */
			bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;						/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, just ride with it			*/
			if (ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lei[1]=lei[0];
				lei[0]=ei;
				ei/=-2;
				}
			}
					
		/* If likelihood has not increased, then reset and change the increment  value	*/
		else	{
			/* if we went from x -> -x, then we want to cut the increment in half */
			if (ei==-1*lei[0] || ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				if (ei==iei)	{
					/* determine whether es[0] or es[2] is the second best - move towards that	*/
					if (es[0]>es[2])	ei/=2;
					else				ei/=-2;
					}
				else	{
					/* determine whether es[0] or es[1] is the second best - move towards that	*/
					if (es[0]>es[1])	ei/=2;
					else				ei/=-2;
					}
				}
			/* if we just divided in half, then we want to reverse the increment */
			else if (ei==-1*iei || (2*ei==lei[0] || -2*ei==lei[0]))	{
				lei[1]=lei[0];
				lei[0]=ei;
				/* if ev+x/2 got closer than ev-x, back up and go to ev+x/4	*/
				if (es[0]>es[2])	ei/=2;
				/* otherwise go to ev-x/2									*/
				else				ei*=-1;
				}
			ev=bep[1];
			}

		/* make sure that ei does not take ev below 1.0 */
		while (ev+ei<= emin)	ei/=2;
		es[2] = es[1];									/* Store last 2 attempts to identify	*/
		es[1] = es[0];									/* when the peak is past 				*/
		}

	supportup[abs(bestR-r)-1]=bep[0];
	evup[abs(bestR-r)-1]=bep[1];
	rub=r;

	}

free_dvector(obsrvd);

sbars=dmatrix(2+(rub-rlb),3);
for (r=rlb; r<=rub; ++r)	{
	sbars[r-rlb][0]=r;
	if (r<bestR)	{
		sbars[r-rlb][1]=evlo[abs(bestR-r)-1];
		sbars[r-rlb][2]=supportlo[abs(bestR-r)-1];
		}
	else if (r>bestR)	{
		sbars[r-rlb][1]=evup[abs(bestR-r)-1];
		sbars[r-rlb][2]=supportup[abs(bestR-r)-1];
		}
	else	{
		sbars[r-rlb][1]=bestM;
		sbars[r-rlb][2]=bestS;
		}
	}

sbars[(rub-rlb)+1][0]=sbars[(rub-rlb)+1][1]=sbars[(rub-rlb)+1][2]=100;
return sbars;
}

/*support_bars_zipf - routine to put error bars on richness given a Zipf distribution
*****************************************************************************************/
double **support_bars_zipf (double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar)
{
int i = 0;					/* LOOP VARIABLE														*/
int r = 0;					/* LOOP RICHNESS														*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
int rlb;					/* lower bound for richness												*/
int rub;					/* upper bound for richness												*/
int mtaxa=1000;				/* maximum number of taxa to consider									*/
int nspec;					/* number of specimens													*/
double ev = 0.000f;			/* LOOP SLOPE 															*/
double emin = 1.0000001f;	/* min slope															*/
double ein = 1.00001f;		/* initial slope														*/
double lei[2];				/* last two evenness incrementers										*/
double iei;					/* initial slope increment at each richness								*/
double ei = 0.000f;			/* how much to increment ev in each loop								*/
double mei=0.000001f;		/* minimum evenness increment											*/
double es[3];				/* previous log likelihoods (cell number = num previous).				*/
double bep[2];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format		*/
double *fitdist;			/* fit distribution														*/
double *expect;				/* expected number of species with 0Émax finds							*/
double *obsrvd;				/* observed number of species with 0Émax finds							*/
double pbes;				/* previous best support for decay rate									*/
double supportlo[1000];		/* support for richness values below the ML value						*/
double supportup[1000];		/* support for richness values above the ML value						*/
double evlo[1000];			/* evenness maximizing likelihood of richness values below ML richness	*/
double evup[1000];			/* evenness maximizing likelihood of richness values above ML richness	*/
double **sbars;				/* support and slopes for richness values within bar units of support	*/

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);
nspec=sumivector(empdist,ntaxa);

ein=bestM;
iei=-1*(2/((double) nspec));
if (iei==0)	iei=-0.00001f;

bep[0]=bestS;
for (r=bestR-1; fabs(bestS-bep[0])<=bar && r>=ntaxa; --r)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;

	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
	lei[0]=lei[1]=0.0f;
	ei=iei;
	while (ei+ein<= emin)	ei/=2;

	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; ((ev>=emin && fabs(ei)>mei) && ((pbes == 0.0f) || (bep[0] > (pbes + (SUPINC/100))))); ev += ei) {	
		
		/* generate zipf distribution with richness r and log-decay of ev */
		fitdist = proportional_zipf_distribution(ev,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,nspec,nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {					/* IF BETTER THAN BEST FIT */
			pbes = bep[0];						/* save last best ssq for evenness */
			bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;						/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, just ride with it			*/
			if (ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lei[1]=lei[0];
				lei[0]=ei;
				ei/=-2;
				}
			}
					
		/* If likelihood has not increased, then reset and change the increment  value	*/
		else	{
			/* if we went from x -> -x, then we want to cut the increment in half */
			if (ei==-1*lei[0] || ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				if (ei==iei)	{
					/* determine whether es[0] or es[2] is the second best - move towards that	*/
					if (es[0]>es[2])	ei/=2;
					else				ei/=-2;
					}
				else	{
					/* determine whether es[0] or es[1] is the second best - move towards that	*/
					if (es[0]>es[1])	ei/=2;
					else				ei/=-2;
					}
				}
			/* if we just divided in half, then we want to reverse the increment */
			else if (ei==-1*iei || (2*ei==lei[0] || -2*ei==lei[0]))	{
				lei[1]=lei[0];
				lei[0]=ei;
				/* if ev+x/2 got closer than ev-x, back up and go to ev+x/4	*/
				if (es[0]>es[2])	ei/=2;
				/* otherwise go to ev-x/2									*/
				else				ei*=-1;
				}
			ev=bep[1];
			}

		/* make sure that ei does not take ev below 1.0 */
		while (ev+ei<= emin)	ei/=2;
		es[2] = es[1];									/* Store last 2 attempts to identify	*/
		es[1] = es[0];									/* when the peak is past 				*/
		}

	supportlo[abs(bestR-r)-1]=bep[0];
	evlo[abs(bestR-r)-1]=bep[1];

	rlb=r;
	}

ein=bestM;
iei=(2/((double) nspec));
if (iei==0)	iei=0.00001f;

bep[0]=bestS;
for (r=bestR+1; fabs(bestS-bep[0])<=bar && r<=mtaxa; ++r)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;

	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
	lei[0]=lei[1]=0.0f;
	ei=iei;
	while (ei+ein<= emin)	ei/=2;

	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; ((ev>=emin && fabs(ei)>mei) && ((pbes == 0.0f) || (bep[0] > (pbes + (SUPINC/100))))); ev += ei) {	
		
		/* generate zipf distribution with richness r and log-decay of ev */
		fitdist = proportional_zipf_distribution(ev,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,nspec,nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {					/* IF BETTER THAN BEST FIT */
			pbes = bep[0];						/* save last best ssq for evenness */
			bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;						/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, just ride with it			*/
			if (ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lei[1]=lei[0];
				lei[0]=ei;
				ei/=-2;
				}
			}
					
		/* If likelihood has not increased, then reset and change the increment  value	*/
		else	{
			/* if we went from x -> -x, then we want to cut the increment in half */
			if (ei==-1*lei[0] || ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				if (ei==iei)	{
					/* determine whether es[0] or es[2] is the second best - move towards that	*/
					if (es[0]>es[2])	ei/=2;
					else				ei/=-2;
					}
				else	{
					/* determine whether es[0] or es[1] is the second best - move towards that	*/
					if (es[0]>es[1])	ei/=2;
					else				ei/=-2;
					}
				}
			/* if we just divided in half, then we want to reverse the increment */
			else if (ei==-1*iei || (2*ei==lei[0] || -2*ei==lei[0]))	{
				lei[1]=lei[0];
				lei[0]=ei;
				/* if ev+x/2 got closer than ev-x, back up and go to ev+x/4	*/
				if (es[0]>es[2])	ei/=2;
				/* otherwise go to ev-x/2									*/
				else				ei*=-1;
				}
			ev=bep[1];
			}

		/* make sure that ei does not take ev below 1.0 */
		while (ev+ei<= emin)	ei/=2;
		es[2] = es[1];									/* Store last 2 attempts to identify	*/
		es[1] = es[0];									/* when the peak is past 				*/
		}

	supportup[abs(bestR-r)-1]=bep[0];
	evup[abs(bestR-r)-1]=bep[1];
	rub=r;

	}

free_dvector(obsrvd);

sbars=dmatrix(2+(rub-rlb),3);
for (r=rlb; r<=rub; ++r)	{
	sbars[r-rlb][0]=r;
	if (r<bestR)	{
		sbars[r-rlb][1]=evlo[abs(bestR-r)-1];
		sbars[r-rlb][2]=supportlo[abs(bestR-r)-1];
		}
	else if (r>bestR)	{
		sbars[r-rlb][1]=evup[abs(bestR-r)-1];
		sbars[r-rlb][2]=supportup[abs(bestR-r)-1];
		}
	else	{
		sbars[r-rlb][1]=bestM;
		sbars[r-rlb][2]=bestS;
		}
	}

sbars[(rub-rlb)+1][0]=sbars[(rub-rlb)+1][1]=sbars[(rub-rlb)+1][2]=100;
return sbars;
}


double **support_bars_lnu (double bestS, int ntaxa, int bestR, double bestM, int *empdist, double bar)
{
int i = 0;					/* LOOP VARIABLE														*/
int r = 0;					/* LOOP RICHNESS														*/
int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/
int rlb;					/* lower bound for richness												*/
int rub;					/* upper bound for richness												*/
int mtaxa=1000;				/* maximum number of taxa to consider									*/
int nspec;					/* number of specimens													*/
double ev = 0.000f;			/* LOOP SLOPE 															*/
double emin = 1.0000001f;	/* min slope															*/
double ein = 1.00001f;		/* initial slope														*/
double lei[2];				/* last two evenness incrementers										*/
double iei;					/* initial slope increment at each richness								*/
double ei = 0.000f;			/* how much to increment ev in each loop								*/
double mei=0.000001f;		/* minimum evenness increment											*/
double es[3];				/* previous log likelihoods (cell number = num previous).				*/
double bep[2];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format		*/
double *fitdist;			/* fit distribution														*/
double *expect;				/* expected number of species with 0Émax finds							*/
double *obsrvd;				/* observed number of species with 0Émax finds							*/
double pbes;				/* previous best support for decay rate									*/
double supportlo[1000];		/* support for richness values below the ML value						*/
double supportup[1000];		/* support for richness values above the ML value						*/
double evlo[1000];			/* evenness maximizing likelihood of richness values below ML richness	*/
double evup[1000];			/* evenness maximizing likelihood of richness values above ML richness	*/
double **sbars;				/* support and slopes for richness values within bar units of support	*/

/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */
obsrvd=idhistogram(empdist,ntaxa);
nspec=sumivector(empdist,ntaxa);

ein=bestM;
iei=-1*(1/((double) nspec));
if (iei==0)	iei=-0.00001f;

bep[0]=bestS;
for (r=bestR-1; fabs(bestS-bep[0])<=bar && r>=ntaxa; --r)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;

	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
	lei[0]=lei[1]=0.0f;
	ei=iei;
	while (ei+ein<= emin)	ei/=2;

	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; ((ev>=emin && fabs(ei)>mei) && ((pbes == 0.0f) || (bep[0] > (pbes + (SUPINC/100))))); ev += ei) {	
		
		/* generate zipf distribution with richness r and log-decay of ev */
		fitdist = proportional_lgn_distribution(ev,0,0,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,nspec,nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {					/* IF BETTER THAN BEST FIT */
			pbes = bep[0];						/* save last best ssq for evenness */
			bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;						/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, just ride with it			*/
			if (ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lei[1]=lei[0];
				lei[0]=ei;
				ei/=-2;
				}
			}
					
		/* If likelihood has not increased, then reset and change the increment  value	*/
		else	{
			/* if we went from x -> -x, then we want to cut the increment in half */
			if (ei==-1*lei[0] || ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				if (ei==iei)	{
					/* determine whether es[0] or es[2] is the second best - move towards that	*/
					if (es[0]>es[2])	ei/=2;
					else				ei/=-2;
					}
				else	{
					/* determine whether es[0] or es[1] is the second best - move towards that	*/
					if (es[0]>es[1])	ei/=2;
					else				ei/=-2;
					}
				}
			/* if we just divided in half, then we want to reverse the increment */
			else if (ei==-1*iei || (2*ei==lei[0] || -2*ei==lei[0]))	{
				lei[1]=lei[0];
				lei[0]=ei;
				/* if ev+x/2 got closer than ev-x, back up and go to ev+x/4	
				if (es[0]>es[2])	ei/=2;*/
				if (fabs(lei[0])==fabs(lei[1]))	ei/=2;
				/* otherwise go to ev-x/2				*/								
				else				ei*=-1;	
				}
			ev=bep[1];
			}

		/* make sure that ei does not take ev below 1.0 */
		while (ev+ei<= emin)	ei/=2;
		es[2] = es[1];									/* Store last 2 attempts to identify	*/
		es[1] = es[0];									/* when the peak is past 				*/
		}

	supportlo[abs(bestR-r)-1]=bep[0];
	evlo[abs(bestR-r)-1]=bep[1];

	rlb=r;
	}

ein=bestM;
iei=(2/((double) nspec));
if (iei==0)	iei=0.00001f;

bep[0]=bestS;
for (r=bestR+1; fabs(bestS-bep[0])<=bar && r<=mtaxa; ++r)	{
	obsrvd[0]=r-ntaxa;
	pbes = 0.0f;

	/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/
	for (i=0; i<2; i++) bep[i] = -1.0*DBL_MAX;
	for (i=0; i<3; i++) es[i] = -1.0*DBL_MAX;
	lei[0]=lei[1]=0.0f;
	ei=iei;
	while (ei+ein<= emin)	ei/=2;

	/* increment slope until that fails to improve likelihood or resolution limit reached */
	for (ev = ein; ((ev>=emin && fabs(ei)>mei) && ((pbes == 0.0f) || (bep[0] > (pbes + (SUPINC/100))))); ev += ei) {	
		
		/* generate zipf distribution with richness r and log-decay of ev */
		fitdist = proportional_lgn_distribution(ev,0,0,r);				/* MAKE DISTRIBUTION */
		/* find the expected proportions of taxa with 0Éx finds */
		expect=expfinds(fitdist,r,nspec,nspec);
		free_dvector(fitdist);
		es[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);
		free_dvector(expect);
		/*Debugging line */
		if (ev<=emin) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);

		if (es[0] >= bep[0]) {					/* IF BETTER THAN BEST FIT */
			pbes = bep[0];						/* save last best ssq for evenness */
			bep[0] = es[0];						/* STORE FIT */
			bep[1] = ev;						/* STORE SLOPE */
			
			/* while we are getting better on the initial increment, just ride with it			*/
			if (ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				}
			/* if we have a later improvement, wander halfway back to the last improvement 	*/
			/* (remember, we always start at the most likely slope up to that point			*/
			else	{
				lei[1]=lei[0];
				lei[0]=ei;
				ei/=-2;
				}
			}
					
		/* If likelihood has not increased, then reset and change the increment  value	*/
		else	{
			/* if we went from x -> -x, then we want to cut the increment in half */
			if (ei==-1*lei[0] || ei==iei)	{
				lei[1]=lei[0];
				lei[0]=ei;
				if (ei==iei)	{
					/* determine whether es[0] or es[2] is the second best - move towards that	*/
					if (es[0]>es[2])	ei/=2;
					else				ei/=-2;
					}
				else	{
					/* determine whether es[0] or es[1] is the second best - move towards that	*/
					if (es[0]>es[1])	ei/=2;
					else				ei/=-2;
					}
				}
			/* if we just divided in half, then we want to reverse the increment */
			else if (ei==-1*iei || (2*ei==lei[0] || -2*ei==lei[0]))	{
				lei[1]=lei[0];
				lei[0]=ei;
				/* if ev+x/2 got closer than ev-x, back up and go to ev+x/4	*/
				if (es[0]>es[2])	ei/=2;
				/* otherwise go to ev-x/2									*/
				else				ei*=-1;
				}
			ev=bep[1];
			}

		/* make sure that ei does not take ev below 1.0 */
		while (ev+ei<= emin)	ei/=2;
		es[2] = es[1];									/* Store last 2 attempts to identify	*/
		es[1] = es[0];									/* when the peak is past 				*/
		}

	supportup[abs(bestR-r)-1]=bep[0];
	evup[abs(bestR-r)-1]=bep[1];
	rub=r;

	}

free_dvector(obsrvd);

sbars=dmatrix(2+(rub-rlb),3);
for (r=rlb; r<=rub; ++r)	{
	sbars[r-rlb][0]=r;
	if (r<bestR)	{
		sbars[r-rlb][1]=evlo[abs(bestR-r)-1];
		sbars[r-rlb][2]=supportlo[abs(bestR-r)-1];
		}
	else if (r>bestR)	{
		sbars[r-rlb][1]=evup[abs(bestR-r)-1];
		sbars[r-rlb][2]=supportup[abs(bestR-r)-1];
		}
	else	{
		sbars[r-rlb][1]=bestM;
		sbars[r-rlb][2]=bestS;
		}
	}
sbars[(rub-rlb)+1][0]=sbars[(rub-rlb)+1][1]=sbars[(rub-rlb)+1][2]=100;

return sbars;
}