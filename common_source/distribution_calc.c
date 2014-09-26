#define distribution_calc
#include "distribution_calc.h"
#include "matrixanalysis.h"
#include "memory.h"
#include "more_math.h"
#include "probability.h"
#include "sort.h"

#define FUNC(a,b,c,d,e,f)	((*func)(a,b,c,d,e,f))	/* used in midpnt - added 2004.12.27	*/
#define EPS			1.0e-6         					/* used in qromo - added 2004.12.27		*/
#define JMAX		14								/* used in qromo - added 2004.12.27		*/
#define JMAXP		(JMAX+1)       					/* used in qromo - added 2004.12.27		*/
#define K			5         						/* used in qromo - added 2004.12.27		*/

/*
CALCULATES GEOMETRIC DISTRIBUTION. This uses a linear decay between R and ln relative abundance where R
	is the rank of the taxon.  The area under this curve from 1ÉS then is used to standardize the height
	so that it represents the proportions of taxon ranked 1ÉS.  
NEEDS:
    - M (slope in log-linear space)
    - min (the minimum abundance value that you wish to consider)
RETURNS:
	- abundance: proportional abundances
***********************************************************************/
double *proportional_geo_distribution(double M, long double min)
{
int i = 0;							/* loop variable */
long	S=1;
double	x;
double *A;							/* array to return */

if (M>1.0f) {
	printf("\nproportional_geo_distribution, illegal slope = %f\n",M);
	exit(1);
	}

x=1+((log(min)-log(1-M))/log(M));
S=x;
if ((x-S)>=0.5)	++S;
/*for (S=1; ((1-M)*pow(M,(S-1)))>min; S=S)	++S;	*/

A = dvector(S);					/* allocate array 									*/

A[0] = 1-M;						/* a[0] = 100; everything else is M^i times less 	*/
for (i=1; i<S; i++)
	A[i] = A[i-1] * M;			/* taxon i = taxon (i-1) / slope 					*/
								/* 		a probability function						*/

return A;
}

/*
CALCULATES GEOMETRIC DISTRIBUTION. This uses a linear decay between R and ln relative abundance where R
	is the rank of the taxon.  The area under this curve from 1ÉS then is used to standardize the height
	so that it represents the proportions of taxon ranked 1ÉS.  
NEEDS:
    - M (slope in log-linear space)
    - S (number of taxa)
RETURNS:
	- abundance: proportional abundances
***********************************************************************/
double *proportional_geos_distribution(double M, int S)
{
int i = 0;							/* loop variable */
double sum = 0.0f;					/* number of "occurrences" */
double *A;							/* array to return */

if (M<1.0f) {
	printf("\nproportional_geos_distribution, illegal slope = %f\n",M);
	exit(1);
	}
if (S<=0)	{
	printf("\nproportional_geos_distribution, illegal number of taxa = %d\n",S);
	exit(1);
	}

A = dvector(S);						/* allocate array */

sum = A[0] = 100;					/* a[0] = 100; everything else is M^i times less */
for (i=1; i<S; i++) {
	A[i] = A[i-1] / M;				/* taxon i = taxon (i-1) / slope */
	sum += A[i];					/* sum relative abundances - dividing by this creates a */
	}								/* 		a probability function							*/

for (i=0; i<S; i++) A[i] /= sum;	/* set the area under the function to 1.0 */

return A;
}

/*
CALCULATES ZIPF DISTRIBUTION. This uses a linear decay between ln R and ln relative abundance
	where R is the rank of the taxon.  The area under this curve from 1ÉS then is used to standardize the 
	height so that it represents the proportions of taxon ranked 1ÉS.  
NEEDS:
    - M (slope in log-log space)
    - S (number of taxa)
RETURNS:
	- abundance: proportional abundances
***********************************************************************/
double *proportional_zipf_distribution(double M, int S)
{
int i = 0;							/* loop variable */
double sum = 0.0f;					/* number of "occurrences" */
double *A;							/* array to return */

if (M<0) { printf("\nproportional_zipf_distribution, illegal slope = %f\n",M); exit(1); }
if (S<=0) { printf("\nproportional_zipf_distribution, illegal number of taxa = %d\n",S); exit(1); }

A = dvector(S);						/* allocate array */
/*check this!	*/
sum = A[0] = 1;						/* a[0] = 100; everything else is M^i times less */
for (i=1; i<S; i++) {
	A[i] = pow(i+1,-1*M);				/* decay is log-log										*/
	sum += A[i];					/* sum relative abundances - dividing by this creates a */
	}								/* 		a probability function							*/

for (i=0; i<S; i++) A[i] /= sum;	/* set the area under the function to 1.0 */

return A;
}

/*
CALCULATES ZIPF-MANDELBROT DISTRIBUTION. This creates a "humped" curvilinear decay between ln R and ln relative abundance
	where R is the rank of the taxon.  The area under this curve from 1ÉS then is used to standardize the 
	height so that it represents the proportions of taxon ranked 1ÉS.  
NEEDS:
    - M (slope in log-log space; represents "predictability" of diversity)
    - beta (creates hump; represents "niche space")
    - S (number of taxa)
RETURNS:
	- abundance: proportional abundances
***********************************************************************/
double *proportional_zm_distribution(double M, double beta, int S)
{
int i = 0;							/* loop variable */
double sum = 0.0f;					/* number of "occurrences" */
double *A;							/* array to return */

if (M<0) { printf("\nproportional_zm_distribution, illegal slope = %f\n",M); exit(1); }
if (S<=0) { printf("\nproportional_zm_distribution, illegal number of taxa = %d\n",S); exit(1); }

A = dvector(S);						/* allocate array */

sum = 0;					/* a[0] = 100; everything else is M^i times less */
for (i=0; i<S; i++) {
/*	A[i] = A[0]/(pow(M,log(i+1+beta)));	/* decay is log-log										*/
	A[i] = pow(i+1+beta,-1*M);				
	sum += A[i];					/* sum relative abundances - dividing by this creates a */
	}								/* 		a probability function							*/

for (i=0; i<S; i++) A[i] /= sum;	/* set the area under the function to 1.0 */

return A;
}

/*
CALCULATES LOG-POWER DISTRIBUTION.  This is based on a delusion of Wagner's that distributions look like this sometimes.  
NEEDS:
    - C (coefficient on power function; a in a*x^b)
    - X (exponent on power function; b in a*x^b)
    - S (number of taxa)
RETURNS:
	- abundance: proportional abundances
***********************************************************************/
double *proportional_lp_distribution(double C, double X, int S) {
int i = 0;							/* loop variable */
double sum = 0.0f;					/* number of "occurrences" */
double x = 0.0f, y = 0.0f;			/* temp variables */
double *A;							/* array to return */

if (C<0.0f) { printf("\nproportional_lp_distribution, illegal C: C = %f, X = %f, S = %d\n",C,X,S); exit(1); }
if (X<0.0f) { printf("\nproportional_lp_distribution, illegal X: C = %f, X = %f, S = %d\n",C,X,S); exit(1); }
if (S<=0) { printf("\nproportional_lp_distribution, illegal S: C = %f, X = %f, S = %d\n",C,X,S); exit(1); }

A = dvector(S);						/* allocate array */

for (i=0; i<S; i++)	{
	x = i+1;
	y = -1*C*pow(x,X);
	A[i] = pow(10,y);				/* taxon i = formula */
	sum += A[i];					/* sum number of occurrences */
	}

for (i=0; i<S; i++) A[i] /= sum;	/* make proportional */

return A;
}


/*
CALCULATES ZERO SUM MULTINOMIAL DISTRIBUTION. This is modified from Olziewski's code.  
NEEDS:
    - theta ()
    - m ()
    - J (population size)
RETURNS:
	- abundance: proportional abundances
***********************************************************************/
double *proportional_zsm_distribution(double theta, double m, int J)
{
long int a, b, n, s1, s0, S;
double *phi;					/* Expected # of species with x specimens 			*/
double *sphi;					/* Sum of expected # of species with x specimens 	*/
double *whit;					/* Whittaker plot derived from phi					*/

phi=zerosum(theta,m,J);

sphi=dvector(J+1);
n=s0=0;
for (a=1; a<=J; ++a)	{
	sphi[a]=sphi[a-1]+phi[a];
	s1=sphi[a];
	if (s1>s0)	n+=(s1-s0)*(a+1);
	s0=s1;
	}

a=J;
S=sphi[a];
whit=dvector(S+1);
whit[S]=-1;
s0=0;
for (a=1; a<=J; ++a)	{
	s1=sphi[a];
	if (s1>s0)	{
		for (b=(S-s1); b<(S-s0); ++b)
			whit[b]=((double) (a+1))/((double) n);
		}
	s0=s1;
	}
free_dvector(phi);
free_dvector(sphi);
return whit;
}

/*
CALCULATES LOG-NORMAL DISTRIBUTION - this is a real lognormal distribution! 
	This based on dividing the area under a normal curve in S units.  It calculates this for one half of 
	the normal curve (using the standard erf equation for estimating area between 0 and X SD's away from 
	the mean).  It then uses that half to plot the positions along the X-axis on the remainder of the curve. 
	If it is untruncated, then it is a mirror image.  If it it truncated, then the untruncated half is
	calculated first and it is used to fill in the truncated part.  
	Note: the position along the X-axis gives the relative abundance as mag^X where M is the magnitude of 
	increase and X is the distance away from th first taxon.  
NEEDS:
    - mag (magnitude of increase between octaves)
    - trunc (1: the distribution is truncated; 0: untruncated)
    - mode (distance between the mode and the truncation - if positive, then you have the whole left side
    			and part of the right; if negative, then the opposite)
    - S (number of taxa)
RETURNS:
	- A: proportional abundances
***********************************************************************/
double *proportional_lgn_distribution (double mag, int S)
{
int		a;
double spar=0.0f;
long double x, m;
double	*A;

A=dvector(S);
for (a=0; a<S; ++a)	{
	x=((long double) (S-a))/((long double) (S+1));
	m=normsinv(x);
	A[a]=pow(mag,m);
	spar+=A[a];
	}
for (a=0; a<S; ++a)	A[a]/=spar;

return A;
}


/*
CALCULATES LOG-NORMAL DISTRIBUTION - this is a real lognormal distribution! 
	This based on dividing the area under a normal curve in S units.  It calculates this for one half of 
	the normal curve (using the standard erf equation for estimating area between 0 and X SD's away from 
	the mean).  It then uses that half to plot the positions along the X-axis on the remainder of the curve. 
	If it is untruncated, then it is a mirror image.  If it it truncated, then the untruncated half is
	calculated first and it is used to fill in the truncated part.  
	Note: the position along the X-axis gives the relative abundance as mag^X where M is the magnitude of 
	increase and X is the distance away from th first taxon.  
NEEDS:
    - mag (magnitude of increase between octaves)
    - trunc (1: the distribution is truncated; 0: untruncated)
    - mode (distance between the mode and the truncation - if positive, then you have the whole left side
    			and part of the right; if negative, then the opposite)
    - S (number of taxa)
RETURNS:
	- A: proportional abundances
***********************************************************************/
double *proportional_tlgn_distribution (double mag, int S, int St)
{
int		a;
double spar=0.0f;
long double x, m;
double	*A;

A=dvector(St);
for (a=0; a<St; ++a)	{
	x=((long double) (S-a))/((long double) (S+1));
	m=normsinv(x);
	A[a]=pow(mag,m);
	spar+=A[a];
	}
for (a=0; a<St; ++a)	A[a]/=spar;

return A;
}


/*
CALCULATES LOG-NORMAL DISTRIBUTION - this is a real lognormal distribution! 
	This based on dividing the area under a normal curve in S units.  It calculates this for one half of 
	the normal curve (using the standard erf equation for estimating area between 0 and X SD's away from 
	the mean).  It then uses that half to plot the positions along the X-axis on the remainder of the curve. 
	If it is untruncated, then it is a mirror image.  If it it truncated, then the untruncated half is
	calculated first and it is used to fill in the truncated part.  
	Note: the position along the X-axis gives the relative abundance as mag^X where M is the magnitude of 
	increase and X is the distance away from th first taxon.  
NEEDS:
    - mag (magnitude of increase between octaves)
    - trunc (1: the distribution is truncated; 0: untruncated)
    - mode (distance between the mode and the truncation - if positive, then you have the whole left side
    			and part of the right; if negative, then the opposite)
    - S (number of taxa)
RETURNS:
	- A: proportional abundances
***********************************************************************/
double *proportional_lgn_distribution_old2 (double mag, int trunc, double mode, int S)
{
int		a=0, b, c=0, d=0, s=0, odd;
double long	spar=0.0f;
double	rar;
double	x, rnd=0.0f, ttl=0.0f, lttl=-1.0f;
double	lx, dx, px;
double	ar=0.0f;
double	at=0.0f;
double	*A;

A=dvector(S);
odd=S%2;

at=1/((double) S);
x=lx=0.0;							/* 0.0 is the middle of the normal curve	*/
dx=0.001;							/* initial increment						*/
/* we start from the middle of the curve and work outwards; because it is symmetrical (or partially, if
	truncated), we can add in both directions until we are done	*/

a=at*1000000;
/* if there are an odd number of taxa, then the median taxon has A[median]=0 for now.  The adjacent	*/
/*		species histogram bars will start ar/2 to the left/right of 0.								*/
c=S/2;
if (odd==1)	{	
	A[(S/2)]=0.0f;
	b=ar=0.0f;							/* area under the curve from start to current x	*/
	rar=0.5*(1+erf(lx/pow(2,0.5)));		/* area to the right of the outset				*/
	while (b < (a/2))	{
		x+=dx;								/* slide right on x-axis							*/
		ar=0.5*(1+erf(x/pow(2,0.5)))-rar;	/* recalculate area between x & lx;					*/
		b=ar*1000000;						/* this is for precision: we never get it exact		*/
		if (b>(a/2))	{					/* we want the middle of the bar: we've overshot!	*/
			x=px;							/* move x back to where it was						*/
			dx/=2;							/* reduce dx										*/
			b=(a/2)-1;						/* reset b											*/						
			}	/* end overshot	*/
		px=x;								/* set prior x: we might need to return				*/
		}	/* keep going until ar ~= at	*/
	c=(S/2)+1;
	}	/* end initial routine for odd numbers of species	*/

for (d=(S/2)-1; d>=0; --d)	{
	dx=at;								/* amount to increment on x-axis				*/
	b=ar=0.0f;							/* area under the curve from start to current x	*/
	rar=0.5*(1+erf(lx/pow(2,0.5)));		/* area to the right of the outset				*/
	while (b < (a/2))	{
		x+=dx;								/* slide right on x-axis							*/
		ar=0.5*(1+erf(x/pow(2,0.5)))-rar;	/* recalculate area between x & lx;					*/
		b=ar*1000000;						/* this is for precision: we never get it exact		*/
		if (b>(a/2))	{					/* we want the middle of the bar: we've overshot!	*/
			x=px;							/* move x back to where it was						*/
			dx/=2;							/* reduce dx										*/
			b=(a/2)-1;						/* reset b											*/						
			}	/* end overshot	*/
		
		else if (b==a/2)	{
			A[d]=x;							/* tally position on X-axis that species d occupies		*/
			A[c]=-x;						/* tally position on X-axis that species c occupies	*/
			dx=at;							/* reset dx												*/
			while (b<a)	{
				x+=dx;
				ar=0.5*(1+erf(x/pow(2,0.5)))-rar;	/* recalculate area between x & lx;					*/
				b=ar*1000000;						/* this is for precision: we never get it exact		*/
				if (b>a && d>0)	{					/* we want the middle of the bar: we've overshot!	*/
					x=px;							/* move x back to where it was						*/
					dx/=2;							/* reduce dx										*/
					b=a-1;							/* reset b											*/						
					}	/* end overshot	*/
				if (d==0)	b=a;					/* if it is the last species, then do not bother	*/
				px=x;								/* set prior x: we might need to return				*/
				}	/* finish loop to get the second half of the bar		*/
			}	/* end search for halfway point	*/
		px=x;								/* set prior x: we might need to return				*/
		}	/* keep going until ar ~= at	*/
	px=lx=x;	/* we need to know the last x					*/
	++c;		/* increment the right half of the log-normal	*/
	}

for (b=0; b<S; ++b)	A[b]=pow(mag,A[b]);	/* now give abundances from octave values	*/

/* if trunc is zero, then this is a standard log-normal */
x=sumdvector(A,S);				/* calculate the area under the function	*/

for (b=0; b<S; ++b)	A[b]/=x;	/* set the area under the functon to 1.0	*/

return A;
}


/*
CALCULATES LOG-NORMAL DISTRIBUTION - this is a real lognormal distribution! 
	This based on dividing the area under a normal curve in S units.  It calculates this for one half of 
	the normal curve (using the standard erf equation for estimating area between 0 and X SD's away from 
	the mean).  It then uses that half to plot the positions along the X-axis on the remainder of the curve. 
	If it is untruncated, then it is a mirror image.  If it it truncated, then the untruncated half is
	calculated first and it is used to fill in the truncated part.  
	Note: the position along the X-axis gives the relative abundance as mag^X where M is the magnitude of 
	increase and X is the distance away from th first taxon.  
NEEDS:
    - mag (magnitude of increase between octaves)
    - trunc (1: the distribution is truncated; 0: untruncated)
    - mode (distance between the mode and the truncation - if positive, then you have the whole left side
    			and part of the right; if negative, then the opposite)
    - S (number of taxa)
RETURNS:
	- A: proportional abundances
***********************************************************************/
double *proportional_lgn_distribution_old (double mag, int trunc, double mode, int S)
{
int		a=0, b, d=0, s=0, odd, taxa, rem;
double long	spar=0.0f;
double	ha, her;
double	x, rnd=0.0f, ttl=0.0f, lttl=-1.0f;
double	lx, dx;
double	ar=0.0f;
double	*T;
double	*A;

A=dvector(S);
odd=S%2;

/* if trunc is zero, then this is a standard log-normal */
if (trunc==0)	{
	ar=1.0f;
	rem=taxa=(S/2);
	}

/* if trunc=1, then we are dealing with a truncated log-normal. 					*/
/* "mode" gives where it is truncated, with 0 being the normal mode,				*/
/* 1 being 1 SD to the right and -1 being 1 SD to the left							*/
/* If mode is negative, then you have the whole right side and part of the left.	*/
else	{
	ha=normheight(0,mode,1);
	if (mode>0)	her=erf(mode);
	else		her=erf(-1*mode);
	ha*=her;
	ar=0.5+ha;
	taxa=(((double) S) * 0.5)/ar;
	rem=S-taxa;
	}

spar=ar/((double) S);	/* area encompassed by any one species */

/* if there are an even number of taxa, then taxon (S/2) starts at spar/2 */
/* if there are an odd number of taxa, then taxon (S/2)+1/2 starts at spar */
T=dvector(taxa);

if (odd==0)
	T[0]=spar/2;
else	{
	T[0]=spar;
	A[taxa]=1;
	}			/**** START HERE! 	****/

for (a=1; a<taxa; ++a)	T[a]=T[a-1]+spar;

x=spar/(0.5*0.39894228);	/* distance left/right of mean of histogram bar occupying 1/S of the area	*/	

if (odd==1)	{
	T[(taxa/2)]=0;
	lx=0;
	}
else	{
	T[(taxa/2)]=x;
	T[(taxa/2)+1]=-x;
	lx=x;
	}

for (s=1; s<(taxa/2); ++s)	{
	
	}

/*	lttl=erf(ocs)*normheight(0,ocs,1);	*/
lx=lttl=0;
dx=0.5;
/* go through the first half of the taxa	*/ 
for (s=0; s<taxa; ++s)	{
	a=1; b=2;
	for (x=lx+dx; a!=b;)	{
//		y=erf(x);
//		z=normheight(0,x,1);
//		ttl=erf(x)*normheight(0,x,1);
//		ttl=y*z;
		ttl=(x-lx)*normheight(0,lx+((x-lx)/2),1);
		
		if ((lttl>T[s] && ttl<T[s]) || (lttl<T[s] && ttl>T[s]))
			dx/=-2;
/*		else
/*			dx/=2;	*/
		a=lttl*100000000;
		b=ttl*100000000;
		lttl=ttl;
		if ((x-lx)>1)	{
			x=x;
			}
		x+=dx;
		}
	A[taxa-(s+1)]=pow(mag,x);
	/* reset search parameters */
	dx=x-lx;
	lx=x;
	}

/* if mode>0, then is truncated on the right, then we have done the rare half first	*/
/* 	Therefore, we have to flip-flop all of this.							*/	
if (trunc==1 && mode>0)	for (s=0; s<taxa; ++s)	A[s]=1/A[s];

/* fill in the other half with the reciprocal of the first half 	*/
/*  Note: if truncated, then only part of the first half is used	*/
for (s=0; taxa+s+odd<S; ++s)	{
	A[taxa+s+odd]=1/A[taxa-(s+1)];
	}

x=sumdvector(A,S);				/* calculate the area under the function	*/

for (s=0; s<S; ++s)	A[s]/=x;	/* set the area under the functon to 1.0	*/

if (trunc==1 && mode>0)	A = dshellsort_dec(A,S);	/* flip-flop order if it began with the rarest */

free_dvector(T);
return A;
}

/*
draw_lgn_octaves - this shows where along a normal curve 1/Sth of the area is used.   
	This based on
NEEDS:
    - mag (magnitude of increase between octaves)
    - trunc (1: the distribution is truncated; 0: untruncated)
    - mode (position of the mode)
    - S (number of taxa)
RETURNS:
	- A: proportional abundances
***********************************************************************/
double *draw_lgn_octaves (int trunc, double mode, int S)
{
int		a=0, b, d=0, s=0, odd, taxa, rem;
double long	spar=0.0f;
double	ha, her;
double	x, rnd=0.0f, ttl=0.0f, lttl=-1.0f;
double	lx, dx;
/*double	lp=0.0000000f, lnc=0.0000000f, y=0.0000000f;	*/
double	ar=0.0f;
double	*T;
double	*ord;


ord=dvector(S);
odd=S%2;

/* if trunc is zero, then this is a standard log-normal */
if (trunc==0)	{
	ar=1.0f;
	rem=taxa=(S/2);
	}

/* if trunc=1, then we are dealing with a truncated log-normal. 					*/
/* "mode" gives where it is truncated, with 0 being the normal mode,				*/
/* 1 being 1 SD to the right and -1 being 1 SD to the left							*/
/* If mode is negative, then you have the whole right side and part of the left.	*/
else	{
	ha=normheight(0,mode,1);
	her=erf(mode);
	ha*=her;
	ar=0.5+ha;
	taxa=(((double) S) * 0.5)/ar;
	rem=S-taxa;
	}

spar=ar/((double) S);	/* area encompassed by any one species */

/* if there are an even number of taxa, then taxon (S/2) starts at spar/2 */
/* if there are an odd number of taxa, then taxon (S/2)+1/2 starts at spar */
T=dvector(taxa);

if (odd==0)
	T[0]=spar/2;
else	{
	T[0]=spar;
	ord[taxa]=1;
	}

for (a=1; a<taxa; ++a)	T[a]=T[a-1]+spar;

/*lttl=erf(ocs)*normheight(0,ocs,1);	*/
lx=lttl=0;
dx=0.5;
/* go through the first half of the taxa	*/ 
for (s=0; s<taxa; ++s)	{
	a=1; b=2;
	for (x=lx+dx; a!=b; x+=dx)	{
		ttl=erf(x)*normheight(0,x,1);
		if ((lttl>T[s] && ttl<T[s]) || (lttl<T[s] && ttl>T[s]))
			dx/=-2;
/*		else
/*			dx/=2;	*/
		a=lttl*100000000;
		b=ttl*100000000;
		lttl=ttl;
		}
	ord[taxa-(s+1)]=x;
	/* reset search parameters */
	dx=x-lx;
	lx=x;
	}

/* if mode>0, then is truncated on the right, then we have done the rare half first	*/
/* 	Therefore, we have to flip-flop all of this.							*/	
if (trunc==1 && mode>0)	for (s=0; s<taxa; ++s)	ord[s]*=-1;

/* fill in the other half with the reciprocal of the first half */
for (s=0; s<rem; ++s)	{
	ord[taxa+s+odd]=-1*ord[taxa-(s+1)];
	}

if (trunc==1 && mode>0)	ord = dshellsort_dec(ord,S);

free_dvector(T);
return ord;
}


/*
CALCULATES LOG-SERIES DISTRIBUTION - 
	This based on Hayek & Buzas' observation that the there is a linear relationship between sampled S
		and log sample size.  
NEEDS:
    - mag (magnitude of increase between specimens before we expect to sample another new species)
    - S (number of taxa)
RETURNS:
	- A: proportional abundances
***********************************************************************/
double *proportional_ls_distribution (double mag, int S)
{
int		a=0, s=0;
double	N, lN;
double	*A;


A=dvector(S);
for (s=0; s<S; ++s)	{
	N=pow(mag,s);
	A[s]=1;
	for (a=s-1; a>=0; --a)	A[a]=(N-1)*(A[a]/lN);
	lN=N;
	}

for (s=0; s<S; ++s)	A[s]/=N;
A = dshellsort_dec(A,S);

return A;
}


/*
CALCULATES BROKEN-STICK DISTRIBUTION - 
	This based on Hayek & Buzas' observation that there is a linear relationship between sampled S and N
		where N is sample size.  
NEEDS:
    - mag (how many times more speciemsn we need to sample the next species)
RETURNS:
	- A: proportional abundances
***********************************************************************/
double *proportional_bs_distribution (double mag, int S)
{
int		a=0, s=0;
double	N, lN;
double	*A;


A=dvector(S);
for (s=0; s<S; ++s)	{
	if (s==0)	N=1;
	else		N+=mag;
	A[s]=1;
	for (a=s-1; a>=0; --a)	A[a]=(N-1)*(A[a]/lN);
	lN=N;
	}

for (s=0; s<S; ++s)	A[s]/=N;

return A;
}


/*
CALCULATES BROKEN-STICK DISTRIBUTION - 
	This based on Hayek & Buzas' observation that there is a linear relationship between ln S and ln N
		where N is sample size.  
NEEDS:
    - mag (how many times more speciemsn we need to sample the next species)
RETURNS:
	- A: proportional abundances
***********************************************************************/
double *proportional_ln_distributionBH (double mag, int S)
{
int		a=0, s=0;
double	y, N, lN;
double	*A;


A=dvector(S);
for (s=0; s<S; ++s)	{
	y=log(s+1);

	N=pow(mag,y);
	A[s]=1;

	for (a=s-1; a>=0; --a)	A[a]=(N-1)*(A[a]/lN);
	lN=N;
	}

for (s=0; s<S; ++s)	A[s]/=N;
A = dshellsort_dec(A,S);

return A;
}


/*
CALCULATES RELATIVE ABUNDANCE DISTRIBUTION FROM EXPECTED NUMBERS OF SPECIES WITH X SPECIMENS DISTRIBUTION - 
NEEDS:
    - expected (expected number of species with x specimens)
    - max (total number of specimens)
RETURNS:
	- A: proportional abundances
***********************************************************************/
double *proportional_from_expt_inds(double *expected, int max)
{
int		a, b, c, r, t=0;
double	sum=0.499999f, S;
double *A;

S=sumdvector(expected,max);
r=S;
if (S-((double) r) >0.5)	++r;

A=dvector(r);

for (a=0; a<max && t<r; ++a)	{
	if (a>0)	expected[a]+=expected[a-1];
	if (expected[a]>sum)	{
		b=roundtoint(expected[a]);
		for (c=t; c<b; ++c)	A[c]=a;
		t=b;
		sum=0.5+t;
		}
	}

sum=0.0f;
sum=sumdvector(A,r);

for (a=0; a<r; ++a)	A[a]/=sum;

A=dshellsort_dec(A,r);
return A;
}

/* 	CALCULATES AN EXPECTED DISTRIBUTION GIVEN A PROPORTIONAL DISTRIBUTION AND A SAMPLE SIZE */
double* ideal_distribution(double *A, int N)	{
int i;

for (i=0; i<N; ++i)	A[i]=A[i]*((double) N);

return A;
}

/*
CALCULATES NORMAL DISTRIBUTION, 
NEEDS:
    - num_stdev - Number of divisions along normal curve to be used
    - oct_per_stdev - Octaves per SD (1 makes 1 Octave = 1 SD; 2 makes 2 Octaves = 1 SD
    - modal_oct - Octave with "mean" (=Octaves to the left of the "mean" on the normal curve.)
		NOTE: I am assuming that if someone enters "5," then they mean the fifth element, which is element 4;
			hence, we used modal_oct-1 now.  
RETURNS:
    - ARRAY GIVING NORMAL DISTRIBUTION BEGINING [START] SD's BEFORE THE MEAN
*************************************************************/
double* normdistodd(int num_stdev, int oct_per_stdev, int modal_oct) {
	int	a;
	int length = oct_per_stdev * num_stdev;
	double y;
	double Oct = -1.0 * (double) (modal_oct-1) / (double) oct_per_stdev;	/* changed to modal_oct-1 to accomodate start at zero */
	double Area = 0.0f;
	double *NormA;
	double p = pow(2*PI,0.5);

	NormA=dvector(length);

	/* Find the height of the histogram for each octave 		*/
	/* This will be used to determine the proportion of species */
	/* that fall into a category								*/
	for (a=0 ; a<length ; a++)	{
		y = exp(-(Oct*Oct)/2);
		y/= p;
		NormA[a] = y;
		Area+= y;
		Oct+= 1.0/ (double) oct_per_stdev;
		}
	
	for (a=0; a<(length); a++)
		NormA[a] /= Area;
	
	return NormA;
}

/*
CALCULATES NORMAL DISTRIBUTION HISTOGRAM WITH SYMMETRICAL BARS AROUND MEAN,
	THIS MEANS THAT THE HEIGHT 
NEEDS:
    - num_stdev - Number of divisions along normal curve to be used
    - oct_per_stdev - Octaves per SD (1 makes 1 Octave = 1 SD; 2 makes 2 Octaves = 1 SD
    - modal_oct - Octave with "mean" (=Octaves to the left of the "mean" on the normal curve.)
		NOTE: I am assuming that if someone enters "5," then they mean the fifth element, which is element 4;
			hence, we used modal_oct-1 now.  
RETURNS:
    - ARRAY GIVING NORMAL DISTRIBUTION BEGINING [START] SD's BEFORE THE MEAN
*************************************************************/
double* normdistevn(int num_stdev, int oct_per_stdev, int modal_oct) {
	int	a;
	int length = oct_per_stdev * num_stdev;
	double y;
	double Oct = 0.5+(-1.0 * (double) (modal_oct-1) / (double) oct_per_stdev);	/* changed to modal_oct-1 to accomodate start at zero */
	double Area = 0.0f;
	double *NormA;
	double p = pow(2*PI,0.5);

	NormA=dvector(length);

	/* Find the height of the histogram for each octave 		*/
	/* This will be used to determine the proportion of species */
	/* that fall into a category								*/
	for (a=0 ; a<length ; a++)	{
		y = exp(-(Oct*Oct)/2);
		y/= p;
		NormA[a] = y;
		Area+= y;
		Oct+= 1.0/ (double) oct_per_stdev;
		}
	
	for (a=0; a<(length); a++)
		NormA[a] /= Area;
	
	return NormA;
}

/***************************************************************************
zerosum: calculates the expected number of species with x specimens.
	requires:
		theta: 
		m: "seed" probability of finding taxa
		J: approximate number of individuals
		mxn: maximum sampled abundance - calculate only twice beyond this to save time
			(added by Wagner 01.03.2005)
	modified from code written by T. Olszewski 2003
***************************************************************************/
double *zerosum(double theta, double m, long int J)
{
double gam,coef,*phi;
double x, y, z;
long int a;

/*n=2*mxn;
if ((2*mxn)>J)	n=J;	
n=J;*/

phi=dvector(J+1);
gam = m*(J-1)/(1-m);
x=gammln(J+1);
y=gammln(gam);
z=gammln(J+gam);

phi[0]=0.0f;
for (a=1; a<=J; ++a)	{
	coef = (log(theta)+x+y-gammln(a+1)-gammln(J-a+1)-z);
	phi[a] = qromo(volkov,theta,gam,a,J,coef,0.0,gam,midpnt);
	
	/* Wagner 2005.0331 - if the remainder cannot expect to yield another species, then stop	*/
	if ((phi[a]*((double)(J-a))) < 0.5)	a=J;
	}
/*fprintf(fpout,"1\t%e\t%e\n",phi,S);	*/

return phi;
}

/***************************************************************************
From Volkov et al. 2003 Nature 424:1035-1037
***************************************************************************/
double volkov(double theta, double gam, int n, long int J, double coef, double y)
{
return exp(gammln(n+y)+gammln(J-n+gam-y)-gammln(1+y)-gammln(gam-y)-y*theta/gam+coef);
}


/***************************************************************************
This routine computes the nth stage refinement of an extended midpoint rule.
	func is input as a pointer to the function to integrated between limits a and b,
	also input.  When called with n=1, thie routine returns the crudest estimate of
	the integral of f(x)dx over the interval a to b.  Subsequent calls with n=2,3,...
	(in that sequential order) will improve the accuracy of s by adding (2/3)x3^(n-1)
	additional interior points.  s should not be modified between sequential calls.

	written by T. Olszewski 2003
***************************************************************************/
double midpnt(double (*func)(double,double,int,long int,double,double), double theta, double gam, int nn, long int J, double coef, double a, double b, int n)
{
double x,tnm,sum,del,ddel;
static double s;
int it,j;

if (n == 1) {
	return (s=(b-a)*FUNC(theta,gam,nn,J,coef,0.5*(a+b)));
	} 
else {
	for (it=1,j=1;j<n-1;j++) it *=3;
	tnm=it;
	del=(b-a)/(3.0*tnm);
	ddel=del+del;
	x=a+0.5*del;
	sum=0.0;
	for (j=1;j<=it;j++) {
		sum += FUNC(theta,gam,nn,J,coef,x);
		x += ddel;
		sum += FUNC(theta,gam,nn,J,coef,x);
		x += del;
		}
	s=(s+(b-a)*sum/tnm)/3.0;
	return s;
	}
}

/***************************************************************************
qromo: Romberg integration on an open interval.  
	Returns the integral of the function func from a to b, using any specified integrating function choosen and 
	Romberg's method.  Normally choose will be an open formula, not evaluating the function at the endpoints.  
	It is assumed that choose triples the number of steps on each call, and that its error series contains only 
	even powers of the number of steps. The routines midpnt, midinf, midsql, midsqu, medexp, are possible choices 
	for choose.  The parameters have the same meaning as in qromb.

	written by T. Olszewski 2003
***************************************************************************/
double qromo(double (*func)(double,double,int,long int,double,double), double theta, double gam, int nn, long int J, double coef, double a, double b, double (*choose)(double(*)(double,double,int,long int,double,double), double, double, int, long int, double, double, double, int))
{
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
int j;
double ss,dss,h[JMAXP+1],s[JMAXP+1];

h[1]=1.0;
for (j=1;j<=JMAX;j++) {
	s[j]=(*choose)(func,theta,gam,nn,J,coef,a,b,j);
	if (j >= K) {
		polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
		if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
	s[j+1]=s[j];
	h[j+1]=h[j]/9.0;	/*this is where the assumption of step tripling and an even error series is used*/
	}
return 0.0;				/*never get here*/
}

/***************************************************************************
Given arrays  xa[1..n] and ya[1..n], and given a value x, this routine returns a value y, and an error estimate dy.  
	If P(x) is the polynomial of degree N-1 such that P(xai)=yai, i=1,...,n, then the returned value y = P(x).

	written by T. Olszewski 2003
***************************************************************************/
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
int i,m,ns=1;
double den,dif,dift,ho,hp,w;
double *c,*d;

dif=fabs(x-xa[1]);
c=dvector(n+1);
d=dvector(n+1);

for (i=1;i<=n;i++) {				/*here we find the index ns of the closest table entry,*/
	if ( (dift=fabs(x-xa[i])) < dif) {
		ns=i;
		dif=dift;
		}
	c[i]=ya[i];        			/*and initialize the tableau of c's and d's.*/
	d[i]=ya[i];
	}
*y=ya[ns--];					/*This is the initial approximation to y.*/
for (m=1;m<n;m++) {				/*For each column of the tableau,*/
	for (i=1;i<=n-m;i++) {		/*we loop over the current c's and d's and update them.*/
		ho=xa[i]-x;
		hp=xa[i+m]-x;
		w=c[i+1]-d[i];
		den=ho-hp;
		den=w/den;
		d[i]=hp*den;
		c[i]=ho*den;
		}
	*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	/*After each column in the tableau is completed, we decide which correction,
	c or d, we want to add to our accumulating value of y, i.e., which path to
	take through the tableau--forking up or down.  We do this in such a way as
	to take the most "straight line" route through the tableau to its apex,
	updating ns accordingly to keep track of where we are.  This route keeps the
	partial approximation centered (insofar as possible) on the target of x.
	The last dy added is thus the error indication.*/
	}
free_dvector(c);
free_dvector(d);
}

/* stateschangefromcharacterchange: given probabilistic numbers of changes per character, determine the
		probabilistic number of transitions to the nth state;
		NOTE: for binary characters, this regurgitates exactly what was entered	
Requires:
	PComp: an nchar x (mxstp+1) matrix giving the probability that a character changed N times
	nstates: number of states per character;
	nchar: total number of characters
	mxstp: maximum number of steps
***************************************************************************************************/
double *stateschangefromcharacterchange(double **Pcomp, int *nstates, int nchars, int mxstp)
{
int		a, c, s;
double	x;
double	*zs;

zs=dvector(mxstp+1);

for (c=0; c<nchars; ++c)	{
	/* at s steps, determine the expected number of states with 1, 2, 3, etc., changes	*/
	for (s=(nstates[c]-1); s<mxstp; ++s)	{
		/* if binary, then it's just the probability of change	*/
//		x=Pcomp[c][s];
		if (nstates[c]<=2)	zs[s]+=Pcomp[c][s];
		/* figure out the distribution of the changes per state for N+1 state characters 	*/
		/* do this by assuming that N steps for N states represents one change per state	*/
		/* assume that each state is sampled with P=1/N for each extra step.				*/
		/* Thus, for a 3-state (N=2) character and any one state:
			P[1 derivation | N steps] = 1.0 x P[N steps]
			P[1 derivation | N+1 steps] = 0.50 x P[N+1 steps]
			P[2 derivation | N+1 steps] = 0.50 x P[N+1 steps]
			P[1 derivation | N+2 steps] = 0.25 x P[N+2 steps]
			P[2 derivation | N+2 steps] = 0.50 x P[N+2 steps]
			P[3 derivation | N+2 steps] = 0.25 x P[N+2 steps]
			
			NOTE: P[1 derivation | N+X steps] = P[0 derivation | X steps]
			
			So, at X steps, the question is, what is P[Y-N | X-N]?

		/* a will be the number of derivations per state */
		/* what we have to do at 3 extra steps is find out the distribution of states
				with 1, 2 and 3 extra steps.  Note that this is a waste of time for
				binary characters, where the "other" step is always activated		*/
		/************************************************************************************/
		else if (nstates[c]>2)	{
			/* this starts at 0 because we need to calculate with no "extra" steps	*/
			for (a=0; a<=s-(nstates[c]-1); ++a)	{
				/* this is the probability of the state being derived a "extra" times	*/
				if (s==(nstates[c]-1))
/*					dummy1[a][s]=dummy1[s][a]=xx[a]=1.0f;	*/
					x=1.0f;
				else
/*					dummy1[a][s]=dummy1[s][a]=xx[a]=binomexact(a,(s-(nstates[c]-1)),((double) 1/(nstates[c]-1)));	*/
					x=binomexact(a,(s-(nstates[c]-1)),((double) 1/(nstates[c]-1)));
				/* state must be observed once, so observations are 1+a							*/
				/* this gives the "fraction" of states from a character with a+1 derivations	*/
/*				dummy2[a][s]=dummy2[s][a]=(nstates[c]-1)*xx[a]*Pcomp[c][s];	*/
/*				zs[a+1]+=(nstates[c]-1)*xx[a]*Pcomp[c][s];	/* note: this is a+1 because they all have at least one change	*/
				zs[a+1]+=(nstates[c]-1)*x*Pcomp[c][s];	/* note: this is a+1 because they all have at least one change	*/
/*				ps[c]+=(nstates[c]-1)*xx[a]*Pcomp[c][s];	*/
				}	/* end P[a+1 state steps | s character steps]	*/
			}	/* end case of multistate character	*/
		}	/* end search through different numbers of steps	*/
	}

return zs;
}