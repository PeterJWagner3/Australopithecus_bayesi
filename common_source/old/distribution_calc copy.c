#define distribution_calc
#include "distribution_calc.h"
#include "memory.h"
#include "more_math.h"
#include "Probability.h"
#include "sort.h"

/*
CALCULATES GEOMETRIC DISTRIBUTION. 
NEEDS:
    - M (slope in log-linear space)
    - S (number of taxa)
RETURNS:
	- abundance: proportional abundances
***********************************************************************/
double *proportional_gs_distribution(double M, int S)
{
int i = 0;							/* loop variable */
double sum = 0.0f;					/* number of "occurrences" */
double *A;							/* array to return */

if (M<1.0f) {
	printf("\nproportional_gs_distribution, illegal slope = %f\n",M);
	exit(1);
	}
if (S<=0)	{
	printf("\nproportional_gs_distribution, illegal number of taxa = %d\n",S);
	exit(1);
	}

A = dvector(S);						/* allocate array */

sum = A[0] = 100;					/* taxon 0 = 100 occurrences */
for (i=1; i<S; i++) {
	A[i] = A[i-1] / M;				/* taxon i = taxon (i-1) / slope */
	sum += A[i];					/* sum number of occurrences */
	}

for (i=0; i<S; i++) A[i] /=  sum;	/* make proportional */

return A;
}


/*
CALCULATES ZIPF-MANDELBROT DISTRIBUTION. 
NEEDS:
    - M (slope in log-log space)
    - S (number of taxa)
RETURNS:
	- abundance: proportional abundances
***********************************************************************/
double *proportional_zm_distribution(double M, int S)
{
int i = 0;							/* loop variable */
double sum = 0.0f;					/* number of "occurrences" */
double *A;							/* array to return */

if (M<1) { printf("\nproportional_zm_distribution, illegal slope = %f\n",M); exit(1); }
if (S<=0) { printf("\nproportional_zm_distribution, illegal number of taxa = %d\n",S); exit(1); }

A = dvector(S);						/* allocate array */

sum = A[0] = 100;					/* taxon 0 = 100 occurrences */
for (i=1; i<S; i++) {
	A[i] = A[0]/(pow(M,log(i+1)));
	sum += A[i];					/* sum number of occurrences */
	}

for (i=0; i<S; i++) A[i] /= sum;	/* make proportional */

return A;
}

/*
CALCULATES LOG-POWER DISTRIBUTION. 
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
CALCULATES LOG-NORMAL DISTRIBUTION. 
NEEDS:
    - SD (number of standard devations of lognormal to use)
    - Oct (number of octaves per standard deviation)
    - S (number of taxa)
RETURNS:
	- abundance: proportional abundances
***********************************************************************/
double *proportional_ln_distribution(int num_stdev, int oct_per_stdev, int modal_octave, double mag_per_oct, int S) {

double *ddistribution;
double *pdistribution;
int *idistribution;

if (num_stdev<=0) { printf("proportional_ln_distribution, illegal num_stdev= %d\n",num_stdev); exit(1); }
if (oct_per_stdev<=0) { printf("proportional_ln_distribution, illegal oct_per_stdev= %d\n",oct_per_stdev); exit(1); }
if (modal_octave<0) { printf("proportional_ln_distribution, illegal modal_octave= %d\n",modal_octave); exit(1); }
if (mag_per_oct<0.0f) { printf("proportional_ln_distribution, illegal mag_per_oct = %f\n",mag_per_oct); exit(1); }
if (S<=0) { printf("proportional_ln_distribution, illegal S= %d",S); exit(1); }

ddistribution = normdistodd(num_stdev, oct_per_stdev, modal_octave);
idistribution = OctaveRichness2(S, ddistribution, num_stdev, oct_per_stdev);
pdistribution = LNAbundCalc(S, num_stdev, oct_per_stdev, mag_per_oct, idistribution);

free_dvector(ddistribution);
free_ivector(idistribution);

return  pdistribution;
}


/*
CALCULATES FAUX LOG-NORMAL DISTRIBUTION - this is a geometric distribution with a shifting decay rate. 
Note: this was initially written for an untruncated log-normal, where the mode is the median rank.  This now allows for truncated log-normals
NEEDS:
    - M (slope at the mode)
    - dM (initial slope)
    - S (number of taxa)
    - mode (position of the mode)
RETURNS:
	- abundance: proportional abundances
***********************************************************************/
double *proportional_fln_distribution(double M, double dM, int S, double mode)
{
int i = 0;							/* loop variable */
double sum = 0.0f;					/* number of "occurrences" */
double x = 0.0f, y = 0.0f;			/* temp variables */
double *A;							/* array to return */
/*double	median;	*/			

A=dvector(S);

sum=A[0]=100.0f;
for (i=0; i<S-1; ++i)	{
	x=(mode-((double) i))/mode;
	if (x<0)	x=-1*x;
	
	y=M+(x*(dM-M));
	A[i+1]=A[i]/y;

	sum=sum+A[i+1];
	}

for (i=0; i<S; ++i)	A[i]=A[i]/sum;

/* if M <1, then the slope will rise at some point */
if (M<1.0)	A = dshellsort_dec(A,S);

return(A);
}

/*
CALCULATES LOG-NORMAL DISTRIBUTION. 
This is an iterative solution, which avoids weird effects of rounding error
	Octave Richness is calculated at lower richness, then higher richness is added to that Octave Richness
	This means that if there are 3 in Octave 2 at 100 taxa, there will be at least 3 in Octave 2 at 105 taxa.
		NOTE:  THIS IS NOT TRUE WHEN USING proportional_ln_distribution!!!!! 
NEEDS:
    - SD (number of standard devations of lognormal to use)
    - Oct (number of octaves per standard deviation)
    - S (number of taxa)
RETURNS:
	- abundance: proportional abundances
***********************************************************************/
double *proportional_ln_distribution2(int num_stdev, int oct_per_stdev, int modal_octave, double mag_per_oct, int S, int *idistribution) {
	
	double *ddistribution;
	double *pdistribution;
	
	if (num_stdev<=0) { printf("proportional_ln_distribution, illegal num_stdev= %d\n",num_stdev); exit(1); }
	if (oct_per_stdev<=0) { printf("proportional_ln_distribution, illegal oct_per_stdev= %d\n",oct_per_stdev); exit(1); }
	if (modal_octave<0) { printf("proportional_ln_distribution, illegal modal_octave= %d\n",modal_octave); exit(1); }
	if (mag_per_oct<0.0f) { printf("proportional_ln_distribution, illegal mag_per_oct = %f\n",mag_per_oct); exit(1); }
	if (S<=0) { printf("proportional_ln_distribution, illegal S= %d",S); exit(1); }

	ddistribution = normdistodd(num_stdev, oct_per_stdev, modal_octave);
	pdistribution = LNAbundCalc(S, num_stdev, oct_per_stdev, mag_per_oct, idistribution);

	free_dvector(ddistribution);
	free_ivector(idistribution);

	return  pdistribution;
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

/*
CALCULATES SPECIES PER LOG-NORMAL OCTAVE
NEEDS:
    - R - richness (R)
    - NormA - Number of divisions along normal curve to be used
    - num_stdev - Number of divisions along normal curve to be used
    - oct_per_stdev - Octaves per SD (1 makes 1 Octave = 1 SD; 2 makes 2 Octaves = 1 SD
RETURNS:
    - ARRAY GIVING NUMBER OF SPECIES PER OCTAVE
****************************************************************************/
/*** NOTE - THIS DOES NOT AVOID ROUNDING ERROR - USE OctaveRichnes2 FOR THAT
****************************************************************************/
int *OctaveRichness(int R, double *NormA, int num_stdev, int oct_per_stdev) {
int	a, mode;
int c = 0;
int length = oct_per_stdev * num_stdev;
int *Rich, *round;
double	*deviat;
double	w=0.0f, x=0.0f, y=R, z=0.0f;

Rich=ivector(length);
deviat=dvector(length);
for (a=length-1; a>=0; --a)	{
	/* do not round yet - we will do that later */
	Rich[a]=x=NormA[a]*y;
	deviat[a]=x-Rich[a];
	c+=Rich[a];
	}
	
if (c<R)	{
	/* find modal octave */
	if (NormA[0]>NormA[1])
		mode=0;
	else	{
		for (a=1; a<length; ++a)	{
			if (NormA[a]>NormA[a-1] && NormA[a]>NormA[a+1])	{
				mode=a;
				a=length;
				}
			}
		}
	/* modal octave found */

	round=ivector(length);
	round=sort_decintbydouble(round,deviat,length);
	
	for (a=0; a<(R-c); ++a)	{
		++Rich[round[a]];
		}			
	
	free_ivector(round);
	}
free_dvector(deviat);
return Rich;
}

/*CALCULATES SPECIES PER LOG-NORMAL OCTAVE
NEEDS:
    - R - richness (R)
    - NormA - Number of divisions along normal curve to be used
    - num_stdev - Number of divisions along normal curve to be used
    - oct_per_stdev - Octaves per SD (1 makes 1 Octave = 1 SD; 2 makes 2 Octaves = 1 SD
RETURNS:
    - ARRAY GIVING NUMBER OF SPECIES PER OCTAVE
*************************************************************/
/* Now, rescale NormA so that it gives the number of species in each octave */
/* this is done from the top down to avoid rounding errors */
int* OctaveRichness2(int R, double *NormA, int num_stdev, int oct_per_stdev) {
int	a, mode;
int c = 0;
int length = oct_per_stdev * num_stdev;
int *Rich, *round;
double	*deviat, *prop;
double	v=100.0f, w=0.0f, x=0.0f, y=R, z=0.0f;

Rich=ivector(length);
deviat=dvector(length);
prop=dvector(length);

for (a=length-1; a>=0; --a)	{
	/* do not round yet - we will do that later */
	Rich[a]=x=NormA[a]*v;
	deviat[a]=x-Rich[a];
	c+=Rich[a];
	}
	
if (c<100)	{
	/* find modal octave */
	if (NormA[0]>NormA[1])
		mode=0;
	else	{
		for (a=1; a<length; ++a)	{
			if (NormA[a]>NormA[a-1] && NormA[a]>NormA[a+1])	{
				mode=a;
				a=length;
				}
			}
		}
	
	round=ivector(length);
	round=sort_decintbydouble(round,deviat,length);
	
	for (a=0; a<(100-c); ++a)	{
		++Rich[round[a]];
		}			
	
	free_ivector(round);
	}

/* Now, add additional taxa to the distribution 																	*/
/* do this by finding out which octaves deviate most from the expectations of a Gaussian curve given new richness	*/
/* This eliminates rounding error and makes the normal curve for X+5 a modification of the curve for X rather		*/
/*		than possibly a very different curve (as happens due to rounding error */
if (R>100)	{
	for (v=105; v<=y; v=v+5)	{
		/* find which octaves deviate most from expectation */
		for (a=0; a<length; ++a)	{
			prop[a] = (x=Rich[a])/v;
			deviat[a] = NormA[a] - prop[a];
			}

		round=ivector(length);
		round=sort_decintbydouble(round,deviat,length);

		for (a=0; a<5; ++a)
			++Rich[round[a]];
		
		free_ivector(round);
		}
	}

free_dvector(deviat);
free_dvector(prop);
return Rich;
}

/*
CALCULATES TRUCNATED LOG NORMAL ABUNDANCE DISTRIBUTION, 
NEEDS:
    - num_stdev - Number of divisions along normal curve to be used
    - oct_per_stdev - Octaves per SD (1 makes 1 Octave = 1 SD; 2 makes 2 Octaves = 1 SD
    - R - richness (R)
    - LS - Log scale on octaves.  A species abundance is LS^Octave
    - Start - Number of octaves to the left of the "mean" on the normal curve.  
RETURNS:
    - ARRAY CONTAINING LOG NORMAL DISTRIBUTION
*************************************************************/
double* LNAbundCalc(int R, int num_stdev, int oct_per_stdev, double LS, int *Rich)	{
int	a, b, c;
int T = 0;
int length = num_stdev * oct_per_stdev;
double	x, z;
double Area = 0.0f;
double	*A;

A=dvector(R);

/* Now, calculate numbers for each species */
/* divide each octave by species in it; abundance is LS^thatnumber */
for (a=0; a<length; ++a)	{
	x=0;
	for (b=0; b<Rich[a]; ++b)	{
		x+= 1.0 / (double) Rich[a];
		z = pow(LS,(double) (a+x));
		Area+= z;
		
		/* sort abundances as you go (initial abundances less than subsequent ones) */
		for (c=T-1; z>A[c] && c>=0; --c)	{
			A[c+1]=A[c];
			}

		A[c+1]=z;
		
		++T;
		if (T>=R)	{
			a = length;
			b = Rich[a];
			}
		}
	if (Rich[a+1]==0 && T>(R/2)) {
		a = length;
		}
	}
	
for (a=0; a<R; ++a)	{
	A[a]/= Area;
	}

return A;
}

/* 	CALCULATES AN EXPECTED DISTRIBUTION GIVEN A PROPORTIONAL DISTRIBUTION AND A SAMPLE SIZE */
double* ideal_distribution(double *A, int N)	{
int i;

for (i=0; i<N; ++i)	A[i]=A[i]*((double) N);

return A;
}


double *proportional_lgn_distribution (double mag, double mode, int S)
{
int		a=0, b, c, d=0, n, spec=0, in=1, pn=0;
double	ra, la, rer, ler;
double	oct, moct, x, rnd=0.0f, ttl=0.0f, lttl=-1.0f, kttl;
/*double	d1, d2, h;*/
double	lp=0.0000000f, lnc=0.0000000f, y=0.0000000f;
double	fn, ar=0.0f;
double	*A;

/* make sure that mxsp>=mxfds */
/* if S<100, then there are 5 octaves that should have species */
if (S<100)	{
	moct=mode+2.5;
	/* if the distribution is truncated, then we need to adjust the area under the curve	*/
	/* dividing each fn below by this increases it to the proportion of the possible area	*/
/*	for (x=0; x<=moct+3.5; x+=0.125)	{
		 y=normheight(x,moct,1);
		 ar+=y/8;
		 }	*/
	/* find the area on the left side of the bell curve */
	/* if the distribution is truncated, then there will be less of this than otherwise 		*/
	/* even if not, then because we begin 2.5 SD's away, we want to eliminate that small amount */
	if (moct>0)	{
		la=normheight(0,moct,1);
		ler=erf(moct);
		la*=ler;
		}
	else	la=0.0f;
	
	if (mode<=0)	{
		ra=normheight(0,2.5,1);
		rer=erf(2.5);
		ra*=rer;
		}
	else	{
		ra=normheight(0,mode+2.5,1);
		rer=erf(mode+2.5);
		ra*=rer;
		}
	ar=ra+la;
	}
/* if S>1000, then there are 6 octaves that should have species */
else	{
	moct=mode+3.0f;
	for (x=0; x<=moct+3.0; x+=0.125)	{
		 y=normheight(x,moct,1);
		 ar+=y/8;
		 }
	}

y=0.0000000f;
A=dvector(S);

in=1;
/* find the probability of species beginning with n specimens*/
for (n=1; a<S; n=n+in)	{
	/* calculate how far along the x-axis we are given that the axis really is in a log-scale */
	oct=log(n)/log(mag);
	x=0.5+n;
	/* find the area within the histogram centred on n - that is the expected frequency of species with n specimens */
/*	h=normheight(oct,moct,1);
	d1=log(x)/log(mag);
	d2=log(x-1)/log(mag);
	fn=h*(d1-d2);
	fn=((double) S)*fn;*/

	if (in==1)	{
/*		fn=(1/ar)*((double) S)*normheight(oct,moct,1)*((log(x)/log(mag))-(log(x-1)/log(mag)));	*/
		/* estimate the total number of expected individuals */
		/* then see if we expect the next species to fall on this point */
		la=log(x-1)/log(mag);
		ra=log(x)/log(mag);
		fn=(1/ar) * ((double) S) * normareabetween(la,ra,moct,1);
		ttl+=fn;
		a=ttl;
		b=lttl;

/*		if (((ttl-a)>=0.5 && (lttl-b)<0.5) || (a-b)>1)	{	*/
		if (ttl>=(((double) d)+0.5) && lttl<(((double) d)+0.5))	{
			if ((ttl-a)>=0.5)	++a;
			if (b<0)			b=0;
			for (c=d; c<a; ++c)	pn=A[c]=n;
			spec+=((a-d)*n);
			d=a;		/* this shouldn't be necessary, but it is..... */
			if (a==1 && A[0]>=10)				in=A[0]/2;
			else if (a>1 && (A[a-1]-A[a-2]>=10))in=(A[a-1]-A[a-2])/2;
			else								in=1;
			}
		}
		
	else	{
		y=(n-in)+0.5;
		x=n+0.5;
		la=log(y)/log(mag);
		ra=log(x)/log(mag);
/*		fn=normareabetween(la,ra,moct,1);
		fn*=(1/ar);
		fn*=((double) S);	*/
		fn=(1/ar) * ((double) S) * normareabetween(la,ra,moct,1);
		ttl+=fn;
/*		d=1+lttl;	*/
		a=ttl;
		b=lttl;
		/* now if we overshot X.5, we need to back up to the n at which it happened */
/*		if (((ttl-a)>=0.5 && (lttl-b)<0.5) || (ttl-lttl >= ((((double) d)+0.5)-lttl)))	{	*/
		if (ttl>=(((double) d)+0.5) && lttl<(((double) d)+0.5))	{
			c=kttl=ttl;					/* kttl will be used to keep track of prior ttls	*/
			ttl-=fn;					/* reset ttl			*/
			in/=-2;						/* we need to go backwards in smaller steps */
			while (n+in<pn)	in/=2;
/*			if (n-1==A[a-1])	in=0;	*/
			if (n-1==pn)	in=0;
			
			while (abs(in)>0)	{
				while ((n+in)<pn)	in/=2;
				if (in==0)	in=1;
				n+=in;
				x=0.5+n;
/*				la=log(y)/log(mag);	*/
				ra=log(x)/log(mag);
/*				fn=normareabetween(la,ra,moct,1);
				fn*=(1/ar);
				fn*=((double) S);	*/
				fn=(1/ar) * ((double) S) * normareabetween(la,ra,moct,1);
				ttl+=fn;
				
				/* whenever the current and previous area are on either side of X.5, reverse a half-step */
				if ((ttl-a)<0.5 && (kttl-c)>0.5 || (ttl-a)>0.5 && (kttl-c)<0.5)	{	
					in/=-2;						/* we need to go backwards in smaller steps */
					if (n==(pn+1) && (ttl-a)>=0.5)	in=0;
					/* once this reaches 0, we are there */
					if (in==0)	{
						/* if ttl is below X.5, then it is n+1 that we want */
						if ((ttl-a)<0.5)	{
							++n;
							ttl=kttl;
							}
						pn=A[a]=n;
						d=a+1;
						spec+=n;
						}
					}
				c=kttl=ttl;
				if (in!=0)	ttl-=fn;
				/* do not let n become lower than the n prior to X.5 being passed */
				}
			
			if (a==1 && A[0]>=10)				in=A[0]/2;
			else if (a>1 && (A[a-1]-A[a-2]>=10))in=(A[a]-A[a-1])/2;
			else								in=1;
			}
		}

	pn=n;
	lttl=ttl;
	}

for (n=0; n<S; ++n)
	A[n]/=((double) spec);

A = dshellsort_dec(A,S);

return A;
}