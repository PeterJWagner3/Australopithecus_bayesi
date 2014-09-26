#define info_theory
#include "infotheory.h"
#include "sort.h"

/* calculates Akaike's Information Criterion (aic)
	L = model's log likelihood
	k = number of parameters.
	2 * L + 2 * k
	
Akaike H. 1973 Information theory and an extension of the maximum likelihood principle. 2nd
International Symposium on Information Theory 267-281, .
***********************************************************************/
double calc_aic( int k, double L)	{
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

/* calculates Bayesian Information Criterion (BIC) for sent data.
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


/* finds best partitions given a matrix of log-likelihoods.
	loglike: cases x hyp matrix of log-likelihoods, where loglike[c][h] gives the loglikelihood of h for case c given the data
	cases: number of individual partitions that might share the same "rate" or other similar hypothesis (e.g., characters in a clade)
	hyp: number of rate hypotheses used for the matrix
 k = number of parameters.
 n = sample size.
 
 Burnham & Anderson 1998 p. 68 suggest not using this one...
 ***********************************************************************/
double **testpartitions(double **loglike, int cases, int hyp, double *hypos, double dp, int n)
{
int	c1, c2, c3, h, p, r1, r2, r3;
int	b, m, mp=1;
double	mlmpr=0.0f;
double	y;
int	*rnkcase, *mlrhyp, *prtcse;
int	**parts, **mlparts;
double	*mlcase, *mlchyp, *prtlgl, *prtaicc;
double	**partitionlgl, **partitionhyp;
double	**product;

mlcase=dvector(cases);						/* most likely hypothesis for each character		*/
mlchyp=dvector(cases);						/* maximum log-likelihood for each character		*/
mlrhyp=ivector(cases);						/* most likely hypothesis rank for each character	*/
prtcse=ivector(cases);						/* vector used to print original character numbers	*/
cleardvector(mlchyp,cases,-1*RAND_MAX);						

parts=imatrix(cases+1,cases+1);				/* matrix giving first ranked character of each partition, with parts[x][0…x] giving the divisions for a 3-parameter hypothesis 	*/
mlparts=imatrix(cases+1,cases+1);			/* matrix giving first ranked character of each ml partition, with mlparts[x][0…x] giving the divisions for the best 3-parameter hypothesis 	*/
/*parts[0][0]=1;							/* first element gives number of partitions								*/
/*parts[0][1]=0;							/* elements 1…partitions gives the first character of each partition	*/
prtlgl=dvector(cases+1);					/* prtlgl[x] gives the log-likelihood of x partitions					*/
prtaicc=dvector(cases+1);					/* prtlgl[x] gives the log-likelihood of x partitions					*/
prtlgl[0]=prtlgl[1]=-1*RAND_MAX;
for (h=0; h<hyp; ++h)	{
	y=0;
	for (c1=0; c1<cases; ++c1)	{
		y+=loglike[c1][h];
		if (loglike[c1][h]>mlchyp[c1])	{
			mlcase[c1]=hypos[h];		/* new most likely rate hypothesis	*/
			mlchyp[c1]=loglike[c1][h];
			mlrhyp[c1]=h;
			}
		}
	if (y>prtlgl[1])	{
		prtlgl[1]=y;
		mlmpr=hypos[h];
		}
	}

prtaicc[1]=calc_aic_c(1,prtlgl[1],n);
printf("Rates\tlgN\t\tAICc\t\tPartitions\n");
printf("%5d\t%3.2f\t%3.2f\t\ta = %2.1f: %d-%d\n",1,prtlgl[1],prtaicc[1],mlmpr,1,cases);
/*printf("The most likely overall 1-rate is %3.2f with a log-likelihood of %3.2f.\n",mlmpr,prtlgl[1]);	*/
	
rnkcase=ivector(cases);
for (c1=0; c1<cases; ++c1)	rnkcase[c1]=c1;
rnkcase=sort_incintbydouble(rnkcase,mlcase,cases);

partitionlgl=dmatrix(cases+1,cases+1);		/* partitionlgl[x][y] gives maximum log-likelihoods for case set x…y	*/
partitionhyp=dmatrix(cases+1,cases+1);		/* partitionhyp[x][y] gives hypothesis maximizing log-likelihoods for case set x…y	*/

for (c1=0; c1<cases; ++c1)	{
	r1=rnkcase[c1];						/* r1 is ranked case, based on most likely rate		*/
	partitionlgl[c1][c1]=mlchyp[r1];	/* loglikelihood if partition is solely of this character	*/
	partitionhyp[c1][c1]=((double) (1+mlrhyp[r1]))*dp;
	for (c2=c1+1; c2<cases; ++c2)	{
		r2=rnkcase[c2];					/* r2 is ranked case, based on most likely rate		*/
		/* go from ml rate of first character to ml rate of r2	*/
		partitionlgl[c1][c2]=-1*RAND_MAX;
		/* take h from ml hypothesis for char r1 to ml hypothesis for char r2	*/
		for (h=mlrhyp[r1]; h<=mlrhyp[r2]; ++h)	{
			y=0;
			/* log likelihood of partition at rate is sum of log-likelihoods of r1…r2	*/
			for (c3=c1; c3<=c2; ++c3)	{
				r3=rnkcase[c3];			/* r3 is the 3rd slowest of this rate set	*/
				y+=loglike[r3][h];		/* add loglikelihood of r3 at this rate		*/
				}
			if (y>partitionlgl[c1][c2])	{
				partitionlgl[c1][c2]=y;						/* new maximum log-likelihood	*/
				partitionhyp[c1][c2]=((double) (1+h))*dp;	/* new ml hypothesis			*/
				}
			else	h=1+mlrhyp[r2];		/* kill search: we've gone past the best	*/
			}	/* finish searching range of hypotheses for most likely hyp. for this set	*/
		}	/* 	*/
	}

mlparts[1][0]=parts[1][0]=0;
mlparts[1][1]=parts[1][1]=cases;
/* p gives parameters-1	*/
for (p=2; p<cases; ++p)	{
	prtlgl[p]=-1*RAND_MAX;		/* set likelihood to minimum						*/
	parts[p][0]=0;				/* always start with the first character			*/
	for (c1=1; c1<p; ++c1)
		parts[p][c1]=c1;		/* seed initial array to (say) p[3][0]=0, p[3][1]=1, p[3][2]=2	*/
	parts[p][p]=cases;			/* always end with last character, which is cases-1	*/
	
	/* modify highest partition until no improvement	*/
	for (b=parts[p][p-1]; b<cases; ++b)	{
		y=0.0f;								/* likelihood of particular partition	*/
		for (h=p; h>0; --h)	{
			r1=parts[p][h-1];				/* gets lower ranked character: remember, cell h gives first ranked character of next partition	*/
			r2=parts[p][h]-1;				/* gets higher ranked character: remember, cell h gives first ranked character of next partition	*/
			y+=partitionlgl[r1][r2];		/* add log-likelihood of this partition	*/
			}
		if (y>prtlgl[p])	{
			prtlgl[p]=y;	/* add new greater likelihood			*/
			/* tally new best set of partitions	*/
			for (m=0; m<=p; ++m)	mlparts[p][m]=parts[p][m];
			}
		/* DELETE TRUNCATION: Rebooting for a new (say) 3-rate set will drop likelihood below best	*/
/*		else	b=cases-1;	*/
		/* end this search */
		/* 2010-01-25: START FROM HERE!!!!!	*/
		if (b==cases-1)	{
			/* count backwards and make sure that lower parts can increment still	*/
			/* trick: we need to go from:
				p[0][0]=0;
				p[0][1]=1;
				p[0][2]=17;
				p[0][3]=18	where cases = 17 to
				p[0][0]=0;
				p[0][1]=2;
				p[0][2]=3;
				p[0][3]=18 and get no further than:
				p[0][0]=0;
				p[0][1]=16;		(cases-[partitions-1])
				p[0][2]=17;		(cases-[partitions-2])
				p[0][3]=18.		(cases-[partitions-3])
			/* increment final partition: but if it runs out, then reset it and increment the prior increment	*/
			/* note: this should reset all increments back to the 2nd if allowed to do so.						*/
			/* go from p[3][0]=0, p[3][1]=1, p[3][2]=17, p[3][3]=18 to p[3][0]=0, p[3][1]=2, p[3][2]=3, p[3][3]=18	*/
			m=p-2;
			while(parts[p][m]==cases-(p-m-1) && m>0)	{
				--m;
				}	/* count back to next partition that can be incremented	*/
			if (m>0)	{
				++parts[p][m];
				for(c1=m+1; c1<p; ++c1)	parts[p][c1]=parts[p][c1-1]+1;
				b=parts[p][p-2];			/* restart increments	*/
				}
			}	/* end case of incrementing lower partition	*/
		else
			++parts[p][p-1];
		}	
	prtaicc[p]=calc_aic_c(p,prtlgl[p],n);
	printf("%5d\t%3.2f\t%3.2f\t\t",p,prtlgl[p],prtaicc[p]);
	for (b=0; b<p; ++b)	{
		mlmpr=partitionhyp[mlparts[p][b]][mlparts[p][b+1]-1];		/* this gives the "rate" maximizing partition likelihood	*/
		printf("a = %2.1f: ",mlmpr);
		/* set up dummy array with original character numbers	*/
		for (m=0; m<(mlparts[p][b+1]-mlparts[p][b]); ++m)		{
			prtcse[m]=rnkcase[mlparts[p][b]+m];
			}
		prtcse=ishellsort_inc(prtcse,(mlparts[p][b+1]-mlparts[p][b]));
		for (m=0; m<(mlparts[p][b+1]-mlparts[p][b]); ++m)		{
			printf("%d",prtcse[m]+1);
			if (m<((mlparts[p][b+1]-mlparts[p][b])-1))	{
				if(prtcse[m]==(prtcse[m+1]-1))	{
					printf("-");
					while (prtcse[m]==(prtcse[m+1]-1) && m<(mlparts[p][b+1]-mlparts[p][b]))
						++m;	
					printf("%d",prtcse[m]+1);
					if (m<((mlparts[p][b+1]-mlparts[p][b])-1))
						printf(", ");
				}	/* end loop to print out x-y	*/
				else printf(", ");
			}	/* finish printing original character numbers	*/
		}
		if (b<(p-1))	printf("; ");
		else			printf("\n");
		}

	if (prtaicc[p]>prtaicc[p-1])	p=cases;
	else							++mp;
	}

/* prepare output	*/
product=dmatrix(cases+2,mp+2);
product[cases][0]=mp+1;				/* gives size of the output matrix	*/ 
for (p=1; p<=mp+1; ++p)	{
	product[cases][p]=prtlgl[p];
	product[cases+1][p]=prtaicc[p];
	}
	
for (p=1; p<=mp+1; ++p)	{
	for (b=1; b<=p; ++b)	{
		for (r1=mlparts[p][b-1]; r1<mlparts[p][b]; ++r1)	{
			c1=rnkcase[r1];
			product[c1][p]=partitionhyp[mlparts[p][b-1]][mlparts[p][b]-1];
			}
		}
	}
for (r1=0; r1<cases; ++r1)	{
	c1=rnkcase[r1];				/* get original character number	*/
	product[c1][0]=mlcase[c1];	/* most likely hypothesis for case	*/
//	for (p=1; p<=mp+1; ++p)	{
//		b=0;
//		while (mlparts[p][b]<=r1 && b<p)	{
//			++b;	/* find partition in which r1/c1 fits	*/
//			}
//		product[c1][p]=partitionhyp[mlparts[p][b-1]][mlparts[p][b]-1];
//		}	/* end search for rate for c1 in the p-rate hypothesis	*/
	}
	
free_ivector(rnkcase);
free_ivector(mlrhyp);
free_imatrix(parts,cases+1,cases+1);
free_imatrix(mlparts,cases+1,cases+1);
free_dvector(mlcase);
free_dvector(mlchyp);
free_dvector(prtaicc);
free_dvector(prtlgl);
free_dmatrix(partitionlgl,cases+1,cases+1);
free_dmatrix(partitionhyp,cases+1,cases+1);
return product;
}