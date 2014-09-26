/* Routines to manipulate matrices/arrays already in computer memory.  
/*	Peter Wagner	05/02
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define matrixanalysis
#include "matrixanalysis.h"
#include "matrixchange.h"
#include "matrixreading.h"
#include "memory.h"
#include "minmax.h"
#include "probability.h"
#include "sort.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

/*
Function returning the similarities among characters given mutual compatiblity (see O’Keefe & Wagner 2001 Syst. Biol.)
 	given a compatibility matrix.
Neads:
	comatrix: binary matrix, such as a char x char compatibility matrix, 
		with 1 = compatible and 0 = incompatible
	n: number of characters
*****************************************************************************/
double **gowertrans(double **matrix, int n)
{
int		i, j;
double	ave;
double **dmat;
double	*average;
double	**gower;

dmat=dmatrix(n,n);
for (i=0; i<n; ++i)
	for (j=0; j<n; ++j)
		dmat[i][j]=-0.5*matrix[i][j];
		
gower=dmatrix(n,n);
average=dvector(n);

for (i=0; i<n; ++i)	average[i]=dvectoraver(dmat[i],n);

ave=dmatrixaver(dmat,n,n);

for (i=0; i<n-1; ++i)	{
	for (j=i+1; j<n; ++j)	{
		gower[i][j]=gower[j][i]=matrix[i][j]-(average[i]+average[j])+ave;
		}
	}

free_dvector(average);
free_dmatrix(dmat,n,n);
return gower;
}

/* dvectoraver - calculates the average value in vector v;
/*
/* Requires:
/*		v - an array of real numbers;
/*		n - the length of the array;
/* Returns:
/*		ave - the average value;
***************************************************************************************/
double dvectoraver(double *v, int n)
{
double denom=0.0f, sum=0.0f, ave=0.0f;

denom=n;
sum=sumdvector(v,n);
ave=sum/denom;

return ave;
}

/* dmatrixaver - calculates the average value in vector v;
/*
/* Requires:
/*		m - an array of real numbers;
/*		r - the number of rows;
/*		c - the number of columns;
/* Returns:
/*		ave - the average value;
***************************************************************************************/
double dmatrixaver(double **m, int r, int c)
{
int a, b;
double denom=0.0f, sum=0.0f, ave=0.0f;

for (a=0; a<r; ++a)	{
	for (b=0; b<c; ++b)	{
		sum=sum+m[a][b];
		++denom;
		}
	}

ave=sum/denom;

return ave;
}

/* sumdarray - sums two arrays of real numbers
/*
/* Requires:
/*		i - the array of sums;
/*		j - an array of real numbers;
/*		k - another array of real numbers;
/*		N - the length of the array;
/* Returns:
/*		i - the sum of the elements;
***************************************************************************************/
void *adddvectors(double *i, double *j, int N)
{
int	a;
for (a=0; a<N; ++a)	i[a]+=j[a];
return i;
}


/* sumdarray - sums two arrays of real numbers
/*
/* Requires:
/*		i - the array of sums;
/*		j - an array of real numbers;
/*		k - another array of real numbers;
/*		N - the length of the array;
/* Returns:
/*		i - the sum of the elements;
***************************************************************************************/
double *sumdvectors(double *i, double *j, double *k, int N)
{
int	a;
for (a=0; a<N; ++a)	i[a]=k[a]+j[a];
return i;
}

/* sumivectosr - sums two arrays of real numbers
/*
/* Requires:
/*		i - the array of sums;
/*		j - an array of real numbers;
/*		k - another array of real numbers;
/*		N - the length of the array;
/* Returns:
/*		i - the sum of the elements;
***************************************************************************************/
int *sumivectors(int *i, int *j, int *k, int N)
{
int	a;
for (a=0; a<N; ++a)	i[a]=k[a]+j[a];
return i;
}

/* sumivector - sums one array of integer
/*
/* Requires:
/*		i - the array of long;
/*		n - the length of l;
/* Returns:
/*		s - the sum of the elements;
***************************************************************************************/
int sumivector(int *i, int n)
{
int	a, s=0;
for (a=0; a<n; ++a)	s=s+i[a];
return s;
}

/* sumdvector - sums one array of double
/*
/* Requires:
/*		d - the array of long;
/*		n - the length of l;
/* Returns:
/*		s - the sum of the elements;
***************************************************************************************/
double sumdvector(double *d, int n)
{
int	a;
double s=0.0f;

for (a=0; a<n; ++a)	s=s+d[a];
return s;
}

/* sumlvector - sums one array of long
/*
/* Requires:
/*		l - the array of long;
/*		n - the length of l;
/* Returns:
/*		s - the sum of the elements;
***************************************************************************************/
long sumlvector(long *l, int n)
{
int	a;
long s=0;

for (a=0; a<n; ++a)	s=s+l[a];
return s;
}

/* sumulvector - sums one array of unsigned long
/*
/* Requires:
/*		l - the array of long;
/*		n - the length of l;
/* Returns:
/*		s - the sum of the elements;
***************************************************************************************/
unsigned long sumulvector(unsigned long *l, int n)
{
int	a;
unsigned long s=0;

for (a=0; a<n; ++a)	s=s+l[a];
return s;
}

/* sumlmatrxcol - sums theh column of a long matrix
/*
/* Requires:
/*		l - the matrix of long;
/*		r - the number of rows;
/*		c - the column to be summed;
/* Returns:
/*		s - the sum of the elements;
***************************************************************************************/
long sumlmatrixcol (long **l, int r, int c)
{
int	a;
long s=0;

for (a=0; a<r; ++a)	s=s+l[a][c];
return s;
}

/* sumlmatrxcol - sums theh column of a long matrix
/*
/* Requires:
/*		l - the matrix of long;
/*		r - the number of rows;
/*		c - the column to be summed;
/* Returns:
/*		s - the sum of the elements;
***************************************************************************************/
unsigned long sumulmatrixcol (unsigned long **l, int r, int c)
{
int	a;
unsigned long s=0;

for (a=0; a<r; ++a)	s=s+l[a][c];
return s;
}

/* countivector - counts the number of nonzero elements in an array.  
/*
/* Requires:
/*		i - the array of sums;
/*		j - an array of real numbers;
/*		k - another array of real numbers;
/*		N - the length of the array;
/* Returns:
/*		i - the sum of the elements;
***************************************************************************************/
int countivector(int *i, int n)
{
int	a, s=0;
for (a=0; a<n; ++a)	if (i[a]>0)	++s;
return s;
}

/* darraytotal - returns the sum of values in an array of real numbers
/*
/* Requires:
/*		d - an array of real numbers;
/*		N - the length of the array;
/* Returns:
/*		T - the sum of the elements;
***************************************************************************************/
double darraytotal(double *d, int N)
{
int		i;
double	T=0;

for (i=0; i<N; ++i)	T=T+d[i];
return T;
}

/* intarraytotal - returns the sum of values in an array of integers
/*
/* Requires:
/*		d - an array of integers;
/*		N - the length of the array;
/* Returns:
/*		T - the sum of the elements;
***************************************************************************************/
int iarraytotal(int *d, int N)
{
int		i, T=0;

for (i=0; i<N; ++i)	T=T+d[i];
return T;
}

/* larraytotal - returns the sum of values in an array of integers
/*
/* Requires:
/*		d - an array of integers;
/*		N - the length of the array;
/* Returns:
/*		T - the sum of the elements;
***************************************************************************************/
long larraytotal(long *d, int N)
{
int		i;
long	T=0;

for (i=0; i<N; ++i)	T=T+d[i];
return T;
}

/* charraydifferences - Returns the number of differences between two species
/* Requires:
/*		sp1: species 1's characters
/*		sp2: species 2's characters
/*		nchars: # characters
/*		UNKNOWN: number for unknown characters
/*		INAP: number for inapplicable characters
/*
/* Returns:
/*		diffs: number of characters that differ
***************************************************************************************/
int charraydifferences(long *sp1, long *sp2, int nchars, int UNKNOWN, int INAP)
{
int		ch, diffs=0;

for (ch=0; ch<nchars; ++ch)	{
	while (ch<nchars && ((sp1[ch]==UNKNOWN || sp1[ch]==INAP) || (sp2[ch]==UNKNOWN || sp2[ch]==INAP)))
		++ch;
	if (ch>=nchars)	break;
	if (sp1[ch]!=sp2[ch])	++diffs;
	}

return diffs;
}

/* charraycomparable - Returns the number of characters comparable between two species
/* Requires:
/*		sp1: species 1's characters
/*		sp2: species 2's characters
/*		nchars: # characters
/*		UNKNOWN: number for unknown characters
/*		INAP: number for inapplicable characters
/*
/* Returns:
/*		comp: proportion of characters that differ
***************************************************************************************/
int charraycomparable(long *sp1, long *sp2, int nchars, int UNKNOWN, int INAP)
{
int		ch, comp=0;

for (ch=0; ch<nchars; ++ch)	{
	while (ch<nchars && ((sp1[ch]==UNKNOWN || sp1[ch]==INAP) || (sp2[ch]==UNKNOWN || sp2[ch]==INAP)))
		++ch;
	if (ch>=nchars)	break;
	++comp;
	}

return comp;
}

/* ihistogram - Returns array h in which h[x] gives the number of times x is observed in array v
/* Requires:
/*		v: data array
/*		n: length of data array
/*
/* Returns:
/*		histo: array where histo[x] gives the number of times x occurs in v
/* NOTE: returned vector is long

/* altered 2012-12-29: for some reason, I was setting it up to not allocate enough space if there were not entities observed only 1 time
***************************************************************************************/
long *ilhistogram (int *v, int n)
{
int i, mx;
//int	mn;
long *histo;

mx=maxiarray(v,n);
//mn=miniarray(v,n);
//histo=lvector((mx-mn)+2);
histo=lvector(mx+2);

for (i=0; i<n; ++i)	++histo[v[i]];

return histo;
}

/* lhistogram - Returns array h in which h[x] gives the number of times x is observed in array v
/* Requires:
/*		v: data array
/*		n: length of data array
/*
/* Returns:
/*		histo: array where histo[x] gives the number of times x occurs in v
/* NOTE: returned vector is long
***************************************************************************************/
long *llhistogram (long *v, int n)
{
int i, mn, mx;
long *histo;

mx=maxlarray(v,n);
mn=minlarray(v,n);
histo=lvector((mx-mn)+2);

for (i=0; i<n; ++i)	++histo[v[i]];

return histo;
}

/* ulhistogram - Returns array h in which h[x] gives the number of times x is observed in array v
/* Requires:
/*		v: data array
/*		n: length of data array
/*
/* Returns:
/*		histo: array where histo[x] gives the number of times x occurs in v
/* NOTE: returned vector is long
***************************************************************************************/
unsigned long *ulhistogram (long *v, int n)
{
int i, mn, mx;
unsigned long *histo;

mx=maxlarray(v,n);
mn=minlarray(v,n);
histo=ulvector((mx-mn)+2);

for (i=0; i<n; ++i)	++histo[v[i]];

return histo;
}

/* uihistogram - Returns array h in which h[x] gives the number of times x is observed in array v
/* Requires:
/*		v: data array
/*		n: length of data array
/*
/* Returns:
/*		histo: array where histo[x] gives the number of times x occurs in v
/* NOTE: returned vector is long
***************************************************************************************/
unsigned long *uihistogram (int *v, int n)
{
int i, mn, mx;
unsigned long *histo;

mx=maxiarray(v,n);
mn=miniarray(v,n);
histo=ulvector((mx-mn)+2);

for (i=0; i<n; ++i)	++histo[v[i]];

return histo;
}


/* ihistogram - Returns array h in which h[x] gives the number of times x is observed in array v
/* Requires:
/*		v: data array
/*		n: length of data array
/*
/* Returns:
/*		histo: array where histo[x] gives the number of times x occurs in v
/* NOTE: returned vector is double
***************************************************************************************/
double *idhistogram (int *v, int n)
{
int i, mn, mx;
double *histo;

mx=maxiarray(v,n);
mn=miniarray(v,n);
histo=dvector(mx+1);

for (i=0; i<n; ++i)	++histo[v[i]];

return histo;
}

/* ihistogramfull - Returns array h in which h[x] gives the number of times x is observed in array v with zeros between max and possible max
/* Requires:
/*		v: data array
/*		n: length of data array
/*		k: maximum possible length
/*
/* Returns:
/*		histo: array where histo[x] gives the number of times x occurs in v
/* NOTE: returned vector is double
***************************************************************************************/
double *idhistogramfull (int *v, int n, int k)
{
int i;
double *histo;

histo=dvector(k+1);
cleardvector(histo,k+1,0.0f);

for (i=0; i<n; ++i)	++histo[v[i]];

return histo;
}

/* lhistogram - Returns array h in which h[x] gives the number of times x is observed in array v
/* Requires:
/*		v: data array
/*		n: length of data array
/*
/* Returns:
/*		histo: array where histo[x] gives the number of times x occurs in v
/* NOTE: returned vector is double
***************************************************************************************/
double *ldhistogram (long *v, int n)
{
int i, mn, mx;
double *histo;

mx=maxlarray(v,n);
mn=minlarray(v,n);
histo=dvector(2*((mx-mn)+2));

for (i=0; i<n; ++i)	++histo[v[i]];

return histo;
}


/* codedcombinations - Returns a long matrix giving observed morphotypes (combinations) and a list of species with that morphotype
/*	Requires:
/*		matrix: character matrix
/*		notu:	number of taxa
/*		nchar:	number of characters
/*		miss:	unknown score
/*		gap:	inapplicable score
/* Returns:
/*		results:
/*			results[0]: combination
/*			results[1]: number of species with combination
/*			results[2…n]: species number
***************************************************************************************/
long **codedcombinations (long **matrix, int notu, int nchar, int miss, int gap)
{
int		b, c, f, sp, t, u;
int 	mxtx=1, unique=0;
/*int		*mxst, *combos;	*/
long	**results;

/*states=numberstates(matrix, notu, nchar, miss, gap);
/*mxst=maxcharstates(matrix, notu, nchar, miss, gap);	/* maximum state values */
/*combos=ivector(nchar);	/* the number of combinations for any state of that character with subsequent characters */
/*for (a=0; a<nchar; ++a)	{
/*	if (a<(nchar-1))	combos[a]=mxst[a+1];
/*	else				combos[a]=1;
/*	for (b=a+2; b<nchar; ++b)
/*		combos[a]=combos[a]*(mxst[b]+1);						/* use mxst[b]+1 to allow for 0…mxst */
/*	}												*/

mxtx=maxduplicates(matrix,notu,nchar);
unique=countunique(matrix,notu,nchar);

results=lmatrix(unique,mxtx+3);	/* results[a][0]: morphotype; results[a][1]: richness; results[a][2…mxtx+1]: taxon #'s */

/*mnst=mincharstates(matrix, notu, nchar, miss, gap);*/

u=0;
/* reduce tab-delimited matrix to a single 'string' of numbers */ 
for (sp=0; sp<notu; ++sp)	{

	c=matrix[sp][0];
	for (b=1; b<nchar; ++b)	{
		c*=10;
		c+=matrix[sp][b];
		}

	f=-1;
	for (b=0; b<u && f==-1; ++b)
		if (results[b][0]==c)	f=b;

	/* if f is not -1, then sp has been found before */
	if (f>=0)	{
		++results[f][1];
		t=1+results[f][1];
		if (t<0 || t>(mxtx+1))	printf("We have an array element %d when the maximum is %d\n",t,mxtx+1);
		results[f][results[f][1]+1]=sp;
		}
	
	/* this marks the first time a morphotype is found */
	else	{
		results[u][0]=matrix[sp][0];
		for (b=1; b<nchar; ++b)	{
			results[u][0]*=10;
			results[u][0]+=matrix[sp][b];
			}	/* finish condensing morphotype to a single numberr */
		results[u][1]=1;
		results[u][results[u][1]+1]=sp;
		++u;
		if (u>unique || u<0)	printf("We have an array element %d when the maximum is %d\n",u,unique-1);
		}	/* end case */
	}	/* end search of species */

/*free_ivector(combos);
/*free_ivector(mxst);	*/

return results;
}

/* maxduplicates - finds the maximum number of times a particular coding appears in a matrix
/*	Requires:
/*		matrix: character matrix
/*		notu:	number of taxa
/*		nchar:	number of characters
/*  Returns:
/*		mxtx:	the greatest number of duplicates
***************************************************************************************/
int maxduplicates(long **matrix, int notu, int nchar)
{
int		a, c, s1, s2;
int 	mx=0, mxtx=1;

for (s1=0; s1<notu-mxtx; ++s1)	{
	mx=1;
	for (s2=s1+1; s2<notu; ++s2)	{
		a=1;
		for (c=0; (c<nchar && a==1); ++c)	{
			if (matrix[s1][c]!=matrix[s2][c])	a=0;
			}	/* returns 1 if identical, 0 if different */
		mx+=a;
		}	/* end search of subsequent taxa for identical codings */
	if (mx>mxtx)	mxtx=mx;
	}	/* end search of taxon for duplicates */
return mxtx;
}

/* maxduplicates - finds the maximum number of times a particular coding appears in a matrix
/*	Requires:
/*		matrix: character matrix
/*		notu:	number of taxa
/*		nchar:	number of characters
/*		miss:	unknown score
/*		gap:	inapplicable score
/*  Returns:
/*		unique:	the greatest number of duplicates
***************************************************************************************/
int countunique(long **matrix, int notu, int nchar)
{
int		a, c, s1, s2;
int 	unique=1;

for (s1=1; s1<notu; ++s1)	{
	a=0;	/* a=0 means that a match has not yet been found */
	for (s2=s1-1; s2>=0 && a==0; --s2)	{
		a=1;	/* a stays at 1 if the taxa are identical */
		for (c=0; (c<nchar && a==1); ++c)	{
			if (matrix[s1][c]!=matrix[s2][c])	a=0;
			}	/* return 1 if identical, 0 if different */
		if (a==1)	s2=0;
		}	/* finish searching prior taxa for identical scores */
	if (a==0)	++unique;
	}	/* finish searching taxa for unique codings */
return unique;
}

int countuniquechars(long **matrix, int notu, int nchar)
{
int	c1, c2, s;
int match=0, count;

count=nchar;			/* number of uniquely coded characters	*/

for (c1=(nchar-1); c1>0; --c1)	{
	for (c2=(c1-1); c2>=0; --c2)	{
		match=1;	/* fixed 2014-09-15	*/
		for (s=0; s<notu && match==1; ++s)	{
			if (matrix[s][c1]!=matrix[s][c2])	match=0;
			}	/* compare all species unti a mismatch is found	*/
		if (match==1)	{
			--count;		/* reduce novel character count	*/
			c2=0;
			}	/* designate character as matching another with a 0	*/
		}	/* end search of other characters	*/
	}
return count;
}

int *countcharreplicates(long **matrix, int notu, int nchar, int annchar)
{
int	c1, c2, s;
int match=0, count;
int *novel, *charreps;

novel=ivector(nchar);	/* number of characters matching this one (0 if duplicate of earlier character(	*/
count=nchar;			/* number of uniquely coded characters	*/

novel[0]=1;
for (c1=(nchar-1); c1>0; --c1)	{
	++novel[c1];
	for (c2=(c1-1); c2>=0; --c2)	{
		match=1;	/* fixed 2014-09-15	*/
		for (s=0; s<notu && match==1; ++s)	{
			if (matrix[s][c1]!=matrix[s][c2])	match=0;
			}	/* compare all species unti a mismatch is found	*/
		if (match==1)	{
			novel[c2]+=novel[c1];	/* this will find total matching characters	*/
			--count;		/* reduce novel character count	*/
			c2=novel[c1]=0;
			}	/* designate character as matching another with a 0	*/
		}	/* end search of other characters	*/
	}

charreps=ivector(count);
/* now condense vector by removing zeros	*/
c1=0;
for (c2=0; c2<nchar; ++c2)	{
	if (novel[c2]>0)	{
		charreps[c1]=novel[c2];
		++c1;
		if (c1==annchar)	c2=nchar;
		}
	}

free_ivector(novel);
return charreps;
}


/* avedvector - finds the average value of a vector
/*	Requires:
/*		v: a vector of real numbers
/*		n: number of values
/*  Returns:
/*		x: the mean
***************************************************************************************/
double avedvector(double *v, int n)
{
int a;
double x=0.0f;

for (a=0; a<n; ++a)	
	x+=v[a];

x/=((double) n);

return x;
}


/* multiplydmatrix: multiplies an r x c and a c x r matrix
	A: r x c matrix of double
	B: c x r matrix of double
	r: rows in A, columns in B
	c: columns in B, rows in A
*******************/
void **multiplydmatrix(double **A, double **B, double **C, int r, int c)
{
int	i, j, k;

for (i=0; i<r; ++i)	{
	for (j=0; j<r; ++j)	{
		C[i][j]=0;
		for (k=0; k<c; ++k)	C[i][j]+=A[i][k]*B[k][j];
		}
	}
}

/* mexp: calculates P=e^tQ = I + ∑((tQ)^n)/n!
	Q: k x k rate matrix
	t: duration over which t occurs
	k: size of matrix
 *******************/
double **mexp(double **Q, int k)
{
int	a, b, i;
double	x=1.0f;
double **P, **Z, **A;

A=dmatrix(k,k);
Z=dmatrix(k,k);
P=dmatrix(k,k);

cleardmatrix(P, k, k, 0);
for (a=0; a<k; ++a)	{
	for (b=0; b<k; ++b)	{
		if (b==a)	P[a][b]=1+Q[a][b];
		else		P[a][b]=Q[a][b];
		}
	}
equaldmatrix(Z, Q, k, k);

for (i=2; i<10; ++i)	{
	x*=((double) i);
	multiplydmatrix(Q,Z,A,k,k);
	for (a=0; a<k; ++a)
		for (b=0; b<k; ++b)
			P[a][b]+=(A[a][b]/x);	
	equaldmatrix(Z,A,k,k);
	}

//for (a=0; a<k; ++a)	P[a][a]+=1.0f;	/* add identity matrix to P matrix	*/

free_dmatrix(Z,k,k);
free_dmatrix(A,k,k);

return P;
}

/* rank elements of a unsigned long vector; returns double because of ties	*/
double *rankivector_dec(int *v, int n)
{
int a, b, c;
int *tv;
double	*ranks;

tv=ivector(n);
equalivector(tv,v,n);
ishellsort_dec_command(tv,n);

ranks=dvector(n);
ranks[n-1]=n;		/* rank the last one for now: it might get erased	*/

for (a=0; a<(n-1); ++a)	{
	b=a;
	while (tv[b+1]==tv[a] && b<n)	++b;
	for (c=a; c<=b; ++c)	ranks[c]=1+((double) (a+b))/2;
	a=b;	/* slide a along	*/
	}

free_ivector(tv);
return ranks;
}


/* rank elements of a unsigned long vector; returns double because of ties	*/
double *rankulvector_dec(unsigned long *v, int n)
{
int a, b, c;
unsigned long *tv;
double	*ranks;

tv=ulvector(n);
equalulvector(tv,v,n);
ulshellsort_dec_command(tv,n);

ranks=dvector(n);
ranks[n-1]=n;		/* rank the last one for now: it might get erased	*/

for (a=0; a<(n-1); ++a)	{
	b=a;
	while (tv[b+1]==tv[a] && b<n)	++b;
	for (c=a; c<=b; ++c)	ranks[c]=1+((double) (a+b))/2;
	a=b;	/* slide a along	*/
	}

free_ulvector(tv);
return ranks;
}


/* rank elements of a unsigned long vector; returns double because of ties	*/
double *rankdvector_dec(double *v, int n)
{
int a, b, c;
double *tv;
double	*ranks;

tv=dvector(n);
equaldvector(tv,v,n);
dshellsort_dec_command(tv,n);

ranks=dvector(n);
ranks[n-1]=n;		/* rank the last one for now: it might get erased	*/

for (a=0; a<(n-1); ++a)	{
	b=a;
	while (tv[b+1]==tv[a] && b<n)	++b;
	for (c=a; c<=b; ++c)	ranks[c]=1+((double) (a+b))/2;
	a=b;	/* slide a along	*/
	}

free_dvector(tv);
return ranks;
}

/* rank elements of a unsigned long vector; returns double because of ties	*/
double *rankielements_dec(int *v, int n)
{
int a, b, c;
int *tv;
double	*dummy, *ranks;

tv=ivector(n);
equalivector(tv,v,n);
ishellsort_dec_command(tv,n);

dummy=dvector(n);
dummy[n-1]=n;		/* rank the last one for now: it might get erased	*/

for (a=0; a<(n-1); ++a)	{
	b=a;
//	while (tv[b+1]==tv[a] && b<n)	++b;
	while (tv[b+1]==tv[a] && b<(n-1))	++b;
	for (c=a; c<=b; ++c)	dummy[c]=1+((double) (a+b))/2;
	a=b;	/* slide a along	*/
	}

ranks=dvector(n);

for (a=0; a<n; ++a)	{
	b=0;
	while (tv[b]>v[a] && b<n)	++b;
	ranks[a]=dummy[b];
	}

free_ivector(tv);
free_dvector(dummy);
return ranks;
}

/* rank elements of a unsigned long vector; returns double because of ties	*/
double *rankdelements_dec(double *v, int n)
{
int a, b, c;
double *tv;
double	*dummy, *ranks;

tv=dvector(n);
equaldvector(tv,v,n);
dshellsort_dec_command(tv,n);

dummy=dvector(n);
dummy[n-1]=n;		/* rank the last one for now: it might get erased	*/

for (a=0; a<(n-1); ++a)	{
	b=a;
//	while (tv[b+1]==tv[a] && b<n)	++b;
	while (tv[b+1]==tv[a] && b<(n-1))	++b;
	for (c=a; c<=b; ++c)	dummy[c]=1+((double) (a+b))/2;
	a=b;	/* slide a along	*/
	}

ranks=dvector(n);

for (a=0; a<n; ++a)	{
	b=0;
	while (tv[b]>v[a] && b<n)	++b;
	ranks[a]=dummy[b];
	}

free_dvector(tv);
free_dvector(dummy);
return ranks;
}

