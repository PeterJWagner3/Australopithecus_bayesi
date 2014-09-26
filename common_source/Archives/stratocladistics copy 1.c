/* Routines to date nodes and calculate the gaps implied by trees     */
/*		Written by Peter Wagner 01/95
/*		Updated:	03/96
/*					04/96
/*					02/97
/*					09/97
/*					08/98
/*					09/02
/*					01/03
/*					06/03 - interger arrays converted to long for dates
/*					03/05 - stratocompatibility added
/*					09/05 - more stratocompatibility added
/**********************************************************************/

#define historicaldiversity
#include "historicaldiversity.h"
#define matrixanalysis
#include "matrixanalysis.h"
#define matrixchange
#include "matrixchange.h"
#define matrixreading
#include "matrixreading.h"
#define memory
#include "memory.h"
#define minmax
#include "minmax.h"
#define Optimization
#include "optimization.h"
#define sort
#include "sort.h"
#define stratocladistics
#include "stratocladistics.h"
#define TreeRead
#include "tree_read.h"

//* dateclade - Finds age of clades based on the oldest member in the clade
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	fa - First Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	F - modified array fa, with elements notu…notu+clades now filled.
****************************************************************************/
long *dateclade(long *fa, long **tree, int clades, int notu)
{
int	node, b, sp, htu;
long *F;

F=lvector(notu+clades);
b=maxlarray(fa,notu);

for (sp=0; sp<notu; ++sp)				F[sp]=fa[sp];

for (node=clades-1; node>=0; --node)	{
	F[htu=node+notu]=F[sp=tree[node][1]];
	for (b=2; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (F[htu]>F[sp])	F[htu]=F[sp];
		}
	}
	
return F;
}

//* dateclade - Finds age of clades based on the oldest member in the clade
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	fa - First Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	F - modified array fa, with elements notu…notu+clades now filled.
****************************************************************************/
double *datecladereal(double *fa, long **tree, int clades, int notu)
{
int	node, b, sp, htu;
double *F;

F=dvector(notu+clades);
b=maxdarray(fa,notu);

for (sp=0; sp<notu; ++sp)				F[sp]=fa[sp];

for (node=clades-1; node>=0; --node)	{
	F[htu=node+notu]=F[sp=tree[node][1]];
	for (b=2; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (F[htu]>F[sp])	F[htu]=F[sp];
		}
	}
	
return F;
}


//* datecladeFile - Finds age of clades based on the oldest member in the clade and an fa file
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	fa - First Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	fa - modified array fa, with elements notu…notu+clades now filled.
****************************************************************************/
long* datecladefile(char *Origins, long **tree, int clades, int notu)
{
int	 node, b, sp, htu;
long *fa;
FILE *fopen();
FILE *FKA;

fa=lvector(notu+clades);
FKA = fopen(Origins,"r");
for (b=0; b<notu; ++b)
	fscanf(FKA,"%i",&fa[b]);
fclose(FKA);

b=maxlarray(fa,notu);

for (node=clades-1; node>=0; --node)	{
	fa[htu=node+notu]=fa[sp=tree[node][1]];
	for (b=2; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (fa[htu]>fa[sp])	fa[htu]=fa[sp];
		}
	}
	
return fa;
}


//* dateclade - Finds age of clades based on the oldest member in the clade
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	fa - First Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	F - modified array fa, with elements notu…notu+clades now filled.
****************************************************************************/
double *datecladerealfile(char *Origins, long **tree, int clades, int notu)
{
int	 node, b, sp, htu;
double *fa;
FILE *fopen();
FILE *FKA;

fa=dvector(notu+clades);
FKA = fopen(Origins,"r");
for (b=0; b<notu; ++b)
	fscanf(FKA,"%lf",&fa[b]);
fclose(FKA);

for (node=clades-1; node>=0; --node)	{
	fa[htu=node+notu]=fa[sp=tree[node][1]];
	for (b=2; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (fa[htu]>fa[sp])	fa[htu]=fa[sp];
		}
	}
	
return fa;
}

//* datecladeAdd - Finds minimum divergence dates assuming monophyletic groups
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	fa - First Appearances (these are modified along the way)
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
****************************************************************************/
void datecladeadd(long *fa, long **tree, int clades, int notu)
{
int	node, b, sp;

for (node=clades-1; node>=0; --node)	{
	fa[node+notu]=100000;
	for (b=1; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (fa[node+notu]>fa[sp])	fa[node+notu]=fa[sp];
		}
	}
}

//* datecladeExtra - Finds minimum divergence dates allowing for ancestors
/* Requires:
/*	Diverg - A matrix of divergence times between taxa
/*	tree - a matrix containing the tree structure;
/*	div2 - diversity of each clade;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	fa - modified array fa, with elements notu…notu+clades now filled.
****************************************************************************/
long *datecladeextra(int **Diverg, long *fa, long **tree, int clades, int notu)
{
int		a, node, sp1, sp2, both;
long	**VennTree;	/* tree listing all taxa in the node */

VennTree=lmatrix(notu, notu);
for (node=0; node<clades; ++node)
	for (a=0; a<=tree[node][0]; ++a)	VennTree[node][a]=tree[node][a];
VennTree=clademember(VennTree, notu, clades);

for (sp1=0; sp1<notu-1; ++sp1)	{
	for (sp2=sp1+1; sp2<notu; ++sp2)	{
		if (Diverg[sp1][sp2]!=0)	{
			/* work down the tree */
			/* the highest node with both sp1 and sp2 now gets fa[node+notu]=Diverg[sp1][sp2]	*/
			both=0;
			for (node=clades-1; both<2; --node)	{
				both=0;
				/* increment both whenever there is a match 	*/
				/* if it equals two, then there are two matches */
				for (a=1; a<=VennTree[node][0]; ++a)	{
					if (VennTree[node][a]==sp1 || VennTree[node][a]==sp2)	++both;
					if (both==2)	{
						a=VennTree[node][0];
						if (fa[node+notu]>Diverg[sp1][sp2])
							fa[node+notu]=Diverg[sp1][sp2];
						}	/* register fa for node, and end search */
					}	/* end search of node */
				}
			}	/* only bother if we have a divergence date of note */
		}
	}
	
free_lmatrix(VennTree, notu, notu);
return fa;
}

long *datecladediverge(long **DT, long *fa, long **tree, int clades, int notu)
{
int sp1, sp2, max, anc;
int **LCA;

max=DT[0][0];
LCA=lastcommonancestormatrix(tree,clades,notu);

for (sp1=0; sp1<notu+clades; ++sp1)	fa[sp1]=max;
for (sp1=0; sp1<notu-1; ++sp1)	{
	for (sp2=sp1+1; sp2<notu; ++sp2)	{
		while (LCA[sp1][sp2]==-1 && sp2<notu)	++sp2;
		if (sp2>=notu)	break;
		anc=LCA[sp1][sp2]+notu;
		if (fa[anc]>DT[sp1][sp2])	fa[anc]=DT[sp1][sp2];
		}
	}
free_imatrix(LCA,notu,notu);
return fa;
}


//* BranchAgeBin - Finds age of branches using a discrete scale assuming no ancestors
/*		NOTE: if you have continuous time, use "BranchAgeCont"
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - First Appearances
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
****************************************************************************/
long *blbinnaive(long **tree, long *fa, int clades, int notu)
{
int	 a, b, sp, anc;
long *bt;

bt=lvector(notu+clades);

for (a=0; a<clades; ++a)	{
	anc=a+notu;
	for (b=1; b<=tree[a][0]; ++b)	{
		sp=tree[a][b];
		bt[sp]=fa[sp]-fa[anc];
		}
	}
	
return bt;
}

//* BranchAgeBin - Finds age of branches using a continuous scale assuming no ancestors
/*		NOTE: if you have time bins (e.g., stages), use "BranchAgeBin"
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - First Appearances
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
****************************************************************************/
double *blmanaive(long **tree, double *fa, int clades, int notu)
{
int	 a, b, sp, anc;
double	 *bt;

bt=dvector(notu+clades);

for (a=0; a<clades; ++a)	{
	anc=a+notu;
	for (b=1; b<=tree[a][0]; ++b)	{
		sp=tree[a][b];
		bt[sp]=fa[sp]-fa[anc];
		}
	}
	
return bt;
}

//* blma - Tallies branch lengths in real time allowing for ancestors.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	la - an array of last appearance dates;
/*	bl - an array of branch lengths;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	bt - range extension on each branch.
****************************************************************************/
double *branchlngma(long **tree, double **ranges, int *bl, int notu)
{
int d, nd, sp, clades, anc;
double *bt;
double *F, *fa;

fa=dvector(notu);
for (d=0; d<notu; ++d)	fa[d]=ranges[d][0];

clades=cladecountbytaxa(tree,notu);
bt=dvector(notu+clades);

F=datecladereal(fa,tree,clades,notu);

for (nd=clades-1; nd>=0; --nd)	{
	anc=-1;
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		if (sp<notu && bl[sp]==0)	{
			anc=sp;
			d=tree[nd][0];
			}
		}
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		/* if ancestor present and ranges overlap or abut, then no gap */
		if (anc>-1)	{
			/* if descendent precedes ancestor, then there is a gap    */
			if (F[sp]<F[anc])		bt[sp]=F[anc]-F[sp];
			/* if ancestor disappears before desc appears, then a gap  */
			else if (F[sp]>ranges[anc][1])	bt[sp]=(F[sp]-ranges[anc][1])-1;
			}
		else
			bt[sp]=F[sp]-F[notu+nd];
		}
	}
free_dvector(F);
free_dvector(fa);
return bt;
}


//* CladeSurvive - Determines how long clades last
/* Requires:
/*	la - Last Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	L - modified array la, with elements notu…notu+clades now filled.
****************************************************************************/
long *CladeSurvive(long *la, long **tree, int clades, int notu)
{
int node, b, sp;
long *L;

L=lvector(notu+clades);

for (sp=0; sp<notu; ++sp)				L[sp]=la[sp];
for (sp=notu; sp<(notu+clades); ++sp)	L[sp]=0;

for (node=clades-1; node>=0; --node)	{
	for (b=1; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (L[node+notu]<L[sp])	L[node+notu]=L[sp];
		}
	}

return L;
}
//* CladeSurviveAdd - Finds minimum divergence dates assuming monophyletic groups
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	la - Last Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	L - modified array la, with elements notu…notu+clades now filled.
****************************************************************************/
void CladeSurviveAdd(long *la, long **tree, int clades, int notu)
{
int node, b, sp;

for (node=clades-1; node>=0; --node)	{
	for (b=1; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (la[node+notu]<la[sp])	la[node+notu]=la[sp];
		}
	}
}

//* temporalbranchlength - Routine to get basic branch lengths in time units.  
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	bt - modified array fa, with elements notu…notu+clades now filled.
****************************************************************************/
double *temporalbranchlength(long **tree, double *fa, int clades, int notu)
{
int	cl, sp, sc;
double	*bt;

bt=dvector(notu+clades);

for (cl=clades-1; cl>=0; --cl)	{
	for (sc=1; sc<=tree[cl][0]; ++sc)	{
		sp=tree[cl][sc];

		bt[sp]=fa[sp]-fa[cl+notu];
		if (bt[sp]<0)	bt[sp]*=-1;

		}
	}

return bt;
}

double *temporalbranchlengthApo(long **tree, long *fa, int *bl, int clades, int notu)
{
int	cl, sp, sc;
int	*anc;
double	*bt;

anc=ivector(clades);
bt=dvector(notu+clades);

for (cl=clades-1; cl>=0; --cl)	{
	for (sc=1; sc<=tree[cl][0]; ++sc)	{
		sp=tree[cl][sc];
		if (bl[sp]==0 && sp<notu)	anc[cl]=1;
		}
	}

for (cl=clades-1; cl>=0; --cl)	{
/* if none of the sister taxa are ancestral, then it is the same as temporalbranchlength */
	if (anc[cl]==0)	{
		for (sc=1; sc<=tree[cl][0]; ++sc)	{
			sp=tree[cl][sc];
			if (sp<notu)
				bt[sp]=0-fa[cl+notu];
			else
				bt[sp]=fa[sp]-fa[cl+notu];
			}
		}
	/* if species is ancestral, then it can acquire apomorphies from the clade's
	origin until the derived taxa first appear */
	else	{
		for (sc=1; sc<=tree[cl][0]; ++sc)	{
			sp=tree[cl][sc];
			if (sp<notu && bl[sp]==0)
				bt[sp]=0;
			else if (sp<notu && bl[sp]>0)
				bt[sp]=0-fa[sp];
			}
		}
	}

return bt;
}
/*
double *AdjustTemporalBLAp(long **tree, long *fa, int *bl, double *bt,double rate, int clades, int notu)
{
int	sp;
double	t;
double *Clock;

Clock=dvector(clades+notu);

for (sp=0; sp<clades+notu; ++sp)
	Clock[sp]=(t=bl[sp])/rate;


free_dvector(Clock);
}	*/

//* naivestratdebt - Tallies stratigraphic gaps assuming no ancestors.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	sd - Number of gaps (in integers - use naiveMIG for real numbers).
****************************************************************************/
int naivestratdebt(long **tree, long *fa, int notu)
{
int d, nd, sp, clades;
int	sd=0;
long *F;

clades=cladecountbytaxa(tree,notu);
F=dateclade(fa,tree,clades,notu);

for (nd=0; nd<clades; ++nd)	{
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		sd=sd+(F[sp]-F[nd+notu]);
		}
	}
return sd;
}

//* stratdebt - Tallies stratigraphic gaps in bins (e.g., stages) allowing for ancestors.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	la - an array of last appearance dates;
/*	bl - an array of branch lengths;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	sd - Number of gaps (in integers - use MIGcalc for real numbers).
****************************************************************************/
int stratdebt(long **tree, long**ranges, int *bl, int notu)
{
int d, nd, sp, clades, anc;
int	sd=0;
long *F, *fa;

fa=lvector(notu);
for (d=0; d<notu; ++d)	fa[d]=ranges[d][0];

clades=cladecountbytaxa(tree,notu);
F=dateclade(fa,tree,clades,notu);

for (nd=clades-1; nd>=0; --nd)	{
	anc=-1;
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		if (sp<notu && bl[sp]==0)	{
			anc=sp;
			d=tree[nd][0];
			}
		}
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		/* if ancestor present and ranges overlap or abut, then no gap */
		if (anc>-1)	{
			/* if descendent precedes ancestor, then there is a gap    */
			if (F[sp]<F[anc])		sd=sd+(F[anc]-F[sp]);
			/* if ancestor disappears before desc appears, then a gap  */
			else if (F[sp]>ranges[anc][1])	sd=sd+((F[sp]-ranges[anc][1])-1);
			}
		else
			sd=sd+(F[sp]-F[notu+nd]);
		}
	}
free_lvector(F);
free_lvector(fa);
return sd;
}

//* stratdebtreal - Tallies stratigraphic gaps in real time allowing for ancestors.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	la - an array of last appearance dates;
/*	bl - an array of branch lengths;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	sd - Number of gaps (in integers - use MIGcalc for real numbers).
****************************************************************************/
double stratdebtreal(long **tree, double **ranges, int *bl, int notu)
{
int d, nd, sp, clades, anc;
double	sd=0.0;
double *F, *fa;

fa=dvector(notu);
for (d=0; d<notu; ++d)	fa[d]=ranges[d][0];

clades=cladecountbytaxa(tree,notu);
F=datecladereal(fa,tree,clades,notu);

for (nd=clades-1; nd>=0; --nd)	{
	anc=-1;
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		if (sp<notu && bl[sp]==0)	{
			anc=sp;
			d=tree[nd][0];
			}
		}
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		/* if ancestor present and ranges overlap or abut, then no gap */
		if (anc>-1)	{
			/* if descendent precedes ancestor, then there is a gap    */
			if (F[sp]<F[anc])		sd=sd+(F[anc]-F[sp]);
			/* if ancestor disappears before desc appears, then a gap  */
			else if (F[sp]>ranges[anc][1])	sd=sd+((F[sp]-ranges[anc][1])-1);
			}
		else
			sd=sd+(F[sp]-F[notu+nd]);
		}
	}
free_dvector(F);
free_dvector(fa);
return sd;
}

//* SCI - Calculates Stratigraphic Consistency Index after Huelsenbeck 1994.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	SCI - consistent nodes / total nodes.
****************************************************************************/
double calcSCI(long **tree, long *fa, int notu)
{
int	d, dcl, sp, nd, d1, d2;
int	clades;
int	clfirst, spfirst, sim;
long *F;
double	con=0.0f, incon=0.0f, SCI=0.0f;

clades=cladecountbytaxa(tree,notu);

F=dateclade(fa,tree,clades,notu);

for (nd=clades-1; nd>=0; --nd)	{
	dcl=0;
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		if (sp>notu)	++dcl;
		}
	if (dcl>0)	{
		/* if it is a bifurcation, then its easy */
		if (tree[nd][0]==2)	{
			if (dcl==2)	{
				d1=tree[nd][1];
				d2=tree[nd][2];
				if (F[d1]==F[d2])	con=con+2.0;
				else	{
					con=con+1.0;
					incon=incon+1.0;
					}
				}
			else if (dcl==1)	{
				if (tree[nd][1]<notu)	{
					d1=tree[nd][1];
					d2=tree[nd][2];
					}
				else	{
					d1=tree[nd][2];
					d2=tree[nd][1];
					}
				if (F[d1]<=F[d2])	con=con+1;
				else				incon=incon+1;
				}
			}	/* end routine for bifurcation */

		/* polytomies are trickier */
		/* if a clade is part of a polytomy and there is any taxon as old	*/
		/*		as or older than that clade, then it is consistent			*/
		/*		NOTE: only one clade can be inconsistent; 			 		*/
		else	{
			clfirst=spfirst=0;
			for (d=1; d<=tree[nd][0]; ++d)
				spfirst=clfirst=clfirst+F[tree[nd][0]];
			
			for (d=1; d<=tree[nd][0]; ++d)	{
				sp=tree[nd][d];
				/* keep track of which TU's appear first */
				if (sp<notu && F[sp]<spfirst)	{
					spfirst=F[sp];
					}
				else if (sp>=notu)	{
					if (F[sp]<clfirst)	{
						clfirst=F[sp];
						sim=1;
						}
					else if (F[sp]==clfirst)	++sim;
					}
				}
			if (clfirst<spfirst && sim==1)	{
				con=con+(dcl-1);			/* at most, one node attached to a polytomy can be inconsistent */
				++incon;				
				}
			else	con=con+dcl;
			}
		}
	}

SCI=con/(con+incon);

free_lvector(F);
return SCI;
}

//* SCIm - Calculates modified Stratigraphic Consistency Index after Habib & Gittelman.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	SCI - consistent nodes / total nodes.
****************************************************************************/
double calcSCIm(long **tree, long *fa, int notu)
{
int	b, c, cl, d, dcl, nd, scl, mxdc;
int	clades, *f1, *dc;
long *F, *fn, **treecl;
double	denom=0.0f, SCIm=0.0f;

clades=cladecountbytaxa(tree,notu);
treecl=cladesinclades(tree,clades,notu);
dc=ivector(clades);
mxdc=0;
for (c=0; c<clades; ++c)	{
	dc[c]=0;
	for (d=1; d<=tree[c][0]; ++d)	if (tree[c][d]>=notu)	++dc[c];
	if (dc[c]>mxdc)	mxdc=dc[c];
	}
f1=ivector(mxdc);

F=dateclade(fa,tree,clades,notu);
fn=lvector(clades);
for (c=0; c<clades; ++c)	fn[c]=F[c+notu];
free_lvector(F);

for (nd=0; nd<clades; ++nd)	{
	/* only examine nodes with sister nodes */
	while (dc[nd]<2 && nd<clades)	++nd;
	if (nd>=clades)	break;
	/* go through node and find daughter clades (converting to clade # as you go) */
	dcl=-1;
	for (d=1; d<=tree[nd][0]; ++d)	if (tree[nd][d]>=notu)	f1[++dcl]=tree[nd][d]-notu;

	/* Now, go through nodes derived from sister node and see how many appear later	*/
	for (d=0; d<dc[nd]; ++d)	{
		cl=f1[d];
		for (c=0; c<dc[nd]; ++c)	{
			if (c!=d)	{
				scl=f1[c];						/* sister clade */
				denom+=treecl[scl][0];
				if (fn[cl]<=fn[scl])	SCIm+=treecl[scl][0];
				else	{
					/* tally all of sister clades descendant nodes that are as old as or younger */
					for (b=1; b<=treecl[scl][0]; ++b)	{
						dcl=treecl[scl][b];		/* daughter clade of sister clade */
						if (fn[dcl]>=fn[cl])	++SCIm;
						}
					}
				/* tally all of sister clades descendant nodes */
				}	/* finish test on condition that c and d are not the same */
			}	/* finish search of sister clades */
		}	/* finish comparisons of clades within node nd */
	}	/* finish search of tree */

SCIm/=denom;

free_ivector(f1);
free_ivector(dc);
free_lvector(fn);
free_lmatrix(treecl,clades,clades);
return SCIm;
}

/* GER - Calculates Gap Excess Ratio after Wills 1998, 1999.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/*		NOTE: this calculates for a CLADOGRAM, not a phylogeny;
/*			Use GERcor to calculate GER for a phylogeny (with branch lengths);
/* Returns:
/*	GER: 1 - (implied - minimum gaps)/(max gaps - min gaps).
****************************************************************************/
double calcGER(long **tree, long *fa, int notu)
{
int		d;
int		clades;
double	mngap=0.0f, mxgap=0.0f, gaps=0.0f;
double	GER=0.0f;

clades=cladecountbytaxa(tree,notu);

mngap=maxlarray(fa,notu)-minlarray(fa,notu);
mxgap=larraytotal(fa,notu)-(minlarray(fa,notu)*notu);

d=naivestratdebt(tree,fa,notu);
gaps=d;

GER=(mxgap-gaps)/(mxgap-mngap);
return GER;
}

//* RCI - Calculates Relative Completeness Index after Benton 1993.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/*		NOTE: this calculates for a CLADOGRAM, not a phylogeny;
/*			Use RCIcor to calculate GER for a phylogeny (with branch lengths);
/* Returns:
/*	RCI - difference bn. observed and minimum gaps over difference bn. max and min.
****************************************************************************/
double calcRCI(long **tree, long**ranges, int notu)
{
int		d, sp;
int		clades;
long	*fa;
double	gaps=0.0f, SRL=0.0f;
double	RCI=0.0f;

fa=lvector(notu);
for (d=0; d<notu; ++d)	fa[d]=ranges[d][0];
clades=cladecountbytaxa(tree,notu);
for (sp=0; sp<notu; ++sp)	SRL=SRL+(1+(ranges[sp][1]-ranges[sp][0]));

d=naivestratdebt(tree,fa,notu);
gaps=d;

RCI=1-(gaps/SRL);
free_lvector(fa);
return RCI;
}

//* MSM - Manhattan Stratigraphic Metric after Siddall ????.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/*		NOTE: this calculates for a CLADOGRAM, not a phylogeny;
/*			Use RCIcor to calculate GER for a phylogeny (with branch lengths);
/* Returns:
/*	MSM - difference bn. observed and minimum gaps over difference bn. max and min.
****************************************************************************/
double calcMSM(long **tree, long *fa, int notu)
{
int		d;
int		clades;
double	gaps=0.0f, mngap=0.0f;
double	MSM=0.0f;

clades=cladecountbytaxa(tree,notu);

mngap=maxlarray(fa,notu)-minlarray(fa,notu);
	
d=naivestratdebt(tree,fa,notu);
gaps=d;
MSM=(mngap/gaps);

return MSM;
}

/* GER - Calculates Gap Excess Ratio after Wills 1998, 1999.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/*		NOTE: this calculates for a PHYLOGENY;
/* Returns:
/*	GER: 1 - (implied - minimum gaps)/(max gaps - min gaps).
****************************************************************************/
double calcGERcor(long **tree, long **ranges, int *bl, int notu)
{
int		d;
int		clades;
int		onset, end;
long	*richness;
double	mngap=0.0f, mxgap=0.0f, gaps=0.0f;
double	GER=0.0f;

clades=cladecountbytaxa(tree,notu);

mxgap=sumlmatrixcol(ranges,notu,0)-(minlmatrixcol(ranges,notu,0)*notu);

onset=minlmatrixcol(ranges,notu,0);
end=maxlmatrixcol(ranges,notu,1);

/* necessary gap exists only if no taxa are sampled over an interval */
richness=richnesstally(ranges, notu);
for (d=onset; d<=end; ++d)
	if (richness[d]==0)	++mngap;

d=stratdebt(tree,ranges,bl,notu);
gaps=d;

free_lvector(richness);

GER=(mxgap-gaps)/(mxgap-mngap);
return GER;
}

//* RCI - Calculates Relative Completeness Index after Benton 1993.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/*		NOTE: this calculates for a CLADOGRAM, not a phylogeny;
/*			Use RCIcor to calculate GER for a phylogeny (with branch lengths);
/* Returns:
/*	RCI - difference bn. observed and minimum gaps over difference bn. max and min.
****************************************************************************/
double calcRCIcor(long **tree, long**ranges, int *bl, int notu)
{
int		d, sp;
int		clades;
double	gaps=0.0f, SRL=0.0f;
double	RCI=0.0f;

clades=cladecountbytaxa(tree,notu);
for (sp=0; sp<notu; ++sp)	SRL=SRL+(1+(ranges[sp][1]-ranges[sp][0]));

d=stratdebt(tree,ranges,bl,notu);
gaps=d;

RCI=1-(gaps/SRL);
return RCI;
}

//* MSM - Manhattan Stratigrahic Metric after Siddall 1997.  
/* NOTE: this metric is completely nonsensical and should not be used.
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/*		NOTE: this calculates for a CLADOGRAM, not a phylogeny;
/*			Use RCIcor to calculate GER for a phylogeny (with branch lengths);
/* Returns:
/*	MSM - difference bn. observed and minimum gaps over difference bn. max and min.
****************************************************************************/
double calcMSMcor(long **tree, long**ranges, int *bl, int notu)
{
int		d;
int		clades;
int		onset, end;
long	*richness;
double	gaps=0.0f, mngap=0.0f;
double	MSM=0.0f;

clades=cladecountbytaxa(tree,notu);

onset=minlmatrixcol(ranges,notu,0);
end=maxlmatrixcol(ranges,notu,1);

richness=richnesstally(ranges, notu);

for (d=onset; d<=end; ++d)
	if (richness[d]==0)	++mngap;
	
d=stratdebt(tree,ranges,bl,notu);
gaps=d;
if (d==0)	MSM=0.0;
else		MSM=(mngap/gaps);

free_lvector(richness);
return MSM;
}


//* stratcompatibilty - determines whether each character in a matrix is compatible with stratigraphy or not.  
/* Requires:
/*	ranges - a matrix of first and last appearances where:
/*		ranges[x][0]=first appearance;
/*		ranges[x][1]=last appearanc;
/*	matrix - cladistic character matrix;
/*	types - type for each character (0: ordered; 1: unordered);
/*	nchars - number of characters;
/*	notu - number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	stratcomp - an array telling you whether the character shows gaps in its  stratigraphic range or not.
****************************************************************************/
unsigned long *stratcompatibilty (long **ranges, long **matrix, int *types, int nchars, int notu, int UNKNOWN, int INAP)
{
int				a, c, q, s, t;
int				fs, ls;
int				mnst, mxst, sts;
int				onset, end;
int				*diverse, **stord;
long			*rich;
unsigned long	*stratcomp;

onset=minlmatrixcol(ranges,notu,0);
end=maxlmatrixcol(ranges,notu,1);

rich=richnesstally(ranges,notu);

diverse=ivector(end+1);
stratcomp=ulvector(nchars);

sts=maxinclmatrix(matrix,notu,nchars,UNKNOWN,INAP);
stord=imatrix(sts+1,3);

for (c=0; c<nchars; ++c)	{
	mnst=colminclmatrix(matrix,notu,c,UNKNOWN,INAP);
	mxst=colmaxclmatrix(matrix,notu,c,UNKNOWN,INAP);
	stratcomp[c]=1;
	for (s=mnst; s<=mxst && stratcomp[c]==1; ++s)	{
		clearivector(diverse,end+1,0);
		for (t=0; t<notu; ++t)	{
			if (matrix[t][c]==s)	for (a=ranges[t][0]; a<=ranges[t][1]; ++a)	++diverse[a];
			}
		q=0;
		for (a=onset; a<=end; ++a)	{
			if (diverse[a]>0 && q==0)						q=1;		/* first find of character 			*/
			else if ((diverse[a]==0 && rich[a]>0) && q==1)	q=2;		/* first absence after being found	*/
			else if (diverse[a]>0 && q==2)	{							/* a gap is demonstrated			*/
				a=end;
				stratcomp[c]=0;
				}
			}	/* search diverse array for gaps between first and last sightings of state s */
		}	/* search each state s in character c for gaps	*/ 
	
	/* one more test for ordered multistates	*/
	if (stratcomp[c]==1 && (types[c]==0 && (mxst-mnst)>1))	{
		/* find the oldest state(s)	*/
		clearimatrix(stord,sts,3,RAND_MAX);
		for (s=0; s<=sts; ++s)	{
			for (t=0; t<notu; ++t)	{
				if (matrix[t][c]==s)	{
					if (stord[s][0]==RAND_MAX)		{
						stord[s][0]=s;
						stord[s][1]=ranges[t][0];
						stord[s][2]=ranges[t][1];
						}
					else {
						if (ranges[t][0]<stord[s][1])	stord[s][1]=ranges[t][0];
						if (ranges[t][1]>stord[s][2])	stord[s][2]=ranges[t][1];
						}
					}	/* end case where taxon has the right state */
				}	/* end search of taxa	*/
			
			if (fs==RAND_MAX && stord[s][1]==onset)	fs=s;
			}	/* end search through states	*/
		
		/* now make sure that everything is adjacent	*/
		/* fs gives the first appearing state; all states above it and below it must appear after that BUT without gap between it and fs	*/
		/*    Also, each state moving away from the starting point has to appear later BUT without gap between it and “prior” state			*/
		ls=fs;	/* because we might skip states, look at the last coded state; this ignores “gaps” in continuous character sequence	*/
		for (s=fs+1; s<=mxst && stratcomp[c]==1; ++s)	{
			while (stord[s][0]==RAND_MAX && s<=mxst)	++s;
			if (s>mxst)	break;
			if (stord[s][1]<stord[ls][1])				stratcomp[c]=0;
			else if (stord[s][1]>(stord[ls][2]+1))		stratcomp[c]=0;
			ls=s;
			}

		ls=fs;	/* because we might skip states, look at the last coded state; this ignores “gaps” in continuous character sequence	*/
		for (s=fs-1; s>=mnst && stratcomp[c]==1; --s)	{
			while (stord[s][0]==RAND_MAX && s>=mnst)	--s;
			if (s<mnst)	break;
			if (stord[s][1]<stord[ls][1])				stratcomp[c]=0;
			else if (stord[s][1]>(stord[ls][2]+1))		stratcomp[c]=0;
			ls=s;
			}
		}	/* end case of ordered multistate	*/
	}	/* test each character for stratigraphic compatibility	*/

free_imatrix(stord,sts,3);
free_ivector(diverse);
free_lvector(rich);
return stratcomp;
}

//* stratcompatfull - determines whether a compatible character pair are compatible with stratigraphy or not.  
/* Requires:
/*	ranges - a matrix of first and last appearances where:
/*		ranges[x][0]=first appearance;
/*		ranges[x][1]=last appearanc;
/*	matrix - cladistic character matrix;
/*	types - type for each character (0: ordered; 1: unordered);
/*	states - number of states for ach character;
/*	ch1 - 1st character being compared;
/*	ch2 - 2nd character being compared;
/*	notu - number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	strcmp - an array telling you whether the character shows gaps in its  stratigraphic range or not.
		strcmp[0]: completely compatible with stratigraphy
		strcmp[1]: first appearances compatible with stratigraphy
		strcmp[2]: all character pairs consistent with stratigraphy
		strcmp[3]: all characters consistent with stratigraphy
****************************************************************************/
unsigned long *stratcompatfull(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int	a, b, c, d, p;
int	cc, sp, max, st1, st2, sc1, sc2;
int	keyst1, keyst2;
unsigned long *strcmp;
unsigned long **combos;
long **pairfnd, **pairrng, **ch1fnd, **ch2fnd, **staterng;
long **pairs;
long **prfl;

max=maxlmatrixcol(ranges, notu, 1);

pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);
pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);
ch1fnd=statefinds(ranges, matrix, states[ch1], ch1, notu, UNKNOWN, INAP);
ch2fnd=statefinds(ranges, matrix, states[ch2], ch2, notu, UNKNOWN, INAP);

clearulvector(strcmp,4,1);

/* determine whether states have consistent ranges	*/
staterng=stateranges(ranges,matrix,states[ch1],ch1,notu,UNKNOWN,INAP);
strcmp[3]=stratcompatchar(staterng,ch1fnd,states[ch1]);
free_lmatrix(staterng,states[ch1],2);
free_lmatrix(ch1fnd,states[ch1],2);

if (strcmp[3]==1)	{
	staterng=stateranges(ranges,matrix,states[ch2],ch2,notu,UNKNOWN,INAP);
	strcmp[3]=stratcompatchar(staterng,ch2fnd,states[ch2]);
	free_lmatrix(staterng,states[ch2],2);
	free_lmatrix(ch2fnd,states[ch2],2);
	}
	
/* determine whether state pairs have consistent ranges	*/
strcmp[2]=stratconsistchpair(pairrng, pairfnd,states[ch1],states[ch2],max);

combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));

for (st1=0; st1<states[ch1]; ++st1)	{
	for (st2=0; st2<states[ch2]; ++st2)	{
		a=st1*states[ch2]+st2;
		if (pairrng[a][0]<max && pairrng[a][0]>=0)	{
			++combos[0][st1];
			++combos[1][st2];
			}
		}
	}

sc1=sc2=0;
for (st1=0; st1<states[ch1]; ++st1)	if (combos[0][st1]>1)	++sc1;
for (st2=0; st2<states[ch2]; ++st2)	if (combos[1][st2]>1)	++sc2;

pairs=lmatrix(p=(states[ch1]+states[ch2]),2);
prfl=lmatrix(p,2);
/*	Find each state paired with multiple states: these are the crux	*/
keyst1=-1;

/* find the state that is paired with the state with multiple	*/
st1=0;
for (c=0; c<sc1; ++c)	{
	clearlmatrix(pairs,states[ch1]+states[ch2],2,-1);
	while (combos[0][st1]<2)	++st1;
	if (st1>=states[ch1])	break;
	keyst1=st1;
	/* find the other matches for keyst1	*/
	p=0;
	for (sp=0; sp<notu; ++sp)	{
		if (matrix[sp][ch1]==keyst1 && (matrix[sp][ch2]!=UNKNOWN && matrix[sp][ch2]!=INAP))	{
			b=0;
			/* make sure that this pair has not been found yet	*/
			for (d=0; d<=p && b==0; ++d)	if (matrix[sp][ch2]==pairs[d][1])	b=1;	
			if (b==0)	{
				pairs[p][0]=keyst1;
				pairs[p][1]=matrix[sp][ch2];
				prfl[p][0]=pairrng[keyst1*states[ch2]+matrix[sp][ch2]][0];
				prfl[p][1]=pairrng[keyst1*states[ch2]+matrix[sp][ch2]][1];
				++p;
				}
			}	/* end search for additional 2nd char state pairs	*/
		}

	keyst2=-1;
	for (st2=0; st2<states[ch2] && keyst2==-1; ++st2)	{
		if (combos[1][st2]>1)	{
			keyst2=st2;
			for (sp=0; sp<notu; ++sp)	{
				if (matrix[sp][ch2]==keyst2 && (matrix[sp][ch1]!=UNKNOWN && matrix[sp][ch1]!=INAP))	{
					b=0;
					/* make sure that we don't add pairs twice	*/ 
					for (d=0; d<p && b==0; ++d)
						if (matrix[sp][ch1]==pairs[d][0])	b=1;

					if (b==0)	{
						pairs[p][0]=matrix[sp][ch1];
						pairs[p][1]=keyst2;
						prfl[p][0]=pairrng[matrix[sp][ch1]*states[ch2]+keyst2][0];
						prfl[p][1]=pairrng[matrix[sp][ch1]*states[ch2]+keyst2][1];
						++p;
						}	/* end pair	*/
					}	/* end case where we've found part of the pair	*/
				}	/*	end search for pairs	*/
			}	/* end search for key state	*/
		}
	/* find the "swing" pair: for 00, 10, 11, this would be 10: it must be either intermediate (00->10->11) or ancestral (00<-10->11)	*/
	/*		thus, only one state pair can be older																						*/
	for (cc=0; pairs[cc][0]!=keyst1 || pairs[cc][1]!=keyst2; cc=cc)	++cc;
	
	/* pair cc can be younger than only one other pair; if younger than two more, then it is incompatible	*/
	b=0;
	for (a=0; a<p && b<2; ++a)
		/* tally state pairs with first appearances older than our "swing" pair	*/
		if (a!=cc && prfl[a][0]<prfl[cc][0])		++b;
	if (b>1)	strcmp[0]=strcmp[1]=0;	/* gap in character state tree	*/
	else	{
		/* tally state pairs with first appearances after the last appearance of our "swing" pair	*/
		b=0;
		for (a=0; a<p && b<1; ++a)	{
			if (a!=cc && prfl[a][0]>(1+prfl[cc][1]))	{
				/* make sure that there is not another state that does not bridge the gap	*/
				for (d=0; d<p; ++d)	{
					b=1;
					/* the next condition will be met only if there are three states paired with keyst1	*/
					if (d!=a && pairs[d][0]==keyst1 && pairs[d][1]!=keyst2)	{
						/* if an intermediate state bridges the gap, then make it compatible again	*/
						if (prfl[d][0]<=(prfl[cc][1]+1) && (prfl[d][1]+1)>=prfl[a][0])	b=0;
						}	/* end look for filling of gap between first and last appearances	*/
					}	/* end look for filling of gap between first and last appearances	*/
				}	/* end case where there is a gap in ranges	*/
			}
		if (b>0)	strcmp[0]=0;		/* gap in character state tree	*/
		}
	++st1;
	}

if (strcmp[0]==1 && (states[ch1]>2 && states[ch2]>2))	{
	st2=0;
	for (c=0; c<sc2; ++c)	{
		clearlmatrix(pairs,states[ch1]+states[ch2],2,-1);
		while (combos[0][st2]<2)	++st2;
		if (st2>=states[ch2])		break;
		keyst2=st2;
		keyst1=-1;
		/* now, find the other matches for character 1	*/
		p=0;
		for (sp=0; sp<notu; ++sp)	{
			if (matrix[sp][ch2]==keyst2 && (matrix[sp][ch1]!=keyst1 && (matrix[sp][ch1]!=UNKNOWN && matrix[sp][ch1]!=INAP)))	{
				b=0;
				/* make sure that this pair has not been found yet	*/
				for (d=0; d<=p && b==0; ++d)	if (matrix[sp][ch1]==pairs[d][1])	b=1;	
				if (b==0)	{
					pairs[p][0]=matrix[sp][ch1];
					pairs[p][1]=keyst2;
					prfl[p][0]=pairrng[matrix[sp][ch1]*states[ch1]+keyst2][0];
					prfl[p][1]=pairrng[matrix[sp][ch1]*states[ch1]+keyst2][0];
					++p;
					}
				}	/* end search for additional 1st char state pairs	*/
			}

		for (st1=0; st1<states[ch1] && keyst1==-1; ++st1)	{
			if (combos[1][st1]>1)	{
				keyst1=st1;
				for (sp=0; sp<notu && p<combos[1][st1]; ++sp)	{
					if (matrix[sp][ch1]==keyst1 && (matrix[sp][ch2]!=UNKNOWN && matrix[sp][ch2]!=INAP))	{
						b=0;
						/* make sure that we don't add pairs twice	*/ 
						for (d=0; d<p && b==0; ++d)
							if (matrix[sp][ch2]==pairs[d][1])	b=1;

						if (b==0)	{
							pairs[p][0]=keyst1;
							pairs[p][1]=matrix[sp][ch2];
							prfl[p][0]=pairrng[(keyst1*states[ch2])+matrix[sp][ch2]][0];
							prfl[p][1]=pairrng[(keyst1*states[ch2])+matrix[sp][ch2]][1];
							++p;
							}	/* end pair	*/
						}	/* end case where we've found part of the pair	*/
					}	/*	end search for pairs	*/
				}	/* end search for key state	*/
			}
		/* find the "swing" pair: for 00, 10, 11, this would be 10: it must be either intermediate (00->10->11) or ancestral (00<-10->11)	*/
		/*		thus, only one state pair can be older																						*/
		for (cc=0; pairs[cc][0]!=keyst1 || pairs[cc][1]!=keyst2; cc=cc)	++cc;
		
		/* pair cc can be younger than only one other pair; if younger than two more, then it is incompatible	*/
		b=0;
		for (a=0; a<p && b<2; ++a)	{
			/* tally state pairs with first appearances older than our "swing" pair	*/
		/* tally state pairs with first appearances older than our "swing" pair	*/
			if (a!=cc && prfl[a][0]<prfl[cc][0])		++b;
			}
		if (b>1)	strcmp[0]=strcmp[1]=0;	/* gap in character state tree	*/
		else	{
			/* tally state pairs with first appearances after the last appearance of our "swing" pair	*/
			b=0;
			for (a=0; a<p && b<1; ++a)	{
				if (a!=cc && prfl[a][0]>(1+prfl[cc][1]))	{
					/* make sure that there is not another state that does not bridge the gap	*/
					for (d=0; d<p; ++d)	{
						b=1;
						/* the next condition will be met only if there are three states paired with keyst2	*/
						if (d!=a && pairs[d][0]==keyst2 && pairs[d][1]!=keyst1)	{
							/* if an intermediate state bridges the gap, then make it compatible again	*/
							if (prfl[d][0]<=(prfl[cc][1]+1) && (prfl[d][1]+1)>=prfl[a][0])	b=0;
							}
						}
					}
				}
			if (b>0)	strcmp[0]=0;		/* gap in character state tree	*/
			}
		}
	++st2;
	}

/* if there is only one state per character that is paired with multiple states, then it is easy	*/
free_lmatrix(pairrng,states[ch1]*states[ch2],2);
free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
free_lmatrix(prfl,states[ch1]+states[ch2],2);
free_lmatrix(pairs,(states[ch1]+states[ch2]),2);
free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));
return strcmp;
}

//* statepairranges - determines the first and last appearance of each pair.  
/* Requires:
/*	ranges - ranges of each taxon:
/*		range[x][y] = 1 if taxon x is found in bin y, 0 if absent;
/*	matrix - matrix of characters for notu taxa
/*	states - the number of states for character ch;
/*	ch - the character examined;
/*	notu - the number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	pairrng - first and last appearance of each pair
****************************************************************************/
long **statepairranges(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int  a, st1, st2, sp, pairs;
long **pairrng;

pairs=1+(states[ch1]*states[ch2])+states[ch2];
pairrng=lmatrix(pairs,2);

for (a=0; a<pairs; ++a)	{
	pairrng[a][0]=RAND_MAX;
	pairrng[a][1]=-1*RAND_MAX;
	}

for (sp=0; sp<notu; ++sp)	{
	if (((st1=matrix[sp][ch1])!=UNKNOWN && (st2=matrix[sp][ch2])!=UNKNOWN) && (matrix[sp][ch1]!=INAP && matrix[sp][ch2]!=INAP))	{
		a=st1*states[ch2]+st2;
		if (ranges[sp][0]<pairrng[a][0])	pairrng[a][0]=ranges[sp][0];
		if (ranges[sp][1]>pairrng[a][1])	pairrng[a][1]=ranges[sp][1];
		}
	}
return	pairrng;
}

//* statepairfinds - determines whether each pair is found in each bin.  
/* Requires:
/*	ranges - ranges of each taxon:
/*		range[x][y] = 1 if taxon x is found in bin y, 0 if absent;
/*	matrix - matrix of characters for notu taxa
/*	states - the number of states for character ch;
/*	ch - the character examined;
/*	notu - the number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	pairfnd - whether each pair is found in each bin
****************************************************************************/
long **statepairfinds (long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int  a, st1, st2, sp, st, pairs;
long **pairfnd;

pairs=1+(states[ch1]*states[ch2])+states[ch2];
a=maxlmatrixcol(ranges,notu,1);
pairfnd=lmatrix(pairs,1+a);

for (sp=0; sp<notu; ++sp)	{
	if (((st1=matrix[sp][ch1])!=UNKNOWN && (st2=matrix[sp][ch2])!=UNKNOWN) && (matrix[sp][ch1]!=INAP && matrix[sp][ch2]!=INAP))	{
		a=st1*states[ch2]+st2;
		for (st=ranges[sp][0]; st<=ranges[sp][1]; ++st)
			if (pairfnd[a][st]==0)			pairfnd[a][st]=1;
		}
	}

return pairfnd;
}


//* stateranges - determines the first and last bin for each state.  
/* Requires:
/*	ranges - ranges of each taxon:
/*		range[x][y] = 1 if taxon x is found in bin y, 0 if absent;
/*	matrix - matrix of characters for notu taxa
/*	states - the number of states for character ch;
/*	ch - the character examined;
/*	notu - the number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	staterng - range of each state
****************************************************************************/
long **stateranges (long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP)
{
int  s, sp;
long **staterng;

staterng=lmatrix(states,2);

for (s=0; s<states; ++s)	{
	staterng[s][0]=RAND_MAX;
	staterng[s][1]=-1*RAND_MAX;
	}

for (sp=0; sp<notu; ++sp)	{
	if ((s=matrix[sp][ch])!=UNKNOWN && matrix[sp][ch]!=INAP)	{
		if (ranges[sp][0]<staterng[s][0])	staterng[s][0]=ranges[sp][0];
		if (ranges[sp][1]>staterng[s][1])	staterng[s][1]=ranges[sp][1];
		}
	}
return	staterng;
}

//* statefinds - determines whether each state occurs in each bin.  
/* Requires:
/*	ranges - ranges of each taxon:
/*		range[x][y] = 1 if taxon x is found in bin y, 0 if absent;
/*	matrix - matrix of characters for notu taxa
/*	states - the number of states for character ch;
/*	ch - the character examined;
/*	notu - the number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	stfnd - whether each state is found in each bin
****************************************************************************/
long **statefinds (long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP)
{
int a, s, sp, st;
long **stfnd;

stfnd=lmatrix(states,1+(a=maxlmatrixcol(ranges,notu,1)));

for (sp=0; sp<notu; ++sp)	{
	if ((s=matrix[sp][ch])!=UNKNOWN && matrix[sp][ch]!=INAP)	{
		for (st=ranges[sp][0]; st<=ranges[sp][1]; ++st)
			if (stfnd[s][st]==0)			stfnd[s][st]=1;
		}
	}

return stfnd;
}


//* stratcompatstate - determines whether all state in a character have continuous ranges.  
/* Requires:
/*	staterng - ranges of each state:
/*		staterng[x][0]=first appearance of state x;
/*		staterng[x][1]=last appearance of state x;
/*	chfnd - whether a state is found in each bin
/*		chfnd[x][y] = 1 if state x is found in bin y, 0 if absent;
/*	s - the state being examined;
/*
/* Returns:
/*	scs - 0 if not continuous, 1 if continuous
****************************************************************************/
unsigned long stratcompatstate(long **staterng, long **chfnd, int s)
{
int st;
unsigned long scs=1;

for (st=staterng[s][0]+1; st<staterng[s][1] && scs==1; ++st)	{
	if (chfnd[s][st]==0)	scs=0;
	}
return scs;
}


//* stratcompatchar - determines whether all state in a character have continuous ranges.  
/* Requires:
/*	staterng - ranges of each state:
/*		staterng[x][0]=first appearance of state x;
/*		staterng[x][1]=last appearance of state x;
/*	states - the number fo states for the character;
/*
/* Returns:
/*	scc - 0 if not all are continuous, 1 if all are continuous
****************************************************************************/
unsigned long stratcompatchar(long **staterng, long **chfnd, int states)
{
int	a, s;
unsigned long scc=1;

for (s=0; s<states && scc==1; ++s)	{
	a=stratcompatstate(staterng, chfnd, s);
	if (a==0)	scc=0;
	}
return scc;
}


//* stratconsistchpair - determines whether all state pairs for two characters have continuous ranges.  
/* Requires:
/*	pairrng - ranges of each pair:
/*		pairrng[x][0]=first appearance of pair x;
/*		pairrng[x][1]=last appearance of pair x;
/*	pr - the pair number, given by (states_char_1 x states_char_2) + state_of_char_2;
/*
/* Returns:
/*	scsp - 0 if incompatible, 1 if compatible
****************************************************************************/
unsigned long stratconsiststpair(long **pairrng, long **pairfnd, int pr)
{
int st;
unsigned long scsp=1;

/* determine whether state pairs have consistent ranges	*/
for (st=pairrng[pr][0]+1; st<pairrng[pr][1]; ++st)	{
	if (pairfnd[pr][st]==0)	scsp=0;
	} 
return scsp;
}

//* stratconsistchpair - determines whether all state pairs for two characters have continuous ranges.  
/* Requires:
/*	pairrng - ranges of each pair:
/*		pairrng[x][0]=first appearance of pair x;
/*		pairrng[x][1]=last appearance of pair x;
/*	st1 - states for 1st character;
/*	st2 - states for 2nd character;
/*	max - oldest stratigraphic bin;
/*
/* Returns:
/*	scp - 0 if incompatible, 1 if compatible
****************************************************************************/
unsigned long stratconsistchpair(long **pairrng, long **pairfnd, int st1, int st2, int max)
{
int pr, cs1, cs2;
unsigned long scp=1;

for (cs1=0; cs1<st1 && scp==1; ++cs1)	{
	for (cs2=0; cs2<st2 && scp==1; ++cs2)	{
		pr=(cs1*st2)+cs2;
		if (pairrng[pr][0]>=0 && pairrng[pr][0]<=max)	
			scp=stratconsiststpair(pairrng, pairfnd,pr);
		}
	}
return scp;
}
