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
/*					02/09 - final assault on stratocompatibility launched!
/**********************************************************************/

#define compatibility_functions
#include "compatibility_functions.h"
#define disparityanalyses
#include "disparityanalyses.h"
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
#define MonteCarloPhylogenyFunctions
#include "MonteCarloPhylogenyFunctions.h"
#define Optimization
#include "optimization.h"
#define sort
#include "sort.h"
#define stratocladistics
#include "stratocladistics.h"
#define TreeRead
#include "tree_read.h"

/** dateclade - Finds age of clades based on the oldest member in the clade
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

/** dateclade - Finds age of clades based on the oldest member in the clade
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	range - First Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	F - modified array fa, with elements notu…notu+clades now filled.
****************************************************************************/
long **datecladefull(long **range, long **tree, int clades, int notu)
{
int	b;
long *fa, *dates, **clrange;

fa=lvector(notu);
for (b=0; b<notu; ++b)	fa[b]=range[b][0];
dates=dateclade(fa,tree,clades,notu);
free_lvector(fa);

clrange=lmatrix(notu+clades,2);
for (b=0; b<notu+clades; ++b)	{
	if (b<notu)	{
		clrange[b][0]=range[b][0];
		clrange[b][1]=range[b][1];
		}
	else
		clrange[b][1]=dates[b];
	}

/* next step: drage clrange[b][0] to clrange[anc][1] for clades, where anc is the ancestral node of b	*/
return clrange;
}

/** datecladesim - Finds age of clades based on the oldest member in the clade
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	tree - a matrix containing the tree structure plus first-last appearance data;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	F - modified array fa, with elements notu…notu+clades now filled.
****************************************************************************/

/** dateclade - Finds age of clades based on the oldest member in the clade
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
double *datecladereal(double *fa, long **tree, int clades, int notu, int sign)
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
		if (F[htu]>F[sp] && sign==-1)		F[htu]=F[sp];
		else if (F[htu]<F[sp] && sign==1)	F[htu]=F[sp];
		}
	}
	
return F;
}


/** datecladeFile - Finds age of clades based on the oldest member in the clade and an fa file
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


/** dateclade - Finds age of clades based on the oldest member in the clade
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

/** datecladeAdd - Finds minimum divergence dates assuming monophyletic groups
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

/** datecladeExtra - Finds minimum divergence dates allowing for ancestors
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


/** BranchAgeBin - Finds age of branches using a discrete scale assuming no ancestors
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

/** BranchAgeBin - Finds age of branches using a continuous scale assuming no ancestors
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

/** blma - Tallies branch lengths in real time allowing for ancestors.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	la - an array of last appearance dates;
/*	bl - an array of branch lengths;
/*	notu - number of taxa;
/*	sign: positive means that 394 is older than 392 (e.g., 394 Ma is befor 392 Ma); negative means -10 is before -9;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	bt - range extension on each branch.
2012-12-29: CHECK THIS!!!!
****************************************************************************/
double *branchlngma(long **tree, double **ranges, int *bl, int notu, int sign)
{
int d, nd, sp, clades, anc;
double *bt;
double *F, *fa;

fa=dvector(notu);
for (d=0; d<notu; ++d)	fa[d]=ranges[d][0];

clades=cladecountbytaxa(tree,notu);
bt=dvector(notu+clades);

F=datecladereal(fa,tree,clades,notu,sign);	/* negative 1 means that 1 is older than 2; positive 1 means that 2 is older than 1	*/

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
		if (sign==-1)	{
			if (anc>-1)	{
				/* if descendent precedes ancestor, then there is a gap    */
				if (F[sp]<F[anc])		bt[sp]=F[anc]-F[sp];
				/* if ancestor disappears before desc appears, then a gap  */
				else if (F[sp]>ranges[anc][1])	bt[sp]=(F[sp]-ranges[anc][1])-1;
				}
			else
				bt[sp]=F[sp]-F[notu+nd];
			}
		else if (sign==1)	{
			if (anc>-1)	{
				/* if descendent precedes ancestor, then there is a gap    */
				if (F[sp]>F[anc])		bt[sp]=F[sp]-F[anc];
				/* if ancestor disappears before desc appears, then a gap  */
				else if (F[sp]<ranges[anc][1])	bt[sp]=(ranges[anc][1]-ranges[sp][0])-1;
				}
			else
				bt[sp]=F[notu+nd]-F[sp];
			}
		}
	}
free_dvector(F);
free_dvector(fa);
return bt;
}

/** blma - Tallies branch lengths in real time allowing for ancestors.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	la - an array of last appearance dates;
/*	bl - an array of branch lengths;
/*	notu - number of taxa;
/*	sign: positive means that 394 is older than 392 (e.g., 394 Ma is befor 392 Ma); negative means -10 is before -9;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	bt - range extension on each branch.
2012-12-29: CHECK THIS!!!!
****************************************************************************/
double **branchlength_ma_to_draw(long **tree, double **ranges, int *bl, int notu, int sign, double min)
{
int d, nd, sp, clades, anc, htu;
double **bt;

clades=cladecountbytaxa(tree,notu);
bt=dmatrix(notu+clades,2);			/* bt[x][0] gives onset of ghost lineage/ghost taxon; bt[x][1] gives end of ghost lineage/taxon			*/
									/* if x is a clade, then make sure that bt[x][0] is always min before bt[x][1] for artistic purposes	*/

for (sp=0; sp<notu; ++sp)	bt[sp][0]=bt[sp][1]=ranges[sp][0];	/* extension ends at FA	*/

/* start at top of tree and work down	*/
/* bt[node][1] = oldest taxon in the clade	*/
/* bt[node][0] = bt[node][1]-min at first: we'll lower that later	*/
/* if there is no ancestor,then bt[sp][0]= min(bt[sp][x]) for all taxa in clade;
/* if there is an ancestor,then bt[sp][0] is either FA OR LA of ancestor LA[anc] precedes FA[descendant]	*/
for (nd=clades-1; nd>=0; --nd)	{
	htu=nd+notu;
	bt[htu][1]=bt[tree[nd][1]][0];		/* use the FA of the first taxon in the clade to get started	*/
	anc=-1;
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		if (sp<notu && bl[sp]==0)	anc=sp;
		if (bt[sp][0]<bt[htu][1])	bt[htu][1]=bt[sp][0];
		}
	/* if no ancestor,then draw everything down to oldest species	*/
	if (anc==-1)	{
		bt[htu][0]=bt[htu][1]+(((double) sign)*min);	/* draw node down some minimum length	*/
		for (d=1; d<=tree[nd][0]; ++d)	{
			sp=tree[nd][d];
			bt[sp][0]=bt[htu][1];						/* ghost lineage/taxon set at clade age	*/
			}
		}
	else	{
		bt[htu][0]=bt[htu][1];							/* do not draw node down unless required later	*/
		bt[anc][0]=bt[htu][1];		/* onset of ancestor must equal end point of node	*/
		
		for (d=1; d<=tree[nd][0]; ++d)	{
			sp=tree[nd][d];
			if (sp!=anc)	{
				/* if descendant appears after last appearance of ancestor, then draw it "down" to the ancestor	*/
				/*  otherwise, leave it alone; note taht we use ranges here because ancestor must be an OTU		*/
				if (bt[sp][0]>ranges[anc][1])	bt[sp][0]=ranges[anc][1];
				}	/* end examination of descendant species	*/
			}	/* end search of clade	*/
		}	/* end case of ancestors	*/
	}
return bt;
}


/** CladeSurvive - Determines how long clades last
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
/** CladeSurviveAdd - Finds minimum divergence dates assuming monophyletic groups
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

/** temporalbranchlength - Routine to get basic branch lengths in time units.  
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

/** temporalbranchlength - Routine to get basic branch lengths in time units.  
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
double **temporalbranchlengthequal(long **tree, long *fba, double *byr, int clades, int notu)
{
int		a,b,c,d,g,s;
int		anc, cl, sp, ties;
int		*lumped;
long	*ancestors, *height, *taxa, *cladebins;
long	**vtree;
double	en, f, l;
double	*bt;
double	**phy;

/* get list of ancestors */
ancestors=((long*) listancestor(tree,clades,notu));
height=((long*) patristicheight(tree,clades,notu));
bt=dvector(notu+clades);
vtree=lmatrix(clades,notu);
equallmatrix(vtree,tree,clades,notu);
vtree=clademember(vtree,notu,clades);

/* initially set birth to middle of bin	*/
/*for (sp=0; sp<notu; ++sp)		bt[sp]=(byr[fba[sp]]+byr[1+fba[sp]])/2;	*/

cladebins=dateclade(fba, tree, clades, notu);
/*cladeages=datecladereal(fa,tree,clades,notu);	*/
for (sp=0; sp<clades+notu; ++sp)		bt[sp]=(byr[cladebins[sp]]+byr[1+cladebins[sp]])/2;

lumped=ivector(clades);

taxa=lvector(notu);
for (sp=0; sp<notu; ++sp)	taxa[sp]=sp;
taxa=sort_declongbylong(taxa,height,notu);
taxa=sort_declongbylong(taxa,fba,notu);

for (s=0; s<notu; ++s)	{
	sp=taxa[s];
	cl=(anc=ancestors[sp])-notu;
	while (cl<0 && s<notu)	{
		sp=taxa[++s];
		cl=(anc=ancestors[sp])-notu;
		}
	if (s>=notu)	break;
	f=bt[sp];
	l=bt[anc];
	
	/* If a taxon's clade already has been adjusted, then the ancestor might be younger		*/
	/* make sure that all descendants of that ancestor are at least as old as that ancestor	*/
	if (f<l)	{
		bt[sp]=f=l;
/*		cl=anc-notu;
/*		for (d=1; d<=vtree[cl][0]; ++d)	{
/*			g=vtree[cl][d];
/*			if (bt[g]<bt[anc]
/*			}	*/
		}
	ties=0;
	lumped[0]=sp;
	while (f==l && cl>0)	{
		lumped[++ties]=anc;
		anc=ancestors[cl+notu];
		cl=anc-notu;
		f=l;

		if (f<l)	{
			f=l;
			for (a=0; a<=ties; ++a)	{
				b=lumped[a];
				bt[b]=l;
				}
			}

		if (anc>=0)	l=bt[anc];
		else		l=0;
		}
	lumped[ties+1]=anc;

	if (ties>0)	{
		/* reset birth to the end of the bin	*/
		bt[sp]=byr[1+fba[sp]];
		
		a=fba[sp]+1;
		en=(byr[a]-bt[anc])/((double) (ties+1));
		
		/* if the first clade already has been adjusted, then start from that; otherwise, start from its ancestor	*/
		a=lumped[ties];
		if (bt[a]!=(byr[fba[a]]+byr[fba[a]+1])/2)	--ties;
		
		for (a=ties; a>=0; --a)	{
			b=lumped[a];
			c=lumped[a+1];
			/* if b is a clade,then make sure all of its descendants have their ages modified accordingly	*/
			if (b>notu && a>1)	{
				cl=b-notu;
				for (d=1; d<=vtree[cl][0]; ++d)	{
					g=vtree[cl][d];
					if (bt[g]==bt[b])	bt[g]=bt[c]+en;
					}
				}
			bt[b]=bt[c]+en;
			}
		}
	}
	
/*for (sp=0; sp<notu+clades; ++sp)	{
	if (ancestors[sp]>=notu)	{
		anc=ancestors[sp];
		bt[1][sp]=bt[0][sp]-bt[0][anc];
		}
	}	*/

phy=dmatrix(notu+clades,2);
for (a=0; a<notu+clades; ++a)	{
	phy[a][0]=bt[a];
	anc=ancestors[a];
	phy[a][1]=bt[a]-bt[anc];
	}

free_dvector(bt);
free_lvector(height);
free_lvector(taxa);
free_lmatrix(vtree,clades,notu);
free_ivector(lumped);
free_lvector(ancestors);
free_lvector(cladebins);
return phy;
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

/** naivestratdebt - Tallies stratigraphic gaps assuming no ancestors.  
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

/** stratdebt - Tallies stratigraphic gaps in bins (e.g., stages) allowing for ancestors.  
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

/** stratdebtnodes - Tallies stratigraphic gaps in bins (e.g., stages) allowing for ancestors.  
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
long *stratdebtnodes(long **tree, long**ranges, int *bl, int notu)
{
int d, nd, sp, clades, anc;
long	sd=0;
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
return fa;		/* CHANGE THIS!!!	*/
}

/** stratdebtreal - Tallies stratigraphic gaps in real time allowing for ancestors.  
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
F=datecladereal(fa,tree,clades,notu,-1);	/* negative 1 means that 1 is older than 2; positive 1 means that 2 is older than 1	*/

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

/** SCI - Calculates Stratigraphic Consistency Index after Huelsenbeck 1994.  
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

/** SCIm - Calculates modified Stratigraphic Consistency Index after Habib & Gittelman.  
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

/** RCI - Calculates Relative Completeness Index after Benton 1993.  
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

/** MSM - Manhattan Stratigraphic Metric after Siddall ????.  
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

/** RCI - Calculates Relative Completeness Index after Benton 1993.  
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

/** MSM - Manhattan Stratigrahic Metric after Siddall 1997.  
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


/** stratcompatibilty - determines whether each character in a matrix is compatible with stratigraphy or not.  
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
unsigned long *stratcompatibility (long **ranges, long **matrix, int *types, int nchars, int notu, int UNKNOWN, int INAP)
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
	/* find minimum and maximum states	*/
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

/** stratcompatfull - determines whether a compatible character pair are compatible with stratigraphy or not.  
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
/*	stratcomp - an array telling you whether the character shows gaps in its  stratigraphic range or not.
		stratcomp[0]: character 1
		stratcomp[1]: character 2
		stratcomp[2]: character 1 states
		stratcomp[3]: character 2 states
		stratcomp[4]: state comparisons (2 for binary; possibly more for multistate)
		stratcomp[5]: first appearances compatible with stratigraphy
		stratcomp[6]: pair synoptic ranges completely compatible with stratigraphy
		stratcomp[7]: taxon synoptic ranges completely compatible with stratigraphy
		stratcomp[8]: divergent (0) vs. hierarchical (1) stratigraphic compatibility
****************************************************************************/
unsigned long *stratcompatfull(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int	a, b, c, d, f;
int	cc, sp, max, st1, st2, sc1, sc2, rc, scf1, scf2, scu1, scu2, prs;
int cl[20], cr[20];
int	keyst1, keyst2, sw1, sw2;				/* these are the "swing" states between characters + important info	*/
unsigned long *stratcomp;
unsigned long **combos;
/*long **pairfnd;*/
long **pairs, **pairrng, **prfl, **pairfnd, *gaps;

max=maxlmatrixcol(ranges, notu, 1);

pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/

stratcomp=ulvector(9);
clearulvector(stratcomp,9,0);

stratcomp[0]=ch1;
stratcomp[1]=ch2;
stratcomp[2]=states[ch1];
stratcomp[3]=states[ch2];

combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));
for (st1=0; st1<states[ch1]; ++st1)	{
	for (st2=0; st2<states[ch2]; ++st2)	{
		a=st1*states[ch2]+st2;
		if (pairrng[a][0]<=max && pairrng[a][0]>=0)	{
			++combos[0][st1];
			++combos[1][st2];
			}
		}
	}

sc1=sc2=0;
for (st1=0; st1<states[ch1]; ++st1)	if (combos[0][st1]>1)	++sc1;
for (st2=0; st2<states[ch2]; ++st2)	if (combos[1][st2]>1)	++sc2;

pairs=lmatrix(prs=(states[ch1]+states[ch2]),2);	/* state pairs between two characters		*/
prfl=lmatrix(prs,2);								/* first/last appearance of character pairs	*/
gaps=lvector(prs);

/*	Find each state paired with multiple states: these are the crux	*/
keyst1=-1;

/* play off of the character with the most swing states	*/
if (sc2>sc1)	{
	free_lmatrix(pairrng,states[ch1]*states[ch2],2);
	free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
	free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));
	c=ch1;
	ch1=ch2;
	ch2=c;
	c=sc1;
	sc1=sc2;
	sc2=c;
	pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
	pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/
	combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));
	for (st1=0; st1<states[ch1]; ++st1)	{
		for (st2=0; st2<states[ch2]; ++st2)	{
			a=st1*states[ch2]+st2;
			if (pairrng[a][0]<=max && pairrng[a][0]>=0)	{
				++combos[0][st1];
				++combos[1][st2];
				}
			}
		}
	}	

/*if (mxcm>1)	{	*/
st1=0;
for (c=0; c<sc1; ++c)	{
	clearlmatrix(pairs,states[ch1]+states[ch2],2,-1);
	while (combos[0][st1]<2)	++st1;
	if (st1>=states[ch1])	break;
	keyst1=st1;			/* the swing state is found!	*/
	/* find the other matches for keyst1	*/
	prs=0;		/* p is the number of pairs	*/
	clearlvector(gaps,(states[ch1]+states[ch2]),0);
	for (sp=0; sp<notu; ++sp)	{
		if (matrix[sp][ch1]==keyst1 && (matrix[sp][ch2]!=UNKNOWN && matrix[sp][ch2]!=INAP))	{
			b=0;
			/* make sure that this pair has not been found yet	*/
			for (d=0; d<=prs && b==0; ++d)	if (matrix[sp][ch2]==pairs[d][1])	b=1;	
			if (b==0)	{
				pairs[prs][0]=keyst1;					/* pair #p st. 1 found	*/
				pairs[prs][1]=matrix[sp][ch2];		/* pair #p st. 2 found	*/
				prfl[prs][0]=pairrng[(keyst1*states[ch2])+pairs[prs][1]][0];	/* FA of pair p	*/
				prfl[prs][1]=pairrng[(keyst1*states[ch2])+pairs[prs][1]][1];	/* LA of pair p	*/
				/* look for gaps in the range of pair prs	*/
				for (f=prfl[prs][0]+1; f<prfl[prs][1]; ++f)
					if (pairfnd[(keyst1*states[ch2])+matrix[sp][ch2]][f]==0)	++gaps[prs];
				++prs;	/* increment pairs	*/
				}
			}	/* end search for additional 2nd char state pairs	*/
		}

	keyst2=-1;
	for (st2=0; st2<states[ch2] && keyst2==-1; ++st2)	{
		if (combos[1][st2]>1)	{
			for (sp=0; sp<notu; ++sp)	{
				/* look at this only if this swing is paired with the swing from character 1: 
					this will not always be true for multistates; 2009.02.10						*/
				if (matrix[sp][ch2]==st2 && matrix[sp][ch1]==keyst1)	{
					keyst2=st2;
					sp=notu;		/* keystate 2 found, so kill the loop	*/
					}
				}

			if (keyst2>-1)	{
				rc=1;	/* remaining combos: when this hits combo[1][keyst2], then end	*/
				for (sp=0; sp<notu & rc<combos[1][keyst2]; ++sp)	{
					if (matrix[sp][ch2]==keyst2 && (matrix[sp][ch1]!=UNKNOWN && matrix[sp][ch1]!=INAP))	{
						b=0;
						/* make sure that we don't add pairs twice	*/ 
						for (d=0; d<prs && b==0; ++d)
							if (matrix[sp][ch1]==pairs[d][0])	b=1;

						if (b==0)	{
							pairs[prs][0]=matrix[sp][ch1];	/* pair #p st. 1 found	*/
							pairs[prs][1]=keyst2;				/* pair #p st. 2 found	*/
							prfl[prs][0]=pairrng[(pairs[prs][0]*states[ch2])+keyst2][0];	/* FA of pair p	*/
							prfl[prs][1]=pairrng[(pairs[prs][0]*states[ch2])+keyst2][1];	/* LA of pair p	*/
							/* look for gaps in the range of pair p	*/
							for (f=prfl[prs][0]+1; f<prfl[prs][1]; ++f)
								if (pairfnd[(pairs[prs][0]*states[ch2])+keyst2][f]==0)	++gaps[prs];
							++prs;	/* increment pairs	*/
							++rc;	/* one more of the "linked" pairs found	*/
							}	/* end pair	*/
						}	/* end case where we've found part of the pair	*/
					}	/*	end search for pairs	*/
				/* find cc, the "swing" pair; for 00, 10, 11, this would be 10: 							*/
				/*		it must be either intermediate (00->10->11 or 11->01->00) or ancestral (00<-10->11)	*/
				/*		thus, only one state pair can be older												*/
				/* ADD cl & cr: the "left" and "right" combos, to which this one will be compared; damn, I'm good..... */
				cc=-1;
				scf1=scf2=sw1=sw2=0;
				for (d=0; d<prs && cc==-1; ++d)	{
					if (pairs[d][0]==keyst1 && pairs[d][1]==keyst2)	cc=d;	/* "swing" pair	*/
					}
				/* now examine other state pairs with keystate 1	*/
				scu1=1;		/* this notes a gap between swing pairs and end pairs: it takes only one to set this to 0	*/
				for (d=0; d<prs; ++d)	{
					if (pairs[d][0]==keyst1 && d!=cc)	{
						cl[sw1]=d;											/* other pairs with key state 1	*/
						++sw1;
						/* if "end" pair appears before "swing" pair	*/
						if (prfl[d][0]<prfl[cc][0])	{
							scf1=1;						/* flag that end pair precedes swing pair	*/
							/* if there is no gap between the "end" pair's LA and the "swing" pair's FA	*/
							/* NOTE: in the case of multistates, only ONE of the "end" pairs need overlap with swing pair	*/
							if (prfl[d][1]>=(prfl[cc][0]-1))	scu1=0;
							/* d=p;	*/
							}	/* end case where "end" pair appears before "swing" pair	*/
						/* if there is not gap between the "swing" pair's LA and the "end" pair's FA, then erase gap	*/
						else if (prfl[cc][1]>=(prfl[d][0]-1))	scu1=0;
						}	/* end examination of "swing" pair involving keystate 1	*/
					}	/* end search for "swing" pairs involving keystate 1	*/
				/* now examine other state pairs with keystate 2	*/
				scu2=1;
				for (d=0; d<prs; ++d)	{
					if (pairs[d][1]==keyst2 && d!=cc)	{
						cr[sw2]=d;											/* other pairs with key state 2	*/
						++sw2;
						/* if "end" pair appears before "swing" pair	*/
						if (prfl[d][0]<prfl[cc][0])	{
							scf2=1;						/* flag that end pair precedes swing pair	*/
							/* if there is no gap between the "end" pair's LA and the "swing" pair's FA	*/
							/* NOTE: in the case of multistates, only ONE of the "end" pairs need overlap with swing pair	*/
							if (prfl[d][1]>=(prfl[cc][0]-1))	scu2=0;
							}	/* end case where "end" pair appears before "swing" pair	*/
						/* if there is not gap between the "swing" pair's LA and the "end" pair's FA, then erase gap	*/
						else if (prfl[cc][1]>=(prfl[d][0]-1))	scu2=0;
						}	/* end examination of "swing" pair involving keystate 2	*/
					}	/* end search for "swing" pairs involving key state 2	*/					
				}	/* end case where 2nd keystate is found	*/
				
			/* now, tally up stratigraphic compatibility	*/
			++stratcomp[4];						/* another comparison made				*/
			if ((a=scf1+scf2)<2)	{
				++stratcomp[5];									/* general stratigraphic compatibility			*/
				if (scu1==0 && scu2==0)			++stratcomp[6];	/* range stratigraphic compatibility			*/
				if ((b=maxlarray(gaps,prs))==0)	++stratcomp[7];	/* no gaps in the relevant character pairs		*/
				if (a==1)						++stratcomp[8];	/* 10->00->01 instead of 10<-00->01				*/
				}
			}	/* end test of whether "swing state" is key state 2 	*/
		}	/* end search for 2nd key state	*/
	++st1;	/* if there are multiple swing states, then move on to the next possible one!	*/
	}	/* end survey of "swing" states matched	*/

free_lmatrix(pairrng,states[ch1]*states[ch2],2);
free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
free_lmatrix(prfl,states[ch1]+states[ch2],2);
free_lmatrix(pairs,(states[ch1]+states[ch2]),2);
free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));
free_lvector(gaps);
return stratcomp;
}

/** stratcompatfullplus - determines whether a compatible character pair are compatible with stratigraphy or not.  
/* Requires:
/*	ranges - a matrix of first and last appearances where:
/*		ranges[x][0]=first appearance;
/*		ranges[x][1]=last appearanc;
/*	matrix - cladistic character matrix;
/*	states - number of states for ach character;
/*	charcomps - compatibility of each character;
/*	ch1 - 1st character being compared;
/*	ch2 - 2nd character being compared;
/*	notu - number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	stratcomp - an array telling you whether the character shows gaps in its  stratigraphic range or not.
		stratcomp[0]: character 1
		stratcomp[1]: character 2
		stratcomp[2]: character 1 states
		stratcomp[3]: character 2 states
		stratcomp[4]: state comparisons (2 for binary; possibly more for multistate)
		stratcomp[5]: first appearances compatible with stratigraphy
		stratcomp[6]: pair synoptic ranges completely compatible with stratigraphy
		stratcomp[7]: taxon synoptic ranges completely compatible with stratigraphy
		stratcomp[8]: divergent (0) vs. hierarchical (1) stratigraphic compatibility
		stratcomp[9]: hierarchical stratigraphic compatibility consistent with anagenesis
		stratcomp[10]: hierarchical stratigraphic compatibility consistent with budding
		stratcomp[11]: compatibility & stratigraphy suggests 00->11->10 transition rather than 10<-00->11
		stratcomp[12]: compatibility & stratigraphy indecisive about 00->11->10 vs 10<-00->11
***************************************************************************************************************************/
unsigned long *stratcompatfullplus(long **ranges, long **matrix, int *states, unsigned long *charcomps, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int	a, b, c, d, f, g, h;
int	cc, sp, max, st1, st2, sc1, sc2, rc, scf1, scf2, scu1, scu2, prs, m1, m2;
int cl[20], cr[20];
int	keyst1, keyst2, sw1, sw2, swing;				/* these are the "swing" states between characters + important info, with "swing" giving pair #	*/
unsigned long *stratcomp;
unsigned long **combos;
/*long **pairfnd;*/
long **pairs, **pairrng, **prfl, **pairfnd, *gaps;

max=maxlmatrixcol(ranges, notu, 1);

pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/

stratcomp=ulvector(13);
clearulvector(stratcomp,13,0);

stratcomp[0]=ch1;
stratcomp[1]=ch2;
stratcomp[2]=states[ch1];
stratcomp[3]=states[ch2];

combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));
for (st1=0; st1<states[ch1]; ++st1)	{
	for (st2=0; st2<states[ch2]; ++st2)	{
		a=st1*states[ch2]+st2;
		if (pairrng[a][0]<=max && pairrng[a][0]>=0)	{
			++combos[0][st1];
			++combos[1][st2];
			}
		}
	}

sc1=sc2=0;
for (st1=0; st1<states[ch1]; ++st1)	if (combos[0][st1]>1)	++sc1;
for (st2=0; st2<states[ch2]; ++st2)	if (combos[1][st2]>1)	++sc2;

pairs=lmatrix(prs=(states[ch1]+states[ch2]),2);	/* state pairs between two characters		*/
prfl=lmatrix(prs,2);							/* first/last appearance of character pairs	*/
gaps=lvector(prs);

/*	Find each state paired with multiple states: these are the crux	*/
keyst1=-1;

/* play off of the character with the most swing states	*/
if (sc2>sc1)	{
	free_lmatrix(pairrng,states[ch1]*states[ch2],2);
	free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
	free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));
	c=ch1;
	ch1=ch2;
	ch2=c;
	c=sc1;
	sc1=sc2;
	sc2=c;
	pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
	pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/
	combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));
	for (st1=0; st1<states[ch1]; ++st1)	{
		for (st2=0; st2<states[ch2]; ++st2)	{
			a=st1*states[ch2]+st2;
			if (pairrng[a][0]<=max && pairrng[a][0]>=0)	{
				++combos[0][st1];
				++combos[1][st2];
				}
			}
		}
	}	

/*if (mxcm>1)	{	*/
st1=0;
for (c=0; c<sc1; ++c)	{
	clearlmatrix(pairs,states[ch1]+states[ch2],2,-1);
	while (combos[0][st1]<2)	++st1;
	if (st1>=states[ch1])	break;
	keyst1=st1;			/* the swing state is found!	*/
	/* find the other matches for keyst1	*/
	prs=0;		/* p is the number of pairs	*/
	clearlvector(gaps,(states[ch1]+states[ch2]),0);
	for (sp=0; sp<notu; ++sp)	{
		if (matrix[sp][ch1]==keyst1 && (matrix[sp][ch2]!=UNKNOWN && matrix[sp][ch2]!=INAP))	{
			b=0;
			/* make sure that this pair has not been found yet	*/
			for (d=0; d<=prs && b==0; ++d)	if (matrix[sp][ch2]==pairs[d][1])	b=1;	
			if (b==0)	{
				pairs[prs][0]=keyst1;					/* pair #p st. 1 found	*/
				pairs[prs][1]=matrix[sp][ch2];		/* pair #p st. 2 found	*/
				prfl[prs][0]=pairrng[(keyst1*states[ch2])+pairs[prs][1]][0];	/* FA of pair p	*/
				prfl[prs][1]=pairrng[(keyst1*states[ch2])+pairs[prs][1]][1];	/* LA of pair p	*/
				/* look for gaps in the range of pair prs	*/
				for (f=prfl[prs][0]+1; f<prfl[prs][1]; ++f)
					if (pairfnd[(keyst1*states[ch2])+matrix[sp][ch2]][f]==0)	++gaps[prs];
				++prs;	/* increment pairs	*/
				}
			}	/* end search for additional 2nd char state pairs	*/
		}

	keyst2=-1;
	for (st2=0; st2<states[ch2] && keyst2==-1; ++st2)	{
		if (combos[1][st2]>1)	{
			for (sp=0; sp<notu; ++sp)	{
				/* look at this only if this swing is paired with the swing from character 1: 
					this will not always be true for multistates; 2009.02.10						*/
				if (matrix[sp][ch2]==st2 && matrix[sp][ch1]==keyst1)	{
					keyst2=st2;
					sp=notu;		/* keystate 2 found, so kill the loop	*/
					}
				}

			if (keyst2>-1)	{
				rc=1;	/* remaining combos: when this hits combo[1][keyst2], then end	*/
				for (sp=0; sp<notu & rc<combos[1][keyst2]; ++sp)	{
					if (matrix[sp][ch2]==keyst2 && (matrix[sp][ch1]!=UNKNOWN && matrix[sp][ch1]!=INAP))	{
						b=0;
						/* make sure that we don't add pairs twice	*/ 
						for (d=0; d<prs && b==0; ++d)
							if (matrix[sp][ch1]==pairs[d][0])	b=1;

						if (b==0)	{
							pairs[prs][0]=matrix[sp][ch1];	/* pair #p st. 1 found	*/
							pairs[prs][1]=keyst2;				/* pair #p st. 2 found	*/
							prfl[prs][0]=pairrng[(pairs[prs][0]*states[ch2])+keyst2][0];	/* FA of pair p	*/
							prfl[prs][1]=pairrng[(pairs[prs][0]*states[ch2])+keyst2][1];	/* LA of pair p	*/
							/* look for gaps in the range of pair p	*/
							for (f=prfl[prs][0]+1; f<prfl[prs][1]; ++f)
								if (pairfnd[(pairs[prs][0]*states[ch2])+keyst2][f]==0)	++gaps[prs];
							++prs;	/* increment pairs	*/
							++rc;	/* one more of the "linked" pairs found	*/
							}	/* end pair	*/
						}	/* end case where we've found part of the pair	*/
					}	/*	end search for pairs	*/
				/* find cc, the "swing" pair; for 00, 10, 11, this would be 10: 							*/
				/*		it must be either intermediate (00->10->11 or 11->01->00) or ancestral (00<-10->11)	*/
				/*		thus, only one state pair can be older												*/
				/* ADD cl & cr: the "left" and "right" combos, to which this one will be compared; damn, I'm good..... */
				cc=-1;
				scf1=scf2=sw1=sw2=0;
				for (d=0; d<prs && cc==-1; ++d)	{
					if (pairs[d][0]==keyst1 && pairs[d][1]==keyst2)	cc=d;	/* "swing" pair	*/
					}
				/* now examine other state pairs with keystate 1	*/
				scu1=1;		/* this notes a gap between swing pairs and end pairs: it takes only one to set this to 0	*/
				for (d=0; d<prs; ++d)	{
					if (pairs[d][0]==keyst1 && d!=cc)	{
						cl[sw1]=d;											/* other pairs with key state 1	*/
						++sw1;
						/* if "end" pair appears before "swing" pair	*/
						if (prfl[d][0]<prfl[cc][0])	{
							scf1=1;						/* flag that end pair precedes swing pair	*/
							/* if there is no gap between the "end" pair's LA and the "swing" pair's FA	*/
							/* NOTE: in the case of multistates, only ONE of the "end" pairs need overlap with swing pair	*/
							if (prfl[d][1]>=(prfl[cc][0]-1))	scu1=0;
							/* d=p;	*/
							}	/* end case where "end" pair appears before "swing" pair	*/
						/* if there is not gap between the "swing" pair's LA and the "end" pair's FA, then erase gap	*/
						else if (prfl[cc][1]>=(prfl[d][0]-1))	scu1=0;
						}	/* end examination of "swing" pair involving keystate 1	*/
					}	/* end search for "swing" pairs involving keystate 1	*/
				/* now examine other state pairs with keystate 2	*/
				scu2=1;
				for (d=0; d<prs; ++d)	{
					if (pairs[d][1]==keyst2 && d!=cc)	{
						cr[sw2]=d;											/* other pairs with key state 2	*/
						++sw2;
						/* if "end" pair appears before "swing" pair	*/
						if (prfl[d][0]<prfl[cc][0])	{
							scf2=1;						/* flag that end pair precedes swing pair	*/
							/* if there is no gap between the "end" pair's LA and the "swing" pair's FA	*/
							/* NOTE: in the case of multistates, only ONE of the "end" pairs need overlap with swing pair	*/
							if (prfl[d][1]>=(prfl[cc][0]-1))	scu2=0;
							}	/* end case where "end" pair appears before "swing" pair	*/
						/* if there is not gap between the "swing" pair's LA and the "end" pair's FA, then erase gap	*/
						else if (prfl[cc][1]>=(prfl[d][0]-1))	scu2=0;
						}	/* end examination of "swing" pair involving keystate 2	*/
					}	/* end search for "swing" pairs involving key state 2	*/					
				}	/* end case where 2nd keystate is found	*/
				
			/* now, tally up stratigraphic compatibility	*/
			++stratcomp[4];						/* another comparison made				*/
			if ((a=scf1+scf2)<2)	{
				++stratcomp[5];									/* general stratigraphic compatibility			*/
				if (scu1==0 && scu2==0)			++stratcomp[6];	/* range stratigraphic compatibility			*/
				if ((b=maxlarray(gaps,prs))==0)	++stratcomp[7];	/* no gaps in the relevant character pairs		*/
				if (a==1 && prs>=3)	{
					++stratcomp[8];	/* 10->00->01 instead of 10<-00->01				*/
					/* 2011-05-16 Determine whether ranges overlap	*/
					/* go through relevant state pairs and determine whether LA[older]<=FA[younger]	*/
					/* this uses prfl for the 3 pairs in question: just need to figure out which 3 they are all of these years later....	*/
					/* OK, try this: 
							1) find swing pair based on keyst1 & keyst2.  IF we are here, then it is intermediate in age.
							2) Find the older pair with one of the swing states: did it disappear at same time or earlier?
								a) repeat if there are multiple pairs with the swing state
							3) Find the younger pair based on keyst1 & keyst2: did it appear at same time or later?	*/
					swing=-1;
					for (d=0; d<prs; ++d)	{
						if (pairs[d][0]==keyst1 && pairs[d][1]==keyst2)	{
							swing=d;
							d=prs;
							}
						}
					if (swing!=-1)	{
						g=h=0;
						for (d=0; d<prs; ++d)	{
							if (d!=swing && (pairs[d][0]==keyst1 || pairs[d][1]==keyst2))	{
								/* if state pair appears after or at same time swing pair disappears, then consistent with anagenesis	*/
								if (prfl[d][0]>=prfl[swing][1])	{
									if (g==0)	++stratcomp[9];	/* consistent with anagenesis			*/
									else		++stratcomp[10];	/* only one case of anagenesis possible	*/
									g=1;
									}
								/* if swing pair appears after or at same time state pair disappears, then consistent with anagenesis	*/
								else if (prfl[swing][0]>=prfl[d][1])	{
									if (h==0)	++stratcomp[9];	/* consistent with anagenesis	*/	
									else		++stratcomp[10];	/* only one case of anagenesis possible	*/
									h=1;
									}
								/* otherwise, we have overlapping ranges: consistent with budding	*/
								else			++stratcomp[10];	/* consistent with budding		*/
								}	/* end comparisons with other state pairs	*/
							}	/* end search of other state pairs	*/
						}	/* double check: but there really has to be a swing pair in here!	*/
					}	/* end case of hierarchically compatible pair	*/
				}	/* end case of stratigraphically compatible characters	*/
			
			else	{
				/* 2011-05-16: Find swing pair and determine which character "reverse" to make this a literal interpretation of the fossil record	*/
				/* Try this: 
					1) find swing pair;
					2) find which of the two states is older	*/
				m1=m2=RAND_MAX;
				for (d=0; d<prs; ++d)	{
					if (pairs[d][0]==keyst1 && prfl[d][0]<m1)	m1=prfl[d][0];
					if (pairs[d][1]==keyst2 && prfl[d][0]<m2)	m2=prfl[d][0];
					}
				/* if the older state is from the less compatible character, then increment	*/
				if ((m1<m2 && charcomps[ch2]>charcomps[ch1]) || (m2<m1 && charcomps[ch1]>charcomps[ch2]))	++stratcomp[11];
				else if (m2==m1 || charcomps[ch2]==charcomps[ch1])											++stratcomp[12];
				}
			}	/* end test of whether "swing state" is key state 2 	*/
		}	/* end search for 2nd key state	*/
	++st1;	/* if there are multiple swing states, then move on to the next possible one!	*/
	}	/* end survey of "swing" states matched	*/

free_lmatrix(pairrng,states[ch1]*states[ch2],2);
free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
free_lmatrix(prfl,states[ch1]+states[ch2],2);
free_lmatrix(pairs,(states[ch1]+states[ch2]),2);
free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));
free_lvector(gaps);
return stratcomp;
}

/** stratcompatfullplusplus - determines whether a compatible character pair are compatible with stratigraphy or not.  
/* switched from unsigned long to double on 2011-10-19: that way, "half" can be credited for ties	*/
/* Requires:
/*	ranges - a matrix of first and last appearances where:
/*		ranges[x][0]=first appearance;
/*		ranges[x][1]=last appearanc;
/*	matrix - cladistic character matrix;
/*	states - number of states for ach character;
/*	charcomps - compatibility of each character;
/*	ch1 - 1st character being compared;
/*	ch2 - 2nd character being compared;
/*	notu - number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	stratcomp - an array telling you whether the character shows gaps in its  stratigraphic range or not.
		stratcomp[0]: character 1
		stratcomp[1]: character 2
		stratcomp[2]: character 1 states
		stratcomp[3]: character 2 states
		stratcomp[4]: state comparisons (2 for binary; possibly more for multistate)
		stratcomp[5]: general stratigraphic compatibility!
		stratcomp[6]: strict stratigraphic compatibility
		stratcomp[7]: super strict stratigraphic compatibility
		stratcomp[8]: divergent (0) vs. hierarchical (1) stratigraphic compatibility
		stratcomp[9]: hierarchical stratigraphic compatibility consistent with anagenesis
		stratcomp[10]: hierarchical stratigraphic compatibility consistent with budding
		stratcomp[11]: compatibility suggests 10<-01->11 rather than 10->01->11
		stratcomp[12]: compatibility suggests 10->01->11 rather than 10<-01->11
		stratcomp[13]: stratigraphy suggests 10<-01->11 rather than 10->01->11
		stratcomp[14]: stratigraphy suggests 10->01->11 rather than 10<-01->11
		stratcomp[15]: relative diversity of younger swing pair in SIN cases
***************************************************************************************************************************/
double *stratcompatfullplusplus(long **ranges, long **matrix, int *states, unsigned long *charcomps, int ch1, int ch2, int notu, int tiebreaker, int UNKNOWN, int INAP)
{
int	a, b, c, d, f, p;
int	maxocc, st1, st2, ttlprs=0, maxst;
//int	cc, sp, maxocc, st1, st2, sc1, sc2, rc, scf1, scf2, scu1, scu2, prs, m1, m2, ttlprs=0, maxst;
//int cl[20], cr[20];
int	swing, end1, end2, cpe1, cpe2;						/* these are the "swing" states between characters + important info, with "swing" giving pair #	*/
/*int	kp1, kp2;										/* key pairs			*/
/*int	prdv1, prdv2;									/* key pair diversities	*/
int cp1, cp2, cp3;
int epd1, epd2;											/* diversity of end pairs for incompatible pairs	*/
unsigned long **combos;
/*long **pairfnd;*/
long **pairs, **allpairs, **pairrng, **pairfnd;
long gap[3];
long prfl[3][2];
double *stratcomp;
int	hier, divr, s1, s2;
//long pairs[3][2];

/* 2011-10-19: NEW APPROACH!  
	1. Put all pairs into a pair x 2 matrix
	2. Put in a 3-tiered for loop to find all combinations of 3 pairs;
	3. If there is a swing-pair, then analyze this pair.
****************************************************/
maxocc=maxlmatrixcol(ranges, notu, 1);
maxst=imax(states[ch1],states[ch2]);
pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));

for (st1=0; st1<states[ch1]; ++st1)	{
	for (st2=0; st2<states[ch2]; ++st2)	{
		a=st1*states[ch2]+st2;
		/* find the number of combinations in which each states of both characters appears	*/
		if (pairrng[a][0]<=maxocc && pairrng[a][0]>=0)	{
			++combos[0][st1];
			++combos[1][st2];
			++ttlprs;			/* total pairs observed: if this is less than 3, then we can quit now!	*/
			}
		}
	}
free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));

stratcomp=dvector(16);
cleardvector(stratcomp,16,0);

stratcomp[0]=ch1;
stratcomp[1]=ch2;
stratcomp[2]=states[ch1];
stratcomp[3]=states[ch2];

if (ttlprs>2)	{
	allpairs=lmatrix(ttlprs,2);
	p=0;
	for (st1=0; st1<states[ch1]; ++st1)	{
		for (st2=0; st2<states[ch2]; ++st2)	{
			a=st1*states[ch2]+st2;
			/* find the number of combinations in which each states of both characters appears	*/
			if (pairrng[a][0]<=maxocc && pairrng[a][0]>=0)	{
				allpairs[p][0]=st1;
				allpairs[p][1]=st2;
				++p;			/* total pairs observed: if this is less than 3, then we can quit now!	*/
				}
			}
		}

	/* make sure that pairfnd uses (st1*states[ch2])+states[ch2] throughout, and not simply the pair number!	*/
	pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/

	/* now, go through all combinations of three pairs, and exampine those that have a linking pair	*/
	pairs=lmatrix(3,2);
	for (a=0; a<(ttlprs-2); ++a)	{
		pairs[0][0]=allpairs[a][0];
		pairs[0][1]=allpairs[a][1];
		/* get first-last appearances of pair a	*/
		prfl[0][0]=pairrng[cp1=(pairs[0][0]*states[ch2])+pairs[0][1]][0];
		prfl[0][1]=pairrng[cp1][1];
		/* determine if this pair has a gap in it	*/
		gap[0]=0;
		for (d=prfl[0][0]; d<=prfl[0][1]; ++d)	{
			if (pairfnd[cp1][d]==0)	{
				gap[0]=1;
				d=prfl[0][1];
				}
			}
		for (b=a+1; b<(ttlprs-1); ++b)	{
			pairs[1][0]=allpairs[b][0];
			pairs[1][1]=allpairs[b][1];
			/* get first-last appearances of pair b	*/
			prfl[1][0]=pairrng[cp2=(pairs[1][0]*states[ch2])+pairs[1][1]][0];
			prfl[1][1]=pairrng[cp2][1];
			/* determine if this pair has a gap in it	*/
			gap[1]=0;
			for (d=prfl[1][0]; d<=prfl[1][1]; ++d)	{
				if (pairfnd[cp2][d]==0)	{
					gap[1]=1;
					d=prfl[1][1];
					}
				}
			for (c=b+1; c<ttlprs; ++c)	{
				pairs[2][0]=allpairs[c][0];
				pairs[2][1]=allpairs[c][1];
				/* get first-last appearances of pair c	*/
				prfl[2][0]=pairrng[cp3=(pairs[2][0]*states[ch2])+pairs[2][1]][0];
				prfl[2][1]=pairrng[cp3][1];
				/* determine if this pair has a gap in it	*/
				gap[2]=0;
				for (d=prfl[2][0]; d<=prfl[2][1]; ++d)	{
					if (pairfnd[cp3][d]==0)	{
						gap[2]=1;
						d=prfl[2][1];
						}
					}
				/* determine whether this set is linked: if so, then swing will be positive	*/
				swing=findswingpair(pairs,3,maxst);
				if (swing>-1)	{
					stratcomp[4]+=1.0f;							/* another comparison made						*/
					end1=0;
					end2=2;
					if (swing==0)
						end1=1;
					else if (swing==2)
						end2=1;
					
					/* make end1 older than end2	*/
					if (prfl[end1][0]>prfl[end2][0])	{
						d=end1;
						end1=end2;
						end2=d;
						}
					/* insert something here so that pairfnd's first number points to ALL pairs, not just the three being considered	*/
					cpe1=(pairs[end1][0]*states[ch2])+pairs[end1][1];
					cpe2=(pairs[end2][0]*states[ch2])+pairs[end2][1];
/*					cpsw=(pairs[swing][0]*states[ch2])+pairs[swing][1];	/* this might not be needed	*/
					
					/* if swing pair appears after both end pairs appear, then this is stratigraphically incompatible	*/
					if ((prfl[swing][0]>prfl[end1][0]) && (prfl[swing][0]>prfl[end2][0]))	{
						/* WORK HERE 2011-10-20	*/
						
						d=prfl[swing][0];	/* first appearance of swing	*/
						/* if second appearing pair is more diverse, then assume that it is ancestral	*/
						epd1=pairfnd[cpe1][d-1]+pairfnd[cpe1][d];
						epd2=pairfnd[cpe2][d-1]+pairfnd[cpe2][d];
						if ((pairfnd[cpe1][d-1]+pairfnd[cpe1][d]) < (pairfnd[cpe2][d-1]+pairfnd[cpe2][d]))
							stratcomp[14]+=1.0f;				/* "hierarchical" stratigraphic incompatibility	*/
						else if ((pairfnd[cpe1][d-1]+pairfnd[cpe1][d]) == (pairfnd[cpe2][d-1]+pairfnd[cpe2][d]))	{
							stratcomp[13]+=0.5f;				/* split between "hierarchical" & "divergent" stratigraphic incompatibility	*/
							stratcomp[14]+=0.5f;				/* split between "hierarchical" & "divergent" stratigraphic incompatibility	*/
							}
						else
							stratcomp[13]+=1.0f;				/* "divergent" stratigraphic incompatibility	*/
						
						/* calculate relative diversity of younger morphotype & add it																*/
						/* if neither is present immediately before or after, then they had equal diversity											*/
						if ((epd1+epd2)>0)	stratcomp[15]+=((double) (pairfnd[cpe2][d-1]+pairfnd[cpe2][d]))/(((double) (pairfnd[cpe1][d-1]+pairfnd[cpe1][d]))+((double) (pairfnd[cpe2][d-1]+pairfnd[cpe2][d])));
						else				stratcomp[15]+=0.5;

						/* determine which character was more apt to change, and thus whether end1 or end2 is more apt to be ancestral	*/
						/* old routine	*/
//						m1=m2=RAND_MAX;
//						for (d=0; d<prs; ++d)	{
//							if (pairs[d][0]==keyst1 && prfl[d][0]<m1)	m1=prfl[d][0];
//							if (pairs[d][1]==keyst2 && prfl[d][0]<m2)	m2=prfl[d][0];
//							}
						/* if the older state is from the less compatible character, then increment	*/
//						if ((m1<m2 && charcomps[ch2]>charcomps[ch1]) || (m2<m1 && charcomps[ch1]>charcomps[ch2]))	++stratcomp[11];
//						else if (m2==m1 || charcomps[ch2]==charcomps[ch1])											++stratcomp[12];

						/* new routine	*/
						/* use the character with greater compatibility. If character 1, then assume that 00 evolved from 01, not 10	*/ 
						/* That means finding whether end1 or end2 matches swing for character 1										*/
						if (charcomps[ch1]>charcomps[ch2])	{
							/* assume that swing evolved from end2: 10->01->11	*/
							if (pairs[end1][0]==pairs[swing][0])
								stratcomp[11]+=1.0f;
							/* assume that swing evolved from end1: 01<-10->11	*/
							else
								stratcomp[12]+=1.0f;
							}
						else if (charcomps[ch1]<charcomps[ch2])	{
							/* assume that swing evolved from end1: 01<-10->11	*/
							if (pairs[end2][0]==pairs[swing][0])
								stratcomp[12]+=1.0f;
							/* assume that swing evolved from end2: 10->01->11	*/
							else
								stratcomp[11]+=1.0f;
							}
						/* split tie	*/
						else	{
							stratcomp[11]+=0.5f;
							stratcomp[12]+=0.5f;
							}
						
						}	/* end case of stratigraphic incompatibility	*/
						
					else	{
						stratcomp[5]+=1.0f;							/* general stratigraphic compatibility						*/
						
						/* divergent vs. hierarchical compatibility	*/
//						if (prfl[swing][0]>prfl[end1][0] && prfl[swing][0]<=prfl[end2][0])	{
						hier=divr=0;
						if (prfl[swing][0]>prfl[end1][0])		hier=1;
						else if (prfl[swing][0]<prfl[end1][0])	divr=1;
						else if (tiebreaker==1)	{
							/* routine if 2 appear at the same time	*/
							if (prfl[end2][0]>prfl[swing][0])	{
								s1=(pairs[end1][0]*states[ch2])+pairs[end1][1];		/* 2013-04-12: get number of state pair for ranges */
								s2=(pairs[swing][0]*states[ch2])+pairs[swing][1];	/* 2013-04-16: get number of state pair for ranges */
								if (pairfnd[s1][prfl[swing][0]] > pairfnd[s2][prfl[end1][0]])
									hier=1;
								else if (pairfnd[s1][prfl[swing][0]] < pairfnd[s2][prfl[end1][0]])
									divr=1;
								}	/* end test to see whether one of the two oldest pairs is more diverse than the other	*/
							}

						/* routine for clearly hierarchical	*/
//						if (prfl[swing][0]>prfl[end1][0])	{
						if (hier==1)	{
							/* second condition should not be necessary	*/
							/*  e1	sw	e2
								•	•	•
								•	•---•
								•---•
								•			*/
							stratcomp[8]+=1.0f;						/*  tally hierarchical compatilibility						*/
							
							/* ask whether hierarchical compatibility is consistent with anagenesis or budding					*/
							if (prfl[swing][0]>=prfl[end1][1])
								stratcomp[9]+=0.5f;					/* consistent with anagenesis								*/
							else
								stratcomp[10]+=0.5f;				/* consistent with budding									*/
							if (prfl[end2][0]>=prfl[swing][1])
								stratcomp[9]+=0.5f;					/* consistent with anagenesis								*/
							else
								stratcomp[10]+=0.5f;				/* consistent with budding									*/

							/* end1 should overlap with swing and swing should overlap with end2		*/ 
							if (((prfl[end1][1]+1)>=prfl[swing][0]) && ((prfl[swing][1]+1)>=prfl[end2][0]))	{
								stratcomp[6]+=1.0f;					/* no gap between morphotypes								*/
								/* super-strict requires no gaps at all!														*/
								f=sumlvector(gap,3);
								if (f==0)		stratcomp[7]+=1.0f;	/* no gaps within or between morphotypes					*/ 
								}
							}	/* end routine for clearly hierarchical															*/

						/* routine for clearly divergent	*/
//						else if (prfl[swing][0]<prfl[end1][0])	{
						else if (divr==1)	{
							/*  e1	sw	e2
								•	•	•
								•	•---•
								•---•
									•		*/
							/* strict stratigraphic compatibility achieved if swing state overlaps with both derived states	*/
							if (((prfl[swing][1]+1)>=prfl[end1][0]) && ((prfl[swing][1]+1)>=prfl[end2][0]))	{
								stratcomp[6]+=1.0f;					/* no gap between morphotypes								*/
								/* super-strict requires no gaps at all!														*/
								f=sumlvector(gap,3);
								if (f==0)	stratcomp[7]+=1.0f;				/* no gaps within or between morphotypes					*/ 
								}	/* end search for strict stratigraphic compatibiilty	*/
							}	/* end check of strict and superstrict stratigraphic compatibility for divergent stratigraphic compatibility	*/

						/* if the first two morphotypes appear at the same time, then it as a split	*/
//						else if (prfl[swing][0]==prfl[end1][0])	{
						else if (hier==0 && divr==0)	{
							/*  e1	sw	e2
								•	•	•
								•	•---•
								•   •
								•-?-•	*/
							/* if end1 or swing disappears before end2 appears OR if end1 or swing disappears first when both disappear before end1, then assume it is not ancestral	*/
/*							if ((prfl[end1][1]+1)<prfl[end2][0] && (prfl[swing][1]+1)>=prfl[end2][0])	*/

							/* otherwise, assume that this could be either hierarchical or divergent									*/
							stratcomp[8]+=0.5f;

							/* now worry about budding and anagenesis.  First, it's 50% that it's divergent, so all values are chopped in half			*/
							/*		2nd, the first two pairs can only fit budding as they are contemporaneous											*/
							stratcomp[10]+=0.25f;					/* end1 & swing consistent with budding IF the pair is hierarchical		    		*/
							/* ask whether possible hierarchical compatibility is consistent with anagenesis or budding									*/
							if (prfl[end2][0]>=prfl[swing][1])
								stratcomp[9]+=0.25f;				/* consistent with anagenesis: but only half because it might not be hierarchical	*/
							else
								stratcomp[10]+=0.25f;				/* one pair definitely consistent with budding, other might be		    			*/
							
							/* we know that end1 and swing overlap; the quesiton is whether swing overlaps with end2	*/
							if ((prfl[swing][1]+1)>=prfl[end2][0])	{
								stratcomp[6]+=1.0f;					/* no gap between morphotypes								*/
								/* super-strict requires no gaps at all!														*/
								f=sumlvector(gap,3);
								if (f==0)		stratcomp[7]+=1.0f;	/* no gaps within or between morphotypes					*/ 
								}
							}	/* end routine for ambivalent															*/
						/* WORK HERE 2011-10-20	*/
						}	/* end case of stratigraphic compatibility	*/
					}	/* end case of compatible trio of character pairs	*/
				}	/* end possible 3rd pairs	*/
			}	/* end possible 2nd pairs	*/
		}	/* end possible 1st pairs	*/
	
	free_lmatrix(pairs,3,2);
	free_lmatrix(allpairs,ttlprs,2);
	free_lmatrix(pairfnd,states[ch1]*states[ch2],1+maxocc);
	}	/* end cases of 2+ pairs of character states	*/

free_lmatrix(pairrng,states[ch1]*states[ch2],2);

return stratcomp;
}


/** stratcompatfullplusplusplus - determines whether a compatible character pair are compatible with stratigraphy or not.  
/* switched from unsigned long to double on 2011-10-19: that way, "half" can be credited for ties	*/
/* Requires:
/*	ranges - a matrix of first and last appearances where:
/*		ranges[x][0]=first appearance;
/*		ranges[x][1]=last appearanc;
/*	matrix - cladistic character matrix;
/*	states - number of states for ach character;
/*	charcomps - compatibility of each character;
/*	ch1 - 1st character being compared;
/*	ch2 - 2nd character being compared;
/*	notu - number of taxa;
/*	onset - when clade starts;
/*	term - when clade ends;
/*  tiebreak - 0 means ignore ties, 1 means use earliest richness as a tie-breaker
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	stratcomp - an array telling you whether the character shows gaps in its  stratigraphic range or not.
		stratcomp[0]: character 1
		stratcomp[1]: character 2
		stratcomp[2]: character 1 states
		stratcomp[3]: character 2 states
		stratcomp[4]: state comparisons (2 for binary; possibly more for multistate)
		stratcomp[5]: general stratigraphic compatibility!
		stratcomp[6]: strict stratigraphic compatibility
		stratcomp[7]: super strict stratigraphic compatibility
		stratcomp[8]: divergent (0) vs. hierarchical (1) stratigraphic compatibility
		stratcomp[9]: hierarchical stratigraphic compatibility consistent with anagenesis
		stratcomp[10]: hierarchical stratigraphic compatibility consistent with budding
		stratcomp[11]: compatibility suggests 10<-01->11 rather than 10->01->11
		stratcomp[12]: compatibility suggests 10->01->11 rather than 10<-01->11
		stratcomp[13]: stratigraphy suggests 10<-01->11 rather than 10->01->11
		stratcomp[14]: stratigraphy suggests 10->01->11 rather than 10<-01->11
		stratcomp[15]: relative diversity of younger swing pair in SIN cases
		stratcomp[16]: center of gravity for oldest state pair
		stratcomp[17]: center of gravity for intermediate state pair
		stratcomp[18]: when hierarchical ancestral pair goes extinct;
		stratcomp[19]: when divergent ancestral pair goes extinct;
		stratcomp[20]: purely hierarchical comparisons
		stratcomp[21]: purely divergent comparisons
		
		MAKE SURE PROPORTION OF "LIVING FOSSILS" RECORDED
****************************************************************************************************************************/
double *stratcompatfullplusplusplus(long **ranges, long **matrix, int *states, unsigned long *charcomps, int ch1, int ch2, int notu, int onset, int term, int tiebreaker, int UNKNOWN, int INAP)
{
int	a, b, c, d, f, p, t;
int	maxocc, st1, st2, ttlprs=0, maxst;
//int	cc, sp, maxocc, st1, st2, sc1, sc2, rc, scf1, scf2, scu1, scu2, prs, m1, m2, ttlprs=0, maxst;
//int cl[20], cr[20];
int	swing, end1, end2, cpe1, cpe2;						/* these are the "swing" states between characters + important info, with "swing" giving pair #	*/
/*int	kp1, kp2;										/* key pairs			*/
/*int	prdv1, prdv2;									/* key pair diversities	*/
int cp1, cp2, cp3;
int epd1, epd2;											/* diversity of end pairs for incompatible pairs	*/
unsigned long **combos;
/*long **pairfnd;*/
long **pairs, **allpairs, **pairrng, **pairfnd;
long gap[3];
long prfl[3][2];
double *stratcomp;
double cgan, cgcl;										/* 2013-04-12: center of gravity	*/
int ancp, desc1, desc2;									/* 2013-04-12: ancestral pair		*/
int hier, divr;
int trioonset, trioterm;
long *S;
//long pairs[3][2];

/* 2011-10-19: NEW APPROACH!  
	1. Put all pairs into a pair x 2 matrix
	2. Put in a 3-tiered for loop to find all combinations of 3 pairs;
	3. If there is a swing-pair, then analyze this pair.
****************************************************/
maxocc=maxlmatrixcol(ranges, notu, 1);
maxst=imax(states[ch1],states[ch2]);
pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));

/* added again on 2013-04-12: find beginning and end of clade	*/
//onset=minlmatrixcol(ranges,notu,0);
//end=maxlmatrixcol(ranges,notu,0);

for (st1=0; st1<states[ch1]; ++st1)	{
	for (st2=0; st2<states[ch2]; ++st2)	{
		a=st1*states[ch2]+st2;
		/* find the number of combinations in which each states of both characters appears	*/
		if (pairrng[a][0]<=maxocc && pairrng[a][0]>=0)	{
			++combos[0][st1];
			++combos[1][st2];
			++ttlprs;			/* total pairs observed: if this is less than 3, then we can quit now!	*/
			}
		}
	}
free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));

stratcomp=dvector(22);
cleardvector(stratcomp,22,0);

stratcomp[0]=ch1;
stratcomp[1]=ch2;
stratcomp[2]=states[ch1];
stratcomp[3]=states[ch2];

S=lvector(term+1);

if (ttlprs>2)	{
	allpairs=lmatrix(ttlprs,2);
	p=0;
	for (st1=0; st1<states[ch1]; ++st1)	{
		for (st2=0; st2<states[ch2]; ++st2)	{
			a=st1*states[ch2]+st2;
			/* find the number of combinations in which each states of both characters appears	*/
			if (pairrng[a][0]<=maxocc && pairrng[a][0]>=0)	{
				allpairs[p][0]=st1;
				allpairs[p][1]=st2;
				++p;			/* total pairs observed: if this is less than 3, then we can quit now!	*/
				}
			}
		}

	/* make sure that pairfnd uses (st1*states[ch2])+states[ch2] throughout, and not simply the pair number!	*/
	pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/

	/* now, go through all combinations of three pairs, and exampine those that have a linking pair	*/
	pairs=lmatrix(3,2);
	for (a=0; a<(ttlprs-2); ++a)	{
		pairs[0][0]=allpairs[a][0];
		pairs[0][1]=allpairs[a][1];
		/* get first-last appearances of pair a	*/
		prfl[0][0]=pairrng[cp1=(pairs[0][0]*states[ch2])+pairs[0][1]][0];
		prfl[0][1]=pairrng[cp1][1];
		/* determine if this pair has a gap in it	*/
		gap[0]=0;
		for (d=prfl[0][0]; d<=prfl[0][1]; ++d)	{
			if (pairfnd[cp1][d]==0)	{
				gap[0]=1;
				d=prfl[0][1];
				}
			}
		for (b=a+1; b<(ttlprs-1); ++b)	{
			pairs[1][0]=allpairs[b][0];
			pairs[1][1]=allpairs[b][1];
			/* get first-last appearances of pair b	*/
			prfl[1][0]=pairrng[cp2=(pairs[1][0]*states[ch2])+pairs[1][1]][0];
			prfl[1][1]=pairrng[cp2][1];
			/* determine if this pair has a gap in it	*/
			gap[1]=0;
			for (d=prfl[1][0]; d<=prfl[1][1]; ++d)	{
				if (pairfnd[cp2][d]==0)	{
					gap[1]=1;
					d=prfl[1][1];
					}	/* gap found	*/
				}
			for (c=b+1; c<ttlprs; ++c)	{
				pairs[2][0]=allpairs[c][0];
				pairs[2][1]=allpairs[c][1];
				/* get first-last appearances of pair c	*/
				prfl[2][0]=pairrng[cp3=(pairs[2][0]*states[ch2])+pairs[2][1]][0];
				prfl[2][1]=pairrng[cp3][1];
				/* determine if this pair has a gap in it	*/
				gap[2]=0;
				for (d=prfl[2][0]; d<=prfl[2][1]; ++d)	{
					if (pairfnd[cp3][d]==0)	{
						gap[2]=1;
						d=prfl[2][1];
						}
					}
				/* determine whether this set is linked: if so, then swing will be positive	*/
				swing=findswingpair(pairs,3,maxst);
				if (swing>-1)	{
					stratcomp[4]+=1.0f;							/* another comparison made						*/
					end1=0;
					end2=2;
					if (swing==0)
						end1=1;
					else if (swing==2)
						end2=1;
					
					/* make end1 older than end2	*/
					if (prfl[end1][0]>prfl[end2][0])	{
						d=end1;
						end1=end2;
						end2=d;
						}
					/* insert something here so that pairfnd's first number points to ALL pairs, not just the three being considered	*/
					cpe1=(pairs[end1][0]*states[ch2])+pairs[end1][1];
					cpe2=(pairs[end2][0]*states[ch2])+pairs[end2][1];
/*					cpsw=(pairs[swing][0]*states[ch2])+pairs[swing][1];	/* this might not be needed	*/
					
					/* if swing pair appears after both end pairs appear, then this is stratigraphically incompatible	*/
					if ((prfl[swing][0]>prfl[end1][0]) && (prfl[swing][0]>prfl[end2][0]))	{
						/* WORK HERE 2011-10-20	*/
						
						d=prfl[swing][0];	/* first appearance of swing	*/
						/* if second appearing pair is more diverse, then assume that it is ancestral	*/
						epd1=pairfnd[cpe1][d-1]+pairfnd[cpe1][d];
						epd2=pairfnd[cpe2][d-1]+pairfnd[cpe2][d];
						if ((pairfnd[cpe1][d-1]+pairfnd[cpe1][d]) < (pairfnd[cpe2][d-1]+pairfnd[cpe2][d]))
							stratcomp[14]+=1.0f;				/* "hierarchical" stratigraphic incompatibility	*/
						else if ((pairfnd[cpe1][d-1]+pairfnd[cpe1][d]) == (pairfnd[cpe2][d-1]+pairfnd[cpe2][d]))	{
							stratcomp[13]+=0.5f;				/* split between "hierarchical" & "divergent" stratigraphic incompatibility	*/
							stratcomp[14]+=0.5f;				/* split between "hierarchical" & "divergent" stratigraphic incompatibility	*/
							}
						else
							stratcomp[13]+=1.0f;				/* "divergent" stratigraphic incompatibility	*/
						
						/* calculate relative diversity of younger morphotype & add it																*/
						/* if neither is present immediately before or after, then they had equal diversity											*/
						if ((epd1+epd2)>0)	stratcomp[15]+=((double) (pairfnd[cpe2][d-1]+pairfnd[cpe2][d]))/(((double) (pairfnd[cpe1][d-1]+pairfnd[cpe1][d]))+((double) (pairfnd[cpe2][d-1]+pairfnd[cpe2][d])));
						else				stratcomp[15]+=0.5;

						/* determine which character was more apt to change, and thus whether end1 or end2 is more apt to be ancestral	*/
						/* old routine	*/
//						m1=m2=RAND_MAX;
//						for (d=0; d<prs; ++d)	{
//							if (pairs[d][0]==keyst1 && prfl[d][0]<m1)	m1=prfl[d][0];
//							if (pairs[d][1]==keyst2 && prfl[d][0]<m2)	m2=prfl[d][0];
//							}
						/* if the older state is from the less compatible character, then increment	*/
//						if ((m1<m2 && charcomps[ch2]>charcomps[ch1]) || (m2<m1 && charcomps[ch1]>charcomps[ch2]))	++stratcomp[11];
//						else if (m2==m1 || charcomps[ch2]==charcomps[ch1])											++stratcomp[12];

						/* new routine	*/
						/* use the character with greater compatibility. If character 1, then assume that 00 evolved from 01, not 10	*/ 
						/* That means finding whether end1 or end2 matches swing for character 1										*/
						if (charcomps[ch1]>charcomps[ch2])	{
							/* assume that swing evolved from end2: 10->01->11	*/
							if (pairs[end1][0]==pairs[swing][0])
								stratcomp[11]+=1.0f;
							/* assume that swing evolved from end1: 01<-10->11	*/
							else
								stratcomp[12]+=1.0f;
							}
						else if (charcomps[ch1]<charcomps[ch2])	{
							/* assume that swing evolved from end1: 01<-10->11	*/
							if (pairs[end2][0]==pairs[swing][0])
								stratcomp[12]+=1.0f;
							/* assume that swing evolved from end2: 10->01->11	*/
							else
								stratcomp[11]+=1.0f;
							}
						/* split tie	*/
						else	{
							stratcomp[11]+=0.5f;
							stratcomp[12]+=0.5f;
							}
						
						}	/* end case of stratigraphic incompatibility	*/
						
					else	{
						stratcomp[5]+=1.0f;							/* general stratigraphic compatibility						*/
						
						/* divergent vs. hierarchical compatibility	*/
//						if (prfl[swing][0]>prfl[end1][0] && prfl[swing][0]<=prfl[end2][0])	{
						hier=divr=0;
						if (prfl[swing][0]>prfl[end1][0])		hier=1;
						else if (prfl[swing][0]<prfl[end1][0])	divr=1;
						else if (tiebreaker==1)	{
							/* routine if 2 appear at the same time	*/
							if (prfl[end2][0]>prfl[swing][0])	{
								desc1=(pairs[end1][0]*states[ch2])+pairs[end1][1];		/* 2013-04-12: get number of state pair for ranges */
								desc2=(pairs[swing][0]*states[ch2])+pairs[swing][1];	/* 2013-04-16: get number of state pair for ranges */
								if (pairfnd[desc1][prfl[swing][0]] > pairfnd[desc2][prfl[end1][0]])
									hier=1;
								else if (pairfnd[desc1][prfl[swing][0]] < pairfnd[desc2][prfl[end1][0]])
									divr=1;
								}	/* end test to see whether one of the two oldest pairs is more diverse than the other	*/
							}

						/* routine for clearly hierarchical	*/
//						if (prfl[swing][0]>prfl[end1][0])	{
						if (hier==1)	{
							/*  e1	sw	e2
								•	•	•
								•	•---•
								•---•
								•			*/
							stratcomp[8]+=1.0f;						/*  tally hierarchical compatilibility						*/
							stratcomp[20]+=1.0f;					/*  2013-04-13: purely hierarchical compatilibility						*/

							/* 2013-04-12: AGAIN (and hope it sticks): find center of gravity for original state pair			*/
							ancp=(pairs[end1][0]*states[ch2])+pairs[end1][1];		/* 2013-04-12: get number of state pair for ranges */
							desc1=(pairs[swing][0]*states[ch2])+pairs[swing][1];	/* 2013-04-16: get number of state pair for ranges */
							desc2=(pairs[end2][0]*states[ch2])+pairs[end2][1];		/* 2013-04-16: get number of state pair for ranges */
							/* 2013-04-16: get clade shape for just these three state-pairs	*/
							trioterm=trioonset=prfl[end1][0];
							for (t=0; t<3; ++t)	if (prfl[t][1]>trioterm)	trioterm=imin(prfl[t][1],term);
							for (t=onset; t<=term; ++t)	S[t]=pairfnd[ancp][t]+pairfnd[desc1][t]+pairfnd[desc2][t];
							
							cgcl=centerofgravityfromrichness(S,trioonset,trioterm);					/* CoG for whole clade, time standardized	*/
							cgan=centerofgravityfromrichness(pairfnd[ancp],trioonset,trioterm);		/* CoG for ancestral paraclade, time standardized	*/
//							stratcomp[16]=cg[1]/ccg[0];
//							stratcomp[16]=cg[1]/ccg[1];
							/* in the end, I decided to use the time standardized equation and make it a proportion of the total clade's time-standardized CoG	*/
							/* this should account for differences in clade durations as well as multistate trios taking up only part of the clade's history	*/
//							stratcomp[16]=(((1+trioterm-trioonset) - (cg[1]/ccg[0]))/(1+trioterm-trioonset)) / (((1+trioterm-trioonset) - (ccg[1]/ccg[0]))/(1+trioterm-trioonset));
							stratcomp[16]=(1-cgan)/(1-cgcl);
							/************ end clade shape part ************/
							
							/* record how long ancestral state pair survives	*/
							t=imin(prfl[end1][1],term);							/* end of ancestral pair: but do not let it exceed clade age	*/
							stratcomp[18]=((double) (1+(t-trioonset)))/((double) (1+(trioterm-trioonset)));

							/* ask whether hierarchical compatibility is consistent with anagenesis or budding					*/
							if (prfl[swing][0]>=prfl[end1][1])
								stratcomp[9]+=0.5f;					/* consistent with anagenesis								*/
							else
								stratcomp[10]+=0.5f;				/* consistent with budding									*/
							if (prfl[end2][0]>=prfl[swing][1])
								stratcomp[9]+=0.5f;					/* consistent with anagenesis								*/
							else
								stratcomp[10]+=0.5f;				/* consistent with budding									*/

							/* end1 should overlap with swing and swing should overlap with end2		*/ 
							if (((prfl[end1][1]+1)>=prfl[swing][0]) && ((prfl[swing][1]+1)>=prfl[end2][0]))	{
								stratcomp[6]+=1.0f;					/* no gap between morphotypes								*/
								/* super-strict requires no gaps at all!														*/
								f=sumlvector(gap,3);
								if (f==0)		stratcomp[7]+=1.0f;	/* no gaps within or between morphotypes					*/ 
								}
							}	/* end routine for clearly hierarchical															*/

						/* routine for clearly divergent	*/
//						else if (prfl[swing][0]<prfl[end1][0])	{
						else if (divr==1)	{
							/*  e1	sw	e2
								•	•	•
								•	•---•
								•---•
									•		*/
							stratcomp[21]+=1.0f;									/* 2013-04-15: purely divergent compatilibility						*/

							/* ask whether hierarchical compatibility is consistent with anagenesis or budding					*/
							/* 2013-04-12: AGAIN (and hope it sticks): find center of gravity for original state pair			*/
							ancp=(pairs[swing][0]*states[ch2])+pairs[swing][1];	/* 2013-04-16: get number of state pair for ranges */
							desc1=(pairs[end1][0]*states[ch2])+pairs[end1][1];		/* 2013-04-12: get number of state pair for ranges */
							desc2=(pairs[end2][0]*states[ch2])+pairs[end2][1];		/* 2013-04-16: get number of state pair for ranges */
							/* 2013-04-16: get clade shape for just these three state-pairs	*/
							trioterm=trioonset=prfl[swing][0];
							for (t=0; t<3; ++t)	if (prfl[t][1]>trioterm)	trioterm=imin(prfl[t][1],term);
							for (t=onset; t<=term; ++t)	S[t]=pairfnd[ancp][t]+pairfnd[desc1][t]+pairfnd[desc2][t];

							cgcl=centerofgravityfromrichness(S,trioonset,trioterm);					/* CoG for whole clade, time standardized	*/
							cgan=centerofgravityfromrichness(pairfnd[ancp],trioonset,trioterm);		/* CoG for ancestral paraclade, time standardized	*/
							stratcomp[17]=(1-cgan)/(1-cgcl);
							
							/* record how long ancestral state pair survives	*/
							t=imin(prfl[swing][1],term);							/* end of ancestral pair: but do not let it exceed clade age	*/
							stratcomp[19]=((double) (1+(t-trioonset)))/((double) (1+(trioterm-trioonset)));

							/* strict stratigraphic compatibility achieved if swing state overlaps with both derived states	*/
							if (((prfl[swing][1]+1)>=prfl[end1][0]) && ((prfl[swing][1]+1)>=prfl[end2][0]))	{
								stratcomp[6]+=1.0f;					/* no gap between morphotypes								*/
								/* super-strict requires no gaps at all!														*/
								f=sumlvector(gap,3);
								if (f==0)	stratcomp[7]+=1.0f;				/* no gaps within or between morphotypes					*/ 
								}	/* end search for strict stratigraphic compatibiilty	*/
							}	/* end check of strict and superstrict stratigraphic compatibility for divergent stratigraphic compatibility	*/

						/* if the first two morphotypes appear at the same time, then it as a split	*/
//						else if (prfl[swing][0]==prfl[end1][0])	{
						else if (hier==0 && divr==0)	{
							/*  e1	sw	e2
								•	•	•
								•	•---•
								•   •
								•-?-•	*/
							/* if end1 or swing disappears before end2 appears OR if end1 or swing disappears first when both disappear before end1, then assume it is not ancestral	*/
/*							if ((prfl[end1][1]+1)<prfl[end2][0] && (prfl[swing][1]+1)>=prfl[end2][0])	*/

							/* otherwise, assume that this could be either hierarchical or divergent									*/
							stratcomp[8]+=0.5f;
						/* skip CoG stuff for 
							/* now worry about budding and anagenesis.  First, it's 50% that it's divergent, so all values are chopped in half			*/
							/*		2nd, the first two pairs can only fit budding as they are contemporaneous											*/
							stratcomp[10]+=0.25f;					/* end1 & swing consistent with budding IF the pair is hierarchical		    		*/
							/* ask whether possible hierarchical compatibility is consistent with anagenesis or budding									*/
							if (prfl[end2][0]>=prfl[swing][1])
								stratcomp[9]+=0.25f;				/* consistent with anagenesis: but only half because it might not be hierarchical	*/
							else
								stratcomp[10]+=0.25f;				/* one pair definitely consistent with budding, other might be		    			*/
							
							/* we know that end1 and swing overlap; the quesiton is whether swing overlaps with end2	*/
							if ((prfl[swing][1]+1)>=prfl[end2][0])	{
								stratcomp[6]+=1.0f;					/* no gap between morphotypes								*/
								/* super-strict requires no gaps at all!														*/
								f=sumlvector(gap,3);
								if (f==0)		stratcomp[7]+=1.0f;	/* no gaps within or between morphotypes					*/ 
								}
							}	/* end routine for ambivalent															*/
						/* WORK HERE 2011-10-20	*/
						}	/* end case of stratigraphic compatibility	*/
					}	/* end case of compatible trio of character pairs	*/
				}	/* end possible 3rd pairs	*/
			}	/* end possible 2nd pairs	*/
		}	/* end possible 1st pairs	*/
	
	free_lmatrix(pairs,3,2);
	free_lmatrix(allpairs,ttlprs,2);
	free_lmatrix(pairfnd,states[ch1]*states[ch2],1+maxocc);
	}	/* end cases of 2+ pairs of character states	*/

free_lvector(S);
free_lmatrix(pairrng,states[ch1]*states[ch2],2);

return stratcomp;
}



/** statepairranges - determines the first and last appearance of each pair.  
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

/** statepairfinds - determines whether each pair is found in each bin.  
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
/*	pairfnd - number of taxa with state pair found in each bin
****************************************************************************/
long **statepairfinds(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int  a, st1, st2, sp, st, pairs;
long **pairfnd;

//pairs=1+(states[ch1]*states[ch2])+states[ch2];
pairs=1+(states[ch1]*states[ch2]);
a=maxlmatrixcol(ranges,notu,1);
pairfnd=lmatrix(pairs,1+a);

for (sp=0; sp<notu; ++sp)	{
	if (((st1=matrix[sp][ch1])!=UNKNOWN && (st2=matrix[sp][ch2])!=UNKNOWN) && (matrix[sp][ch1]!=INAP && matrix[sp][ch2]!=INAP))	{
		a=(st1*states[ch2])+st2;			/* for 2 binary: 0<-00; 1<-01; 2<-10; 3<-11	*/
		for (st=ranges[sp][0]; st<=ranges[sp][1]; ++st)
/*			if (pairfnd[a][st]==0)			pairfnd[a][st]=1;	*/
			++pairfnd[a][st];
		}
	}

return pairfnd;
}


/** stateranges - determines the first and last bin for each state.  
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
int  a, s, sp;
long **staterng;

staterng=lmatrix(states,2);

for (s=0; s<states; ++s)	{
	staterng[s][0]=RAND_MAX;
	staterng[s][1]=-1*RAND_MAX;
	}

for (sp=0; sp<notu; ++sp)	{
	if ((s=matrix[sp][ch])!=UNKNOWN && matrix[sp][ch]!=INAP)	{
		if (s<0)	{
			s/=-1;
			while (s>0)	{
				a=s%10;
				if (ranges[sp][0]<staterng[a][0])	staterng[a][0]=ranges[sp][0];
				if (ranges[sp][1]>staterng[a][1])	staterng[a][1]=ranges[sp][1];
				s/=10;
				}
			}
		else	{
			if (ranges[sp][0]<staterng[s][0])	staterng[s][0]=ranges[sp][0];
			if (ranges[sp][1]>staterng[s][1])	staterng[s][1]=ranges[sp][1];
			}
		}
	}
return	staterng;
}

/** statefinds - determines whether each state occurs in each bin.  
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
int a, b, s, sp, st;
long **stfnd;

/*stfnd=lmatrix(states,1+(a=maxlmatrixcol(ranges,notu,1)));		sometimes state numbers >= states	*/
stfnd=lmatrix(1+(b=maxlmatrixcol(matrix,notu,ch)),1+(a=maxlmatrixcol(ranges,notu,1)));

a=states;
for (sp=0; sp<notu; ++sp)	{
	if ((s=matrix[sp][ch])!=UNKNOWN && matrix[sp][ch]!=INAP)	{
		/* if polymorphic	*/
		if (s<0)	{
			s/=-1;
			while (s>0)	{
				a=s%10;
				for (st=ranges[sp][0]; st<=ranges[sp][1]; ++st)
					if (stfnd[a][st]==0)			stfnd[a][st]=1;
				s/=10;
				}
			}
		/* if not polymorphic	*/
		else	{
			for (st=ranges[sp][0]; st<=ranges[sp][1]; ++st)
				if (stfnd[s][st]==0)			stfnd[s][st]=1;
			}
		}
	}

return stfnd;
}


/** stratconsiststate - determines whether all state in a character have continuous ranges.  
/* Requires:
/*	staterng - ranges of each state:
/*		staterng[x][0]=first appearance of state x;
/*		staterng[x][1]=last appearance of state x;
/*	stfnd - whether a state is found in each bin
/*		stfnd[x][y] = 1 if state x is found in bin y, 0 if absent;
/*	s - the state being examined;
/*
/* Returns:
/*	scs - 0 if not continuous, 1 if continuous
****************************************************************************/
unsigned long stratconsiststate(long **staterng, long **stfnd, int s)
{
int a, st;
unsigned long scs=1;

for (st=staterng[s][0]+1; st<staterng[s][1] && scs==1; ++st)	{
	a=stfnd[s][st];
	if (stfnd[s][st]==0)	scs=0;
	}
return scs;
}


/** stratconsistchar - determines whether all state in a character have continuous ranges.  
/* Requires:
/*	staterng - ranges of each state:
/*		staterng[x][0]=first appearance of state x;
/*		staterng[x][1]=last appearance of state x;
/*	states - the number fo states for the character;
/*
/* Returns:
/*	scc - 0 if not all are continuous, 1 if all are continuous
****************************************************************************/
unsigned long stratconsistchar(long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP)
{
int	a, s;
unsigned long scc=1;
long **staterng, **stfnd;

stfnd=statefinds(ranges, matrix, states, ch, notu, UNKNOWN, INAP);		/* pres/abs of states within strat bins	*/
staterng=stateranges(ranges, matrix, states, ch, notu, UNKNOWN, INAP);	/* find ranges of individual states	*/

for (s=0; s<states && scc==1; ++s)	{
	a=stratconsiststate(staterng, stfnd, s);
	if (a==0)	scc=0;
	}
free_lmatrix(staterng,states,2);
free_lmatrix(stfnd,states,2);
return scc;
}


/** stratconsiststpair - determines whether a state pair for two characters has a continuous ranges.  
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

/** stratconsistchpair - determines whether all state pairs for two characters have continuous ranges.  
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
unsigned long stratconsistchpair(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int pr, cs1, cs2;
unsigned long scp=1;
long max;
long **pairrng, **pairfnd;

pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/
max=maxlmatrixcol(ranges, notu, 1);

for (cs1=0; cs1<states[ch1] && scp==1; ++cs1)	{
	for (cs2=0; cs2<states[ch2] && scp==1; ++cs2)	{
		pr=(cs1*states[ch2])+cs2;
		if (pairrng[pr][0]>=0 && pairrng[pr][0]<=max)	
			scp=stratconsiststpair(pairrng, pairfnd,pr);
		}
	}
free_lmatrix(pairrng,states[ch1]*states[ch2],2);
free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
return scp;
}


/** stratfrconsistchpair - determines the proportion of state pairs for two characters have continuous ranges.  
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
double stratfrconsistchpair(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int a=0, pr, cs1, cs2;
long max;
long **pairrng, **pairfnd;
double y=0.0f, scp=0.0f;

pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/
max=maxlmatrixcol(ranges, notu, 1);

for (cs1=0; cs1<states[ch1]; ++cs1)	{
	for (cs2=0; cs2<states[ch2]; ++cs2)	{
		pr=(cs1*states[ch2])+cs2;
		if (pairrng[pr][0]>=0 && pairrng[pr][0]<=max)	{
			a=a+stratconsiststpair(pairrng, pairfnd,pr);
			++y;
			}
		}
	}
scp=((double) a)/y;

free_lmatrix(pairrng,states[ch1]*states[ch2],2);
free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
return scp;
}

/** stratconsistchar - determines the proportion of  all states in a character have continuous ranges.  
/* Requires:
/*	staterng - ranges of each state:
/*		staterng[x][0]=first appearance of state x;
/*		staterng[x][1]=last appearance of state x;
/*	states - the number fo states for the character;
/*
/* Returns:
/*	scc - 0 if not all are continuous, 1 if all are continuous
****************************************************************************/
double stratfrconsistchar(long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP)
{
int	a=0, s;
double scc=0.0f;
long **staterng, **stfnd;

stfnd=statefinds(ranges, matrix, states, ch, notu, UNKNOWN, INAP);		/* pres/abs of states within strat bins	*/
staterng=stateranges(ranges, matrix, states, ch, notu, UNKNOWN, INAP);	/* find ranges of individual states	*/

for (s=0; s<states; ++s)	{
	a=a+stratconsiststate(staterng, stfnd, s);
	scc+=a;
	}

scc=((double) a)/((double) states);
free_lmatrix(staterng,states,2);
free_lmatrix(stfnd,states,2);
return scc;
}

/* cleancladerangedata: remove gaps from a clade's range
****************************************************************************/
void cleancladerangedata(long **ranges, int notu)
{
int	b, c, sp;
int	lfa, lla;
//int good;
long	*rich;

lfa=maxlmatrixcol(ranges,notu,0);
lla=maxlmatrixcol(ranges,notu,1);
rich=richnesstally(ranges,notu);		/* Sepkoski style synoptic richness per bin	*/
/* 2013-04-22: weird memory error appearing here	*/

for (c=0; c<=lla; ++c)	{
//	good=0;
/*	for (sp=0; sp<notu; ++sp)	{
		if (ranges[sp][0]<=c && ranges[sp][1]>=c)	{
			good=1;
			sp=notu;
			}
		}			kluge to avoid memory problem: but it comes back later, anyway	*/
	if (rich[c]==0)	{
//	if (good==0)	{
		for (sp=0; sp<notu; ++sp)	{
			if (ranges[sp][0]>c)	{
				--ranges[sp][0];
				--ranges[sp][1];
				}
			else if (ranges[sp][1]>c)
				--ranges[sp][1];
			}
		for (b=c; b<lla; ++b)	rich[b]=rich[b+1];
		--c;
		--lla;						/* decrement the final last appearance								*/
		if (c<lfa)	--lfa;			/* decrement the final first appearance only if it is after bin c	*/
		}	/* end case of gap in clade range */
	}	/* end search for gaps in clade range	*/
free_lvector(rich);
/*return (0);*/
}

/* nu_stratcompat: number of stratigraphic compatibilities for each character
Requires:
	charmat: character matrix;
	compmat: compatibilty matrix;
	ranges: stratigraphic ranges of taxa;
	nstats: number of states for each character
	notu: number of taxa
	strong: whether stratigraphic compatibilty is "weak" (0=based only on first appearances) or "strong" (1=using whole ranges)
	hier: whether we count all (0=divergent+hierarchical) or only hierarchical (1=00->01->11) stratigraphic compatibilities
	UNKNOWNN: code for unknown states
	INAP: code for inapplicable state
****************************************************************************
unsigned long *nu_stratcompat(unsigned long **charmat, unsigned long **compmat, long **ranges, int nchars, int notu, int strong, int hier, int UNKNOWN, int INAP)
{
int c1, c2;


for (c1=0; c1<(nchars-1); ++c1)	{
	for (c2=(c1+1); c2<nchars; ++c2)	{
		}
	}
}	*/


/* stratcompatparambstest: simulations for stratigraphic compatibility with parameters determined by best fit for stratigraphic data
Requires:
	taxonname: group analyzed
	citation: source of data
	summary: summary of observed data
	mbl: simulation parameters
	omatrix: observed matrix
	ctype: character types
	nstates: states per character
	bias: character evolution bias
	maxch: maximum change per character
	depend: character dependencies
	notu: number of taxa
	nchars: number of characters
	compat: compatibility type (general or hierarchical)
	RUNS: simulations
	debug: 1 to seed based on replications, 0 for time
	excl: exclude autapomorphies
	UNKNOWN: code for unknown states
	INAP: code for inapplicable states
Returns:
	simmary[0]: p observed general stratigraphic compatibility
	simmary[1]: p observed strict stratigraphic compatibility
	simmary[2]: p observed hierarchical stratigraphic compatibility
	simmary[3]: expected general stratigraphic compatibility
	simmary[4]: expected strict stratigraphic compatibility
	simmary[5]: expected hierarchical stratigraphic compatibility
	simmary[6]: p proportional hierarchical stratigraphic compatibility
	simmary[7]: expected proportional hierarchical stratigraphic compatibility
	simmary[8]: p hierarchical compatibility consistent with anagenesis
	simmary[9]:	expected hierarchical compatibility consistent with anagenesis
	simmary[10]: p observed “reversal” character being less compatible
	simmary[11]: expected number of stratigraphically incompatible pairs in which “reversed” character is less compatible
	simmary[12]: p observed average compatibility of stratigraphically compatible pairs
	simmary[13]: expected compatibility of stratigraphically compatible pairs
	simmary[14]: p observed average compatibility of stratigraphically incompatible pairs
	simmary[15]: expected compatibility of stratigraphically incompatible pairs
NOTE: Should I add expected compatibility of SIN & GSC characters?  Ummm.... sure, why not?

TO REDO:
	simmary[0]: P[observed general stratigraphic compatibility]
	simmary[1]: E[general stratigraphic compatibility]
	simmary[2]: P[observed strict stratigraphic compatibility]
	simmary[3]: E[strict stratigraphic compatibility]
	simmary[4]: P[observed super strict compatible character pairs]
	simmary[5]: E[super strict compatible character pairs]
	simmary[6]: P[observed hierarchical stratigraphic compatibility]
	simmary[7]: E[hierarchical stratigraphic compatibility]
	simmary[8]: P[observed proportion hierarchical stratigraphic compatibility;]
	simmary[9]: E[proportion hierarchical stratigraphic compatibility;]
	simmary[10]: P[observed hierarchical stratigraphic compatible pairs consistent with anagenesis]
	simmary[11]: E[hierarchical stratigraphic compatible pairs consistent with anagenesis]
	simmary[12]: P[observed hierarchical stratigraphic compatible pairs consistent with budding]
	simmary[13]: E[hierarchical stratigraphic compatible pairs consistent with budding]
	simmary[14]: P[observed compatibilty implies divergent stratigraphically incompatibilty (SIN)]
	simmary[15]: E[compatibilty implies divergent SIN]
	simmary[16]: P[observed compatibilty implies hierarchical stratigraphically incompatibilty SIN]
	simmary[17]: E[stratigraphy implies hierarchical SIN]
	simmary[18]: P[observed stratigraphy implies divergent SIN]
	simmary[19]: E[stratigraphy implies divergent SIN]
	simmary[20]: P[observed stratigraphy implies hierarchical SIN]
	simmary[21]: E[stratigraphy implies hierarchical SIN]
	simmary[22]: P[observed proportion of younger end state pair in SIN characters]
	simmary[23]: E[proportion of younger end state pair in SIN characters]
	simmary[24]: P[compatibility of stratigraphically compatible pairs]
	simmary[25]: E[compatibility of stratigraphically compatible pairs]
	simmary[26]: P[compatibility of SIN pairs]
	simmary[27]: E[compatibility of SIN pairs]
	simmary[28]: P[diff. in compatibility between SC & SIN pairs]
	simmary[29]: E[diff. in compatibility between SC & SIN pairs]
	simmary[30]: E[prop. clade shape ancestral pairs when hierarchical]
	simmary[31]: P[prop. clade shape ancestral pairs when hierarchical]
	simmary[32]: E[prop. clade shape ancestral pairs when divergent]
	simmary[33]: P[prop. clade shape ancestral pairs when divergent]
	simmary[34]: E[overall clade shape]
	simmary[35]: P[overall clade shape]
double	*genstratcomp, *strstratcomp, *hierstratcmp, *hierstratprp;
double	*supstratcomp, *gapfreepairs, *hierstratana, *hierstratbud;
double	*compdivSIN, *comphierSIN, *stratdivSIN, *strathierSIN;
double	*propend2cp, *compSCP, *compSIN;
****************************************************************************/
double *stratcompatparambstest(char taxonname[90], char citation[90], double *summary, double *mbl, long **omatrix, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, char excl, int tiebreaker, int debug, int UNKNOWN, int INAP)
{
/* part one - matrix properties */
int		d, r, s;
int		c, c1, c2;
int		keep;
int 	clades;
int		hierflop=0;
int		*simauts;
//int	uninf, out=0, uninfcompat, posscompat, auts, ttlcompat, comptype=0;
long  	**simatrix, **ranges;
long	**tree, **vtree;
//long	**dbtree;
double	*sstrcmp;
unsigned long **compmat;
unsigned long	*charcomps;
//double	*genstratcomp, *strstratcomp, *hierstratcmp, *hierstratprp,*scha,*sih;
double	*genstratcomp, *strstratcomp, *hierstratcmp, *hierstratprp;
double	*supstratcomp, *hierstratana, *hierstratbud;
double	*compdivSIN, *comphierSIN, *stratdivSIN, *strathierSIN;
double	*propend2cp, *avecompscp, *avecompsin, *avecompdif;
double	*simmary;
double	cp;
double	tiegs=0.0f, tiess=0.0f, tiehs=0.0f, tiehp=0.0f, tiescha=0.0f, tiesih=0.0f, tieavecscp=0.0f, tieavecsin=0.0f, tieavecdif=0.0f;
double	tiexs=0.0f, tieha=0.0f, tiehb=0.0f, tiesincd=0.0f, tiesinch=0.0f, tiesinsd=0.0f, tiesinsh=0.0f, tiesind2=0.0f;
//double	**genrstcmpcor, **strcstcmpcor, **hierstcmpcor;
//long	**histories;
long	secs;
char	outfile[120];
FILE	*output, *debugoutput;
int		end, onset;				/* 2013-04-13: for center of gravity madness	*/
double	*CG, *HAG, *DAG, *HAD, *DAD;
double	*sC;					/* raw clade shape	*/

if (RUNS%2==0)	++RUNS;			/* it's easier to use an odd number....	*/

time(&secs);
srand((unsigned int) secs);

genstratcomp=dvector(RUNS);		/* tallies sstrcomp[5], for simmary[0,1]				*/
strstratcomp=dvector(RUNS);		/* tallies sstrcomp[6], for simmary[2,3]				*/	
supstratcomp=dvector(RUNS);		/* tallies sstrcomp[7], for simmary[8,9]				*/
hierstratcmp=dvector(RUNS);		/* tallies sstrcomp[8], for simmary[4,5]				*/
hierstratprp=dvector(RUNS);		/* tallies sstrcomp[8]/sstrcomp[5], for simmary[6,7]	*/
hierstratana=dvector(RUNS);		/* tallies sstrcomp[9], for simmary[10,11]				*/
hierstratbud=dvector(RUNS);		/* tallies sstrcomp[10], for simmary[12,13]				*/
compdivSIN=dvector(RUNS);		/* tallies sstrcomp[11], for simmary[14,15]				*/
comphierSIN=dvector(RUNS);		/* tallies sstrcomp[12], for simmary[16,17]				*/
stratdivSIN=dvector(RUNS);		/* tallies sstrcomp[13], for simmary[18,19]				*/
strathierSIN=dvector(RUNS);		/* tallies sstrcomp[14], for simmary[20,21]				*/
propend2cp=dvector(RUNS);		/* tallies sstrcomp[15], for simmary[22,23]				*/
avecompscp=dvector(RUNS);		/* tallies for simmary[24,25]							*/
avecompsin=dvector(RUNS);		/* tallies for simmary[26.27]							*/
avecompdif=dvector(RUNS);		/* tallies for simmary[28,29]							*/
HAG=dvector(RUNS);				/* tallies for simmary[30]								*/
HAD=dvector(RUNS);				/* tallies for simmary[31]								*/
DAG=dvector(RUNS);				/* tallies for simmary[32]								*/
DAD=dvector(RUNS);				/* tallies for simmary[33]								*/
CG=dvector(RUNS);				/* tallies for simmary[34]								*/

ranges=lmatrix(notu,2);
simatrix=lmatrix(notu,nchars);
simmary=dvector(35);

d=0;							/* for debugging	*/
for (r=d; r<RUNS; ++r)	{
	/* 2011-03-28: BLOWOUT at r=5	*/
	genstratcomp[r]=strstratcomp[r]=hierstratcmp[r]=hierstratprp[r]=0.0f;
	supstratcomp[r]=hierstratana[r]=hierstratbud[r]=compdivSIN[r]=comphierSIN[r]=stratdivSIN[r]=strathierSIN[r]=propend2cp[r]=avecompscp[r]=avecompsin[r]=0.0f;

	for (s=0; s<notu; ++s)	{
		for (c=0; c<nchars; ++c)	{
			if (omatrix[s][c]==UNKNOWN)
				simatrix[s][c]=UNKNOWN;
			else if (omatrix[s][c]==INAP)
				simatrix[s][c]=INAP;
			else
				simatrix[s][c]=0;
			}
		}
	if (debug==1)	{
		srand((unsigned int) (r+1)*notu*nchars);
		}
	
	/* 	trees[notu-1]: branch length of otus	
		trees[notu]:   branch length of clades
		trees[notu+1]: first appearances of otus
		trees[notu+2]: last appearances of otus
	**********************************************/
	tree=evolvetree(notu,mbl,1);
	clades=cladecountbytaxa(tree,notu);
	/* pull range data out of back of tree matrix	*/
	for (s=0; s<notu; ++s)	{
		ranges[s][0]=tree[notu+1][s];
		ranges[s][1]=tree[notu+2][s];
		}
	cleancladerangedata(ranges,notu);
	
	/* set of center of gravity madness	*/
	onset=minlmatrixcol(ranges,notu,0);
	end=minlmatrixcol(ranges,notu,0);
	CG[r]=centerofgravityfromranges(ranges,notu,0,onset,end);
	sC=cladeshapefromranges(ranges,notu,0);

	vtree=VennTreePlus(tree,clades,notu,notu);
	/* evolve character matrix, simatrix	*/
//	long **evolvetocompat(long **tree, int tcomp, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int *depend, int comptype, int UNKNOWN, int INAP)
//	simatrix=evolvetocompat(tree,empcompat,notu,matrix,nchars,nstates,ctype,bias,chmax,depend,comptype,UNKNOWN,INAP);
	simatrix=evolvetocompatibility(vtree, compat, notu, clades, simatrix, nchars, nstates, ctype, bias, maxch, depend, 0, UNKNOWN, INAP);
	
	if (d>0)	{
		debugoutput=fopen("SuspectMatrix.txt","w");
		for (s=0; s<notu; ++s)	{
			for (c=0; c<nchars; ++c)	{
				if (simatrix[s][c]==UNKNOWN)	fprintf(debugoutput,"?\t");
				else if (simatrix[s][c]==INAP)	fprintf(debugoutput,"—\t");
				else							fprintf(debugoutput,"%d\t",simatrix[s][c]);
				}
			fprintf(debugoutput,"%d\t%d\n",ranges[s][0],ranges[s][1]);
			}
		fclose(debugoutput);
		}
//	simatrix=evolvematrix(vtree, notu, clades, simatrix, nchars, nstates, ctype, bias, maxch, 2*nchars, UNKNOWN, INAP);
/*		sstates=numberstates(simatrix,notu,nchars,UNKNOWN,INAP);	*/
//	compat=nu_comp(nstates, notu, simatrix, ctype, nchars, 0, 0, UNKNOWN, INAP);
	compmat=compatible(nstates,notu,simatrix,ctype,nchars,0,0,UNKNOWN,INAP);
	charcomps=char_comp(nstates,notu,simatrix,ctype,nchars,0,0,UNKNOWN,INAP);
	simauts=autapomorphies(simatrix,nstates,notu,nchars,UNKNOWN,INAP);	/* get the minimum number of derived states	*/

/*	for (c1=0; c1<nchars; ++c1)	++examples[compmat[c1][c1]];	*/
	cp=0.0f;
	for (c1=0; c1<(nchars-1); ++c1)	{
		/* reinstated 2012-02-29	*/
		if (excl=='y')		while ((simauts[c1]<2 || simauts[c1]>(notu-2)) && c1<(nchars-1))	++c1;
		if (c1>=nchars)		break;
		for (c2=(c1+1); c2<nchars; ++c2)	{
			/* reinstated 2012-02-29	*/
			if (excl=='y')		while ((simauts[c2]<2 || simauts[c2]>(notu-2)) && c2<(nchars-1))	++c2;
			if (c2>=nchars)		break;

			if (compmat[c1][c2]==1)	{
				/*sstrcmp=stratcompatfull(ranges, simatrix, nstates, c1, c2, notu, UNKNOWN, INAP);	*/
				/* 2011-10-14: changed to "plus"													*/
				sstrcmp=stratcompatfullplusplusplus(ranges, simatrix, nstates, charcomps, c1, c2, notu, onset, end, tiebreaker, UNKNOWN, INAP);
				/* rewrite this so that we simply take the median values from sstrcomp (simulated strat comp) for expectations	*/
				/* and the distribution for p-values when examining the real results											*/
				cp+=sstrcmp[4];										/* compatible comparisons						*/
				genstratcomp[r]+=sstrcmp[5];						/* general stratigraphic consistency			*/
				strstratcomp[r]+=sstrcmp[6];						/* strict stratigraphic consistency				*/
				supstratcomp[r]+=sstrcmp[7];						/* superstrict stratigraphic consistency		*/
				hierstratcmp[r]+=sstrcmp[8];						/* hierarchical stratigraphic consistency		*/
				hierstratana[r]+=sstrcmp[9];						/* HSC consistent with anagenesis				*/
				hierstratbud[r]+=sstrcmp[10];						/* HSC consistent with budding					*/
				compdivSIN[r]+=sstrcmp[11];							/* compatibility-supported divergent SIN		*/
				comphierSIN[r]+=sstrcmp[12];						/* compatibility-supported hierarchical SIN		*/
				stratdivSIN[r]+=sstrcmp[13];						/* stratigraphy-supported divergent SIN			*/
				strathierSIN[r]+=sstrcmp[14];						/* stratigraphy-supported hierarchical SIN		*/
				propend2cp[r]+=sstrcmp[15];							/* relative diversity of younger SIN end pair	*/
				if (sstrcmp[5]>0)				avecompscp[r]+=(((double) charcomps[c1])+((double) charcomps[c2]))/2;
				if (sstrcmp[5]<sstrcmp[4])		avecompsin[r]+=(((double) charcomps[c1])+((double) charcomps[c2]))/2;


/*				scha[r]+=sstrcmp[9];						/* hierarchical compatible pairs consistent with anagenesis			2011-10-14	*/
/*				sih[r]+=sstrcmp[11];						/* stratigraphically incompatible pairs consistent with hierarchy	2011-10-14	*/
				/* 2011-10-14: tally compatibilities of stratigraphically compatible & incompatible character pairs	*/
/*				if (sstrcmp[5]==1)	{
					avecompscp[r]+=((double) (charcomps[c1]+charcomps[c2]))/2;
					scgw+=1.0;													/* this will tally characters, not state pairs						*/
/*					}
				else	{
					avecompsin[r]+=((double) (charcomps[c1]+charcomps[c2]))/2;
					sci+=1.0;
					}
	
				/* this will give the proportion of compatibilities that are stratigraphically compatible by compatibility	*/
/*				genrstcmpcor[compmat[c1][c1]]+=sstrcmp[5];			/* general stratigraphic consistency		*/
/*				strcstcmpcor[compmat[c1][c1]]+=sstrcmp[6];			/* strict stratigraphic consistency			*/
/*				hierstcmpcor[compmat[c1][c1]]+=sstrcmp[8];			/* hierarchical stratigraphic consistency	*/
/*				genrstcmpcor[compmat[c2][c2]]+=sstrcmp[5];			/* general stratigraphic consistency		*/
/*				strcstcmpcor[compmat[c2][c2]]+=sstrcmp[6];			/* strict stratigraphic consistency			*/
/*				hierstcmpcor[compmat[c2][c2]]+=sstrcmp[8];			/* hierarchical stratigraphic consistency	*/
				free_dvector(sstrcmp);
				}
			}
		}
	
	/* get proportions of hierarchical stratigraphic characters relative to total								*/
	/* Also, get proportions of particular types of hierarchical stratigraphic characters						*/
	hierstratprp[r]+=(hierstratcmp[r]/genstratcomp[r]);						/* proportional hierarchical stratigraphic compatibility	*/
	/* if there is no hierarchical compatibility, then modify the routine: 2011-11-01	*/
	if (hierstratcmp[r]>0)	{
		hierstratana[r]/=hierstratcmp[r];										/* proportion of HSC matching anagenesis					*/
		hierstratbud[r]/=hierstratcmp[r];										/* proportion of HSC demanding budding						*/
		}
	else	{
		++hierflop;																/* tally runs in which there was no hierarchical compatibility	*/
		hierstratana[r]=hierstratbud[r]=MAXRAND;
		}
	/* get average compatibilities of stratigraphically compatible and incompatible character pairs					*/
	avecompscp[r]/=genstratcomp[r];											/* added 2011-10-14	*/
	avecompsin[r]/=(cp-genstratcomp[r]);									/* added 2011-10-14	*/
	avecompdif[r]=avecompscp[r]-avecompsin[r];								/* added 2011-11-01: because I sort these later, this is otherwise lost */
	/* get proportion of SIN character pairs matching divergent or hierarchical patterns				*/
	compdivSIN[r]/=(cp-genstratcomp[r]);									/* compatibility-supported divergent SIN		*/
	comphierSIN[r]/=(cp-genstratcomp[r]);									/* compatibility-supported hierarchical SIN		*/
	stratdivSIN[r]/=(cp-genstratcomp[r]);									/* stratigraphy-supported divergent SIN			*/
	strathierSIN[r]/=(cp-genstratcomp[r]);									/* stratigraphy-supported hierarchical SIN		*/
	propend2cp[r]/=(cp-genstratcomp[r]);									/* relative diversity of younger SIN end pair	*/
	/* cp will vary a little from simulation to simulation because it is based on state numbers, so report things in proportions	*/
	genstratcomp[r]/=((double) cp);											/* general stratigraphic consistency						*/
	strstratcomp[r]/=((double) cp);											/* strict stratigraphic consistency							*/
	supstratcomp[r]/=((double) cp);											/* super strict stratigraphic consistency					*/
	hierstratcmp[r]/=((double) cp);											/* hierarchical stratigraphic consistency					*/

	/* tally cases where simulated values are lower than observed values		*/
	if (genstratcomp[r]<summary[0])		simmary[0]+=(1/((double) RUNS));
	if (strstratcomp[r]<summary[1])		simmary[2]+=(1/((double) RUNS));
	if (supstratcomp[r]<summary[2])		simmary[4]+=(1/((double) RUNS));
	if (hierstratcmp[r]<summary[3])		simmary[6]+=(1/((double) RUNS));
	if (hierstratprp[r]<summary[4])		simmary[8]+=(1/((double) RUNS));
	if (hierstratana[r]<summary[5])		simmary[10]+=(1/((double) RUNS));
	if (hierstratbud[r]<summary[6])		simmary[12]+=(1/((double) RUNS));
	if (compdivSIN[r]<summary[7])		simmary[14]+=(1/((double) RUNS));
	if (comphierSIN[r]<summary[8])		simmary[16]+=(1/((double) RUNS));
	if (stratdivSIN[r]<summary[9])		simmary[18]+=(1/((double) RUNS));
	if (strathierSIN[r]<summary[10])	simmary[20]+=(1/((double) RUNS));
	if (propend2cp[r]<summary[11])		simmary[22]+=(1/((double) RUNS));
	if (avecompscp[r]<summary[12])		simmary[24]+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompsin[r]<summary[13])		simmary[26]+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompdif[r]<(summary[12]-summary[13]))
										simmary[27]+=(1/((double) RUNS));					/* added 2011-11-01: because I sort these later, this is otherwise lost */

	/* tally cases where simulated values are exactly equal to observed values	*/
	if (genstratcomp[r]==summary[0])	tiegs+=(1/((double) RUNS));
	if (strstratcomp[r]==summary[1])	tiess+=(1/((double) RUNS));
	if (supstratcomp[r]==summary[2])	tiexs+=(1/((double) RUNS));
	if (hierstratcmp[r]==summary[3])	tiehs+=(1/((double) RUNS));
	if (hierstratprp[r]==summary[4])	tiehp+=(1/((double) RUNS));
	if (hierstratana[r]==summary[5])	tieha+=(1/((double) RUNS));
	if (hierstratbud[r]==summary[6])	tiehb+=(1/((double) RUNS));
	if (compdivSIN[r]==summary[7])		tiesincd+=(1/((double) RUNS));
	if (comphierSIN[r]==summary[8])		tiesinch+=(1/((double) RUNS));
	if (stratdivSIN[r]==summary[9])		tiesinsd+=(1/((double) RUNS));
	if (strathierSIN[r]==summary[10])	tiesinsh+=(1/((double) RUNS));
	if (propend2cp[r]==summary[11])		tiesind2+=(1/((double) RUNS));
	if (avecompscp[r]==summary[12])		tieavecscp+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompsin[r]==summary[13])		tieavecsin+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompdif[r]==(summary[12]-summary[13]))
										tieavecdif+=(1/((double) RUNS));					/* added 2011-11-01: because I sort these later, this is otherwise lost */

	if ((r%10)==9)	{
//			printf("Doing Rate = %5.4f",mbl[2]);
		if (r>10)	{
			printf("\b\b\b\b\b");					/* clear ", Tree %d"		*/
			}
		if (r>18 && r<100) 			printf("\b\b\b");
		else if (r>99 && r<1000)	printf("\b\b\b\b");
		else if (r>999)				printf("\b\b\b\b\b");
		printf("Tree %d\n", r+1);
		}
	
	free_ulmatrix(compmat,nchars,nchars);
	free_lmatrix(vtree,clades+1,notu);
	free_lmatrix(tree,clades+3,notu);
	}

/* sort for medians	*/
genstratcomp=dshellsort_inc(genstratcomp,RUNS); 	/* 0:  sort stratigraphic compatiblity						*/
strstratcomp=dshellsort_inc(strstratcomp,RUNS); 	/* 1:  sort strict stratigraphic compatiblity				*/
supstratcomp=dshellsort_inc(supstratcomp,RUNS);		/* 2:  sort super strict stratigraphic compatibility		*/
hierstratcmp=dshellsort_inc(hierstratcmp,RUNS);		/* 3:  sort hierarchical stratigraphic compatiblity (HSC)	*/
hierstratprp=dshellsort_inc(hierstratprp,RUNS);		/* 4:  sort proportion HSC									*/
hierstratana=dshellsort_inc(hierstratana,RUNS);		/* 5:  sort HSC consistent with anagenesis					*/
hierstratbud=dshellsort_inc(hierstratbud,RUNS);		/* 6:  sort HSC demanding budding							*/
compdivSIN=dshellsort_inc(compdivSIN,RUNS);			/* 7:  sort compatibility suggests divergent SIN			*/
comphierSIN=dshellsort_inc(comphierSIN,RUNS);		/* 8:  sort compatibility suggests hierarchical SIN			*/
stratdivSIN=dshellsort_inc(stratdivSIN,RUNS);		/* 9:  sort stratigraphy suggests divergent SIN				*/
strathierSIN=dshellsort_inc(strathierSIN,RUNS);		/* 10: sort stratigraphy suggests hierarchical SIN			*/
propend2cp=dshellsort_inc(propend2cp,RUNS);			/* 11: sort relative diversity of derived end pair			*/
avecompscp=dshellsort_inc(avecompscp,RUNS);			/* 12: sort average compatibilities of SC character pairs	*/
avecompsin=dshellsort_inc(avecompsin,RUNS);			/* 13: sort average compatibilities of SIN character pairs	*/
avecompdif=dshellsort_inc(avecompdif,RUNS);			/* 14: sort average differences between SC & SIN pairs		*/

/* get alpha values: if observed values are "low" then we want p[observed or lower]: so, add ties to total	*/
/* for "high" values, we do not need to do this: 1-simmary gives p[observed or more extreme] as we tallied only sims less than observed	*/
/* we want the probability of the observed OR MORE EXTREME; because we tallied every time the observed was lower than what we saw,
		we don't need to do this when p> 0.5	*/
if (simmary[0]<0.5)		simmary[0]+=tiegs;			/* general stratigraphic compatibility	*/
if (simmary[2]<0.5)		simmary[2]+=tiess;			/* strict stratigraphic compatibility	*/
if (simmary[4]<0.5)		simmary[4]+=tiexs;
if (simmary[6]<0.5)		simmary[6]+=tiehs;
if (simmary[8]<0.5)		simmary[8]+=tiehp;
if (simmary[10]<0.5)	simmary[10]+=tieha;
if (simmary[12]<0.5)	simmary[12]+=tiehb;
if (simmary[14]<0.5)	simmary[14]+=tiesincd;
if (simmary[16]<0.5)	simmary[16]+=tiesinch;
if (simmary[18]<0.5)	simmary[18]+=tiesinsd;
if (simmary[20]<0.5)	simmary[20]+=tiesinsh;
if (simmary[22]<0.5)	simmary[22]+=tiesind2;
if (simmary[24]<0.5)	simmary[24]+=tieavecscp;
if (simmary[26]<0.5)	simmary[26]+=tieavecsin;
if (simmary[28]<0.5)	simmary[28]+=tieavecsin;

/* adjust for runs with no hierarchical compatibility	*/
simmary[10]*=(((double) (RUNS-hierflop))/((double) RUNS));		/* adjust significance of anagenetic hierarchical compatibility	*/
simmary[12]*=(((double) (RUNS-hierflop))/((double) RUNS));		/* adjust significance of budding hierarchical compatibility	*/


/* tally median expected stratigraphic compatibilities	*/
simmary[1]=genstratcomp[RUNS/2];
simmary[3]=strstratcomp[RUNS/2];
simmary[5]=supstratcomp[RUNS/2];
simmary[7]=hierstratcmp[RUNS/2];
simmary[9]=hierstratprp[RUNS/2];
simmary[11]=hierstratana[(RUNS-hierflop)/2];			/* accommodate runs with no hierarchical compatibility	*/
simmary[13]=hierstratbud[(RUNS-hierflop)/2];			/* accommodate runs with no hierarchical compatibility	*/
simmary[15]=compdivSIN[RUNS/2];
simmary[17]=comphierSIN[RUNS/2];
simmary[19]=stratdivSIN[RUNS/2];
simmary[21]=strathierSIN[RUNS/2];
simmary[23]=propend2cp[RUNS/2];
simmary[25]=avecompscp[RUNS/2];
simmary[27]=avecompsin[RUNS/2];
simmary[29]=avecompdif[RUNS/2];

/* return output to here!	*/
/* TEMPORARY!  UNTIL I'VE DEBUGGED THIS.....	*/
keep=1;
//keep=0;
//printf("Enter '1' if you want output for this taxon: ");
//scanf("%i",&keep);	

if (keep==1)	{
	strcpy(outfile,taxonname);
	strcat(outfile,"_");
	strcat(outfile,citation);
	strcat(outfile,"_Expectations");
	if (mbl[4]==0)	strcat(outfile,"_Bud.xls");
	else			strcat(outfile,"_Bif.xls");
	
	output=fopen(outfile,"w");
	fprintf(output,"Run");
	fprintf(output,"\tGeneral Stratigraphic Compatibility (SC)");
	fprintf(output,"\tStrict SC");
	fprintf(output,"\tSuper Strict SC");
	fprintf(output,"\tHierarchical SC");
	fprintf(output,"\tProportional HSC");
	fprintf(output,"\tAnagenetic HSC");
	fprintf(output,"\tBudding HSC");
	fprintf(output,"\tDivergent SIN (Comp)");
	fprintf(output,"\tHierarchical SIN (Comp)");
	fprintf(output,"\tDivergent SIN (Strat)");
	fprintf(output,"\tHierarchical SIN (Strat)");
	fprintf(output,"\tProp Diversity Younger SIN Pair");
	fprintf(output,"\tµ Character Compabibility (SC)");					/* added 2011-10-14	*/
	fprintf(output,"\tµ Character Compabibility (SIN)");					/* added 2011-10-14	*/
	fprintf(output,"\tµ Compabibility Diff (SC-SIN)\n");				/* added 2011-11-01	*/
	
	for (r=0; r<RUNS; ++r)	{
		fprintf(output,"%d",r+1);
		fprintf(output,"\t%4.3f",genstratcomp[r]);
		fprintf(output,"\t%4.3f",strstratcomp[r]);
		fprintf(output,"\t%4.3f",supstratcomp[r]);
		fprintf(output,"\t%4.3f",hierstratcmp[r]);
		fprintf(output,"\t%4.3f",hierstratprp[r]);
		if (hierstratana[r]<MAXRAND)	fprintf(output,"\t%4.3f",hierstratana[r]);
		else							fprintf(output,"\t•");
		if (hierstratbud[r]<MAXRAND)	fprintf(output,"\t%4.3f",hierstratbud[r]);
		else							fprintf(output,"\t•");
		fprintf(output,"\t%4.3f",compdivSIN[r]);
		fprintf(output,"\t%4.3f",comphierSIN[r]);
		fprintf(output,"\t%4.3f",stratdivSIN[r]);
		fprintf(output,"\t%4.3f",strathierSIN[r]);
		fprintf(output,"\t%4.3f",propend2cp[r]);
		fprintf(output,"\t%4.3f",avecompscp[r]);					/* added 2011-10-14	*/
		fprintf(output,"\t%4.3f",avecompsin[r]);					/* added 2011-10-14	*/
		fprintf(output,"\t%4.3f\n",avecompdif[r]);					/* added 2011-11-01	*/
		}
	fclose(output);
	}


free_dvector(genstratcomp);
free_dvector(strstratcomp);
free_dvector(supstratcomp);		/* tallies sstrcomp[7], for simmary[8,9]				*/
free_dvector(hierstratcmp);
free_dvector(hierstratprp);
free_dvector(hierstratana);		/* tallies sstrcomp[9], for simmary[10,11]				*/
free_dvector(hierstratbud);		/* tallies sstrcomp[10], for simmary[12,13]				*/
free_dvector(compdivSIN);		/* tallies sstrcomp[11], for simmary[14,15]				*/
free_dvector(comphierSIN);		/* tallies sstrcomp[12], for simmary[16,17]				*/
free_dvector(stratdivSIN);		/* tallies sstrcomp[13], for simmary[18,19]				*/
free_dvector(strathierSIN);		/* tallies sstrcomp[14], for simmary[20,21]				*/
free_dvector(propend2cp);		/* tallies sstrcomp[15], for simmary[22,23]				*/
free_dvector(avecompscp);
free_dvector(avecompsin);
free_dvector(avecompdif);
free_lmatrix(ranges,notu,2);

return simmary;
}

/* stratcompatparambstest2: simulations for stratigraphic compatibility with parameters determined by best fit for stratigraphic data
Requires:
	taxonname: group analyzed
	citation: source of data
	summary: summary of observed data
/*		summary[0]:  prop. of generally compatible pairs										2011-10-15	*/
/*		summary[1]:  prop. of strictly compatible pairs											2011-10-15	*/
/*		summary[2]:  prop. of super strictly compatible pairs									2011-10-15	*/
/*		summary[3]:  prop. of hierarchically compatible pairs									2011-10-15	*/
/*		summary[4]:  prop. of hierarchically compatible pairs given strat. consistency			2011-10-15	*/
/*		summary[5]:  prop. of hierarchically compatible pairs consistent with anagenesis		2011-10-15	*/
/*		summary[6]:  prop. of hierarchically compatible pairs consistent with budding			2011-10-15	*/
/*		summary[7]:  prop. of SIN pairs consistent with divergence given compatibility			2011-10-15	*/
/*		summary[8]:  prop. of SIN pairs consistent with hierarchy given compatibility			2011-10-15	*/
/*		summary[9]:  prop. of SIN pairs consistent with divergence given stratigraphy			2011-10-15	*/
/*		summary[10]: prop. of SIN pairs consistent with hierarchy given stratigraphy			2011-10-15	*/
/*		summary[11]: average relative diversity of derived end-pair for SIN pairs				2011-10-15	*/
/*		summary[12]: average compatibility of stratigraphically compatible pairs				2011-10-15	*/
/*		summary[13]: average compatibility of stratigraphically incompatible pairs				2011-10-15	*/
/*		summary[14]: average center of gravity for hierarchical 00								2013-04-14	*/
/*		summary[15]: average center of gravity for divergent 00									2013-04-14	*/
/*		summary[16]: Center of gravity for whole clade											2013-04-14	*/
/*		summary[17]: average duration of hierarchial 00's										2013-04-14	*/
/*		summary[18]: average duration of divergent 00's											2013-04-14	*/
/*		summary[19]: "living fossil" hierarchial 00's											2013-04-15	*/
/*		summary[20]: "living fossil" divergent 00's												2013-04-15	*/
/*	mbl: simulation parameters
	omatrix: observed matrix
	ctype: character types
	nstates: states per character
	bias: character evolution bias
	maxch: maximum change per character
	depend: character dependencies
	notu: number of taxa
	nchars: number of characters
	compat: compatibility type (general or hierarchical)
	RUNS: simulations
	debug: 1 to seed based on replications, 0 for time
	excl: exclude autapomorphies
	UNKNOWN: code for unknown states
	INAP: code for inapplicable states
Returns:
	simmary[0]: P[observed general stratigraphic compatibility]
	simmary[1]: E[general stratigraphic compatibility]
	simmary[2]: P[observed strict stratigraphic compatibility]
	simmary[3]: E[strict stratigraphic compatibility]
	simmary[4]: P[observed super strict compatible character pairs]
	simmary[5]: E[super strict compatible character pairs]
	simmary[6]: P[observed hierarchical stratigraphic compatibility]
	simmary[7]: E[hierarchical stratigraphic compatibility]
	simmary[8]: P[observed proportion hierarchical stratigraphic compatibility;]
	simmary[9]: E[proportion hierarchical stratigraphic compatibility;]
	simmary[10]: P[observed hierarchical stratigraphic compatible pairs consistent with anagenesis]
	simmary[11]: E[hierarchical stratigraphic compatible pairs consistent with anagenesis]
	simmary[12]: P[observed hierarchical stratigraphic compatible pairs consistent with budding]
	simmary[13]: E[hierarchical stratigraphic compatible pairs consistent with budding]
	simmary[14]: P[observed compatibilty implies divergent stratigraphically incompatibilty (SIN)]
	simmary[15]: E[compatibilty implies divergent SIN]
	simmary[16]: P[observed compatibilty implies hierarchical stratigraphically incompatibilty SIN]
	simmary[17]: E[stratigraphy implies hierarchical SIN]
	simmary[18]: P[observed stratigraphy implies divergent SIN]
	simmary[19]: E[stratigraphy implies divergent SIN]
	simmary[20]: P[observed stratigraphy implies hierarchical SIN]
	simmary[21]: E[stratigraphy implies hierarchical SIN]
	simmary[22]: P[observed proportion of younger end state pair in SIN characters]
	simmary[23]: E[proportion of younger end state pair in SIN characters]
	simmary[24]: P[compatibility of stratigraphically compatible pairs]
	simmary[25]: E[compatibility of stratigraphically compatible pairs]
	simmary[26]: P[compatibility of SIN pairs]
	simmary[27]: E[compatibility of SIN pairs]
	simmary[28]: P[diff. in compatibility between SC & SIN pairs]
	simmary[29]: E[diff. in compatibility between SC & SIN pairs]
	simmary[30]: P[clade shape 00 when hierarchical]
	simmary[31]: E[clade shape 00 when hierarchical]
	simmary[32]: P[clade shape 00 when divergent]
	simmary[33]: E[clade shape 00 when divergent]
	simmary[34]: P[overall clade shape]
	simmary[35]: E[overall clade shape]
	simmary[36]: P[hierarchical 00 duration]
	simmary[37]: E[hierarchical 00 duration]
	simmary[38]: P[divergent 00 duration]
	simmary[39]: E[divergent 00 duration]
	simmary[40]: P[hierarchical "living fossils"]
	simmary[41]: E[hierarchical "living fossils"]
	simmary[42]: P[divergent "living fossils"]
	simmary[43]: E[divergent "living fossils"]
****************************************************************************/
double *stratcompatparambstest2(char taxonname[90], char citation[90], double *summary, double *mbl, long **omatrix, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, char excl, int tiebreaker, int debug, int UNKNOWN, int INAP)
{
/* part one - matrix properties */
int		d, r, s;
int		c, c1, c2;
int		keep;
int 	clades;
int		hierflop=0;
int		*simauts;
long  	**simatrix, **ranges;
long	**tree, **vtree;
double	*sstrcmp;
unsigned long **compmat;
unsigned long	*charcomps;
double	*genstratcomp, *strstratcomp, *hierstratcmp, *hierstratprp;
double	*supstratcomp, *hierstratana, *hierstratbud;
double	*compdivSIN, *comphierSIN, *stratdivSIN, *strathierSIN;
double	*propend2cp, *avecompscp, *avecompsin, *avecompdif;
double	*simmary;
double	cp;
double	tiegs=0.0f, tiess=0.0f, tiehs=0.0f, tiehp=0.0f, tiescha=0.0f, tiesih=0.0f, tieavecscp=0.0f, tieavecsin=0.0f, tieavecdif=0.0f;
double	tiexs=0.0f, tieha=0.0f, tiehb=0.0f, tiesincd=0.0f, tiesinch=0.0f, tiesinsd=0.0f, tiesinsh=0.0f, tiesind2=0.0f;
long	secs;
char	outfile[120];
FILE	*output, *debugoutput;
/* 2013-04-13: for center of gravity madness	*/
int		term, onset;				
double	*CG, *aveCGH, *aveCGD, *aveADH, *aveADD, *freqLFH, *freqLFD;
double	tieavecgh=0.0f, tieavecgd=0.0f, tieaveadh=0.0f, tieaveadd=0.0f, tierc1c0=0.0f, tiefreqlfh=0.0f, tiefreqlfd=0.0f;
int	hiercounter, divcounter;	/* cases of hierarchy or divergence clearly	*/

if (RUNS%2==0)	++RUNS;			/* it's easier to use an odd number....	*/

time(&secs);
srand((unsigned int) secs);

genstratcomp=dvector(RUNS);		/* tallies sstrcomp[5], for simmary[0,1]				*/
strstratcomp=dvector(RUNS);		/* tallies sstrcomp[6], for simmary[2,3]				*/	
supstratcomp=dvector(RUNS);		/* tallies sstrcomp[7], for simmary[8,9]				*/
hierstratcmp=dvector(RUNS);		/* tallies sstrcomp[8], for simmary[4,5]				*/
hierstratprp=dvector(RUNS);		/* tallies sstrcomp[8]/sstrcomp[5], for simmary[6,7]	*/
hierstratana=dvector(RUNS);		/* tallies sstrcomp[9], for simmary[10,11]				*/
hierstratbud=dvector(RUNS);		/* tallies sstrcomp[10], for simmary[12,13]				*/
compdivSIN=dvector(RUNS);		/* tallies sstrcomp[11], for simmary[14,15]				*/
comphierSIN=dvector(RUNS);		/* tallies sstrcomp[12], for simmary[16,17]				*/
stratdivSIN=dvector(RUNS);		/* tallies sstrcomp[13], for simmary[18,19]				*/
strathierSIN=dvector(RUNS);		/* tallies sstrcomp[14], for simmary[20,21]				*/
propend2cp=dvector(RUNS);		/* tallies sstrcomp[15], for simmary[22,23]				*/
avecompscp=dvector(RUNS);		/* tallies for simmary[24,25]							*/
avecompsin=dvector(RUNS);		/* tallies for simmary[26.27]							*/
avecompdif=dvector(RUNS);		/* tallies for simmary[28,29]							*/
aveCGH=dvector(RUNS);			/* tallies for simmary[30,31]							*/
aveCGD=dvector(RUNS);			/* tallies for simmary[32,33]							*/
aveADH=dvector(RUNS);			/* tallies for simmary[30,31]							*/
aveADD=dvector(RUNS);			/* tallies for simmary[32,33]							*/
freqLFH=dvector(RUNS);			/* tallies for simmary[34,35]							*/
freqLFD=dvector(RUNS);			/* tallies for simmary[36,37]							*/
CG=dvector(RUNS);				/* tallies for simmary[38,39]							*/

ranges=lmatrix(notu,2);			/* simulated ranges										*/
simatrix=lmatrix(notu,nchars);	/* simulated morphologies								*/
simmary=dvector(45);			/* OUTPUT												*/

d=0;							/* for debugging	*/
for (r=d; r<RUNS; ++r)	{
	/* 2011-03-28: BLOWOUT at r=5	*/
	genstratcomp[r]=strstratcomp[r]=hierstratcmp[r]=hierstratprp[r]=0.0f;
	supstratcomp[r]=hierstratana[r]=hierstratbud[r]=compdivSIN[r]=comphierSIN[r]=stratdivSIN[r]=strathierSIN[r]=propend2cp[r]=avecompscp[r]=avecompsin[r]=0.0f;
	CG[r]=aveCGH[r]=aveCGD[r]=aveADH[r]=aveADD[r]=0.0f;
	
	for (s=0; s<notu; ++s)	{
		for (c=0; c<nchars; ++c)	{
			if (omatrix[s][c]==UNKNOWN)
				simatrix[s][c]=UNKNOWN;
			else if (omatrix[s][c]==INAP)
				simatrix[s][c]=INAP;
			else
				simatrix[s][c]=0;
			}
		}
	if (debug==1)	{
		srand((unsigned int) (r+1)*notu*nchars);
		}
	
	/* 	trees[notu-1]: branch length of otus	
		trees[notu]:   branch length of clades
		trees[notu+1]: first appearances of otus
		trees[notu+2]: last appearances of otus
	**********************************************/
	tree=evolvetree(notu,mbl,1);
	clades=cladecountbytaxa(tree,notu);
	/* pull range data out of back of tree matrix	*/
	for (s=0; s<notu; ++s)	{
		ranges[s][0]=tree[notu+1][s];
		ranges[s][1]=tree[notu+2][s];
		}
	cleancladerangedata(ranges,notu);
	
	/* set of center of gravity madness	*/
	onset=minlmatrixcol(ranges,notu,0);							/* first stage in which something appears						*/
	term=maxlmatrixcol(ranges,notu,0);							/* last stage in which something appears						*/
//	sC=cladeshapefromranges(ranges,notu,0);						/* vector with c0, c1, c2 & c3 for simulated clade				*/
//	CG[r]=sC[1]/sC[0];											/* scalar = c1/c0 for simulated clade (corrected 2013-04-16		*/
	CG[r]=centerofgravityfromranges(ranges,notu,0,onset,term);	/* scalar = c1/c0 for simulated clade (corrected 2013-04-16		*/
	if (CG[r]==summary[16])	++tierc1c0;							/* observed center of gravity is matched						*/

	vtree=VennTreePlus(tree,clades,notu,notu);
	/* evolve character matrix, simatrix	*/
//	long **evolvetocompat(long **tree, int tcomp, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int *depend, int comptype, int UNKNOWN, int INAP)
//	simatrix=evolvetocompat(tree,empcompat,notu,matrix,nchars,nstates,ctype,bias,chmax,depend,comptype,UNKNOWN,INAP);
	simatrix=evolvetocompatibility(vtree, compat, notu, clades, simatrix, nchars, nstates, ctype, bias, maxch, depend, 0, UNKNOWN, INAP);
	
	if (d>0)	{
		debugoutput=fopen("SuspectMatrix.txt","w");
		for (s=0; s<notu; ++s)	{
			for (c=0; c<nchars; ++c)	{
				if (simatrix[s][c]==UNKNOWN)	fprintf(debugoutput,"?\t");
				else if (simatrix[s][c]==INAP)	fprintf(debugoutput,"—\t");
				else							fprintf(debugoutput,"%d\t",simatrix[s][c]);
				}
			fprintf(debugoutput,"%d\t%d\n",ranges[s][0],ranges[s][1]);
			}
		fclose(debugoutput);
		}
//	simatrix=evolvematrix(vtree, notu, clades, simatrix, nchars, nstates, ctype, bias, maxch, 2*nchars, UNKNOWN, INAP);
/*		sstates=numberstates(simatrix,notu,nchars,UNKNOWN,INAP);	*/
//	compat=nu_comp(nstates, notu, simatrix, ctype, nchars, 0, 0, UNKNOWN, INAP);
	compmat=compatible(nstates,notu,simatrix,ctype,nchars,0,0,UNKNOWN,INAP);
	charcomps=char_comp(nstates,notu,simatrix,ctype,nchars,0,0,UNKNOWN,INAP);
	simauts=autapomorphies(simatrix,nstates,notu,nchars,UNKNOWN,INAP);	/* get the minimum number of derived states	*/

/*	for (c1=0; c1<nchars; ++c1)	++examples[compmat[c1][c1]];	*/
	cp=0.0f;
	hiercounter=divcounter=0;
	for (c1=0; c1<(nchars-1); ++c1)	{
		/* reinstated 2012-02-29	*/
		if (excl=='y')		while ((simauts[c1]<2 || simauts[c1]>(notu-2)) && c1<(nchars-1))	++c1;
		if (c1>=nchars)		break;
		for (c2=(c1+1); c2<nchars; ++c2)	{
			/* reinstated 2012-02-29	*/
			if (excl=='y')		while ((simauts[c2]<2 || simauts[c2]>(notu-2)) && c2<(nchars-1))	++c2;
			if (c2>=nchars)		break;

			if (compmat[c1][c2]==1)	{
				/*sstrcmp=stratcompatfull(ranges, simatrix, nstates, c1, c2, notu, UNKNOWN, INAP);	*/
				/* 2011-10-14: changed to "plus"													*/
				sstrcmp=stratcompatfullplusplusplus(ranges, simatrix, nstates, charcomps, c1, c2, notu, onset, term, tiebreaker, UNKNOWN, INAP);
				/*************************************************
					sstrcmp[0]:  character 1
					sstrcmp[1]:  character 2
					sstrcmp[2]:  character 1 states
					sstrcmp[3]:  character 2 states
					sstrcmp[4]:  state comparisons (2 for binary; possibly more for multistate)
					sstrcmp[5]:  general stratigraphic compatibility!
					sstrcmp[6]:  strict stratigraphic compatibility
					sstrcmp[7]:  super strict stratigraphic compatibility
					sstrcmp[8]:  divergent (0) vs. hierarchical (1) stratigraphic compatibility
					sstrcmp[9]:  hierarchical stratigraphic compatibility consistent with anagenesis
					sstrcmp[10]: hierarchical stratigraphic compatibility consistent with budding
					sstrcmp[11]: compatibility suggests 10<-01->11 rather than 10->01->11
					sstrcmp[12]: compatibility suggests 10->01->11 rather than 10<-01->11
					sstrcmp[13]: stratigraphy suggests 10<-01->11 rather than 10->01->11
					sstrcmp[14]: stratigraphy suggests 10->01->11 rather than 10<-01->11
					sstrcmp[15]: relative diversity of younger swing pair in SIN cases
					sstrcmp[16]: center of gravity for oldest state pair
					sstrcmp[17]: center of gravity for intermediate state pair
					sstrcmp[18]: when hierarchical ancestral pair goes extinct;
					sstrcmp[19]: when divergent ancestral pair goes extinct;
					sstrcmp[20]: purely hierarchical comparisons;
					sstrcmp[21]: purely divergent comparisons;			*******************/
				/* rewrite this so that we simply take the median values from sstrcomp (simulated strat comp) for expectations	*/
				/* and the distribution for p-values when examining the real results											*/
				cp+=sstrcmp[4];										/* compatible comparisons						*/
				genstratcomp[r]+=sstrcmp[5];						/* general stratigraphic consistency			*/
				strstratcomp[r]+=sstrcmp[6];						/* strict stratigraphic consistency				*/
				supstratcomp[r]+=sstrcmp[7];						/* superstrict stratigraphic consistency		*/
				hierstratcmp[r]+=sstrcmp[8];						/* hierarchical stratigraphic consistency		*/
				hierstratana[r]+=sstrcmp[9];						/* HSC consistent with anagenesis				*/
				hierstratbud[r]+=sstrcmp[10];						/* HSC consistent with budding					*/
				compdivSIN[r]+=sstrcmp[11];							/* compatibility-supported divergent SIN		*/
				comphierSIN[r]+=sstrcmp[12];						/* compatibility-supported hierarchical SIN		*/
				stratdivSIN[r]+=sstrcmp[13];						/* stratigraphy-supported divergent SIN			*/
				strathierSIN[r]+=sstrcmp[14];						/* stratigraphy-supported hierarchical SIN		*/
				propend2cp[r]+=sstrcmp[15];							/* relative diversity of younger SIN end pair	*/
				if (sstrcmp[5]>0)				avecompscp[r]+=(((double) charcomps[c1])+((double) charcomps[c2]))/2;
				if (sstrcmp[5]<sstrcmp[4])		avecompsin[r]+=(((double) charcomps[c1])+((double) charcomps[c2]))/2;

				/* 2013-04-14: Add stuff COG & surivorship here!  2013-04-16: just use binaries	*/
				if (sstrcmp[20]>0)	{
					aveCGH[r]+=sstrcmp[16];
					aveADH[r]+=sstrcmp[18];
					if (sstrcmp[18]==1)	++freqLFH[r];						/* tally another living fossil				*/
					hiercounter+=sstrcmp[20];								/* number of purely hierarchical tallies	*/
					}
				if (sstrcmp[21]>0)	{
					aveCGD[r]+=sstrcmp[17];
					aveADD[r]+=sstrcmp[19];
					if (sstrcmp[19]==1)	++freqLFD[r];						/* tally another living fossil	*/
					divcounter+=sstrcmp[21];								/* number of purely divergent tallies	*/
					}			/* 2013-04-15: RESUME HERE!!!!	*/

/*				scha[r]+=sstrcmp[9];						/* hierarchical compatible pairs consistent with anagenesis			2011-10-14	*/
/*				sih[r]+=sstrcmp[11];						/* stratigraphically incompatible pairs consistent with hierarchy	2011-10-14	*/
				/* 2011-10-14: tally compatibilities of stratigraphically compatible & incompatible character pairs	*/
/*				if (sstrcmp[5]==1)	{
					avecompscp[r]+=((double) (charcomps[c1]+charcomps[c2]))/2;
					scgw+=1.0;													/* this will tally characters, not state pairs						*/
/*					}
				else	{
					avecompsin[r]+=((double) (charcomps[c1]+charcomps[c2]))/2;
					sci+=1.0;
					}
	
				/* this will give the proportion of compatibilities that are stratigraphically compatible by compatibility	*/
/*				genrstcmpcor[compmat[c1][c1]]+=sstrcmp[5];			/* general stratigraphic consistency		*/
/*				strcstcmpcor[compmat[c1][c1]]+=sstrcmp[6];			/* strict stratigraphic consistency			*/
/*				hierstcmpcor[compmat[c1][c1]]+=sstrcmp[8];			/* hierarchical stratigraphic consistency	*/
/*				genrstcmpcor[compmat[c2][c2]]+=sstrcmp[5];			/* general stratigraphic consistency		*/
/*				strcstcmpcor[compmat[c2][c2]]+=sstrcmp[6];			/* strict stratigraphic consistency			*/
/*				hierstcmpcor[compmat[c2][c2]]+=sstrcmp[8];			/* hierarchical stratigraphic consistency	*/
				free_dvector(sstrcmp);
				}	/* end case of simulated compatible pair	*/
			}	/* end search of remaining simulated compatible pairs for possible compatibilities	*/
		}	/* search characters for possible compatibilities	*/
	
	/* get proportions of hierarchical stratigraphic characters relative to total								*/
	/* Also, get proportions of particular types of hierarchical stratigraphic characters						*/
	hierstratprp[r]+=(hierstratcmp[r]/genstratcomp[r]);						/* proportional hierarchical stratigraphic compatibility	*/
	/* if there is no hierarchical compatibility, then modify the routine: 2011-11-01	*/
	if (hierstratcmp[r]>0)	{
		hierstratana[r]/=hierstratcmp[r];										/* proportion of HSC matching anagenesis					*/
		hierstratbud[r]/=hierstratcmp[r];										/* proportion of HSC demanding budding						*/
		/* stick COG stuff in here???	*/
		}
	else	{
		++hierflop;																/* tally runs in which there was no hierarchical compatibility	*/
		hierstratana[r]=hierstratbud[r]=MAXRAND;
		}
	/* get average compatibilities of stratigraphically compatible and incompatible character pairs					*/
	avecompscp[r]/=genstratcomp[r];											/* added 2011-10-14	*/
	avecompsin[r]/=(cp-genstratcomp[r]);									/* added 2011-10-14	*/
	avecompdif[r]=avecompscp[r]-avecompsin[r];								/* added 2011-11-01: because I sort these later, this is otherwise lost */
	/* get proportion of SIN character pairs matching divergent or hierarchical patterns				*/
	compdivSIN[r]/=(cp-genstratcomp[r]);									/* compatibility-supported divergent SIN		*/
	comphierSIN[r]/=(cp-genstratcomp[r]);									/* compatibility-supported hierarchical SIN		*/
	stratdivSIN[r]/=(cp-genstratcomp[r]);									/* stratigraphy-supported divergent SIN			*/
	strathierSIN[r]/=(cp-genstratcomp[r]);									/* stratigraphy-supported hierarchical SIN		*/
	propend2cp[r]/=(cp-genstratcomp[r]);									/* relative diversity of younger SIN end pair	*/
	/* cp will vary a little from simulation to simulation because it is based on state numbers, so report things in proportions	*/
	genstratcomp[r]/=((double) cp);											/* general stratigraphic consistency						*/
	strstratcomp[r]/=((double) cp);											/* strict stratigraphic consistency							*/
	supstratcomp[r]/=((double) cp);											/* super strict stratigraphic consistency					*/
	hierstratcmp[r]/=((double) cp);											/* hierarchical stratigraphic consistency					*/
	
	/* 2013-04-14: Add stuff COG & surivorship here!	*/
	aveCGH[r]/=((double) hiercounter);
/*	aveCGH[r]/=sC[r];					/* standardize first by opportunities and second by total clade shape	*/
	aveADH[r]/=((double) hiercounter);
	freqLFH[r]/=((double) hiercounter);
	aveCGD[r]/=((double) divcounter);
/*	aveCGD[r]/=sC[r];					/* standardize first by opportunities and second by total clade shape	*/
	aveADD[r]/=((double) divcounter);
	freqLFD[r]/=((double) divcounter);
	
	/* tally cases where simulated values are lower than observed values		*/
	if (genstratcomp[r]<summary[0])		simmary[0]+=(1/((double) RUNS));
	if (strstratcomp[r]<summary[1])		simmary[2]+=(1/((double) RUNS));
	if (supstratcomp[r]<summary[2])		simmary[4]+=(1/((double) RUNS));
	if (hierstratcmp[r]<summary[3])		simmary[6]+=(1/((double) RUNS));
	if (hierstratprp[r]<summary[4])		simmary[8]+=(1/((double) RUNS));
	if (hierstratana[r]<summary[5])		simmary[10]+=(1/((double) RUNS));
	if (hierstratbud[r]<summary[6])		simmary[12]+=(1/((double) RUNS));
	if (compdivSIN[r]<summary[7])		simmary[14]+=(1/((double) RUNS));
	if (comphierSIN[r]<summary[8])		simmary[16]+=(1/((double) RUNS));
	if (stratdivSIN[r]<summary[9])		simmary[18]+=(1/((double) RUNS));
	if (strathierSIN[r]<summary[10])	simmary[20]+=(1/((double) RUNS));
	if (propend2cp[r]<summary[11])		simmary[22]+=(1/((double) RUNS));
	if (avecompscp[r]<summary[12])		simmary[24]+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompsin[r]<summary[13])		simmary[26]+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompdif[r]<(summary[12]-summary[13]))
										simmary[27]+=(1/((double) RUNS));					/* added 2011-11-01: because I sort these later, this is otherwise lost */
	/* 2013-04-14: Add stuff COG & surivorship here!	*/
	if (aveCGH[r]<summary[14])			simmary[30]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (aveCGD[r]<summary[15])			simmary[32]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (CG[r]<summary[16])				simmary[34]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (aveADH[r]<summary[17])			simmary[36]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (aveADD[r]<summary[18])			simmary[38]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (freqLFH[r]<summary[19])			simmary[40]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (freqLFD[r]<summary[20])			simmary[42]+=(1/((double) RUNS));					/* added 2013-04-15 */

	/* tally cases where simulated values are exactly equal to observed values	*/
	if (genstratcomp[r]==summary[0])	tiegs+=(1/((double) RUNS));
	if (strstratcomp[r]==summary[1])	tiess+=(1/((double) RUNS));
	if (supstratcomp[r]==summary[2])	tiexs+=(1/((double) RUNS));
	if (hierstratcmp[r]==summary[3])	tiehs+=(1/((double) RUNS));
	if (hierstratprp[r]==summary[4])	tiehp+=(1/((double) RUNS));
	if (hierstratana[r]==summary[5])	tieha+=(1/((double) RUNS));
	if (hierstratbud[r]==summary[6])	tiehb+=(1/((double) RUNS));
	if (compdivSIN[r]==summary[7])		tiesincd+=(1/((double) RUNS));
	if (comphierSIN[r]==summary[8])		tiesinch+=(1/((double) RUNS));
	if (stratdivSIN[r]==summary[9])		tiesinsd+=(1/((double) RUNS));
	if (strathierSIN[r]==summary[10])	tiesinsh+=(1/((double) RUNS));
	if (propend2cp[r]==summary[11])		tiesind2+=(1/((double) RUNS));
	if (avecompscp[r]==summary[12])		tieavecscp+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompsin[r]==summary[13])		tieavecsin+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompdif[r]==(summary[12]-summary[13]))
										tieavecdif+=(1/((double) RUNS));					/* added 2011-11-01: because I sort these later, this is otherwise lost */
	/* 2013-04-14: Add stuff COG & surivorship here!	*/
	if (aveCGH[r]==summary[14])			tieavecgh+=(1/((double) RUNS));						/* added 2013-04-15 */
	if (aveCGD[r]==summary[15])			tieavecgd+=(1/((double) RUNS));						/* added 2013-04-15 */
	if (CG[r]==summary[16])				tierc1c0+=(1/((double) RUNS));						/* added 2013-04-15 */
	if (aveADH[r]==summary[17])			tieaveadh+=(1/((double) RUNS));						/* added 2013-04-15 */
	if (aveADD[r]==summary[18])			tieaveadd+=(1/((double) RUNS));						/* added 2013-04-15 */
	if (freqLFH[r]==summary[19])		tiefreqlfh+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (freqLFD[r]==summary[20])		tiefreqlfd+=(1/((double) RUNS));					/* added 2013-04-15 */
	
	if ((r%10)==9)	{
//			printf("Doing Rate = %5.4f",mbl[2]);
		if (r>10)	{
			printf("\b\b\b\b\b");					/* clear ", Tree %d"		*/
			}
		if (r>18 && r<100) 			printf("\b\b\b");
		else if (r>99 && r<1000)	printf("\b\b\b\b");
		else if (r>999)				printf("\b\b\b\b\b");
		printf("Tree %d\n", r+1);
		}
	
/*	free_dvector(sC);						/* 2013-04-15		*/
	free_ulmatrix(compmat,nchars,nchars);
	free_lmatrix(vtree,clades+1,notu);
	free_lmatrix(tree,clades+3,notu);
	free_ivector(simauts);				/* found 2013-05-14	*/
	free_ulvector(charcomps);			/* found 2013-05-14	*/
	}

/* sort for medians	*/
genstratcomp=dshellsort_inc(genstratcomp,RUNS); 	/* 0:  sort stratigraphic compatiblity						*/
strstratcomp=dshellsort_inc(strstratcomp,RUNS); 	/* 1:  sort strict stratigraphic compatiblity				*/
supstratcomp=dshellsort_inc(supstratcomp,RUNS);		/* 2:  sort super strict stratigraphic compatibility		*/
hierstratcmp=dshellsort_inc(hierstratcmp,RUNS);		/* 3:  sort hierarchical stratigraphic compatiblity (HSC)	*/
hierstratprp=dshellsort_inc(hierstratprp,RUNS);		/* 4:  sort proportion HSC									*/
hierstratana=dshellsort_inc(hierstratana,RUNS);		/* 5:  sort HSC consistent with anagenesis					*/
hierstratbud=dshellsort_inc(hierstratbud,RUNS);		/* 6:  sort HSC demanding budding							*/
compdivSIN=dshellsort_inc(compdivSIN,RUNS);			/* 7:  sort compatibility suggests divergent SIN			*/
comphierSIN=dshellsort_inc(comphierSIN,RUNS);		/* 8:  sort compatibility suggests hierarchical SIN			*/
stratdivSIN=dshellsort_inc(stratdivSIN,RUNS);		/* 9:  sort stratigraphy suggests divergent SIN				*/
strathierSIN=dshellsort_inc(strathierSIN,RUNS);		/* 10: sort stratigraphy suggests hierarchical SIN			*/
propend2cp=dshellsort_inc(propend2cp,RUNS);			/* 11: sort relative diversity of derived end pair			*/
avecompscp=dshellsort_inc(avecompscp,RUNS);			/* 12: sort average compatibilities of SC character pairs	*/
avecompsin=dshellsort_inc(avecompsin,RUNS);			/* 13: sort average compatibilities of SIN character pairs	*/
avecompdif=dshellsort_inc(avecompdif,RUNS);			/* 14: sort average differences between SC & SIN pairs		*/
aveCGH=dshellsort_inc(aveCGH,RUNS);					/* 15: sort average 00 CoGs for hierarchical pairs		*/
aveCGD=dshellsort_inc(aveCGD,RUNS);					/* 16: sort average 00 CoGs for divergent pairs			*/
aveADH=dshellsort_inc(aveADH,RUNS);					/* 17: sort average 00 durations for hierarchical pairs	*/
aveADD=dshellsort_inc(aveADD,RUNS);					/* 18: sort average 00 durations for divergent pairs	*/
CG=dshellsort_inc(CG,RUNS);							/* 19: sort average CoG for simulated clades				*/
freqLFH=dshellsort_inc(freqLFH,RUNS);				/* 20: sort hierarchical living fossils						*/
freqLFD=dshellsort_inc(freqLFD,RUNS);				/* 20: sort hierarchical living fossils						*/

/* get alpha values: if observed values are "low" then we want p[observed or lower]: so, add ties to total	*/
/* for "high" values, we do not need to do this: 1-simmary gives p[observed or more extreme] as we tallied only sims less than observed	*/
/* we want the probability of the observed OR MORE EXTREME; because we tallied every time the observed was lower than what we saw,
		we don't need to do this when p> 0.5	*/
if (simmary[0]<0.5)		simmary[0]+=tiegs;			/* general stratigraphic compatibility	*/
if (simmary[2]<0.5)		simmary[2]+=tiess;			/* strict stratigraphic compatibility	*/
if (simmary[4]<0.5)		simmary[4]+=tiexs;
if (simmary[6]<0.5)		simmary[6]+=tiehs;
if (simmary[8]<0.5)		simmary[8]+=tiehp;
if (simmary[10]<0.5)	simmary[10]+=tieha;
if (simmary[12]<0.5)	simmary[12]+=tiehb;
if (simmary[14]<0.5)	simmary[14]+=tiesincd;
if (simmary[16]<0.5)	simmary[16]+=tiesinch;
if (simmary[18]<0.5)	simmary[18]+=tiesinsd;
if (simmary[20]<0.5)	simmary[20]+=tiesinsh;
if (simmary[22]<0.5)	simmary[22]+=tiesind2;
if (simmary[24]<0.5)	simmary[24]+=tieavecscp;
if (simmary[26]<0.5)	simmary[26]+=tieavecsin;
if (simmary[28]<0.5)	simmary[28]+=tieavecsin;
if (simmary[30]<0.5)	simmary[30]+=tieavecgh;
if (simmary[32]<0.5)	simmary[32]+=tieavecgd;
if (simmary[34]<0.5)	simmary[34]+=tieaveadh;
if (simmary[36]<0.5)	simmary[36]+=tieaveadd;
if (simmary[38]<0.5)	simmary[38]+=tierc1c0;
if (simmary[40]<0.5)	simmary[40]+=tiefreqlfh;
if (simmary[42]<0.5)	simmary[42]+=tiefreqlfd;

/* adjust for runs with no hierarchical compatibility	*/
simmary[10]*=(((double) (RUNS-hierflop))/((double) RUNS));		/* adjust significance of anagenetic hierarchical compatibility	*/
simmary[12]*=(((double) (RUNS-hierflop))/((double) RUNS));		/* adjust significance of budding hierarchical compatibility	*/


/* tally median expected stratigraphic compatibilities	*/
simmary[1]=genstratcomp[RUNS/2];
simmary[3]=strstratcomp[RUNS/2];
simmary[5]=supstratcomp[RUNS/2];
simmary[7]=hierstratcmp[RUNS/2];
simmary[9]=hierstratprp[RUNS/2];
simmary[11]=hierstratana[(RUNS-hierflop)/2];			/* accommodate runs with no hierarchical compatibility	*/
simmary[13]=hierstratbud[(RUNS-hierflop)/2];			/* accommodate runs with no hierarchical compatibility	*/
simmary[15]=compdivSIN[RUNS/2];
simmary[17]=comphierSIN[RUNS/2];
simmary[19]=stratdivSIN[RUNS/2];
simmary[21]=strathierSIN[RUNS/2];
simmary[23]=propend2cp[RUNS/2];
simmary[25]=avecompscp[RUNS/2];
simmary[27]=avecompsin[RUNS/2];
simmary[29]=avecompdif[RUNS/2];
simmary[31]=aveCGH[RUNS/2];
simmary[33]=aveCGD[RUNS/2];
simmary[35]=CG[RUNS/2];
simmary[37]=aveADH[RUNS/2];
simmary[39]=aveADD[RUNS/2];
simmary[41]=freqLFH[RUNS/2];
simmary[43]=freqLFD[RUNS/2];

/* return output to here!	*/
/* TEMPORARY!  UNTIL I'VE DEBUGGED THIS.....	*/
keep=1;
//keep=0;
//printf("Enter '1' if you want output for this taxon: ");
//scanf("%i",&keep);	

if (keep==1)	{
	strcpy(outfile,taxonname);
	strcat(outfile,"_");
	strcat(outfile,citation);
	strcat(outfile,"_Expectations");
	if (mbl[4]==0)	strcat(outfile,"_Bud.xls");
	else			strcat(outfile,"_Bif.xls");
	
	output=fopen(outfile,"w");
	fprintf(output,"Run");
	fprintf(output,"\tGeneral Stratigraphic Compatibility (SC)");
	fprintf(output,"\tStrict SC");
	fprintf(output,"\tSuper Strict SC");
	fprintf(output,"\tHierarchical SC");
	fprintf(output,"\tProportional HSC");
	fprintf(output,"\tAnagenetic HSC");
	fprintf(output,"\tBudding HSC");
	fprintf(output,"\tDivergent SIN (Comp)");
	fprintf(output,"\tHierarchical SIN (Comp)");
	fprintf(output,"\tDivergent SIN (Strat)");
	fprintf(output,"\tHierarchical SIN (Strat)");
	fprintf(output,"\tProp Diversity Younger SIN Pair");
	fprintf(output,"\tµ Character Compabibility (SC)");					/* added 2011-10-14	*/
	fprintf(output,"\tµ Character Compabibility (SIN)");				/* added 2011-10-14	*/
	fprintf(output,"\tµ Compabibility Diff (SC-SIN)");					/* added 2011-11-01	*/
	fprintf(output,"\tµ H00 CoG");										/* added 2013-04-16	*/
	fprintf(output,"\tµ D00 CoG");										/* added 2013-04-16	*/
	fprintf(output,"\tµ Clade CoG");									/* added 2013-04-16	*/
	
	for (r=0; r<RUNS; ++r)	{
		fprintf(output,"%d",r+1);
		fprintf(output,"\t%4.3f",genstratcomp[r]);
		fprintf(output,"\t%4.3f",strstratcomp[r]);
		fprintf(output,"\t%4.3f",supstratcomp[r]);
		fprintf(output,"\t%4.3f",hierstratcmp[r]);
		fprintf(output,"\t%4.3f",hierstratprp[r]);
		if (hierstratana[r]<MAXRAND)	fprintf(output,"\t%4.3f",hierstratana[r]);
		else							fprintf(output,"\t•");
		if (hierstratbud[r]<MAXRAND)	fprintf(output,"\t%4.3f",hierstratbud[r]);
		else							fprintf(output,"\t•");
		fprintf(output,"\t%4.3f",compdivSIN[r]);
		fprintf(output,"\t%4.3f",comphierSIN[r]);
		fprintf(output,"\t%4.3f",stratdivSIN[r]);
		fprintf(output,"\t%4.3f",strathierSIN[r]);
		fprintf(output,"\t%4.3f",propend2cp[r]);
		fprintf(output,"\t%4.3f",avecompscp[r]);					/* added 2011-10-14	*/
		fprintf(output,"\t%4.3f",avecompsin[r]);					/* added 2011-10-14	*/
		fprintf(output,"\t%4.3f",avecompdif[r]);					/* added 2011-11-01	*/
		fprintf(output,"\t%4.3f",aveCGH[r]);						/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f",aveCGD[r]);						/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f",CG[r]);							/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f",aveADH[r]);						/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f",aveADD[r]);						/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f",aveADD[r]);						/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f",freqLFH[r]);						/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f\n",freqLFD[r]);						/* added 2013-04-16	*/
		}
	fclose(output);
	}


free_dvector(aveADD);				/* tallies for simmary[32,33]							*/	
free_dvector(aveADH);				/* tallies for simmary[30,31]							*/	
free_dvector(aveCGD);				/* tallies for simmary[32,33]							*/	
free_dvector(aveCGH);				/* tallies for simmary[30,31]							*/	
free_dvector(avecompdif);			/* tallies for simmary[28,29]							*/		
free_dvector(avecompscp);			/* tallies for simmary[24,25]							*/		
free_dvector(avecompsin);			/* tallies for simmary[26.27]							*/		
free_dvector(CG);					/* tallies for simmary[38,39]							*/
free_dvector(compdivSIN);			/* tallies sstrcomp[11], for simmary[14,15]				*/					
free_dvector(comphierSIN);			/* tallies sstrcomp[12], for simmary[16,17]				*/					
free_dvector(freqLFD);				/* tallies for simmary[36,37]							*/	
free_dvector(freqLFH);				/* tallies for simmary[34,35]							*/	
free_dvector(genstratcomp);			/* tallies sstrcomp[5], for simmary[0,1]				*/					
free_dvector(hierstratana);			/* tallies sstrcomp[9], for simmary[10,11]				*/					
free_dvector(hierstratbud);			/* tallies sstrcomp[10], for simmary[12,13]				*/					
free_dvector(hierstratcmp);			/* tallies sstrcomp[8], for simmary[4,5]				*/					
free_dvector(hierstratprp);			/* tallies sstrcomp[8]/sstrcomp[5], for simmary[6,7]	*/								
free_dvector(propend2cp);			/* tallies sstrcomp[15], for simmary[22,23]				*/					
free_dvector(stratdivSIN);			/* tallies sstrcomp[13], for simmary[18,19]				*/					
free_dvector(strathierSIN);			/* tallies sstrcomp[14], for simmary[20,21]				*/					
free_dvector(strstratcomp);			/* tallies sstrcomp[6], for simmary[2,3]				*/					
free_dvector(supstratcomp);			/* tallies sstrcomp[7], for simmary[8,9]				*/					

free_lmatrix(ranges,notu,2);		/* simulated ranges										*/
free_lmatrix(simatrix,notu,nchars);	/* simulated morphologies								*/


return simmary;
}

	
/* 	fa: observed first appearance
	origin: hypothesized origin
	r: vector giving preservation rates in different intervals
	vr: vector giving lognormal variation in preservation rates in different intervals
******************************************************************************************/	
double stratlikelihoodbranch(double fa, double origin, double *r, double *vr)
{
int a, b, c, d;
double	t;
double	x, y, z;
double	divisions[4];
double L=0.0f, lnL=0.0f;

a=1+fa;			/* bin of observed first appearance; note that bins go from .00 to 0.9999…		*/
b=1+origin;		/* bin of hypothesized origin; note that bins go from .00 to 0.9999…			*/

/* these are the divisions of a normal curve giving midpoints of 4 equal area units	*/
divisions[0]=-1.150349380;
divisions[1]=-0.318639364;
divisions[2]=0.318639364;
divisions[3]=1.150349380;

/* go through each bin (c) from first appearance (a) to hypothesized origin (b)	*/
for (c=a; c>=b; --c)	{
	if (c==a && c>b)
		t=fa-((double) (c-1));		/* take into account partial presence in last "stage"	*/
	else if (c==b && c<a)
		t=((double) c)-origin;		/* take into account partial presence in first "stage"	*/
	else if (a==b)
		t=fa-origin;				/* take into account partial presence in only "stage"	*/
	else
		t=1.0f;						/* presence in whole "stage"	*/

	if (t>0)	{
		L=0.0f;
		for (d=0; d<4; ++d)	{
			x=pow(vr[c],divisions[d]);	/* rate modifier at this point of lognormal	*/
			y=x*r[c];					/* modified rate							*/
			z=y*t;						/* expected finds							*/
			L+=pow(e,-z);				/* probability of 0 finds at this rate		*/
			}
		L/=4;			/* take into account our 1 in 4 priors on the different preservation rates	*/
		lnL+=log(L);	/* log it	*/
		}
	}
return lnL;
}

/* 2012-12-29: PICKUP HERE!!!	*/
double gaplikelihoodgeog (double top, double bottom, double **preservation, double *timescale, double *geogpr, int stages, int areas, int rates)
{
int	a, r;
int up, lo, st;
double x, t;
double lnlgap=0.0f;
double m, n, p;

/* determine where gaps are	in time	*/
/* make sure that these handle tops and bottoms at boundaries	*/
/* time scale gives the end of the stage.  So, if the top = timescale[st],
	then st is the upper end.
	If bottom = timescale[st+1], then st is the bottom				*/
for (st=0; st<stages; ++st)	{
	if (top>=timescale[st] && top<timescale[st+1])
		up=st;
	if (bottom>timescale[st] && bottom<=timescale[st+1])	{
		lo=st;
		st=stages;
		}
	}

for (st=up; st<=lo; ++st)	{
	/* three possible cases:
		1. lower part of the uppermost stage
		2. an entire stage
		3. upper part of the lowermost stage */
	x=0.0f;
	if (st==up)	{
		if (st==lo)		t=bottom-top;
		else			t=timescale[st+1]-top;
		}
	else if (st==lo)	t=bottom-timescale[st];
	else				t=timescale[st+1]-timescale[st];

	/* now tally probability of the gaps over different areas	*/
	for (a=0; a<areas; ++a)	{
		for (r=0; r<rates; ++r)	{
			m=preservation[st][(rates*a)+r];
			n=1-preservation[st][(rates*a)+r];
			p=pow((1-preservation[st][(rates*a)+r]),t);
			x+=(geogpr[a]*pow((1-preservation[st][(rates*a)+r]),t)/((double) rates));
			}	/* go throuh all recovery rates, with rates on a per time basis	*/
		}	/* go through all areas that might be covered in a gap	*/
	lnlgap+=log(x);
	}	/* go through all stages in a gap	*/

return (lnlgap);
} 

/* 2013-01-13: PICKUP HERE!!!	*/
void gaplikelihoodgeogsep (double **geogstrlnl, double top, double bottom, double **preservation, double *timescale, int stages, int areas, int rates, int taxon)
{
int	a, r;
int up, lo, st;
double x, t;
//double *lnlgap;
double m, n, p;

//lnlgap=dvector(areas);

/* determine where gaps are	in time	*/
/* make sure that these handle tops and bottoms at boundaries	*/
/* time scale gives the end of the stage.  So, if the top = timescale[st],
	then st is the upper end.
	If bottom = timescale[st+1], then st is the bottom				*/
for (st=0; st<stages; ++st)	{
	if (top>=timescale[st] && top<timescale[st+1])
		up=st;
	if (bottom>timescale[st] && bottom<=timescale[st+1])	{
		lo=st;
		st=stages;
		}
	}

for (a=0; a<areas; ++a) geogstrlnl[taxon][a]=0.0f;

for (st=up; st<=lo; ++st)	{
	/* three possible cases:
		1. lower part of the uppermost stage
		2. an entire stage
		3. upper part of the lowermost stage */
	if (st==up)	{
		if (st==lo)		t=bottom-top;
		else			t=timescale[st+1]-top;
		}
	else if (st==lo)	t=bottom-timescale[st];
	else				t=timescale[st+1]-timescale[st];

	/* now tally probability of the gaps over different areas	*/
	for (a=0; a<areas; ++a)	{
		x=0.0f;			/* this will sum the posterior probabilities */
		for (r=0; r<rates; ++r)	{
			m=preservation[st][(rates*a)+r];
			n=1-preservation[st][(rates*a)+r];
			p=pow((1-preservation[st][(rates*a)+r]),t);
			x+=(pow((1-preservation[st][(rates*a)+r]),t)/((double) rates));
			}	/* go throuh all recovery rates, with rates on a per time basis	*/
		geogstrlnl[taxon][a]+=log(x);
		}	/* go through all areas that might be covered in a gap	*/
	}	/* go through all stages in a gap	*/

/*return (0);	*/
}


/* stratcompatparambstest3: simulations for stratigraphic compatibility with parameters determined by best fit for stratigraphic data
Requires:
	taxonname: group analyzed
	citation: source of data
	summary: summary of observed data
/*		summary[0]:  prop. of generally compatible pairs										2011-10-15	*/
/*		summary[1]:  prop. of strictly compatible pairs											2011-10-15	*/
/*		summary[2]:  prop. of super strictly compatible pairs									2011-10-15	*/
/*		summary[3]:  prop. of hierarchically compatible pairs									2011-10-15	*/
/*		summary[4]:  prop. of hierarchically compatible pairs given strat. consistency			2011-10-15	*/
/*		summary[5]:  prop. of hierarchically compatible pairs consistent with anagenesis		2011-10-15	*/
/*		summary[6]:  prop. of hierarchically compatible pairs consistent with budding			2011-10-15	*/
/*		summary[7]:  prop. of SIN pairs consistent with divergence given compatibility			2011-10-15	*/
/*		summary[8]:  prop. of SIN pairs consistent with hierarchy given compatibility			2011-10-15	*/
/*		summary[9]:  prop. of SIN pairs consistent with divergence given stratigraphy			2011-10-15	*/
/*		summary[10]: prop. of SIN pairs consistent with hierarchy given stratigraphy			2011-10-15	*/
/*		summary[11]: average relative diversity of derived end-pair for SIN pairs				2011-10-15	*/
/*		summary[12]: average compatibility of stratigraphically compatible pairs				2011-10-15	*/
/*		summary[13]: average compatibility of stratigraphically incompatible pairs				2011-10-15	*/
/*		summary[14]: average center of gravity for hierarchical 00								2013-04-14	*/
/*		summary[15]: average center of gravity for divergent 00									2013-04-14	*/
/*		summary[16]: Center of gravity for whole clade											2013-04-14	*/
/*		summary[17]: average duration of hierarchial 00's										2013-04-14	*/
/*		summary[18]: average duration of divergent 00's											2013-04-14	*/
/*		summary[19]: "living fossil" hierarchial 00's											2013-04-15	*/
/*		summary[20]: "living fossil" divergent 00's												2013-04-15	*/
/*		summary[21]: "living fossil" divergent 00's												2013-04-15	*/
/*		summary[22]: disparity at 25% of species												2014-01-14	*/
/*		summary[23]: disparity at 50% of species												2014-01-14	*/
/*		summary[24]: disparity at 75% of species												2014-01-14	*/
/*		summary[25]: disparity at all of species												2014-01-14	*/

/*	mbl: simulation parameters
	omatrix: observed matrix
	ctype: character types
	nstates: states per character
	bias: character evolution bias
	maxch: maximum change per character
	depend: character dependencies
	notu: number of taxa
	nchars: number of characters
	compat: compatibility type (general or hierarchical)
	RUNS: simulations
	debug: 1 to seed based on replications, 0 for time
	excl: exclude autapomorphies
	UNKNOWN: code for unknown states
	INAP: code for inapplicable states
Returns:
	simmary[0]: P[observed general stratigraphic compatibility]
	simmary[1]: E[general stratigraphic compatibility]
	simmary[2]: P[observed strict stratigraphic compatibility]
	simmary[3]: E[strict stratigraphic compatibility]
	simmary[4]: P[observed super strict compatible character pairs]
	simmary[5]: E[super strict compatible character pairs]
	simmary[6]: P[observed hierarchical stratigraphic compatibility]
	simmary[7]: E[hierarchical stratigraphic compatibility]
	simmary[8]: P[observed proportion hierarchical stratigraphic compatibility;]
	simmary[9]: E[proportion hierarchical stratigraphic compatibility;]
	simmary[10]: P[observed hierarchical stratigraphic compatible pairs consistent with anagenesis]
	simmary[11]: E[hierarchical stratigraphic compatible pairs consistent with anagenesis]
	simmary[12]: P[observed hierarchical stratigraphic compatible pairs consistent with budding]
	simmary[13]: E[hierarchical stratigraphic compatible pairs consistent with budding]
	simmary[14]: P[observed compatibilty implies divergent stratigraphically incompatibilty (SIN)]
	simmary[15]: E[compatibilty implies divergent SIN]
	simmary[16]: P[observed compatibilty implies hierarchical stratigraphically incompatibilty SIN]
	simmary[17]: E[stratigraphy implies hierarchical SIN]
	simmary[18]: P[observed stratigraphy implies divergent SIN]
	simmary[19]: E[stratigraphy implies divergent SIN]
	simmary[20]: P[observed stratigraphy implies hierarchical SIN]
	simmary[21]: E[stratigraphy implies hierarchical SIN]
	simmary[22]: P[observed proportion of younger end state pair in SIN characters]
	simmary[23]: E[proportion of younger end state pair in SIN characters]
	simmary[24]: P[compatibility of stratigraphically compatible pairs]
	simmary[25]: E[compatibility of stratigraphically compatible pairs]
	simmary[26]: P[compatibility of SIN pairs]
	simmary[27]: E[compatibility of SIN pairs]
	simmary[28]: P[diff. in compatibility between SC & SIN pairs]
	simmary[29]: E[diff. in compatibility between SC & SIN pairs]
	simmary[30]: P[clade shape 00 when hierarchical]
	simmary[31]: E[clade shape 00 when hierarchical]
	simmary[32]: P[clade shape 00 when divergent]
	simmary[33]: E[clade shape 00 when divergent]
	simmary[34]: P[overall clade shape]
	simmary[35]: E[overall clade shape]
	simmary[36]: P[hierarchical 00 duration]
	simmary[37]: E[hierarchical 00 duration]
	simmary[38]: P[divergent 00 duration]
	simmary[39]: E[divergent 00 duration]
	simmary[40]: P[hierarchical "living fossils"]
	simmary[41]: E[hierarchical "living fossils"]
	simmary[42]: P[divergent "living fossils"]
	simmary[43]: E[divergent "living fossils"]
	simmary[44]: P[observed time to qet 1/4 diversity]
	simmary[45]: E[observed time to qet 1/4 diversity]
	simmary[46]: P[observed time to qet 1/2 diversity]
	simmary[47]: E[observed time to qet 1/2 diversity]
	simmary[48]: P[observed time to qet 3/4 diversity]
	simmary[48]: E[observed time to qet 3/4 diversity]
	simmary[50]: P[observed disparity at S/2]
	simmary[51]: E[observed disparity at S/2]
**********************************************************************************************************/
double *stratcompatparambstest3(char taxonname[90], char citation[90], double *summary, double *mbl, long **omatrix, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, char excl, int tiebreaker, int debug, int UNKNOWN, int INAP)
{
/* part one - matrix properties */
int		d, r, s;
int		c, c1, c2;
int		keep;
int 	clades;
int		hierflop=0;
int		*simauts;
long  	**simatrix, **ranges;
long	**tree, **vtree;
double	*sstrcmp;
unsigned long **compmat;
unsigned long	*charcomps;
double	*genstratcomp, *strstratcomp, *hierstratcmp, *hierstratprp;
double	*supstratcomp, *hierstratana, *hierstratbud;
double	*compdivSIN, *comphierSIN, *stratdivSIN, *strathierSIN;
double	*propend2cp, *avecompscp, *avecompsin, *avecompdif;
double	*simmary;
double	cp;
double	tiegs=0.0f, tiess=0.0f, tiehs=0.0f, tiehp=0.0f, tiescha=0.0f, tiesih=0.0f, tieavecscp=0.0f, tieavecsin=0.0f, tieavecdif=0.0f;
double	tiexs=0.0f, tieha=0.0f, tiehb=0.0f, tiesincd=0.0f, tiesinch=0.0f, tiesinsd=0.0f, tiesinsh=0.0f, tiesind2=0.0f;
double	tieq1=0.0f, tieq2=0.0f, tieq3=0.0f, tiehalftaxa=0.0f;
long	secs;
char	outfile[120];
FILE	*output, *debugoutput;
/* 2013-04-13: for center of gravity madness	*/
int		term, onset;				
double	*CG, *aveCGH, *aveCGD, *aveADH, *aveADD, *freqLFH, *freqLFD;
double **PWDis, *sdispsummary;
double	tieavecgh=0.0f, tieavecgd=0.0f, tieaveadh=0.0f, tieaveadd=0.0f, tierc1c0=0.0f, tiefreqlfh=0.0f, tiefreqlfd=0.0f;
int	hiercounter, divcounter;	/* cases of hierarchy or divergence clearly	*/
double *quart1, *quart2, *quart3, *halftaxa;


if (RUNS%2==0)	++RUNS;			/* it's easier to use an odd number....	*/

time(&secs);
srand((unsigned int) secs);

genstratcomp=dvector(RUNS);		/* tallies sstrcomp[5], for simmary[0,1]				*/
strstratcomp=dvector(RUNS);		/* tallies sstrcomp[6], for simmary[2,3]				*/	
supstratcomp=dvector(RUNS);		/* tallies sstrcomp[7], for simmary[8,9]				*/
hierstratcmp=dvector(RUNS);		/* tallies sstrcomp[8], for simmary[4,5]				*/
hierstratprp=dvector(RUNS);		/* tallies sstrcomp[8]/sstrcomp[5], for simmary[6,7]	*/
hierstratana=dvector(RUNS);		/* tallies sstrcomp[9], for simmary[10,11]				*/
hierstratbud=dvector(RUNS);		/* tallies sstrcomp[10], for simmary[12,13]				*/
compdivSIN=dvector(RUNS);		/* tallies sstrcomp[11], for simmary[14,15]				*/
comphierSIN=dvector(RUNS);		/* tallies sstrcomp[12], for simmary[16,17]				*/
stratdivSIN=dvector(RUNS);		/* tallies sstrcomp[13], for simmary[18,19]				*/
strathierSIN=dvector(RUNS);		/* tallies sstrcomp[14], for simmary[20,21]				*/
propend2cp=dvector(RUNS);		/* tallies sstrcomp[15], for simmary[22,23]				*/
avecompscp=dvector(RUNS);		/* tallies for simmary[24,25]							*/
avecompsin=dvector(RUNS);		/* tallies for simmary[26.27]							*/
avecompdif=dvector(RUNS);		/* tallies for simmary[28,29]							*/
aveCGH=dvector(RUNS);			/* tallies for simmary[30,31]							*/
aveCGD=dvector(RUNS);			/* tallies for simmary[32,33]							*/
aveADH=dvector(RUNS);			/* tallies for simmary[30,31]							*/
aveADD=dvector(RUNS);			/* tallies for simmary[32,33]							*/
freqLFH=dvector(RUNS);			/* tallies for simmary[34,35]							*/
freqLFD=dvector(RUNS);			/* tallies for simmary[36,37]							*/
CG=dvector(RUNS);				/* tallies for simmary[38,39]							*/
quart1=dvector(RUNS);			/* tallies for simmary[40,41]							*/
quart2=dvector(RUNS);			/* tallies for simmary[42,43]							*/
quart3=dvector(RUNS);			/* tallies for simmary[44,45]							*/
halftaxa=dvector(RUNS);			/* tallies for simmary[46,47]							*/

ranges=lmatrix(notu,2);			/* simulated ranges										*/
simatrix=lmatrix(notu,nchars);	/* simulated morphologies								*/
simmary=dvector(52);			/* OUTPUT												*/
	
d=0;							/* for debugging	*/
for (r=d; r<RUNS; ++r)	{
	/* 2011-03-28: BLOWOUT at r=5	*/
	genstratcomp[r]=strstratcomp[r]=hierstratcmp[r]=hierstratprp[r]=0.0f;
	supstratcomp[r]=hierstratana[r]=hierstratbud[r]=compdivSIN[r]=comphierSIN[r]=stratdivSIN[r]=strathierSIN[r]=propend2cp[r]=avecompscp[r]=avecompsin[r]=0.0f;
	CG[r]=aveCGH[r]=aveCGD[r]=aveADH[r]=aveADD[r]=0.0f;
	
	for (s=0; s<notu; ++s)	{
		for (c=0; c<nchars; ++c)	{
			if (omatrix[s][c]==UNKNOWN)
				simatrix[s][c]=UNKNOWN;
			else if (omatrix[s][c]==INAP)
				simatrix[s][c]=INAP;
			else
				simatrix[s][c]=0;
			}
		}
	if (debug==1)	{
		srand((unsigned int) (r+1)*notu*nchars);
		}
	
	/* 	trees[notu-1]: branch length of otus	
		trees[notu]:   branch length of clades
		trees[notu+1]: first appearances of otus
		trees[notu+2]: last appearances of otus
	**********************************************/
	tree=evolvetree(notu,mbl,1);
	clades=cladecountbytaxa(tree,notu);
	/* pull range data out of back of tree matrix	*/
	for (s=0; s<notu; ++s)	{
		ranges[s][0]=tree[notu+1][s];
		ranges[s][1]=tree[notu+2][s];
		}
	cleancladerangedata(ranges,notu);
	
	/* set of center of gravity madness	*/
	onset=minlmatrixcol(ranges,notu,0);							/* first stage in which something appears						*/
	term=maxlmatrixcol(ranges,notu,0);							/* last stage in which something appears						*/
//	sC=cladeshapefromranges(ranges,notu,0);						/* vector with c0, c1, c2 & c3 for simulated clade				*/
//	CG[r]=sC[1]/sC[0];											/* scalar = c1/c0 for simulated clade (corrected 2013-04-16		*/
	CG[r]=centerofgravityfromranges(ranges,notu,0,onset,term);	/* scalar = c1/c0 for simulated clade (corrected 2013-04-16		*/
	if (CG[r]==summary[16])	++tierc1c0;							/* observed center of gravity is matched						*/

	vtree=VennTreePlus(tree,clades,notu,notu);
	/* evolve character matrix, simatrix	*/
//	long **evolvetocompat(long **tree, int tcomp, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int *depend, int comptype, int UNKNOWN, int INAP)
//	simatrix=evolvetocompat(tree,empcompat,notu,matrix,nchars,nstates,ctype,bias,chmax,depend,comptype,UNKNOWN,INAP);
	simatrix=evolvetocompatibility(vtree, compat, notu, clades, simatrix, nchars, nstates, ctype, bias, maxch, depend, 0, UNKNOWN, INAP);
	
	if (d>0)	{
		debugoutput=fopen("SuspectMatrix.txt","w");
		for (s=0; s<notu; ++s)	{
			for (c=0; c<nchars; ++c)	{
				if (simatrix[s][c]==UNKNOWN)	fprintf(debugoutput,"?\t");
				else if (simatrix[s][c]==INAP)	fprintf(debugoutput,"—\t");
				else							fprintf(debugoutput,"%d\t",simatrix[s][c]);
				}
			fprintf(debugoutput,"%d\t%d\n",ranges[s][0],ranges[s][1]);
			}
		fclose(debugoutput);
		}
//	simatrix=evolvematrix(vtree, notu, clades, simatrix, nchars, nstates, ctype, bias, maxch, 2*nchars, UNKNOWN, INAP);
/*		sstates=numberstates(simatrix,notu,nchars,UNKNOWN,INAP);	*/
//	compat=nu_comp(nstates, notu, simatrix, ctype, nchars, 0, 0, UNKNOWN, INAP);
	compmat=compatible(nstates,notu,simatrix,ctype,nchars,0,0,UNKNOWN,INAP);
	charcomps=char_comp(nstates,notu,simatrix,ctype,nchars,0,0,UNKNOWN,INAP);
	simauts=autapomorphies(simatrix,nstates,notu,nchars,UNKNOWN,INAP);	/* get the minimum number of derived states	*/

/*	for (c1=0; c1<nchars; ++c1)	++examples[compmat[c1][c1]];	*/
	cp=0.0f;
	hiercounter=divcounter=0;
	for (c1=0; c1<(nchars-1); ++c1)	{
		/* reinstated 2012-02-29	*/
		if (excl=='y')		while ((simauts[c1]<2 || simauts[c1]>(notu-2)) && c1<(nchars-1))	++c1;
		if (c1>=nchars)		break;
		for (c2=(c1+1); c2<nchars; ++c2)	{
			/* reinstated 2012-02-29	*/
			if (excl=='y')		while ((simauts[c2]<2 || simauts[c2]>(notu-2)) && c2<(nchars-1))	++c2;
			if (c2>=nchars)		break;

			if (compmat[c1][c2]==1)	{
				/*sstrcmp=stratcompatfull(ranges, simatrix, nstates, c1, c2, notu, UNKNOWN, INAP);	*/
				/* 2011-10-14: changed to "plus"													*/
				sstrcmp=stratcompatfullplusplusplus(ranges, simatrix, nstates, charcomps, c1, c2, notu, onset, term, tiebreaker, UNKNOWN, INAP);
				/*************************************************
					sstrcmp[0]:  character 1
					sstrcmp[1]:  character 2
					sstrcmp[2]:  character 1 states
					sstrcmp[3]:  character 2 states
					sstrcmp[4]:  state comparisons (2 for binary; possibly more for multistate)
					sstrcmp[5]:  general stratigraphic compatibility!
					sstrcmp[6]:  strict stratigraphic compatibility
					sstrcmp[7]:  super strict stratigraphic compatibility
					sstrcmp[8]:  divergent (0) vs. hierarchical (1) stratigraphic compatibility
					sstrcmp[9]:  hierarchical stratigraphic compatibility consistent with anagenesis
					sstrcmp[10]: hierarchical stratigraphic compatibility consistent with budding
					sstrcmp[11]: compatibility suggests 10<-01->11 rather than 10->01->11
					sstrcmp[12]: compatibility suggests 10->01->11 rather than 10<-01->11
					sstrcmp[13]: stratigraphy suggests 10<-01->11 rather than 10->01->11
					sstrcmp[14]: stratigraphy suggests 10->01->11 rather than 10<-01->11
					sstrcmp[15]: relative diversity of younger swing pair in SIN cases
					sstrcmp[16]: center of gravity for oldest state pair
					sstrcmp[17]: center of gravity for intermediate state pair
					sstrcmp[18]: when hierarchical ancestral pair goes extinct;
					sstrcmp[19]: when divergent ancestral pair goes extinct;
					sstrcmp[20]: purely hierarchical comparisons;
					sstrcmp[21]: purely divergent comparisons;			*******************/
				/* rewrite this so that we simply take the median values from sstrcomp (simulated strat comp) for expectations	*/
				/* and the distribution for p-values when examining the real results											*/
				cp+=sstrcmp[4];										/* compatible comparisons						*/
				genstratcomp[r]+=sstrcmp[5];						/* general stratigraphic consistency			*/
				strstratcomp[r]+=sstrcmp[6];						/* strict stratigraphic consistency				*/
				supstratcomp[r]+=sstrcmp[7];						/* superstrict stratigraphic consistency		*/
				hierstratcmp[r]+=sstrcmp[8];						/* hierarchical stratigraphic consistency		*/
				hierstratana[r]+=sstrcmp[9];						/* HSC consistent with anagenesis				*/
				hierstratbud[r]+=sstrcmp[10];						/* HSC consistent with budding					*/
				compdivSIN[r]+=sstrcmp[11];							/* compatibility-supported divergent SIN		*/
				comphierSIN[r]+=sstrcmp[12];						/* compatibility-supported hierarchical SIN		*/
				stratdivSIN[r]+=sstrcmp[13];						/* stratigraphy-supported divergent SIN			*/
				strathierSIN[r]+=sstrcmp[14];						/* stratigraphy-supported hierarchical SIN		*/
				propend2cp[r]+=sstrcmp[15];							/* relative diversity of younger SIN end pair	*/
				if (sstrcmp[5]>0)				avecompscp[r]+=(((double) charcomps[c1])+((double) charcomps[c2]))/2;
				if (sstrcmp[5]<sstrcmp[4])		avecompsin[r]+=(((double) charcomps[c1])+((double) charcomps[c2]))/2;

				/* 2013-04-14: Add stuff COG & surivorship here!  2013-04-16: just use binaries	*/
				if (sstrcmp[20]>0)	{
					aveCGH[r]+=sstrcmp[16];
					aveADH[r]+=sstrcmp[18];
					if (sstrcmp[18]==1)	++freqLFH[r];						/* tally another living fossil				*/
					hiercounter+=sstrcmp[20];								/* number of purely hierarchical tallies	*/
					}
				if (sstrcmp[21]>0)	{
					aveCGD[r]+=sstrcmp[17];
					aveADD[r]+=sstrcmp[19];
					if (sstrcmp[19]==1)	++freqLFD[r];						/* tally another living fossil	*/
					divcounter+=sstrcmp[21];								/* number of purely divergent tallies	*/
					}			/* 2013-04-15: RESUME HERE!!!!	*/

/*				scha[r]+=sstrcmp[9];						/* hierarchical compatible pairs consistent with anagenesis			2011-10-14	*/
/*				sih[r]+=sstrcmp[11];						/* stratigraphically incompatible pairs consistent with hierarchy	2011-10-14	*/
				/* 2011-10-14: tally compatibilities of stratigraphically compatible & incompatible character pairs	*/
/*				if (sstrcmp[5]==1)	{
					avecompscp[r]+=((double) (charcomps[c1]+charcomps[c2]))/2;
					scgw+=1.0;													/* this will tally characters, not state pairs						*/
/*					}
				else	{
					avecompsin[r]+=((double) (charcomps[c1]+charcomps[c2]))/2;
					sci+=1.0;
					}
	
				/* this will give the proportion of compatibilities that are stratigraphically compatible by compatibility	*/
/*				genrstcmpcor[compmat[c1][c1]]+=sstrcmp[5];			/* general stratigraphic consistency		*/
/*				strcstcmpcor[compmat[c1][c1]]+=sstrcmp[6];			/* strict stratigraphic consistency			*/
/*				hierstcmpcor[compmat[c1][c1]]+=sstrcmp[8];			/* hierarchical stratigraphic consistency	*/
/*				genrstcmpcor[compmat[c2][c2]]+=sstrcmp[5];			/* general stratigraphic consistency		*/
/*				strcstcmpcor[compmat[c2][c2]]+=sstrcmp[6];			/* strict stratigraphic consistency			*/
/*				hierstcmpcor[compmat[c2][c2]]+=sstrcmp[8];			/* hierarchical stratigraphic consistency	*/
				free_dvector(sstrcmp);
				}	/* end case of simulated compatible pair	*/
			}	/* end search of remaining simulated compatible pairs for possible compatibilities	*/
		}	/* search characters for possible compatibilities	*/
	
	/* get proportions of hierarchical stratigraphic characters relative to total								*/
	/* Also, get proportions of particular types of hierarchical stratigraphic characters						*/
	hierstratprp[r]+=(hierstratcmp[r]/genstratcomp[r]);						/* proportional hierarchical stratigraphic compatibility	*/
	/* if there is no hierarchical compatibility, then modify the routine: 2011-11-01	*/
	if (hierstratcmp[r]>0)	{
		hierstratana[r]/=hierstratcmp[r];										/* proportion of HSC matching anagenesis					*/
		hierstratbud[r]/=hierstratcmp[r];										/* proportion of HSC demanding budding						*/
		/* stick COG stuff in here???	*/
		}
	else	{
		++hierflop;																/* tally runs in which there was no hierarchical compatibility	*/
		hierstratana[r]=hierstratbud[r]=MAXRAND;
		}
	/* get average compatibilities of stratigraphically compatible and incompatible character pairs					*/
	avecompscp[r]/=genstratcomp[r];											/* added 2011-10-14	*/
	avecompsin[r]/=(cp-genstratcomp[r]);									/* added 2011-10-14	*/
	avecompdif[r]=avecompscp[r]-avecompsin[r];								/* added 2011-11-01: because I sort these later, this is otherwise lost */
	/* get proportion of SIN character pairs matching divergent or hierarchical patterns				*/
	compdivSIN[r]/=(cp-genstratcomp[r]);									/* compatibility-supported divergent SIN		*/
	comphierSIN[r]/=(cp-genstratcomp[r]);									/* compatibility-supported hierarchical SIN		*/
	stratdivSIN[r]/=(cp-genstratcomp[r]);									/* stratigraphy-supported divergent SIN			*/
	strathierSIN[r]/=(cp-genstratcomp[r]);									/* stratigraphy-supported hierarchical SIN		*/
	propend2cp[r]/=(cp-genstratcomp[r]);									/* relative diversity of younger SIN end pair	*/
	/* cp will vary a little from simulation to simulation because it is based on state numbers, so report things in proportions	*/
	genstratcomp[r]/=((double) cp);											/* general stratigraphic consistency						*/
	strstratcomp[r]/=((double) cp);											/* strict stratigraphic consistency							*/
	supstratcomp[r]/=((double) cp);											/* super strict stratigraphic consistency					*/
	hierstratcmp[r]/=((double) cp);											/* hierarchical stratigraphic consistency					*/
	
	/* 2013-04-14: Add COG, surivorship, disparity here!	*/
	aveCGH[r]/=((double) hiercounter);
/*	aveCGH[r]/=sC[r];					/* standardize first by opportunities and second by total clade shape	*/
	aveADH[r]/=((double) hiercounter);
	freqLFH[r]/=((double) hiercounter);
	aveCGD[r]/=((double) divcounter);
/*	aveCGD[r]/=sC[r];					/* standardize first by opportunities and second by total clade shape	*/
	aveADD[r]/=((double) divcounter);
	freqLFD[r]/=((double) divcounter);
	/* disparity stuff added 2014-01-14	*/
	PWDis=PhDfromData(simatrix,notu,nchars,UNKNOWN,INAP);
	sdispsummary=wholecladedisparitytaxa(PWDis, notu, ranges, onset, term);
	quart1[r]=sdispsummary[0];
	quart2[r]=sdispsummary[1];
	quart3[r]=sdispsummary[2];
	halftaxa[r]=sdispsummary[3];
	free_dmatrix(PWDis,notu,notu);
	free_dvector(sdispsummary);
	/* compare these to end of 

	
	/*	simmary[0]: P[observed time to qet 1/4 diversity]
		simmary[1]: E[observed time to qet 1/4 diversity]
		simmary[2]: P[observed time to qet 1/2 diversity]
		simmary[3]: E[observed time to qet 1/2 diversity]
		simmary[4]: P[observed time to qet 3/4 diversity]
		simmary[5]: E[observed time to qet 3/4 diversity]
		simmary[6]: P[observed disparity at S/2]
		simmary[7]: E[observed disparity at S/2]
		******************************/
	
	
	/* tally cases where simulated values are lower than observed values		*/
	if (genstratcomp[r]<summary[0])		simmary[0]+=(1/((double) RUNS));
	if (strstratcomp[r]<summary[1])		simmary[2]+=(1/((double) RUNS));
	if (supstratcomp[r]<summary[2])		simmary[4]+=(1/((double) RUNS));
	if (hierstratcmp[r]<summary[3])		simmary[6]+=(1/((double) RUNS));
	if (hierstratprp[r]<summary[4])		simmary[8]+=(1/((double) RUNS));
	if (hierstratana[r]<summary[5])		simmary[10]+=(1/((double) RUNS));
	if (hierstratbud[r]<summary[6])		simmary[12]+=(1/((double) RUNS));
	if (compdivSIN[r]<summary[7])		simmary[14]+=(1/((double) RUNS));
	if (comphierSIN[r]<summary[8])		simmary[16]+=(1/((double) RUNS));
	if (stratdivSIN[r]<summary[9])		simmary[18]+=(1/((double) RUNS));
	if (strathierSIN[r]<summary[10])	simmary[20]+=(1/((double) RUNS));
	if (propend2cp[r]<summary[11])		simmary[22]+=(1/((double) RUNS));
	if (avecompscp[r]<summary[12])		simmary[24]+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompsin[r]<summary[13])		simmary[26]+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompdif[r]<(summary[12]-summary[13]))
										simmary[27]+=(1/((double) RUNS));					/* added 2011-11-01: because I sort these later, this is otherwise lost */
	/* 2013-04-14: Add stuff COG & surivorship here!	*/
	if (aveCGH[r]<summary[14])			simmary[30]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (aveCGD[r]<summary[15])			simmary[32]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (CG[r]<summary[16])				simmary[34]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (aveADH[r]<summary[17])			simmary[36]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (aveADD[r]<summary[18])			simmary[38]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (freqLFH[r]<summary[19])			simmary[40]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (freqLFD[r]<summary[20])			simmary[42]+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (quart1[r]<summary[21])			simmary[44]+=(1/((double) RUNS));					/* added 2014-01-14 */
	if (quart2[r]<summary[22])			simmary[46]+=(1/((double) RUNS));					/* added 2014-01-14 */
	if (quart3[r]<summary[23])			simmary[48]+=(1/((double) RUNS));					/* added 2014-01-14 */
	if (halftaxa[r]<summary[24])		simmary[50]+=(1/((double) RUNS));					/* added 2014-01-14 */

	/* tally cases where simulated values are exactly equal to observed values	*/
	if (genstratcomp[r]==summary[0])	tiegs+=(1/((double) RUNS));
	if (strstratcomp[r]==summary[1])	tiess+=(1/((double) RUNS));
	if (supstratcomp[r]==summary[2])	tiexs+=(1/((double) RUNS));
	if (hierstratcmp[r]==summary[3])	tiehs+=(1/((double) RUNS));
	if (hierstratprp[r]==summary[4])	tiehp+=(1/((double) RUNS));
	if (hierstratana[r]==summary[5])	tieha+=(1/((double) RUNS));
	if (hierstratbud[r]==summary[6])	tiehb+=(1/((double) RUNS));
	if (compdivSIN[r]==summary[7])		tiesincd+=(1/((double) RUNS));
	if (comphierSIN[r]==summary[8])		tiesinch+=(1/((double) RUNS));
	if (stratdivSIN[r]==summary[9])		tiesinsd+=(1/((double) RUNS));
	if (strathierSIN[r]==summary[10])	tiesinsh+=(1/((double) RUNS));
	if (propend2cp[r]==summary[11])		tiesind2+=(1/((double) RUNS));
	if (avecompscp[r]==summary[12])		tieavecscp+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompsin[r]==summary[13])		tieavecsin+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompdif[r]==(summary[12]-summary[13]))
										tieavecdif+=(1/((double) RUNS));					/* added 2011-11-01: because I sort these later, this is otherwise lost */
	/* 2013-04-14: Add stuff COG & surivorship here!	*/
	if (aveCGH[r]==summary[14])			tieavecgh+=(1/((double) RUNS));						/* added 2013-04-15 */
	if (aveCGD[r]==summary[15])			tieavecgd+=(1/((double) RUNS));						/* added 2013-04-15 */
	if (CG[r]==summary[16])				tierc1c0+=(1/((double) RUNS));						/* added 2013-04-15 */
	if (aveADH[r]==summary[17])			tieaveadh+=(1/((double) RUNS));						/* added 2013-04-15 */
	if (aveADD[r]==summary[18])			tieaveadd+=(1/((double) RUNS));						/* added 2013-04-15 */
	if (freqLFH[r]==summary[19])		tiefreqlfh+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (freqLFD[r]==summary[20])		tiefreqlfd+=(1/((double) RUNS));					/* added 2013-04-15 */
	if (quart1[r]==summary[21])			tieq1+=(1/((double) RUNS));							/* added 2014-01-14 */
	if (quart2[r]==summary[22])			tieq2+=(1/((double) RUNS));							/* added 2014-01-14 */
	if (quart3[r]==summary[23])			tieq3+=(1/((double) RUNS));							/* added 2014-01-14 */
	if (halftaxa[r]==summary[24])		tiehalftaxa+=(1/((double) RUNS));					/* added 2014-01-14 */
	
	
	if ((r%10)==9)	{
//			printf("Doing Rate = %5.4f",mbl[2]);
		if (r>10)	{
			printf("\b\b\b\b\b");					/* clear ", Tree %d"		*/
			}
		if (r>18 && r<100) 			printf("\b\b\b");
		else if (r>99 && r<1000)	printf("\b\b\b\b");
		else if (r>999)				printf("\b\b\b\b\b");
		printf("Tree %d\n", r+1);
		}
	
/*	free_dvector(sC);						/* 2013-04-15		*/
	free_ulmatrix(compmat,nchars,nchars);
	free_lmatrix(vtree,clades+1,notu);
	free_lmatrix(tree,clades+3,notu);
	free_ivector(simauts);				/* found 2013-05-14	*/
	free_ulvector(charcomps);			/* found 2013-05-14	*/
	}

/* sort for medians	*/
dshellsort_inc_command(genstratcomp,RUNS); 		/* 0:  sort stratigraphic compatiblity						*/
dshellsort_inc_command(strstratcomp,RUNS); 		/* 1:  sort strict stratigraphic compatiblity				*/
dshellsort_inc_command(supstratcomp,RUNS);		/* 2:  sort super strict stratigraphic compatibility		*/
dshellsort_inc_command(hierstratcmp,RUNS);		/* 3:  sort hierarchical stratigraphic compatiblity (HSC)	*/
dshellsort_inc_command(hierstratprp,RUNS);		/* 4:  sort proportion HSC									*/
dshellsort_inc_command(hierstratana,RUNS);		/* 5:  sort HSC consistent with anagenesis					*/
dshellsort_inc_command(hierstratbud,RUNS);		/* 6:  sort HSC demanding budding							*/
dshellsort_inc_command(compdivSIN,RUNS);		/* 7:  sort compatibility suggests divergent SIN			*/
dshellsort_inc_command(comphierSIN,RUNS);		/* 8:  sort compatibility suggests hierarchical SIN			*/
dshellsort_inc_command(stratdivSIN,RUNS);		/* 9:  sort stratigraphy suggests divergent SIN				*/
dshellsort_inc_command(strathierSIN,RUNS);		/* 10: sort stratigraphy suggests hierarchical SIN			*/
dshellsort_inc_command(propend2cp,RUNS);		/* 11: sort relative diversity of derived end pair			*/
dshellsort_inc_command(avecompscp,RUNS);		/* 12: sort average compatibilities of SC character pairs	*/
dshellsort_inc_command(avecompsin,RUNS);		/* 13: sort average compatibilities of SIN character pairs	*/
dshellsort_inc_command(avecompdif,RUNS);		/* 14: sort average differences between SC & SIN pairs		*/
dshellsort_inc_command(aveCGH,RUNS);			/* 15: sort average 00 CoGs for hierarchical pairs		*/
dshellsort_inc_command(aveCGD,RUNS);			/* 16: sort average 00 CoGs for divergent pairs			*/
dshellsort_inc_command(aveADH,RUNS);			/* 17: sort average 00 durations for hierarchical pairs	*/
dshellsort_inc_command(aveADD,RUNS);			/* 18: sort average 00 durations for divergent pairs	*/
dshellsort_inc_command(CG,RUNS);				/* 19: sort average CoG for simulated clades				*/
dshellsort_inc_command(freqLFH,RUNS);			/* 20: sort hierarchical living fossils						*/
dshellsort_inc_command(freqLFD,RUNS);			/* 21: sort hierarchical living fossils						*/
dshellsort_inc_command(quart1,RUNS);
dshellsort_inc_command(quart2,RUNS);
dshellsort_inc_command(quart3,RUNS);
dshellsort_inc_command(halftaxa,RUNS);

/* get alpha values: if observed values are "low" then we want p[observed or lower]: so, add ties to total	*/
/* for "high" values, we do not need to do this: 1-simmary gives p[observed or more extreme] as we tallied only sims less than observed	*/
/* we want the probability of the observed OR MORE EXTREME; because we tallied every time the observed was lower than what we saw,
		we don't need to do this when p> 0.5	*/
if (simmary[0]<0.5)		simmary[0]+=tiegs;			/* general stratigraphic compatibility	*/
if (simmary[2]<0.5)		simmary[2]+=tiess;			/* strict stratigraphic compatibility	*/
if (simmary[4]<0.5)		simmary[4]+=tiexs;
if (simmary[6]<0.5)		simmary[6]+=tiehs;
if (simmary[8]<0.5)		simmary[8]+=tiehp;
if (simmary[10]<0.5)	simmary[10]+=tieha;
if (simmary[12]<0.5)	simmary[12]+=tiehb;
if (simmary[14]<0.5)	simmary[14]+=tiesincd;
if (simmary[16]<0.5)	simmary[16]+=tiesinch;
if (simmary[18]<0.5)	simmary[18]+=tiesinsd;
if (simmary[20]<0.5)	simmary[20]+=tiesinsh;
if (simmary[22]<0.5)	simmary[22]+=tiesind2;
if (simmary[24]<0.5)	simmary[24]+=tieavecscp;
if (simmary[26]<0.5)	simmary[26]+=tieavecsin;
if (simmary[28]<0.5)	simmary[28]+=tieavecsin;
if (simmary[30]<0.5)	simmary[30]+=tieavecgh;
if (simmary[32]<0.5)	simmary[32]+=tieavecgd;
if (simmary[34]<0.5)	simmary[34]+=tieaveadh;
if (simmary[36]<0.5)	simmary[36]+=tieaveadd;
if (simmary[38]<0.5)	simmary[38]+=tierc1c0;
if (simmary[40]<0.5)	simmary[40]+=tiefreqlfh;
if (simmary[42]<0.5)	simmary[42]+=tiefreqlfd;
if (simmary[44]<0.5)	simmary[44]+=tieq1;
if (simmary[46]<0.5)	simmary[46]+=tieq2;
if (simmary[48]<0.5)	simmary[48]+=tieq3;
if (simmary[50]<0.5)	simmary[50]+=tiehalftaxa;

/* adjust for runs with no hierarchical compatibility	*/
simmary[10]*=(((double) (RUNS-hierflop))/((double) RUNS));		/* adjust significance of anagenetic hierarchical compatibility	*/
simmary[12]*=(((double) (RUNS-hierflop))/((double) RUNS));		/* adjust significance of budding hierarchical compatibility	*/


/* tally median expected stratigraphic compatibilities	*/
simmary[1]=genstratcomp[RUNS/2];
simmary[3]=strstratcomp[RUNS/2];
simmary[5]=supstratcomp[RUNS/2];
simmary[7]=hierstratcmp[RUNS/2];
simmary[9]=hierstratprp[RUNS/2];
simmary[11]=hierstratana[(RUNS-hierflop)/2];			/* accommodate runs with no hierarchical compatibility	*/
simmary[13]=hierstratbud[(RUNS-hierflop)/2];			/* accommodate runs with no hierarchical compatibility	*/
simmary[15]=compdivSIN[RUNS/2];
simmary[17]=comphierSIN[RUNS/2];
simmary[19]=stratdivSIN[RUNS/2];
simmary[21]=strathierSIN[RUNS/2];
simmary[23]=propend2cp[RUNS/2];
simmary[25]=avecompscp[RUNS/2];
simmary[27]=avecompsin[RUNS/2];
simmary[29]=avecompdif[RUNS/2];
simmary[31]=aveCGH[RUNS/2];
simmary[33]=aveCGD[RUNS/2];
simmary[35]=CG[RUNS/2];
simmary[37]=aveADH[RUNS/2];
simmary[39]=aveADD[RUNS/2];
simmary[41]=freqLFH[RUNS/2];
simmary[43]=freqLFD[RUNS/2];
simmary[45]=quart1[RUNS/2];
simmary[47]=quart2[RUNS/2];
simmary[49]=quart3[RUNS/2];
simmary[51]=halftaxa[RUNS/2];

/* return output to here!	*/
/* TEMPORARY!  UNTIL I'VE DEBUGGED THIS.....	*/
keep=1;
//keep=0;
//printf("Enter '1' if you want output for this taxon: ");
//scanf("%i",&keep);	

if (keep==1)	{
	strcpy(outfile,taxonname);
	strcat(outfile,"_");
	strcat(outfile,citation);
	strcat(outfile,"_Expectations");
	if (mbl[4]==0)	strcat(outfile,"_Bud.xls");
	else			strcat(outfile,"_Bif.xls");
	
	output=fopen(outfile,"w");
	fprintf(output,"Run");
	fprintf(output,"\tGeneral Stratigraphic Compatibility (SC)");
	fprintf(output,"\tStrict SC");
	fprintf(output,"\tSuper Strict SC");
	fprintf(output,"\tHierarchical SC");
	fprintf(output,"\tProportional HSC");
	fprintf(output,"\tAnagenetic HSC");
	fprintf(output,"\tBudding HSC");
	fprintf(output,"\tDivergent SIN (Comp)");
	fprintf(output,"\tHierarchical SIN (Comp)");
	fprintf(output,"\tDivergent SIN (Strat)");
	fprintf(output,"\tHierarchical SIN (Strat)");
	fprintf(output,"\tProp Diversity Younger SIN Pair");
	fprintf(output,"\tµ Character Compabibility (SC)");					/* added 2011-10-14	*/
	fprintf(output,"\tµ Character Compabibility (SIN)");				/* added 2011-10-14	*/
	fprintf(output,"\tµ Compabibility Diff (SC-SIN)");					/* added 2011-11-01	*/
	fprintf(output,"\tµ H00 CoG");										/* added 2013-04-16	*/
	fprintf(output,"\tµ D00 CoG");										/* added 2013-04-16	*/
	fprintf(output,"\tµ Clade CoG");									/* added 2013-04-16	*/
	
	for (r=0; r<RUNS; ++r)	{
		fprintf(output,"%d",r+1);
		fprintf(output,"\t%4.3f",genstratcomp[r]);
		fprintf(output,"\t%4.3f",strstratcomp[r]);
		fprintf(output,"\t%4.3f",supstratcomp[r]);
		fprintf(output,"\t%4.3f",hierstratcmp[r]);
		fprintf(output,"\t%4.3f",hierstratprp[r]);
		if (hierstratana[r]<MAXRAND)	fprintf(output,"\t%4.3f",hierstratana[r]);
		else							fprintf(output,"\t•");
		if (hierstratbud[r]<MAXRAND)	fprintf(output,"\t%4.3f",hierstratbud[r]);
		else							fprintf(output,"\t•");
		fprintf(output,"\t%4.3f",compdivSIN[r]);
		fprintf(output,"\t%4.3f",comphierSIN[r]);
		fprintf(output,"\t%4.3f",stratdivSIN[r]);
		fprintf(output,"\t%4.3f",strathierSIN[r]);
		fprintf(output,"\t%4.3f",propend2cp[r]);
		fprintf(output,"\t%4.3f",avecompscp[r]);					/* added 2011-10-14	*/
		fprintf(output,"\t%4.3f",avecompsin[r]);					/* added 2011-10-14	*/
		fprintf(output,"\t%4.3f",avecompdif[r]);					/* added 2011-11-01	*/
		fprintf(output,"\t%4.3f",aveCGH[r]);						/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f",aveCGD[r]);						/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f",CG[r]);							/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f",aveADH[r]);						/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f",aveADD[r]);						/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f",aveADD[r]);						/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f",freqLFH[r]);						/* added 2013-04-16	*/
		fprintf(output,"\t%4.3f\n",freqLFD[r]);						/* added 2013-04-16	*/
		}
	fclose(output);
	}


free_dvector(aveADD);				/* tallies for simmary[32,33]							*/	
free_dvector(aveADH);				/* tallies for simmary[30,31]							*/	
free_dvector(aveCGD);				/* tallies for simmary[32,33]							*/	
free_dvector(aveCGH);				/* tallies for simmary[30,31]							*/	
free_dvector(avecompdif);			/* tallies for simmary[28,29]							*/		
free_dvector(avecompscp);			/* tallies for simmary[24,25]							*/		
free_dvector(avecompsin);			/* tallies for simmary[26.27]							*/		
free_dvector(CG);					/* tallies for simmary[38,39]							*/
free_dvector(compdivSIN);			/* tallies sstrcomp[11], for simmary[14,15]				*/					
free_dvector(comphierSIN);			/* tallies sstrcomp[12], for simmary[16,17]				*/					
free_dvector(freqLFD);				/* tallies for simmary[36,37]							*/	
free_dvector(freqLFH);				/* tallies for simmary[34,35]							*/	
free_dvector(genstratcomp);			/* tallies sstrcomp[5], for simmary[0,1]				*/					
free_dvector(hierstratana);			/* tallies sstrcomp[9], for simmary[10,11]				*/					
free_dvector(hierstratbud);			/* tallies sstrcomp[10], for simmary[12,13]				*/					
free_dvector(hierstratcmp);			/* tallies sstrcomp[8], for simmary[4,5]				*/					
free_dvector(hierstratprp);			/* tallies sstrcomp[8]/sstrcomp[5], for simmary[6,7]	*/								
free_dvector(propend2cp);			/* tallies sstrcomp[15], for simmary[22,23]				*/					
free_dvector(stratdivSIN);			/* tallies sstrcomp[13], for simmary[18,19]				*/					
free_dvector(strathierSIN);			/* tallies sstrcomp[14], for simmary[20,21]				*/					
free_dvector(strstratcomp);			/* tallies sstrcomp[6], for simmary[2,3]				*/					
free_dvector(supstratcomp);			/* tallies sstrcomp[7], for simmary[8,9]				*/					
free_dvector(quart1);
free_dvector(quart2);
free_dvector(quart3);
free_dvector(halftaxa);
free_lmatrix(ranges,notu,2);		/* simulated ranges										*/
free_lmatrix(simatrix,notu,nchars);	/* simulated morphologies								*/

return simmary;
}


/** dateclade - Finds age of clades based on the oldest member in the clade
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
long **dateallbranchesfull(long **ranges, long **tree, long *bl, int clades, int notu)
{
int	a, b, anc, sp, htu, node;
long **ofa;

a=notu+clades;
ofa=lmatrix(a,2);			/* ofa[0]: inferred origin; ofa[1]: "sighting"	*/

for (sp=0; sp<notu; ++sp)				ofa[sp][1]=ranges[sp][0];		/* note sighting of species	*/

/* date each clade	*/
for (node=clades-1; node>=0; --node)	{
	htu=node+notu;
	ofa[htu][1]=ofa[a=tree[node][1]][1];	/* first appearance of taxon within node is the "sighting" of the node	*/
	for (b=2; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (ofa[htu][1]>ofa[sp][1])	ofa[htu][1]=ofa[sp][1];
		}	/* end finding "sightings" of each clade	*/
/*	a=tree[node][1];				/* possible ancestor	*/
	anc=-1;
	if (a<notu && bl[a]==0)	{
		anc=a;
		ofa[anc][0]=ofa[htu][1];	/* inferred origin of ancestor must match that of clade (usually will anyway)	*/
		}

	for (b=2; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (anc==-1)	ofa[sp][0]=ofa[htu][1];
		else			ofa[sp][0]=lmin(ofa[sp][1],(1+ranges[anc][1]));
		}
	}	
return ofa;
}


/* routine to compare any given character pair and assess how many comparisons can be made (0 for incompatible, 1 for binary compatible with 3 pairs, possible more for multistate)
 	and how many are generall stratigraphically compatible		*/
int	*genstratcompat(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int	a, b, c, d, p;
int	maxocc, st1, st2, ttlprs=0, maxst;
//int	cc, sp, maxocc, st1, st2, sc1, sc2, rc, scf1, scf2, scu1, scu2, prs, m1, m2, ttlprs=0, maxst;
//int cl[20], cr[20];
int	swing, end1, end2;									/* these are the "swing" states between characters + important info, with "swing" giving pair #	*/
/*int	kp1, kp2;										/* key pairs			*/
/*int	prdv1, prdv2;									/* key pair diversities	*/
int cp1, cp2, cp3;
unsigned long **combos;
/*long **pairfnd;*/
long **pairs, **allpairs, **pairrng, **pairfnd;
long prfl[3][2];
int	*stratcompat;
int	comparisons=0, genstratcomp=1;

stratcompat=ivector(2);
//long pairs[3][2];

/* 2011-10-19: NEW APPROACH!  
	1. Put all pairs into a pair x 2 matrix
	2. Put in a 3-tiered for loop to find all combinations of 3 pairs;
	3. If there is a swing-pair, then analyze this pair.
****************************************************/
maxocc=maxlmatrixcol(ranges, notu, 1);
maxst=imax(states[ch1],states[ch2]);
pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));

for (st1=0; st1<states[ch1]; ++st1)	{
	for (st2=0; st2<states[ch2]; ++st2)	{
		a=st1*states[ch2]+st2;
		/* find the number of combinations in which each states of both characters appears	*/
		if (pairrng[a][0]<=maxocc && pairrng[a][0]>=0)	{
			++combos[0][st1];
			++combos[1][st2];
			++ttlprs;			/* total pairs observed: if this is less than 3, then we can quit now!	*/
			}
		}
	}
free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));

if (ttlprs>2)	{
	allpairs=lmatrix(ttlprs,2);
	p=0;
	for (st1=0; st1<states[ch1]; ++st1)	{
		for (st2=0; st2<states[ch2]; ++st2)	{
			a=st1*states[ch2]+st2;
			/* find the number of combinations in which each states of both characters appears	*/
			if (pairrng[a][0]<=maxocc && pairrng[a][0]>=0)	{
				allpairs[p][0]=st1;
				allpairs[p][1]=st2;
				++p;			/* total pairs observed: if this is less than 3, then we can quit now!	*/
				}
			}
		}

	/* make sure that pairfnd uses (st1*states[ch2])+states[ch2] throughout, and not simply the pair number!	*/
	pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/

	/* now, go through all combinations of three pairs, and exampine those that have a linking pair	*/
	pairs=lmatrix(3,2);
	for (a=0; a<(ttlprs-2); ++a)	{
		pairs[0][0]=allpairs[a][0];
		pairs[0][1]=allpairs[a][1];
		/* get first-last appearances of pair a	*/
		prfl[0][0]=pairrng[cp1=(pairs[0][0]*states[ch2])+pairs[0][1]][0];
		prfl[0][1]=pairrng[cp1][1];
		/* determine if this pair has a gap in it	*/
		for (b=a+1; b<(ttlprs-1); ++b)	{
			pairs[1][0]=allpairs[b][0];
			pairs[1][1]=allpairs[b][1];
			/* get first-last appearances of pair b	*/
			prfl[1][0]=pairrng[cp2=(pairs[1][0]*states[ch2])+pairs[1][1]][0];
			prfl[1][1]=pairrng[cp2][1];
			/* determine if this pair has a gap in it	*/
			for (c=b+1; c<ttlprs; ++c)	{
				pairs[2][0]=allpairs[c][0];
				pairs[2][1]=allpairs[c][1];
				/* get first-last appearances of pair c	*/
				prfl[2][0]=pairrng[cp3=(pairs[2][0]*states[ch2])+pairs[2][1]][0];
				prfl[2][1]=pairrng[cp3][1];
				/* determine if this pair has a gap in it	*/
				/* determine whether this set is linked: if so, then swing will be positive	*/
				swing=findswingpair(pairs,3,maxst);
				if (swing>-1)	{
					stratcompat[comparisons]+=1.0f;							/* another comparison made						*/
					end1=0;
					end2=2;
					if (swing==0)
						end1=1;
					else if (swing==2)
						end2=1;
					
					/* make end1 older than end2	*/
					if (prfl[end1][0]>prfl[end2][0])	{
						d=end1;
						end1=end2;
						end2=d;
						}
					/* if swing pair appears after both end pairs appear, then this is stratigraphically incompatible	*/
					if ((prfl[swing][0]<=prfl[end1][0]) || (prfl[swing][0]<=prfl[end2][0]))	{
						stratcompat[genstratcomp]+=1.0f;							/* general stratigraphic compatibility						*/
						}	/* end case of stratigraphic compatibility	*/
					}	/* end case of compatible trio of character pairs	*/
				}	/* end possible 3rd pairs	*/
			}	/* end possible 2nd pairs	*/
		}	/* end possible 1st pairs	*/
	
	free_lmatrix(pairs,3,2);
	free_lmatrix(allpairs,ttlprs,2);
	free_lmatrix(pairfnd,states[ch1]*states[ch2],1+maxocc);
	}	/* end cases of 2+ pairs of character states	*/

free_lmatrix(pairrng,states[ch1]*states[ch2],2);

return stratcompat;
}


/* routine to compare any given character pair and assess how many comparisons can be made (0 for incompatible, 1 for binary compatible with 3 pairs, possible more for multistate)
 	and how many are generall stratigraphically compatible		*/
double	propgenstratcompat(long **ranges, long **matrix, int *states, unsigned long **compmatrix, int notu, int nchars, int UNKNOWN, int INAP)
{
int	a, b, c, d, p;
int ch1, ch2;
int	maxocc, st1, st2, ttlprs=0, maxst;
//int	cc, sp, maxocc, st1, st2, sc1, sc2, rc, scf1, scf2, scu1, scu2, prs, m1, m2, ttlprs=0, maxst;
//int cl[20], cr[20];
int	swing, end1, end2;									/* these are the "swing" states between characters + important info, with "swing" giving pair #	*/
/*int	kp1, kp2;										/* key pairs			*/
/*int	prdv1, prdv2;									/* key pair diversities	*/
int cp1, cp2, cp3;
//unsigned long **combos;
/*long **pairfnd;*/
long **pairrng, **pairs;//, **allpairs,, **pairfnd;
long prfl[3][2], allpairs[100][2];//, pairs[3][2];
double	comparisons=0.0f, genstratcomp=0.0f;

/* 2011-10-19: NEW APPROACH!  
	1. Put all pairs into a pair x 2 matrix
	2. Put in a 3-tiered for loop to find all combinations of 3 pairs;
	3. If there is a swing-pair, then analyze this pair.
****************************************************/
maxocc=maxlmatrixcol(ranges, notu, 1);
pairs=lmatrix(3,2);

for (ch1=0; ch1<(nchars-1); ++ch1)	{
	for (ch2=(ch1+1); ch2<nchars; ++ch2)	{
		if (compmatrix[ch1][ch2]==1)	{
			ttlprs=0;							/* this must be reset at each comparison	*/
			maxst=imax(states[ch1],states[ch2]);
			pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
//			combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));

			for (st1=0; st1<states[ch1]; ++st1)	{
				for (st2=0; st2<states[ch2]; ++st2)	{
					a=st1*states[ch2]+st2;
					/* find the number of combinations in which each states of both characters appears	*/
					if (pairrng[a][0]<=maxocc && pairrng[a][0]>=0)	{
//						++combos[0][st1];
//						++combos[1][st2];
						++ttlprs;			/* total pairs observed: if this is less than 3, then we can quit now!	*/
						}
					}	/* end states from character 2	*/
				}	/* end states from character 1 */
//			free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));

			if (ttlprs>2)	{
//				allpairs=lmatrix(ttlprs,2);
				for (a=0; a<ttlprs; ++a)	for (b=0; b<2; ++b)	allpairs[a][b]=0;
				p=0;
				for (st1=0; st1<states[ch1]; ++st1)	{
					for (st2=0; st2<states[ch2]; ++st2)	{
						a=st1*states[ch2]+st2;	/* this puts a number on the particular combination	*/
						/* find the number of combinations in which each states of both characters appears	*/
						if (pairrng[a][0]<=maxocc && pairrng[a][0]>=0)	{
							allpairs[p][0]=st1;
							allpairs[p][1]=st2;
							++p;			/* total pairs observed: if this is less than 3, then we can quit now!	*/
							}	/*	allpairs[x][0] gives the state of character 1, allpairs[x][1] gives the state of character 2	*/
						}
					}

				/* make sure that pairfnd uses (st1*states[ch2])+states[ch2] throughout, and not simply the pair number!	*/
				/* pairfnd no longer seems to be necessary	*/
//				pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/

				/* now, go through all combinations of three pairs, and exampine those that have a linking pair	*/				
				for (a=0; a<(ttlprs-2); ++a)	{
					pairs[0][0]=allpairs[a][0];
					pairs[0][1]=allpairs[a][1];
					/* get first-last appearances of pair a	*/
					prfl[0][0]=pairrng[cp1=(pairs[0][0]*states[ch2])+pairs[0][1]][0];
					prfl[0][1]=pairrng[cp1][1];
					/* determine if this pair has a gap in it	*/
					for (b=a+1; b<(ttlprs-1); ++b)	{
						pairs[1][0]=allpairs[b][0];
						pairs[1][1]=allpairs[b][1];
						/* get first-last appearances of pair b	*/
						prfl[1][0]=pairrng[cp2=(pairs[1][0]*states[ch2])+pairs[1][1]][0];
						prfl[1][1]=pairrng[cp2][1];
						/* determine if this pair has a gap in it	*/
						for (c=b+1; c<ttlprs; ++c)	{
							pairs[2][0]=allpairs[c][0];
							pairs[2][1]=allpairs[c][1];
							/* get first-last appearances of pair c	*/
							prfl[2][0]=pairrng[cp3=(pairs[2][0]*states[ch2])+pairs[2][1]][0];
							prfl[2][1]=pairrng[cp3][1];
							/* determine if this pair has a gap in it	*/
							/* determine whether this set is linked: if so, then swing will be positive	*/
							swing=findswingpair(pairs,3,maxst);
							if (swing>-1)	{
								comparisons+=1.0f;							/* another comparison made	*/
								end1=0;
								end2=2;
								if (swing==0)
									end1=1;
								else if (swing==2)
									end2=1;
								
								/* make end1 older than end2	*/
								if (prfl[end1][0]>prfl[end2][0])	{
									d=end1;
									end1=end2;
									end2=d;
									}
								/* if swing pair appears after both end pairs appear, then this is stratigraphically incompatible	*/
								if ((prfl[swing][0]<=prfl[end1][0]) || (prfl[swing][0]<=prfl[end2][0]))	{
									genstratcomp+=1.0f;							/* general stratigraphic compatibility						*/
									}	/* end case of stratigraphic compatibility	*/
								}	/* end case of compatible trio of character pairs	*/
							}	/* end possible 3rd pairs	*/
						}	/* end possible 2nd pairs	*/
					}	/* end possible 1st pairs	*/
				}	/* end cases of 2+ pairs of character states	*/
			/* pairs=1+(states[ch1]*states[ch2])+states[ch2];
				pairrng=lmatrix(pairs,2);	*/
			free_lmatrix(pairrng,(1+(states[ch1]*states[ch2])+states[ch2]),2);			/* make sure that this deallocation is correct	*/
/*			free_lmatrix(pairrng,states[ch1]*states[ch2],2);			/* make sure that this deallocation is correct	*/
			}	/* end case of compatible character pair	*/
		}	/* cend comparisons with other characters for ch1	*/
	}	/* end examination of character pairs	*/

//free_lmatrix(allpairs,ttlprs,2);
//free_lmatrix(pairfnd,states[ch1]*states[ch2],1+maxocc);
free_lmatrix(pairs,3,2);
return (genstratcomp/comparisons);
}