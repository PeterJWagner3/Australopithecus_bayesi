#define compatibility_functions
#include <time.h>
#include "compatibility_functions.h"
#include "matrixanalysis.h"
#include "matrixchange.h"
#include "matrixreading.h"
#include "memory.h"
#include "minmax.h"
#include "MonteCarloPhylogenyFunctions.h"
#include "tree_read.h"
#define COMPLIKE		"_L[Steps|C,D]"


/* charcombs - counts the number of combinations between two characters, ch1 & ch2
Requires:
	matrix: notu x nchars array of character states
	states: array giving number of states for character x
	ch1: the 1st character
	ch2: the 2nd character
	notu: number of taxa
	UNKNOWN: value for "?"
	INAP: value for inapplicable
Returns:
	incompatible: 0 for compatible, 1 for incompatible
**********************************************************************************************************************************************************/
int charcombs(long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int a, b, c, cc=0, sp=0, st1=0;
long **combos;

combos=lmatrix(2,1+(a=imax(states[ch1],states[ch2])));
clearlmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])),-1);
for (st1<=0; st1<states[ch1]; ++st1)	{
	for (sp=0; sp<notu; ++sp)	{
		if (matrix[sp][ch1]==st1)	{
			b=0;
			for (c=0; c<cc && b==0; ++c)	{
				if (combos[st1][c]==matrix[sp][ch2] || (matrix[sp][ch2]!=UNKNOWN && matrix[sp][ch2]!=INAP))
					b=1;										/* combo already found	*/
				else if (combos[st1][c]==-1)	{
					combos[st1][c]=matrix[sp][ch2];
					+cc;										/* new combo found		*/
					b=1;										/* end loop				*/
					}
				}
			}
		}
	}
free_lmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));
}
/* charcombos - calculates the number of combinations
Requires:
	t1: sorted character array, with character having the greatest number of states
	t2: second character array sorted by character one
	st1: number of states for the first character
	st2: number of states for the first character
	notu: number of taxa
	UNKNOWN: value for "?"
	INAP: value for inapplicable
Returns:
	incompatible: 0 for compatible, 1 for incompatible
**********************************************************************************************************************************************************/
int *charcombos(long *t1, long *t2, int st, int notu, int UNKNOWN, int INAP)
{
int sp=0, s=0;
int *combos;

combos=ivector(st);

sp=0;
for (s=0; s<st; ++s)	{
	if (t2[sp]!=UNKNOWN && t2[sp]!=INAP)	combos[s]=1;
	for (sp=sp; t1[sp]==s && sp<notu; ++sp)	{
		while (sp<notu && (t2[sp+1]==t2[sp] && t1[sp+1]==t1[sp]))	++sp;
		if (t1[sp+1]==s && (t2[sp+1]!=UNKNOWN && t2[sp+1]!=INAP))	++combos[s];
		}
	}
return combos;
}
/* unordcompatible - determines whether a character pair involving an unordered multistate are compatible
Requires:
	t1: sorted character array, with character having the greatest number of states
	t2: second character array sorted by character one
	st1: number of states for the first character
	st2: number of states for the first character
	notu: number of taxa
	UNKNOWN: value for "?"
	INAP: value for inapplicable
Returns:
	incompatible: 0 for compatible, 1 for incompatible
**********************************************************************************************************************************************************/
int unordcompatible(int *t1, int *t2, int st1, int st2, int notu, int UNKNOWN, int INAP)
{
int	a, b, c, d, f, s, sp, sp2, maxst;
int	total, incompatible=0;
int link[20][20], shared[20][20], links[20];
int *combos;

maxst=st1;
if (st2>maxst)	maxst=st2;
combos=charcombos(((long *) t1),((long *) t2), st1, notu, UNKNOWN, INAP);

for (a=0; a<st1-1; ++a)	{
	while (combos[a]<=1 && a<st1-1)	++a;
	if (a>=st1-1)	break;
	for (b=a+1; b<st1; ++b)	{
		while (combos[b]<=1 && b<st1)	++b;
		if (b>=st1)	break;
		for (c=0; c<st2-1; ++c)	{
			for (d=c+1; d<st2; ++d)	{
				total=0;
				for (sp=0; sp<notu && total<4; ++sp)	{
					if (t1[sp]==a)	{
						if (t2[sp]==c)	{
							++total;		/* note that this requires t1 to be sorted */
							while (t2[sp+1]==c && sp<notu)	++sp;	/* get through all ac taxa */
							}	/* end search for ac taxa */
						else if (t2[sp]==d)	{
							++total;
							while (t2[sp+1]==d && sp<notu)	++sp;	/* get through all ad taxa */
							}	/* end search for ad taxa */
						/* if we have passed d, then there is no point at looking at these taxa anymore */
						else if (t2[sp]>d || (t2[sp]==UNKNOWN || t2[sp]==INAP))	
							while (t1[sp+1]==a)	++sp;				/* get through a¥ taxa */
						else
							while (t2[sp+1]<c || (t2[sp+1]>c && t2[sp+1]<d))	++sp;	/* skip states between c and d */
						}	/* end search for a¥ taxa */
					else if (t1[sp]==b)	{
						if (t2[sp]==c)	{
							++total;
							while (t2[sp+1]==c && sp<notu)	++sp;	/* get through all bc taxa */
							}	/* end search for bc taxa */
						else if (t2[sp]==d)	{
							++total;
							while (t2[sp+1]==d && sp<notu)	++sp;	/* get through all bd taxa */
							}	/* end search for bd taxa */
						/* if we have passed d, then there is no point at looking at these taxa anymore */
						else if (t2[sp]>d || (t2[sp]==UNKNOWN || t2[sp]==INAP))	sp=notu;
						else
							while (t2[sp+1]<c || (t2[sp+1]>c && t2[sp+1]<d))	++sp;	/* skip states between c and d */
						}	/* end search for b¥ taxa */
					else if (t1[sp]>b || (t1[sp]==UNKNOWN ||t1[sp]==INAP))	sp=notu;	/* skip ahead if past b */
					else
						while (t1[sp+1]<a || (t1[sp+1]>a && t1[sp+1]<b))	++sp;		/* skip states between a and b */
					}	/* end examination of taxa */
				if (total>=4)	{
					incompatible=1;
					a=b=c=d=maxst;		/* sets all loop variables to end of loop or beyond, ending routine */
					}	/* if incompatible, then there is no point in searching further */
				}	/* end examination of combinations for ¥d */
			}	/* end examination of combinations for ¥c */
		}	/* end examination of combinations for b¥ */
	}	/* end examination of combinations for a¥ */

if (st2>2 && incompatible==0)	{

	for (a=0; a<20; ++a)	{
		links[a]=0;
		for (b=0; b<20; ++b)	shared[a][b]=link[a][b]=-1;
		}
	a=0;
	for (sp=0; sp<notu-1; ++sp)	{
		for (sp2=sp+1; sp2<notu; ++sp2)	{
			if (t2[sp2]==t2[sp] && t1[sp2]!=t1[sp])	{
				link[t1[sp]][links[a]]=t1[sp2];
				++links[a];
				shared[t1[sp]][t1[sp2]]=shared[t1[sp2]][t1[sp]]=t2[sp];
				while ((t2[sp2+1]==t2[sp2] && t1[sp2+1]==t1[sp2]) && sp2<notu-1)	++sp2;
				}
			}
		while ((t2[sp+1]==t2[sp] && t1[sp+1]==t1[sp]) && sp<notu-1)	++sp;
		if (t1[sp]!=t1[sp+1])	++a;
		}
		
	/* now, make sure that there are no "triangles" among compatible pairs; for example, if state A is linked to states B and C, then B and C cannot be linked	*/ 
	/* example:		00
					02
					11
					12
					20
					21	*/
	/* first, go through all of character 1's states that are paired with the same states of character 2; e.g., char 1 state 0 above is paired with char 1 state 1
		through char 2 state 2 and char 1 state 2 through char 2 state 0; char 1 state 1 is paired with car 1 state 2 through char 2 state 1					*/
	for (a=0; a<st1-2 && incompatible==0; ++a)	{
		for (b=0; b<links[a]-1 && incompatible==0; ++b)	{
			s=link[a][b];	/* s is a state that is linked to char 1 state a by char 2 state X	*/
			/* sometimes missing or inapplicable data will result in one state being omitted	*/
			/* 		this will skip past that to the next state									*/
			while (s==-1 && b<links[a]-1)	{
				++b;
				s=link[a][b];
				}
			if (b>=links[a]-1)	break;
			/* now search through states to which car 1 state a is linked						*/
			for (c=0; c<links[s] && incompatible==0; ++c)	{
				for (d=b+1; d<links[a] && incompatible==0; ++d)	{
					if ((f=link[a][d])==link[s][c])	{
						/* both a & s are linked to a third character	*/
						/* if linked by different states, then the characters are incompatible	*/
						if (shared[a][s]!=shared[s][f])	incompatible=1;
						}
					}
				}	/* search states linked to state s and make sure none also are linked to state a	*/ 
			}
		}	/* finish making sure that no triangles exist	*/
	}

free_ivector(combos);
return incompatible;
}

/* ordcompatible - determines whether a character pair involving an ordered multistate are compatible
Requires:
	t1: sorted character array, with character having the greatest number of states
	t2: second character array sorted by character one
	notu: number of taxa
	UNKNOWN: value for "?"
	INAP: value for inapplicable
Returns:
	incompatible: 0 for compatible, 1 for incompatible
**********************************************************************************************************************************************************/
int ordcompatible(int *t1, int *t2, int st1, int st2, int notu, int UNKNOWN, int INAP)
{
int	drop, rise, sp, rev;
int incompatible=0;

drop=rise=0;

incompatible=unordcompatible(t1, t2, st1, st2, notu, UNKNOWN, INAP);

if (incompatible==0)	{
	/* make sure that there are coded species	*/
	drop=rise=rev=0;

	for (sp=1; sp<notu && rev<2; ++sp)	{
		if (t2[sp]>t2[sp-1])	{
			if (drop==1)	++rev;
			rise=1;
			}
		else if (t2[sp]<t2[sp-1])	{
			if (rise==1)	++rev;
			drop=1;
			}
		
		if (rev==2)	incompatible=1;
		}		/* work this out based on drops and rises; states cannot zig-zag when ordered	*/
	}
return incompatible;
}	/* end routine for ordered multistates	*/

/* binarycompatible - determines whether a pair of binary characters are compatible
Requires:
	t1: sorted character array, with character having the greatest number of states
	t2: second character array sorted by character one
	notu: number of taxa
	UNKNOWN: value for "?"
	INAP: value for inapplicable
Returns:
	incompatible: 0 for compatible, 1 for incompatible
**********************************************************************************************************************************************************/
int binarycompatible (int *t1, int *t2, int notu, int UNKNOWN, int INAP)
{

int sp, ttl=0, incompatible=0;

for (sp=1; sp<notu && ttl<4; ++sp)	{
	if (t1[sp]==UNKNOWN || t1[sp]==INAP)			sp = notu;
	else if (t2[sp]>t2[sp-1] && (t2[sp]!=UNKNOWN && t2[sp]!=INAP))	++ttl;
	else if (t1[sp]>t1[sp-1])						++ttl;
	if (ttl>=3)	{
		sp = notu;
		incompatible=1;
		}	
	}
return incompatible;
}

/* Function returning compability matrix.  O means incompatible, 1 means compatible.
Neads:
	states: #states per character
	notu: number of taxa
	chmatrix: character matrix
	type: character type (0 = ordered, 1 = unordered)
	nchars: number of characters
	comptype: 0: general compatibility; 1: hierarchical compatibilty
	OUTGROUP: number of outgroup taxon
	UNKNOWN: value for "?"
	INAP: value for inapplicable
Returns:
	comatrix: compatibility matrix
*****************************************************************************/
unsigned long **compatible(int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP)
{
int a, b, d, ch1, ch2, sp, placer, incompatible, inhiercompat, mxst, hi, lo, compsp, maxst;
int	mn1, mn2, mx2, rev;
int *autap, *character1, *character2, *mnst, *mtch, *tallied, *taxst, *test1, *test2;
unsigned long	**comatrix;

a=outgroup;		/* delete this if you do not restore hierarchical compatibility	*/

comatrix=ulmatrix(nchars,nchars);
mnst=ivector(nchars);
character1=ivector(notu);
character2=ivector(notu);
test1=ivector(notu);
test2=ivector(notu);

maxst=maxiarray(nstates,nchars);
for (a=0; a<nchars; ++a)	{
	mnst[a]=100;
	for (b=0; b<notu; ++b)	{
		if (chmatrix[b][a]<mnst[a] && (chmatrix[b][a]!=UNKNOWN && chmatrix[b][a]!=INAP))
			mnst[a]=chmatrix[b][a];
		}
	}

autap=autapomorphies(chmatrix, nstates, notu, nchars, UNKNOWN, INAP);

for (a=0; a<notu; ++a)
	character1[a]=character2[a]=test1[a]=test2[a]=0;

tallied=ivector(maxst+1);
taxst=ivector(maxst+1);
mtch=ivector(maxst+1);

for (ch1=0; ch1<nchars; ch1++)	{
	/* if the character is an autapomorphy or invariant, do not waste your time */
	/*    the character must be compatible with all other characters	*/
	for (ch1=ch1; autap[ch1]<2 && ch1<nchars; ++ch1)	{
		comatrix[ch1][ch1]=0;
		for (ch2=ch1+1; ch2<nchars; ++ch2)	{
			if (comptype==0)	comatrix[ch1][ch2]=comatrix[ch2][ch1]=1;
			else				comatrix[ch1][ch2]=comatrix[ch2][ch1]=0;
			}	
		}	

	if (ch1>=nchars)	break;
	
	comatrix[ch1][ch1]=0;

	/* Make sure that outgroup state is appropriately coded	
	if (nstates[ch1]==2 && chmatrix[OUTGROUP][ch1]==1)	{
		for (sp=0; sp<notu; ++sp)	{
			if (chmatrix[sp][ch1]==0)			chmatrix[sp][ch1]=1;
			else if (chmatrix[sp][ch1]==1)		chmatrix[sp][ch1]=0;
			}
		}					*/
	if (ch1>=nchars)	break;	
	for (ch2=ch1+1; ch2<nchars; ++ch2)	{
		mxst = nstates[ch1];
		for (sp=0; sp<notu; ++sp)			character1[sp]=chmatrix[sp][ch1];
		/* characters are compatible if second is autapomorphic */
		for (ch2=ch2; (autap[ch2]<2 && ch2<nchars); ++ch2)	{
			if (comptype==0)	comatrix[ch1][ch2]=comatrix[ch2][ch1]=1;
			else				comatrix[ch1][ch2]=comatrix[ch2][ch1]=0;
			}	
		if (ch2>=nchars)	break;
		comatrix[ch1][ch2]=comatrix[ch2][ch1]=0;
		clearivector(test1,notu,RAND_MAX);
		clearivector(test2,notu,RAND_MAX);
		for (sp=0; sp<notu; ++sp)
			character2[sp]=chmatrix[sp][ch2];

		/*** Determine compatibility ***/		
		incompatible = 0;

		if ((nstates[ch1]>nstates[ch2] && type[ch1]==1) || (nstates[ch2]>nstates[ch1] && type[ch2]==0))	{
			for (sp=0; sp<notu; ++sp)	{
				b = character2[sp];
				character2[sp]=character1[sp];
				character1[sp]=b;
				mxst = nstates[ch2];
				}
			hi=ch2;
			lo=ch1;
			}
		else	{
			hi=ch1;
			lo=ch2;
			}

		/* sort on placers for comparisons */
		/* rewrite to ignore unknowns & inapplicables */
		mx2=-1*RAND_MAX;
		mn2=RAND_MAX;
		mn1=RAND_MAX;
		for (sp=0; sp<notu; ++sp)	{
			if ((character1[sp]!=UNKNOWN && character1[sp]!=INAP) && character1[sp]<mn1)	mn1=character1[sp];
			if ((character2[sp]!=UNKNOWN && character2[sp]!=INAP) && character2[sp]>mx2)	mx2=character2[sp];
			if ((character2[sp]!=UNKNOWN && character2[sp]!=INAP) && character2[sp]<mn2)	mn2=character2[sp];
			}
		rev=0;
		for (sp=0; rev==0 && sp<notu; ++sp)	if (character1[sp]==mn1 && character2[sp]==mx2)	rev=1;
		for (sp=0; rev==1 && sp<notu; ++sp)	if (character1[sp]==mn1 && character2[sp]==mn2)	rev=0;
		
		
		compsp=0;	/* this counts the number of comparable species - i.e., those scored for both characters */
		for (sp=0; sp<notu; ++sp)	{
			placer = 0;
			while (sp<notu && ((character1[sp]==UNKNOWN || character1[sp]==INAP) || (character2[sp]==UNKNOWN || character2[sp]==INAP)))
				++sp;
			if (sp>=notu)	break;
			++compsp;
			/**** Sort on 1st character ****/
			for (placer=0; (character1[sp]>test1[placer] && placer< notu); placer=placer)	++placer;
			/**** Sort on 2nd character ****/
			if (rev==0)	for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]>test2[placer]) && placer<notu); placer=placer)	++placer;
			else		for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]<test2[placer]) && placer<notu); placer=placer)	++placer;

			for (b=sp; b>placer; --b)	{
				test1[b]=test1[b-1];
				test2[b]=test2[b-1];
				}
			test1[placer] = character1[sp];
			test2[placer] = character2[sp];
			}
			
		for (sp=compsp; sp<notu; ++sp)	test1[sp]=test2[sp]=UNKNOWN;

/*		ms=0;
		for (sp=0; sp<notu-ms; ++sp)	{
			if ((test1[sp]==UNKNOWN || test1[sp]==INAP) || (test2[sp]==UNKNOWN || test2[sp]==INAP))	{
				for (b=sp; b<(notu-ms)-1;	++b)	{
					test1[b]=test1[b+1];
					test2[b]=test2[b+1];
					}
				--sp;
				++ms;
				}
			}	*/
		
		if (compsp>0)	{	
			/* for ease of computation, set lowest state to 0 */
			if (test1[0]>0 && (test1[0]!=INAP && test1[0]!=UNKNOWN))	{
				mn1=test1[0];
				if (mn1>0)	{
					for (sp=0; sp<compsp; ++sp)	{
						if (test1[sp]!=INAP && test1[sp]!=UNKNOWN)	{
							test1[sp]=test1[sp]-mn1;
							}
						else	sp=compsp;
						}
					}
				}

			/**** Routine for Binary Characters ****/
			if (nstates[ch1]==2 && nstates[ch2]==2)
				incompatible=binarycompatible(test1, test2, compsp, UNKNOWN, INAP);
			
			/**** Routine for Multistate Characters ****/
			else	{
				/* if ordered multistate	*/
				if (type[hi]==0)
					incompatible=ordcompatible(test1, test2,nstates[hi], nstates[lo], compsp, UNKNOWN, INAP);
					/* end routine for ordered multistates	*/
				else
					incompatible=unordcompatible(test1, test2, nstates[hi], nstates[lo], compsp, UNKNOWN, INAP);
					/* end compatibility test for unordered characters */
				}	/* end routine for multstates	*/

			inhiercompat = 1;
			if (incompatible==0)	{
				if (comptype==0)	{
					comatrix[ch1][ch2]=comatrix[ch2][ch1]=1;
					}
				else	{						
					a=0;
					/* simple test if both are binary: */
					if (nstates[ch1]==2 && nstates[ch2]==2)	{
						a=0;
						/* first look for 00 & 01 */
						for (a=0; test1[a]==mnst[ch1] && a<compsp; ++a)	a=a;
						for (sp=0; sp<=a; ++sp)	{
							/* 2nd state found?	*/
							if (test2[sp]!=test2[0])	{
								inhiercompat = 0;
								sp = compsp;
								}				
							}
						/* if only 00's found (or 01's), make sure that two states are found with 1- */
						if (inhiercompat==1)	{
							b=a;
							for (b=a; (test1[b]==1+mnst[ch1] && b<compsp); ++b)	b=b;
							for (sp=a; sp<=b; ++sp)	{
								if (test2[sp]!=test2[a])	{
									inhiercompat = 0;
									sp = compsp;
									}
								} 
							}
						}
					/* if two multistates with different # states OR different numbers of taxa
					 with derived conditions, then they must be HC if compatible at all 
					 (unless they are "complements")				*/
					
					else	{
						/* check to see if any derived states have multiple counterparts	*/
						for (d=0; d<compsp; ++d)	{
							if (test1[d]>0 && test1[d]<mxst)	{
								if (test2[d]>0 && test2[d]<mxst)	{
									inhiercompat = 0;
									d=compsp;
									}
								}
							}
						/* even if all deriveds for one state paired with single derived, this might not be true */
						if (inhiercompat==1)	{
							for (d=0; d<compsp; ++d)	{
								if (test2[d]>0 && test2[d]<mxst)	{
									if (test1[d]>0 && test1[d]<mxst)	{
										inhiercompat = 0;
										d=compsp;
										}
									}
								}
							}
						}
					}
				if (inhiercompat==0)	{
					comatrix[ch1][ch2]=comatrix[ch2][ch1]=1;
					/* 1 for hierarchically compatible, 0 for hierarchically incompatible	*/
					}
				}
			else	comatrix[ch1][ch2]=comatrix[ch2][ch1]=0;
			}	/* only look for compatibility if there is anything to compare */
		}	/* end comparison between ch1 & ch2 */
	}

for (ch1=0; ch1<nchars; ++ch1)	{
	comatrix[ch1][ch1]=0;
	for (ch2=0; ch2<nchars; ++ch2)	{
		if (ch2!=ch1)	comatrix[ch1][ch1]=comatrix[ch1][ch1]+comatrix[ch1][ch2];
		}
	}

free_ivector(autap);
free_ivector(character1);
free_ivector(character2);
free_ivector(mnst);
free_ivector(mtch);
free_ivector(tallied);
free_ivector(taxst);
free_ivector(test1);
free_ivector(test2);

return comatrix;
}


/*
Function returning the number of compatible characters WITHOUT requiring a compatibility matrix.
	IF you have a compatiblity matrix, use int count_comp instead!
Neads:
	nstates: #states per character
	notu: number of taxa
	chmatrix: character matrix
	type: character type (0 = ordered, 1 = unordered)
	nchars: number of characters
	comptype: 0: general compatibility; 1: hierarchical compatibilty
	OUTGROUP: number of outgroup taxon
	UNKNOWN: value for "?"
	INAP: value for inapplicable
*****************************************************************************/
int nu_comp(int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP)
{
int a, b, d, ch1, ch2, sp, placer, incompatible, inhiercompat, mxst, hi, lo, compsp, maxst;
int	mn1, mn2, mx2, rev;
int *autap, *character1, *character2, *mnst, *mtch, *tallied, *taxst, *test1, *test2;
int CP=0;

a=outgroup;		/* delete this if you do not restore hierarchical compatibility	*/

mnst=ivector(nchars);
character1=ivector(2*notu);		/* this should not be necessary	*/
character2=ivector(2*notu);		/* this should not be necessary	*/
test1=ivector(2*notu);			/* this should not be necessary	*/
test2=ivector(2*notu);			/* this should not be necessary	*/

maxst=maxiarray(nstates,nchars);
for (a=0; a<nchars; ++a)	{
	mnst[a]=100;
	for (b=0; b<notu; ++b)	{
		if (chmatrix[b][a]<mnst[a] && (chmatrix[b][a]!=UNKNOWN && chmatrix[b][a]!=INAP))
			mnst[a]=chmatrix[b][a];
		}
	}

autap=autapomorphies(chmatrix, nstates, notu, nchars, UNKNOWN, INAP);

for (a=0; a<notu; ++a)
	character1[a]=character2[a]=test1[a]=test2[a]=0;

tallied=ivector(maxst+1);
taxst=ivector(maxst+1);
mtch=ivector(maxst+1);

for (ch1=0; ch1<nchars; ch1++)	{
	/* if the character is an autapomorphy or invariant, do not waste your time */
	/*    the character must be compatible with all other characters	*/
	for (ch1=ch1; autap[ch1]<2 && ch1<nchars; ++ch1)	{
		for (ch2=ch1+1; ch2<nchars; ++ch2)	{
			if (comptype==0)	++CP;
			}	
		}	

	if (ch1>=nchars)	break;
	
	/* Make sure that outgroup state is appropriately coded	
	if (nstates[ch1]==2 && chmatrix[OUTGROUP][ch1]==1)	{
		for (sp=0; sp<notu; ++sp)	{
			if (chmatrix[sp][ch1]==0)			chmatrix[sp][ch1]=1;
			else if (chmatrix[sp][ch1]==1)		chmatrix[sp][ch1]=0;
			}
		}					*/
	if (ch1>=nchars)	break;	
	for (ch2=ch1+1; ch2<nchars; ++ch2)	{
		mxst = nstates[ch1];
		for (sp=0; sp<notu; ++sp)			character1[sp]=chmatrix[sp][ch1];
		/* characters are compatible if second is autapomorphic */
		for (ch2=ch2; (autap[ch2]<2 && ch2<nchars); ++ch2)	{
			if (comptype==0)	++CP;
			}	
		if (ch2>=nchars)	break;
		clearivector(test1,notu,RAND_MAX);
		clearivector(test2,notu,RAND_MAX);
		for (sp=0; sp<notu; ++sp)
			character2[sp]=chmatrix[sp][ch2];

		/*** Determine compatibility ***/		
		incompatible = 0;

		if ((nstates[ch1]>nstates[ch2] && type[ch1]==1) || (nstates[ch2]>nstates[ch1] && type[ch2]==0))	{
			for (sp=0; sp<notu; ++sp)	{
				b = character2[sp];
				character2[sp]=character1[sp];
				character1[sp]=b;
				mxst = nstates[ch2];
				}
			hi=ch2;
			lo=ch1;
			}
		else	{
			hi=ch1;
			lo=ch2;
			}

		/* sort on placers for comparisons */
		/* rewrite to ignore unknowns & inapplicables */
		mx2=-1*RAND_MAX;
		mn2=RAND_MAX;
		mn1=RAND_MAX;
		for (sp=0; sp<notu; ++sp)	{
			if ((character1[sp]!=UNKNOWN && character1[sp]!=INAP) && character1[sp]<mn1)	mn1=character1[sp];
			if ((character2[sp]!=UNKNOWN && character2[sp]!=INAP) && character2[sp]>mx2)	mx2=character2[sp];
			if ((character2[sp]!=UNKNOWN && character2[sp]!=INAP) && character2[sp]<mn2)	mn2=character2[sp];
			}
		rev=0;
		for (sp=0; rev==0 && sp<notu; ++sp)	if (character1[sp]==mn1 && character2[sp]==mx2)	rev=1;
		for (sp=0; rev==1 && sp<notu; ++sp)	if (character1[sp]==mn1 && character2[sp]==mn2)	rev=0;
		
		
		compsp=0;	/* this counts the number of comparable species - i.e., those scored for both characters */
		for (sp=0; sp<notu; ++sp)	{
			placer = 0;
			while (sp<notu && ((character1[sp]==UNKNOWN || character1[sp]==INAP) || (character2[sp]==UNKNOWN || character2[sp]==INAP)))
				++sp;
			if (sp>=notu)	break;
			++compsp;
			/**** Sort on 1st character ****/
			for (placer=0; (character1[sp]>test1[placer] && placer< notu); placer=placer)	++placer;
			/**** Sort on 2nd character ****/
			if (rev==0)	for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]>test2[placer]) && placer<notu); placer=placer)	++placer;
			else		for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]<test2[placer]) && placer<notu); placer=placer)	++placer;

			for (b=sp; b>placer; --b)	{
				test1[b]=test1[b-1];
				test2[b]=test2[b-1];
				}
			test1[placer] = character1[sp];
			test2[placer] = character2[sp];
			}
			
		for (sp=compsp; sp<notu; ++sp)	test1[sp]=test2[sp]=UNKNOWN;

		if (compsp>0)	{	
			/* for ease of computation, set lowest state to 0 */
			if (test1[0]>0 && (test1[0]!=INAP && test1[0]!=UNKNOWN))	{
				mn1=test1[0];
				if (mn1>0)	{
					for (sp=0; sp<compsp; ++sp)	{
						if (test1[sp]!=INAP && test1[sp]!=UNKNOWN)	{
							test1[sp]=test1[sp]-mn1;
							}
						else	sp=compsp;
						}
					}
				}

			/**** Routine for Binary Characters ****/
			if (nstates[ch1]==2 && nstates[ch2]==2)
				incompatible=binarycompatible(test1, test2, compsp, UNKNOWN, INAP);
			
			/**** Routine for Multistate Characters ****/
			else	{
				/* if ordered multistate	*/
				if (type[hi]==0)
					incompatible=ordcompatible(test1, test2,nstates[hi], nstates[lo], compsp, UNKNOWN, INAP);
					/* end routine for ordered multistates	*/
				else
					incompatible=unordcompatible(test1, test2, nstates[hi], nstates[lo], compsp, UNKNOWN, INAP);
					/* end compatibility test for unordered characters */
				}	/* end routine for multstates	*/

			inhiercompat = 1;
			if (incompatible==0)	{
				if (comptype==0)	++CP;
				else	{						
					a=0;
					/* simple test if both are binary: */
					if (nstates[ch1]==2 && nstates[ch2]==2)	{
						a=0;
						/* first look for 00 & 01 */
						for (a=0; test1[a]==mnst[ch1] && a<compsp; ++a)	a=a;
						for (sp=0; sp<=a; ++sp)	{
							/* 2nd state found?	*/
							if (test2[sp]!=test2[0])	{
								inhiercompat = 0;
								sp = compsp;
								}				
							}
						/* if only 00's found (or 01's), make sure that two states are found with 1- */
						if (inhiercompat==1)	{
							b=a;
							for (b=a; (test1[b]==1+mnst[ch1] && b<compsp); ++b)	b=b;
							for (sp=a; sp<=b; ++sp)	{
								if (test2[sp]!=test2[a])	{
									inhiercompat = 0;
									sp = compsp;
									}
								} 
							}
						}
					/* if two multistates with different # states OR different numbers of taxa
					 with derived conditions, then they must be HC if compatible at all 
					 (unless they are "complements")				*/
					
					else	{
						/* check to see if any derived states have multiple counterparts	*/
						for (d=0; d<compsp; ++d)	{
							if (test1[d]>0 && test1[d]<mxst)	{
								if (test2[d]>0 && test2[d]<mxst)	{
									inhiercompat = 0;
									d=compsp;
									}
								}
							}
						/* even if all deriveds for one state paired with single derived, this might not be true */
						if (inhiercompat==1)	{
							for (d=0; d<compsp; ++d)	{
								if (test2[d]>0 && test2[d]<mxst)	{
									if (test1[d]>0 && test1[d]<mxst)	{
										inhiercompat = 0;
										d=compsp;
										}
									}
								}
							}
						}
					}
				if (inhiercompat==0)	++CP;
				}
			}	/* only look for compatibility if there is anything to compare */
		}	/* end comparison between ch1 & ch2 */
	}

free_ivector(autap);
free_ivector(character1);
free_ivector(character2);
free_ivector(mnst);
free_ivector(mtch);
free_ivector(tallied);
free_ivector(taxst);
free_ivector(test1);
free_ivector(test2);

return CP;
}

/*
Function returning the number of compatiblities per character WITHOUT requiring a compatibility matrix.
Neads:
	nstates: #states per character
	notu: number of taxa
	chmatrix: character matrix
	type: character type (0 = ordered, 1 = unordered)
	nchars: number of characters
	comptype: 0: general compatibility; 1: hierarchical compatibilty
	OUTGROUP: number of outgroup taxon
	UNKNOWN: value for "?"
	INAP: value for inapplicable
	charcomps: the array giving the number of compatibilities per character; this is allocated in advance
		in order to discourage memory errors from repeated allocations/deallocations.
*****************************************************************************/
unsigned long *char_comp(int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP)
{
int a, b, d, ch1, ch2, sp, placer, incompatible, inhiercompat, mxst, hi, lo, compsp, maxst;
int	mn1, mn2, mx2, rev;
int *autap, *character1, *character2, *mnst, *mtch, *tallied, *taxst, *test1, *test2;
unsigned long *charcomps;

a=outgroup;		/* delete this if you do not restore hierarchical compatibility	*/

charcomps=ulvector(nchars);
mnst=ivector(nchars);
character1=ivector(notu);
character2=ivector(notu);
test1=ivector(notu);
test2=ivector(notu);

maxst=maxiarray(nstates,nchars);
for (a=0; a<nchars; ++a)	{
	mnst[a]=100;
	for (b=0; b<notu; ++b)	{
		if (chmatrix[b][a]<mnst[a] && (chmatrix[b][a]!=UNKNOWN && chmatrix[b][a]!=INAP))
			mnst[a]=chmatrix[b][a];
		}
	}

autap=autapomorphies(chmatrix, nstates, notu, nchars, UNKNOWN, INAP);

for (a=0; a<notu; ++a)
	character1[a]=character2[a]=test1[a]=test2[a]=0;

tallied=ivector(maxst+1);
taxst=ivector(maxst+1);
mtch=ivector(maxst+1);

for (ch1=0; ch1<nchars; ch1++)	{
	/* if the character is an autapomorphy or invariant, do not waste your time */
	/*    the character must be compatible with all other characters	*/
	for (ch1=ch1; autap[ch1]<2 && ch1<nchars; ++ch1)	{
		for (ch2=ch1+1; ch2<nchars; ++ch2)	{
			if (comptype==0)	{
				++charcomps[ch1];
				++charcomps[ch2];
				}
			}	
		}	

	if (ch1>=nchars)	break;
	
	/* Make sure that outgroup state is appropriately coded	
	if (nstates[ch1]==2 && chmatrix[OUTGROUP][ch1]==1)	{
		for (sp=0; sp<notu; ++sp)	{
			if (chmatrix[sp][ch1]==0)			chmatrix[sp][ch1]=1;
			else if (chmatrix[sp][ch1]==1)		chmatrix[sp][ch1]=0;
			}
		}					*/
	if (ch1>=nchars)	break;	
	for (ch2=ch1+1; ch2<nchars; ++ch2)	{
		mxst = nstates[ch1];
		for (sp=0; sp<notu; ++sp)			character1[sp]=chmatrix[sp][ch1];
		/* characters are compatible if second is autapomorphic */
		for (ch2=ch2; (autap[ch2]<2 && ch2<nchars); ++ch2)	{
			if (comptype==0)	{
				++charcomps[ch1];
				++charcomps[ch2];
				}
			}	
		if (ch2>=nchars)	break;
		clearivector(test1,notu,RAND_MAX);
		clearivector(test2,notu,RAND_MAX);
		for (sp=0; sp<notu; ++sp)
			character2[sp]=chmatrix[sp][ch2];

		/*** Determine compatibility ***/		
		incompatible = 0;

		if ((nstates[ch1]>nstates[ch2] && type[ch1]==1) || (nstates[ch2]>nstates[ch1] && type[ch2]==0))	{
			for (sp=0; sp<notu; ++sp)	{
				b = character2[sp];
				character2[sp]=character1[sp];
				character1[sp]=b;
				mxst = nstates[ch2];
				}
			hi=ch2;
			lo=ch1;
			}
		else	{
			hi=ch1;
			lo=ch2;
			}

		/* sort on placers for comparisons */
		/* rewrite to ignore unknowns & inapplicables */
		mx2=-1*RAND_MAX;
		mn2=RAND_MAX;
		mn1=RAND_MAX;
		for (sp=0; sp<notu; ++sp)	{
			if ((character1[sp]!=UNKNOWN && character1[sp]!=INAP) && character1[sp]<mn1)	mn1=character1[sp];
			if ((character2[sp]!=UNKNOWN && character2[sp]!=INAP) && character2[sp]>mx2)	mx2=character2[sp];
			if ((character2[sp]!=UNKNOWN && character2[sp]!=INAP) && character2[sp]<mn2)	mn2=character2[sp];
			}
		rev=0;
		for (sp=0; rev==0 && sp<notu; ++sp)	if (character1[sp]==mn1 && character2[sp]==mx2)	rev=1;
		for (sp=0; rev==1 && sp<notu; ++sp)	if (character1[sp]==mn1 && character2[sp]==mn2)	rev=0;
		
		
		compsp=0;	/* this counts the number of comparable species - i.e., those scored for both characters */
		for (sp=0; sp<notu; ++sp)	{
			placer = 0;
			while (sp<notu && ((character1[sp]==UNKNOWN || character1[sp]==INAP) || (character2[sp]==UNKNOWN || character2[sp]==INAP)))
				++sp;
			if (sp>=notu)	break;
			++compsp;
			/**** Sort on 1st character ****/
			for (placer=0; (character1[sp]>test1[placer] && placer< notu); placer=placer)	++placer;
			/**** Sort on 2nd character ****/
			if (rev==0)	for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]>test2[placer]) && placer<notu); placer=placer)	++placer;
			else		for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]<test2[placer]) && placer<notu); placer=placer)	++placer;

			for (b=sp; b>placer; --b)	{
				test1[b]=test1[b-1];
				test2[b]=test2[b-1];
				}
			test1[placer] = character1[sp];
			test2[placer] = character2[sp];
			}
			
		for (sp=compsp; sp<notu; ++sp)	test1[sp]=test2[sp]=UNKNOWN;

		if (compsp>0)	{	
			/* for ease of computation, set lowest state to 0 */
			if (test1[0]>0 && (test1[0]!=INAP && test1[0]!=UNKNOWN))	{
				mn1=test1[0];
				if (mn1>0)	{
					for (sp=0; sp<compsp; ++sp)	{
						if (test1[sp]!=INAP && test1[sp]!=UNKNOWN)	{
							test1[sp]=test1[sp]-mn1;
							}
						else	sp=compsp;
						}
					}
				}

			/**** Routine for Binary Characters ****/
			if (nstates[ch1]==2 && nstates[ch2]==2)
				incompatible=binarycompatible(test1, test2, compsp, UNKNOWN, INAP);
			
			/**** Routine for Multistate Characters ****/
			else	{
				/* if ordered multistate	*/
				if (type[hi]==0)
					incompatible=ordcompatible(test1, test2,nstates[hi], nstates[lo], compsp, UNKNOWN, INAP);
					/* end routine for ordered multistates	*/
				else
					incompatible=unordcompatible(test1, test2, nstates[hi], nstates[lo], compsp, UNKNOWN, INAP);
					/* end compatibility test for unordered characters */
				}	/* end routine for multstates	*/

			inhiercompat = 1;
			if (incompatible==0)	{
				if (comptype==0)	{
					++charcomps[ch1];
					++charcomps[ch2];
					}
				else	{						
					a=0;
					/* simple test if both are binary: */
					if (nstates[ch1]==2 && nstates[ch2]==2)	{
						a=0;
						/* first look for 00 & 01 */
						for (a=0; test1[a]==mnst[ch1] && a<compsp; ++a)	a=a;
						for (sp=0; sp<=a; ++sp)	{
							/* 2nd state found?	*/
							if (test2[sp]!=test2[0])	{
								inhiercompat = 0;
								sp = compsp;
								}				
							}
						/* if only 00's found (or 01's), make sure that two states are found with 1- */
						if (inhiercompat==1)	{
							b=a;
							for (b=a; (test1[b]==1+mnst[ch1] && b<compsp); ++b)	b=b;
							for (sp=a; sp<=b; ++sp)	{
								if (test2[sp]!=test2[a])	{
									inhiercompat = 0;
									sp = compsp;
									}
								} 
							}
						}
					/* if two multistates with different # states OR different numbers of taxa
					 with derived conditions, then they must be HC if compatible at all 
					 (unless they are "complements")				*/
					
					else	{
						/* check to see if any derived states have multiple counterparts	*/
						for (d=0; d<compsp; ++d)	{
							if (test1[d]>0 && test1[d]<mxst)	{
								if (test2[d]>0 && test2[d]<mxst)	{
									inhiercompat = 0;
									d=compsp;
									}
								}
							}
						/* even if all deriveds for one state paired with single derived, this might not be true */
						if (inhiercompat==1)	{
							for (d=0; d<compsp; ++d)	{
								if (test2[d]>0 && test2[d]<mxst)	{
									if (test1[d]>0 && test1[d]<mxst)	{
										inhiercompat = 0;
										d=compsp;
										}
									}
								}
							}
						}
					}
				if (inhiercompat==0)	{
					++charcomps[ch1];
					++charcomps[ch2];
					}
				}
			}	/* only look for compatibility if there is anything to compare */
		}	/* end comparison between ch1 & ch2 */
	}

free_ivector(autap);
free_ivector(character1);
free_ivector(character2);
free_ivector(mnst);
free_ivector(mtch);
free_ivector(tallied);
free_ivector(taxst);
free_ivector(test1);
free_ivector(test2);

return charcomps;
}

/*
Function returning the number of compatibilities given a compatibility matrix.
Neads:
	comatrix: char x char compatibility matrix, with 1 = compatible and 0 = incompatible
	nchars: number of characters
*****************************************************************************/
int countcomp(unsigned long **comatrix, int nchars)
{
int	c1, c2, CP=0;

for (c1=0; c1<nchars-1; ++c1)	{
	for (c2=c1+1; c2<nchars; ++c2)	{
		if (comatrix[c1][c2]==1)
			++CP;
		}
	}
return CP;
}

/*
Function returning the number of compatibilities per character given a compatibility matrix.
Neads:
	comatrix: char x char compatibility matrix, with 1 = compatible and 0 = incompatible
	nchars: number of characters
	charcomps: the array giving compatibilities per character (this is allocated in advance
		in order to discourage memory errors from repeated allocations / deallocations.)
*****************************************************************************/
unsigned long *countcharcomps(unsigned long **comatrix, int nchars, unsigned long *charcomps)
{
int	c1, c2;

for (c1=0; c1<nchars; ++c1)	charcomps[c1]=0;

for (c1=0; c1<nchars-1; ++c1)	{
	for (c2=c1+1; c2<nchars; ++c2)	{
		if (comatrix[c1][c2]==1)	{
			++charcomps[c1];
			++charcomps[c2];
			}
		}
	}
return charcomps;
}

/*
Function returning the number of mutual compatibilities (see OÕKeefe & Wagner 2001 Syst. Biol.
 	given a compatibility matrix.
Neads:
	comatrix: char x char compatibility matrix, with 1 = compatible and 0 = incompatible
	nchars: number of characters
*****************************************************************************/
unsigned long **mutualcomp(unsigned long **comatrix, int nchars)
{
int	c1, c2, c3;
unsigned long **m;

m=ulmatrix(nchars,nchars);

for (c1=0; c1<nchars; ++c1)	{
	m[c1][c1]=comatrix[c1][c1];
	for (c2=c1+1; c2<nchars; ++c2)	{
		m[c1][c2]=0;
		for (c3=0; c3<nchars; ++c3)	{
			while (c3==c1 || c3==c2)	++c3;
			if (c3>=nchars)	break;
			if (comatrix[c1][c3]==1 && comatrix[c2][c3]==1)
				++m[c1][c2];
			}
		m[c2][c1]=m[c1][c2];
		}
	}
return m;
}

/*
Function returning the number of compatible characters for a particular character
	WITHOUT requiring a compatibility matrix in advance.

Neads:
	ch1: character being examined
	nstates: #states per character
	notu: number of taxa
	chmatrix: character matrix
	type: character type (0 = ordered, 1 = unordered)
	nchars: number of characters
	comptype: 0: general compatibility; 1: hierarchical compatibilty
	OUTGROUP: number of outgroup taxon
	UNKNOWN: value for "?"
	INAP: value for inapplicable
*****************************************************************************/
int char_nu_comp(int ch1, int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP)
{
int a, b, d, ch2, sp, placer, incompatible, inhiercompat, mxst, hi, lo, compsp, maxst;
int	mn1, mn2, mx2, rev;
int *autap, *character1, *character2, *mnst, *mtch, *tallied, *taxst, *test1, *test2;
int numcmp=0;

a=outgroup;		/* delete this if you do not restore hierarchical compatibility	*/

mnst=ivector(nchars);
character1=ivector(2*notu);	/* should be only notu	*/
character2=ivector(2*notu);	/* should be only notu	*/
test1=ivector(2*notu);		/* should be only notu	*/
test2=ivector(2*notu);		/* should be only notu	*/

maxst=maxiarray(nstates,nchars);
for (a=0; a<nchars; ++a)	{
	mnst[a]=100;
	for (b=0; b<notu; ++b)	{
		if (chmatrix[b][a]<mnst[a] && (chmatrix[b][a]!=UNKNOWN && chmatrix[b][a]!=INAP))
			mnst[a]=chmatrix[b][a];
		}
	}

autap=autapomorphies(chmatrix, nstates, notu, nchars, UNKNOWN, INAP);

for (a=0; a<notu; ++a)
	character1[a]=character2[a]=test1[a]=test2[a]=0;

tallied=ivector(maxst+1);
taxst=ivector(maxst+1);
mtch=ivector(maxst+1);

/* Make sure that outgroup state is appropriately coded	
if (nstates[ch1]==2 && chmatrix[OUTGROUP][ch1]==1)	{
	for (sp=0; sp<notu; ++sp)	{
		if (chmatrix[sp][ch1]==0)			chmatrix[sp][ch1]=1;
		else if (chmatrix[sp][ch1]==1)		chmatrix[sp][ch1]=0;
		}
	}	*/

if (autap[ch1]==nstates[ch1]-1)	numcmp=nchars-1;
else	{
	for (ch2=0; ch2<nchars; ++ch2)	{
		if (ch2==ch1)		++ch2;
		if (ch2>=nchars)	break;
		
		mxst = nstates[ch1];
		for (sp=0; sp<notu; ++sp)			character1[sp]=chmatrix[sp][ch1];
		/* characters are compatible if second is autapomorphic */
		for (ch2=ch2; (autap[ch2]<2 && ch2<nchars); ++ch2)	{
			if (comptype==0)	++numcmp;
			}	
		if (ch2>=nchars)	break;
		clearivector(test1,notu,RAND_MAX);
		clearivector(test2,notu,RAND_MAX);
		for (sp=0; sp<notu; ++sp)
			character2[sp]=chmatrix[sp][ch2];

		/*** Determine compatibility ***/		
		incompatible = 0;

		if ((nstates[ch1]>nstates[ch2] && type[ch1]==1) || (nstates[ch2]>nstates[ch1] && type[ch2]==0))	{
			for (sp=0; sp<notu; ++sp)	{
				b = character2[sp];
				character2[sp]=character1[sp];
				character1[sp]=b;
				mxst = nstates[ch2];
				}
			hi=ch2;
			lo=ch1;
			}
		else	{
			hi=ch1;
			lo=ch2;
			}

		/* sort on placers for comparisons */
		/* rewrite to ignore unknowns & inapplicables */
		mx2=-1*RAND_MAX;
		mn2=RAND_MAX;
		mn1=RAND_MAX;
		for (sp=0; sp<notu; ++sp)	{
			if ((character1[sp]!=UNKNOWN && character1[sp]!=INAP) && character1[sp]<mn1)	mn1=character1[sp];
			if ((character2[sp]!=UNKNOWN && character2[sp]!=INAP) && character2[sp]>mx2)	mx2=character2[sp];
			if ((character2[sp]!=UNKNOWN && character2[sp]!=INAP) && character2[sp]<mn2)	mn2=character2[sp];
			}
		rev=0;
		for (sp=0; rev==0 && sp<notu; ++sp)	if (character1[sp]==mn1 && character2[sp]==mx2)	rev=1;
		for (sp=0; rev==1 && sp<notu; ++sp)	if (character1[sp]==mn1 && character2[sp]==mn2)	rev=0;
		
		
		compsp=0;	/* this counts the number of comparable species - i.e., those scored for both characters */
		for (sp=0; sp<notu; ++sp)	{
			placer = 0;
			while (sp<notu && ((character1[sp]==UNKNOWN || character1[sp]==INAP) || (character2[sp]==UNKNOWN || character2[sp]==INAP)))
				++sp;
			if (sp>=notu)	break;
			++compsp;
			/**** Sort on 1st character ****/
			for (placer=0; (character1[sp]>test1[placer] && placer< notu); placer=placer)	++placer;
			/**** Sort on 2nd character ****/
			if (rev==0)	for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]>test2[placer]) && placer<notu); placer=placer)	++placer;
			else		for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]<test2[placer]) && placer<notu); placer=placer)	++placer;

			for (b=sp; b>placer; --b)	{
				test1[b]=test1[b-1];
				test2[b]=test2[b-1];
				}
			test1[placer] = character1[sp];
			test2[placer] = character2[sp];
			}
			
		for (sp=compsp; sp<notu; ++sp)	test1[sp]=test2[sp]=UNKNOWN;

		if (compsp>0)	{	
			/* for ease of computation, set lowest state to 0 */
			if (test1[0]>0 && (test1[0]!=INAP && test1[0]!=UNKNOWN))	{
				mn1=test1[0];
				if (mn1>0)	{
					for (sp=0; sp<compsp; ++sp)	{
						if (test1[sp]!=INAP && test1[sp]!=UNKNOWN)	{
							test1[sp]=test1[sp]-mn1;
							}
						else	sp=compsp;
						}
					}
				}

			/**** Routine for Binary Characters ****/
			if (nstates[ch1]==2 && nstates[ch2]==2)
				incompatible=binarycompatible(test1, test2, compsp, UNKNOWN, INAP);
			
			/**** Routine for Multistate Characters ****/
			else	{
				/* if ordered multistate	*/
				if (type[hi]==0)
					incompatible=ordcompatible(test1, test2,nstates[hi], nstates[lo], compsp, UNKNOWN, INAP);
					/* end routine for ordered multistates	*/
				else
					incompatible=unordcompatible(test1, test2, nstates[hi], nstates[lo], compsp, UNKNOWN, INAP);
					/* end compatibility test for unordered characters */
				}	/* end routine for multstates	*/

			inhiercompat = 1;
			if (incompatible==0)	{
				if (comptype==0)		++numcmp;

				else	{						
					a=0;
					/* simple test if both are binary: */
					if (nstates[ch1]==2 && nstates[ch2]==2)	{
						a=0;
						/* first look for 00 & 01 */
						for (a=0; test1[a]==mnst[ch1] && a<compsp; ++a)	a=a;
						for (sp=0; sp<=a; ++sp)	{
							/* 2nd state found?	*/
							if (test2[sp]!=test2[0])	{
								inhiercompat = 0;
								sp = compsp;
								}				
							}
						/* if only 00's found (or 01's), make sure that two states are found with 1- */
						if (inhiercompat==1)	{
							b=a;
							for (b=a; (test1[b]==1+mnst[ch1] && b<compsp); ++b)	b=b;
							for (sp=a; sp<=b; ++sp)	{
								if (test2[sp]!=test2[a])	{
									inhiercompat = 0;
									sp = compsp;
									}
								} 
							}
						}
					/* if two multistates with different # states OR different numbers of taxa
					 with derived conditions, then they must be HC if compatible at all 
					 (unless they are "complements")				*/
					
					else	{
						/* check to see if any derived states have multiple counterparts	*/
						for (d=0; d<compsp; ++d)	{
							if (test1[d]>0 && test1[d]<mxst)	{
								if (test2[d]>0 && test2[d]<mxst)	{
									inhiercompat = 0;
									d=compsp;
									}
								}
							}
						/* even if all deriveds for one state paired with single derived, this might not be true */
						if (inhiercompat==1)	{
							for (d=0; d<compsp; ++d)	{
								if (test2[d]>0 && test2[d]<mxst)	{
									if (test1[d]>0 && test1[d]<mxst)	{
										inhiercompat = 0;
										d=compsp;
										}
									}
								}
							}
						}
					}
				if (inhiercompat==0)	++numcmp;
				}
			}	/* only look for compatibility if there is anything to compare */
		}	/* end comparison between ch1 & ch2 */
	}	/* end examination of compatibility with other characters	*/
free_ivector(autap);
free_ivector(character1);
free_ivector(character2);
free_ivector(mnst);
free_ivector(mtch);
free_ivector(tallied);
free_ivector(taxst);
free_ivector(test1);
free_ivector(test2);

return numcmp;
}

/*
Function returning the number of compatible characters for a particular character
	WITHOUT requiring a compatibility matrix in advance.

Neads:
	ch1: the 1st character
	ch2: the 2nd character
	nstates: #states per character
	notu: number of taxa
	chmatrix: character matrix
	type: character type (0 = ordered, 1 = unordered)
	nchars: number of characters
	comptype: 0: general compatibility; 1: hierarchical compatibilty
	OUTGROUP: number of outgroup taxon
	UNKNOWN: value for "?"
	INAP: value for inapplicable
*****************************************************************************/
int pair_comp(int ch1, int ch2, int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP)
{
int a, b, d, sp, placer, incompatible, inhiercompat, mxst, hi, lo, compsp, maxst;
int	mn1, mn2, mx2, rev;
int *autap, *character1, *character2, *mnst, *mtch, *tallied, *taxst, *test1, *test2;
int numcmp=0;

a=outgroup;		/* delete this if you do not restore hierarchical compatibility	*/

mnst=ivector(nchars);
character1=ivector(notu);
character2=ivector(notu);
test1=ivector(notu);
test2=ivector(notu);

maxst=maxiarray(nstates,nchars);
for (a=0; a<nchars; ++a)	{
	mnst[a]=100;
	for (b=0; b<notu; ++b)	{
		if (chmatrix[b][a]<mnst[a] && (chmatrix[b][a]!=UNKNOWN && chmatrix[b][a]!=INAP))
			mnst[a]=chmatrix[b][a];
		}
	}

autap=autapomorphies(chmatrix, nstates, notu, nchars, UNKNOWN, INAP);

for (a=0; a<notu; ++a)
	character1[a]=character2[a]=test1[a]=test2[a]=0;

tallied=ivector(maxst+1);
taxst=ivector(maxst+1);
mtch=ivector(maxst+1);

/* Make sure that outgroup state is appropriately coded	
if (nstates[ch1]==2 && chmatrix[OUTGROUP][ch1]==1)	{
	for (sp=0; sp<notu; ++sp)	{
		if (chmatrix[sp][ch1]==0)			chmatrix[sp][ch1]=1;
		else if (chmatrix[sp][ch1]==1)		chmatrix[sp][ch1]=0;
		}
	}	*/

if (autap[ch1]==nstates[ch1]-1)	numcmp=1;

else	{
	mxst = nstates[ch1];
	for (sp=0; sp<notu; ++sp)			character1[sp]=chmatrix[sp][ch1];
	/* characters are compatible if second is autapomorphic */
	for (ch2=ch2; (autap[ch2]<2 && ch2<nchars); ++ch2)	{
		if (comptype==0 && autap[ch2]<2)	++numcmp;
		}
/*	if (numcmp==1)	break;	*/
	clearivector(test1,notu,RAND_MAX);
	clearivector(test2,notu,RAND_MAX);
	for (sp=0; sp<notu; ++sp)
		character2[sp]=chmatrix[sp][ch2];

	/*** Determine compatibility ***/		
	incompatible = 0;

	if ((nstates[ch1]>nstates[ch2] && type[ch1]==1) || (nstates[ch2]>nstates[ch1] && type[ch2]==0))	{
		for (sp=0; sp<notu; ++sp)	{
			b = character2[sp];
			character2[sp]=character1[sp];
			character1[sp]=b;
			mxst = nstates[ch2];
			}
		hi=ch2;
		lo=ch1;
		}
	else	{
		hi=ch1;
		lo=ch2;
		}

	/* sort on placers for comparisons */
	/* rewrite to ignore unknowns & inapplicables */
	mx2=-1*RAND_MAX;
	mn2=RAND_MAX;
	mn1=RAND_MAX;
	for (sp=0; sp<notu; ++sp)	{
		if ((character1[sp]!=UNKNOWN && character1[sp]!=INAP) && character1[sp]<mn1)	mn1=character1[sp];
		if ((character2[sp]!=UNKNOWN && character2[sp]!=INAP) && character2[sp]>mx2)	mx2=character2[sp];
		if ((character2[sp]!=UNKNOWN && character2[sp]!=INAP) && character2[sp]<mn2)	mn2=character2[sp];
		}
	rev=0;
	for (sp=0; rev==0 && sp<notu; ++sp)	if (character1[sp]==mn1 && character2[sp]==mx2)	rev=1;
	for (sp=0; rev==1 && sp<notu; ++sp)	if (character1[sp]==mn1 && character2[sp]==mn2)	rev=0;
	
	
	compsp=0;	/* this counts the number of comparable species - i.e., those scored for both characters */
	for (sp=0; sp<notu; ++sp)	{
		placer = 0;
		while (sp<notu && ((character1[sp]==UNKNOWN || character1[sp]==INAP) || (character2[sp]==UNKNOWN || character2[sp]==INAP)))
			++sp;
		if (sp>=notu)	break;
		++compsp;
		/**** Sort on 1st character ****/
		for (placer=0; (character1[sp]>test1[placer] && placer< notu); placer=placer)	++placer;
		/**** Sort on 2nd character ****/
		if (rev==0)	for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]>test2[placer]) && placer<notu); placer=placer)	++placer;
		else		for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]<test2[placer]) && placer<notu); placer=placer)	++placer;

		for (b=sp; b>placer; --b)	{
			test1[b]=test1[b-1];
			test2[b]=test2[b-1];
			}
		test1[placer] = character1[sp];
		test2[placer] = character2[sp];
		}
		
	for (sp=compsp; sp<notu; ++sp)	test1[sp]=test2[sp]=UNKNOWN;

	if (compsp>0)	{	
		/* for ease of computation, set lowest state to 0 */
		if (test1[0]>0 && (test1[0]!=INAP && test1[0]!=UNKNOWN))	{
			mn1=test1[0];
			if (mn1>0)	{
				for (sp=0; sp<compsp; ++sp)	{
					if (test1[sp]!=INAP && test1[sp]!=UNKNOWN)	{
						test1[sp]=test1[sp]-mn1;
						}
					else	sp=compsp;
					}
				}
			}

		/**** Routine for Binary Characters ****/
		if (nstates[ch1]==2 && nstates[ch2]==2)
			incompatible=binarycompatible(test1, test2, compsp, UNKNOWN, INAP);
		
		/**** Routine for Multistate Characters ****/
		else	{
			/* if ordered multistate	*/
			if (type[hi]==0)
				incompatible=ordcompatible(test1, test2,nstates[hi], nstates[lo], compsp, UNKNOWN, INAP);
				/* end routine for ordered multistates	*/
			else
				incompatible=unordcompatible(test1, test2, nstates[hi], nstates[lo], compsp, UNKNOWN, INAP);
				/* end compatibility test for unordered characters */
			}	/* end routine for multstates	*/

		inhiercompat = 1;
		if (incompatible==0)	{
			if (comptype==0)		++numcmp;

			else	{						
				a=0;
				/* simple test if both are binary: */
				if (nstates[ch1]==2 && nstates[ch2]==2)	{
					a=0;
					/* first look for 00 & 01 */
					for (a=0; test1[a]==mnst[ch1] && a<compsp; ++a)	a=a;
					for (sp=0; sp<=a; ++sp)	{
						/* 2nd state found?	*/
						if (test2[sp]!=test2[0])	{
							inhiercompat = 0;
							sp = compsp;
							}				
						}
					/* if only 00's found (or 01's), make sure that two states are found with 1- */
					if (inhiercompat==1)	{
						b=a;
						for (b=a; (test1[b]==1+mnst[ch1] && b<compsp); ++b)	b=b;
						for (sp=a; sp<=b; ++sp)	{
							if (test2[sp]!=test2[a])	{
								inhiercompat = 0;
								sp = compsp;
								}
							} 
						}
					}
				/* if two multistates with different # states OR different numbers of taxa
				 with derived conditions, then they must be HC if compatible at all 
				 (unless they are "complements")				*/
				
				else	{
					/* check to see if any derived states have multiple counterparts	*/
					for (d=0; d<compsp; ++d)	{
						if (test1[d]>0 && test1[d]<mxst)	{
							if (test2[d]>0 && test2[d]<mxst)	{
								inhiercompat = 0;
								d=compsp;
								}
							}
						}
					/* even if all deriveds for one state paired with single derived, this might not be true */
					if (inhiercompat==1)	{
						for (d=0; d<compsp; ++d)	{
							if (test2[d]>0 && test2[d]<mxst)	{
								if (test1[d]>0 && test1[d]<mxst)	{
									inhiercompat = 0;
									d=compsp;
									}
								}
							}
						}
					}
				}
			if (inhiercompat==0)	++numcmp;
			}
		}	/* only look for compatibility if there is anything to compare */
	}	/* end comparison between ch1 & ch2 */
free_ivector(autap);
free_ivector(character1);
free_ivector(character2);
free_ivector(mnst);
free_ivector(mtch);
free_ivector(tallied);
free_ivector(taxst);
free_ivector(test1);
free_ivector(test2);

return numcmp;
}
/*
Function returning possible numbers of combinations per character.
Neads:
	states: #states per character
	notu: number of taxa
	chmatrix: character matrix
	nchars: number of characters
	UNKNOWN: value for "?"
	INAP: value for inapplicable
*****************************************************************************/
long *PossibleCharCombinations(int *nstates, int nchars)
{
int	 ch1, ch2;
long *combos;

combos=lvector(nchars);
for (ch1=0; ch1<nchars; ++ch1)	{
	combos[ch1]=0;
	for (ch2=0; ch2<nchars; ++ch2)	{
		if (ch2==ch1)	++ch2;
		if (ch2>=nchars)	break;
		combos[ch1]=combos[ch1]+nstates[ch1]*nstates[ch2];
		}	/* end addition of combinations on ch1 */
	}	/* end search through characters	*/

return combos;

}
/*
Function returning observed numbers of combinations per character.
Neads:
	states: #states per character
	notu: number of taxa
	chmatrix: character matrix
	nchars: number of characters
	UNKNOWN: value for "?"
	INAP: value for inapplicable
	
*****************************************************************************/
long **charcombinations(long **combomatrix, int *nstates, int notu, long **chmatrix, int nchars, int UNKNOWN, int INAP)
{
int 	a, b, d, ch1, ch2, sp, placer, mxst, maxst, hi, lo;
int 	*character1, *character2, *mtch, *test1, *test2, *taxst;
int 	**pair;

character1=ivector(notu);
character2=ivector(notu);
test1=ivector(notu);
test2=ivector(notu);
/*dependent=ivector(nchars);	*/

maxst=0;
for (a=0; a<nchars; ++a)	{
	for (b=0; b<notu; ++b)	{
		if (chmatrix[b][a]>maxst && (chmatrix[b][a]!=UNKNOWN && chmatrix[b][a]!=INAP))	
			maxst=chmatrix[b][a];
		}
	}

for (a=0; a<notu; ++a)
	character1[a]=character2[a]=test1[a]=test2[a]=0;

mtch=ivector(maxst+1);
taxst=ivector(maxst+1);
pair=imatrix(maxst+1,maxst+1);

for (ch1=0; ch1<nchars; ch1++)	{

	combomatrix[ch1][ch1]=0;

	for (ch2=ch1+1; ch2<nchars; ++ch2)	{
	
		combomatrix[ch1][ch2]=0;

		clearimatrix(pair,maxst+1,maxst+1,-1);

		mxst = nstates[ch1];
		for (sp=0; sp<notu; ++sp)	character1[sp]=chmatrix[sp][ch1];

		clearivector(test1,notu,RAND_MAX);
		clearivector(test2,notu,RAND_MAX);	 
		for (sp=0; sp<notu; ++sp)
			character2[sp]=chmatrix[sp][ch2];

		/*** Determine combinations ***/		
		if (nstates[ch1]>nstates[ch2])	{
			for (sp=0; sp<notu; ++sp)	{
				b = character2[sp];
				character2[sp]=character1[sp];
				character1[sp]=b;
				mxst = nstates[ch2];
				}
			hi=ch2;
			lo=ch1;
			}
		else	{
			hi=ch1;
			lo=ch2;
			}

		/*sort on placers for comparisons */
		for (sp=0; sp<notu; ++sp)	{
			placer = 0;
			/**** Sort on 1st character ****/
			for (placer=0; (character1[sp]>test1[placer] && placer< notu); placer=placer)	++placer;
			/**** Sort on 2nd character ****/
			for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]>test2[placer]) && placer<notu); placer=placer)	++placer;

			for (b=sp; b>placer; --b)	{
				test1[b]=test1[b-1];
				test2[b]=test2[b-1];
				}
			test1[placer] = character1[sp];
			test2[placer] = character2[sp];
			}
			
		/* for ease of computation, set lowest state to 0 */
		if (test1[0]>0 && (test1[0]!=INAP && test1[0]!=UNKNOWN))	{
			a=test1[0];
			for (sp=0; sp<notu; ++sp)	{
				if (test1[sp]!=INAP && test1[sp]!=UNKNOWN)	{
					test1[sp]=test1[sp]-a;
					}
				else	sp=notu;
				}
			}

		/* find where each character starts	*/
		if (test1[0]<=mxst)	sp=taxst[test1[0]]=0;
		else if (test1[0]!=INAP && test1[0]!=UNKNOWN)
			printf("Too many states for chars %d & %d on line 1173\n",ch1,ch2);
		for (d=0; d<mxst-1; ++d)	{
			taxst[d+1]=0;
			for (sp=sp; test1[sp]==d && sp<notu; sp=sp)	++sp;
			if (test1[sp]==UNKNOWN || test1[sp]==INAP)
					d=mxst;
			else	taxst[d+1]=sp;
			}

		for (a=0; a<mxst; ++a)
			mtch[a]=-1;

		mtch[a=0]=0;
		pair[0][a]=test2[0];

		/* find all realized character pairs */
		/* look at first n-1 characters; tally each state from char 2 that 
		    is paired with state a from character 1if state d+1 != state d, then 
			begin anew with the next state; repeat. */
		/* check at 17 sp., ch=0 & 20 */
		for (d=0; d<notu-1; ++d)	{
			if (test1[d]==UNKNOWN || test1[d]==INAP)	{
				d=notu;
				break;
				}
			else if (test1[d]!=test1[d+1])	{
				/* if the first state is unknown or inapplicable, then this state has no pairs */
				if ((test1[d+1]!=UNKNOWN && test1[d+1]!=INAP) && (test2[d+1]!=UNKNOWN && test2[d+1]!=INAP))	{
					++a;
					mtch[a]=b=0;
					if (a>maxst)	printf("Too many states for chars %d & %d on line 1232\n",ch1,ch2);
					pair[a][b]=test2[d+1];
					}
				}
			else if (test2[d]!=test2[d+1])	{
				if (test2[d+1]!=UNKNOWN && test2[d+1]!=INAP)	{
					++mtch[a];
					if (a>maxst)	printf("Too many states for chars %d & %d on line 1239\n",ch1,ch2);
					b=mtch[a];
					pair[a][b]=test2[d+1];
					}
				}	/* end test to see if b is new in a¥b comparison */
			}	/* end search through taxa */
			
		for (a=0; (a<=maxst /*&& pair[a][0]!=-1*/); ++a)	{
			for (a=a; pair[a][0]==-1 & a<maxst; a=a)	++a;
			if (a>maxst)	break;
			for (b=0; (b<=maxst /*&& pair[a][b]!=-1*/); ++b)	{
				for (b=b; pair[a][b]==-1; b=b)	++b;
				if (b>maxst)	break;
				if (pair[a][b]!=UNKNOWN && pair[a][b]!=INAP)	++combomatrix[ch1][ch2];
				}
			}
		combomatrix[ch2][ch1]=combomatrix[ch1][ch2];
		}	/* end search of combinations	*/
	}

free_ivector(character1);
free_ivector(character2);
free_ivector(test1);
free_ivector(test2);
free_ivector(mtch);
free_ivector(taxst);
free_imatrix(pair,maxst+1,maxst+1);

return combomatrix;
}

/*
likeallstepsgivencomp - calculates the probability of X compatibilities given Y steps and observed matrix.
Neads:
	matrix: empirical character matrix
	ctype: array giving character types (0: ordered; 1: unordered)
	nstates: array giving the number of states per character
	nchars: number of characters
	empcomp: compatibility of matrix
	bias: array giving biases in gains / loses for character, with 60 meaning P[increase]=0.6
	maxd: maximum number of changes per character
	notu: number of taxa
	fossil: 0 if no fossils, 1 if fossils included.  
	mbl: simulation parameters (origination, extinction, sampling, speciation mode)
	debug: 0 if using a random number, 1 if using a number generated by the replicate for debugging purposes
	UNKNOWN: value for "?"
	INAP: value for inapplicable
**********************************************************************************************************************************************************/
double *likeallstepsgivencomp(long **matrix, int *ctype, int *nstates, int nchars, int notu, int empcomp, int comptype, int fossils, double *mbl, int *bias, int *maxd, int pars, int debug, int UNKNOWN, int INAP)
{
int		a, b,st,max,reps,truns;
int		simcomp;
long	**simat, **tree;
double	*Lsteps;

max=iarraytotal(nstates,nchars)-nchars;
Lsteps=dvector(1+(max-pars));

printf("Enter the number of runs: ");
scanf("%s",&reps);
printf("\n");
printf("The program will examine %d trees from %d to %d steps\n",reps,pars,max);
printf("\tto find the best tree lengths for this test.\n");

printf("Doing tree #");

tree=lmatrix(notu+2,notu);		/* also will give clade diversity in first cell */
								/* finally, gives branch lengths in final two lines */
/* set the simulated matrix equal to the true matrix initially - this will make things easier later */
simat=lmatrix(notu,nchars);
equallmatrix(simat,matrix,notu,nchars);

for (truns=0; truns<reps; ++truns)	{

	if (debug==1)	srand((unsigned int) (truns));
	/* This is easier because you do not know how many nodes you'll get when sampling over time */
	if (fossils==1)	tree=evolvetreeVenn(notu, mbl, fossils);
	/* if no fossil taxa, then just make a cladogram and a Venn tree from it */
	else	{
		tree=evolvecladogram(notu,tree);
		tree=clademember(tree, notu, notu-1);
		}
	printf("%d",reps+1);
	
	for (st=pars; st<=max; ++st)	{
		a=st-pars;
		printf(" %d st\n",st);
		if (debug==1)	srand((unsigned int) (truns+st));
		simat=evolvematrix(tree, notu, simat, nchars, nstates, ctype, bias, maxd, st, UNKNOWN, INAP);
		simcomp=nu_comp(nstates, notu, simat, ctype, nchars, comptype, 0, UNKNOWN, INAP);
		if (simcomp==empcomp)	++Lsteps[a];
		printf("\b");
		for (b=1; b<=st; b=b*10) {
			printf("\b");
			}
		printf("\b\b\b\b");
		if (simcomp<(2*empcomp/3) && reps==0)	max=st-1;
		}

	for (b=(reps+1); b>=1; b=b/10)	printf("\b");
	}

free_lmatrix(simat,notu,nchars);
free_lmatrix(tree,notu+2,notu);

return Lsteps;
}

/*
probcompgivenallsteps - calculates the probability of X compatibilities given Y steps and observed matrix.
Neads:
	matrix: empirical character matrix
	ctype: array giving character types (0: ordered; 1: unordered)
	nstates: array giving the number of states per character
	nchars: number of characters
	empcomp: compatibility of matrix
	bias: array giving biases in gains / loses for character, with 60 meaning P[increase]=0.6
	maxd: maximum number of changes per character
	notu: number of taxa
	fossil: 0 if no fossils, 1 if fossils included.  
	mbl: simulation parameters (origination, extinction, sampling, speciation mode)
	debug: 0 if using a random number, 1 if using a number generated by the replicate for debugging purposes
	UNKNOWN: value for "?"
	INAP: value for inapplicable
Returns:
	PComp: each P[a][b] give the probability of a compatibilities given b steps. This is formatted for reading likelihoods.....
**********************************************************************************************************************************************************/
double **probcompgivenallsteps(long **matrix, int *ctype, int *nstates, int nchars, int notu, int empcomp, int comptype, int fossils, double *mbl, int *bias, int *maxd, int pars, int max, int debug, int UNKNOWN, int INAP)
{
int		a, b, c, s, st,mreps,treps,truns, mruns;
long	*compats;
long	**simat, **tree;
double	x;
double	**Pcomp;

a=(nchars*(nchars-1))/2;
Pcomp=dmatrix(1+(max-pars),a);

printf("Enter the number of trees to evolve: ");
scanf("%i",&treps);
printf("\n");
printf("Enter the number of matrices to evolve per tree: ");
scanf("%i",&mreps);
printf("\n");
printf("The program will examine %d trees, generating %d matrices for each tree with\n",treps,mreps);
printf("\t%d to %d steps to find the P[%d compatible pairs | steps].\n",pars,max,empcomp);
printf("Note that for each matrix, it will check compatibility after each step from %d to %d steps\n",pars,max);
printf("Doing tree #");

tree=lmatrix(notu+2,notu);		/* also will give clade diversity in first cell */
								/* finally, gives branch lengths in final two lines */
/* set the simulated matrix equal to the true matrix initially - this will make things easier later */
simat=lmatrix(notu,nchars);

for (truns=0; truns<treps; ++truns)	{

	if (debug==1)	srand((unsigned int) (truns));
	/* This is easier because you do not know how many nodes you'll get when sampling over time */
	if (fossils==1)	tree=evolvetreeVenn(notu, mbl, fossils);
	/* if no fossil taxa, then just make a cladogram and a Venn tree from it */
	else	{
		tree=evolvecladogram(notu,tree);
		tree=clademember(tree, notu, notu-1);
		}
	if ((truns+1)%10==0)	printf("%d matrix ",truns+1);
	
	for (mruns=0; mruns<mreps; ++mruns)	{

		if ((truns+1)%10==0 && (mruns+1)%5==0)	printf("%d\n",mruns+1);
		equallmatrix(simat,matrix,notu,nchars);

		compats=evolvecompat(tree, notu, simat, nchars, nstates, ctype, bias, maxd, pars, max, comptype, UNKNOWN, INAP);

		st=pars;
		while (compats[st-pars]==0)	++st;

		for (st=st; st<=max; ++st)	{
			if ((c=compats[s=st-pars])>0)
				++Pcomp[s][c];
			}
		if ((truns+1)%10==0 && (mruns+1)%5==0)	{
			for (b=(mruns+1); b>=1; b=b/10)	printf("\b");
			printf("\b");
			}
		}

	if ((truns+1)%10==0)	{
		printf("\b\b\b\b\b\b\b");
		for (b=(truns+1); b>=1; b=b/10)	printf("\b");
		printf("\b");
		}
	
	free_lvector(compats);
	}

for (s=0; s<(1+(max-pars)); ++s)	{
	x=0;
	for (c=0; c<(nchars*(nchars-1))/2; ++c)
		x+=Pcomp[s][c];
	for (c=0; c<(nchars*(nchars-1))/2; ++c)
		Pcomp[s][c]/=x;
	}

free_lmatrix(simat,notu,nchars);
free_lmatrix(tree,notu+2,notu);

return Pcomp;
}

/*
likelystepsperch - calculates the probability of X compatibilities given Y steps and observed matrix.
Requires:
	matrix: empirical character matrix
	nstates: array giving the number of states per character
	ctype: array giving character types (0: ordered; 1: unordered)
	maxd: maximum numbers of steps per character
	compat: array giving number of compatibilties per character character types (0: ordered; 1: unordered)
	nchars: number of characters
	empcomp: compatibility of matrix
	bias: array giving biases in gains / loses for character, with 60 meaning P[increase]=0.6
	maxd: maximum number of changes per character
	notu: number of taxa
	fossil: 0 if no fossils, 1 if fossils included.  
	mbl: simulation parameters (origination, extinction, sampling, speciation mode)
	debug: 0 if using a random number, 1 if using a number generated by the replicate for debugging purposes
	UNKNOWN: value for "?"
	INAP: value for inapplicable
Returns:
	Lsteps: each L[a][b] give the probability of a compatibilities given b steps and the rest of the matrix. This is formatted for reading likelihoods.....
**********************************************************************************************************************************************************/
double **likelystepsperch(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP)
{
int	a, b, c, d, ch, sp, steps, m, t, s;
int mxdl, apo;
int	truns, mruns, cruns;
int *chapos, *chuns, *chmax;
int	*dstates, *dtypes, *done;
long *chvector;
long **tree, **smatrix;
unsigned long *compat;
double x;
double **Lsteps, **psteps;
char outfile[60];
FILE *output;

strcpy(outfile,taxonname);
strcat(outfile,"_L[Steps|CharComp]");

/* get the real compatibilities for each character */
compat=char_comp(nstates,notu,matrix,ctype,nchars,comptype,0,UNKNOWN,INAP);

/* see if there are dependent inapplicables */
a=0;
for (ch=0; ch<nchars && a==0; ++ch)	if (depend[ch]!=ch)	a=1;

/* allocate memory for "dummy" arrays & matrices - these simply jackknife the real data to make a nchars-1 matrix & character info */
dstates=ivector(nchars+1);
dtypes=ivector(nchars+1);
done=ivector(nchars);
chapos=autapomorphies(matrix,nstates,notu,nchars,UNKNOWN,INAP);
chuns=unknownstates(matrix,notu,nchars,UNKNOWN,INAP);
chmax=ivector(nchars);
smatrix=lmatrix(notu, nchars+1);

chvector=lvector(notu);

mxdl=maxiarray(maxd,nchars);	/* maximum number of changes possible for any character */
if (mxdl>notu)	mxdl=notu-1;
Lsteps=dmatrix(nchars,mxdl+1);

clearivector(chmax,nchars,notu-1);
/*for (ch=0; ch<nchars; ++ch)
/*	if (nstates[ch]>=4 && chmax[ch]<(notu-1))	chmax[ch]=notu-1;	*/

equalivector(dstates,nstates,nchars);
equalivector(dtypes,ctype,nchars);

printf("\nFor each character, this routine will generate X trees and Y matrices per tree.\n");
printf("  Each matrix is of N-1 characters, with N being the observed number.  The compatibility.\n");
printf("  of each matrix will be the same as the compatibility of the other N-1 characters in\n");
printf("  the observed matrix.  The program then will evolve the character K, K+1, K+2É steps Z times\n");
printf("  (with K being one less than the number of states).  The program then calculate the\n");
printf("  compatibility of the character.  It also will determine whether the proper number of states\n");
printf("  and the proper number of derived taxa evolve.  Over X*Y*Z runs, this will determine the\n");
printf("  probability of the observed structure of each character given K, K+1, etc., steps.  \n\n");

printf("Enter the number of trees to use for each character (X above): ");
scanf("%i",&truns);
printf("\n");
printf("Enter the number of matrices to evolve per tree (Y above): ");
scanf("%i",&mruns);
printf("\n");
printf("Enter the number of times to replicate each number of steps (Z above): ");
scanf("%i",&cruns);
printf("\n");
printf("The program will examine %d trees, generating %d matrices for each tree with\n",truns,mruns);
printf("\teach number of steps replicated %d times to find L[setps | CPs, states, apomorphies].\n",cruns);

tree=lmatrix(notu+2,notu);			/* also will give clade diversity in first cell 			*/
									/* 		finally, gives branch lengths in final two lines 	*/
psteps=dmatrix(mxdl+1,nchars*notu);	/* matrix that is maxsteps X (nchars*notu); emulates a cube rather than a matrix	*/

printf("Doing tree ");

for (t=0; t<truns; ++t)	{
	/* evolve a tree to evolve characters over */
	if (debug==1)	srand((unsigned int) t);
	/* This is easier because you do not know how many nodes you'll get when sampling over time */
	if (fossil==1)	tree=evolvetreeVenn(notu, mbl, fossil);
	/* if no fossil taxa, then just make a cladogram and a Venn tree from it */

	else	{
		tree=evolvecladogram(notu,tree);
		tree=clademember(tree, notu, notu-1);
		}
	printf("%d, matrix ",t+1);
	for (m=0; m<mruns; ++m)	{

		if (debug==1)	srand((unsigned int) (t*mruns)+m);

		/* evolve a matrix matching the compatibility of the original matrix	*/
		smatrix=evolvetocompat(tree,empcompat,notu,matrix,nchars,nstates,ctype,bias,chmax,depend,comptype,UNKNOWN,INAP);
		
		clearivector(done,nchars,0);
		for (ch=0; ch<nchars; ++ch)	{
			while (done[ch]==1 && ch<nchars)	++ch;
			if (ch>=nchars)	break;
			
			dtypes[nchars]=ctype[ch];
			dstates[nchars]=nstates[ch];
			for (c=0; c<cruns; ++c)	{
				if (debug==1)	srand((unsigned int) (t*mruns*cruns)+(m*cruns)+c);
/*				fprintf(debugout,"tree %d matrix %d character run %d character %d\n",t,m,c,ch);
/*				fflush(stdout);	*/
				for (steps=nstates[ch]-1; steps<=mxdl; ++steps)	{
					chvector=evolvecharacterNsteps(tree,steps,notu,chvector,nstates[ch],ctype[ch],bias[ch],UNKNOWN,INAP);
					/* next, add it to the rest of the matrix */
					for (sp=0; sp<notu; ++sp)	smatrix[sp][nchars]=chvector[sp];
					/* find its compatibility */
					d=char_nu_comp(nchars, dstates, notu, smatrix, dtypes, nchars, comptype, 0, UNKNOWN, INAP);
					apo=autapo_char(smatrix,notu,c,UNKNOWN,INAP);
					/* increment the probability of the observed combination of compatibilities and taxa with derived condition by 1/cruns	*/
					x=readdcube(psteps,steps,apo,d,nchars);			/*	read value from */
					assigndcube(psteps,steps,apo,d,nchars,x+1);		/*	*/
					clearlvector(chvector,notu,0);					/*	*/
					}
				}
			
			
			for (steps=nstates[ch]-1; steps<=mxdl; ++steps)
				Lsteps[ch][steps]=(readdcube(psteps,steps,chapos[ch],compat[ch],nchars)+readdcube(psteps,steps,chapos[ch],compat[ch]+1,nchars))/2;
			done[ch]=1;
			for (c=ch+1; c<nchars; ++c)	{
				while (done[c]==1 && c<nchars)	++c;
				if (chuns[c]==chuns[ch] && (nstates[c]==nstates[ch] && ctype[c]==ctype[ch]))	{
					for (steps=nstates[ch]-1; steps<=mxdl; ++steps)
						Lsteps[c][steps]=(readdcube(psteps,steps,chapos[c],compat[c],nchars)+readdcube(psteps,steps,chapos[c],compat[c]+1,nchars))/2;
					done[c]=1;
					}	/* tally likelihoods for similar characters	*/
				}	/* search for other characters with the same number of states and unknowns	*/
			}	/* try the next character now	*/

		printf("%d\n",m+1);
		if (m<9)		 			printf("\b\b");
		else if (m>8 && m<99) 		printf("\b\b\b");
		else if (m>98 && m<999)		printf("\b\b\b\b");
		}	/* end analysis with this matrix	*/
	printf("\b\b\b\b\b\b\b\b");
	for (b=(t+1); b>=1; b=b/10)	printf("\b");
	printf("\b");
	}	/* end simulation of trees	*/


x=truns*mruns*cruns;
for (ch=0; ch<nchars; ++ch)	{
	for (steps=nstates[ch]-1; steps<mxdl; ++steps)	{
		Lsteps[ch][steps]/=x;	
		}
	}

output=fopen(outfile,"w");
fprintf(output,"Steps");
for (ch=0; ch<nchars; ++ch)	{
	if (ch<9)	{
		if (nchars>9 && nchars<100)	fprintf(output,"\tChar_0%d",ch+1);
		else if (nchars>99)			fprintf(output,"\tChar_00%d",ch+1);
		else						fprintf(output,"\tChar_%d",ch+1);
		}
	else	{
		if (nchars>99)				fprintf(output,"\tChar_0%d",ch+1);
		else						fprintf(output,"\tChar_%d",ch+1);
		}
	}
fprintf(output,"\n");

m=maxiarray(maxd,nchars);
for (s=1; s<m; ++s)	{
	fprintf(output,"%d",s);
	for (ch=0; ch<nchars; ++ch)	{
		if (s>=(nstates[ch]-1))	fprintf(output,"\t%9.8f",Lsteps[ch][s]);
		else					fprintf(output,"\t¥");
		}
	fprintf(output,"\n");
	} 
fclose(output);

free_ivector(dstates);
free_ivector(dtypes);
free_ivector(done);
free_ivector(chapos);
free_ivector(chuns);
free_lmatrix(smatrix,notu, nchars+1);
free_lmatrix(tree,notu+2,notu);		/* also will give clade diversity in first cell */
free_dmatrix(psteps,mxdl+1,nchars*notu);

return Lsteps;
}


/*likelystepsperch - calculates the probability of X compatibilities given Y steps and observed matrix.
Requires:
	matrix: empirical character matrix
	nstates: array giving the number of states per character
	ctype: array giving character types (0: ordered; 1: unordered)
	maxd: maximum numbers of steps per character
	compat: array giving number of compatibilties per character character types (0: ordered; 1: unordered)
	nchars: number of characters
	empcomp: compatibility of matrix
	bias: array giving biases in gains / loses for character, with 60 meaning P[increase]=0.6
	maxd: maximum number of changes per character
	notu: number of taxa
	fossil: 0 if no fossils, 1 if fossils included.  
	mbl: simulation parameters (origination, extinction, sampling, speciation mode)
	debug: 0 if using a random number, 1 if using a number generated by the replicate for debugging purposes
	UNKNOWN: value for "?"
	INAP: value for inapplicable
Returns:
	Lsteps: each L[a][b] give the probability of a compatibilities given b steps and the rest of the matrix. This is formatted for reading likelihoods.....
**********************************************************************************************************************************************************/
double **likelystepsperchar(long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP)
{
int	a, b, c, d, f, g, ch, sp, steps, m, t;
int dcomp, mxst, mxdl, redch, nxt, nxtcmpa, nxtcmpb, nxtcmpc, nxtok;
int	truns, mruns, cruns;
int	apo, dapo, unsc, logdp, mxdp;
int	*dstates, *dtypes, *dmaxd, *dbias, *ddpnd, *sychos;
long *chvector, *nxtch, **chmat;
long **tree, **dsmatrix;
/*double runs;	*/
unsigned long *compat, *dcompat;
double **Lsteps, **denom;

/* get the real compatibilities for each character */
compat=char_comp(nstates,notu,matrix,ctype,nchars,comptype,0,UNKNOWN,INAP);

/* see if there are dependent inapplicables */
a=0;
for (ch=0; ch<nchars && a==0; ++ch)	if (depend[ch]!=ch)	a=1;
sychos=ivector(nchars);				/* sychos[x] gives the number of characters depending of character x */
if (a==1)
	for (ch=0; ch<nchars && a==0; ++ch)	++sychos[depend[ch]]; 
else
	clearivector(sychos,nchars,1);

/* allocate memory for "dummy" arrays & matrices - these simply jackknife the real data to make a nchars-1 matrix & character info */
dstates=ivector(nchars);
dtypes=ivector(nchars);
dmaxd=ivector(nchars-1);
dbias=ivector(nchars-1);
ddpnd=ivector(nchars-1);
dsmatrix=lmatrix(notu, nchars);

chvector=lvector(notu);
nxtch=lvector(notu);
mxdp=maxiarray(sychos,nchars);
if (mxdp>1)
	chmat=lmatrix(notu,mxdp);

mxdl=maxiarray(maxd,nchars);	/* maximum number of changes possible for any character */
if (mxdl>notu)	mxdl=notu-1;
Lsteps=dmatrix(nchars,3*(mxdl+1));
denom=dmatrix(nchars,mxdl+1);

printf("\nFor each character, this routine will generate X trees and Y matrices per tree.\n");
printf("  Each matrix is of N-1 characters, with N being the observed number.  The compatibility.\n");
printf("  of each matrix will be the same as the compatibility of the other N-1 characters in\n");
printf("  the observed matrix.  The program then will evolve the character K, K+1, K+2É steps Z times\n");
printf("  (with K being one less than the number of states).  The program then calculate the\n");
printf("  compatibility of the character.  It also will determine whether the proper number of states\n");
printf("  and the proper number of derived taxa evolve.  Over X*Y*Z runs, this will determine the\n");
printf("  probability of the observed structure of each character given K, K+1, etc., steps.  \n\n");

printf("Enter the number of trees to use for each character (X above): ");
scanf("%i",&truns);
printf("\n");
printf("Enter the number of matrices to evolve per tree (Y above): ");
scanf("%i",&mruns);
printf("\n");
printf("Enter the number of times to replicate each number of steps (Z above): ");
scanf("%i",&cruns);
printf("\n");
printf("The program will examine %d trees, generating %d matrices for each tree with\n",truns,mruns);
printf("\teach number of steps replicated %d times to find L[setps | CPs, states, apomorphies].\n",cruns);

tree=lmatrix(notu+2,notu);		/* also will give clade diversity in first cell */
								/* finally, gives branch lengths in final two lines */

printf("Doing tree ");
for (t=0; t<truns; ++t)	{
	/* evolve a tree to evolve characters over */
	if (debug==1)	srand((unsigned int) t);
	/* This is easier because you do not know how many nodes you'll get when sampling over time */
	if (fossil==1)	tree=evolvetreeVenn(notu, mbl, fossil);
	/* if no fossil taxa, then just make a cladogram and a Venn tree from it */

	else	{
		tree=evolvecladogram(notu,tree);
		tree=clademember(tree, notu, notu-1);
		}
	printf("%d, matrix ",t+1);
	for (m=0; m<mruns; ++m)	{
		printf("%d, character ",m+1);
		for (ch=0; ch<nchars; ++ch)	{

			nxt=0;					/* this becomes 1 when we 
			printf("%d\n",ch+1);
			/* rewrite array of dependents to take into account that the characters are migrating */


			/* count the unknowns and inapplicables */
			unsc=0;
			for (sp=0; sp<notu; ++sp)	if (matrix[sp][ch]==UNKNOWN || matrix[sp][ch]==INAP)	++unsc;
			
			/* if we are on the first character, then evolve the remaining characters; we'll replace each one as we go	*/
			if (ch==0 || nxt==0)	{
				redch=0;	/* reduced number of characters */
				logdp=0;

				for (c=0; c<nchars; ++c)	{
					/************************************************************************************/
					/* logdp counts the logical dependents; it will be at least 1 when we reach ch		*/
					/*    if there are multiple dependents, then it will be >1;	 The reduced matrix		*/
					/*    removed all characters dependent on ch (including ch); thus, the properies of	*/
					/*    the remaining characters must be moved up; the logical dependence must be		*/
					/*    changed because the character numbers change.  For example, if chars. 4 & 5	*/
					/*    depend on char. 3 and we remove char. 1, then char. 3 becomes char. 2 and 	*/
					/*	  chars. 4&5 become 3&4 and dependent on 2.  Note that when we remove 3, chars	*/
					/*    1 & 2 stay put, but char. 6 is moved up to char. 3 because 4 & 5 also are		*/
					/*    removed.  Those will be evolved along with char. 3 after the matrix evolves.  */
					/************************************************************************************/
					if (depend[c]==ch)	++logdp;
					else	{
						dstates[c-logdp]=nstates[c];
						dtypes[c-logdp]=ctype[c];
						dmaxd[c-logdp]=maxd[c];
						dbias[c-logdp]=bias[c];
						ddpnd[c-logdp]=depend[c]-logdp;
						for (sp=0; sp<notu; ++sp)	dsmatrix[sp][c-logdp]=matrix[sp][c];
						++redch;
						}
					}

				/* clear vector for character change */
				if (logdp==1)	for (sp=0; sp<notu; ++sp)	chvector[sp]=matrix[sp][ch];

				dstates[nchars-1]=nstates[ch];
				dtypes[nchars-1]=ctype[ch];

				/* create reduced matrix and character info arrays excluding character ch and any characters logically dependent upon it*/
				apo=autapo_char(matrix,notu,ch,UNKNOWN,INAP);
				dcomp=empcompat-compat[ch];
				for (c=0; c<nchars; ++c)	{
					/* remove compatibilities of logically dependent characters */
					if (depend[c]==ch && c!=ch)
						dcomp=dcomp-compat[c];
					}

				/* if character is apomorphic, then skip this part - it has a perfect compatibility regardless */
				dsmatrix=evolvetocompat(tree,dcomp,notu,dsmatrix,redch,dstates,dtypes,dbias,dmaxd,ddpnd,comptype,UNKNOWN,INAP);
				dcompat=char_comp(dstates,notu,dsmatrix,dtypes,redch,comptype,0,UNKNOWN,INAP);
				
				}
			/* otherwise, put the last character into the slot where this character was; note, however, that it will be scored	*/
			/* so that the remaining compatibility matches that of the rest of the matrix										*/
			else	{
				dstates[ch]=nstates[ch];
				dtypes[ch]=ctype[ch];
				dmaxd[ch]=maxd[ch];
				dbias[ch]=bias[ch];
				ddpnd[ch]=depend[ch]-logdp;
				for (sp=0; sp<notu; ++sp)	dsmatrix[sp][ch]=nxtch[sp];
				}

			/* this is easier if the simulated character has compatibility similar to the real character 		*/
			/* 	so, try to find another that does (with same ordering and number of states, and flip-flop them.	*/
			if (ch<nchars-1 && dcompat[ch]!=compat[ch+1])	{
				nxtok=0;
				for (c=ch+1; c<nchars-1; ++c)	{
					if (dcompat[c]==compat[ch+1] && (nstates[ch+1]==dstates[c] && ctype[ch+1]==dtypes[c]))	{
						for (sp=0; sp<notu; ++sp)	{
							chvector[sp]=dsmatrix[sp][c];
							dsmatrix[sp][c]=dsmatrix[sp][ch];
							dsmatrix[sp][ch]=chvector[sp];
							}
						d=dmaxd[c];
						dmaxd[c]=dmaxd[ch];
						dmaxd[ch]=d;
						
						d=dcompat[c];
						dcompat[c]=dcompat[ch];
						dcompat[ch]=d;
						
						c=nchars;
						
						nxtok=1;
						}	/* end the finding and swapping with a similar character	*/
					}	/* search remaining characters for one looks like the real character	*/
				}	/* end case where simulated character looks different from real one	*/

			if (ch<nchars-1)	{
				/* find compatibility among the N-2 characters other than ch & ch+1 - this will help us determine how compatible ch */
				/*		needs to be when evaluating ch+1 coming up; this way, we don't have to re-evolve the matrix!				*/
/*				nxtcomp=0;
				for (c=0; c<nchars-1; ++c)	{
					if (c!=ch && c!=ch+1)	{
						for (d=c+1; d<nchars; ++d)	{
							if (d!=ch && d!=ch+1)	nxtcomp+=comatrix[c][d];
							}
						}
					}	*/
				nxtcmpa=empcompat-dcompat[ch];	/* this gives compatibility among the N-2 other characters	*/
				nxtcmpb=empcompat-compat[ch+1];	/* this is where we need to be								*/
				nxtcmpc=empcompat-compat[ch];	/* this is where we need to be								*/
				}
			
			if (debug==1)	srand((unsigned int) c+pow((t+2),(ch+1)));
		/*		printf("%d matrix ",t+1);	*/

			if (apo>1 || (apo+unsc)<notu-1)	{

				printf("%d\n",m+1);
				/* now evolve the reduced matrix until it matches the compatibility of the other nchar-1 characters 		*/
				/* What we will do is seed the random number generator with a particular value and retain that value;		*/
				/* We then use the evolvecompat routine to generate an array of compatibilities per step; if the desired	*/
				/* compatibility is in that array, then we will use evolvematrix to evolve that matrix AFTER reseeding the	*/
				/* random number generator with the same seed with the number of steps known to produce the desired 		*/
				/* compatibility.  This means that we need to keep track of our random numbers!								*/
				
				/* now, find the probability of getting observed compatibility AND observed apomorphy distribution given rates */
					/* easy routine for characters with no dependents */
				if (logdp==1)	{
					for (steps=nstates[ch]-1; steps<maxd[ch]+4 && steps<(notu-1); ++steps)	{
						for (c=0; c<cruns; ++c)	{
							/* now, evolve character */
							chvector=evolvecharacterNsteps(tree,steps,notu,chvector,nstates[ch],ctype[ch],bias[ch],UNKNOWN,INAP);
							/* next, add it to the rest of the matrix */
							for (sp=0; sp<notu; ++sp)	dsmatrix[sp][nchars-1]=chvector[sp];
							/* find its compatibility */
							d=char_nu_comp(nchars-1, dstates, notu, dsmatrix, dtypes, nchars, comptype, 0, UNKNOWN, INAP);
							/* now see if compatibilities match */
							if (nstates[ch]>2 && ctype[ch]==0)	{
								mxst=0;
								for (sp=0; sp<notu; ++sp)
									if (chvector[sp]!=UNKNOWN && chvector[sp]!=INAP)
										if (chvector[sp]>mxst)	mxst=chvector[sp];
								}
							else	mxst=nstates[ch]-1;

							if (d==compat[ch] && mxst==(nstates[ch]-1))	{
								++Lsteps[ch][steps];
								dapo=autapo_char(dsmatrix,notu,nchars-1,UNKNOWN,INAP);	/* check this: should it be nchars-1 or ch?	*/

								if (dapo==apo || (notu-dapo)==apo)	++Lsteps[ch][notu+steps];
								}	/* end examination of whether we duplicated compatibility */
							if ((nxt==0 && nxtok==1) && (d>=compat[ch]-1 && d<=compat[ch]+1))	{
								f=pair_comp(ch, nchars-1, dstates, notu, dsmatrix, dtypes, nchars, comptype, 0, UNKNOWN, INAP);
								g=pair_comp(ch, ch+1, nstates, notu, matrix, ctype, nchars, comptype, 0, UNKNOWN, INAP);
								if (f==g && d==compat[ch])	{
									equallvector(nxtch,chvector,sp);
									nxt=1;
									}
								else if ((f==1 && g==0) && d==compat[ch]+1)	{
									equallvector(nxtch,chvector,sp);
									nxt=1;
									}
								else if ((f==0 && g==1) && d==compat[ch]-1)	{
									equallvector(nxtch,chvector,sp);
									nxt=1;
									}
								}
/*							if (nxt==0 && (d+dcomp-dcompat[ch])==nxtcmpb || (d+dcomp-dcompat[ch])==nxtcmpb)	{
/*								twostates[0]=nstates[ch];
/*								twostates[1]=dstates[ch];
/*								twotypes[0]=ctypes[ch];
/*								twotypes[1]=dtypes[ch];
/*								pair_comp(ch, nchars-1, int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int outgroup, int UNKNOWN, int INAP);
/*								f=pair_comp(ch, nchars-1, dstates, notu, dsmatrix, dtypes, nchars, comptype, UNKNOWN, INAP);
/*								if ((d==nxtcmpa+1 && f==1) || (d==nxtcmpa && f==0))	{
/*									equallvector(nxtch,chvector,sp);
/*									nxt=1;
/*									}
/*								}	*/
							++denom[ch][steps];	/* increment number of runs */
							}	/* end run of possible steps */
						}	/* end case where character has no logical dependents */
					}

				/* more difficult routine for those with depdendents */
				else	{
					for (steps=nstates[ch]-1; steps<maxd[ch]+4 && steps<(notu-1); ++steps)	{
						/* now, evolve character */
						chvector=evolvecharacterNsteps(tree,steps,notu,chvector,nstates[ch],ctype[ch],bias[ch],UNKNOWN,INAP);
						/* next, add it to the rest of the matrix */
						for (sp=0; sp<notu; ++sp)	dsmatrix[sp][nchars-logdp]=chvector[sp];
						/* find the characters compatibility compatibility 						*/
						/* note: you need to omit the dependent characters in two ways; 		*/
						/*   1) simulated charact is #chars-dependents (not #chars-1)			*/
						/*   2) because all dependents are compatible with their independent,	*/
						/*      we want to look for #compatibilites - dependents				*/
						d=char_nu_comp(nchars-logdp, dstates, notu, dsmatrix, dtypes, (nchars-logdp+1), comptype, 0, UNKNOWN, INAP);
						
						/* now see if compatibilities match */
						/* first, determine the maximum number of st */
						if (nstates[ch]>2)	{
							mxst=0;
							for (sp=0; sp<notu; ++sp)
								if (chvector[sp]!=UNKNOWN && chvector[sp]!=INAP)
									if (chvector[sp]>mxst)	mxst=chvector[sp];
							}
						else	mxst=nstates[ch]-1;

						if (d==compat[ch]-logdp && mxst==(nstates[ch]-1))	{
							++Lsteps[ch][steps];
							dapo=autapo_char(dsmatrix,notu,nchars-1,UNKNOWN,INAP);	/* check this: should it be nchars-1 or ch?	*/

							/* if we replicate the number of taxa with derived conditions, continue */
							if (dapo==apo || (notu-dapo)==apo)	{
								++Lsteps[ch][mxdl+steps+1];
								/* if we have the right number of apomorphies, then evolve the dependent characters */
								if (dapo==apo)	{
									/* the minimum number of steps could be 0 IF the number of steps for the independent matches 	*/
									/* the number of states for the dependent character												*/
/*										for ()	*/
									}
								}	/* tally cases where everything matches */
							}	/* end examination of whether we duplicated compatibility */
						
						if (dapo==apo || (notu-dapo)==apo)
							++Lsteps[ch][2*(mxdl+1)+steps];
						
						}	/* end run of possible steps */
					}	/* end case where character has logical dependents */
				for (b=(m+1); b>=1; b=b/10)	printf("\b");
				printf("\b");
				}	/* end case for non-apomorphic taxa */
			/* if apomorphic, then just determine the probability of apomorphy given K, K+1, K+2, etc. steps. */
			else	{
				for (steps=nstates[ch]-1; steps<maxd[ch]+3 && steps<notu; ++steps)	{
					chvector=evolvecharacterNsteps(tree,steps,notu,chvector,nstates[ch],ctype[ch],bias[ch],UNKNOWN,INAP);
					dapo=0;
					/* if only one taxon has '1' or if all but one have '1' then it is apomorphic */
					for (sp=0; sp<notu; ++sp)	if (chvector[ch]==1)	++dapo;
					if (dapo==1 || dapo==(notu-unsc-1))	{
						++Lsteps[ch][steps];
						++Lsteps[ch][notu+steps];
						}	/* tally apomorphic successes */
					}	/* end X steps */
				}	/* end routine for apomorphic states */
			}	/* end evolved trees */
		}
	dstates[ch]=nstates[ch];
	dtypes[ch]=dtypes[ch];
	dmaxd[ch]=maxd[ch];

	printf("\b\b\b\b\b\b\b");
	for (b=(t+1); b>=1; b=b/10)	printf("\b");
	printf("\b");
	}
}


/*
statepairs - counts the number of state pairs in a matrix.
Requires:
	matrix: empirical character matrix
	notu: number of taxa
	nchar: number of characters
	maxstate: maximum number of states for any character
	UNKNOWN: value for "?"
	INAP: value for inapplicable
Returns:
	pairs: Number of state pairs.....
**********************************************************************************************************************************************************/
int statepairs(long **mat, int notu, int nchar, int maxstate/*, int UNKNOWN, int INAP*/)
{
int	ch1, ch2, s1, s2, sp;
int	a, b;
int pairs=0;
int	**prmt;
/*int	prmt[50][50];
/*for (s1=0; s1<50; ++s1)	for (s2=0; s2<50; ++s2)	prmt[s1][s2]=0;	*/
prmt=imatrix(nchar*(maxstate+1),nchar*(maxstate+1));

for (sp=0; sp<notu; ++sp)	{
	for (ch1=0; ch1<nchar-1; ++ch1)	{
		s1=mat[sp][ch1];
		a=(ch1*maxstate)+s1;
		for (ch2=ch1+1; ch2<nchar; ++ch2)	{
			s2=mat[sp][ch2];
			b=(ch2*maxstate)+s2;
			
			if (prmt[a][b]==0)	{
				++prmt[a][b];
				++pairs;
				}
			}
		}
	}

free_imatrix(prmt,nchar*(maxstate+1),nchar*(maxstate+1));

return pairs;
}

/* obspairs - lists the pairs of states observed together on two character vectors.
Requires:
	t1: 1st character vector
	t2: 2nd character vector
	combos: matrix giving the number of combinations for each state of character 1
	st: the number of states for character 1 (vector t1)
	notu: number of taxa (length of vector t1)
Returns:
	pairs: observed pairings of characters 1 and 2.....
**********************************************************************************************************************************************************/
int **obsstatepairs(int *t1, int *t2, int *combos, int st, int notu)
{
int	mx, st1=0, st2=0, sp;
int **pairs;

mx=maxiarray(combos,st);
pairs=imatrix(st,mx);

pairs[0][0]=t2[0];
for (sp=1; sp<notu; ++sp)	{
	if (t1[sp]!=t1[sp-1])
		pairs[++st1][st2=0]=t2[sp];

	else if (t2[sp]!=t2[sp-1])	pairs[st1][++st2]=t2[sp];
	}

return pairs;
}


/*
likelystepspercharacter - calculates the probability of X compatibilities given Y steps and observed matrix.
Requires:
	taxonname: to be used for output files
	matrix: empirical character matrix
	nstates: array giving the number of states per character
	ctype: array giving character types (0: ordered; 1: unordered)
	maxd: maximum numbers of steps per character
	compat: array giving number of compatibilties per character
	ctype: array giving character types (0: ordered; 1: unordered)
	nchars: number of characters
	empcomp: compatibility of matrix
	bias: array giving biases in gains / loses for character, with 60 meaning P[increase]=0.6
	depend: depend[x] gives on which character that character x depends; if char 50 is feather and char 51 is feather color, then depend[50]=50, depend[51]=50
	maxd: maximum number of changes per character
	notu: number of taxa
	fossil: 0 if no fossils, 1 if fossils included.  
	mbl: simulation parameters (origination, extinction, sampling, speciation mode)
	debug: 0 if using a random number, 1 if using a number generated by the replicate for debugging purposes
	UNKNOWN: value for "?"
	INAP: value for inapplicable
Returns:
	Lsteps: each L[a][b] give the probability of a compatibilities given b steps and the rest of the matrix. This is formatted for reading likelihoods.....
**********************************************************************************************************************************************************/
void likelystepspercharacter(char *taxonname, long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP)
{
int	a, b, c, d, ch, m, n, t, s;
int mxdl, apo, combos, smcm;
int	truns, mruns;
int *chapos, *chuns, *chmax, *chstcmb, *chsteps;
/*int	*dstates, *dtypes, *done;	*/
long **tree, **smatrix, **chstcls;
unsigned long *compat;
double x, y, mr;
double **Lsteps, **psteps, **examples, **lnL;
char outfile1[60], outfile2[60], outfile3[60], outfile4[60];
FILE *output1, *output2, *output3, *output4;

/* get the real compatibilities for each character */
compat=char_comp(nstates,notu,matrix,ctype,nchars,comptype,0,UNKNOWN,INAP);

/* see if there are dependent inapplicables */
a=0;
for (ch=0; ch<nchars && a==0; ++ch)	if (depend[ch]!=ch)	a=1;

/* allocate memory for "dummy" arrays & matrices - these simply jackknife the real data to make a nchars-1 matrix & character info */
/*dstates=ivector(nchars+1);
/*dtypes=ivector(nchars+1);
/*done=ivector(nchars);	*/
chapos=autapomorphies(matrix,nstates,notu,nchars,UNKNOWN,INAP);		/* number of otus with "minority" condition		*/
chuns=unknownstates(matrix,notu,nchars,UNKNOWN,INAP);				/* number of unknown conditions per character	*/
chstcls=charstateunkncombos(nstates,chuns,nchars);				/* combinations of states and unknowns observed	*/
chmax=ivector(nchars);
chsteps=ivector(nchars);
smatrix=lmatrix(notu, nchars+1);									/* matrix to be simulated						*/

chstcmb=ivector(nchars);											/* state+unknown combination# for each character*/

for (combos=0; chstcls[combos][0]>-1; combos=combos)	++combos;

for (ch=0; ch<nchars; ++ch)	{
	for (c=0; c<combos; ++c)	{
		/* if character ch has chstcls[c][0] states and chstcls[c][1] unknowns, then it is combo # c */
		if (chstcls[c][0]==nstates[ch] && chstcls[c][1]==chuns[ch])	{
			chstcmb[ch]=c;
			c=combos;
			}	/* terminate search	*/
		}	/* look through all combos	*/
	}	/* find proper combo for each character	*/

mxdl=maxiarray(maxd,nchars);	/* maximum number of changes possible for any character */
if (mxdl>notu)	mxdl=notu-1;	/* do not bother making it greater than the # of otus */
Lsteps=dmatrix(nchars,mxdl+1);

clearivector(chmax,nchars,notu-1);
/*for (ch=0; ch<nchars; ++ch)
/*	if (nstates[ch]>=4 && chmax[ch]<(notu-1))	chmax[ch]=notu-1;	

/*equalivector(dstates,nstates,nchars);
/*equalivector(dtypes,ctype,nchars);	*/

printf("\nFor each character, this routine will generate X trees and Y matrices per tree.\n");
printf("\nThis routine will generate X trees and Y matrices per tree.  The number of steps\n");
printf("  the compatibility, and the number of ÒminorityÓ conditions for each character.\n");
printf("  are tallied when the whole simulated matrix has the same compatibility as the.\n");
printf("  empirical matrix.  This is repeated, generating:\n");
printf("  P[Compatibility, Distribution | Steps, Overall Compatibility]\n");
printf("  = L[Steps, Overall Compatibility | Compatibility, Distribution].\n");

printf("Enter the number of trees to use (X above): ");
scanf("%i",&truns);
printf("\n");
printf("Enter the number of matrices to evolve (Y above): ");
scanf("%i",&mruns);
printf("\n");
printf("The program will examine %d trees, generating %d matrices for each tree.\n",truns,mruns);

for (combos=0; chstcls[combos][0]>-1; ++combos)	combos=combos;

tree=lmatrix(notu+2,notu);			/* also will give clade diversity in first cell 			*/
									/* 		finally, gives branch lengths in final two lines 	*/
psteps=dmatrix(mxdl+1,combos*notu);	/* matrix that is maxsteps X (nchars*notu); emulates a cube rather than a matrix	*/
examples=dmatrix(combos,mxdl);		/* this will give the number of times that a simulation has X steps given the combination of Y States & Z unknowns	*/

printf("Doing tree ");

for (t=0; t<truns; ++t)	{
	/* evolve a tree to evolve characters over */
	if (debug==1)	srand((unsigned int) t);
	/* This is easier because you do not know how many nodes you'll get when sampling over time */
	if (fossil==1)	tree=evolvetreeVenn(notu, mbl, fossil);
	/* if no fossil taxa, then just make a cladogram and a Venn tree from it */

	else	{
		tree=evolvecladogram(notu,tree);
		tree=clademember(tree, notu, notu-1);
		}
	printf("%d, matrix ",t+1);
	for (m=0; m<mruns; ++m)	{

		if (debug==1)	srand((unsigned int) (t*mruns)+m);

		/* evolve a matrix matching the compatibility of the original matrix	*/
		clearivector(chsteps,nchars,0);
		smatrix=evolvetocompatsavessteps(tree,empcompat,notu,matrix,nchars,nstates,ctype,bias,chmax,depend,chsteps,comptype,UNKNOWN,INAP);
		for (ch=0; ch<nchars; ++ch)	{
			smcm=char_nu_comp(ch, nstates, notu, smatrix, ctype, nchars, comptype, 0, UNKNOWN, INAP);	/* get the compatibility of character ch									*/
			apo=autapo_char(smatrix,notu,ch,UNKNOWN,INAP);												/*  get "minority" number (=1 for autapomorphy or only outgroup primitive)	*/
			x=readdcube(psteps,chsteps[ch],apo,smcm,combos);											/*	read value from psteps 													*/
			assigndcube(psteps,chsteps[ch],apo,smcm,combos,x+1);										/*	find the appropriate cell number for the "3D" matrix					*/
			++examples[chstcmb[ch]][chsteps[ch]];														/* another example of stateunkn combo X with Y steps						*/
			}
		printf("%d\n",m+1);
		if (m<9)		 			printf("\b\b");
		else if (m>8 && m<99) 		printf("\b\b\b");
		else if (m>98 && m<999)		printf("\b\b\b\b");
		}	/* end analysis with this matrix	*/
	printf("\b\b\b\b\b\b\b\b");
	for (b=(t+1); b>=1; b=b/10)	printf("\b");
	printf("\b");
	}

/* tally likelihoods from finished results	*/	
for (ch=0; ch<nchars; ++ch)	{
	c=compat[ch];								/* get the real compatibility						*/
	a=chapos[ch];								/* get the real "apomorphies"						*/
	b=chstcmb[ch];								/* get the real combination of states and unknowns	*/
	for (s=1; s<mxdl; ++s)	{
		x=readdcube(psteps,s,a,c,b);			/*	# times obs. compatibility & apomorphies seen given steps, unknowns									*/
		y=examples[b][s];						/*	# times given steps, unknowns happened for this states+unknowns combo								*/
		Lsteps[ch][s]=x/y;						/*	p[compatibility, apos | steps, states, unknowns] = L[steps, states, unknowns | compatibility, apos	*/  
		}
	}

strcpy(outfile1,taxonname);
strcat(outfile1,"_L[Steps|CharComp].xls");
output1=fopen(outfile1,"w");
fprintf(output1,"Steps");
for (ch=0; ch<nchars; ++ch)	{
	if (ch<9)	{
		if (nchars>9 && nchars<100)	fprintf(output1,"\tChar_0%d",ch+1);
		else if (nchars>99)			fprintf(output1,"\tChar_00%d",ch+1);
		else						fprintf(output1,"\tChar_%d",ch+1);
		}
	else	{
		if (nchars>99)				fprintf(output1,"\tChar_0%d",ch+1);
		else						fprintf(output1,"\tChar_%d",ch+1);
		}
	}
fprintf(output1,"\n");

m=maxiarray(maxd,nchars);
for (s=1; s<m; ++s)	{
	fprintf(output1,"%d",s);
	for (ch=0; ch<nchars; ++ch)	{
		if (s>=(nstates[ch]-1))	fprintf(output1,"\t%9.8f",Lsteps[ch][s]);
		else					fprintf(output1,"\t¥");
		}
	fprintf(output1,"\n");
	} 
fclose(output1);

strcpy(outfile2,taxonname);
strcat(outfile2,"_L[Net Rate|CharComp].xls");
strcpy(outfile3,taxonname);
strcat(outfile3,"_lnL[Net Rate|CharComp].xls");
strcpy(outfile4,taxonname);
strcat(outfile4,"_S[Net Rate|CharComp].xls");

output2=fopen(outfile1,"w");
fprintf(output2,"Net Rate");
output3=fopen(outfile1,"w");
fprintf(output3,"Net Rate");
output4=fopen(outfile1,"w");
fprintf(output4,"Net Rate");
for (ch=0; ch<nchars; ++ch)	{
	if (ch<9)	{
		if (nchars>9 && nchars<100)	{
			fprintf(output2,"\tChar_0%d",ch+1);
			fprintf(output3,"\tChar_0%d",ch+1);
			fprintf(output4,"\tChar_0%d",ch+1);
			}
		else if (nchars>99)	{
			fprintf(output2,"\tChar_00%d",ch+1);
			fprintf(output3,"\tChar_00%d",ch+1);
			fprintf(output4,"\tChar_00%d",ch+1);
			}
		else	{
			fprintf(output2,"\tChar_%d",ch+1);
			fprintf(output3,"\tChar_%d",ch+1);
			fprintf(output4,"\tChar_%d",ch+1);
			}
		}
	else	{
		if (nchars>99)	{
			fprintf(output2,"\tChar_0%d",ch+1);
			fprintf(output3,"\tChar_0%d",ch+1);
			fprintf(output4,"\tChar_0%d",ch+1);
			}
		else	{
			fprintf(output2,"\tChar_%d",ch+1);
			fprintf(output3,"\tChar_%d",ch+1);
			fprintf(output4,"\tChar_%d",ch+1);
			}
		}
	}
fprintf(output2,"\n");
fprintf(output3,"\n");
fprintf(output4,"\n");

mr=0.05;	/* put in a way to read this!	*/
n=((int) (((double) (mxdl+1))/mr));
lnL=dmatrix(nchars,n);

d=0;
for (x=0.05; x<(mxdl+1); x=x+0.05)	{
	fprintf(output2,"4.3f",x);
	fprintf(output2,"4.3f",x);
	for (ch=0; ch<nchars; ++ch)	{
		y=0;
		for (s=1; s<=mxdl; ++s)	{
			y=y+(Poisson((double) s, x, 1) * Lsteps[ch][s]);
			}
		fprintf(output2,"\t%9.8f",y);
		fprintf(output2,"\t%9.8f",(lnL[ch][d]=log(y)));
		}
	++d;
	fprintf(output2,"\n");
	fprintf(output2,"\n");
	}
fclose(output2);
fclose(output3);

/* recalculate lnL as Support	*/
for (ch=0; ch<nchars; ++ch)	{
	y=maxdmatrixcol(lnL,n,ch);
	d=0;
	for (x=0.05; x<m; x=x+0.05)	{
		lnL[ch][d]-=y;
		++d;
		}
	}

for (x=0.05; x<m; x=x+0.05)	{
	fprintf(output4,"4.3f",x);
	d=0;
	for (ch=0; ch<nchars; ++ch)	{
		y=0;
		fprintf(output4,"\t%9.8f",lnL[ch][d]);
		}
	++d;
	fprintf(output4,"\n");
	}
fclose(output4);

free_ivector(chapos);
free_ivector(chuns);
free_ivector(chstcmb);
free_ulvector(compat);
free_lmatrix(smatrix,notu, nchars+1);
free_lmatrix(tree,notu+2,notu);
free_dmatrix(psteps,mxdl+1,nchars*notu);
free_dmatrix(Lsteps,nchars,mxdl+1);
free_dmatrix(lnL,nchars,n);
}

