#define CompatibilityFunctions
#include <time.h>
#include "CompatibilityFunctions.h"
#include "matrixanalysis.h"
#include "matrixchange.h"
#include "matrixreading.h"
#include "memory.h"
#include "minmax.h"
#include "MonteCarloPhylogenyFunctions.h"
#include "TreeRead.h"

int *charcombos(int *t1, int *t2, int st, int notu, int UNKNOWN, int INAP)
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
/* unordecompatible - determines whether a character pair involving an unordered multistate are compatible
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
int	a, b, c, d, sp, maxst;
int	total, incompatible=0;
int *combos;

maxst=st1;
if (st2>maxst)	maxst=st2;
combos=charcombos(t1,t2, st1, notu, UNKNOWN, INAP);

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
							while (t1[sp+1]==a)	++sp;				/* get through a� taxa */
						else
							while (t2[sp+1]<c || (t2[sp+1]>c && t2[sp+1]<d))	++sp;	/* skip states between c and d */
						}	/* end search for a� taxa */
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
						}	/* end search for b� taxa */
					else if (t1[sp]>b || (t1[sp]==UNKNOWN ||t1[sp]==INAP))	sp=notu;	/* skip ahead if past b */
					else
						while (t1[sp+1]<a || (t1[sp+1]>a && t1[sp+1]<b))	++sp;		/* skip states between a and b */
					}	/* end examination of taxa */
				if (total>=4)	{
					incompatible=1;
					a=b=c=d=maxst;		/* sets all loop variables to end of loop or beyond, ending routine */
					}	/* if incompatible, then there is no point in searching further */
				}	/* end examination of combinations for �d */
			}	/* end examination of combinations for �c */
		}	/* end examination of combinations for b� */
	}	/* end examination of combinations for a� */

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
int ordcompatible(int *t1, int *t2, int type2, int st1, int st2, int notu, int UNKNOWN, int INAP)
{
int	c, d, drop, rise, last;
int incompatible=0;

drop=rise=0;
c=0;
while ((t2[c]==UNKNOWN || t2[c]==INAP) && c<notu) ++c;

if (c<notu)	{
	last = t2[c];
	for (d=c+1; d<notu; ++d)	{
		if (t2[d]<st1 && t1[c]<st2)	{
			if (t2[d]<last)	++drop;
			if (t2[d]>last)	++rise;
			last = t2[d];
			
			if (type2==0 && (drop+rise)>=st2)	{
				incompatible=1;
				d=notu;
				}
			else if ((drop+rise)>st2)	{
				incompatible=1;
				d=notu;
				}
			}
		}
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
unsigned long **compatible(int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int OUTGROUP, int UNKNOWN, int INAP)
{
int a, b, d, ch1, ch2, sp, placer, incompatible, inhiercompat, mxst, hi, lo, compsp, maxst;
int *autap, *character1, *character2, *mnst, *mtch, *tallied, *taxst, *test1, *test2;
unsigned long	**comatrix;

comatrix=ulmatrix(nchars,nchars);
mnst=ivector(nchars);
character1=ivector(notu);
character2=ivector(notu);
test1=ivector(notu);
test2=ivector(notu);

maxst=maxiarray(nstates,notu);
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

	/* Make sure that outgroup state is appropriately coded	*/
	if (nstates[ch1]==2 && chmatrix[OUTGROUP][ch1]==1)	{
		for (sp=0; sp<notu; ++sp)	{
			if (chmatrix[sp][ch1]==0)			chmatrix[sp][ch1]=1;
			else if (chmatrix[sp][ch1]==1)		chmatrix[sp][ch1]=0;
			}
		}					
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
		compsp=0;	/* this counts the number of comparable species - i.e., those scored for both characters */
		for (sp=0; sp<notu; ++sp)	{
			placer = 0;
			while (sp<notu && ((character1[sp]==UNKNOWN || character1[sp]==INAP) || (character2[sp]==UNKNOWN || character2[sp]==INAP)))
				++sp;
			if (sp>=notu)	break;
			++compsp;
			/**** Sort on 1st character ****/
			for (placer=0; (character1[sp]>test1[placer] && placer< notu); ++placer);
			/**** Sort on 2nd character ****/
			for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]>test2[placer]) && placer<notu); ++placer);
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
				a=test1[0];
				for (sp=0; sp<compsp; ++sp)	{
					if (test1[sp]!=INAP && test1[sp]!=UNKNOWN)	{
						test1[sp]=test1[sp]-a;
						}
					else	sp=compsp;
					}
				}

			/**** Routine for Binary Characters ****/
			if (nstates[ch1]==2 && nstates[ch2]==2)
				incompatible=binarycompatible(test1, test2, notu, UNKNOWN, INAP);
			
			/**** Routine for Multistate Characters ****/
			else	{
				/* if ordered multistate	*/
				if (type[hi]==0)
					incompatible=ordcompatible(test1, test2, type[lo],nstates[hi], nstates[lo], notu, UNKNOWN, INAP);
					/* end routine for ordered multistates	*/
				else
					incompatible=unordcompatible(test1, test2, nstates[hi], nstates[lo], notu, UNKNOWN, INAP);
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
						for (a=0; test1[a]==mnst[ch1] && a<compsp; ++a);
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
							for (b=a; (test1[b]==1+mnst[ch1] && b<compsp); ++b);
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
int nu_comp(int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int OUTGROUP, int UNKNOWN, int INAP)
{
int a, b, d, ch1, ch2, sp, placer, incompatible, inhiercompat, mxst, hi, lo, CP=0, compsp, maxst;
int *character1, *character2, *test1, *test2, *mnst, *autap;
int *tallied, *taxst;

mnst=ivector(nchars);
character1=ivector(notu);
character2=ivector(notu);
test1=ivector(notu);
test2=ivector(notu);

autap=autapomorphies(chmatrix, nstates, notu, nchars, UNKNOWN, INAP);

maxst=maxiarray(nstates,nchars);
for (a=0; a<nchars; ++a)	{
	mnst[a]=100;
	for (b=0; b<notu; ++b)
		if (chmatrix[b][a]<mnst[a] && (chmatrix[b][a]!=UNKNOWN && chmatrix[b][a]!=INAP))	
			mnst[a]=chmatrix[b][a];
	}

tallied=ivector(maxst+1);
taxst=ivector(maxst+1);

for (ch1=0; ch1<nchars; ch1++)	{
	/* if the character is an autapomorphy or invariant, do not waste your time */
	/*    the character must be compatible with all other characters	*/
	for (ch1=ch1; (autap[ch1]<2 && ch1<nchars); ++ch1)	{
		for (ch2=ch1+1; ch2<nchars; ++ch2)	{
			if (comptype==0)	++CP;
			}	
		}	

/* Make sure that outgroup state is appropriately coded	*/
	if (nstates[ch1]==2 && chmatrix[OUTGROUP][ch1]==1)	{
		for (sp=0; sp<notu; ++sp)	{
			if (chmatrix[sp][ch1]==0)			chmatrix[sp][ch1]=1;
			else if (chmatrix[sp][ch1]==1)		chmatrix[sp][ch1]=0;
			}
		}					
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
		compsp=0;	/* this counts the number of comparable species - i.e., those scored for both characters */
		for (sp=0; sp<notu; ++sp)	{
			placer = 0;
			while (sp<notu && ((character1[sp]==UNKNOWN || character1[sp]==INAP) || (character2[sp]==UNKNOWN || character2[sp]==INAP)))
				++sp;
			if (sp>=notu)	break;
			++compsp;
			/**** Sort on 1st character ****/
			for (placer=0; (character1[sp]>test1[placer] && placer< notu); ++placer);
			/**** Sort on 2nd character ****/
			for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]>test2[placer]) && placer<notu); ++placer);
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
				a=test1[0];
				for (sp=0; sp<compsp; ++sp)	{
					if (test1[sp]!=INAP && test1[sp]!=UNKNOWN)	{
						test1[sp]=test1[sp]-a;
						}
					else	sp=compsp;
					}
				}

			/**** Routine for Binary Characters ****/
			if (nstates[ch1]==2 && nstates[ch2]==2)
				incompatible=binarycompatible(test1,test2,notu,UNKNOWN,INAP);
			
			/**** Routine for Multistate Characters ****/
			else	{
				/* if ordered multistate	*/
				if (type[hi]==0)
					incompatible=ordcompatible(test1,test2,type[lo],nstates[hi],nstates[lo],notu,UNKNOWN,INAP);

				/**** Routine for Unorderd Multistates ****/
				/* this basically breaks multistates into multiple binary combinations and finds whether any state pair */
				/*		has four combinations and thus is incombatible */
				else
					incompatible=unordcompatible(test1, test2, nstates[hi], nstates[lo], notu, UNKNOWN, INAP);
				}	/* end routine for multstates	*/

			inhiercompat = 1;

			if (incompatible==0 && comptype==0)	++CP;
			else if (incompatible==0)	{						
				a=0;
				/* simple test if both are binary: */
				if (nstates[ch1]==2 && nstates[ch2]==2)	{
					a=0;
					/* first look for 00 & 01 */
					for (a=0; (test1[a]==mnst[ch1] && a<compsp); ++a);
					for (sp=0; sp<=a; ++sp)	{
						/* 2nd state found?	*/
						if (test2[sp]!=test2[0])	{
							inhiercompat = 0;
							sp = compsp;
							}				
						}
					/* if only 00's found (or 01's), make sure that two nstates are found with 1- */
					if (inhiercompat==1)	{
						b=a;
						for (b=a; (test1[b]==1+mnst[ch1] && b<compsp); ++b);
						for (sp=a; sp<=b; ++sp)	{
							if (test2[sp]!=test2[a])	{
								inhiercompat = 0;
								sp = compsp;
								}
							} 
						}
					}
				/* if two multistates with different # nstates OR different numbers of taxa
				 with derived conditions, then they must be HC if compatible at all 
				 (unless they are "complements")				*/
				
				else	{
					/* check to see if any derived nstates have multiple counterparts	*/
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

				/* 1 for hierarchically compatible, 0 for hierarchically incompatible	*/
				if (inhiercompat==0)
					++CP;						
				}	/* end routine for hierachically compatible */
			}	/* end test of whether they are compatible */
		}	/* end comparison between ch1 and ch2 */
	}	/* end ch1 */
		
free_ivector(character1);
free_ivector(character2);
free_ivector(test1);
free_ivector(test2);
free_ivector(mnst);
free_ivector(autap);
free_ivector(tallied);
free_ivector(taxst);

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
unsigned long *char_comp(int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int OUTGROUP, int UNKNOWN, int INAP)
{
int a, b, d, ch1, ch2, sp, placer, incompatible, inhiercompat, mxst, hi, lo, maxst;
int *character1, *character2, *test1, *test2, *mnst, *autap;
int *tallied, *taxst;
unsigned long *charcomps;

charcomps=ulvector(nchars);
mnst=ivector(nchars);
character1=ivector(notu);
character2=ivector(notu);
test1=ivector(notu);
test2=ivector(notu);
autap=autapomorphies(chmatrix, nstates, notu, nchars, UNKNOWN, INAP);

maxst=maxiarray(nstates,nchars);
for (a=0; a<nchars; ++a)	{
	charcomps[a]=0;
	mnst[a]=100;
	for (b=0; b<notu; ++b)	{
		if (chmatrix[b][a]<mnst[a] && (chmatrix[b][a]!=UNKNOWN && chmatrix[b][a]!=INAP))	
			mnst[a]=chmatrix[b][a];
		}
	}

tallied=ivector(maxst+1);
taxst=ivector(maxst+1);

for (ch1=0; ch1<nchars; ch1++)	{
	/* if the character is an autapomorphy or invariant, do not waste your time */
	/*    the character must be compatible with all other characters	*/
	for (ch1=ch1; (autap[ch1]<2 && ch1<nchars); ++ch1)	{
		for (ch2=ch1+1; ch2<nchars; ++ch2)	{
			if (comptype==0)	{
				++charcomps[ch1];
				++charcomps[ch2];
				}
			}	
		}	

/* Make sure that outgroup state is appropriately coded	*/
	if (nstates[ch1]==2 && chmatrix[OUTGROUP][ch1]==1)	{
		for (sp=0; sp<notu; ++sp)	{
			if (chmatrix[sp][ch1]==0)			chmatrix[sp][ch1]=1;
			else if (chmatrix[sp][ch1]==1)		chmatrix[sp][ch1]=0;
			}
		}					
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

		/*sort on placers for comparisons */
		for (sp=0; sp<notu; ++sp)	{
			placer = 0;
			/**** Sort on 1st character ****/
			for (placer=0; (character1[sp]>test1[placer] && placer< notu); ++placer);
			/**** Sort on 2nd character ****/
			for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]>test2[placer]) && placer<notu); ++placer);
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

		/**** Routine for Binary Characters ****/
		if (nstates[ch1]==2 && nstates[ch2]==2)
			incompatible=binarycompatible(test1,test2,notu,UNKNOWN,INAP);
		
		/**** Routine for Multistate Characters ****/
		else	{
			/* if ordered multistate	*/
			if (type[hi]==0)
				incompatible=ordcompatible(test1,test2,type[lo],nstates[hi],nstates[lo],notu,UNKNOWN,INAP);

			/**** Routine for Unorderd Multistates ****/
			/* this basically breaks multistates into multiple binary combinations and finds whether any state pair */
			/*		has four combinations and thus is incombatible */
			else
				incompatible=unordcompatible(test1, test2, nstates[hi], nstates[lo], notu, UNKNOWN, INAP);

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
					for (a=0; (test1[a]==mnst[ch1] && a<notu); ++a);
					for (sp=0; sp<=a; ++sp)	{
						/* 2nd state found?	*/
						if (test2[sp]!=test2[0])	{
							inhiercompat = 0;
							sp = notu;
							}				
						}
					/* if only 00's found (or 01's), make sure that two states are found with 1- */
					if (inhiercompat==1)	{
						b=a;
						for (b=a; (test1[b]==1+mnst[ch1] && b<notu); ++b);
						for (sp=a; sp<=b; ++sp)	{
							if (test2[sp]!=test2[a])	{
								inhiercompat = 0;
								sp = notu;
								}
							} 
						}
					}
				/* if two multistates with different # states OR different numbers of taxa
				 with derived conditions, then they must be HC if compatible at all 
				 (unless they are "complements")				*/
				
				else	{
					/* check to see if any derived states have multiple counterparts	*/
					for (d=0; d<notu; ++d)	{
						if (test1[d]>0 && test1[d]<mxst)	{
							if (test2[d]>0 && test2[d]<mxst)	{
								inhiercompat = 0;
								d=notu;
								}
							}
						}
					/* even if all deriveds for one state paired with single derived, this might not be true */
					if (inhiercompat==1)	{
						for (d=0; d<notu; ++d)	{
							if (test2[d]>0 && test2[d]<mxst)	{
								if (test1[d]>0 && test1[d]<mxst)	{
									inhiercompat = 0;
									d=notu;
									}
								}
							}
						}
					}
				}
			if (inhiercompat==0)	{
				++charcomps[ch1];
				++charcomps[ch2];
				/* 1 for hierarchically compatible, 0 for hierarchically incompatible	*/
				}
			}
		}
	}
		
free_ivector(character1);
free_ivector(character2);
free_ivector(test1);
free_ivector(test2);
free_ivector(mnst);
free_ivector(autap);
free_ivector(tallied);
free_ivector(taxst);

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
Function returning the number of mutual compatibilities (see O�Keefe & Wagner 2001 Syst. Biol.
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
int char_nu_comp(int CH, int *nstates, int notu, long **chmatrix, int *type, int nchars, int comptype, int OUTGROUP, int UNKNOWN, int INAP)
{
int a, b, d, ch2, sp, placer, incompatible, inhiercompat, mxst, hi, lo, NC=0, maxst;
int *character1, *character2, *test1, *test2, *mnst, *autap;
int *tallied, *taxst;

mnst=ivector(nchars);
character1=ivector(notu);
character2=ivector(notu);
test1=ivector(notu);
test2=ivector(notu);

autap=autapomorphies(chmatrix, nstates, notu, nchars, UNKNOWN, INAP);

maxst=maxiarray(nstates,nchars)+1;
for (a=0; a<nchars; ++a)	{
	mnst[a]=100;
	for (b=0; b<notu; ++b)
		if (chmatrix[b][a]<mnst[a] && (chmatrix[b][a]!=UNKNOWN && chmatrix[b][a]!=INAP))	mnst[a]=chmatrix[b][a];
	}

tallied=ivector(maxst+1);
taxst=ivector(maxst+1);

/* if the character is an autapomorphy or invariant, do not waste your time */
/*    the character must be compatible with all other characters	*/
if (autap[CH]==1)	NC=nchars-1;

else	{
	/* Make sure that outgroup state is appropriately coded	*/
	if (nstates[CH]==2 && chmatrix[OUTGROUP][CH]==1)	{
		for (sp=0; sp<notu; ++sp)	{
			if (chmatrix[sp][CH]==0)			chmatrix[sp][CH]=1;
			else if (chmatrix[sp][CH]==1)		chmatrix[sp][CH]=0;
			}
		}					

	for (ch2=0; ch2<nchars; ++ch2)	{
		if (ch2==CH)	++ch2;
		if (ch2>=nchars)	break;

		mxst = nstates[CH];
		for (sp=0; sp<notu; ++sp)			character1[sp]=chmatrix[sp][CH];
		/* characters are compatible if second is autapomorphic */
		for (ch2=ch2; (autap[ch2]<2 && ch2<nchars); ++ch2)	{
			if (comptype==0)	++NC;
			}	
		if (ch2>=nchars)	break;

		/* set up arrays for testing compatibility */
		for (a=0; a<notu; ++a)
			test1[a]=test2[a]=100;	 
		for (sp=0; sp<notu; ++sp)
			character2[sp]=chmatrix[sp][ch2];

		/*** Determine compatibility ***/		
		incompatible = 0;

		if ((nstates[CH]>nstates[ch2] && type[CH]==1) || (nstates[ch2]>nstates[CH] && type[ch2]==0))	{
			for (sp=0; sp<notu; ++sp)	{
				b = character2[sp];
				character2[sp]=character1[sp];
				character1[sp]=b;
				mxst = nstates[ch2];
				}
			hi=ch2;
			lo=CH;
			}
		else	{
			hi=CH;
			lo=ch2;
			}

		/*sort on placers for comparisons */
		for (sp=0; sp<notu; ++sp)	{
			placer = 0;
			/**** Sort on 1st character ****/
			for (placer=0; (character1[sp]>test1[placer] && placer< notu); ++placer);
			/**** Sort on 2nd character ****/
			for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]>test2[placer]) && placer<notu); ++placer);
			for (b=sp; b>placer; --b)	{
				test1[b]=test1[b-1];
				test2[b]=test2[b-1];
				}
			test1[placer] = character1[sp];
			test2[placer] = character2[sp];
			}

		/**** Routine for Binary Characters ****/
		if (nstates[CH]==2 && nstates[ch2]==2)
			incompatible=binarycompatible(test1,test2,notu,UNKNOWN,INAP);
		
		/**** Routine for Multistate Characters ****/
		else	{
			/* if ordered multistate	*/
			if (type[hi]==0)
				incompatible=ordcompatible(test1,test2,type[lo],nstates[hi],nstates[lo],notu,UNKNOWN,INAP);

			/**** Routine for Unorderd Multistates ****/
			/* this basically breaks multistates into multiple binary combinations and finds whether any state pair */
			/*		has four combinations and thus is incombatible */
			else
				incompatible=unordcompatible(test1, test2, nstates[hi], nstates[lo], notu, UNKNOWN, INAP);

			}	/* end routine for multstates	*/

		inhiercompat = 1;
		if (incompatible==0)	{
			if (comptype==0)	++NC;
			else	{						
				a=0;
				/* simple test if both are binary: */
				if (nstates[CH]==2 && nstates[ch2]==2)	{
					a=0;
					/* first look for 00 & 01 */
					for (a=0; (test1[a]==mnst[CH] && a<notu); ++a);
					for (sp=0; sp<=a; ++sp)	{
						/* 2nd state found?	*/
						if (test2[sp]!=test2[0])	{
							inhiercompat = 0;
							sp = notu;
							}				
						}
					/* if only 00's found (or 01's), make sure that two states are found with 1- */
					if (inhiercompat==1)	{
						b=a;
						for (b=a; (test1[b]==1+mnst[CH] && b<notu); ++b);
						for (sp=a; sp<=b; ++sp)	{
							if (test2[sp]!=test2[a])	{
								inhiercompat = 0;
								sp = notu;
								}
							} 
						}
					}
				/* if two multistates with different # states OR different numbers of taxa
				 with derived conditions, then they must be HC if compatible at all 
				 (unless they are "complements")				*/
				
				else	{
					/* check to see if any derived states have multiple counterparts	*/
					for (d=0; d<notu; ++d)	{
						if (test1[d]>0 && test1[d]<mxst)	{
							if (test2[d]>0 && test2[d]<mxst)	{
								inhiercompat = 0;
								d=notu;
								}
							}
						}
					/* even if all deriveds for one state paired with single derived, this might not be true */
					if (inhiercompat==1)	{
						for (d=0; d<notu; ++d)	{
							if (test2[d]>0 && test2[d]<mxst)	{
								if (test1[d]>0 && test1[d]<mxst)	{
									inhiercompat = 0;
									d=notu;
									}
								}
							}
						}
					}
				}
			if (inhiercompat==0)	{
				++NC;
				/* 1 for hierarchically compatible, 0 for hierarchically incompatible	*/
				}
			}
		}
	}
free_ivector(character1);
free_ivector(character2);
free_ivector(test1);
free_ivector(test2);
free_ivector(mnst);
free_ivector(autap);
free_ivector(tallied);
free_ivector(taxst);

return NC;
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
//dependent=ivector(nchars);

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
			for (placer=0; (character1[sp]>test1[placer] && placer< notu); ++placer);
			/**** Sort on 2nd character ****/
			for (placer=placer; ((character1[sp]==test1[placer]&&character2[sp]>test2[placer]) && placer<notu); ++placer);

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
			for (sp=sp; test1[sp]==d && sp<notu; ++sp);
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
				}	/* end test to see if b is new in a�b comparison */
			}	/* end search through taxa */
			
		for (a=0; (a<=maxst /*&& pair[a][0]!=-1*/); ++a)	{
			for (a=a; pair[a][0]==-1 & a<maxst; ++a);
			if (a>maxst)	break;
			for (b=0; (b<=maxst /*&& pair[a][b]!=-1*/); ++b)	{
				for (b=b; pair[a][b]==-1; ++b);
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
double *likeallstepsgivencomp(long **matrix, int *ctype, int *nstates, int nchars, int notu, int empcomp, int comptype, int out, int fossils, double *mbl, int *bias, int *maxd, int pars, int debug, int UNKNOWN, int INAP)
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
		simcomp=nu_comp(nstates, notu, simat, ctype, nchars, comptype, out, UNKNOWN, INAP);
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
int		a, b,st,mreps,treps,truns, mruns;
long	*compats;
long	**simat, **tree;
double	x;
double	**Pcomp;

Pcomp=dmatrix((nchars*(nchars-1))/2,1+(max-pars));

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
	printf("%d matrix ",truns+1);
	
	for (mruns=0; mruns<mreps; ++mruns)	{

		if (debug==1)	srand((unsigned int) (truns+1)*(mruns+1));
		printf("%d\n",mruns+1);
		equallmatrix(simat,matrix,notu,nchars);

		compats=evolvecompat(tree, notu, simat, nchars, nstates, ctype, bias, maxd, pars, max, comptype, UNKNOWN, INAP);

		for (st=pars; st<=max; ++st)	{
			a=compats[st-pars];
			++Pcomp[a][st-pars];
			}
		for (b=(mruns+1); b>=1; b=b/10)	printf("\b");
		printf("\b");
		}
	printf("\b\b\b\b\b\b\b");
	for (b=(truns+1); b>=1; b=b/10)	printf("\b");
	printf("\b");
	
	free_lvector(compats);
	}

x=treps*mreps;
for (b=0; b<(1+(max-pars)); ++b)	{
	x=0;
	for (a=0; a<(nchars*(nchars-1))/2; ++a)
		x=Pcomp[a][b]+x;
	for (a=0; a<(nchars*(nchars-1))/2; ++a)
		Pcomp[a][b]=Pcomp[a][b]/x;
	}

free_lmatrix(simat,notu,nchars);
free_lmatrix(tree,notu+2,notu);

return Pcomp;
}

/*
likelystepsperchar - calculates the probability of X compatibilities given Y steps and observed matrix.
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
int	a, b, c, d, ch, sp, steps;
int dcompat, mxst, mxdl, redch;
int	truns, mruns, cruns;
int	apo, dapo, unsc, logdp, mxdp;
int	*dstates, *dtypes, *dmaxd, *dbias, *ddpnd, *sychos;
long *chvector, **chmat;
long **tree, **dsmatrix;
double runs;
unsigned long *compat;
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
mxdp=maxiarray(sychos,nchars);
if (mxdp>1)
	chmat=lmatrix(notu,mxdp);

mxdl=maxiarray(maxd,nchars);	/* maximum number of changes possible for any character */
if (mxdl>notu)	mxdl=notu-1;
Lsteps=dmatrix(nchars,2*mxdl);
denom=dmatrix(nchars,mxdl);

printf("For each character, this routine will generate X trees and Y matrices per tree.\n");
printf("  Each matrix is of N-1 characters, with N being the observed number.  The compatibility.\n");
printf("  of each matrix will be the same as the compatibility of the other N-1 characters in\n");
printf("  the observed matrix.  The program then will evolve the character K, K+1, K+2� steps Z times\n");
printf("  (with K being one less than the number of states).  The program then calculate the\n");
printf("  compatibility of the character.  It also will determine whether the proper number of states\n");
printf("  and the proper number of derived taxa evolve.  Over X*Y*Z runs, this will determine the\n");
printf("  probability of the observed structure of each character given K, K+1, etc., steps.  \n");

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


for (ch=0; ch<nchars; ++ch)	{

	redch=0;	/* reduced number of characters */
	logdp=0;

	/* rewrite array of dependents to take into account that the characters are migrating */
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
			++redch;
			}
		}

	/* clear vector for character change */
	if (logdp==1)	for (sp=0; sp<notu; ++sp)	chvector[sp]=matrix[sp][ch];
	else
	

	/* count the unknowns and inapplicables */
	unsc=0;
	for (sp=0; sp<notu; ++sp)	if (matrix[sp][ch]==UNKNOWN || matrix[sp][ch]==INAP)	++unsc;
	
	dstates[nchars-1]=nstates[ch];
	dtypes[nchars-1]=ctype[ch];

	apo=autapo_char(matrix,notu,ch,UNKNOWN,INAP);
	dcompat=empcompat-compat[ch];
	for (c=0; c<nchars; ++c)	{
		/* remove compatibilities of logically dependent characters */
		if (depend[c]==ch && c!=ch)
			dcompat=dcompat-compat[c];
		}
	/* create reduced matrix and character info arrays excluding character ch and any characters logically dependent upon it*/

	c=0;
	for (a=ch; a<nchars; ++a)
		if (depend[a]!=ch)		for (b=0; b<notu; ++b)	dsmatrix[b][c]=matrix[b][a];
	
	/* evolve a tree to evolve characters over */
	if (debug==1)	srand((unsigned int) ch);
	for (a=0; a<truns; ++a)	{
		/* This is easier because you do not know how many nodes you'll get when sampling over time */
		if (fossil==1)	tree=evolvetreeVenn(notu, mbl, fossil);
		/* if no fossil taxa, then just make a cladogram and a Venn tree from it */
		else	{
			tree=evolvecladogram(notu,tree);
			tree=clademember(tree, notu, notu-1);
			}

		if (debug==1)	srand((unsigned int) ch);
		/* if character is apomorphic, then skip this part - it has a perfect compatibility regardless */
		if (apo>1 || (apo+unsc)<notu-1)	{
			for (b=0; b<mruns; ++b)	{
				/* now evolve the reduced matrix until it matches the compatibility of the other nchar-1 characters 		*/
				/* What we will do is seed the random number generator with a particular value and retain that value;		*/
				/* We then use the evolvecompat routine to generate an array of compatibilities per step; if the desired	*/
				/* compatibility is in that array, then we will use evolvematrix to evolve that matrix AFTER reseeding the	*/
				/* random number generator with the same seed with the number of steps known to produce the desired 		*/
				/* compatibility.  This means that we need to keep track of our random numbers!								*/

				dsmatrix=evolvetocompat(tree,dcompat,notu,dsmatrix,redch,dstates,dtypes,dbias,dmaxd,ddpnd,comptype,UNKNOWN,INAP);
				
				/* now, find the probability of getting observed compatibility AND observed apomorphy distribution given rates */
				for (c=0; c<cruns; ++c)	{
					/* easy routine for characters with no dependents */
					if (logdp==1)	{
						for (steps=nstates[ch]-1; steps<maxd[ch]+4 && steps<(notu-1); ++steps)	{
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
								dapo=autapo_char(dsmatrix,notu,nchars-1,UNKNOWN,INAP);

								if (dapo==apo || (notu-dapo)==apo)	++Lsteps[ch][notu+steps];
								}	/* end examination of whether we duplicated compatibility */
							++denom[ch][steps];	/* increment number of runs */
							}	/* end run of possible steps */
						}	/* end case where character has no logical dependents */
					
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
								dapo=autapo_char(dsmatrix,notu,nchars-1,UNKNOWN,INAP);

								/* if we replicate the number of taxa with derived conditions, continue */
								if (dapo==apo || (notu-dapo)==apo)	{
									++Lsteps[ch][notu+steps];
									/* if we have the right number of apomorphies, then evolve the dependent characters */
									if (dapo==apo)	{
										/* the minimum number of steps could be 0 IF the number of steps for the independent matches 	*/
										/* the number of states for the dependent character												*/
//										for ()
										}
									}	/* tally cases where everything matches */
								}	/* end examination of whether we duplicated compatibility */
							
							}	/* end run of possible steps */
						}	/* end case where character has logical dependents */
					}	/* end runs for evolving character */
				}	/* end evolved matrices */
			}	/* end case for non-apomorphic taxa */
		/* if apomorphic, then just determine the probability of apomorphy given K, K+1, K+2, etc. steps. */
		else	{
			for (b=0; b<mruns*cruns; ++b)	{
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
				}	/* do this mruns*cruns times just to make it as common as the other runs */
			}	/* end routine for apomorphic states */
		}	/* end evolved trees */
	
	dstates[ch]=nstates[ch];
	dtypes[ch]=dtypes[ch];
	dmaxd[ch]=maxd[ch];
	}

runs=truns*mruns*cruns;

free_ulvector(compat);
free_ivector(dstates);
free_ivector(dtypes);
free_ivector(dmaxd);
free_ivector(dbias);
free_lmatrix(dsmatrix,notu, nchars);
free_lvector(chvector);
if (mxdp>1)	free_lmatrix(chmat,notu,mxdp);
return Lsteps;
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
int statepairs(long **mat, int notu, int nchar, int maxstate, int UNKNOWN, int INAP)
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

