double **likelystepsperchar(long **matrix, int *nstates, int *ctype, int *maxd, int *bias, int *depend, int nchars, int empcompat, int comptype, int notu, int fossil, double *mbl, int debug, int UNKNOWN, int INAP){int	a, b, c, d, ch, sp, steps, m, t;int dcomp, mxst, mxdl, redch;int	truns, mruns, cruns;int	apo, dapo, unsc, logdp, mxdp, avbr;int	*dstates, *dtypes, *dmaxd, *dbias, *ddpnd, *sychos;long *chvector, **chmat;long **tree, **dsmatrix;double runs;unsigned long *compat, *dcompat;double **Lsteps, **denom;/* get the real compatibilities for each character */compat=char_comp(nstates,notu,matrix,ctype,nchars,comptype,0,UNKNOWN,INAP);/* see if there are dependent inapplicables */a=0;for (ch=0; ch<nchars && a==0; ++ch)	if (depend[ch]!=ch)	a=1;sychos=ivector(nchars);				/* sychos[x] gives the number of characters depending of character x */if (a==1)	for (ch=0; ch<nchars && a==0; ++ch)	++sychos[depend[ch]]; else	clearivector(sychos,nchars,1);/* allocate memory for "dummy" arrays & matrices - these simply jackknife the real data to make a nchars-1 matrix & character info */dstates=ivector(nchars);dtypes=ivector(nchars);dmaxd=ivector(nchars-1);dbias=ivector(nchars-1);ddpnd=ivector(nchars-1);dsmatrix=lmatrix(notu, nchars);chvector=lvector(notu);mxdp=maxiarray(sychos,nchars);if (mxdp>1)	chmat=lmatrix(notu,mxdp);mxdl=maxiarray(maxd,nchars);	/* maximum number of changes possible for any character */if (mxdl>notu)	mxdl=notu-1;Lsteps=dmatrix(nchars,3*(mxdl+1));denom=dmatrix(nchars,mxdl+1);printf("\nFor each character, this routine will generate X trees and Y matrices per tree.\n");printf("  Each matrix is of N-1 characters, with N being the observed number.  The compatibility.\n");printf("  of each matrix will be the same as the compatibility of the other N-1 characters in\n");printf("  the observed matrix.  The program then will evolve the character K, K+1, K+2� steps Z times\n");printf("  (with K being one less than the number of states).  The program then calculate the\n");printf("  compatibility of the character.  It also will determine whether the proper number of states\n");printf("  and the proper number of derived taxa evolve.  Over X*Y*Z runs, this will determine the\n");printf("  probability of the observed structure of each character given K, K+1, etc., steps.  \n\n");printf("Enter the number of trees to use for each character (X above): ");scanf("%i",&truns);printf("\n");printf("Enter the number of matrices to evolve per tree (Y above): ");scanf("%i",&mruns);printf("\n");printf("Enter the number of times to replicate each number of steps (Z above): ");scanf("%i",&cruns);printf("\n");printf("The program will examine %d trees, generating %d matrices for each tree with\n",truns,mruns);printf("\teach number of steps replicated %d times to find L[setps | CPs, states, apomorphies].\n",cruns);tree=lmatrix(notu+2,notu);		/* also will give clade diversity in first cell */								/* finally, gives branch lengths in final two lines */printf("Doing tree ");for (t=0; t<truns; ++t)	{	/* evolve a tree to evolve characters over */	if (debug==1)	srand((unsigned int) t);	/* This is easier because you do not know how many nodes you'll get when sampling over time */	if (fossil==1)	tree=evolvetreeVenn(notu, mbl, fossil);	/* if no fossil taxa, then just make a cladogram and a Venn tree from it */	else	{		tree=evolvecladogram(notu,tree);		tree=clademember(tree, notu, notu-1);		}	printf("%d, matrix ",t+1);	for (m=0; m<mruns; ++m)	{		printf("%d, character ",m+1);		for (ch=0; ch<nchars; ++ch)	{			printf("%d\n",ch+1);			/* rewrite array of dependents to take into account that the characters are migrating */			/* count the unknowns and inapplicables */			unsc=0;			for (sp=0; sp<notu; ++sp)	if (matrix[sp][ch]==UNKNOWN || matrix[sp][ch]==INAP)	++unsc;						/* if we are on the first character, then evolve the remaining characters; we'll replace each one as we go	*/			if (ch==0)	{				redch=0;	/* reduced number of characters */				logdp=0;				for (c=0; c<nchars; ++c)	{					/************************************************************************************/					/* logdp counts the logical dependents; it will be at least 1 when we reach ch		*/					/*    if there are multiple dependents, then it will be >1;	 The reduced matrix		*/					/*    removed all characters dependent on ch (including ch); thus, the properies of	*/					/*    the remaining characters must be moved up; the logical dependence must be		*/					/*    changed because the character numbers change.  For example, if chars. 4 & 5	*/					/*    depend on char. 3 and we remove char. 1, then char. 3 becomes char. 2 and 	*/					/*	  chars. 4&5 become 3&4 and dependent on 2.  Note that when we remove 3, chars	*/					/*    1 & 2 stay put, but char. 6 is moved up to char. 3 because 4 & 5 also are		*/					/*    removed.  Those will be evolved along with char. 3 after the matrix evolves.  */					/************************************************************************************/					if (depend[c]==ch)	++logdp;					else	{						dstates[c-logdp]=nstates[c];						dtypes[c-logdp]=ctype[c];						dmaxd[c-logdp]=maxd[c];						dbias[c-logdp]=bias[c];						ddpnd[c-logdp]=depend[c]-logdp;						for (sp=0; sp<notu; ++sp)	dsmatrix[sp][c-logdp]=matrix[sp][c];						++redch;						}					}				/* clear vector for character change */				if (logdp==1)	for (sp=0; sp<notu; ++sp)	chvector[sp]=matrix[sp][ch];				dstates[nchars-1]=nstates[ch];				dtypes[nchars-1]=ctype[ch];				/* create reduced matrix and character info arrays excluding character ch and any characters logically dependent upon it*/				apo=autapo_char(matrix,notu,ch,UNKNOWN,INAP);				dcomp=empcompat-compat[ch];				for (c=0; c<nchars; ++c)	{					/* remove compatibilities of logically dependent characters */					if (depend[c]==ch && c!=ch)						dcomp=dcomp-compat[c];					}				/* if character is apomorphic, then skip this part - it has a perfect compatibility regardless */				dsmatrix=evolvetocompat(tree,dcomp,notu,dsmatrix,redch,dstates,dtypes,dbias,dmaxd,ddpnd,comptype,UNKNOWN,INAP);				dcompat=char_comp(dstates,notu,dsmatrix,dtypes,redch,comptype,0,UNKNOWN,INAP);								}			/* otherwise, put the last character into the slot where this character was; note, however, that it will be scored	*/			/* so that the remaining compatibility matches that of the rest of the matrix										*/			else	{				dstates[ch]=nstates[ch];				dtypes[ch]=ctype[ch];				dmaxd[ch]=maxd[ch];				dbias[ch]=bias[ch];				ddpnd[ch]=depend[ch]-logdp;				for (sp=0; sp<notu; ++sp)	dsmatrix[sp][ch]=nxtchar[sp];				}			if (ch<nchars-1)	nxtcomp=dcompat[ch+1];						if (debug==1)	srand((unsigned int) c+pow((t+2),(ch+1)));		//		printf("%d matrix ",t+1);			if (apo>1 || (apo+unsc)<notu-1)	{				printf("%d\n",m+1);				/* now evolve the reduced matrix until it matches the compatibility of the other nchar-1 characters 		*/				/* What we will do is seed the random number generator with a particular value and retain that value;		*/				/* We then use the evolvecompat routine to generate an array of compatibilities per step; if the desired	*/				/* compatibility is in that array, then we will use evolvematrix to evolve that matrix AFTER reseeding the	*/				/* random number generator with the same seed with the number of steps known to produce the desired 		*/				/* compatibility.  This means that we need to keep track of our random numbers!								*/								/* now, find the probability of getting observed compatibility AND observed apomorphy distribution given rates */					/* easy routine for characters with no dependents */				if (logdp==1)	{					for (steps=nstates[ch]-1; steps<maxd[ch]+4 && steps<(notu-1); ++steps)	{						for (c=0; c<cruns; ++c)	{							/* now, evolve character */							chvector=evolvecharacterNsteps(tree,steps,notu,chvector,nstates[ch],ctype[ch],bias[ch],UNKNOWN,INAP);							/* next, add it to the rest of the matrix */							for (sp=0; sp<notu; ++sp)	dsmatrix[sp][nchars-1]=chvector[sp];							/* find its compatibility */							d=char_nu_comp(nchars-1, dstates, notu, dsmatrix, dtypes, nchars, comptype, 0, UNKNOWN, INAP);							/* now see if compatibilities match */							if (nstates[ch]>2 && ctype[ch]==0)	{								mxst=0;								for (sp=0; sp<notu; ++sp)									if (chvector[sp]!=UNKNOWN && chvector[sp]!=INAP)										if (chvector[sp]>mxst)	mxst=chvector[sp];								}							else	mxst=nstates[ch]-1;							if (d==compat[ch] && mxst==(nstates[ch]-1))	{								++Lsteps[ch][steps];								dapo=autapo_char(dsmatrix,notu,nchars-1,UNKNOWN,INAP);								if (dapo==apo || (notu-dapo)==apo)	++Lsteps[ch][notu+steps];								}	/* end examination of whether we duplicated compatibility */							if (nxt==0 && d==nxtcomp)	{								equalivector(nxtchar,chvector,sp);								nxt=1;								}							++denom[ch][steps];	/* increment number of runs */							}	/* end run of possible steps */						}	/* end case where character has no logical dependents */					}				/* more difficult routine for those with depdendents */				else	{					for (steps=nstates[ch]-1; steps<maxd[ch]+4 && steps<(notu-1); ++steps)	{						/* now, evolve character */						chvector=evolvecharacterNsteps(tree,steps,notu,chvector,nstates[ch],ctype[ch],bias[ch],UNKNOWN,INAP);						/* next, add it to the rest of the matrix */						for (sp=0; sp<notu; ++sp)	dsmatrix[sp][nchars-logdp]=chvector[sp];						/* find the characters compatibility compatibility 						*/						/* note: you need to omit the dependent characters in two ways; 		*/						/*   1) simulated charact is #chars-dependents (not #chars-1)			*/						/*   2) because all dependents are compatible with their independent,	*/						/*      we want to look for #compatibilites - dependents				*/						d=char_nu_comp(nchars-logdp, dstates, notu, dsmatrix, dtypes, (nchars-logdp+1), comptype, 0, UNKNOWN, INAP);												/* now see if compatibilities match */						/* first, determine the maximum number of st */						if (nstates[ch]>2)	{							mxst=0;							for (sp=0; sp<notu; ++sp)								if (chvector[sp]!=UNKNOWN && chvector[sp]!=INAP)									if (chvector[sp]>mxst)	mxst=chvector[sp];							}						else	mxst=nstates[ch]-1;						if (d==compat[ch]-logdp && mxst==(nstates[ch]-1))	{							++Lsteps[ch][steps];							dapo=autapo_char(dsmatrix,notu,nchars-1,UNKNOWN,INAP);							/* if we replicate the number of taxa with derived conditions, continue */							if (dapo==apo || (notu-dapo)==apo)	{								++Lsteps[ch][mxdl+steps+1];								/* if we have the right number of apomorphies, then evolve the dependent characters */								if (dapo==apo)	{									/* the minimum number of steps could be 0 IF the number of steps for the independent matches 	*/									/* the number of states for the dependent character												*///										for ()									}								}	/* tally cases where everything matches */							}	/* end examination of whether we duplicated compatibility */						if (dapo==apo || (notu-dapo)==apo)							++Lsteps[ch][2*(mxdl+1)+steps];												}	/* end run of possible steps */					}	/* end case where character has logical dependents */				for (b=(m+1); b>=1; b=b/10)	printf("\b");				printf("\b");				}	/* end case for non-apomorphic taxa */			/* if apomorphic, then just determine the probability of apomorphy given K, K+1, K+2, etc. steps. */			else	{				for (steps=nstates[ch]-1; steps<maxd[ch]+3 && steps<notu; ++steps)	{					chvector=evolvecharacterNsteps(tree,steps,notu,chvector,nstates[ch],ctype[ch],bias[ch],UNKNOWN,INAP);					dapo=0;					/* if only one taxon has '1' or if all but one have '1' then it is apomorphic */					for (sp=0; sp<notu; ++sp)	if (chvector[ch]==1)	++dapo;					if (dapo==1 || dapo==(notu-unsc-1))	{						++Lsteps[ch][steps];						++Lsteps[ch][notu+steps];						}	/* tally apomorphic successes */					}	/* end X steps */				}	/* end routine for apomorphic states */			}	/* end evolved trees */		}	dstates[ch]=nstates[ch];	dtypes[ch]=dtypes[ch];	dmaxd[ch]=maxd[ch];	printf("\b\b\b\b\b\b\b");	for (b=(t+1); b>=1; b=b/10)	printf("\b");	printf("\b");	}}