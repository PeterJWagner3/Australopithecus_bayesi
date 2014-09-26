/* CALCULATE LIKELIHOOD THAT A GIVEN DISTRIBUTION FITS THE GEOMETRIC SERIES. NEEDS:	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest	- ntaxa (size of array)	- nspec (sum of array)RETURNS:	- result[0]: log likelihood	- result[1]: slope of geometric series	- result[2]: optimal richness			COMMENTS:	- SUPINC, FITINC defined in header file***********************************************************************/double *o_fit_gs(int *empdist, int ntaxa, int nspec) {int i = 0, j;				/* LOOP VARIABLE	*/int r = 0;					/* LOOP RICHNESS	*/int ri;						/* richness increment	*/int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*//*int mxspc;					/* maximum number of possible finds to consider when calculating P */double ev = 0.000f;			/* LOOP SLOPE 		*/double emin = 1.0f;			/* min slope		*/double ein = 0.000f;		/* initial slope	*/double ei = 0.000f;			/* how much to increment ev in each loop							*/double es[3];				/* previous log likelihoods (cell number = num previous).		*/double rs[3];				/* previous log likelihoods (cell number = num previous).		*/double bep[3];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/double *brp;				/* BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array	*/double *fitdist;			/* fit distribution												*/double *expect;				/* expected number of species with 0�max finds					*///double *exp;				/* expected number of species with 0�max finds					*/long *obsrvd;				/* observed number of species with 0�max finds					*/long *obsfnd;				/* array of unique values for numbers of finds					*/double pbes;				/* PREVIOUS BEST SUPPORT (E VALUE)								*/brp=dvector(3);for (i = 0 ; i < 3 ; i++) rs[i] = 0.0f;for (i = 0 ; i < 3 ; i++) brp[i]= -1.0*DBL_MAX;/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */obsrvd=ilhistogram(empdist,ntaxa);/* we want one plus the number of different sample sizes because some species have a sample size of 0 (probably) */unqfnd=1+countuniqivector(empdist,ntaxa);obsfnd=lvector(unqfnd);obsfnd[j=0]=0;for (i=ntaxa-1; i>=0; --i)	if (empdist[i]!=empdist[i-1])	obsfnd[++j]=empdist[i];/************************************************************************************************************************//* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */ein=pow(empdist[0],(((double) empdist[ntaxa-1])/((double) ntaxa)));ri=ntaxa/10;if (ri<2)	ri=2;// increment true richness until that fails to improve likelihoodfor (r=ntaxa; abs(ri)>0; r+=ri) {		obsrvd[0]=r-ntaxa;	pbes = 0.0f;	ei = (double) FITINC / 10;	// increment starting at 0.1 since we are almost never going to find slopes of 2+		for (i = 0 ; i < 3 ; i++) bep[i] = -1.0*DBL_MAX;	for (i = 0 ; i < 3 ; i++) es[i] = 0.0f;	// increment slope until that fails to improve likelihood or resolution limit reached	for (ev = ein ; ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))) ; ev += ei) {					/* generate geometric distribution with richnes r and decay of ev */		fitdist = proportional_gs_distribution(ev,r);				// MAKE DISTRIBUTION		/* find the expected proportions of taxa with 0�x finds */		expect=expfindspart(fitdist,r,obsfnd,unqfnd,nspec);		for (i=0; i<unqfnd; ++i)	expect[i]=expect[i]/((float) r);		free_dvector(fitdist);		es[0] = calc_likelihood_exp(expect,obsfnd,obsrvd,unqfnd);		free_dvector(expect);		//Debugging line		if (ev<=(emin-SUPINC)) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);		if (es[0] >= bep[0]) {							// IF BETTER THAN BEST FIT			rs[0]=pbes = bep[0];			bep[0] = es[0];								// STORE FIT			bep[1] = ev;									// STORE SLOPE			}		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED			ev -= ei;									// STEP BACK ONE UNIT			ei *= -1;									// STEP BACKWARD TO FIND PEAK			es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP			}		else {									// NOT IMPROVING TRY SMALLER INCREMENT			ev -= ei;									// STEP BACK ONE UNIT			ei /= 10;									// SET SMALLER UNIT			}		es[2] = es[1];									// Store last 2 attempts to identify		es[1] = es[0];									// when the peak is past		}	/* if this richness is better than the last */	if (bep[0] >= brp[0]) {								// IF BETTER THAN BEST FIT		brp[2] = r;		for (i = 0 ; i < 2 ; i++)			brp[i] = bep[i];		ein=bep[1];										/* set initial decay to best decay found so far */ 		}	/* optimal richness is overshot */	else if ((rs[2]>rs[1] && (rs[0]>rs[1])) && (abs(ri)>1))	{		r-=ri;				/* step back one unit	 */		ri*=-1;				/* step backwards to peak */		rs[1]=rs[2];		/* set to prior ln L to ignore over step */		/* there is no point in letting r drop below ntaxa */		if ((r+ri)<ntaxa)	{			ri*=-1;				/* step backwards to peak */			ri=ri/2;			}		}	else	{		r-=ri;		if (abs(ri)>1)	ri/=2;		else			ri=0;		}	rs[2] = rs[1];	rs[1] = bep[0];//	printf("\nGS R=%f, ev=%1.10f, S=%4.14f\n",brp[2],brp[1],brp[0]);	}return brp;}/* CALCULATE LIKELIHOOD THAT A GIVEN DISTRIBUTION FITS THE GEOMETRIC SERIES. NEEDS:	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest	- ntaxa (size of array)	- nspec (sum of array)RETURNS:	- result[0]: log likelihood	- result[1]: slope of geometric series	- result[2]: optimal richness			COMMENTS:	- SUPINC, FITINC defined in header file***********************************************************************/double *x_fit_gs(int *empdist, int ntaxa, int nspec) {int i = 0, j;				/* LOOP VARIABLE	*/int r = 0;					/* LOOP RICHNESS	*/int ri;						/* richness increment	*/int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*//*int mxspc;					/* maximum number of possible finds to consider when calculating P */double ev = 0.000f;			/* LOOP SLOPE 		*/double emin = 1.0f;			/* min slope		*/double ein = 0.000f;		/* initial slope	*/double ei = 0.000f;			/* how much to increment ev in each loop							*/double es[3];				/* previous log likelihoods (cell number = num previous).		*/double rs[3];				/* previous log likelihoods (cell number = num previous).		*/double bep[3];				/* BEST ev PARAMETERS (DISTRIBUTION SLOPE) - return array format	*/double *brp;				/* BEST r PARAMETERS (DISTRIBUTION RICHNESS) - returned array	*/double *fitdist;			/* fit distribution												*/double *expect;				/* expected number of species with 0�max finds					*///double *exp;				/* expected number of species with 0�max finds					*/double *obsrvd;				/* observed number of species with 0�max finds					*/long *obsfnd;				/* array of unique values for numbers of finds					*/double pbes;				/* PREVIOUS BEST SUPPORT (E VALUE)								*/brp=dvector(3);for (i = 0 ; i < 3 ; i++) rs[i] = 0.0f;for (i = 0 ; i < 3 ; i++) brp[i]= -1.0*DBL_MAX;/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */obsrvd=idhistogram(empdist,ntaxa);/* we want one plus the number of different sample sizes because some species have a sample size of 0 (probably) *///unqfnd=1+countuniqivector(empdist,ntaxa);j=maxiarray(empdist,ntaxa);for (i=1; i<=j; ++i)	if (empdist[i]>0)	++unqfnd;obsfnd=lvector(unqfnd);obsfnd[j=0]=0;for (i=ntaxa-1; i>=0; --i)	if (empdist[i]!=empdist[i-1])	obsfnd[++j]=empdist[i];/************************************************************************************************************************//* find a good seed value for the geometric, based on the slope that would go from max to 1 in ntaxa species */ein=pow(empdist[0],(((double) empdist[ntaxa-1])/((double) ntaxa)));ri=ntaxa/10;if (ri<2)	ri=2;// increment true richness until that fails to improve likelihoodfor (r=ntaxa; abs(ri)>0; r+=ri) {		obsrvd[0]=r-ntaxa;	pbes = 0.0f;	ei = (double) FITINC / 10;	// increment starting at 0.1 since we are almost never going to find slopes of 2+		for (i = 0 ; i < 3 ; i++) bep[i] = -1.0*DBL_MAX;	for (i = 0 ; i < 3 ; i++) es[i] = 0.0f;	// increment slope until that fails to improve likelihood or resolution limit reached	for (ev = ein ; ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))) ; ev += ei) {					/* generate geometric distribution with richnes r and decay of ev */		fitdist = proportional_gs_distribution(ev,r);				// MAKE DISTRIBUTION		/* find the expected proportions of taxa with 0�x finds */		expect=expfindspart(fitdist,r,obsfnd,unqfnd,nspec);		for (i=0; i<unqfnd; ++i)	expect[i]=expect[i]/((float) r);		free_dvector(fitdist);		es[0] = -1*dsumsqdiffs(expect,obsrvd,1+empdist[0]);		free_dvector(expect);		//Debugging line		if (ev<=(emin-SUPINC)) printf("\nDANGER: GS R=%d, ev=%f, S=%f ",r,ev,es[0]);		if (es[0] >= bep[0]) {							// IF BETTER THAN BEST FIT			rs[0]=pbes = bep[0];			bep[0] = es[0];								// STORE FIT			bep[1] = ev;									// STORE SLOPE			}		else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {// TOO FAR: THE PEAK HAS BEEN PASSED			ev -= ei;									// STEP BACK ONE UNIT			ei *= -1;									// STEP BACKWARD TO FIND PEAK			es[1] = es[2];								// SET PREVIOUS S TO IGNORE OVER STEP			}		else {									// NOT IMPROVING TRY SMALLER INCREMENT			ev -= ei;									// STEP BACK ONE UNIT			ei /= 10;									// SET SMALLER UNIT			}		es[2] = es[1];									// Store last 2 attempts to identify		es[1] = es[0];									// when the peak is past		}	/* if this richness is better than the last */	if (bep[0] >= brp[0]) {								// IF BETTER THAN BEST FIT		brp[2] = r;		for (i = 0 ; i < 2 ; i++)			brp[i] = bep[i];		ein=bep[1];										/* set initial decay to best decay found so far */ 		}	/* optimal richness is overshot */	else if ((rs[2]>rs[1] && (rs[0]>rs[1])) && (abs(ri)>1))	{		r-=ri;				/* step back one unit	 */		ri*=-1;				/* step backwards to peak */		rs[1]=rs[2];		/* set to prior ln L to ignore over step */		/* there is no point in letting r drop below ntaxa */		if ((r+ri)<ntaxa)	{			ri*=-1;				/* step backwards to peak */			ri=ri/2;			}		}	else	{		r-=ri;		if (abs(ri)>1)	ri/=2;		else			ri=0;		}	rs[2] = rs[1];	rs[1] = bep[0];//	printf("\nGS R=%f, ev=%1.10f, S=%4.14f\n",brp[2],brp[1],brp[0]);	}return brp;}