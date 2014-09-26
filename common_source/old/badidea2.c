/* FIND THE MOST-LIKELY FAUX-LOG-NORMAL SERIES, VARYING INTIAL AND MODAL SLOPE AND THE MODE. NEEDS:	- empdist (array with sorted absolute abundance)... MUST be sorted from highest to lowest	- ntaxa (size of array)	- nspec (sum of array)RETURNS:	- result[0]: log likelihood	- result[1]: initial slope of shifting geometric series	- result[2]: modal slope of shifting geometric series	- result[3]: optimal richness				- result[4]: mode (in #taxa from median taxon rank)			COMMENTS:	- SUPINC, FITINC defined in header file*****************************************************************************************************/double *poi_fit_fln_3b(int *empdist, int ntaxa, int nspec, double ge, int gr){int i = 0;					/* LOOP VARIABLE	*/int r = 0;					/* LOOP RICHNESS	*/int ri;						/* richness increment													*/int rin;					/* initial richness to use in each search (begins as ntaxa)				*/double mode;				/* mode of log-normal, given as a taxon number							*/int mi=((double) ntaxa)/4;	/* how much to increment the mode										*/double imode;				/* initial mode (set to the median of taxa ranks)						*/int unqfnd=0;				/* number of different counts (=nspec if all counts unique, =1 if all species have X finds )*/double dev = 0.000f;		/* LOOP SHIFT IN SLOPE 													*/double dmin = 1.0f;			/* min change in ev slope												*/double din = 0.000f;		/* initial change ev slope												*/double di = 0.000f;			/* how much to increment change in ev in each loop						*/double ev = 0.000f;			/* LOOP SLOPE 															*/double emin = 0.75f;		/* min modal slope														*/double ein = 0.000f;		/* initial slope														*/double ei = 0.000f;			/* how much to increment ev in each loop								*/double ds[3];				/* previous initial decay log likelihoods (cell number = num previous).	*/double es[3];				/* previous modal decay log likelihoods (cell number = num previous).	*/double rs[3];				/* previous richness log likelihoods (cell number = num previous).		*/double ms[3];				/* previous mode log likelihoods (cell number = num previous).			*/double bdp[3];				/* BEST initial decay parameters - return array format 					*/double bep[4];				/* BEST modal decay parameters - return array format 					*/double brp[5];				/* BEST richness parameters - return array								*/double *bmp;				/* Best mode parameters - return array format							*/double *fitdist;			/* fit distribution														*/double *expect;				/* expected number of species with 0�max finds							*/double *obsrvd;				/* observed number of species with 0�max finds							*/double pbds;				/* previous best support for initial decay								*/double pbes;				/* previous best support for modal decay								*/double pbrs;				/* previous best support for richness									*/double pbms;				/* previous best support for mode location								*/bmp=dvector(5);for (i=0; i<5; i++) bmp[i]= -1.0*DBL_MAX;/* we need arrays giving both the number of taxa with x finds and all x finds with 1+ species */obsrvd=idhistogram(empdist,ntaxa);/* use previously calculated geometric for seed values */rin=chao2(empdist,ntaxa);din=ein = ge;/*for (r = ntaxa; ((pbrs == 0.0f) || (brp[0] > (pbrs + SUPINC))); r += ri) {	*/pbms=0.0f;/* adjust mode until that fails to improve likelihood	*/for (mode=0; abs(mi)>0 && ((pbms == 0.0f) || (bmp[0] > (pbms + SUPINC))); mode+=mi)	{	pbrs=0.0f;	for (i=0; i<4; i++) brp[i]= -1.0*DBL_MAX;	for (i=0; i<3; ++i)	rs[i]=0;	/* adjust true richness until that fails to improve likelihood	*/	ri=ntaxa/2;	if (ri<2)	ri=2;	for (r=rin; abs(ri)>0 && r>=ntaxa; r+=ri)	{		/* make sure that mode starts between beginning and end */		while (abs(mode/2)>r)	++r;				if (r%2==1)	imode = r/2;		else		imode = (r/2)-0.5;		obsrvd[0]=r-ntaxa;		pbes = 0.0f;		ei = -1.0 * (double) FITINC / 10;		/* increment starting at 0.1 since we are almost never going to find slopes of 2+	*/				for (i=0; i<4; i++) bep[i] = -1.0*DBL_MAX;		for (i=0; i<3; i++) es[i]=0.0f;		/* increment slope until that fails to improve likelihood or resolution limit reached */		for (ev=ein; ev>=emin && ((pbes == 0.0f) || (bep[0] > (pbes + SUPINC))); ev += ei) {							for (i=0; i<3; i++) bdp[i]=-1.0*DBL_MAX;			for (i=0; i<3; i++) ds[i]=0.0f;			di = (double) FITINC / 10;			/*Debugging line */	/*		if (ev<=emin) printf("\nDANGER: fln R=%d, ev=%f, S=%f ",r,ev,es[0]);	*/			pbds = 0.0f;			for (dev=din; ((pbds==0.0f) || (bdp[0]>(pbds+SUPINC))); dev+=di)	{				/* generate geometric distribution with richnes r and decay of ev */				fitdist = proportional_fln_distribution(ev,dev,r,imode+mode);				/* MAKE DISTRIBUTION */				/* find the expected proportions of taxa with 0�x finds */				expect=expfinds(fitdist,r,2*empdist[0],nspec);				free_dvector(fitdist);				ds[0] = lnPoisson_vector_part(expect,obsrvd,empdist[0],1);				free_dvector(expect);				if (ds[0] >= bdp[0]) {						/* IF BETTER THAN BEST FIT */					pbds = bdp[0];									/* save last best ssq for evenness */					es[0]=bdp[0]=ds[0];								/* STORE FIT */					bdp[1] = dev;									/* STORE shift in slope	*/					}				else if ((ds[2] > ds[1]) && (ds[0] > ds[1]) && ((dev-di)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */					dev -= di;									/* STEP BACK ONE UNIT */					di *= -1;									/* STEP BACKWARD TO FIND PEAK */					ds[1] = ds[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */					}				else {											/* NOT IMPROVING TRY SMALLER INCREMENT */					dev -= di;									/* STEP BACK ONE UNIT */					di /= 10;									/* SET SMALLER UNIT */					}				/* it might go on and on forever with miniscule increases when it is a very poor fit */				if ((ds[0]>=ds[1] && ds[0]<(ds[1]+SUPINC)) && (ds[1]>=ds[2] && ds[1]<(ds[2]+SUPINC)))	break;				ds[2] = ds[1];									/* Store last 2 attempts to identify */				ds[1] = ds[0];									/* when the peak is past */				}			if (es[0] >= bep[0]) {							/* IF BETTER THAN BEST FIT */				pbes = bep[0];								/* save last best ssq for evenness */				rs[0] = bep[0] = es[0];						/* STORE FIT */				din = bep[1] = bdp[1];				bep[2] = ev;								/* STORE SLOPE */				}			else if ((es[2] > es[1]) && (es[0] > es[1]) && ((ev-ei)>emin)) {/* TOO FAR: THE PEAK HAS BEEN PASSED */				ev -= ei;									/* STEP BACK ONE UNIT */				ei *= -1;									/* STEP BACKWARD TO FIND PEAK */				es[1] = es[2];								/* SET PREVIOUS S TO IGNORE OVER STEP */				}			else {											/* NOT IMPROVING TRY SMALLER INCREMENT */				ev -= ei;									/* STEP BACK ONE UNIT */				ei /= 10;									/* SET SMALLER UNIT */				}			while ((ev+ei)<emin && ei>0.00001)	ei /= 10;						/* it might go on and on forever with miniscule increases when it is a very poor fit */			if ((es[0]>=es[1] && es[0]<(es[1]+SUPINC)) && (es[1]>=es[2] && es[1]<(es[2]+SUPINC)))	break;			es[2] = es[1];									/* Store last 2 attempts to identify */			es[1] = es[0];									/* when the peak is past */			}		/* if this richness is better than the last */		if (bep[0] >= brp[0]) {								/* IF BETTER THAN BEST FIT */			brp[4] = r;			for (i=0; i<3; i++)				brp[i] = bep[i];			ms[0]=brp[0];			din=bep[1];										/* set initial decay to best decay found so far */ 			ein=bep[2];										/* set initial decay to best decay found so far */ 			if (ein<1.0)	ein+=0.1;			}		/* optimal richness is overshot */		else if ((rs[2]>rs[1] && (rs[0]>rs[1])) && (abs(ri)>1))	{			r-=ri;				/* step back one unit	 */			ri*=-1;				/* step backwards to peak */			rs[1]=rs[2];		/* set to prior ln L to ignore over step */			}		else	{			r-=ri;			if (abs(ri)>1)	ri/=2;			else			ri=0;			}		/* it might go on and on forever with miniscule increases when it is a very poor fit */		if ((rs[0]>=rs[1] && rs[0]<(rs[1]+SUPINC)) && (rs[1]>=rs[2] && rs[1]<(rs[2]+SUPINC)))	ri=0;		rs[2] = rs[1];		rs[1] = bep[0];		}	/* if this mode is better than the last */	if (brp[0] >= bmp[0]) {								/* IF BETTER THAN BEST FIT */		for (i=0; i<5; i++)			bmp[i] = brp[i];		bmp[3] = mode;		din=bmp[1];										/* set initial initial slope to best initial slope found so far */ 		ein=bmp[2];										/* set initial modal decay to best modal decay found so far 	*/ 		rin=bmp[4];										/* set initial richness to best richness found so far			*/		}	/* optimal richness is overshot */	else if ((ms[2]>ms[1] && (ms[0]>ms[1])) && (abs(mi)>1))	{		mode-=mi;				/* step back one unit	 */		mi*=-1;					/* step backwards to peak */		ms[1]=ms[2];		/* set to prior ln L to ignore over step */		}	else	{		mode-=mi;		if (abs(mi)>1)	mi/=2;		else			mi=0;		}	/* it might go on and on forever with miniscule increases when it is a very poor fit */	if ((ms[0]>=ms[1] && ms[0]<(ms[1]+SUPINC)) && (ms[1]>=ms[2] && ms[1]<(ms[2]+SUPINC)))	mi=0;	ms[2] = ms[1];	ms[1] = brp[0];		}free_dvector(obsrvd);return bmp;}