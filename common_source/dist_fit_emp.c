/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*	Matthew Kosnik: mkosnik@uchicago.edu
/*
/*	This file is copyright (C) 2002 Matthew Kosnik
/*
/*	This program is free software; you can redistribute it and/or modify it 
/*	under the terms of version 2 the GNU General Public License as published 
/*	by the Free Software Foundation.
/*
/*	This program is distributed in the hope that it will be useful, but WITHOUT
/*	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
/*	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
/*	more details.
/*
/*	To view a copy of the license go to:
/*	http://www.fsf.org/copyleft/gpl.html
/*	To receive a copy of the GNU General Public License write the Free Software
/* 	Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
/*
/*	Copies of this source code are available without cost from:
/*	http://geosci.uchicago.edu/paleo/csource/
/*	
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define fit_emperical_distributions
#include "dist_fit_emp.h"

int main(void) {

unsigned long **RawData;				/* raw data matrix */
unsigned long **s_data;					/* specimen data matrix (list by taxa abundance) */
unsigned long **o_data;					/* occurrence data matrix (list by taxa) */
unsigned long *list_nums;				/* list numbers */
unsigned long *taxa_nums;				/* taxon numbers */
unsigned long *list_occrs;				/* number of occurrences per list */
unsigned long *list_specs;				/* number of specimens per list */
int noccr=1;							/* number of occurrences */
int ntaxa=1;							/* number of taxa */
int nlists=1;							/* number of lists */
char input[NAME], output[NAME];			/* file names */

int distsumB[N_DIST+1];					// number of best fits
int distsumG[N_DIST+1];					// number of good fits
int distsumP[N_DIST+1];					// number of poor fits
int DFIT[N_DIST+1];						// WHICH DISTRIBUTIONS TO FIT

int	a, b, i, k, n;
int *Abund;								// abundance vector for current list
double *RelAbund;						// relative abundance vector for current list
double EvenKW;							// average deviation from uniform for current list.
double EvenBG;							// Buzas-Gibson evenness for current list
double *GSRes;							// fit Geometric series parameters
double *ZMRes;							// fit Ziph-Mandelbrot parameters
double *LNRes;							// fit log-normal parameters
double *LPRes;							// fit log-power parameters
double *LSRes;							// fit log-series parameters
double *BestDist;						// actual parameters from best lists.
double *aic;							// AIC valves for each distribution, minimum in 0
double *bic;							// BIC valves for each distribution, minimum in 0
FILE *fopen();
FILE *Outfile;
FILE *Outfile1;

// Set random number generator seed to clock.
srand((unsigned int)time((time_t *)NULL));

for (i=0 ; i<N_DIST+1 ; i++) {
	distsumB[i]=0;
	distsumG[i]=0;
	distsumP[i]=0;
}

if (LN==1) DFIT[1]=1; else DFIT[1]=0;
if (GS==1) DFIT[2]=1; else DFIT[2]=0;
if (ZM==1) DFIT[3]=1; else DFIT[3]=0;
if (LS==1) DFIT[4]=1; else DFIT[4]=0;
if (LP==1) DFIT[5]=1; else DFIT[5]=0;

// BASIC INFORMATION ABOUT PROGRAM TO SCREEN
//**************************************************************************
printf("DISTRIBUTION FIT");
printf("\nversion: %s  released: %s\n", VERSION, DATE);
printf("--------------------------------------------------------------------------------------------------------\n");
printf(" ... calculates evenness metrics and fits to distributions for lists of taxon abundance.\n");
printf(" ... is copyright 2002 under the terms of version 2 the GNU General Public License (it is free as in speech).\n");
printf(" ... the source code can found at http://geosci.uchicago.edu/paleo/csource/ (it is free as in beer).\n");
printf(" ... all the numbers in this intro are defined in the header file and can be changed.\n");
printf(" ... was written by M. Kosnik (mkosnik@uchicago.edu).\n");
printf("\nKnown limitations:\n");
printf("--------------------------------------------------------------------------------------------------------\n");
printf(" ... only lists averaging %d or more specimens per taxon and at least %d taxa will be fit.\n",MINMULT,MINOCCRS);
printf(" ... distribution parameters will be fit until the change is support is less than %f.\n",SUPINC);
printf(" ... evenness will be calculated for all samples.\n");
printf("\nInput file format:\n");
printf("--------------------------------------------------------------------------------------------------------\n");
printf(" ... tab delimited (%d column) text files only, with no headers.\n",C_TOTAL);
printf(" ... file should only contain interger numbers, tabs, and return characters.\n");
printf(" ... Blank spaces or non-numeral characters will cause the data file to not load.\n");
printf(" ... Memory is dynamically allocated so file length sould be limited only by your computer.\n");
printf(" column %d: list number (unsigned long)\n",C_LIST+1);
printf(" column %d: taxon number (unsigned long)\n",C_SPEC+1);
printf(" column %d: abundance (unsigned long greater than 0)\n",C_ABUN+1);

// GET INPUT FILE NAME AND READ DATA
//**************************************************************************
printf("\nEnter the name of the datafile (%d character limit): ",NAME);
scanf("%s",&input);
noccr = getfilelength(input, C_TOTAL);						/*calc file length (memory) */
RawData = ul_readdata(input, noccr, C_TOTAL);				/*read input*/
nlists = ul_countcolumn(RawData, noccr, C_LIST);  			/*number of lists*/
ntaxa = ul_countcolumn(RawData, noccr, C_SPEC);  			/*number of taxa*/

list_nums = ul_pullcolumn(RawData, noccr, C_LIST, C_LIST);	/*list numbers*/
taxa_nums = ul_pullcolumn(RawData, noccr, C_SPEC, C_SPEC);	/*taxon numbers*/

s_data = makerectmatrix(RawData, noccr, ntaxa, nlists, C_LIST, C_SPEC, C_ABUN);	/*specimen data*/	
o_data = ul_convert_data_to_pa(s_data, ntaxa, nlists);		/*occurance data*/

list_specs = ul_sumrows(s_data, ntaxa, nlists);				/* number of specimens per list */
list_occrs = ul_sumrows(o_data, ntaxa, nlists);				/* number of occrrences per list */
free_ulmatrix(o_data, nlists, ntaxa);

free_ulmatrix(RawData, noccr, C_TOTAL);

// OUTPUT SUMMARY OF INPUT FILE
//**************************************************************************
printf("\n\"%s\" contains:\n",input);
printf("---------------------------------------------------------------------------------------------------------\n");
printf(" %5d lists\n",nlists);
printf(" %5d taxa\n",ntaxa);
printf(" %5d occurrences (file length)\n",noccr);

// OUTPUT SUMMARY OF OUTPUT
//**************************************************************************
printf("\nThe program returns:\n");
printf("---------------------------------------------------------------------------------------------------------\n");
printf(" For all lists:\n");
printf(" - list number.\n");
printf(" - taxon richness (number of taxa in the list).\n");
printf(" - sample size (number of specimens in the list).\n");
printf(" - average log 10 deviation from uniform.\n");
printf(" - Buzas-Gibson (1969) evenness.\n");
printf(" For all lists averaging %d or more specimens per taxon and %d taxa:\n",MINMULT,MINOCCRS);
printf(" - log likelihood for each distribution.\n");
printf(" - delta AICc for each distribution.\n");

printf("------------------------------------------------");
if (LN==1) printf("------------------------------");
if (GS==1) printf("------------------------------");
if (ZM==1) printf("------------------------------");
if (LS==1) printf("------------------------------");
if (LP==1) printf("------------------------------");
printf("-\n");
printf("List      Taxa      Size       E[ALD]   E[BG]   ");
if (LN==1) printf("     L(LN)");
if (GS==1) printf("     L(GS)");
if (ZM==1) printf("     L(ZM)");
if (LS==1) printf("     L(LS)");
if (LP==1) printf("     L(LP)");
if (LN==1) printf("  AICc(LN)   BIC(LN)");
if (GS==1) printf("  AICc(GS)   BIC(GS)");
if (ZM==1) printf("  AICc(ZM)   BIC(ZM)");
if (LS==1) printf("  AICc(LS)   BIC(LS)");
if (LP==1) printf("  AICc(LP)   BIC(LP)");
printf("\n");
printf("------------------------------------------------");
if (LN==1) printf("------------------------------");
if (GS==1) printf("------------------------------");
if (ZM==1) printf("------------------------------");
if (LS==1) printf("------------------------------");
if (LP==1) printf("------------------------------");
printf("-\n");

// NAME FILE TO CONTAIN DISTRIBUTIONS (EMPERICAL AND FIT DISTRIBUTIONS)
strcpy(output,input);
strcat(output,".distributions");
Outfile1=fopen(output,"w");

// NAME FILE TO CONTAIN SUMMARY INFORMATION (LIST INFORMATION AND FIT INFORMATION)
strcpy(output,input);
strcat(output,".summary");
Outfile=fopen(output,"w");

// SUMMARY FILE HEADER ROW
fprintf(Outfile,"Locality\tTaxa\tSpecimens\tALD\tBG89\t");
if (LN==1) fprintf(Outfile,"ln.support\tln.mag\tln.startoct\tln.num_oct\tln.rich\t");
if (GS==1) fprintf(Outfile,"gs.support\tgs.slope\tgs.rich\t");
if (ZM==1) fprintf(Outfile,"zm.support\tzm.slope\tzm.rich\t");
if (LS==1) fprintf(Outfile,"ls.support\tls.slope\tls.rich\t");
if (LP==1) fprintf(Outfile,"lp.support\tlp.coef\tlp.exp\tlp.rich\t");
fprintf(Outfile,"\n");

// PERFORM ANALYSES FOR EACH LIST AND OUTPUT INFORMATION
//**************************************************************************
aic = dvector(N_DIST+1);			// aicc values for each distibution
bic = dvector(N_DIST+1);			// bic values for each distibution
for (a=0 ; a<nlists ; a++)	{
	
	// print basic list information to screen and file
	//**************************************************************************
	printf("%4ld\t%5ld\t%8ld\t",list_nums[a],list_occrs[a],list_specs[a]);
	fprintf(Outfile,"%ld\t%ld\t%ld\t",list_nums[a],list_occrs[a],list_specs[a]);
	fflush(stdout);
	
	// create abundance arrays
	//**************************************************************************
	Abund = ivector(list_occrs[a]);
	RelAbund = dvector(list_occrs[a]);
	i=0;
	for (b=0 ; b<ntaxa ; b++) {
		if (s_data[a][b]>0) {
			Abund[i] = s_data[a][b];
			RelAbund[i] = (float)s_data[a][b]/(float)list_specs[a];
			i++;
		}
	}
	/* sort abundance vector */
	Abund=ishellsort_dec(Abund, list_occrs[a]);

	// Calculate evenness metrics
	//**************************************************************************
	EvenKW = calc_kw_evenness(RelAbund,list_occrs[a]);			/* calculate Avg. Dev. from uniform*/
	EvenBG = calc_bg_evenness(RelAbund,list_occrs[a]);			/* calculate Buzas-Gibson evenness */
	printf("%4.3f\t%4.3f\t",EvenKW,EvenBG);						/* output evenness metrics */
	fprintf(Outfile,"%f\t%f\t",EvenKW,EvenBG);
	fflush(stdout);   
	free_dvector(RelAbund);

	// FIT DISTRIBUTIONS IF THE LIST HAS ENOUGH TAXA AND SPECIMENS
	//**************************************************************************
	if (((list_specs[a]/list_occrs[a]) >= MINMULT) && (list_occrs[a] >= MINOCCRS)){

	// CALCULATE BEST FIT DISTRIBUTION FOR EACH DISTRIBUTION
	//**************************************************************************
		
		// LOG-NORMAL
		if (LN==1) {
			LNRes = fit_ln(Abund, list_occrs[a], list_specs[a]);
			printf("%9.1f ",LNRes[0]);
			fprintf(Outfile,"%f\t%f\t%d\t%d\t%d\t",LNRes[0],LNRes[1],(int) LNRes[2],(int) LNRes[3], (int) LNRes[4]);
		}

		// GEOMETRIC SERIES
		if (GS==1) {
			GSRes = fit_gs(Abund, list_occrs[a], list_specs[a]);
			printf("%9.1f ",GSRes[0]);
			fprintf(Outfile,"%f\t%f\t%d\t",GSRes[0],GSRes[1],(int) GSRes[2]);
		}	

		// ZIPH-MANDELBROT
		if (ZM==1) {
			ZMRes = fit_zm(Abund, list_occrs[a], list_specs[a]);
			printf("%9.1f ",ZMRes[0]);
			fprintf(Outfile,"%f\t%f\t%d\t",ZMRes[0],ZMRes[1],(int) ZMRes[2]);
		}
		
		// LOG-SERIES
		if (LS==1) {
			LSRes = fit_ls(Abund, list_occrs[a], list_specs[a]);
			printf("%9.1f ",LSRes[0]);
			fprintf(Outfile,"%f\t%f\t%d\t",LSRes[0],LSRes[1],(int) LSRes[2]);
		}
		
		// LOG-POWER
		if (LP==1) {
			LPRes = fit_lp(Abund, list_occrs[a], list_specs[a]);
			printf("%9.1f ",LPRes[0]);
			fprintf(Outfile,"%f\t%f\t%f\t%d\t",LPRes[0],LPRes[1],LPRes[2],(int) LPRes[3]);
		}

	// CALCULATE / OUTPUT BEST FIT DISTRIBUTION OF EACH DISTRIBUTION
	//**************************************************************************
	
		// OUTPUT REAL DISTRIBUTION TO DISTRIBUTIONS FILE
		fprintf(Outfile1,"\nList %ld",list_nums[a]);
		fprintf(Outfile1,"\nReal");
		for (i=0 ; i<list_occrs[a] ; i++) fprintf(Outfile1,"\t%f", (float)Abund[i]/(float)list_specs[a]);

		// LOG-NORMAL
		if (LN==1) {
			BestDist = proportional_ln_distribution(NSTD, OSTD, LNRes[2], LNRes[1], LNRes[4]);
			fprintf(Outfile1,"\nLN");
			for (i=0 ; i<LNRes[4] ; i++) fprintf(Outfile1,"\t%f",BestDist[i]);
			free_dvector(BestDist);
		}			
		// GEOMETRIC SERIES
		if (GS==1) {
			BestDist = proportional_gs_distribution(GSRes[1],GSRes[2]);
			fprintf(Outfile1,"\nGS");
			for (i=0 ; i<GSRes[2] ; i++) fprintf(Outfile1,"\t%f",BestDist[i]);
			free_dvector(BestDist);
		}
		// ZIPH-MANDELBROT
		if (ZM==1) {
			BestDist = proportional_zm_distribution(ZMRes[1],ZMRes[2]);
			fprintf(Outfile1,"\nZM");
			for (i=0 ; i<ZMRes[2] ; i++) fprintf(Outfile1,"\t%f",BestDist[i]);
			free_dvector(BestDist);
		}
		// LOG-SERIES
		if (LS==1) {
// Function pulled by P. Wagner, not working well enough.
//			BestDist = proportional_ls_distribution(LSRes[1],LSRes[2]);
			fprintf(Outfile1,"\nLS");
			for (i=0 ; i<LSRes[2] ; i++) fprintf(Outfile1,"\t%f",BestDist[i]);
			free_dvector(BestDist);
		}
		// LOG-POWER
		if (LP==1) {
			BestDist = proportional_lp_distribution(LPRes[1],LPRes[2],LPRes[3]);
			fprintf(Outfile1,"\nLP");
			for (i=0 ; i<LPRes[3] ; i++) fprintf(Outfile1,"\t%f",BestDist[i]);
			free_dvector(BestDist);
		}
	// IC calculations and output
	//**************************************************************************

		// CALCULATE AICc
		if (LN==1) aic[1] = calc_aic_c( 3, LNRes[0], list_specs[a]);
		if (GS==1) aic[2] = calc_aic_c( 2, GSRes[0], list_specs[a]);
		if (ZM==1) aic[3] = calc_aic_c( 2, ZMRes[0], list_specs[a]);
		if (LS==1) aic[4] = calc_aic_c( 2, LSRes[0], list_specs[a]);
		if (LP==1) aic[5] = calc_aic_c( 3, LPRes[0], list_specs[a]);

		// DETERMINE BEST AICc
		aic[0] = DBL_MAX;
		for (i=1 ; i<=N_DIST ; i++) 
			if ((DFIT[i]==1) && (aic[0]>aic[i])) aic[0] = aic[i];

		// CALCULATE BIC
		if (LN==1) bic[1] = calc_bic( 3, LNRes[0], list_specs[a]);
		if (GS==1) bic[2] = calc_bic( 2, GSRes[0], list_specs[a]);
		if (ZM==1) bic[3] = calc_bic( 2, ZMRes[0], list_specs[a]);
		if (LS==1) bic[4] = calc_bic( 2, LSRes[0], list_specs[a]);
		if (LP==1) bic[5] = calc_bic( 3, LPRes[0], list_specs[a]);

		// DETERMINE BEST BIC
		bic[0] = DBL_MAX;
		for (i=1 ; i<=N_DIST ; i++) 
			if ((DFIT[i]==1) && (bic[0]>bic[i])) bic[0] = bic[i];

		// OUTPUT IC
		for (i=1 ; i<=N_DIST ; i++) 
			if (DFIT[i]==1) { 
				printf ("%10.1f",aic[i]-aic[0]);
				printf ("%10.1f",bic[i]-bic[0]);
			}
		fflush(stdout);
		
		// CLEAN UP MEMORY
		if (LN==1) free_dvector(LNRes);
		if (GS==1) free_dvector(GSRes);
		if (ZM==1) free_dvector(ZMRes);
		if (LS==1) free_dvector(LSRes);
		if (LP==1) free_dvector(LPRes);
		
	// TABULATE SUMMARY INFORMATION
	//**************************************************************************
		n=0;
		k=0;
		for (i=1 ; i< (N_DIST+1) ; i++) {
		
			if ((int) ((aic[i] / AIC_PREC +.5) * AIC_PREC) == (int) ((aic[0] / AIC_PREC +.5) * AIC_PREC)) {
				// if best
				distsumB[i]++;
				k = i;
				n++;
			} else if ((int) ((aic[i] / AIC_PREC +.5) * AIC_PREC) < (int) (((aic[0]+AIC_GOOD) / AIC_PREC +.5) * AIC_PREC)) {
				// if good
				distsumG[i]++;
				n++;
			} else {
				// if bad
				distsumP[i]++;
			}
		}
		if (n > 1) {
			// if no other distribution within good range then best is only good
			distsumG[k]++;
			distsumB[k]--;
		}
	}
	printf("\n");
	fprintf(Outfile,"\n");
	free_ivector(Abund);
}

// PRINT SUMMARY INFORMATION TO SCREEN
//**************************************************************************
printf("\n\nFit Summary");
printf("\n Best = delta AICc = 0");
//printf("\n Best = delta AICc = 0 & others >=%d",AIC_GOOD);
printf("\n Good = delta AICc < %d, (but not best)",AIC_GOOD);
printf("\n Not  = delta AICc >= %d",AIC_GOOD);
printf("\n------------------------");
printf("\nDist  Best   Good    Not");
if (LN==1) printf("\nLN %6d %6d %6d",distsumB[1],distsumG[1],distsumP[1]);
if (GS==1) printf("\nGS %6d %6d %6d",distsumB[2],distsumG[2],distsumP[2]);
if (ZM==1) printf("\nZM %6d %6d %6d",distsumB[3],distsumG[3],distsumP[3]);
if (LS==1) printf("\nLS %6d %6d %6d",distsumB[4],distsumG[4],distsumP[4]);
if (LP==1) printf("\nLP %6d %6d %6d",distsumB[5],distsumG[5],distsumP[5]);
printf("\n");

// CLEAN UP MEMORY / FILES
//**************************************************************************
fclose(Outfile);
fclose(Outfile1);

free_dvector(aic);
free_dvector(bic);
free_ulvector(list_nums);
free_ulvector(taxa_nums);
free_ulvector(list_specs);
free_ulvector(list_occrs);
free_ulmatrix(s_data, nlists, ntaxa);
return 0;
}