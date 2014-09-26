/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
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

#define fit_simulated_distributions
#include "dist_fit_sim.h"

int main(void) {

// KEY PARAMETER VARIABLES
int richness = 0;						// current richness
float slope = 0.0f;					// current slope
double coefficent = 0.0f;				// current lp coefficent
double exponent = 0.0f;					// current lp exponent
int median_octave = 0.0f;				// current ln median octave
double octave_mag = 0.0f;				// current ln octave magnitude
int num_stdev = NSTDDEV;				// current ln number of standard deviations
int oct_per_stdev = OSTDDEV;				// current ln number of octaves per standard deviations

// MAIN SIMULATION VARIABLES
double *SimDist;						// CONTAINS PROPORTIONAL ABUNDANCE
int *SimSamp;						// CONTAINS SIMULATED SAMPLE
int rep = 0;						// main loop var for main replicate loop
int SimSampRich = 0;					// number of taxa in original distribution
double SumSampRich = 0.0f;				// total number of taxa
int SimSampSize = 0;					// number of taxa in original distribution
double SumSampSize = 0.0f;				// total number of taxa
int sample_mult = 0;					// current sample size multiple
int min_sample = 0;					// minimum sample size
int distsumB[N_DIST+1];					// number of Best fits
int distsumG[N_DIST+1];					// number of Good fits
int distsumN[N_DIST+1];					// number of Not fits
int DFIT[N_DIST+1];					// number of Not fits
double daic[N_DIST+1];
double dbic[N_DIST+1];
double aic[N_DIST+1];					// Akaike Information Criterion support for each distribution.
double bic[N_DIST+1];					// Bayes Information Criterion support for each distribution.
double *GSRes, *ZMRes, *LNRes, *LSRes, *LPRes;	// distribution parameters for best distribtion
double fitln[5], fitlp[4];
double fitgs[3], fitzm[3], fitls[3];		// stores values for mean fit values output.
int n_fit;

// MISC VARIABLES
int i=0;							// loop var
char outputname[65];					// name of main output file

FILE	*fopen();	
FILE 	*Out_summary;
FILE 	*Out_fits;
FILE 	*Out_dist;
FILE 	*Out_samp;

// START ACTUAL PROGRAM
//**************************************************************************
srand((unsigned int)time((time_t *)NULL));		// Set random number generator seed to clock
//srand(1);

if (LN==1) DFIT[1]=1; else  DFIT[1]=0;
if (GS==1) DFIT[2]=1; else  DFIT[2]=0;
if (ZM==1) DFIT[3]=1; else  DFIT[3]=0;
if (LS==1) DFIT[4]=1; else  DFIT[4]=0;
if (LP==1) DFIT[5]=1; else  DFIT[5]=0;

// OUTPUT - SCREEN AND FILE
//**************************************************************************
printf("EVALUATION DISTRIBUTION FITS\n");
printf("v. %s %s\n", VERSION, DATE);
printf("by M. Kosnik (mkosnik@uchicago.edu)\n");

// SET FILE NAMES
printf("\nOutputfile names:");

strcpy(outputname,DISTRIB);
strcat(outputname,".summary.out");
Out_summary=fopen(outputname,"w");
printf("\n\t%s\tRun Summary",outputname);

strcpy(outputname,DISTRIB);
strcat(outputname,".fits.out");
Out_fits=fopen(outputname,"w");
printf("\n\t%s\tFit distributions for each rep",outputname);

strcpy(outputname,DISTRIB);
strcat(outputname,".distrib.out");
Out_dist=fopen(outputname,"w");
printf("\n\t%s\tProportional distributions",outputname);

strcpy(outputname,DISTRIB);
strcat(outputname,".samples.out");
Out_samp=fopen(outputname,"w");
printf("\n\t%s\tAbsolute Samples",outputname);

// SET FILE HEADER
if (DISTRIB=="LN") {
	printf("\n-------------------------------------------------------------------------------------------------------------------------------------------");
	printf("\nTrDstr\tTrOcMg\tTrMdOc\tTrNOct\tOct/StD");
	fprintf(Out_summary,"TrDstr\tTrOcMg\tTrMdOc\tTrNOct\ttrue");
	fprintf(Out_samp,"TrDstr\tTrOcMg\tTrMdOc\tTrNOct\tTrOct/TrStDev");
	fprintf(Out_fits,"TrDstr\tTrOcMg\tTrMdOc\tTrNOct\tTrOct/TrStDev");
} else if ( (DISTRIB=="GS") || (DISTRIB=="ZM") || (DISTRIB=="LS") ) {
	printf("\n-------------------------------------------------------------------------------------------------------------------");
	printf("\nDist\tSlope");
	fprintf(Out_summary,"Dist\tTrue.Slope");
	fprintf(Out_samp,"Dist\tTrue.Slope");
	fprintf(Out_fits,"Dist\tTrue.Slope");
} else if (DISTRIB=="LP") {
	printf("\n-------------------------------------------------------------------------------------------------------------------");
	printf("\nDist\tCoef\tExpon");
	fprintf(Out_summary,"Dist\tTrue.Coef\tTrue.Expon");
	fprintf(Out_samp,"Dist\tTrue.Coef\tTrue.Expon");
	fprintf(Out_fits,"Dist\tTrue.Coef\tTrue.Expon");
}

printf("\tTrueS\tSampN\tSampS\t NF");
if (LN==1) printf("\tbLN\tgLN\tnLN");
if (GS==1) printf("\tbGS\tgGS\tnGS");
if (ZM==1) printf("\tbZM\tgZM\tnZM");
if (LS==1) printf("\tbLS\tgLS\tnLS");
if (LP==1) printf("\tbLP\tgLP\tnLP");

if (DISTRIB=="LN") {
	printf("\n-------------------------------------------------------------------------------------------------------------------------------------------");
} else if ( (DISTRIB=="GS") || (DISTRIB=="ZM") || (DISTRIB=="LS") ) {
	printf("\n-------------------------------------------------------------------------------------------------------------------");
} else if (DISTRIB=="LP") {
	printf("\n-------------------------------------------------------------------------------------------------------------------");
}

fprintf(Out_summary,"\tTrRich\tsample.min\tsample.mult\tsample.size\tsample.richness\tNF");
if (LN==1) fprintf(Out_summary,"\taicLN\tbicLN");
if (GS==1) fprintf(Out_summary,"\taicGS\tbicGS");
if (ZM==1) fprintf(Out_summary,"\taicZM\tbicZM");
if (LS==1) fprintf(Out_summary,"\taicLS\tbicLS");
if (LP==1) fprintf(Out_summary,"\taicLP\tbicLP");
if (LN==1) fprintf(Out_summary,"\tA_LN\tln.support\tln.inc/oct\tln.st_oct\tln.num_oct\tln.richness");
if (GS==1) fprintf(Out_summary,"\tA_GS\tgs.support\tgs.slope\tgs.richness");
if (ZM==1) fprintf(Out_summary,"\tA_ZM\tzm.support\tzm.slope\tzm.richness");
if (LS==1) fprintf(Out_summary,"\tA_LS\tls.support\tls.slope\tls.richness");
if (LP==1) fprintf(Out_summary,"\tA_LP\tlp.support\tlp.coef\tlp.expon\tlp.richness");

fprintf(Out_fits,"\tTrue.Rich\tSamp.Size\tSamp.Rich");
if (LN==1) fprintf(Out_fits,"\tLN\tSupport\tFit.inc/oct\tFit.st_oct\tFit.num_oct\tFit.richness");
if (GS==1) fprintf(Out_fits,"\tGS\tSupport\tFit.slope\tFit.richness");
if (ZM==1) fprintf(Out_fits,"\tZM\tSupport\tFit.slope\tFit.richness");
if (LS==1) fprintf(Out_fits,"\tLS\tSupport\tFit.slope\tFit.richness");
if (LP==1) fprintf(Out_fits,"\tLP\tSupport\tFit.coef\tFit.expon\tFit.richness");

if (LN==1) fprintf(Out_fits,"\taic[ln]");
if (GS==1) fprintf(Out_fits,"\taic[gs]");
if (ZM==1) fprintf(Out_fits,"\taic[zm]");
if (LS==1) fprintf(Out_fits,"\taic[ls]");
if (LP==1) fprintf(Out_fits,"\taic[lp]");
if (LN==1) fprintf(Out_fits,"\tdaic[ln]");
if (GS==1) fprintf(Out_fits,"\tdaic[gs]");
if (ZM==1) fprintf(Out_fits,"\tdaic[zm]");
if (LS==1) fprintf(Out_fits,"\tdaic[ls]");
if (LP==1) fprintf(Out_fits,"\tdaic[lp]");
if (LN==1) fprintf(Out_fits,"\tbic[ln]");
if (GS==1) fprintf(Out_fits,"\tbic[gs]");
if (ZM==1) fprintf(Out_fits,"\tbic[zm]");
if (LS==1) fprintf(Out_fits,"\tbic[ls]");
if (LP==1) fprintf(Out_fits,"\tbic[lp]");
if (LN==1) fprintf(Out_fits,"\tdbic[ln]");
if (GS==1) fprintf(Out_fits,"\tdbic[gs]");
if (ZM==1) fprintf(Out_fits,"\tdbic[zm]");
if (LS==1) fprintf(Out_fits,"\tdbic[ls]");
if (LP==1) fprintf(Out_fits,"\tdbic[lp]");

fprintf(Out_samp,"\tSamMin\tSamMul\tTrueS");

// RUN THROUGH SIMULATION PARAMETERS 
//**************************************************************************
for (min_sample=SAMPLE_MIN ; min_sample<=SAMPLE_MAX ; min_sample+=SAMPLE_INC)
for (sample_mult=SMULT_MIN ; sample_mult<=SMULT_MAX ; sample_mult*=SMULT_INC)

/*
//LOG-NORMAL INCREMENT
for (octave_mag=OMAG_START ; octave_mag<=OMAG_STOP ; octave_mag+=OMAG_INC) 
for (median_octave=OCT_START ; median_octave<=OCT_STOP ; median_octave+=OCT_INC) 
*/
//GEOMETRIC / ZIPF-MANDELBROT / LOG-SERIES INCREMENT
for (slope=SLOPE_START ; slope<=SLOPE_STOP ; slope+=SLOPE_INC)
 
/*
//LOG-POWER INCREMENT
for (coefficent=COEF_START ; coefficent<=COEF_STOP ; coefficent+=COEF_INC) 
for (exponent=EXP_START ; exponent<=EXP_STOP ; exponent+=EXP_INC) 
*/

// SIMULATE 
//**************************************************************************
for (richness=RICH_START ; richness<=RICH_STOP ; richness+=RICH_INC) {
	
	// INITALIZE VARS FOR EACH SET OF PARAMETERS
	n_fit=0;
	SumSampRich=0.0f;
	SumSampSize=0.0f;
	for (i=0 ; i<=N_DIST ; i++) {				//initialize result array
		distsumB[i] = 0;
		distsumG[i] = 0;
		distsumN[i] = 0;
		daic[i] = 0;
		dbic[i] = 0;
	 }
	for (i=0 ; i<4 ; i++)						//initialize result array
		fitlp[i] = 0.00000f;
	for (i=0 ; i<5 ; i++)						//initialize result array
		fitln[i] = 0.00000f;
	for (i=0 ; i<3 ; i++) {						//initialize result array
		fitgs[i] = 0.00000f;
		fitzm[i] = 0.00000f;
		fitls[i] = 0.00000f;
	}

	// OUTPUT SIMULATED DISTRIBUTION PARAMETERS TO FILES
	if (DISTRIB=="LN") {
		printf("\nLN\t%3.2f\t%2d\t%d\t%d",octave_mag,median_octave,num_stdev,oct_per_stdev);
		fprintf(Out_summary,"\nLN\t%f\t%d\t%d\t%d",octave_mag,median_octave,num_stdev,oct_per_stdev);
	//	fprintf(Out_fits,"\nLN\t%f\t%d\t%d\t%d",octave_mag,median_octave,num_stdev,oct_per_stdev);
		fprintf(Out_dist,"\nLN\t%f\t%d\t%d\t%d",octave_mag,median_octave,num_stdev,oct_per_stdev);
	} else if (DISTRIB=="GS") {
		printf("\nGS\t%3.2f",slope);
		fprintf(Out_summary,"\nGS\t%3.2f",slope);
	//	fprintf(Out_fits,"\nGS\t%3.2f",slope);
		fprintf(Out_dist,"\nGS\t%3.2f",slope);
	} else if (DISTRIB=="ZM") {
		printf("\nZM\t%3.2f",slope);
		fprintf(Out_summary,"\nZM\t%3.2f",slope);
	//	fprintf(Out_fits,"\nZM\t%3.2f",slope);
		fprintf(Out_dist,"\nZM\t%3.2f",slope);
	} else if (DISTRIB=="LS") {
		printf("\nLS\t%3.2f",slope);
		fprintf(Out_summary,"\nLS\t%3.2f",slope);
	//	fprintf(Out_fits,"\nLS\t%3.2f",slope);
		fprintf(Out_dist,"\nLS\t%3.2f",slope);
	} else if (DISTRIB=="LP") {
		printf("\nLP\t%3.2f\t%3.2f",coefficent,exponent);
		fprintf(Out_summary,"\nLP\t%f\t%f",coefficent,exponent);
	//	fprintf(Out_fits,"\nLP\t%f\t%f",coefficent,exponent);
		fprintf(Out_dist,"\nLP\t%f\t%f",coefficent,exponent);
	}
	printf("\t%3d",richness);
	fprintf(Out_summary,"\t%d",richness);
	fprintf(Out_summary,"\t%d\t%d",min_sample, sample_mult);
//	fprintf(Out_fits,"\t%d",richness);
	fprintf(Out_dist,"\t%d",richness);
	
	// MAKE SIMULATED DISTRIBUTION
	if (DISTRIB=="LN") {
		SimDist = proportional_ln_distribution(num_stdev, oct_per_stdev, median_octave, octave_mag, richness);
	} else if (DISTRIB=="GS") {
		SimDist = proportional_gs_distribution(slope,richness);
	} else if (DISTRIB=="ZM") {
		SimDist = proportional_zm_distribution(slope,richness);
	} else if (DISTRIB=="LS") {
//		SimDist = proportional_ls_distribution(slope,richness);
	} else if (DISTRIB=="LP") {
		SimDist = proportional_lp_distribution(coefficent,exponent,richness);
	} 
	
	// PRINT PROPORTIONAL DISTRIBUTION TO FILE FOR REFERENCE
	for (i=0 ; i<richness ; i++)
		fprintf(Out_dist,"\t%f",SimDist[i]);	
	fflush(stdout);		// THIS MAKES SURE THAT THE CURRENT DISTRIBUTION HITS THE FILE
	
	// DO REPLICATE SAMPLES OF SIMULATED DISTRIBUTION
	//**************************************************************************
	for (rep=0 ; rep<REPS ; rep++) {

//		printf("A"); fflush(stdout);  // DEBUG

		// SAMPLE DISTRIBUTION
		SimSamp = sample_sizexrich(SimDist, richness, min_sample, sample_mult);		// SAMPLE THE SIMULATED DISTRIBUTION
		SimSamp = ishellsort_dec(SimSamp, richness);						// sort abundance vector

		// PRINT SAMPLED DISTRIBUTION TO FILE FOR REFERENCE
		if (DISTRIB=="LN") {
			fprintf(Out_samp,"\nLN\t%f\t%d\t%d\t%d\t%d",octave_mag,median_octave,num_stdev,oct_per_stdev,richness);
		} else if (DISTRIB=="GS") {
			fprintf(Out_samp,"\nGS\t%3.2f\t%d",slope,richness);
		} else if (DISTRIB=="ZM") {
			fprintf(Out_samp,"\nZM\t%3.2f\t%d",slope,richness);
		} else if (DISTRIB=="LS") {
			fprintf(Out_samp,"\nLS\t%3.2f\t%d",slope,richness);
		} else if (DISTRIB=="LP") {
			fprintf(Out_samp,"\nLP\t%f\t%f\t%d",coefficent,exponent,richness);
		}
		fprintf(Out_samp,"\t%d\t%d",min_sample,sample_mult);
		for (i=0 ; i<richness ; i++)
			fprintf(Out_samp,"\t%d",SimSamp[i]);
	
		// CACULATE SAMPLED RICHNESS AND SAMPLE SIZE
		SimSampRich = 0;
		SimSampSize = 0;
		for (i=0 ; i<richness ; i++) {
			SimSampSize += SimSamp[i];
			if (SimSamp[i]>0) SimSampRich++;
		}
		SumSampRich += SimSampRich;
		SumSampSize += SimSampSize;

//		printf("B"); fflush(stdout);  // DEBUG

//		fprintf(Out_fits,"\n\t"); //PRINT THIS SO THAT THE ROW LENGTHS REMAIN MANAGEABLE
		if (DISTRIB=="LN") fprintf(Out_fits,"\nLN\t%f\t%d\t%d\t%d",octave_mag,median_octave,num_stdev,oct_per_stdev);
		else if (DISTRIB=="GS") fprintf(Out_fits,"\nGS\t%3.2f",slope);
		else if (DISTRIB=="ZM") fprintf(Out_fits,"\nZM\t%3.2f",slope);
		else if (DISTRIB=="LS") fprintf(Out_fits,"\nLS\t%3.2f",slope);
		else if (DISTRIB=="LP") fprintf(Out_fits,"\nLP\t%f\t%f",coefficent,exponent);
		fprintf(Out_fits,"\t%d\t%d\t%d",richness,SimSampSize,SimSampRich);

		// FIT DISTRIBUTION (IF THERE IS SUFFICENT TAXA / SPECIMENS)
		//**************************************************************************
		if (((SumSampSize/SimSampRich) >= MINMULT) && (SimSampRich >= MINOCCRS)) {
	
//			printf("C"); fflush(stdout);  // DEBUG

			n_fit++;
	
			// CALCULATE BEST FITS
			if (LN==1) LNRes = fit_ln(SimSamp, SimSampRich, SimSampSize);
			if (GS==1) GSRes = fit_gs(SimSamp, SimSampRich, SimSampSize);
			if (ZM==1) ZMRes = fit_zm(SimSamp, SimSampRich, SimSampSize);
			if (LS==1) LSRes = fit_ls(SimSamp, SimSampRich, SimSampSize);
			if (LP==1) LPRes = fit_lp(SimSamp, SimSampRich, SimSampSize);
	
//			printf("D"); fflush(stdout);  // DEBUG
				
			// SUM PARAMETERS FOR AVERAGES
			if (LN==1) for (i=0 ; i<5 ; i++) fitln[i]+=LNRes[i];
			if (GS==1) for (i=0 ; i<3 ; i++) fitgs[i]+=GSRes[i];
			if (ZM==1) for (i=0 ; i<3 ; i++) fitzm[i]+=ZMRes[i];
			if (LS==1) for (i=0 ; i<3 ; i++) fitls[i]+=LSRes[i];
			if (LP==1) for (i=0 ; i<4 ; i++) fitlp[i]+=LPRes[i];
	
//			printf("E"); fflush(stdout);  // DEBUG
			
			// OUTPUT FITS
			if (LN==1) fprintf(Out_fits,"\tLN\t%f\t%f\t%d\t%d\t%d",LNRes[0],LNRes[1],(int) LNRes[2],(int) LNRes[3],(int) LNRes[4]);
			if (GS==1) fprintf(Out_fits,"\tGS\t%f\t%f\t%d",GSRes[0],GSRes[1],(int) GSRes[2]);
			if (ZM==1) fprintf(Out_fits,"\tZM\t%f\t%f\t%d",ZMRes[0],ZMRes[1],(int) ZMRes[2]);
			if (LS==1) fprintf(Out_fits,"\tLS\t%f\t%f\t%d",LSRes[0],LSRes[1],(int) LSRes[2]);
			if (LP==1) fprintf(Out_fits,"\tLP\t%f\t%f\t%f\t%d",LPRes[0],LPRes[1],LPRes[2],(int) LPRes[3]);

			// CALCULATE AICc
			if (LN==1) aic[1] = calc_aic_c( 3, LNRes[0], SimSampRich);
			if (GS==1) aic[2] = calc_aic_c( 2, GSRes[0], SimSampRich);
			if (ZM==1) aic[3] = calc_aic_c( 2, ZMRes[0], SimSampRich);
			if (LS==1) aic[4] = calc_aic_c( 2, LSRes[0], SimSampRich);
			if (LP==1) aic[5] = calc_aic_c( 3, LPRes[0], SimSampRich);

			// DETERMINE BEST AICc
			aic[0] = DBL_MAX;
			for (i=1 ; i<=N_DIST ; i++) 
				if ((DFIT[i]==1) && (aic[0]>aic[i])) aic[0] = aic[i];

			// CALCULATE BIC
			if (LN==1) bic[1] = calc_bic( 3, LNRes[0], SimSampRich);
			if (GS==1) bic[2] = calc_bic( 2, GSRes[0], SimSampRich);
			if (ZM==1) bic[3] = calc_bic( 2, ZMRes[0], SimSampRich);
			if (LS==1) bic[4] = calc_bic( 2, LSRes[0], SimSampRich);
			if (LP==1) bic[5] = calc_bic( 3, LPRes[0], SimSampRich);

			// DETERMINE BEST BIC
			bic[0] = DBL_MAX;
			for (i=1 ; i<=N_DIST ; i++) 
				if ((DFIT[i]==1) && (bic[0]>bic[i])) bic[0] = bic[i];

			// CLEAN UP MEMORY
			if (LN==1) free_dvector(LNRes);
			if (GS==1) free_dvector(GSRes);
			if (ZM==1) free_dvector(ZMRes);
			if (LS==1) free_dvector(LSRes);
			if (LP==1) free_dvector(LPRes);
			
			// OUTPUT INFORMATION CRITERIA TO FIT FILE.
//			fprintf(Out_fits,"\ta0 %f",aic[0]);
//			fprintf(Out_fits,"\tb0 %f",bic[0]);
			for (i=1 ; i<=N_DIST ; i++) {
				if (DFIT[i]==1) fprintf(Out_fits,"\t%f",aic[i]);
				if (DFIT[i]==1) fprintf(Out_fits,"\t%f",aic[i]-aic[0]);
				if (DFIT[i]==1) fprintf(Out_fits,"\t%f",bic[i]);
				if (DFIT[i]==1) fprintf(Out_fits,"\t%f",bic[i]-bic[0]);
			}

			// TABULATE SUMMARY INFORMATION
			for (i=1 ; i<=N_DIST ; i++) {		// SUM Delta IC's
				if (DFIT[i]==1) daic[i]+=aic[i]-aic[0];
				if (DFIT[i]==1) dbic[i]+=bic[i]-bic[0];
			}
			for (i=1 ; i<=N_DIST ; i++) {		// COUNT HOW GOOD
				if ((DFIT[i]==1)&&(aic[i] < (aic[0]+AIC_GOOD))) {
					distsumB[i]++;
				} else if ((DFIT[i]==1)&&(aic[i] < (aic[0]+AIC_POOR))) {
					distsumG[i]++;
				} else if (DFIT[i]==1) {
					distsumN[i]++;
				}
			}
	
		// OUTPUT RESULTS - IF THERE WAS NOT SUIFFICENT TAXA / SPECIMENS TO FIT THEN...
		//**************************************************************************
		} else {
			distsumB[0]++;
			if (LN==1) fprintf(Out_fits,"\tLN\t.\t.\t.\t.\t.");
			if (GS==1) fprintf(Out_fits,"\tGS\t.\t.\t.");
			if (ZM==1) fprintf(Out_fits,"\tZM\t.\t.\t.");
			if (LS==1) fprintf(Out_fits,"\tLS\t.\t.\t.");
			if (LP==1) fprintf(Out_fits,"\tLP\t.\t.\t.\t.");
			fprintf(Out_fits,"\tNF\t.");
		}
		free_ivector(SimSamp);

//		printf("G"); fflush(stdout);  // DEBUG
 
 		// progress report to screen
		if (PROGRESS>0) print_progress(REPS, rep, PROGRESS);
	}

	//OUTPUT RESULTS - END OF REPS SUMMARY
	//**************************************************************************
	printf("\t%5.1f",SumSampSize/(float) REPS);						// output average sample size
	fprintf(Out_summary,"\t%5.1f",SumSampSize/(float) REPS);
	printf("\t%5.1f",SumSampRich/(float) REPS);						// output average sampled richness
	fprintf(Out_summary,"\t%5.1f",SumSampRich/(float) REPS);

	// number not fit (not enough taxa)
	printf("\t%3d",distsumB[0]);
	fprintf(Out_summary,"\t%d",distsumB[0]);
	
	// OUTPUT THE NUMBER OF TIMES EACH OPTION WON...
	if (LN==1) printf("\t%d\t%d\t%d",distsumB[1],distsumG[1],distsumN[1]);
	if (GS==1) printf("\t%d\t%d\t%d",distsumB[2],distsumG[2],distsumN[2]);
	if (ZM==1) printf("\t%d\t%d\t%d",distsumB[3],distsumG[3],distsumN[3]);
	if (LS==1) printf("\t%d\t%d\t%d",distsumB[4],distsumG[4],distsumN[4]);
	if (LP==1) printf("\t%d\t%d\t%d",distsumB[5],distsumG[5],distsumN[5]);

	for (i=1 ; i<=N_DIST ; i++) {
		if (DFIT[i]==1) fprintf(Out_summary,"\t%3.2f",daic[i]/(REPS-distsumB[0]));
		if (DFIT[i]==1) fprintf(Out_summary,"\t%3.2f",dbic[i]/(REPS-distsumB[0]));
	}

	// OUTPUT MEAN FIT PARAMETERS
	if (n_fit>0) {
		if (LN==1) {
			fprintf(Out_summary,"\tLN");
			for (i=0 ; i<5 ; i++)				//output result array
				fprintf(Out_summary,"\t%f",fitln[i]/(float) n_fit);
		}
	
		// OUTPUT MEAN GS PARAMETERS
		if (GS==1) {
			fprintf(Out_summary,"\tGS");
			for (i=0 ; i<3 ; i++)				//output result array
				fprintf(Out_summary,"\t%f",fitgs[i]/(float) n_fit);
		}
		
		// OUTPUT MEAN ZM PARAMETERS
		if (ZM==1) {
			fprintf(Out_summary,"\tZM");
			for (i=0 ; i<3 ; i++)				//output result array
				fprintf(Out_summary,"\t%f",fitzm[i]/(float) n_fit);
		}
		
		// OUTPUT MEAN LS PARAMETERS
		if (LS==1) {
			fprintf(Out_summary,"\tLS");
			for (i=0 ; i<3 ; i++)				//output result array
				fprintf(Out_summary,"\t%f",fitls[i]/(float) n_fit);
		}
		
		// OUTPUT MEAN LP PARAMETERS
		if (LP==1) {
			fprintf(Out_summary,"\tLP");
			for (i=0 ; i<4 ; i++)				//output result array
				fprintf(Out_summary,"\t%f",fitlp[i]/(float) n_fit);
		}
	} else {
		fprintf(Out_summary,"\t\t\t\t");
	}
	
	free_dvector(SimDist);
	fflush(stdout);	//MAKES SURE THAT ALL THE INFO FROM A PARAMETER SET GET TO FILE
}

//CLEAN UP AND FINISH
//**************************************************************************
fclose(Out_summary);
fclose(Out_fits);
fclose(Out_samp);
fclose(Out_dist);
printf("\n\nDONE!\n");
return 0;
}