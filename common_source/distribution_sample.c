/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
/*	Matthew Kosnik: mkosnik@uchicago.edu
/*
/*	This file is copyright (C) 2001 Matthew Kosnik
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

#define distribution_sample
#include "distribution_sample.h"

/* returns a vector containing taxa sampled from a proportional vector
/*
/* Pro_Dist = proportional abundance array
/* S = richness (length of array)
/* samplesize = how many samples to take from the distribution
***********************************************************************/
int *sample_fixedsize(double *ProDist, int S, int samplesize) {
	
	double *CumDist;					// Cumulative distribution
	int *Sample;						// Sampled distribution
	int i, j;							// loop variable
	double x;							// randomly chosen taxon
	
	CumDist = dvector(S);
	Sample = ivector(S);
	
	// CREATE CUMULATIVE DISTRIBUTION
	CumDist[0]=ProDist[0];
	for (i=1 ; i<S ; i++) {
		CumDist[i] = ProDist[i] + CumDist[i-1];
		Sample[i] = 0;
	}
	
	// SAMPLE DISTRIBUTION
	for (i=0 ; i<samplesize ; i++) {
		x = (double) rand() / ((double) RAND_MAX + 1);
		for (j=0 ; j<S ; j++) {
			if (x < CumDist[j]) {
				Sample[j]++;
				j=S+100;
			}
		}
	}
	
	free_dvector(CumDist);
	return Sample;
}

/* returns a vector containing taxa sampled from a proportional vector
    This function samples a proportional vector to a minimum size.
    Then increases the sample size such that the sample size is (x) times as large as the richness. 
 /* Pro_Dist = proportional abundance array
 /* S = richness (length of array)
 /* minsize = how many samples to take from the distribution
 /* multsize = sample size = multsize * number of taxa sampled
 ***********************************************************************/
 int *sample_sizexrich(double *ProDist, int trueS, int minsize, int multsize) {

     double *CumDist;					// Cumulative distribution
     int *Sample;						// Sampled distribution
     int samplesize=minsize;
     int sampS = 0;
     int i = 0, j = 0;					// loop variable
     double x = 0.0f;					// randomly chosen taxon

     CumDist = dvector(trueS);
     Sample = ivector(trueS);

     // CREATE CUMULATIVE DISTRIBUTION
     CumDist[0]=ProDist[0];
     for (i=1 ; i<trueS ; i++) {
	   CumDist[i] = ProDist[i] + CumDist[i-1];
	   Sample[i] = 0;
     }

     // SAMPLE DISTRIBUTION
     for (i=0 ; i < samplesize ; i++) {
	   x = (double) rand() / ((double) RAND_MAX + 1);
	   for (j=0 ; j<trueS ; j++) {
		 if (x < CumDist[j]) {
		     Sample[j]++;
		     j=trueS+100;
		 }
	   }
	   if ((i == (samplesize-1)) && (multsize!=0)) {		//set new sample size
		 sampS = 0;
		 for (j=0 ; j<trueS ; j++)
		     if (Sample[j]>0) sampS++;
		 samplesize = sampS * multsize;
	   }
     }

     free_dvector(CumDist);
     return Sample;
 }
 