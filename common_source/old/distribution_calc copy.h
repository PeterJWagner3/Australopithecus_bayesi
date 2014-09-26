/*	Matthew Kosnik: mkosnik@uchicago.edu
/*
/*	This file is copyright (C) 2001 Peter J. Wagner & Matthew Kosnik
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

/***********************************************************************
v0.0 	2001.05.04 BY P.J.WAGNER III
	  CODE WRITTEN - CODEWARRIOR 6.1 / MacOS 9.04
v0.1	2001.12.20 BY M. KOSNIK
	  substantively rewritten...
	  also renamed functions
v0.2	2001.12.21 BY M. KOSNIK
            fixed bug in ln generation that was causing memory write at A[-2]
v0.3	2002.01.09 BY M. KOSNIK
			wagner code added:
			-	double* LNAbundCalc(int R, int num_stdev, int oct_per_stdev, double LS, int *Rich);
			-	double* normdistodd(int num_stdev, int oct_per_stdev, int starting_oct);
			-	int* OctaveRichness(int R, double *NormA, int num_stdev, int oct_per_stdev);
			wagner revised:
			-	to make it easier to understand code...
			-	some optimizations made
v0.4	2002.02.09 BY M. KOSNIK
			revised OctaveRichness
				Rich[a] = (int)( (double) NormA[a] * (double) R );
				changed to:
				Rich[a] = (int) ceil( (double) NormA[a] * (double) R );
				removed:
				Rich[0]=Rich[0]+(R-c);

 ***********************************************************************/

#ifdef distribution_calc
    #include <math.h>
    #include <float.h>
    #include <stdlib.h>
    #include <stdio.h>
    #include "memory.h"
    #include "sort.h"

	#define PI 3.141592654f				// constant pi
    #define EX 2.718281828f				// constant e


    double *proportional_gs_distribution(double M, int S); 
    double *proportional_zm_distribution(double M, int S); 
	double *proportional_lp_distribution(double C, double X, int S);
	double *proportional_ln_distribution(int num_stdev, int oct_per_stdev, int starting_octave, double mag_per_oct, int S);
	double *proportional_ln_distribution2(int num_stdev, int oct_per_stdev, int modal_octave, double mag_per_oct, int S, int *idistribution);
	double *proportional_fln_distribution(double M, double dM, int S, double mode);
	double *proportional_lgn_distribution (double mag, int trunc, double mode, int S);
	double *proportional_lgn_distributionN(double mag, double mode, double ocs, int S);
	double *proportional_lgn_distributionA(double mag, double mode, double ocs, int S);
	double *proportional_ls_distribution(double mag, int S);
	double *proportional_bs_distribution(double mag, int S);
	double *proportional_ln_distributionBH(double mag, int S);

	double *LNAbundCalc(int R, int num_stdev, int oct_per_stdev, double LS, int *Rich);
	double *normdistodd(int num_stdev, int oct_per_stdev, int starting_oct);
	double *normdistevn(int num_stdev, int oct_per_stdev, int starting_oct);
	int *OctaveRichness(int R, double *NormA, int num_stdev, int oct_per_stdev);
	int *OctaveRichness2(int R, double *NormA, int num_stdev, int oct_per_stdev);

	double *ideal_distribution(double *A, int N);
	
	double *draw_lgn_octaves(int trunc, double mode, int S);
	
#else

    extern double *proportional_gs_distribution(double M, int S); 
    extern double *proportional_zm_distribution(double M, int S); 
	extern double *proportional_lp_distribution(double C, double X, int S);
	extern double *proportional_ln_distribution(int num_stdev, int oct_per_stdev, int starting_octave, double mag_per_oct, int S);
	extern double *proportional_ln_distribution2(int num_stdev, int oct_per_stdev, int modal_octave, double mag_per_oct, int S, int *idistribution);
	extern double *proportional_fln_distribution(double M, double dM, int S, double mode);
	extern double *proportional_lgn_distribution (double mag, int trunc, double mode, int S);
	extern double *proportional_lgn_distributionN(double mag, double mode, double ocs, int S);
	extern double *proportional_lgn_distributionA(double mag, double mode, double ocs, int S);
	extern double *proportional_ls_distribution(double mag, int S);
	extern double *proportional_bs_distribution(double mag, int S);
	extern double *proportional_ln_distributionBH(double mag, int S);

	extern double *LNAbundCalc(int R, int num_stdev, int oct_per_stdev, double LS, int *Rich);
	extern double *normdistodd(int num_stdev, int oct_per_stdev, int starting_oct);
	extern double *normdistevn(int num_stdev, int oct_per_stdev, int starting_oct);
	extern int *OctaveRichness(int R, double *NormA, int num_stdev, int oct_per_stdev);
	extern int *OctaveRichness2(int R, double *NormA, int num_stdev, int oct_per_stdev);

	extern double *ideal_distribution(double *A, int N);

	extern double *draw_lgn_octaves(int trunc, double mode, int S);

#endif