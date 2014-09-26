/*	Peter Wagner: pwagner@fmnh.org
*	Matthew Kosnik: mkosnik@uchicago.edu
*
*	This file is copyright (C) 2001 Peter J. Wagner & Matthew Kosnik
*
*	This program is free software; you can redistribute it and/or modify it 
*	under the terms of version 2 the GNU General Public License as published 
*	by the Free Software Foundation.
*
*	This program is distributed in the hope that it will be useful, but WITHOUT
*	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
*	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
*	more details.
*
*	To view a copy of the license go to:
*	http://www.fsf.org/copyleft/gpl.html
*	To receive a copy of the GNU General Public License write the Free Software
* 	Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*
*	Copies of this source code are available without cost from:
*	http://geosci.uchicago.edu/paleo/csource/
*	
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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
v0.5	2004.03.15 BY P. WAGNER
			redid proportional_lgn_distribution to use 1/Sth of the area under a normal curve, eliminating the need
				to round-off estimates of richness withing octaves.  
			added draw_lgn_octaves to illustrate where 1/Sth of the area is under a normal curve.  
v0.6	2004.04.21 BY P. WAGNER
			added proportional_ls_distribution and proportional_bs_distribution given Hayek & Buzas formulations.
v0.7	2008.09.03 BY P. WAGNER
			surgery performed on log-normal to accommodate change to erf(x) function.
***********************************************************************/

#ifdef distribution_calc
    #include <float.h>
    #include <limits.h>					/* added 2004.12.27 by P. Wagner */
    #include <math.h>					/* added 2004.12.27 by P. Wagner */
    #include <stdlib.h>
    #include <stdio.h>
	#include <string.h>					/* added 2004.12.27 by P. Wagner */
	#include <time.h>					/* added 2004.12.27 by P. Wagner */
    #include "memory.h"
    #include "sort.h"
    #include "probability.h"			/* added 2004.12.27 by P. Wagner */

	#define PI 3.141592654f				/* constant pi					 */
    #define EX 2.718281828f				/* constant e					 */


    double *proportional_geo_distribution(double M, long double min); 
    double *proportional_geos_distribution(double M, int S); 
    double *proportional_zipf_distribution(double M, int S); 
    double *proportional_zm_distribution(double M, double beta, int S); 
	double *proportional_lp_distribution(double C, double X, int S);
	double *proportional_lgn_distribution(double mag, int S);
	double *proportional_tlgn_distribution(double mag, int S, int St);
	double *proportional_zsm_distribution(double theta, double m, int J);
	double *draw_lgn_octaves(int trunc, double mode, int S);
	double *proportional_ls_distribution(double mag, int S);
	double *proportional_bs_distribution(double mag, int S);
	double *proportional_ln_distributionBH(double mag, int S);
	double *proportional_from_expt_inds(double *expected, int max);

	double *normdistodd(int num_stdev, int oct_per_stdev, int starting_oct);
	double *normdistevn(int num_stdev, int oct_per_stdev, int starting_oct);

	double *ideal_distribution(double *A, int N);
	
	double *zerosum(double theta, double m, long int J);
	double volkov(double theta, double gam, int n, long int J, double coef, double y);
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	double qromo(double (*func)(double,double,int,long int,double,double), double theta, double gam, int nn, long int J, double coef, double a, double b, double (*choose)(double(*)(double,double,int,long int,double,double), double, double, int, long int, double, double, double, int));
	double midpnt(double (*func)(double,double,int,long int,double,double), double theta, double gam, int nn, long int J, double coef, double a, double b, int n);
	
	double *stateschangefromcharacterchange(double **PComp, int *nstates, int nchars, int mxstp);
	
	
#else

    extern double *proportional_geo_distribution(double M, long double min); 
    extern double *proportional_geos_distribution(double M, int S); 
    extern double *proportional_zipf_distribution(double M, int S); 
    extern double *proportional_zm_distribution(double M, double beta, int S); 
	extern double *proportional_lp_distribution(double C, double X, int S);
	extern double *proportional_lgn_distribution(double mag, int S);
	extern double *proportional_tlgn_distribution(double mag, int S, int St);
	extern double *proportional_zsm_distribution(double theta, double m, int J);
	extern double *draw_lgn_octaves(int trunc, double mode, int S);
	extern double *proportional_ls_distribution(double mag, int S);
	extern double *proportional_bs_distribution(double mag, int S);
	extern double *proportional_ln_distributionBH(double mag, int S);
	extern double *proportional_from_expt_inds(double *expected, int max);

	extern double *normdistodd(int num_stdev, int oct_per_stdev, int starting_oct);
	extern double *normdistevn(int num_stdev, int oct_per_stdev, int starting_oct);

	extern double *ideal_distribution(double *A, int N);
	
	extern double *zerosum(double theta, double m, long int J);
	extern double volkov(double theta, double gam, int n, long int J, double coef, double y);
	extern void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	extern double qromo(double (*func)(double,double,int,long int,double,double), double theta, double gam, int nn, long int J, double coef, double a, double b, double (*choose)(double(*)(double,double,int,long int,double,double), double, double, int, long int, double, double, double, int));
	extern double midpnt(double (*func)(double,double,int,long int,double,double), double theta, double gam, int nn, long int J, double coef, double a, double b, int n);

	extern double *stateschangefromcharacterchange(double **PComp, int *nstates, int nchars, int mxstp);

#endif