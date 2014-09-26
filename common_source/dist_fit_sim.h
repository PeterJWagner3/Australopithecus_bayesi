/*	Written by M. Kosnik, posted 2002.xx.xx
/*
/*	Matthew A. Kosnik
/*	Department of Geophysical Sciences
/*	The University of Chicago
/*	5734 South Ellis Avenue
/*	Chicago, Illinois 60637
/*	mkosnik@uchicago.edu
/*
/*	This file is copyright (C) 2002 M. Kosnik
/*
/*	This program is free software; you can redistribute it and/or modify it
/*	under the terms of version 2 of the GNU General Public License as published
/*	by the Free Software Foundation.
/*
/*	This program is distributed in the hope that it will be useful, but WITHOUT
/*	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
/*	FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
/*	more details.
/*
/*	To view a copy of the license go to: http://www.fsf.org/copyleft/gpl.html
/*
/*	To receive a copy of the GNU General Public License write the Free Software
/*	Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
/*
/*	COPIES OF THIS SOURCE CODE ARE AVAILABLE WITHOUT COST FROM:
/*	http://geosci.uchicago.edu/paleo/csource/
/*
/*	v 1.4 - 2001.12
/*		- added (this) main header file
/*		- moved version information to this file.
/*		- moved fitting parameters to tests header file
/*		- changed the way the slope value is incremented to use divisions instead of powers
/*	v 1.6 - 2001.12
/*		- revised code to deal with wagner's revised lognormal generation code.
/*	v 2.0 2002.01
/*		- dropped monte-carlo in favor of Akaike Information Criterion.
/*		- revised header file to include all the simulation parameters
/*	v 2.1 - 2002.02
/*		- changed significance comparsion to only use integer values.
/*	v 2.2 - 2002.02
/*		- changes to fitting functions.
/*	v 2.3 - 2002.02.10
/*		- changes to LN creation functions.
/*	v 2.3.2 - 2002.02.14 - Valentine's Day Special Release
/*		- changes to progress meter.
/*		- average parameters divide by number fit instead of number of reps.
/*		- change to ALL CAPS for constants
/*	v 2.3.3 - 2002.02.17
/*		- minor changes to summary output 
/*			- (no longer outputs NaN when n_fit == 0).
/*			- (no longer outputs extra ln parameter in column list).
/*
/*	v 2.4 - 2002.03.24
/*		- I have no idea what changes wagner made...
/*		- some day I'll train him to work on change logs.
/*		- added Bayes Information Criterion
/*
/*	v 2.4.1 - 2002.04.07
/*		- turned random number generator seed back on
/*		- changed sample size multiples to increment by *= instead of +=
/*		- restored \n in ln output file, prevents lines that are too long
/*		KNOWN ISSUES:
/*		- hang on ln simulations:
/*			LN      2.50     0      4       4       120	100,1x sample size
/*			LN      2.50     4      4       4       180	100,2x sample size
/*
/*	v 2.5.0 - 2002.04.27
/*		- moved to cvs system, did a substantial amount of renaming of files.
/*		- no substantive code changes.
/*
/*	v 2.6 - 2002.07.14
/*		- output clean up (uh, we don't do p-values anymore).
/*		- use aic / bic functions.
/*		- output revised, best fit, good fit, poor fit...
/*
/*	v 2.7 - 2002.07.29
/*		- add log-series support.
/*		- major reorganization.
/*		- add elegible distributions to fit constants.
/*
/*	This code works using gcc 2.96 (MacOSX)
/*
/*	v 2.8 - 2002.09.17
/*		- revised to incompass wagner changes.
/*		- now using: Apple Computer, Inc. GCC version 1161, based on gcc version 3.1 20020420 (prerelease)

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef fit_simulated_distributions

// VERSION INFORMATION
//**************************************************************************
	#define VERSION "2.8"
	#define DATE "2002.09.17"
	#define NAME "evaluate distribution fits"

// DISTRIBUTIONS TO FIT
//**************************************************************************
	#define N_DIST 5			// number of distributions being tested
	#define LN 1			// LN is an option (1 = yes / 0 = no)
	#define GS 1			// GS is an option (1 = yes / 0 = no)
	#define ZM 1			// ZM is an option (1 = yes / 0 = no)
	#define LP 1			// LP is an option (1 = yes / 0 = no)
	// Removed by P. Wagner - DO NOT SET LS = 1
	#define LS 0			// LS is an option (1 = yes / 0 = no)

// RICHNESS PARAMETERS
//**************************************************************************
	#define RICH_START 25		// richness to start at
	#define RICH_STOP 500		// richness to stop at
	#define RICH_INC 5		// amount to increment richness at each step

// SAMPLING PARAMETERS
//**************************************************************************
	#define SAMPLE_MIN 100		// minimum sample size
	#define SAMPLE_MAX 100		// minimum sample size
	#define SAMPLE_INC 100		// minimum sample size
	#define SMULT_MIN 2			// min sample multiple
	#define SMULT_MAX 2			// max sample multiple
	#define SMULT_INC 2			// sample multiple increment value

// BASIC SIMULATION PARAMETERS
//**************************************************************************
	#define MINOCCRS 5		// minimum number of taxa in a list to fit distribution
	#define MINMULT 2			// minimum number of taxa in a list to fit distribution
	#define PROGRESS 0			// how many reps between updating progress report
								// progress = 0 means no progress report.
	#define REPS 100			// number of replicate distributions to sample

// AIC PARAMETERS
//**************************************************************************
	#define AIC_GOOD 10			// value needed to be considered a good fit
	#define AIC_POOR 20			// value needed to be considered a poor fit

// LOG-NORMAL PARAMETERS
//**************************************************************************
//	#define DISTRIB "LN"		// evaluate LN (log-normal)
	#define NSTDDEV 4			// Number of standard deviations
	#define OSTDDEV 4			// Number of octaves per standard deviation
	#define OCT_START 0			// Minimum octave
	#define OCT_STOP 8			// Maximum octave
	#define OCT_INC 2			// Octave increment
	#define OMAG_START 1.25f	// Minimum magnitude of change between octaves
	#define OMAG_STOP 4.25f		// Maximum magnitude of change between octaves
	#define OMAG_INC 0.25f		// Increment magnitude of change between octaves

// ZIPF-MANDELBROT PARAMETERS
//**************************************************************************
	#define DISTRIB "ZM"		// evaluate ZM (Zipf-Mandelbrot)
	#define SLOPE_START 1.00f	// slope to start at
	#define SLOPE_STOP 18.00f	// slope to stop at
	#define SLOPE_INC 0.20f		// amount to increment slope at each step

// GEOMETRIC PARAMETERS
//**************************************************************************
/*	#define DISTRIB "GS"		// evaluate GS (Geometric)
	#define SLOPE_START 1.00f	// slope to start at
	#define SLOPE_STOP 4.00f	// slope to stop at
	#define SLOPE_INC 0.20f		// amount to increment slope at each step
*/
// Removed by P. Wagner - DO NOT USE
// LOG SERIES PARAMETERS
//**************************************************************************
/*	#define DISTRIB "LS"		// evaluate LS (log-series)
	#define SLOPE_START 0.10f	// slope to start at
	#define SLOPE_STOP 0.90f	// slope to stop at
	#define SLOPE_INC 0.10f		// amount to increment slope at each step
*/
// LOG-POWER PARAMETERS
//**************************************************************************
/*	#define DISTRIB "LP"		// evaluate GS (Geometric)
	#define COEF_START 0.2f
	#define COEF_STOP 0.8f
	#define COEF_INC 0.2f
	#define EXP_START 0.2f
	#define EXP_STOP 2.0f
	#define EXP_INC 0.2f
*/
	
// INCLUDES
//**************************************************************************
	#include <stdlib.h>			// just needed.
	#include <stdio.h>			// needed for display
	#include <limits.h>			// needed for RAND_MAX calls
	#include <float.h>			// needed for DBL_MAX calls
	#include <time.h>			// needed for random number seed
	#include <math.h>			// needed for power function, etc.
	#include <string.h>			// needed to make file names

//	THE FOLLOWING FILES ARE NEEDED TO COMPILE THIS PROGRAM:

	#include "distribution_fits.h"		// TEST/FIT DISTRIBUTIONS
	#include "distribution_calc.h"		// CALCULATE DISTRIBUTIONS
	#include "distribution_sample.h"	// RETURNS SAMPLED DISTRIBUTIONS
	#include "memory.h"					// DYNAMIC MEMORY ALLOCATION FUNCTIONS
	#include "sort.h"					// SORTING ROUTINES
	#include "progress.h"				// REPORTS PROGRESS TO SCREEN

#else
#endif
