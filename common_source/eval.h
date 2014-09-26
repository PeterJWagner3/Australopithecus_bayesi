/*	Written by M. Kosnik, posted 2002.xx.xx
/*
/*	Matthew A. Kosnik
/*	Department of Geophysical Sciences
/*	The University of Chicago
/*	5734 South Ellis Avenue
/*	Chicago, Illinois 60637
/*	m-kosnik@uchicago.edu
/*
/*	This code works using gcc 2.96 (MacOSX)
/*
/*	This file is copyright (C) 2002 M. Kosnik
/*
/*	This program is free software; you can redistribute it and/or modify it
/*	under the terms of version 2 of the GNU General Public License as published
/*	by the Free Software Foundation.
/*
/*	This program is distributed in the hope that it will be useful, but WITHOUT
/*	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
/*	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
/*	more details.
/*
/*	To view a copy of the license go to: http://www.fsf.org/copyleft/gpl.html
/*
/*	To receive a copy of the GNU General Public License write the Free Software
/*	Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
/*
/*	COPIES OF THIS SOURCE CODE ARE AVAILABLE WITHOUT COST FROM:
/*
/*  	http://geosci.uchicago.edu/paleo/csource/
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
/*	v 2.1  - 2002.02
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
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef eval_main
// VERSION INFORMATION
//**************************************************************************
    #define VERSION "2.3.3"
    #define DATE "2002.02.17"
    #define NAME "evaluate distributions"

// WHICH DISTRIBUTION TO SIMULATE
//**************************************************************************
    #define DISTRIB "LN"		// evaluate LN (log-normal)
//    #define DISTRIB "GS"		// evaluate GS (Geometric)
//    #define DISTRIB "ZM"		// evaluate ZM (Zipf-Mandelbrot)

// RICHNESS PARAMETERS
//**************************************************************************
    #define RICH_START 20		// richness to start at
    #define RICH_STOP 20		// richness to stop at
    #define RICH_INC 50			// amount to increment richness at each step

// SAMPLING PARAMETERS
//**************************************************************************
    #define SAMPLE_MIN 100		// minimum sample size
    #define SAMPLE_MAX 100000		// maximum sample size
    #define SAMPLE_INC 100000		// minimum sample size
    #define SMULT_MIN 1			// min sample multiple
    #define SMULT_MAX 1			// max sample multiple
    #define SMULT_INC 2			// sample multiple increment value

// BASIC SIMULATION PARAMETERS
//**************************************************************************
    #define MINOCCRS 10				// minimum number of taxa in a list to fit distribution
	// progress_inc = 0 means no progress.
    #define PROGRESS 0				// how many reps between updating progress report
    #define REPS 1				// number of replicate distributions to sample

// LOG-NORMAL PARAMETERS
//**************************************************************************
    #define NSTDDEV 4			// Number of standard deviations
    #define OSTDDEV 4			// Number of octaves per standard deviation
    #define OCT_START 0			// Starting octave
    #define OCT_STOP 8
    #define OCT_INC 2
    #define OMAG_START 1.5f		// Magnitude of change between octaves
    #define OMAG_STOP 3.5f
    #define OMAG_INC 0.5f

// GEOMETRIC AND ZIPF-MANDELBROT PARAMETERS
//**************************************************************************
    #define SLOPE_START 1.00f		// slope to start at
    #define SLOPE_STOP 5.00f		// slope to stop at
    #define SLOPE_INC 0.10f		// amount to increment slope at each step

// LOG-POWER PARAMETERS
//**************************************************************************
    #define LP 0				// log power is an option (1 = yes / 0 = no)
    #define starting_coef 1.0f
    #define stopping_coef 3.5f
    #define coef_increment 0.5f
    #define starting_exp 1.0f
    #define stopping_exp 3.5f
    #define exp_increment 0.5f
    
// INCLUDES
//**************************************************************************
    #include <stdlib.h>				// just needed.
    #include <stdio.h>				// needed for display
    #include <limits.h>				// needed for RAND_MAX calls
    #include <float.h>				// needed for DOUBLE_MAX calls
    #include <time.h>				// needed for random number seed
    #include <string.h>				// needed to make file names
    #include "distribution_fits.h"		// TEST/FIT DISTRIBUTIONS
    #include "distribution_calc.h"		// CALCULATE DISTRIBUTIONS
    #include "sample.h"				// RETURNS SAMPLED DISTRIBUTIONS
    #include "memory.h"				// DYNAMIC MEMORY ALLOCATION FUNCTIONS
    #include "sort.h"				// SORTING ROUTINES
    #include "progress.h"			// REPORTS PROGRESS TO SCREEN
    #include "calculus.h"

#else
#endif
 					  					  