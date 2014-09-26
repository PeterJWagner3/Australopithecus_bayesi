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

/*
v b3.0 - 2002.02.14 - M. Kosnik
	updated to used revised fitting functions.
	updated constants to ALL CAPS.
		
v b3.1 - 2002.02.22 - M. Kosnik
	removed monte carlo pvalue stuff (it was soo 2.x, like totally yesterday's idea).
	added AIC calculation, revised output format.
	cleaned up code, removing unused variables.
	created header file.
	rewrote introductory splash.
	
v b3.2 - 2002.02.27 - M. Kosnik
	replaced AIC with AICc (small number of oberservations correction)
	added summary at end.

v b3.3 - 2002.04.27 - M. Kosnik
	moved code to CVS tree, had to do some renaming... no significant changes.

v b3.4 - 2002.07.13 - M. Kosnik
	use average number of specimens per taxon to decide to fit.
	minor cosmetic changes (they are taxa not species).
	move aic calculation to function.
	increase (maximized?) use of header defined constants.
	add comments in code.
	fix output file header rows (previously richness columns not labeled).
	minor changes to distribution output file (only output real if fitting distribution)
	removed unneeded header files.

v b3.5 - 2002.07.14 - M. Kosnik
	bug fix - now correctly displays reasoning behind best AICc determination
			- code fix, bug not causing errors.

v b3.6 - 2002.07.31 - M. Kosnik
	added support for log-series and log-power.

v b3.7 - 2002.09.17 - M. Kosnik
	changes to keep up with wagner changes.
	added string.h, time.h
	
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifdef fit_emperical_distributions

	#include <stdlib.h>
	#include <stdio.h>
	#include <float.h>
	#include <string.h>
	#include <time.h>
	
	//THE FOLLOWING FILES (AND CORRESPONDING .c FILE) ARE NEEDED & SHOULD BE DISTRIBUTED WITH THIS FILE.
	
	#include "file.h"					// functions to read the data from a file
	#include "memory.h"					// functions deal with memory allocation
	#include "matrix.h"					// functions deal with get summary information from data matrix
	#include "sort.h"					// sort functions
	#include "distribution_evenness.h"
	#include "distribution_calc.h"
	#include "distribution_fits.h"
	
	// VERSION INFORMATION
	#define VERSION "b3.7"
	#define DATE "2002.09.17"
	
	// INPUT FILE FORMAT (numbers refer to columns in data file (starting with 0))
	#define NAME 255			// maximum size of file name
	#define C_LIST 0			// list number
	#define C_SPEC 1			// species number
	#define C_ABUN 2			// species abundance
	#define C_TOTAL 3			// total number of columns
	
	// FITTING LIMITS
	#define MINMULT 2			// minimum average number of specimens per taxon to fit distribution
	#define MINOCCRS 5			// minimum number of taxa in a list to fit distribution
								// this must be at least 1 more than the max number of parameters.

// DISTRIBUTIONS TO FIT
//**************************************************************************
	#define N_DIST 5			// number of distributions being tested
	#define LN 1				// LN is an option (1 = yes / 0 = no)
	#define GS 1				// GS is an option (1 = yes / 0 = no)
	#define ZM 1				// ZM is an option (1 = yes / 0 = no)
	#define LS 0				// LS is an option (1 = yes / 0 = no)
	#define LP 1				// LP is an option (1 = yes / 0 = no)

	// Delta AIC MEANINGS
	#define AIC_GOOD 10			// value needed to be considered a good fit
	#define AIC_PREC 5			// prevsion of aic comparison
		
	// FIXED LOG-NORMAL PARAMETERS
	#define NSTD 5				// number of standard deviations of log-normal to use
	#define OSTD 5				// number of octaves per standard devation to use
	
#else

#endif