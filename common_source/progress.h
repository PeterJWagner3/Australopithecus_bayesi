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

/***********************************************************************
v0.0	2001.12.20 BY M. KOSNIK
			Works with:
            Apple Computer, Inc. version gcc-926, based on gcc version 2.95.2 19991024 (release)

		2003.01.29 - M. KOSNIK
			Works with:
			Apple Computer, Inc. GCC version 1175, based on gcc version 3.1 20020420 (prerelease)
***********************************************************************/

#ifdef mak_progress
    #include <math.h>
    #include <float.h>
    #include <stdlib.h>
    #include <stdio.h>
    
    int print_progress(int n_reps, int rep, int progress_inc);
     
#else
 
    extern int print_progress(int n_reps, int rep, int progress_inc);
    
#endif