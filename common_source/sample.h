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
***********************************************************************/

#ifdef sample
    #include <math.h>
    #include <float.h>
    #include <stdlib.h>
    #include "memory.h"
    
    int *sample_fixedsize(double *Pro_Dist, int S, int samplesize);
    int *sample_sizexrich(double *Pro_Dist, int S, int minsample, int multsample);
    
#else
 
    extern int *sample_fixedsize(double *Pro_Dist, int S, int samplesize);
    extern int *sample_sizexrich(double *Pro_Dist, int S, int minsample, int multsample);

#endif