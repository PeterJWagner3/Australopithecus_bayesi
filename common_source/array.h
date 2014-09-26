/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
/*	Matthew Kosnik: mkosnik@uchicago.edu
/*
/*	This file is copyright (C) 2000 Matthew Kosnik
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

#ifdef mak_array_functions

	#include <math.h>
	#include <limits.h>
	#include "memory.h"				/* includes my memory functions */
	#include "sort.h"				/* includes my sort functions */

	unsigned long i_arraysum(int *, int);
	unsigned long ul_arraysum(unsigned long *, int);
	unsigned long ul_arraycount(unsigned long *, int);
	long i_arraymin(int *, int);
	long i_arraymax(int *, int);
	int i_arraycell(int *, int, int);
	unsigned long ul_arraymin(unsigned long *, int);
	unsigned long ul_arraymax(unsigned long *, int);
	unsigned long *ul_arrayconfidence(unsigned long *, int, float);
	float ul_arraymean(unsigned long *, int);
	unsigned long *ul_random_to_sequence(unsigned long *, int);
	unsigned long *ul_convert_to_sequential(unsigned long *, int, int);

#else

	extern unsigned long i_arraysum(int *, int);
	extern unsigned long ul_arraysum(unsigned long *, int);
	extern unsigned long ul_arraycount(unsigned long *, int);
	extern long i_arraymin(int *, int);
	extern unsigned long ul_arraymin(unsigned long *, int);
	extern long i_arraymax(int *, int);
	extern unsigned long ul_arraymax(unsigned long *, int);
	extern int i_arraycell(int *, int, int);
	extern unsigned long *ul_arrayconfidence(unsigned long *, int, float);
	extern float ul_arraymean(unsigned long *, int);
	extern unsigned long *ul_random_to_sequence(unsigned long *, int);
	extern unsigned long *ul_convert_to_sequential(unsigned long *, int, int);

#endif