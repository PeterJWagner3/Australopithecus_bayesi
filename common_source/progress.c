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

#define mak_progress
#include "progress.h"

/*	Prints progress report to screen
/*	prints current "rep" every "progress_inc"
/*	written 2002 by M. Kosnik
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int print_progress(n_reps, rep, progress_inc) {

	 if (rep==(n_reps-1)) {
	    if (n_reps<1000)	{
			printf("\b\b\b");
			printf("   ");
		} else {
			printf("\b\b\b\b");
			printf("    ");
		}
		fflush(stdout);
	} else if ((rep%progress_inc)==(progress_inc-1))	{
	
	    if (n_reps<1000)	{
		  if (rep>progress_inc)	{
			printf("\b\b\b");
			fflush(stdout);
		  }
		  printf("%3d",rep+1);
		  fflush(stdout);
		  if (rep>(n_reps-progress_inc))	{
			printf("\b\b\b");
			fflush(stdout);
		  }
	    } else {
		  if (rep>progress_inc)	{
			printf("\b\b\b\b");
			fflush(stdout);
		  }
		  printf("%4d",rep+1);
		  fflush(stdout);
		  if (rep>(n_reps-progress_inc))	{
			printf("\b\b\b\b");
			fflush(stdout);
		  }
	    }
	}
return 0;
}