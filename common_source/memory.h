/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//*/*	Matthew Kosnik: mkosnik@uchicago.edu/*	Peter Wagner: wagnerpj@si.edu/*/*	This file is copyright (C) 2000 Matthew Kosnik/*                             2003 Peter J. Wagner/*/*	This program is free software; you can redistribute it and/or modify it /*	under the terms of version 2 the GNU General Public License as published /*	by the Free Software Foundation./*/*	This program is distributed in the hope that it will be useful, but WITHOUT/*	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or /*	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for /*	more details./*/*	To view a copy of the license go to:/*	http://www.fsf.org/copyleft/gpl.html/*	To receive a copy of the GNU General Public License write the Free Software/* 	Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA./*/*	Copies of this source code are available without cost from:/*	http://geosci.uchicago.edu/paleo/csource//*/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//* I have finally decided that NR's memory functions were too difficult./*/* written by M. Kosnik 2000.02.10/* 2002.09.17 M. Kosnik/*	- major revisons, added wagner's clear* functions./* 2003.01.25 P. Wagner took over this! BWAH HA HAAAA!/*  - further updates by Wagner 2009.0127/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */#ifdef memory	#include <stdlib.h>	#include <stdio.h>	#include <time.h>	#include <math.h>	#include <string.h>	/* boolean	*///	bool **bmatrix(int , int );//	bool ***bcube(int , int , int );//	void free_bmatrix(bool **, int);//	void free_bcube(bool ***, int , int );	/* integer	*/	int *ivector(int length);	void free_ivector(int *v);	int **imatrix(int nrows, int ncolumns);	void free_imatrix(int **m, int nrows, int ncolumns);	int ***icube(int length, int width, int height);	void free_icube(int ***i_cube, int length, int width);	/* double	*/	double *dvector(int length);	void free_dvector(double *v);	double **dmatrix(int nrows, int ncolumns);	void free_dmatrix(double **m, int nrows, int ncolumns);	double ***dcube(int length, int width, int height);	double ****dhypcube(int length, int width, int height, int time);	void free_dcube(double ***d_cube, int length, int width);	/* long	*/	long *lvector(long length);	void free_lvector(long *v);	long **lmatrix(int nrows, int ncolumns);	void free_lmatrix(long **m, int nrows, int ncolumns);	long ***lcube(int length, int width, int height);	void free_lcube(long ***l_cube, int length, int width);	/* unsigned long	*/	unsigned long *ulvector(long length);	void free_ulvector(unsigned long *v);	unsigned long **ulmatrix(int nrows, int ncolumns);	void free_ulmatrix(unsigned long **m, int nrows, int ncolumns);	unsigned long ***ulcube(int length, int width, int height);	void free_ulcube(unsigned long ***ul_cube, int length, int width);	/* float	*/	float *fvector(int length);	void free_fvector(float *v);	float **fmatrix(int nrows, int ncolumns);	void free_fmatrix(float **m, int nrows, int ncolumns);	float ***fcube(int length, int width, int height);	/* character	*/	char *cvector(int length);	void free_cvector(char *v);	char **cmatrix(int nrows, int ncolumns);	void free_cmatrix(char **m, int nrows, int ncolumns);	char ***chcube(int length, int width, int height);	void free_fcube(float ***f_cube, int length, int width);	/* unsigned long longs */	unsigned long long **ullmatrix(int , int );	void free_ullmatrix(unsigned long long **, int );	int *ibigvector(int length);	double *dbigvector(int length);#else	extern int *ivector(int length);	extern void free_ivector(int *v);	extern int **imatrix(int nrows, int ncolumns);	extern void free_imatrix(int **m, int nrows, int ncolumns);	extern double *dvector(int length);	extern void free_dvector(double *v);	extern double **dmatrix(int nrows, int ncolumns);	extern void free_dmatrix(double **m, int nrows, int ncolumns);	extern long *lvector(long length);	extern void free_lvector(long *v);	extern long **lmatrix(int nrows, int ncolumns);	extern void free_lmatrix(long **m, int nrows, int ncolumns);	extern unsigned long *ulvector(long length);	extern void free_ulvector(unsigned long *v);	extern unsigned long **ulmatrix(int nrows, int ncolumns);	extern void free_ulmatrix(unsigned long **m, int nrows, int ncolumns);	extern float *fvector(int length);	extern void free_fvector(float *v);	extern float **fmatrix(int nrows, int ncolumns);	extern void free_fmatrix(float **m, int nrows, int ncolumns);	extern char *cvector(int length);	extern void free_cvector(char *v);	extern char **cmatrix(int nrows, int ncolumns);	extern void free_cmatrix(char **m, int nrows, int ncolumns);	extern int *ibigvector(int length);	extern double *dbigvector(int length);	extern float ***fcube(int length, int width, int height);	extern double ***dcube(int length, int width, int height);	extern double ****dhypcube(int length, int width, int height, int time);	extern int ***icube(int length, int width, int height);	extern long ***lcube(int length, int width, int height);	extern unsigned long ***ulcube(int length, int width, int height);	extern char ***chcube(int length, int width, int height);	extern void free_fcube(float ***f_cube, int length, int width);	extern void free_dcube(double ***d_cube, int length, int width);	extern void free_icube(int ***i_cube, int length, int width);	extern void free_lcube(long ***l_cube, int length, int width);	extern void free_ulcube(unsigned long ***ul_cube, int length, int width);#endif