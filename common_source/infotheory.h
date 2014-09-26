#define SUPINC 0.05f				/* minimum change in support 	*/#define FITINC 1					/* increment to fit by			*/#include "memory.h"#include "matrixchange.h"#ifdef info_theory	#include <math.h>	#include <float.h>	#include <stdlib.h>	#include <stdio.h>	double calc_aic_c(int parameters, double fit, int measurements); 	double calc_aic(int parameters, double fit); 	double calc_bic(int parameters, double fit, int measurements);	double **testpartitions(double **loglike, int cases, int hyp, double *hypos, double dp, int n);#else	extern double calc_aic_c(int parameters, double fit, int measurements); 	extern double calc_aic(int parameters, double fit); 	extern double calc_bic(int parameters, double fit, int measurements); 	extern double **testpartitions(double **loglike, int cases, int hyp, double *hypos, double dp, int n);#endif