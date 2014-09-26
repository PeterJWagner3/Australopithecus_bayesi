/************************************************************v0.0 - 2003.01.28 BY WRITTEN BY P.J.WAGNER III ************************************************************/#ifdef historicaldiversity	#include <stdlib.h>	#include <stdio.h>	#include <math.h>	#include "memory.h"	long *richnesstally(long **ranges, int otus);	double *rangethru(long *fa, long *la, int otus);	long *standing(long **ranges, int notu, int S, int st);	long **extant(long **ranges, int notu);	double FreqRat(long *fa, long *la, int otus);	double FreqRatMatrix(long **strat, int notu);	int samefas(long *fa, int otus);	double **gapprescalc(double *gaps, double *rthru, int stg);	double *gapextinction(double *gaps, double *stgaps, double *rng, int stg, int otus);	double completeness(double mu, double r);	double expectedpropsampled(double mu, double r);	unsigned long *rarifiedbylists(unsigned long *sample, int N, int s);	/* richness estimators */	int chao2(int *abundance, int ntaxa);	double chao2real(int *abundance, int ntaxa);	int chao2bc(int *abundance, int ntaxa, int finds);	int jack1(int *abundance, int ntaxa);	int jack2(int *abundance, int ntaxa);	int jack3(int *abundance, int ntaxa);	int jack4(int *abundance, int ntaxa);	int jack5(int *abundance, int ntaxa);		double *chaosharedUV(unsigned long *tad1, unsigned long *tad2, unsigned long *spid1, unsigned long *spid2, int S1, int S2, int w, int z);	double chaosharedrichness(unsigned long *tad1, unsigned long *tad2, unsigned long *spid1, unsigned long *spid2, int S1, int S2, int w, int z);		double *cladeshapefromranges(long **ranges, int otus, int basis);	double *cladeshapefromrichness(long *S, int onset, int end);	double centerofgravityfromranges(long **ranges, int otus, int basis, int onset, int end);	double centerofgravityfromrichness(long *S, int onset, int end);	double rawcenterofgravityfromranges(long **ranges, int otus, int basis);	double rawcenterofgravityfromrichness(long *S, int onset, int end);		double JaccardSim(double A, double B, double C);	double SorensenSim(double A, double B, double C);	double SimpsonSim(double A, double B, double C);	double NestedResultantSim(double A, double B, double C);	double LennonSim(double A, double B, double C);	double JaccardDis(double A, double B, double C);	double SorensenDis(double A, double B, double C);	double SimpsonDis(double A, double B, double C);	double NestedResultantDis(double A, double B, double C);	double LennonDis(double A, double B, double C);		int countshared(unsigned long *taxa1, unsigned long *taxa2, int S1, int S2);	unsigned long *sharedabundance(unsigned long *tad1, unsigned long *taxa1, unsigned long *taxa2, int S1, int S2);	double **sharedrarefaction(unsigned long *finds0, unsigned long *finds1, int n0, int n1, int maxtaxa, int runs);	#else	extern long *richnesstally(long **ranges, int otus);	extern double *rangethru(long *fa, long *la, int otus);	extern long *standing(long **ranges, int notu, int S, int st);	extern long **extant(long **ranges, int notu);	extern double FreqRat(long *fa, long *la, int otus);	extern double FreqRatMatrix(long **strat, int notu);	extern int samefas(long *fa, int otus);	extern double **gapprescalc(double *gaps, double *rthru, int stg);	extern double *gapextinction(double *gaps, double *stgaps, double *rng, int stg, int otus);	extern double completeness(double mu, double r);	extern double expectedpropsampled(double mu, double r);	extern unsigned long *rarifiedbylists(unsigned long *sample, int N, int s);		/* richness estimators */	extern int chao2(int *abundance, int ntaxa);	extern double chao2real(int *abundance, int ntaxa);	extern int chao2bc(int *abundance, int ntaxa, int finds);	extern int jack1(int *abundance, int ntaxa);	extern int jack2(int *abundance, int ntaxa);	extern int jack3(int *abundance, int ntaxa);	extern int jack4(int *abundance, int ntaxa);	extern int jack5(int *abundance, int ntaxa);	extern double *chaosharedUV(unsigned long *tad1, unsigned long *tad2, unsigned long *spid1, unsigned long *spid2, int S1, int S2, int w, int z);	extern double chaosharedrichness(unsigned long *tad1, unsigned long *tad2, unsigned long *spid1, unsigned long *spid2, int S1, int S2, int w, int z);//	extern double *cladeshape(long **ranges, int otus, int basis);	extern double *cladeshapefromranges(long **ranges, int otus, int basis);	extern double *cladeshapefromrichness(long *S, int onset, int end);	extern double centerofgravityfromranges(long **ranges, int otus, int basis, int onset, int end);	extern double centerofgravityfromrichness(long *S, int onset, int end);	extern double rawcenterofgravityfromranges(long **ranges, int otus, int basis);	extern double rawcenterofgravityfromrichness(long *S, int onset, int end);	extern double JaccardSim(double A, double B, double C);	extern double SorensenSim(double A, double B, double C);	extern double SimpsonSim(double A, double B, double C);	extern double NestedResultantSim(double A, double B, double C);	extern double LennonSim(double A, double B, double C);	extern double JaccardDis(double A, double B, double C);	extern double SorensenDis(double A, double B, double C);	extern double SimpsonDis(double A, double B, double C);	extern double NestedResultantDis(double A, double B, double C);	extern double LennonDis(double A, double B, double C);	extern int countshared(unsigned long *taxa1, unsigned long *taxa2, int S1, int S2);	extern unsigned long *sharedabundance(unsigned long *tad1, unsigned long *taxa1, unsigned long *taxa2, int S1, int S2);	extern double **sharedrarefaction(unsigned long *finds0, unsigned long *finds1, int n0, int n1, int maxtaxa, int runs);#endif