/* I have finally decided that NR's memory functions were too difficult.
/*
/* written by M. Kosnik 2000.02.10
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef dynamic_memory
#define dynamic_memory

	#include <iostream>
	
	/* boolean */
	bool **bmatrix(int , int );
	bool ***bcube(int , int , int );
	void free_matrix(bool **, int);
	void free_cube(bool ***, int , int );

	/* integer */
	int ***icube(int , int , int );
	void free_vector(int *);
	void free_matrix(int **, int );
	void free_cube(int ***, int , int );

	/* unsigned shorts */
	unsigned short *usvector(int );
	unsigned short **usmatrix(int , int );
	unsigned short ***uscube(int , int , int );
	void free_vector(unsigned short *);
	void free_matrix(unsigned short **, int );
	void free_cube(unsigned short int ***, int , int );

	/* unsigned longs */
	unsigned long *ulvector(int );
	void free_vector(unsigned long *);
	void free_matrix(unsigned long **, int );

	/* unsigned long longs */
	unsigned long long **ullmatrix(int , int );
	void free_matrix(unsigned long long **, int );

	/* float */
	float ***fcube(int , int , int );
	void free_vector(float *);
	void free_matrix(float **, int );
	void free_cube(float ***, int , int );

	/* double */
	void free_vector(double *);
	void free_matrix(double **, int );

	/* string */
	char *svector(int );
	char **smatrix(int , int );
	void free_vector(char *);
	void free_matrix(char **, int );

#endif

