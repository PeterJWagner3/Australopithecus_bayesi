#include "memory_koz.h"

using namespace std;

bool **bmatrix(int nrows, int ncolumns)
{
int i,j;

/* allocate pointers to rows */
bool **b_matrix = new bool * [nrows];
if (!b_matrix) cout << "allocation error in bmatrix - primary";

/* allocate rows and set pointers to them */
for (i=0 ; i<nrows ; i++) {
	b_matrix[i] = new bool [ncolumns];
	if (!b_matrix[i]) cout << "allocation error in bmatrix - row " << i;
	}

/* initialize to 0 */
for (i=0 ; i<nrows ; i++) {
	for (j=0 ; j<ncolumns ; j++) {
		b_matrix[i][j]=false;
		}
	}

/* return pointer to array of pointers to rows */
return b_matrix;		
}

void free_matrix(bool **b_matrix, int nrows)
{
	int i;
	
	for (i=0 ; i<nrows ; i++) 
		{
		delete [] b_matrix[i];
		b_matrix[i]=0;
		}
	delete [] b_matrix;
	b_matrix=0;
	
	return;
}

bool ***bcube(int length, int width, int height)
{
	int i, j, k;
	
	bool ***b_cube = new bool ** [length];
	if (!b_cube) cout << "allocation error in b_cube - length";

	/* allocate rows and set pointers to them */
	for (i=0 ; i<length ; i++) {
		b_cube[i] = new bool * [width];
		if (!b_cube[i]) cout << "allocation error in b_cube - width " << i;
	}

	for (i=0 ; i<length ; i++) {
		for (j=0 ; j<width ; j++) {
			b_cube[i][j]=new bool [height];
			if (!b_cube[i][j]) cout << "allocation error in b_cube - height " << i<< " x " << j;
		}
	}

	for (i=0 ; i<length ; i++) {
		for (j=0 ; j<width ; j++) {
			for (k=0; k<height ; k++) b_cube[i][j][k]=false;
		}
	}
	
	return b_cube;
}

void free_cube(bool ***b_cube, int length, int width)
{
	int i, j;
	
	for (i=0 ; i<length ; i++) {
		for (j=0 ; j<width ; j++) 
			{
			delete [] b_cube[i][j];
			b_cube[i][j]=0;
			}
	}
	for (i=0 ; i<length ; i++) 
		{
		delete [] b_cube[i];
		b_cube[i]=0;
		}
	delete [] b_cube;
	b_cube=0;
	
	return;
}


void free_vector(int *v)
{
	delete [] v;
	v=0;
	return;		
}

int ***icube(int length, int width, int height)
{
	int i, j, k;
	
	int ***int_cube = new int ** [length];
	if (!int_cube) cout << "allocation error in int_cube - length";

	/* allocate rows and set pointers to them */
	for (i=0 ; i<length ; i++) {
		int_cube[i] = new int * [width];
		if (!int_cube[i]) cout << "allocation error in int_cube - width " << i;
	}

	for (i=0 ; i<length ; i++) {
		for (j=0 ; j<width ; j++) {
			int_cube[i][j]=new int [height];
			if (!int_cube[i]) cout << "allocation error in int_cube - height " << i << " x " << j;
		}
	}

	for (i=0 ; i<length ; i++) {
		for (j=0 ; j<width ; j++) {
			for (k=0; k<height ; k++) int_cube[i][j][k]=0;
		}
	}

	return int_cube;
}

void free_matrix(int **int_matrix, int nrows)
{
	int i;
	
	for (i=0 ; i<nrows ; i++) 
		{
		delete [] int_matrix[i];
		int_matrix[i]=0;
		}
	delete [] int_matrix;
	int_matrix=0;
	
	return;
}

void free_cube(int ***int_cube, int length, int width)
{
	int i, j;
	
	for (i=0 ; i<length ; i++) {
		for (j=0 ; j<width ; j++) 
			{
			delete [] int_cube[i][j];
			int_cube[i][j]=0;
			}
	}
	for (i=0 ; i<length ; i++) 
		{
		delete [] int_cube[i];
		int_cube[i]=0;
		}
	delete [] int_cube;
	int_cube=0;
	
	return;
}

/* Unsigned short allocations
/* written by M. Kosnik 2000.02.10
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* vector
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
unsigned short *usvector(int length)
{
	int i;
	
	unsigned short *v = new unsigned short [length];
	if (!v) cout << "allocation error in ulvector";
	
	/* initialize to 0 */
	for (i=0 ; i<length ; i++) {
		v[i]=0;
	}
	
	/* return pointer to array */
	return v;
}

void free_vector(unsigned short *v)
{
	delete [] v;
	v=0;
	
	return;		
}

/* matrix
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
unsigned short **usmatrix(int nrows, int ncolumns)
{
	int i,j;

	/* allocate pointers to rows */
	unsigned short **m = new unsigned short * [nrows];
	if (!m) cout << "allocation error in ulongmatrix - primary";
	
	/* allocate rows and set pointers to them */
	for (i=0 ; i<nrows ; i++) {
		m[i] = new unsigned short [ncolumns];
		if (!m[i]) cout << "allocation error in ulongmatrix - row " << i;
	}

	/* initialize to 0 */
	for (i=0 ; i<nrows ; i++) {
		for (j=0 ; j<ncolumns ; j++) {
			m[i][j]=0;
		}
	}

	/* return pointer to array of pointers to rows */
	return m;		
}
void free_matrix(unsigned short **m, int nrows)
{
	int i;
	
	for (i=0 ; i<nrows ; i++) delete [] m[i];
	delete [] m;
	m=0;
	
	return;
}

unsigned short ***uscube(int length, int width, int height)
{
	int i, j, k;
	
	unsigned short ***us_cube = new unsigned short ** [length];
	if (!us_cube) cout << "allocation error in us_cube - length";

	/* allocate rows and set pointers to them */
	for (i=0 ; i<length ; i++) {
		us_cube[i] = new unsigned short * [width];
		if (!us_cube[i]) cout << "allocation error in us_cube - width " << i;
	}

	for (i=0 ; i<length ; i++) {
		for (j=0 ; j<width ; j++) {
			us_cube[i][j]=new unsigned short [height];
			if (!us_cube[i]) cout << "allocation error in us_cube - height " << i << " x " <<j;
		}
	}

	for (i=0 ; i<length ; i++) {
		for (j=0 ; j<width ; j++) {
			for (k=0; k<height ; k++) us_cube[i][j][k]=0;
		}
	}
	
	return us_cube;
}

void free_cube(unsigned short ***us_cube, int length, int width)
{
	int i, j;
	
	for (i=0 ; i<length ; i++) {
		for (j=0 ; j<width ; j++) delete [] us_cube[i][j];
	}
	for (i=0 ; i<length ; i++) delete [] us_cube[i];
	delete [] us_cube;
	us_cube=0;
	
	return;
}

/* Unsigned long allocations
/* written by M. Kosnik 2000.02.10
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* vector
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
unsigned long *ulvector(int length)
{
	long i;
	
	unsigned long *v = new unsigned long [length];
	if (!v) cout <<"allocation error in ulvector";
	
	/* initialize to 0 */
	for (i=0 ; i<length ; i++) {
		v[i]=0;
	}
	
	/* return pointer to array */
	return v;
}

void free_vector(unsigned long *v)
{
	delete [] v;
	v=0;
	
	return;		
}

void free_matrix(unsigned long **m, int nrows)
{
	int i;
	
	for (i=0 ; i<nrows ; i++) delete [] m[i];
	delete [] m;
	m=0;
	
	return;
}

/* Unsigned long long allocations
/* written by M. Kosnik 2000.02.10
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* matrix
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
unsigned long long **ullmatrix(int nrows, int ncolumns)
{
	int i,j;

	/* allocate pointers to rows */
	unsigned long long **m = new unsigned long long * [nrows];
	if (!m) cout << "allocation error in ulongmatrix - primary";
	
	/* allocate rows and set pointers to them */
	for (i=0 ; i<nrows ; i++) {
		m[i] = new unsigned long long [ncolumns];
		if (!m[i]) cout << "allocation error in ulongmatrix - row " << i;
	}

	/* initialize to 0 */
	for (i=0 ; i<nrows ; i++) {
		for (j=0 ; j<ncolumns ; j++) {
			m[i][j]=0;
		}
	}

	/* return pointer to array of pointers to rows */
	return m;		
}

void free_matrix(unsigned long long **m, int nrows)
{
	int i;
	
	for (i=0 ; i<nrows ; i++) delete [] m[i];
	delete [] m;
	m=0;
	
	return;
}

void free_vector(float *v)
{
	delete [] v;
	v=0;
	
	return;		
}


void free_matrix(float **m, int nrows)
{
	int i;
	
	for (i=0 ; i<nrows ; i++) delete [] m[i];
	delete [] m;
	m=0;
	
	return;
}

float ***fcube(int length, int width, int height)
{
	int i, j, k;
	
	float ***f_cube = new float ** [length];
	if (!f_cube) cout << "allocation error in f_cube - length";

	/* allocate rows and set pointers to them */
	for (i=0 ; i<length ; i++) {
		f_cube[i] = new float * [width];
		if (!f_cube[i]) cout << "allocation error in f_cube - width " << i;
	}

	for (i=0 ; i<length ; i++) {
		for (j=0 ; j<width ; j++) {
			f_cube[i][j]=new float [height];
			if (!f_cube[i][j]) cout << "allocation error in f_cube - height " << i << " x  " << j;
		}
	}

	for (i=0 ; i<length ; i++) {
		for (j=0 ; j<width ; j++) {
			for (k=0; k<height ; k++) f_cube[i][j][k]=0;
		}
	}
	
	return f_cube;
}

void free_cube(float ***f_cube, int length, int width)
{
	int i, j;
	
	for (i=0 ; i<length ; i++) {
		for (j=0 ; j<width ; j++) 
			{
			delete [] f_cube[i][j];
			f_cube[i][j]=0;
			}
	}
	for (i=0 ; i<length ; i++) 
		{
		delete [] f_cube[i];
		f_cube[i]=0;
		}
	delete [] f_cube;
	f_cube=0;
	
	return;
}

void free_vector(double *v)
{
	delete [] v;
	v=0;
	
	return;		
}


void free_matrix(double **m, int nrows)
{
	int i;
	
	for (i=0 ; i<nrows ; i++) delete [] m[i];
	delete [] m;
	m=0;
	
	return;
}

/* string allocations
/* written by J. Marcot 2002.09.02
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* vector
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
char *svector(int length)
{
	int i;
        
    char *v = new char[length];
	if (!v) cout << "allocation error in fvector";
	
	/* initialize to 0 */
	for (i=0 ; i<length ; i++) {
		v[i]=0;
	}
	
	/* return pointer to array */
	return v;		
}

void free_vector(char *v)
{
	delete [] v;
	v=0;
	
	return;		
}


/* matrix
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
char **smatrix(int nrows, int ncolumns)
{
	int i;
	char **m;

	/* allocate pointers to rows */
	m = new char*[nrows];
	if (!m) cout << "allocation error in dmatrix - primary";
	
	/* allocate rows and set pointers to them */
	for (i=0 ; i<nrows ; i++) {
		m[i] = new char[ncolumns];
		if (!m[i]) cout << "allocation error in dmatrix - row " << i;
	}

	/* initialize to 0 */
//	for (i=0 ; i<nrows ; i++) {
//		for (j=0 ; j<ncolumns ; j++) {
//			m[i][j]="0";
//		}
//	}

	/* return pointer to array of pointers to rows */
	return m;		
}

void free_matrix(char **m, int nrows)
{
	int i;
	
	for (i=0 ; i<nrows ; i++) delete [] m[i];
	delete [] m;
	m=0;
	
	return;
}

