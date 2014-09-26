/* Routines to summarize and matrices/arrays already in files not yet in computer memory
/* 		These are useful for determining how much memory should be allocated for arrays and matrices
/*		Use these routines first, then allocate memory based on results, then re-open files and read
/*			the files "formally."
/*						Peter Wagner	05/02
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define filereading
#include "filereading.h"
#include "memory.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

//* getmatrixinfo - Finds dimensions of a matrix without tabs stored as a file
/* Requires:
/*	filename - the file holding the matrix
/* Returns:
/*	info[0] - number of columns (e.g., characters in a cladistic matrix)
/*	info[1] - number of rows (e.g., taxa in a cladistic matrix)
/* pre 9/10/2003: MatrixInfo
****************************************************************************/
int *getmatrixinfo(char *filename)
{
int		col;
int		*info;
char	i, li;
FILE	*Data;

info=ivector(2);
Data = fopen(filename,"r");

li=i='s';
col=0;
/*while ((i!='\n' && i!='\r') && i!='\m')	{	*/
while (i!='\n' && i!='\r')	{
	i = fgetc(Data);
	if (i=='\n')	break;
	else if (i=='\t' || i==' ')	col=col;
	else if (i=='(')	{
		++col;
		while (i!=')')		i = fgetc(Data);
		}
	else
		++col;
	}
/*++col;*/		/* sometimes this is necessary - it is unpredictable */

fclose(Data);

/*Data = fopen(filename,"r");
row=0;
while (!feof(Data))	{
	i='s';
	for (c=0; i!='ˇ' && (i!='\n' && i!='\r'); c=c)
		i = fgetc(Data);
	++row;
	}
fclose(Data);	*/
	
info[0]=col;
info[1]=getfilerows(filename);
return info;
}

//* gettabdeldim - Finds dimensions of a tab-delimited matrix stored as a file
/* Requires:
/*	filename - the file holding the matrix
/* Returns:
/*	info[0] - number of columns (e.g., characters in a cladistic matrix)
/*	info[1] - number of rows (e.g., taxa in a cladistic matrix)
/* pre 9/10/2003: MatrixInfo
****************************************************************************/
int *gettabdeldim(char *filename)
{
int		col, row;
int		*info;
char	i;
FILE	*Data;

info=ivector(2);
Data = fopen(filename,"r");

i = fgetc(Data);
col=0;
for (i=i; i!='\n' && i!='\r'; i=i)	{
	i = fgetc(Data);
	if (i=='\t')	++col;
	}
++col;
fclose(Data);

Data = fopen(filename,"r");
row=0;
i = fgetc(Data);
while (!feof(Data))	{
	i = fgetc(Data);
	if (i=='\n' || i=='\r')	{
		++row;
		}
	}
++row;
fclose(Data);
	
info[0]=col;
info[1]=row;
return info;
}


//* getcladisticmatrixinfo - Finds dimensions of a matrix stored as a file
/* Requires:
/*	filename - the file holding the matrix
/* Returns:
/*	info[0] - number of columns (e.g., characters in a cladistic matrix)
/*	info[1] - number of rows (e.g., taxa in a cladistic matrix)
/* pre 9/10/2003: cladisticmatrixInfo
****************************************************************************/
int *getcldmatrixinfo(char *filename, char inapcode, char unknown)
{
int		col;
int		*info;
char	i, li;
FILE	*Data;

info=ivector(2);
Data = fopen(filename,"r");

li=i='s';
col=0;
/*while ((i!='\n' && i!='\r') && i!='\m')	{	*/
while (i!='\n' && i!='\r')	{
	i = fgetc(Data);
	if (i=='\t' || i==' ')
		++col;
/*	else if ((i=='\n' || i=='\r') || i=='\m')	*/
	else if (i=='\n' || i=='\r')
		++col;
	}
fclose(Data);
if (col<=1)	{
	Data = fopen(filename,"r");

	li=i='s';
	col=0;
	/*while ((i!='\n' && i!='\r') && i!='\m')	{	*/
	while (i!='\n' && i!='\r')	{
		i = fgetc(Data);
		if (i=='\n')	break;
		else if (i=='&')
			i = fgetc(Data);
		else if (i=='('/*)*/)	{
			while (i!=/*(*/')')	i=fgetc(Data);
			++col;
			}
		/* do not count if 2nd, 3rd, 4th, etc. digit is given */
/*		else if ((i>='0' && i<='9') && (li>='0' && li<='9'))	*/
		else if ((i>='0' && i<='9') || ((i>='a' && li<='z') || (i>='A' && li<='Z')))
			++col;
		else if (i==inapcode || i==unknown)
			++col;
		li=i;
		}
	fclose(Data);
	}

/*Data = fopen(filename,"r");
row=0;
i = ' ';
while (!feof(Data))	
{
	i=fgetc(Data);
	if ( i=='\n' || i=='\r' )
		++row;
}
fclose(Data);	*/
	
info[0]=col;
/*info[1]=row;*/
info[1]=getfilerows(filename);
return info;
}


//* getimatrixinfo - Finds dimensions of a matrix stored as a file
/*		NOTE: if there are non-integers (e.g., question marks or gaps)
/*			then use MatrixInfo
/* Requires:
/*	filename - the file holding the matrix
/* Returns:
/*	info[0] - number of columns (e.g., characters in a cladistic matrix)
/*	info[1] - number of rows (e.g., taxa in a cladistic matrix)
/* pre 9/10/2003: imatrixInfo
****************************************************************************/
int *getimatrixinfo(char *filename, int C)
{
int		row, col;
int		*info;
int		i;
FILE	*Data;

info=ivector(2);

Data = fopen(filename,"r");
row=0;
while (!feof(Data))	{
	for (col=0; col<C; col=col)	{
		fscanf(Data,"%i",&i);
		}
	++row;
	}
fclose(Data);
	
info[0]=col;
info[1]=row;
return info;
}

//* Readimatrix - Reads and store a matrix OF INTEGER VALUES
/*		NOTE: if there are non-integers (e.g., question marks or gaps)
/*			then use Readcmatrix
/* Requires:
/*	filename - the file holding the matrix
/*	R - number of rows (e.g., number of taxa)
/*	C - number of columens (e.g., number of characters)
/* Returns:
/*	matrix - the RxC matrix
****************************************************************************/
int	**readimatrix(char *filename, int R, int C)
{
int		row, col;
int		**matrix;
FILE	*Data;

matrix=imatrix(R,C);
Data = fopen(filename,"r");

for (row=0; row<R; ++row)	{
	for (col=0; col<C; ++col)	{
		fscanf(Data,"%i",&matrix[row][col]);
		}
	}

fclose(Data);
return matrix;
}

//* Readimatrix - Reads and store a matrix OF INTEGER VALUES
/*		NOTE: if there are non-integers (e.g., question marks or gaps)
/*			then use Readcmatrix
/* Requires:
/*	filename - the file holding the matrix
/*	R - number of rows (e.g., number of taxa)
/*	C - number of columens (e.g., number of characters)
/* Returns:
/*	matrix - the RxC matrix
****************************************************************************/
double	**readdmatrix(char *filename, int R, int C)
{
int		row, col;
double	**matrix;
FILE	*Data;

matrix=dmatrix(R,C);
Data = fopen(filename,"r");

for (row=0; row<R; ++row)	{
	for (col=0; col<C; ++col)	{
		fscanf(Data,"%lf",&matrix[row][col]);
		}
	}

fclose(Data);
return matrix;
}

//* Readlmatrix - Reads and store a matrix OF LONG VALUES
/*		NOTE: if there are non-integers (e.g., question marks or gaps)
/*			then use Readcmatrix
/* Requires:
/*	filename - the file holding the matrix
/*	R - number of rows (e.g., number of taxa)
/*	C - number of columens (e.g., number of characters)
/* Returns:
/*	matrix - the RxC matrix
****************************************************************************/
long **readlmatrix(char *filename, int R, int C)
{
int		row, col;
long	**matrix;
FILE	*Data;

matrix=lmatrix(R,C);
Data = fopen(filename,"r");

for (row=0; row<R; ++row)	{
	for (col=0; col<C; ++col)	{
		fscanf(Data,"%i",&matrix[row][col]);
		}
	}

fclose(Data);
return matrix;
}

//* Readclmatrix - Reads and store a matrix that includes question marks and inapplicables
/*		NOTE: if there are non-integers (e.g., question marks or gaps)
/*			then use Readcmatrix
/* Requires:
/*	filename - the file holding the matrix
/*	R - number of rows (e.g., number of taxa)
/*	C - number of columens (e.g., number of characters)
/*	inapcode - the code used for inapplicable characters (e.g., the feather color of a crocodile)
/*	UNKNOWN - the number assigned to question marks (e.g., unknowns such as feather color of a fossil bird)
/*	INAP - the number assigned to inapplicable characters
/* Returns:
/*	matrix - the RxC matrix
****************************************************************************/
long **readclmatrix(char *filename, int R, int C, char inapcode, int UNKNOWN, int INAP)
{
int		row, col;
int		*info;
char	i;
long	**matrix;
FILE	*Data;

info=getcldmatrixinfo(filename, inapcode, '?');

matrix=lmatrix(R,C);
Data = fopen(filename,"r");
for (row=0; row<info[1]; ++row)	{
	i = fgetc(Data);
	while (i=='\n')	i = fgetc(Data);
	for (col=0; col<info[0]; col=col)	{
		if (i=='('/*)*/)	{
			while (i!=/*(*/')')	i=fgetc(Data);
			matrix[row][col]=UNKNOWN;
			++col;
			i=fgetc(Data);
			}	/* read through polymorphics labeled by parens	*/
		else if (i=='?')	{
			matrix[row][col]=UNKNOWN;
			++col;
			i=fgetc(Data);
			}
		else if (i==inapcode)	{
			matrix[row][col]=INAP;
			++col;
			i=fgetc(Data);
			}
		else if (i!='\t' && i!=' ')	{
			if (i>='0' && i<= '9')
				matrix[row][col]=i-48;
			else
				matrix[row][col]=10+(i-'A');
			++col;
			i=fgetc(Data);
			if (i=='&' || i=='/')	{
//				matrix[row][col-1]=UNKNOWN;
				i=fgetc(Data);
				while ((i!='\t' && i!=' ') && i!='\n')	{
					if (i!='&' && i!='/')	{
						if (matrix[row][col-1]!=0)	{
							matrix[row][col-1]=-10*abs(matrix[row][col-1]);
							matrix[row][col-1]-=(i-'0');		/* set to -12 if 1&2 */
							}
						else	{
							matrix[row][col-1]=-10*(i-'0');		/* set to -10 if 0&1 */
							}
						}
					i = fgetc(Data);
					}	/* tab, space or end of line will mark end of character	*/
				}	/* read through polymorphics	*/
			}
		else	i=fgetc(Data);
		}
//	i = fgetc(Data);	/* read end of line marker */
	}
fclose(Data);

free_ivector(info);
return matrix;
}

//* Readclmatrix - Reads and store a matrix that includes question marks and inapplicables
/*		NOTE: if there are non-integers (e.g., question marks or gaps)
/*			then use Readcmatrix
/* Requires:
/*	filename - the file holding the matrix
/*	R - number of rows (e.g., number of taxa)
/*	C - number of columens (e.g., number of characters)
/*	inapcode - the code used for inapplicable characters (e.g., the feather color of a crocodile)
/*	UNKNOWN - the number assigned to question marks (e.g., unknowns such as feather color of a fossil bird)
/*	INAP - the number assigned to inapplicable characters
/* Returns:
/*	matrix - the RxC matrix
****************************************************************************/
long **readclmatrixpoly(char *filename, int R, int C, char inapcode, int UNKNOWN, int INAP)
{
int		a, row, col;
int		*info;
char	i;
long	**matrix;
FILE	*Data;

info=getcldmatrixinfo(filename, inapcode, '?');

matrix=lmatrix(R,C);
Data = fopen(filename,"r");
for (row=0; row<info[1]; ++row)	{
	for (col=0; col<info[0]; col=col)	{
		i = fgetc(Data);
		if (i=='&')	{
			i = fgetc(Data);
			matrix[row][col]=matrix[row][col]-(a*(i-'0'));
			}
		else if (i=='('/*)*/)	{
			i=fgetc(Data);
			matrix[row][col]=-1*(i-'0');
			i=fgetc(Data);
			a=10;
			while (i!=/*(*/')')	{
				if (i!='&' && i!=',')
					matrix[row][col]=matrix[row][col]-(a*(i-'0'));
				i=fgetc(Data);
				a*=10;
				}
			++col;
			}
		else if (i=='?')	{
			matrix[row][col]=UNKNOWN;
			++col;
			}
		else if (i==inapcode)	{
			matrix[row][col]=INAP;
			++col;
			}
		else if (i!='\t' && i!=' ')	{
			if (i>='0' && i<= '9')
				matrix[row][col]=i-48;
			else
				matrix[row][col]=10+(i-'A');
			++col;
			}
		}
	i = fgetc(Data);	/* read end of line marker */
	}
fclose(Data);

free_ivector(info);
return matrix;
}

//* readivectorcol - Reads and stores a array of integers
/* Requires:
/*	filename - the file holding the array
/*	N - number of cells (i.e, array length)
/* Returns:
/*	A - the array
****************************************************************************/
int *readivectorcol(char *filename, int N)
{
int		i, *A;
FILE	*Data;

A=ivector(N);
Data = fopen(filename,"r");

/*while (!feof(Data))	{*/
for (i=0; i<N; ++i)	{
	fscanf(Data,"%i",&A[i]);
	}
fclose(Data);

return A;
}

//* readlvectorcol - Reads and stores a array of long
/* Requires:
/*	filename - the file holding the array
/*	N - number of cells (i.e, array length)
/* Returns:
/*	A - the array
****************************************************************************/
long *readlvectorcol(char *filename, int N)
{
int		i;
long	*A;
FILE	*Data;

A=lvector(N);
Data = fopen(filename,"r");

/*while (!feof(Data))	{*/
for (i=0; i<N; ++i)	{
	fscanf(Data,"%i",&A[i]);
	}
fclose(Data);

return A;
}

//* readlvectorcol - Reads and stores a array of long
/* Requires:
/*	filename - the file holding the array
/*	N - number of cells (i.e, array length)
/* Returns:
/*	A - the array
****************************************************************************/
unsigned long *readulvectorcol(char *filename, int N)
{
int		i;
unsigned long	*A;

FILE	*Data;

A=ulvector(N);
Data = fopen(filename,"r");

/*while (!feof(Data))	{*/
for (i=0; i<N; ++i)	{
	fscanf(Data,"%i",&A[i]);
	}
fclose(Data);

return A;
}


//* readdvectorcol - Reads and stores a array of double
/* Requires:
/*	filename - the file holding the array
/*	N - number of cells (i.e, array length)
/* Returns:
/*	A - the array
****************************************************************************/
double *readdvectorcol(char *filename, int N)
{
int		i;
double	*A;
FILE	*Data;

A=dvector(N);
Data = fopen(filename,"r");

/*while (!feof(Data))	{*/
for (i=0; i<N; ++i)	{
	fscanf(Data,"%lf",&A[i]);
	}
fclose(Data);

return A;
}

//* getfilecols - Finds the number of columns (length) of an matrix/array stored as a file
/* Requires:
/*	filename - the file holding the array
/* Returns:
/*	N - number of cells (i.e, array length)
/* Pre 9/10/2003: "FileLength"
****************************************************************************/
int getfilecols(char *filename)
{
int	i, lng;
FILE *Data;

Data = fopen(filename,"r");
lng=i=0;
while (!feof(Data))	{
	fscanf(Data,"%i",&i);
	++lng;
	}
fclose(Data);
return lng;
}

//* getfilecolsfromtext - Finds the number of columns (length) of an matrix/array stored as a file
/* Requires:
/*	filename - the file holding the array
/* Returns:
/*	N - number of cells (i.e, array length)
/* Pre 9/10/2003: "FileLength"
****************************************************************************/
int getfilecolsfromtext(char *filename)
{
int	lng;
char i;
FILE *Data;

Data = fopen(filename,"r");
lng=1;
while (!feof(Data))	{
	fscanf(Data,"%c",&i);
	if (i=='\t')		++lng;
	else if (i=='\n')	break;
	}
fclose(Data);
return lng;
}

int getlongestnamefromtabdelimitedtext(char *filename)
{
int j=0, lng=0;
char	i;
FILE *Data;

Data = fopen(filename,"r");

while (!feof(Data))	{
	fscanf(Data,"%c",&i);
	if (i!='\t' && i!='\n')	++j;
	else	{
		if (j>lng)	lng=j;
		j=0;
		}
	}
fclose(Data);
return lng;

}

//* getfilerows - Finds the number of rows in a matrix file
/* Requires:
/*	filename - the file holding the array
/* Returns:
/*	N - number of cells (i.e, array length)
/* Pre 9/10/2003: "FileRows"
****************************************************************************/
int getfilerows(char *filename)
{
int		row;
char	i;
FILE	*Data;

system("pwd");
	
/*char buffer5[512];	*/
/*printf("current working directory = %s\n", getcwd(NULL));		FOR XCODE ONLY!	*/

Data = fopen(filename,"r");

row=1;		/* it does not read end-of-file marker, or so it seems	*/
i = ' ';
while (!feof(Data))	{
	i=fgetc(Data);
	if ( i=='\n' || i=='\r' )
		++row;
	}
fclose(Data);

return row;
}
//* getfilemin - Finds the minimum value of an array stored as a file
/* Requires:
/*	filename - the file holding the array
/* Returns:
/*	min - smallest value in file
/* Pre 9/10/2003: "FileMin"
****************************************************************************/
int getfilemin(char *filename)
{
int	i, min;
FILE *Data;

Data = fopen(filename,"r");
min=RAND_MAX;
while (!feof(Data))	{
	fscanf(Data,"%i",&i);
	if (i<min)	min=i;
	}
fclose(Data);
return min;
}

//* getfilemax - Finds the maximum value of an array stored as a file
/* Requires:
/*	filename - the file holding the array
/* Returns:
/*	max - largest number in file
/* Pre 9/10/2003: "FileMax"
****************************************************************************/
int getfilemax(char *filename)
{
int	i, max;
FILE *Data;

Data = fopen(filename,"r");
max=i=0;
while (!feof(Data))	{
	fscanf(Data,"%i",&i);
	if (i>max)	max=i;
	}
fclose(Data);
return max;
}

/* getlongestnameinfile - Reads a list of names and finds the longest
/* Requires:
/*		filename - the file holding the array
/*		names - the number of names
/* Returns:
/*		max - length of the longest name
/* Was LongestNameinFile
****************************************************************************/
int getlongestnameinfile(char *filename, int names)
{
int	n, max, m;
char c;
FILE *Data;

Data = fopen(filename,"r");

max=0;
c=' ';
for (n=0; n<names && !feof(Data); ++n)	{
	m=1;
	while (c!='\n' && !feof(Data))	{
		c=fgetc(Data);
		++m;
		}
	if (m>max)	max=m;
	c=' ';
	}
++max;
fclose(Data);
return max;
}

/* getnames - Reads a list of names and stores it as a matrix
/* Requires:
/*		filename - the file holding the array
/*		names - the number of names
/* Returns:
/*		NameArray - an "array" of names (OK, really a matrix, but....)
/* was ReadnameList
****************************************************************************/
char **getnames(char *filename, int names, int max)
{
int	a,n;
char **NameArray,c;
FILE *Data;

Data = fopen(filename,"r");
NameArray=cmatrix(names,max+1);
for (n=0; n<names; ++n)	{
	for (a=0; a<max && (c!='\n' && !feof(Data)); ++a)	{
		c=fgetc(Data);
		if (c!='\n')	NameArray[n][a]=c; 
		}
	NameArray[n][a-1]='\0';
	c=' ';
/*	fgets(NameArray[n],max,Data);
	fscanf(Data,"%[^\n]",NameArray[n]);*/
	}
fclose(Data);
return NameArray;
}

char ***getinstructions(char *filename)
{
int a=0, b=0, c=0;
int x, y, z;
char d;
FILE *Data;
char ***instructions;

x=getfilerows(filename);
y=getfilecolsfromtext(filename);
z=getlongestnamefromtabdelimitedtext(filename);
z+=10;		/* do this to make sure that there is a final marker at the end of each name	*/

instructions=chcube(x,y,z);

Data = fopen(filename,"r");
while (!feof(Data))	{
	d=fgetc(Data);
	if (d=='\t')	{
		++b;
		c=0;
		}
	else if (d=='\n')	{
		++a;
		b=0;
		c=0;
		}
	else	{
		instructions[a][b][c]=d;
		++c;
		}
	if(a==x)	break;
	}
fclose(Data);
return instructions;
}
/**/


/* 2013-07-11
//* readlvectorcol - Reads and stores a array of long
/* Requires:
/*	filename - the file holding the array
/*	N - number of cells (i.e, array length)
/* Returns:
/*	A - the array
****************************************************************************/
unsigned long *getuldataforvectorA(char *filename, int A, int N)
{
int		i;
unsigned long	*B;
FILE	*Data;

B=ulvector(N);
Data = fopen(filename,"r");

/*while (!feof(Data))	{*/
for (i=0; i<N; ++i)	{
	fscanf(Data,"%i",&B[i]);
	}
fclose(Data);

return B;
}
