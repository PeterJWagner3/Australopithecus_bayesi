#ifdef more_math_functions#include <stdlib.h>#include <math.h>#include <limits.h>int kozrand(int);long combin(int N, int n);double logfact(double );  double lnfact(double );double dvectorproduct(double *v, int l);double dvectorsummation(double *v, int l);unsigned long lngfactorial(int);unsigned long ulngfact(int);int	roundtoint(double);double quadratic_hi(double a, double b, double c);double quadratic_lo(double a, double b, double c);double harmonic_mean(double *data, int N);double geometric_mean(double *data, int N);#elseextern int kozrand(int);										extern long combin(int N, int n);extern double logfact(double );  extern double lnfact(double );extern double dvectorproduct(double *v, int l);extern double dvectorsummation(double *v, int l);extern unsigned long lngfactorial(int);extern unsigned long ulngfact(int);extern int	roundtoint(double);extern double quadratic_hi(double a, double b, double c);extern double quadratic_lo(double a, double b, double c);extern double harmonic_mean(double *data, int N);extern double geometric_mean(double *data, int N);#endif