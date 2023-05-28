#include "seccion.h"
#ifndef SOLVER_H
#define SOLVER_H
//=========================================================================
// METODOS DE SOLUCION PARA ENCONTRAR LAS EXCENTRICIDADES
//.........................................................................

void solver( seccion *ss );

void miSolver( seccion *ss );

void busca_ex_ey( seccion *ss, const int excent );

double busca_ey( double theta, seccion *ss, int *sol );

void solCorrecta( double *y, int *sol, const int der1, const int der2,
                  const double y1, const double y2, seccion *ss );

double busca_ex( double theta, seccion *ss, int *sol );

double funcPn( double theta, double desp, seccion *ss );

double funcEy( double theta, double desp, seccion *ss );

double funcEx( double theta, double desp, seccion *ss );

double funcEtotal( double theta, double desp, seccion *ss );

double despFlexPura( double theta, seccion *ss );

int newtonRaphson( double *xx, double ( *f ) ( double,
                   double, seccion* ), const double xIni,
                   const double dx, double theta, seccion *ss );

// Escribe en un archivo el error en las excentricidades
void escribeError( const double deltaT, const double deltaD, const double limInf,
                   const double limSup, seccion *S, const char * const archivo_ex,
                   const char * const archivo_ey, const char * const archivo_e );

// Checa si un flujo a un archivo se abrio
void checa_flujo( FILE *salida, const char * const archivo );

void escribeError_GiroFijo( const double theta, const double deltaD, seccion *S,
                   const char * const archivo_ex, const char * const archivo_ey,
                   const char * const archivo_e );

// Algoritmo del Nelder-Mead
void nelder_mead( const int n, double *x, const double delta,
                  const double alpha, const double gamma,
                  const double beta, const double tau, seccion *S,
                  double (*f)( double, double, seccion* ) );

// Ordena de menor a mayor un conjunto de valores en funcion de F
void quickSort( const int ini, const int fin, const int n, double *F, double **X );
int particion( const int ini, const int fin, const int n, double *F, double **X );

// Calcula un error total
double error_tipo1( double theta, double desp, seccion *ss );
double error_tipo2( double theta, double desp, seccion *ss );

void imprime_iter( const int n, double *f, double **X, char *mensaje );

// Genera la superficie de intaraccion de una seccion
void sup_interaccion( const double deltaT, const double deltaD,
                      const double desp_inf, const double desp_sup,
                      seccion *S, const char * const archivo );

#endif