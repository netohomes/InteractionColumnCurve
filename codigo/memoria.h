#include "headers.h"

#ifndef MEMORIA_H
#define MEMORIA_H
//=========================================================================
// FUNCIONES PARA MANEJAR MEMORIA
//.........................................................................

// Pide memoria para un vector de doubles
double *mem_array1D_D( const size_t n, const char *nombre );

// Libera la memoria de un vector de doubles
void lib_array1D_D( double *v );

// Pide memoria para un vector de enteros
int *mem_array1D_I( const size_t n, const char *nombre );

// Libera la memoria de un vector de enteros
void lib_array1D_I( int *v );

// Pide memoria para una matriz de m x n de doubles
double **mem_array2D_D( const size_t m, const size_t n,
                        const char *nombre );

// Libera la memoria de una matriz de m x n
void lib_array2D_D( const size_t m, double **A );

//=========================================================================
// FUNCIONES PARA IMPRIMIR ARREGLOS
//.........................................................................

// Imprime los valores de un vector de enteros
void imp_array1D_I( const int n, const int *v, const char *mensaje );

// Imprime los valores de un vector de doubles
void imp_array1D_D( const int n, const double *v, const char *mensaje );

// Imprime los valores de una matriz de doubles
void imp_array2D_D( const int m, const int n, double **A,
                    const char *mensaje );
#endif // MEMORIA_H
