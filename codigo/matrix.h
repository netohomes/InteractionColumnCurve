#include "headers.h"
#include "memoria.h"
#ifndef MATRIX_H
#define MATRIX_H
//=========================================================================
// FUNCIONES PARA OPERAR MATRICES Y/O VECTORES
//.........................................................................

typedef struct{
	int n;       // Numero de elementos definidos en el vector
	int n_mem;   // Numero de elementos ya alojados en memoria
	int *val;    // Valores de los elementos
}miDeque_I;

// Inicia un vector tipo deque de C++
void inicia_miDeque_I( miDeque_I *v );

// Memoria inicial
void mem_miDeque_I( miDeque_I *v, const int n );

// Agrega un elemento al final. La idea es alojar el doble de memoria una
// vez que se llena el vector.
void add_miDeque_I( miDeque_I *v, const int a );

// Libera memoria de miDeque
void lib_mem_miDeque_I( miDeque_I *v );

typedef struct{
	int n;       // Numero de elementos definidos en el vector
	int n_mem;   // Numero de elementos ya alojados en memoria
	double *val; // Valores de los elementos
}miDeque_D;

// Inicia un vector tipo deque de C++
void inicia_miDeque_D( miDeque_D *v );

// Memoria inicial
void mem_miDeque_D( miDeque_D *v, const int n );

// Agrega un elemento al final. La idea es alojar el doble de memoria una
// vez que se llena el vector.
void add_miDeque_D( miDeque_D *v, const double a );

// Libera memoria de miDeque
void lib_mem_miDeque_D( miDeque_D *v );

// Imprime valores de miDeque
void imp_miDeque_D( miDeque_D *v, const char *mensaje );

// Copia un vector en otro {a} = {b}
void copia_vec_D( const int n, double *a, double *b );

// Obtienen el valor maximo y minimo de un vector
int max_vec_D( const int n, double *v, double *max );
int min_vec_D( const int n, double *v, double *min );

// Obtienen el valor maximo y minimo de un vector
int max_vec_I( const int n, int *v, int *max );
int min_vec_I( const int n, int *v, int *min );

// Multiplicacion de matrices [C] = [A][B]
void mult_mat_mat_D( const int m, const int n, const int p, double **C,
                         double **A, double **B );

// {b} = [A]{v}
void mult_mat_vec_D( const int m, const int n, double *b, double **A, double *v );

// [B] = [A]^T
void transpuesta_D( const int m, const int n, double **B, double **A );

// [B] = [A]^T[A]
void mult_matT_mat_D( const int m, const int n, double **A, double **B );

// {c} = [A]^T {b}. [A] es de tamano m x n, {b} es de tam. n x 1
void mult_matT_vec_D( const int m, const int n, double *c, double **A, double *b );

// Calcula la matriz de rotacion en 2D para el angulo de giro dado
void mat_rotacion_2D( const double theta, double **R );

// Invierte una matriz de 2 x 2
void inv_mat_2x2( double **R, double **RI );

// Obtiene el maximo o minimo de una columna de una matriz
double min_col_j_mat( const int j, const int m, double **A );
double max_col_j_mat( const int j, const int m, double **A );


#endif // MATRIX_H
