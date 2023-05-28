#include "memoria.h"
//=========================================================================
// FUNCIONES PARA MANEJAR MEMORIA
//.........................................................................

//.........................................................................
// Arreglos 1D
//.........................................................................

// Pide memoria para un vector de doubles
double *mem_array1D_D( const size_t n, const char *nombre ){
	double *v;
	if( n == 0 ) return NULL;
	v = ( double* ) malloc( n * sizeof( double ) );
	if( v == NULL ){
		printf( "\nMem. insuficiente para el vector %s", nombre );
		exit( EXIT_FAILURE );
	}
	return v;
}

// Libera la memoria de un vector de doubles
void lib_array1D_D( double *v ){
	if( v != NULL ) free( v );
	v = NULL;
}

// Pide memoria para un vector de enteros
int *mem_array1D_I( const size_t n, const char *nombre ){
	int *v;
	if( n <= 0 ) return NULL;
	v = ( int* ) malloc( n * sizeof( int ) );
	if( v == NULL ){
		printf( "\nMem. insuficiente para el vector %s", nombre );
		exit( EXIT_FAILURE );
	}
	return v;
}

// Libera la memoria de un vector de enteros
void lib_array1D_I( int *v ){
	if( v != NULL ) free( v );
	v = NULL;
}

//.........................................................................
// Arreglos 2D
//.........................................................................

// Pide memoria para una matriz de m x n de doubles
double **mem_array2D_D( const size_t m, const size_t n, const char *nombre ){
	double **ptr;
	// Cuando "m" es menor o igual a cero
	if( m == 0 ) return NULL;
	ptr = ( double** ) malloc( m * sizeof( double* ) );
	if( ptr == NULL ){
		printf ( "\nMem. insuficiente para la matriz %s", nombre );
		exit( EXIT_FAILURE );
	}
	else{
		// Cuando "n" es menor o igual a cero
		if( n == 0 ){
			for( size_t i = 0; i < m; i++ ) ptr[i]	= NULL;
		}
		else{
			for( size_t i = 0; i < m; i++ ){
				ptr[i] = ( double* ) malloc( n * sizeof ( double ) );
				if ( ptr[i] == NULL ){
					printf( "\nMem. insuficiente para la matriz %s", nombre );
					exit( EXIT_FAILURE );
				}
			}
		}
	}
	return ptr;
}

// Libera la memoria de una matriz de m x n
void lib_array2D_D( const size_t m, double **A ){
	if( A != NULL ){
		for( size_t i = 0; i < m; i++ ){
			if( A[i] != NULL ) free( A[i] );
			A[i] = NULL;
		}
	}
	free( A );
	A = NULL;
}

//=========================================================================
// FUNCIONES PARA IMPRIMIR ARREGLOS
//.........................................................................

// Imprime los valores de un vector de enteros
void imp_array1D_I( const int n, const int *v, const char *mensaje ){
	printf( "%s", mensaje );
	for( int i = 0; i < n; i++ ) printf( "%d   ", v[i] );
	printf( "\n" );
}

// Imprime los valores de un vector de doubles
void imp_array1D_D( const int n, const double *v, const char *mensaje ){
	printf( "%s", mensaje );
	for( int i = 0; i < n; i++ ) printf( "%lf\n", v[i] );
	printf( "\n" );
}

// Imprime los valores de una matriz de doubles
void imp_array2D_D( const int m, const int n, double **A,
                         const char *mensaje ){
	printf( "%s", mensaje );
	for( int i = 0; i < m; i++ ){
		for( int j = 0; j < n; j++ ){
			printf( "%.9lf   ", A[i][j] );
		}
		printf( "\n" );
	}
	printf( "\n" );
}

