#include "matrix.h"
//=========================================================================
// FUNCIONES PARA OPERAR MATRICES Y/O VECTORES
//.........................................................................

// Inicia un vector tipo deque de C++
void inicia_miDeque_D( miDeque_D *v ){
	v->n     = 0;
	v->n_mem = 0;
	v->val   = NULL;
}

// Memoria inicial
void mem_miDeque_D( miDeque_D *v, const int n ){
	if( n > 0 ){
		v->n = n;
		v->n_mem = n;
		v->val = (double*) malloc( v->n_mem * sizeof(double) );
		if( v->val == NULL ){
			printf( "\nError: no se pudo pedir memoria miDeque" );
			v->n = 0;
			v->n_mem = 0;
			v->val = NULL;
		}
	}
}

// Agrega un elemento al final. La idea es alojar el doble de memoria una
// vez que se llena el vector.
void add_miDeque_D( miDeque_D *v, const double a ){
	double *ptr;
	// Cuando se llena
	if( v->n == v->n_mem ){
		v->n += 1;
		v->n_mem = 2 * v->n;
		ptr = (double*) realloc( v->val, v->n_mem * sizeof(double) );
		if( ptr == NULL ){
			printf( "\nError: no se pudo pedir memoria miDeque" );
			if( v->val != NULL ) free( v->val );
			v->n     = 0;
			v->n_mem = 0;
			v->val   = NULL;
		}
		else{
			v->val = ptr;
			v->val[v->n-1] = a;
		}
	}
	else{
		v->n += 1;
		v->val[v->n-1] = a;
	}
}

// Libera memoria de miDeque
void lib_mem_miDeque_D( miDeque_D *v ){
	if( v->val != NULL ){
		free( v->val );
		v->val = NULL;
		v->n = 0;
		v->n_mem = 0;
	}
}

// Imprime valores de miDeque
void imp_miDeque_D( miDeque_D *v, const char *mensaje ){
	printf( "%s", mensaje );
	printf( "\nn = %d", v->n );
	imp_array1D_D( v->n, v->val, "\nValores:\n" );
}

// Copia un vector en otro {a} = {b}
void copia_vec_D( const int n, double *a, double *b ){
	for( int i = 0; i < n; i++ ) a[i] = b[i];
}

// Inicia un vector tipo deque de C++
void inicia_miDeque_I( miDeque_I *v ){
	v->n     = 0;
	v->n_mem = 0;
	v->val   = NULL;
}

// Memoria inicial
void mem_miDeque_I( miDeque_I *v, const int n ){
	if( n > 0 ){
		v->n = n;
		v->n_mem = n;
		v->val = (int*) malloc( v->n_mem * sizeof(int) );
		if( v->val == NULL ){
			printf( "\nError: no se pudo pedir memoria miDeque" );
			v->n = 0;
			v->n_mem = 0;
			v->val = NULL;
		}
	}
}

// Agrega un elemento al final. La idea es alojar el doble de memoria una
// vez que se llena el vector.
void add_miDeque_I( miDeque_I *v, const int a ){
	int *ptr;
	// Cuando se llena
	if( v->n == v->n_mem ){
		v->n += 1;
		v->n_mem = 2 * v->n;
		ptr = (int*) realloc( v->val, v->n_mem * sizeof(int) );
		if( ptr == NULL ){
			printf( "\nError: no se pudo pedir memoria miDeque" );
			if( v->val != NULL ) free( v->val );
			v->n     = 0;
			v->n_mem = 0;
			v->val   = NULL;
		}
		else{
			v->val = ptr;
			v->val[v->n-1] = a;
		}
	}
	else{
		v->n += 1;
		v->val[v->n-1] = a;
	}
}

// Libera memoria de miDeque
void lib_mem_miDeque_I( miDeque_I *v ){
	if( v->val != NULL ){
		free( v->val );
		v->val = NULL;
		v->n = 0;
		v->n_mem = 0;		
	}
}

// Obtienen el valor maximo y minimo de un vector
int max_vec_D( const int n, double *v, double *max ){
	int encontrado = 0;
	if( n > 0 ){
		*max = v[0];
		for( int i = 0; i < n; i++ ) if( *max < v[i] ) *max = v[i];
		encontrado = 1;
	}
	return encontrado;
}
int min_vec_D( const int n, double *v, double *min ){
	int encontrado = 0;
	if( n > 0 ){
		*min = v[0];
		for( int i = 0; i < n; i++ ) if( *min > v[i] ) *min = v[i];
		encontrado = 1;
	}
	return encontrado;
}

// Obtienen el valor maximo y minimo de un vector
int max_vec_I( const int n, int *v, int *max ){
	int encontrado = 0;
	if( n > 0 ){
		*max = v[0];
		for( int i = 0; i < n; i++ ) if( *max < v[i] ) *max = v[i];
		encontrado = 1;
	}
	return encontrado;
}
int min_vec_I( const int n, int *v, int *min ){
	int encontrado = 0;
	if( n > 0 ){
		*min = v[0];
		for( int i = 0; i < n; i++ ) if( *min > v[i] ) *min = v[i];
		encontrado = 1;
	}
	return encontrado;
}


// Multiplicacion de matrices [C] = [A][B]
void mult_mat_mat_D( const int m, const int n, const int p, double **C,
                         double **A, double **B ){
    double sum;
    for( int i = 0; i < m; i++ ){
        for( int j = 0; j < p; j++ ){
            sum = 0.0;
            for( int k = 0; k < n; k++ ) sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
    }
}

// {b} = [A]{v}
void mult_mat_vec_D( const int m, const int n, double *b, double **A, double *v ){

    double sum;
    for( int i = 0; i < m; i++ ){
        sum = 0.0;
        for( int j = 0; j < n; j++ ) sum += A[i][j] * v[j];
        b[i] = sum;
    }
}

// [B] = [A]^T
void transpuesta_D( const int m, const int n, double **B, double **A ){
	for( int i = 0; i < m; i++ )
		for( int j = 0; j < n; j++ ) B[j][i] = A[i][j];
}

// [B] = [A]^T[A]
void mult_matT_mat_D( const int m, const int n, double **A, double **B ){

	double sum;
	for( int i = 0; i < n; i++ ){
		for( int j = 0; j < n; j++ ){
			sum = 0.0;
			for( int k = 0; k < m; k++ ) sum += A[k][i] * A[k][j];
			B[i][j] = sum;
		}
	}
}

// {c} = [A]^T {b}. [A] es de tamano m x n, {b} es de tam. n x 1
void mult_matT_vec_D( const int m, const int n, double *c, double **A, double *b ){

	double sum;
	for( int i = 1; i < n; i++ ){
		sum = 0.0;
		for( int k = 0; k < m; k++ ) sum += A[k][i] * b[k];
		c[i] = sum;
	}
}

// Calcula la matriz de rotacion en 2D para el angulo de giro dado
void mat_rotacion_2D( const double theta, double **R ){
	R[0][0] =  cos( theta );
	R[0][1] = -sin( theta );
	R[1][0] =  sin( theta );
	R[1][1] =  cos( theta );
}

// Invierte una matriz de 2 x 2
void inv_mat_2x2( double **R, double **RI ){
	double det;
    det =   R[0][0] * R[1][1]
          - R[0][1] * R[1][0];
    RI[0][0] =  R[1][1] / det;
    RI[1][1] =  R[0][0] / det;
    RI[0][1] = -R[0][1] / det;
    RI[1][0] = -R[1][0] / det;	
}

// Obtiene el maximo o minimo de una columna j de una matriz
double max_col_j_mat( const int j, const int m, double **A ){
	double max = A[0][j];
    for( int i = 0; i < m; i++ ) if( max < A[i][j] ) max = A[i][j];
    return max;
}
double min_col_j_mat( const int j, const int m, double **A ){
	double min = A[0][j];
    for( int i = 0; i < m; i++ ) if( min > A[i][j] ) min = A[i][j];
    return min;
}

//=========================================================================
