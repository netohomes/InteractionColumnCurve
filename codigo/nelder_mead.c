#include "memoria.h"
#include "matrix.h"

void imprime_iter( const int n, double *f, double **X, char *mensaje ){
    printf( "%s", mensaje );
    for( int i = 0; i < n; i++ ){
        printf( "\n%f", f[i] );
        for( int j = 0; j < n - 1; j++ ) printf( "   %f", X[i][j] );
        printf( "\n" );
    }
}

// Funcion a usar en el Nelder-Mead
double rosenbrock( const int n, double *x ){
    double s = 0.0;
    for( int i = 0; i < n-1; i++ )
        s += 100.0*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]) + (x[i]-1.0)*(x[i]-1.0);
    return s;
}

int particion( const int ini, const int fin, const int n, double *F, double **X ){
    // Hace la particion para el quickSort usando el ultimo elemento como
    // pivote
    int i, j;
    double a, b, c;
    i = ini - 1;
    j = fin;
    a = F[fin];
    while( 1 ){
        while( F[++i] < a );
        while( a < F[--j] ) if( j == ini ) break;
        if( i >= j ) break;
        b    = F[i];
        F[i] = F[j];
        F[j] = b;
        for( int k = 0; k < n; k++ ){
            c = X[i][k];
            X[i][k] = X[j][k];
            X[j][k] = c;
        }
    }
    b      = F[i];
    F[i]   = F[fin];
    F[fin] = b;
    for( int k = 0; k < n; k++ ){
        c = X[i][k];
        X[i][k]   = X[fin][k];
        X[fin][k] = c;
    }
    return i;
}
// Ordena de menor a mayor un conjunto de valores en funcion de F
void quickSort( const int ini, const int fin, const int n, double *F, double **X ){
    int p;
    if( fin <= ini ) return;
    p = particion( ini, fin, n, F, X );
    quickSort( ini, p-1, n, F, X );
    quickSort( p+1, fin, n, F, X );
}

// Algoritmo del Nelder-Mead
void nelder_mead( const int n, double *x, const double delta,
                  const double alpha, const double gamma,
                  const double beta, const double tau,
                  double (*f)( const int n, double *x ) ){
    int conta;
    double **X, *fX, *xR, fxR, *xE, fxE, *xC, fxC, cte, *x_prom;
    // Nota: {x} ya debe contener el punto inicial
    // Pide memoria
    fX = mem_array1D_D( n + 1, "Fitness" );
    X  = mem_array2D_D( n + 1, n, "Puntos del simplejo" );
    xR = mem_array1D_D( n, "{xR}" );
    xE = mem_array1D_D( n, "{xE}" );
    xC = mem_array1D_D( n, "{xC}" );
    x_prom = mem_array1D_D( n, "{x_prom}" );
    // Incluye el {x} inicial en los puntos del simplejo
    for( int j = 0; j < n; j++ ) X[0][j] = x[j];
    // Inicia los demas puntos del simplejo
    for( int i = 1; i < n + 1; i++ ){
        copia_vec_D( n, X[i], x );
        X[i][i-1] += delta;
    }
    conta = 0;
    while( conta < 10000 ){
        // Evalua la funcion en los puntos del simplejo y en {x}
        for( int i = 0; i < n + 1; i++ ) fX[i] = f( n, X[i] );
            //imprime_iter( n+1, fX, X, "\nAntes de ordenar:" );
        // Ordena los puntos de menor a mayor valor de aptitud
        quickSort( 0, n, n, fX, X );
            //imprime_iter( n+1, fX, X, "\nDespues de ordenar:" );
            //getchar();
        // Promedio excluyendo el peor de {Xi} (en este caso el ultimo)
        cte = 1.0 / (double)(n);
        for( int j = 0; j < n; j++ ) x_prom[j] = 0.0;
        for( int i = 0; i < n; i++ )
            for( int j = 0; j < n; j++ ) x_prom[j] += X[i][j];
        for( int j = 0; j < n; j++ ) x_prom[j] *= cte;
        // Calculo del punto reflejado
        for( int j = 0; j < n; j++ )
            xR[j] = x_prom[j] + alpha*(x_prom[j]-X[n][j]);
        fxR = f( n, xR );
        // Actualiza el simplejo:
        // Cuando la reflexion es mejor que el elite actual
        if( fxR < fX[0] ){
            // Se vuelve a reflejar
            for( int j = 0; j < n; j++ )
                xE[j] = xR[j] + gamma*(xR[j]-x_prom[j]);
            fxE = f( n, xE );
            // Se escoge la mejor de las dos reflexiones
            if( fxE < fX[0] ) copia_vec_D( n, X[n], xE );
            else              copia_vec_D( n, X[n], xR );
        }
        // Cuando la reflexion es al menos mejor el penultimo peor
        else if( fxR < fX[n-1] ){
            copia_vec_D( n, X[n], xR );
        }
        // Cuando la reflexion no mejora (contraccion)
        else{
            if( fxR < fX[n] ) copia_vec_D( n, x, xR );
            else              copia_vec_D( n, x, X[n] );
            for( int j = 0; j < n; j++ )
                xC[j] = x_prom[j] + beta*(x[j]-x_prom[j]);
            fxC = f( n, xC );
            if( fxC < fX[n] ) copia_vec_D( n, X[n], xC );
            else{
                for( int i = 1; i < n + 1; i++ ){
                    for( int j = 0; j < n; j++ )
                        X[i][j] = tau*X[0][j] + (1.0-tau)*X[i][j];
                }
            }
        }
        conta++;
    }
    // Copia la solucion obtenida
    copia_vec_D( n, x, X[0] );
    // Libera memoria
    lib_array1D_D( fX );
    lib_array2D_D( n + 1, X );
    lib_array1D_D( xR );
    lib_array1D_D( xE );
    lib_array1D_D( xC );
    lib_array1D_D( x_prom );
}

int main( void ){
    int n = 10;
    double (*f)(  const int, double *x ), *x;

    x = mem_array1D_D( n, "{x}" );
    f = &rosenbrock; 

    /*void nelder_mead( const int n, double *x, const double delta,
                  const double alpha, const double gamma,
                  const double beta, const double tau,
                  double (*f)( const int n, double *x ) ){*/
    
    for( int i = 0; i < n; i++ ) x[i] = -4.0;
    nelder_mead( n, x, 0.1, 1.0, 1.0, 0.5, 0.5, f );
    imp_array1D_D( n, x, "\nSolucion:\n" );
    printf( "\n" );
    
    lib_array1D_D( x );

    return 0;
}