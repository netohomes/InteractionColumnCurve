#include "geometrias.h"
//=========================================================================
// GEOMETRIAS
//.........................................................................

// Inicia un poligono
void inicia_poligono( poligono *P ){
    P->nNodos = 0;
    P->a  = 0.0;
    P->cx = 0.0;
    P->cy = 0.0;
    P->x  = NULL;
    P->y  = NULL;
}

// Pide memoria para un poligono
void mem_poligono( poligono *P, const int n ){
    P->nNodos = n;
    P->x = mem_array1D_D( n, "{x}" );
    P->y = mem_array1D_D( n, "{y}" );
}

// Libera memoria de un poligono
void lib_mem_poligono( poligono *P ){
    lib_array1D_D( P->x );
    lib_array1D_D( P->y );
    P->nNodos = 0;
}

// Copia los datos de un poligono
void copia_poligono( poligono *p, poligono *q ){
    // Copia propiedades
    p->nNodos = q->nNodos;
    p->a  = q->a;
    p->cx = q->cx;
    p->cy = q->cy;
    // Copia los vertices del poligono
    for( int i = 0; i < q->nNodos; i++ ){
        p->x[i] = q->x[i];
        p->y[i] = q->y[i];
    }
}

// Inicia una zona
void inicia_zona( zona *Z ){
    Z->nAreas = 0;
    Z->areas  = NULL;
}

// Pide memoria para alojar cierta cantidad de areas en una zona
void mem_ini_zona( zona *Z, const int n ){
    Z->nAreas = n;
    Z->areas  = (poligono*) malloc( n * sizeof( poligono ) );
}

// Libera memoria de una zona
void lib_mem_zona( zona *Z ){
    if( Z->areas != NULL ){
        for( int i = 0; i < Z->nAreas; i++ ) lib_mem_poligono( &Z->areas[i] );
        free( Z->areas );
        Z->nAreas = 0;
        Z->areas = NULL;
    }
}

// Copia los valores de una zona a otra
void copia_zona( zona *a, zona *b ){
    // Numero de poligonos
    a->nAreas = b->nAreas;
    // Copia los poligonos
    for( int i = 0; i < b->nAreas; i++ )
        copia_poligono( &a->areas[i], &b->areas[i] );
}

// Imprime los datos de una zona de un material
void imp_zona( zona *Z, const char * const mensaje, const int sangria ){
    char margen[sangria+1];
    for( int i = 0; i < sangria; i++ ) margen[i] = ' ';
    margen[sangria] = '\0';
    printf( "\n%s%s", margen, mensaje );
    printf( "\n%sNo. areas. =  %d", margen, Z->nAreas );
    for( int i = 0; i < Z->nAreas; i++ ){
        printf( "\n%s.. Area No.  %d", margen,  i + 1 );
        for( int j = 0; j < Z->areas[i].nNodos; j++ ){
            printf( "\n%s%.3f", margen, Z->areas[i].x[j] );
            printf( "\t%.3f", Z->areas[i].y[j] );
        }
    }
}

// Obtiene la coord. en "y" mayor de una zona
int yMax_zona( zona *Z, double *ymax ){
    // Regresa 1 si lo encontro, 0 sino
    int cont, encontrado = 0;
    double aux, max = 0.0;
    // Busca un valor inicial valido
    if( Z->nAreas > 0 ){
        cont = 0;
        while( encontrado == 0 ){
            if( Z->areas[cont].nNodos > 0 ){
                max = Z->areas[cont].y[0];
                encontrado = 1;
                break;
            }
            cont++;
        }
    }
    // Busca el mayor
    if( encontrado == 1 ){
        for( int i = cont; i < Z->nAreas; i++ ){
            if( Z->areas[i].nNodos > 0 ){
                max_vec_D( Z->areas[i].nNodos, Z->areas[i].y, & aux );
                if( aux > max ) max = aux;
            }
        }
    }
    *ymax = max;
    return encontrado;
}

// Obtiene la coord. en "y" menor de una zona
int yMin_zona( zona *Z, double *ymin ){
    // Regresa 1 si lo encontro, 0 sino
    int cont, encontrado = 0;
    double aux, min = 0.0;
    // Busca un valor inicial valido
    if( Z->nAreas > 0 ){
        cont = 0;
        while( encontrado == 0 ){
            if( Z->areas[cont].nNodos > 0 ){
                min = Z->areas[cont].y[0];
                encontrado = 1;
                break;
            }
            cont++;
        }
    }
    // Busca el menor
    if( encontrado == 1 ){
        for( int i = cont; i < Z->nAreas; i++ ){
            if( Z->areas[i].nNodos > 0 ){
                min_vec_D( Z->areas[i].nNodos, Z->areas[i].y, &aux );
                if( aux < min ) min = aux;
            }
        }
    }
    *ymin = min;
    return encontrado;
}