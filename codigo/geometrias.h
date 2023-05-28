#include "memoria.h"
#include "headers.h"
#include "matrix.h"

#ifndef GEOMETRIAS_H
#define GEOMETRIAS_H
//=========================================================================
// GEOMETRIAS
//.........................................................................

typedef struct{
    int nNodos;       // Numero de vertices
    double a;         // Area
    double cx;        // Centroide en eje x
    double cy;        // Centroide en eje y
    double *x;        // Coordenadas en x
    double *y;        // Coordenadas en y
}poligono;

// Inicia un poligono
void inicia_poligono( poligono *P );

// Pide memoria para un poligono
void mem_poligono( poligono *P, const int n );

// Libera memoria de un poligono
void lib_mem_poligono( poligono *P );

// Copia los datos de un poligono
void copia_poligono( poligono *p, poligono *q );

typedef struct{
    int nAreas;       // Numero de areas de un tipo
    poligono *areas;  // Poligonos que definen a dichas areas
}zona;

// Inicia una zona
void inicia_zona( zona *Z );

// Pide memoria para alojar cierta cantidad de areas en una zona
void mem_ini_zona( zona *Z, const int n );

// Libera memoria de una zona
void lib_mem_zona( zona *Z );

// Copia los valores de una zona a otra
void copia_zona( zona *a, zona *b );

// Imprime los datos de una zona de un material
void imp_zona( zona *Z, const char * const mensaje, const int sangria );

// Obtiene la coord. en "y" mayor de una zona
int yMax_zona( zona *Z, double *ymax );

// Obtiene la coord. en "y" menor de una zona
int yMin_zona( zona *Z, double *ymax );

#endif // GEOMETRIAS_H
