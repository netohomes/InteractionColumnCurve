#include "geometrias.h"

#ifndef MATERIALES_H
#define MATERIALES_H
//=========================================================================
// MATERIALES Y SUS PROPIEDADES
//.........................................................................

typedef struct{
    double E;
    double Fc;
    zona parte_neg;  // Areas negativas
    zona parte_pos;  // Areas positivas
}concreto;

// Inicia valores para un material tipo concreto
void inicia_mat_concr( concreto *M );

// Inicia la memoria para un material tipo concreto
void mem_ini_mat_concr( concreto *M, const int nAreas_neg, const int nAreas_pos );

// Libera la memoria de un material tipo concreto
void lib_mem_mat_concr( concreto *M );

// Hace una copia de un material de concreto
void copia_mat_concr( concreto *a, concreto *b );

// Imprime los datos de un material tipo concreto
void imp_mat_concr( concreto *M, const char * const mensaje,
                    const int sangria );

typedef struct{
    double E;
    double Fy;
    zona parte_neg;  // Areas negativas
    zona parte_pos;  // Areas positivas
}acero;

// Inicia valores para un material tipo acero
void inicia_mat_acero( acero *M );

// Inicia la memoria para un material tipo acero
void mem_ini_mat_acero( acero *M, const int nAreas_neg, const int nAreas_pos );

// Libera la memoria de un material tipo acero
void lib_mem_mat_acero( acero *M );

// Hace una copia de un material de acero
void copia_mat_acero( acero *a, acero *b );

// Imprime los datos de un material tipo acero
void imp_mat_acero( acero *M, const char * const mensaje,
                    const int sangria );

typedef struct{
    int n;       // Numero de puntos que definen la grafica
    double *def; // Deformacion (eje x)
    double *esf; // Esfuerzo (eje y)
}graf_EsfDef;

// Inicia una grafuca de esfuerzo-deformacion
void inicia_graf_EsfDef( graf_EsfDef *G );

// Pide memoria para una grafica esf-def
void mem_graf_EsfDef( graf_EsfDef *G, const int n );

// Libera memoria de una grafica de esf-def
void lib_mem_graf_EsfDef( graf_EsfDef *G );

void imp_graf_EsfDef( graf_EsfDef *G, const char * const mensaje );

#endif // MATERIALES_H
