#include "headers.h"
#include "materiales.h"

#ifndef DATOS_H
#define DATOS_H
//=========================================================================
// LECTURA Y ESCRITURA DE DATOS
//.........................................................................

// Lee un material tipo concreto
void lee_mat_concr( concreto *M, const char * const palabra_clave,
                   const char * const archivo );

// Lee un material tipo acero
void lee_mat_acero( acero *M, const char * const palabra_clave,
                   const char * const archivo );

// Lee una zona de un material
void lee_zona( zona *Z, FILE *entrada, char *cad );

// Imprime los datos de un material tipo concreto
void imp_mat_concr( concreto *M, const char * const mensaje,
                    const int sangria );

// Imprime los datos de un material tipo acero
void imp_mat_acero( acero *M, const char * const mensaje,
                    const int sangria );

// Imprime los datos de una zona de un material
void imp_zona( zona *Z, const char * const mensaje, const int sangria );

// Lee acciones mecanicas
void lee_acciones_mecanicas( double *P, double *Mx, double *My,
                             const char * const palabra_clave,
                             const char * const archivo );

// Lee una grafica de esfuerzo-deformacion
void lee_graf_EsfDef( graf_EsfDef *G, const char * const palabra_clave,
                      const char * const archivo );

#endif // DATOS_H
