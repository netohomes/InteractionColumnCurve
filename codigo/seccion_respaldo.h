#include "materiales.h"
#include "matrix.h"
#include "memoria.h"
#ifndef SECCION_H
#define SECCION_H
//=========================================================================
// ELEMENTOS SECCION Y FUNCIONES PARA DETERMINAR SU RESISTENCIA
//.........................................................................

typedef struct{
    int nFranjas;   // Numero de franjas para un area (poligono) dado
    double *a;  // Arreglo con las areas de cada franja
    double *cx; // Arreglo con el centroide en x de cada franja
    double *cy; // Arreglo con el centroide en y de cada franja
}discr_pol;

// Inicia la discretizacion de un area (poligono)
void inicia_discr_pol( discr_pol *D );

// Pide memoria para la discretizacion
void mem_discr_pol( discr_pol *D, const int n );

// Libera memoria de discretizacion
void lib_mem_discr_pol( discr_pol *D );

typedef struct{
    int nAreas;     // Numero de areas (poligonos) que han sido discretizados
    discr_pol *pol; // Discretizacion de cada area
}discretizacion;

// Inicia memoria para discretizacion
void inicia_mem_discretizacion( discretizacion *D, const int n );

// Libera memoria de discretizacion
void lib_mem_discretizacion( discretizacion *D );

typedef struct{
    int nDiv;
    double *a;   // La suma de areas pos. y neg.
    double *cx;  // El centroide en "x"
    double *cy;  // El centroide en "y"
    double *def; // La deformacion dado cy con respecto al EN
    double *esf; // El esfuerzo de cada franja en cy
    double *P;   // Contiene la suma de fuerzas de compresion y/o tension
    double *Mx;  // Suma de momentos en eje "x"
    double *My;  // Suma de momentos en eje "y"
}tabla_result;

// Inicia tabla de resultados
void inicia_tabla_result( tabla_result *T );

// Memoria para tabla de resultados
void mem_tabla_result( tabla_result *T, const int n );

// Inicia en ceros todos los valores de la tabla
void llenar_ceros( tabla_result *T );

// Libera memoria de tabla de resultados
void lib_mem_tabla_result( tabla_result *T );

// Imprime tabla de resutlados
void imp_tabla_result( tabla_result *T, const char * const mensaje );

typedef struct{
	// Propiedades de la seccion
	double cx;  // Centroide de la seccion en X
	double cy;  // Centroide de la seccion en Y
	double a;   // Area de la seccion
    // Acciones mecanicas externas
    double P;   // Fuerza axial
    double Mx;  // Momento sobre el eje x
    double My;  // Momento sobre el eje y
    // Resistencia nominal
    double Pn;  // Resistencia a la compresion o tension
    double Mnx; // Resistencia alrededor del eje x de la seccion
    double Mny; // Resistencia alrededor del eje y de la seccion
    // Ubicacion del Eje Neutro
    double theta;// Angulo de giro del EN
    double desp; // Desplazamiento del EN
    // Excentricidades
    double exB, eyB; // Las que se estan buscando
    double ex, ey;   // Las actuales
    // Tamano de discretizacion
    double dh;
    // 
    // Nota: se considera simepre que primero se rota el EN y luego se desplaza.
    // El desplazamiento es en direccion perpendicular al EN.
    concreto region_concr; // Define las areas llenas o huecas de concreto y sus props.
    acero    region_acero; // Define las areas llenas o huecas de acero y sus props.
    concreto region_concr_desp; // Elementos auxiliares donde se rotan y desplazan
    acero    region_acero_desp; // las coordenadas
    // Discretizacion zona de compresion
    discretizacion discr_concr_pos;
    discretizacion discr_concr_neg;
    discretizacion discr_acero_pos;
    discretizacion discr_acero_neg;
    // Discretizacion zona de tension
    discretizacion discr_acero_posT;
    discretizacion discr_acero_negT;
    // Graficas de esfuerzo-deformacion de los materiales
    graf_EsfDef aceroED, concrED;
}seccion;

// Inicia una seccion
void inicia_seccion( seccion *S );

// Libera memoria de una seccion
void lib_mem_seccion( seccion *S );

// Imprime los datos de una seccion
void imp_seccion( seccion *S, const char * const mensaje );

// Inicia valores de la seccion
void inicia_val( seccion *S );

// Evalua la resistencia de una seccion de concreto y acero en funcion
// de la posicion del EN( theta, desp )
void eval_resist_seccion( const double theta, const double desp, seccion *S );

// Rota y luego desplaza la seccion para que el EN quede horizontal
// y pase por el origen (alineado con el eje de las X)
void rota_desplaza( const double theta,const double desp, concreto *cr,
                    acero *ac );

// Rota y desplaza la geometria de un material de tal forma que el EN
// quede horizontal y pase por el origen, esto en funcion de las
// variables theta y desp representan el angulo de giro y el
// desplazamiento del EN.
void rota_despl_poligono( const double theta, const double desp, poligono *geom );

// Rota una serie de coordenadas dado un angulo de giro en 2D
void rota( const double theta, const int n, double *gx, double *gy );

// Desplaza una serie de coordenadas en 2D
void desplaza( const double desp_x, const double desp_y, const int n,
               double *x, double *y );

// Pide la memoria faltante para una seccion
void mem_adicional( seccion *S );

// Aloja nueva memoria para guardar geometrias del mismo tamano
void copia_mem_geom( zona *a,  zona *b );

// Calcula el centroide y area de la seccion
void calc_props_seccion( seccion *S );

// Calcula los centroides y areas de una zona de la seccion
void calc_props_geom( zona *z );

// Suma las areas y primeros momentos de una zona de la seccion
void sum_props_geom( zona *z, double *are, double *mx, double *my,
                     const double f );

// Calcula el area de un poligono definido por sus coordenadas
double area_pol2D( const int n, double *x, double *y );

// Calcula el centroide de un poligono definido por sus coords.
double centroide_x_pol2D( const double area, const int n, double *x, double *y );
double centroide_y_pol2D( const double area, const int n, double *x, double *y );

// Discretiza la zona a compresion de la seccion
void discretiza_comp( const double dh, concreto *concr, acero *ace,
                      discretizacion *discr_concr_pos,
                      discretizacion *discr_concr_neg,
                      discretizacion *discr_acero_pos,
                      discretizacion *discr_acero_neg );

// Discretiza la zona a tension de la seccion
void discretiza_tension( const double dh, acero *ace,
                         discretizacion *discr_acero_pos,
                         discretizacion *discr_acero_neg );

// Discretiza una zona en tension
void discretiza_zona_tension( const double dh, const double yini, zona *Z, discretizacion *discr );

// Discretiza una zona en compresion
void discretiza_zona_comp( const double dh, const double yini, zona *Z, discretizacion *discr );

// Determina el area y centroide de una franja sobre un poligono
void area_centro_franja( const double y_inf, const double y_sup,
                         poligono *P, double *a, double *cx, double *cy );

// Determina los nodos que definen una franja (conectividad)
void conectividad_franja( const double y_inf, const double y_sup, 
                          poligono *P, miDeque_D *xF, miDeque_D *yF );

// Agrega interseccion
void agrega_inter( double x, double y, miDeque_D *xF, miDeque_D *yF );

// Verifica si un lado y sus intersecciones estan dentro de la franja
void agrega_lado_inter( const double y_inf, const double y_sup, const int i,
                        double xp, double yp, double xi, double yi, int li,
                        double xs, double ys, int ls, int *contI, int *contS,
                        miDeque_D *xF, miDeque_D *yF );

// Calcula la distancia euclideana entre dos puntos en 2D
double dist( double P1x, double P1y, double P2x, double P2y );

// Calcula las intersecciones entre un poligono y una recta horizonal
void intersecciones( const double a, const int np, double *px, double *py,
                     miDeque_D *inter_x, miDeque_D *inter_y, miDeque_I *inter_l );

// Determina el punto de interseccion entre dos rectas definidas
// por dos puntos cada una
void intersecta( const double L1x, const double L1y,
                 const double L2x, const double L2y,
                 const double P1x, const double P1y,
                 const double P2x, const double P2y,
                 int *intersecta, double *x, double *y );

// Determina si un punto esta contenido en un segmentos de recta
int punto_segmento( const double L1x, const double L1y,
                    const double L2x, const double L2y,
                   double x, double y );

// Calcula la pendiente de la recta y determina si la recta es
// una horizontal o una vertical
double pendiente_recta( const double x1, const double y1,
                        const double x2, const double y2,
                        int *hor, int *vert );

// Calcula fuerza y momento resistentes
void resist_nominal( const double desp,
                     const double defMaxConcr,
                     const double defMaxAcero,
                     concreto *concr_desp,
                     acero *acero_desp,
                     discretizacion *discr_concr_pos,
                     discretizacion *discr_concr_neg,
                     discretizacion *discr_acero_pos,
                     discretizacion *discr_acero_neg,
                     discretizacion *discr_acero_posT,
                     discretizacion *discr_acero_negT,
                     graf_EsfDef *concrED,
                     graf_EsfDef *aceroED,
                     double *Pn, double *Mnx, double *Mny );

int max_nFranjas( discretizacion *discr );

void llena_tabla( discretizacion *discr, tabla_result *T, const double f );

// Calcula el esfuerzo en base a las graficas de esfuerzo-defor-
// macion dada un deformacion unitaria
double esfuerzo( const double def, graf_EsfDef *G );

// Calcula las excentricidades de acuerdo a los momentos y la fuerza axial
void calcExcent( double *ex, double *ey, double P, double Mx, double My );

// Devuelve a los ejes originales los momentos
void ejes_ini( const double theta, double *Mnx, double *Mny );

void imp_conectividad( miDeque_D *x, miDeque_D *y, const char * const mensaje );

// Fija el valor del tamano de la discretizacion
void fija_dh( const double dh, seccion *S );

// Fija acciones mecanicas externas
void fija_cargas( const double P, const double Mx, const double My, seccion *S );

// Calcula la deformacion a una altura dada
double deformacion( const double y, const double c, const double a,
                    const double defMaxConcr, const double defMaxAcero );

#endif // SECCION_H

