#include "datos.h"
#include "seccion.h"
#include "materiales.h"
#include "solver.h"
const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481;
int main( int argc, char *argv[] ){
    
    double P, Mx, My, dh;
    seccion S;

    // Inicia problema
    inicia_seccion( &S );
    lee_mat_concr( &S.region_concr, "CONCRETO", argv[5] );
    lee_mat_acero( &S.region_acero, "ACERO", argv[5] );
    lee_graf_EsfDef( &S.concrED, "VALORES", "../Graf_esf_def/EsfDefConcResult250.dat" );
    lee_graf_EsfDef( &S.aceroED, "VALORES", "../Graf_esf_def/EsfDefAceroResult.dat" );
    inicia_val( &S );
        //imp_graf_EsfDef( &S.concrED, "\nEsf-Def Concreto" );
        //imp_graf_EsfDef( &S.aceroED, "\nEsf-Def Acero" );

    // Acciones mecanicas externas
    P  = atof( argv[1] );
    Mx = atof( argv[2] );
    My = atof( argv[3] );
    // Tamano de la discretizacion
    dh = atof( argv[4] );

    // Fija los valores
    fija_dh( dh, &S );
    fija_cargas( P, Mx, My, &S );

    /*
    escribeError( 0.2, 0.2, -40.0, 40.0, &S, "resultados/errores_ex.txt", 
                 "resultados/errores_ey.txt", "resultados/errores_e.txt" );
    //*/
    
    /*
    escribeError_GiroFijo( atof(argv[1])/6.0*pi, 0.1, &S, "resultados/errores_ex_GiroFijo.txt",
                           "resultados/errores_ey_GiroFijo.txt",
                           "resultados/errores_e_GiroFijo.txt" );
    */

    //sup_interaccion( 0.2, 0.2, -50.0, 50.0, &S, "resultados/SI_seccion_compuesta_circular.txt" );
    
    // Evalua la resistencia de la seccion dados una rotacion y un desplazamiento
    //eval_resist_seccion( 0.0, 0.0, &S );

    // Corre el solver
    solver( &S );
    //imp_seccion( &S, "\nDatos de la seccion:" ); printf( "\n" );
    
    // Libera la memoria
    lib_mem_seccion( &S );
    
    return 0;
 
 }
