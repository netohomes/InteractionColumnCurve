#include "materiales.h"
//=========================================================================
// MATERIALES Y SUS PROPIEDADES
//.........................................................................

// Inicia valores para un material tipo concreto
void inicia_mat_concr( concreto *M ){
    M->E  = 0.0;
    M->Fc = 0.0;
    inicia_zona( &M->parte_neg );
    inicia_zona( &M->parte_pos );
}

// Inicia la memoria para un material tipo concreto
void mem_ini_mat_concr( concreto *M, const int nAreas_neg, const int nAreas_pos ){
    mem_ini_zona( &M->parte_neg, nAreas_neg );
    mem_ini_zona( &M->parte_pos, nAreas_pos );
}

// Libera la memoria de un material tipo concreto
void lib_mem_mat_concr( concreto *M ){
    lib_mem_zona( &M->parte_neg );  
    lib_mem_zona( &M->parte_pos );
}

// Hace una copia de un material de concreto
void copia_mat_concr( concreto *a, concreto *b ){
    // Propiedades
    a->E  = b->E;
    a->Fc = b->Fc;
    // Areas negativas y positivas
    copia_zona( &a->parte_pos, &b->parte_pos );
    copia_zona( &a->parte_neg, &b->parte_neg );
}

// Imprime los datos de un material tipo concreto
void imp_mat_concr( concreto *M, const char * const mensaje,
                    const int sangria ){
    char margen[sangria+1];
    for( int i = 0; i < sangria; i++ ) margen[i] = ' ';
    margen[sangria] = '\0';
    printf( "%s", mensaje );
    // Propiedades
    printf( "\n%sE  = %f", margen, M->E );
    printf( "\n%sFc = %f", margen, M->Fc );
    // Zona positiva
    imp_zona( &M->parte_pos, "..... Zona positiva", sangria );
    // Zona negativa
    imp_zona( &M->parte_neg, "..... Zona negativa", sangria );
}

// Inicia valores para un material tipo acero
void inicia_mat_acero( acero *M ){
    M->E  = 0.0;
    M->Fy = 0.0;
    inicia_zona( &M->parte_neg );
    inicia_zona( &M->parte_pos );
}

// Inicia la memoria para un material tipo acero
void mem_ini_mat_acero( acero *M, const int nAreas_neg, const int nAreas_pos ){
    mem_ini_zona( &M->parte_neg, nAreas_neg );
    mem_ini_zona( &M->parte_pos, nAreas_pos );
}

// Libera la memoria de un material tipo acero
void lib_mem_mat_acero( acero *M ){
    lib_mem_zona( &M->parte_neg );
    lib_mem_zona( &M->parte_pos );
}

// Hace una copia de un material de acero
void copia_mat_acero( acero *a, acero *b ){
    // Propiedades
    a->E  = b->E;
    a->Fy = b->Fy;
    // Areas negativas y positivas
    copia_zona( &a->parte_pos, &b->parte_pos );
    copia_zona( &a->parte_neg, &b->parte_neg );
}

// Imprime los datos de un material tipo acero
void imp_mat_acero( acero *M, const char * const mensaje,
                    const int sangria ){
    char margen[sangria+1];
    for( int i = 0; i < sangria; i++ ) margen[i] = ' ';
    margen[sangria] = '\0';
    printf( "%s", mensaje );
    // Propiedades
    printf( "\n%sE  = %f", margen, M->E );
    printf( "\n%sFy = %f", margen, M->Fy );
    // Zona positiva
    imp_zona( &M->parte_pos, "..... Zona positiva", sangria );
    // Zona negativa
    imp_zona( &M->parte_neg, "..... Zona negativa", sangria );
}

// Inicia una grafica de esfuerzo-deformacion
void inicia_graf_EsfDef( graf_EsfDef *G ){
    G->n = 0;
    G->def = NULL;
    G->esf = NULL;
}

// Pide memoria para una grafica esf-def
void mem_graf_EsfDef( graf_EsfDef *G, const int n ){
    G->def = mem_array1D_D( n, "deformaciones" );
    G->esf = mem_array1D_D( n, "esfuerzos" );
    G->n = n;
}

// Libera memoria de una grafica de esf-def
void lib_mem_graf_EsfDef( graf_EsfDef *G ){
    lib_array1D_D( G->def );
    lib_array1D_D( G->esf );
    G->n = 0;
}

// Imprime grafica de esfuerzo-deformacion
void imp_graf_EsfDef( graf_EsfDef *G, const char * const mensaje ){
    printf( "%s", mensaje );
    printf( "\nNo. puntos = %d", G->n );
    for( int i = 0; i < G->n; i++ )
        printf( "\n%f\t%f", G->def[i], G->esf[i] );
}


