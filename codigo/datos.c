#include "datos.h"
//=========================================================================
// LECTURA Y ESCRITURA
//.........................................................................

// Lee un material tipo concreto
void lee_mat_concr( concreto *M, const char * const palabra_clave,
                    const char * const archivo ){
    int encontrado;
    char cad[300];
    FILE *entrada;
    // Trata de abrir el flujo al archivo
    entrada = fopen( archivo, "r" );
    if( entrada == NULL ){
        printf( "\nError: no se pudo abrir el archivo %s", archivo );
        exit( EXIT_FAILURE );
    }
    // Intenta encontrar el inicio de la definicion del concreto
    encontrado = 0;
    while( fscanf( entrada, "%s", cad ) != EOF ){
        if( strcmp( palabra_clave, cad ) == 0 ){
            encontrado = 1;
            break;
        }
    }
    // Lee la geometria y propiedades
    if( encontrado ){
        inicia_mat_concr( M );
        // Propiedades
        fscanf( entrada, "%s", cad ); fscanf( entrada, "%lf", &M->E );
        fscanf( entrada, "%s", cad ); fscanf( entrada, "%lf", &M->Fc );
        // Zona de areas positivas
        lee_zona( &M->parte_pos, entrada, cad );
        // Zona de areas negativas
        lee_zona( &M->parte_neg, entrada, cad );
    }
    else{
        printf( "\nError: no se encontro la geometria %s", palabra_clave );
    }
    // Imprime en pantalla la lectura
    //strcpy( cad, palabra_clave );
    //strcat( cad, "\0" );
    // Cierra el flujo
    fclose( entrada );
}

// Lee un material tipo acero
void lee_mat_acero( acero *M, const char * const palabra_clave,
                   const char * const archivo ){
    int encontrado;
    char cad[300];
    FILE *entrada;
    // Trata de abrir el flujo al archivo
    entrada = fopen( archivo, "r" );
    if( entrada == NULL ){
        printf( "\nError: no se pudo abrir el archivo %s", archivo );
        exit( EXIT_FAILURE );
    }
    // Intenta encontrar el inicio de la definicion del concreto
    encontrado = 0;
    while( fscanf( entrada, "%s", cad ) != EOF ){
        if( strcmp( palabra_clave, cad ) == 0 ){
            encontrado = 1;
            break;
        }
    }
    // Lee la geometria y propiedades
    if( encontrado ){
        inicia_mat_acero( M );
        // Propiedades
        fscanf( entrada, "%s", cad ); fscanf( entrada, "%lf", &M->E );
        fscanf( entrada, "%s", cad ); fscanf( entrada, "%lf", &M->Fy );
        // Zona de areas positivas
        lee_zona( &M->parte_pos, entrada, cad );
        // Zona de areas negativas
        lee_zona( &M->parte_neg, entrada, cad );
    }
    else{
        printf( "\nError: no se encontro la geometria %s", palabra_clave );
    }
    // Imprime en pantalla la lectura
    //strcpy( cad, palabra_clave );
    //strcat( cad, "\0" );
    // Cierra el flujo
    fclose( entrada );
}

// Lee una zona de un material
void lee_zona( zona *Z, FILE *entrada, char *cad ){
    int aux;
    // Letrero
    fscanf( entrada, "%s", cad );
    // Num. de areas de la zona
    fscanf( entrada, "%d", &aux );
    // Memoria inicial para la zona (cuantas areas o poligonos va a tener)
    mem_ini_zona( Z, aux );
    for( int i = 0; i < Z->nAreas; i++ ){
        // Letrero de area i
        fscanf( entrada, "%s", cad );
        // Numero de nodos
        fscanf( entrada, "%d", &aux );
        // Memoria para definir el area o poligono
        mem_poligono( &Z->areas[i], aux );
        // Lee las coordenadas que definen el area (o poligono)
        for( int j = 0; j < Z->areas[i].nNodos; j++ ){
            fscanf( entrada, "%lf", &Z->areas[i].x[j] ); // En X
            fscanf( entrada, "%lf", &Z->areas[i].y[j] ); // En Y
        }
    }
}

// Lee acciones mecanicas
void lee_acciones_mecanicas( double *P, double *Mx, double *My,
                             const char * const palabra_clave,
                             const char * const archivo ){
    int encontrado;
    char cad[300];
    FILE *entrada;
    // Trata de abrir el flujo al archivo
    entrada = fopen( archivo, "r" );
    if( entrada == NULL ){
        printf( "\nError: no se pudo abrir el archivo %s", archivo );
        exit( EXIT_FAILURE );
    }
    // Intenta encontrar el inicio donde se definen las cargas
    encontrado = 0;
    while( fscanf( entrada, "%s", cad ) != EOF ){
        if( strcmp( palabra_clave, cad ) == 0 ){
            encontrado = 1;
            break;
        }
    }
    // Lee los parametros de la discretizacion
    if( encontrado ){
        fscanf( entrada, "%s", cad ); fscanf( entrada, "%lf", P );
        fscanf( entrada, "%s", cad ); fscanf( entrada, "%lf", Mx );
        fscanf( entrada, "%s", cad ); fscanf( entrada, "%lf", My );
    }
    // Cierra el flujo
    fclose( entrada );
}

// Lee una grafica de esfuerzo-deformacion
void lee_graf_EsfDef( graf_EsfDef *G, const char * const palabra_clave,
                      const char * const archivo ){
    int encontrado, n;
    char cad[300];
    FILE *entrada;
    // Trata de abrir el flujo al archivo
    entrada = fopen( archivo, "r" );
    if( entrada == NULL ){
        printf( "\nError: no se pudo abrir el archivo %s", archivo );
        exit( EXIT_FAILURE );
    }
    // Intenta encontrar el inicio donde se definen los valores
    encontrado = 0;
    while( fscanf( entrada, "%s", cad ) != EOF ){
        if( strcmp( palabra_clave, cad ) == 0 ){
            encontrado = 1;
            break;
        }
    }
    // Lee puntos que definen la grafica
    if( encontrado ){
        fscanf( entrada, "%d", &n ); // Numero de puntos
        mem_graf_EsfDef( G, n );
        for( int i = 0; i < n; i++ ){
            fscanf( entrada, "%lf", &G->def[i] ); // Deformacion
            fscanf( entrada, "%lf", &G->esf[i] ); // Esfuerzo
        }
    }
    // Cierra el flujo
    fclose( entrada );
}
