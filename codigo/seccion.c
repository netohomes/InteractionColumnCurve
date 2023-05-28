#include "seccion.h"

//=========================================================================
// FUNCIONES PARA DETERMINAR LA RESISTENCIA DE UNA SECCION
//.........................................................................

// Inicia la discretizacion de un area (poligono)
void inicia_discr_pol( discr_pol *D ){
    D->nFranjas = 0;
    D->a  = NULL;
    D->cx = NULL;
    D->cy = NULL;
}

// Pide memoria para la discretizacion
void mem_discr_pol( discr_pol *D, const int n ){
    D->a  = mem_array1D_D( n, "{areas de dicr.}" );
    D->cx = mem_array1D_D( n, "{cx de dicr.}" );
    D->cy = mem_array1D_D( n, "{cy de dicr.}" );
    D->nFranjas = n;
}

// Libera memoria de discretizacion
void lib_mem_discr_pol( discr_pol *D ){
    lib_array1D_D( D->a );
    lib_array1D_D( D->cx );
    lib_array1D_D( D->cy );
    D->nFranjas = 0;
}

// Inicia discretizacion
void inicia_discretizacion( discretizacion *D ){
    D->nAreas = 0;
    D->pol = NULL;
}

// Inicia memoria para discretizacion
void inicia_mem_discretizacion( discretizacion *D, const int n ){
    D->pol = (discr_pol*) malloc( n * sizeof(discr_pol) );
    if( D->pol == NULL ){
        printf( "\nError: no se pudo alojar memoria para discretizacion" );
        D->nAreas = 0;
    }
    else{
        D->nAreas = n;
        for( int i = 0; i < n; i++ ) inicia_discr_pol( &D->pol[i] );
    }
}

// Libera memoria de discretizacion
void lib_mem_discretizacion( discretizacion *D ){
    if( D->pol != NULL ){
        for( int i = 0; i < D->nAreas; i++ )
            if( D->pol[i].nFranjas > 0 ) lib_mem_discr_pol( &D->pol[i] );
        free( D->pol );
        D->pol  = NULL;
        D->nAreas = 0;
    }
}

// Inicia tabla de resultados
void inicia_tabla_result( tabla_result *T ){
    T->nDiv = 0;
    T->a  = NULL;
    T->cx = NULL;
    T->cy = NULL;
    T->def = NULL;
    T->esf = NULL;
    T->P   = NULL;
    T->Mx  = NULL;
    T->My  = NULL;
}

// Memoria para tabla de resultados
void mem_tabla_result( tabla_result *T, const int n ){
    T->a  = mem_array1D_D( n, "result. {a}" );
    T->cx = mem_array1D_D( n, "result. {cx}" );
    T->cy = mem_array1D_D( n, "result. {cy}" );
    T->def = mem_array1D_D( n, "result. {def}" );
    T->esf = mem_array1D_D( n, "result. {esf}" );
    T->P   = mem_array1D_D( n, "result. {P}" );
    T->Mx  = mem_array1D_D( n, "result. {Mx}" );
    T->My  = mem_array1D_D( n, "result. {My}" );
    T->nDiv = n;
}

// Inicia en ceros todos los valores de la tabla
void llenar_ceros( tabla_result *T ){
    for( int i = 0; i < T->nDiv; i++ ){
        T->a[i]   = 0.0;
        T->cx[i]  = 0.0;
        T->cy[i]  = 0.0;
        T->def[i] = 0.0;
        T->esf[i] = 0.0;
        T->P[i]   = 0.0;
        T->Mx[i]  = 0.0;
        T->My[i]  = 0.0;
    }
}

// Libera memoria de tabla de resultados
void lib_mem_tabla_result( tabla_result *T ){
    lib_array1D_D( T->a );
    lib_array1D_D( T->cx );
    lib_array1D_D( T->cy );
    lib_array1D_D( T->def );
    lib_array1D_D( T->esf );
    lib_array1D_D( T->P );
    lib_array1D_D( T->Mx );
    lib_array1D_D( T->My );
    T->nDiv = 0;
}

// Imprime tabla de resutlados
void imp_tabla_result( tabla_result *T, const char * const mensaje ){
    printf( "%s", mensaje );
    printf( "\nNo.   Area   Cx   Cy   Def.   Esf.   P   Mx   My" );
    for( int i = 0; i < T->nDiv; i++ ){
        printf( "\n%d   %f   %f   %f   %f   %f   %f   %f   %f", i+1,
                T->a[i], T->cx[i], T->cy[i], T->def[i], T->esf[i],
                T->P[i], T->Mx[i], T->My[i] );
    }
}

// Inicia una seccion
void inicia_seccion( seccion *S ){
    S->cx = 0.0;
    S->cy = 0.0;
    S->a  = 0.0;
    S->P  = 0.0;
    S->Mx = 0.0;
    S->My = 0.0;
    S->Pn = 0.0;
    S->Mnx = 0.0;
    S->Mny = 0.0;
    S->theta = 0.0;
    S->desp  = 0.0;
    S->ex  = 0.0;
    S->ey  = 0.0;
    S->exB = 0.0;
    S->eyB = 0.0;
    S->dh  = 0.0;
    inicia_mat_concr( &S->region_concr );
    inicia_mat_acero( &S->region_acero );
    inicia_mat_concr( &S->region_concr_desp );
    inicia_mat_acero( &S->region_acero_desp );
    inicia_discretizacion( &S->discr_concr_pos );
    inicia_discretizacion( &S->discr_concr_neg );
    inicia_discretizacion( &S->discr_acero_pos );
    inicia_discretizacion( &S->discr_acero_neg );
    inicia_discretizacion( &S->discr_acero_posT );
    inicia_discretizacion( &S->discr_acero_negT );
    inicia_graf_EsfDef( &S->aceroED );
    inicia_graf_EsfDef( &S->concrED );
}

// Libera memoria de una seccion
void lib_mem_seccion( seccion *S ){
    // Areas de acero y concreto
    lib_mem_mat_concr( &S->region_concr );
    lib_mem_mat_acero( &S->region_acero );
    lib_mem_mat_concr( &S->region_concr_desp );
    lib_mem_mat_acero( &S->region_acero_desp );
    lib_mem_graf_EsfDef( &S->aceroED );
    lib_mem_graf_EsfDef( &S->concrED );
}

// Imprime los datos de una seccion
void imp_seccion( seccion *S, const char * const mensaje ){
    printf( "%s", mensaje );
    printf( "\na  = %f", S->a );
    printf( "\ncx = %f", S->cx );
    printf( "\ncy = %f", S->cy );
    imp_mat_concr( &S->region_concr_desp, "\nRegion de concreto desplazada:", 4 );
    imp_mat_acero( &S->region_acero_desp, "\nRegion de acero desplazada:", 4 );
}

// Rota y luego desplaza la geometria para que el EN quede horizontal
// y pase por el origen (alineado con el eje de las X)
void rota_desplaza( const double theta,const double desp, concreto *cr,
                    acero *ac ){
    // Concreto
    // Areas positivas
    for( int i = 0; i < cr->parte_pos.nAreas; i++ )
        rota_despl_poligono( theta, desp, &cr->parte_pos.areas[i] );
    // Areas negativas
    for( int i = 0; i < cr->parte_neg.nAreas; i++ )
        rota_despl_poligono( theta, desp, &cr->parte_neg.areas[i] );
    // Acero
    // Areas positivas
    for( int i = 0; i < ac->parte_pos.nAreas; i++ )
        rota_despl_poligono( theta, desp, &ac->parte_pos.areas[i] );
    // Areas negativas
    for( int i = 0; i < ac->parte_neg.nAreas; i++ )
        rota_despl_poligono( theta, desp, &ac->parte_neg.areas[i] );
}

// Rota y desplaza la geometria de un material de tal forma que el EN
// quede horizontal y pase por el origen, esto en funcion de las
// variables theta y desp representan el angulo de giro y el
// desplazamiento del EN.
void rota_despl_poligono( const double theta, const double desp, poligono *geom ){
    // Rota el poligono con respecto al centroide de la seccion completa
    rota( -theta, geom->nNodos, geom->x, geom->y );
    // Desplaza el poligono de tal forma que el EN pasa por el origen
    desplaza( 0.0, -desp, geom->nNodos, geom->x, geom->y );
}

// Rota una serie de coordenadas dado un angulo de giro en 2D
void rota( const double theta, const int n, double *gx, double *gy ){
    double x, y, cos_theta, sin_theta;
    // Gira cada punto
    cos_theta = cos( theta );
    sin_theta = sin( theta );
    for( int i = 0; i < n; i++ ){
    	// {PR} = [RI]({P}-{C}), {C} = centroide de toda la seccion
    	// [RI] = [R]^(-1)
    	x = cos_theta * gx[i] - sin_theta * gy[i];
    	y = sin_theta * gx[i] + cos_theta * gy[i];
    	// Guarda el punto rotado
    	gx[i] = x;
    	gy[i] = y;
    }
}

// Desplaza una serie de coordenadas en 2D
void desplaza( const double desp_x, const double desp_y, const int n,
               double *x, double *y ){
    for( int i = 0; i < n; i++ ){
        x[i] += desp_x;
        y[i] += desp_y;
    }
}

// Pide la memoria faltante para una seccion
void mem_adicional( seccion *S ){
    // Concreto
    copia_mem_geom( &S->region_concr_desp.parte_pos,
                    &S->region_concr.parte_pos );
    copia_mem_geom( &S->region_concr_desp.parte_neg,
                    &S->region_concr.parte_neg );
    // Acero
    copia_mem_geom( &S->region_acero_desp.parte_pos,
                    &S->region_acero.parte_pos );
    copia_mem_geom( &S->region_acero_desp.parte_neg,
                    &S->region_acero.parte_neg );
}

// Aloja nueva memoria para guardar geometrias del mismo tamano
void copia_mem_geom( zona *a,  zona *b ){
    // Memoria para alojar "nAreas" poligonos
    mem_ini_zona( a, b->nAreas );
    // Memoria de los poligonos (numero de vertices)
    for( int i = 0; i < b->nAreas; i++ )
        mem_poligono( &a->areas[i], b->areas[i].nNodos );
}

// Calcula el centroide y area de la seccion
void calc_props_seccion( seccion *S ){
    double mx, my, are;
    // Calcula el area y centroide de cada poligono
    // Concreto
    calc_props_geom( &S->region_concr.parte_pos );
    calc_props_geom( &S->region_concr.parte_neg );
    // Acero
    calc_props_geom( &S->region_acero.parte_pos );
    calc_props_geom( &S->region_acero.parte_neg );
    // Calcula el centroide de la seccion
    are = 0.0;
    mx  = 0.0;
    my  = 0.0;
    // Concreto
    sum_props_geom( &S->region_concr.parte_pos, &are, &mx, &my,  1.0 );
    sum_props_geom( &S->region_concr.parte_neg, &are, &mx, &my, -1.0 );
    // Acero
    sum_props_geom( &S->region_acero.parte_pos, &are, &mx, &my,  1.0 );
    sum_props_geom( &S->region_acero.parte_neg, &are, &mx, &my, -1.0 );
    // Centroide
    S->a  = are;
    S->cx = my / are;
    S->cy = mx / are;
        /*
        printf("\nare = %f", are );
        printf("\nmx  = %f", mx );
        printf("\nmy  = %f", my );
        printf("\ncx  = %f", S->cx );
        printf("\ncy  = %f", S->cy );
        //*/

}

// Calcula los centroides y areas de una zona de la seccion
void calc_props_geom( zona *z ){
    int n;
    double *x, *y, a_i;
    for( int i = 0; i < z->nAreas; i++ ){
        n = z->areas[i].nNodos;
        x = z->areas[i].x;
        y = z->areas[i].y;
        a_i = area_pol2D( n, x, y );
        z->areas[i].cx = centroide_x_pol2D( a_i, n, x, y );
        z->areas[i].cy = centroide_y_pol2D( a_i, n, x, y );
        z->areas[i].a  = a_i;
    }
}

// Suma las areas y primeros momentos de una zona de la seccion
void sum_props_geom( zona *z, double *are, double *mx, double *my,
                     const double f ){
    for( int i = 0; i < z->nAreas; i++ ){
        *are += z->areas[i].a * f;
        *mx  += z->areas[i].a * f * z->areas[i].cy;
        *my  += z->areas[i].a * f * z->areas[i].cx;
    }
}

// Calcula el area de un poligono definido por sus coordenadas
double area_pol2D( const int n, double *x, double *y ){
    // D Definido como D = x1*y2 + x2*y3 + ... + xn*y1
    // I Definido como D = y1*x2 + y2*x3 + ... + yn*x1
    double area = 0.0, D = 0.0, I = 0.0;
    if( n >= 2 ){
        for ( int i = 0; i < n - 1; i++ ){
            D += x[i] * y[i+1];
            I += y[i] * x[i+1];
        }
        D += x[n-1] * y[0];
        I += y[n-1] * x[0];
        area = 0.5 * fabs( D - I );
    }
    return area;
}

// Calcula el centroide de un poligono definido por sus coords.
double centroide_x_pol2D( const double area, const int n,
                          double *x, double *y ){
    double x1, x2, y1, y2, c = 0.0, sum = 0.0;
    if( n <= 0 ){
        //printf( "\nCuidado: n = 0" );
        c = 0.0;
    }
    else if( n == 1 || n == 2 ){
        for( int i = 0; i < n; i++ ) sum += x[i];
        c = sum / (double)( n );
        printf("\nCuidado: se esta calculando el centroide en x de una linea o punto");
    }
    else{
        // Para cx
        for ( int i = 0; i < n - 1; i++ ){
            x1 = x[i];
            y1 = y[i];
            x2 = x[i+1];
            y2 = y[i+1];
            sum += ( y2 - y1 ) * ( x1 * x1 + x1 * x2 + x2 * x2 );
        }
        x1 = x[n-1];
        y1 = y[n-1];
        x2 = x[0];
        y2 = y[0];
        sum += ( y2 - y1 ) * ( x1 * x1 + x1 * x2 + x2 * x2 );
        if( area < 1e-10 ) c = 0.0;
        else               c = sum / ( 6.0 * area );
    }
    return c;
}
double centroide_y_pol2D( const double area, const int n,
                          double *x, double *y ){
    double x1, x2, y1, y2, c = 0.0, sum = 0.0;
    if( n <= 0 ){
        //printf( "\nCuidado: n = 0" );
        c = 0.0;
    }    
    else if( n == 1 || n == 2 ){
        for( int i = 0; i < n; i++ ) sum += y[i];
        c = sum / (double)( n );
        printf("\nCuidado: se esta calculando el centroide en y de una linea o punto");
    }
    else{
        // Para cy
        for ( int i = 0; i < n - 1; i++ ){
            x1 = x[i];
            y1 = y[i];
            x2 = x[i+1];
            y2 = y[i+1];
            sum += ( x1 - x2 ) * ( y1 * y1 + y1 * y2 + y2 * y2 );
        }
        x1 = x[n-1];
        y1 = y[n-1];
        x2 = x[0];
        y2 = y[0];
        sum += ( x1 - x2 ) * ( y1 * y1 + y1 * y2 + y2 * y2 );
        if( area < 1e-10 ) c = 0.0;
        else               c = sum / ( 6.0 * area );
    }
    return c;
}

// Discretiza la zona a compresion de la seccion
void discretiza_comp( const double dh, concreto *concr, acero *ace,
                      discretizacion *discr_concr_pos,
                      discretizacion *discr_concr_neg,
                      discretizacion *discr_acero_pos,
                      discretizacion *discr_acero_neg ){
    int existe1, existe2;
    double yini1, yini2, yini;
    // Coordenada "y" a partir de la cual se empezara a discretizar
    existe1 = yMin_zona( &concr->parte_pos, &yini1 );
    existe2 = yMin_zona( &ace->parte_pos,   &yini2 );
    if( existe1 ) yini = yini1;
    if( yini > yini2 && existe2 ) yini = yini2;
    if( yini < 0.0 ) yini = 0.0;
    if( !existe1 ) printf( "\nCuidado: no hay areas pos. de concreto!!" );
    if( !existe2 ) printf( "\nCuidado: no hay areas pos. de acero!!" );
    // Discretizacion del concreto
        //printf( "\n\n------  CONCRETO COMPRESION  ------------------\n\n" );
    discretiza_zona_comp( dh, yini, &concr->parte_pos, discr_concr_pos );
    discretiza_zona_comp( dh, yini, &concr->parte_neg, discr_concr_neg );
    // Discretizacion del acero
        //printf( "\n\n------  ACERO COMPRESION   --------------------\n\n" );
    discretiza_zona_comp( dh, yini, &ace->parte_pos, discr_acero_pos );
    discretiza_zona_comp( dh, yini, &ace->parte_neg, discr_acero_neg );
}

// Discretiza la zona a tension de la seccion
void discretiza_tension( const double dh, acero *ace,
                         discretizacion *discr_acero_pos,
                         discretizacion *discr_acero_neg ){
    int existe;
    double yini;
    // Coordenada "y" a partir de la cual se empezara a discretizar
    existe = yMax_zona( &ace->parte_pos, &yini );
    if( yini > 0.0 ) yini = 0.0;
    if( !existe ) printf( "\nCuidado: no hay areas pos. de acero!!" );
    // Discretizacion del acero
        //printf( "\n\n------  ACERO TENSION   -----------------------\n\n" );
    discretiza_zona_tension( dh, yini, &ace->parte_pos, discr_acero_pos );
    discretiza_zona_tension( dh, yini, &ace->parte_neg, discr_acero_neg );
}

/*
// Discretiza una zona en tension
void discretiza_zona_tension( const double dh, const double yini,
                              zona *Z, discretizacion *discr ){
    int nDiv;
    double ymin, y_inf, y_sup, a, cx, cy;
    // Discretizacion dee la zona
    for( int i = 0; i < Z->nAreas; i++ ){
        // Coordenada y menor de un poligono
        ymin = min_vec_D( Z->areas[i].nNodos, Z->areas[i].y );
        // Numero de divisiones
        nDiv = 0;
        if( ymin < 0.0 ) nDiv = (int)( fabs(ymin)/dh ) + 1;
        // Calcula el centroide y area de cada franja
        if( nDiv > 0 ){
            mem_discr_pol( &discr->pol[i], nDiv );
            y_inf = -dh;
            for( int k = 0; k < nDiv; k++ ){
                y_sup = y_inf + dh;
                area_centro_franja( y_inf, y_sup, &Z->areas[i], &a, &cx, &cy );
                discr->pol[i].a[k]  = a;
                discr->pol[i].cx[k] = cx;
                discr->pol[i].cy[k] = cy;
                y_inf -= dh;
            }
        }
    }
}
//*/

///*
// Discretiza una zona en tension a partir de una "y" inicial
void discretiza_zona_tension( const double dh, const double yini, zona *Z,
                              discretizacion *discr ){
    int nDiv, encontrado;
    double ymin, y_inf, y_sup, a, cx, cy;
    // NOTA: 0.0 >= yini >= ymin
    // Discretizacion dee la zona
    for( int i = 0; i < Z->nAreas; i++ ){
        // Coordenada y menor de un poligono
        encontrado = min_vec_D( Z->areas[i].nNodos, Z->areas[i].y, &ymin );
        // Numero de divisiones
        nDiv = 0;
        if( ymin < 0.0 && encontrado ) nDiv = (int)( fabs(ymin-yini)/dh ) + 1;
        // Calcula el centroide y area de cada franja
        mem_discr_pol( &discr->pol[i], nDiv );
        y_inf = yini - dh;
        for( int k = 0; k < nDiv; k++ ){
            y_sup = y_inf + dh;
            area_centro_franja( y_inf, y_sup, &Z->areas[i], &a, &cx, &cy );
            discr->pol[i].a[k]  = a;
            discr->pol[i].cx[k] = cx;
            discr->pol[i].cy[k] = cy;
            y_inf -= dh;
        }
    }
}
//*/


/*
// Discretiza una zona en compresion
void discretiza_zona_comp( const double dh, const double yini,
                           zona *Z, discretizacion *discr ){
    int nDiv;
    double ymax, y_inf, y_sup, a, cx, cy;
    // Discretizacion dee la zona
    for( int i = 0; i < Z->nAreas; i++ ){
        // Coordenada y mayor de un poligono
        //ymax = max_vec_D( Z->areas[i].nNodos, Z->areas[i].y );
        max_vec_D( Z->areas[i].nNodos, Z->areas[i].y, &ymax );
        // Numero de divisiones
        nDiv = 0;
        if( ymax > 0.0 ) nDiv = (int)( ymax/dh ) + 1;
        // Calcula el centroide y area de cada franja
        if( nDiv > 0 ){
            mem_discr_pol( &discr->pol[i], nDiv );
            y_inf = 0.0;
            for( int k = 0; k < nDiv; k++ ){
                y_sup = y_inf + dh;
                area_centro_franja( y_inf, y_sup, &Z->areas[i], &a, &cx, &cy );
                discr->pol[i].a[k]  = a;
                discr->pol[i].cx[k] = cx;
                discr->pol[i].cy[k] = cy;
                y_inf = y_sup;
            }
        }
    }
}
//*/

///* LA BUENA
// Discretiza una zona en compresion a partir de una "y" inicial
void discretiza_zona_comp( const double dh, const double yini, zona *Z,
                           discretizacion *discr ){
    int nDiv, encontrado;
    double ymax, y_inf, y_sup, a, cx, cy;
    // NOTA: 0.0 <= yini <= ymax
    // Discretizacion de la zona
    for( int i = 0; i < Z->nAreas; i++ ){
        // Coordenada y mayor de un poligono
        encontrado = max_vec_D( Z->areas[i].nNodos, Z->areas[i].y, &ymax );
            //printf("\nyini = %f  ymax = %f", yini, ymax );
            //getchar();
        // Numero de divisiones
        nDiv = 0;
        if( ymax > 0.0 && encontrado ) nDiv = (int)( (ymax-yini)/dh ) + 1;
        // Calcula el centroide y area de cada franja
        mem_discr_pol( &discr->pol[i], nDiv );
        y_inf = yini;
        for( int k = 0; k < nDiv; k++ ){
            y_sup = y_inf + dh;
            area_centro_franja( y_inf, y_sup, &Z->areas[i], &a, &cx, &cy );
            discr->pol[i].a[k]  = a;
            discr->pol[i].cx[k] = cx;
            discr->pol[i].cy[k] = cy;
            y_inf = y_sup;
        }
    }
}
//*/

// Determina el area y centroide de una franja sobre un poligono
void area_centro_franja( const double y_inf, const double y_sup,
                         poligono *P, double *a, double *cx, double *cy ){
    miDeque_D xF, yF;
    inicia_miDeque_D( &xF );
    inicia_miDeque_D( &yF );
    // Determina la conectividad de la franja
    conectividad_franja( y_inf, y_sup, P, &xF, &yF );
    // Calcula el area y los centroides
    *a  = area_pol2D( xF.n, xF.val, yF.val );
    *cx = centroide_x_pol2D( *a, xF.n, xF.val, yF.val );
    *cy = centroide_y_pol2D( *a, xF.n, xF.val, yF.val );
    // Libera memoria
    lib_mem_miDeque_D( &xF );
    lib_mem_miDeque_D( &yF );
}

// Determina los nodos que definen una franja (conectividad)
void conectividad_franja( const double y_inf, const double y_sup, 
                          poligono *P, miDeque_D *xF, miDeque_D *yF ){
    int contI, contS;
    double dist_inf, dist_sup;
    // Variables para guardar las intersecciones
    miDeque_D xI, yI, xS, yS;
    miDeque_I lI, lS;
    inicia_miDeque_D( &xI ); inicia_miDeque_D( &yI ); inicia_miDeque_I( &lI );
    inicia_miDeque_D( &xS ); inicia_miDeque_D( &yS ); inicia_miDeque_I( &lS );
    // Calcula las intersecciones del poligono P con la horizontal y_inf
    intersecciones( y_inf, P->nNodos, P->x, P->y, &xI, &yI, &lI );
    // Calcula las intersecciones del poligono P con la horizontal y_sup
    intersecciones( y_sup, P->nNodos, P->x, P->y, &xS, &yS, &lS );
        //printf( "\nxI.n_mem = %d", xI.n_mem );
        //printf( "\nxS.n_mem = %d", xS.n_mem );
        //getchar( );
    // Determina la conectividad del poligono que define la franja
    contI = 0;
    contS = 0;
    for( int i = 0; i < P->nNodos; i++ ){
        // Agrega el vertice en caso que cumpla
        if( P->y[i] >= y_inf && P->y[i] <= y_sup ){
            agrega_inter( P->x[i], P->y[i], xF, yF );
        }
        // Agrega las intersecciones, en caso de que haya
        if( contI < xI.n && contS < xS.n ){
            if( lI.val[contI] == i && lS.val[contS] == i  ){
                dist_inf = dist( P->x[i], P->y[i], xI.val[contI], yI.val[contI] );
                dist_sup = dist( P->x[i], P->y[i], xS.val[contS], yS.val[contS] );
                if( dist_inf < dist_sup ){
                    agrega_inter( xI.val[contI], yI.val[contI], xF, yF );
                    agrega_inter( xS.val[contS], yS.val[contS], xF, yF );
                }
                else{
                    agrega_inter( xS.val[contS], yS.val[contS], xF, yF);
                    agrega_inter( xI.val[contI], yI.val[contI], xF, yF );
                }
                contI++;
                contS++;
            }
            else if( lI.val[contI] == i ){
                agrega_inter( xI.val[contI], yI.val[contI], xF, yF );
                contI++;
            }
            else if( lS.val[contS] == i ){
                agrega_inter( xS.val[contS], yS.val[contS], xF, yF );
                contS++;
            }
            else{ }            
        }
        else{
            if( contI < xI.n ){
                if( lI.val[contI] == i ){
                    agrega_inter( xI.val[contI], yI.val[contI], xF, yF );
                    contI++;
                }            
            }
            else if( contS < xS.n ){
                if( lS.val[contS] == i ){
                    agrega_inter( xS.val[contS], yS.val[contS], xF, yF );
                    contS++;
                }            
            }
            else{ }
        }
    }
        //imp_miDeque_D( xF, "\nCoord. x de franja:" );
        //imp_miDeque_D( yF, "\nCoord. y de franja:" );
        /*
        if( y_inf >= 0.0 )
            imp_conectividad( xF, yF, "\nConectividad de franja COMPRESION:" );
        else
            imp_conectividad( xF, yF, "\nConectividad de franja TENSION:" );
        printf("\n");
        */
        //getchar( );
    // Libera memoria
    lib_mem_miDeque_I( &lI );
    lib_mem_miDeque_D( &xI );
    lib_mem_miDeque_D( &yI );
    lib_mem_miDeque_I( &lS );
    lib_mem_miDeque_D( &xS );
    lib_mem_miDeque_D( &yS );
}

//.........................................................................

// Agrega interseccion
void agrega_inter( double x, double y, miDeque_D *xF, miDeque_D *yF ){
    add_miDeque_D( xF, x );
    add_miDeque_D( yF, y );
}

// Verifica si un lado y sus intersecciones estan dentro de la franja
void agrega_lado_inter( const double y_inf, const double y_sup, const int i,
                        double xp, double yp, double xi, double yi, int li,
                        double xs, double ys, int ls, int *contI, int *contS,
                        miDeque_D *xF, miDeque_D *yF ){
    double dist_inf, dist_sup;
    // Agrega el vertice en caso que cumpla
    if( yp >= y_inf && yp <= y_sup ){
        add_miDeque_D( xF, xp );
        add_miDeque_D( yF, yp );
    }
    // Agrega las intersecciones, en caso de que haya
    if( li == i && ls == i  ){
        dist_inf = dist( xp, yp, xi, yi );
        dist_sup = dist( xp, yp, xs, ys );
        if( dist_inf < dist_sup ){
            add_miDeque_D( xF, xi );
            add_miDeque_D( yF, yi );
            add_miDeque_D( xF, xs );
            add_miDeque_D( yF, ys );
        }
        else{
            add_miDeque_D( xF, xs );
            add_miDeque_D( yF, ys );
            add_miDeque_D( xF, xi );
            add_miDeque_D( yF, yi );
        }
        (*contI)++;
        (*contS)++;
    }
    else if( li == i ){
        add_miDeque_D( xF, xi );
        add_miDeque_D( yF, yi );
        (*contI)++;
    }
    else if( ls == i ){
        add_miDeque_D( xF, xs );
        add_miDeque_D( yF, ys );
        (*contS)++;
    }
    else{ }
}

// Calcula la distancia euclideana entre dos puntos en 2D
double dist( double P1x, double P1y, double P2x, double P2y ){
    return (P2x-P1x)*(P2x-P1x) + (P2y-P1y)*(P2y-P1y);
}

// Calcula las intersecciones entre un poligono y una recta horizonal
void intersecciones( const double a, const int np, double *px, double *py,
                     miDeque_D *inter_x, miDeque_D *inter_y, miDeque_I *inter_l ){
    double x, y;
    int inter;
    for( int i = 0; i < np - 1; i++ ){
        // Calcula el punto de interseccion entre la recta que
        // contiene al lado y la horizontal
        intersecta( px[i], py[i], px[i+1], py[i+1], 0.0, a, 1.0, a,
                    &inter, &x, &y );
        // Verifica si el punto de interseccion esta contenido
        // en el lado del poligono
        if( inter ){
            if( punto_segmento( px[i], py[i], px[i+1], py[i+1], x, y ) ){
                add_miDeque_D( inter_x, x );
                add_miDeque_D( inter_y, y );
                add_miDeque_I( inter_l, i );
            }
        }
    }
    // Ultimo lado
    intersecta( px[np-1], py[np-1], px[0], py[0], 0.0, a, 1.0, a,
                &inter, &x, &y );
    if( inter ){
        if( punto_segmento( px[np-1], py[np-1], px[0], py[0], x, y ) ){
            add_miDeque_D( inter_x, x );
            add_miDeque_D( inter_y, y );
            add_miDeque_I( inter_l, np-1 );
        }
    }
}

// Determina el punto de interseccion entre dos rectas definidas
// por dos puntos cada una
void intersecta( const double L1x, const double L1y,
                 const double L2x, const double L2y,
                 const double P1x, const double P1y,
                 const double P2x, const double P2y,
                 int *inter, double *x, double *y ){
    int horL, vertL, horP, vertP; // Para saber si son verticales
                                  // u horizontales
    double mL, mP, bL, bP; // Pendientes y constantes
    mL = pendiente_recta( L1x, L1y, L2x, L2y, &horL, &vertL );
    mP = pendiente_recta( P1x, P1y, P2x, P2y, &horP, &vertP );
    // Calcula el punto de interseccion
    *inter = 1;
    // Cuando ambas rectas son verticales
    if( vertL == 1 && vertP == 1 ){
        *inter = 0;
    }
    // Cuando ambas rectas son horizontales
    else if( horL == 1 && horP == 1 ){
        *inter = 0;
    }
    // La primera vertical
    else if( vertL == 1 && vertP == 0 ){
        // Constanes
        bP = P1y - mP * P1x;
        // Punto de interseccion
        *x = L1x;
        *y = mP * L1x + bP;
    }
    // La segunda vertical
    else if( vertL == 0 && vertP == 1 ){
        // Constantes
        bL = L1y - mL * L1x;
        // Punto de intersección
        *x = P1x;
        *y = mL * P1x + bL;
    }
    // No hay verticales pero si puede haber una horizontal
    else{
        bP = P1y - mP * P1x;
        bL = L1y - mL * L1x;
        *x = ( bP - bL ) / ( mL - mP );
        *y = mL * (*x) + bL;
    }
}

// Determina si un punto esta contenido en un segmentos de recta
int punto_segmento( const double L1x, const double L1y,
                    const double L2x, const double L2y,
                   double x, double y ){
    int horL, vertL, valida;
    double lim1_x, lim2_x, lim1_y, lim2_y;
    // Identifica si el segmento es vertical, horizonal o inclinado
    pendiente_recta( L1x, L1y, L2x, L2y, &horL, &vertL );
    // Ordena los puntos de menor a mayor
    lim1_x = L1x;
    lim2_x = L2x;
    if( L2x < L1x ){
        lim1_x = L2x;
        lim2_x = L1x;
    }
    lim1_y = L1y;
    lim2_y = L2y;
    if( L2y < L1y ){
        lim1_y = L2y;
        lim2_y = L1y;
    }
    valida = 0;
    // Cuando no es vertical se emplea el eje "x"
    if( vertL == 0 ){
        if( x >= lim1_x && x <= lim2_x ) valida = 1;
    }
    // Cuando es vertical se emplea el eje "y"
    else{
        if( y >= lim1_y && y <= lim2_y ) valida = 1;
    }
    return valida;
}

// Calcula la pendiente de la recta y determina si la recta es
// una horizontal o una vertical
double pendiente_recta( const double x1, const double y1,
                        const double x2, const double y2,
                        int *hor, int *vert ){
    double m;  // Pendiente
    double pi = 3.14159265358979;
    double t  = 1e-10;   // Factor para determinar si es hor. o vert.
    double lx = x2 - x1; // Longitud en "x"
    double ly = y2 - y1; // Longitud en "y"
    double theta = atan( ly / lx );
    // Cuando se trata de una vertical
    if( fabs( 0.5*pi - fabs(theta) ) < t ){
        *vert = 1;
        *hor  = 0;
        m    = 0.0;
    }
    // Cuando se trata de una horizotal
    else if( fabs(theta) < t ){
        *vert = 0;
        *hor  = 1;
        m    = 0.0;
    }
    else{
        *vert = 0;
        *hor  = 0;
        m = ly / lx;
    }
    return m;
}

// Inicia valores de la seccion
void inicia_val( seccion *S ){
    mem_adicional( S );
    // Calcula propiedades de la seccion
    calc_props_seccion( S );
    ///*
    // Hace coindidir el origen con el centroide de la seccion
    // Concreto
    for( int i = 0; i < S->region_concr.parte_pos.nAreas; i++ ){
        desplaza( -S->cx, -S->cy, S->region_concr.parte_pos.areas[i].nNodos,
                   S->region_concr.parte_pos.areas[i].x,
                   S->region_concr.parte_pos.areas[i].y );
    }
    for( int i = 0; i < S->region_concr.parte_neg.nAreas; i++ ){
        desplaza( -S->cx, -S->cy, S->region_concr.parte_neg.areas[i].nNodos,
                   S->region_concr.parte_neg.areas[i].x,
                   S->region_concr.parte_neg.areas[i].y );
    }
    // Acero
    for( int i = 0; i < S->region_acero.parte_pos.nAreas; i++ ){
        desplaza( -S->cx, -S->cy, S->region_acero.parte_pos.areas[i].nNodos,
                   S->region_acero.parte_pos.areas[i].x,
                   S->region_acero.parte_pos.areas[i].y );
    }
    for( int i = 0; i < S->region_acero.parte_neg.nAreas; i++ ){
        desplaza( -S->cx, -S->cy, S->region_acero.parte_neg.areas[i].nNodos,
                   S->region_acero.parte_neg.areas[i].x,
                   S->region_acero.parte_neg.areas[i].y );
    }
    //*/
}

// Evalua la resistencia de una seccion de concreto y acero en funcion
// de la posicion del EN( theta, desp )
void eval_resist_seccion( const double theta, const double desp, seccion *S ){
    S->theta = theta;
    S->desp  = desp;
    inicia_mem_discretizacion( &S->discr_concr_pos, S->region_concr.parte_pos.nAreas );
    inicia_mem_discretizacion( &S->discr_concr_neg, S->region_concr.parte_neg.nAreas );
    inicia_mem_discretizacion( &S->discr_acero_pos, S->region_acero.parte_pos.nAreas );
    inicia_mem_discretizacion( &S->discr_acero_neg, S->region_acero.parte_neg.nAreas );
    inicia_mem_discretizacion( &S->discr_acero_posT, S->region_acero.parte_pos.nAreas );
    inicia_mem_discretizacion( &S->discr_acero_negT, S->region_acero.parte_neg.nAreas );
    // Hace una copia de los materiales y geometrias
    copia_mat_concr( &S->region_concr_desp, &S->region_concr );
    copia_mat_acero( &S->region_acero_desp, &S->region_acero );
    // Rota y mueve la secc. de forma que el EN se alinea con el eje X
    rota_desplaza( theta, desp, &S->region_concr_desp, &S->region_acero_desp );
    // Discretiza la zona a compresion
    discretiza_comp( S->dh, &S->region_concr_desp, &S->region_acero_desp,
                     &S->discr_concr_pos, &S->discr_concr_neg,
                     &S->discr_acero_pos, &S->discr_acero_neg );
    // Discretiza zona a tension
    discretiza_tension( S->dh, &S->region_acero_desp, &S->discr_acero_posT,
                        &S->discr_acero_negT );
    // Calcula fuerza y momento resistentes
    resist_nominal( desp, 0.003, 50.0*0.003,
                    &S->region_concr_desp, &S->region_acero_desp, 
                    &S->discr_concr_pos, &S->discr_concr_neg,
                    &S->discr_acero_pos, &S->discr_acero_neg,
                    &S->discr_acero_posT, &S->discr_acero_negT,
                    &S->concrED, &S->aceroED, &S->Pn, &S->Mnx,
                    &S->Mny );
    // Devuelve a los ejes originales los momentos
    ejes_ini( S->theta, &S->Mnx, &S->Mny );
    // Calcula excentricidades
    calcExcent( &S->ex, &S->ey, S->Pn, S->Mnx, S->Mny );
        /*
        printf( "\nResistencia nominal de la seccion (ejes originales)" );
        printf( "\n\tPn  = %f", S->Pn );
        printf( "\n\tMnx = %f", S->Mnx );
        printf( "\n\tMny = %f", S->Mny );
        printf( "\n\texB = %f", S->exB );
        printf( "\n\teyB = %f", S->eyB );
        printf( "\n\tex  = %f", S->ex );
        printf( "\n\tey  = %f", S->ey );
        printf( "\n" );
        //*/    
    // Libera memoria
    lib_mem_discretizacion( &S->discr_concr_pos );
    lib_mem_discretizacion( &S->discr_concr_neg );
    lib_mem_discretizacion( &S->discr_acero_pos );
    lib_mem_discretizacion( &S->discr_acero_neg );
    lib_mem_discretizacion( &S->discr_acero_posT );
    lib_mem_discretizacion( &S->discr_acero_negT );
}

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
                     double *Pn, double *Mnx, double *Mny ){
    int nFC, *F, nFT, aux1, aux2;
    double c, a;
    tabla_result TconcrComp, TaceroComp, TaceroTens;
    inicia_tabla_result( &TconcrComp );
    inicia_tabla_result( &TaceroComp );
    inicia_tabla_result( &TaceroTens );
    // Determina el mayor numero de franjas
    // Compresion
        F = mem_array1D_I( 4, "{nF}" );
        F[0] = max_nFranjas( discr_concr_pos );
        F[1] = max_nFranjas( discr_concr_neg );
        F[2] = max_nFranjas( discr_acero_pos );
        F[3] = max_nFranjas( discr_acero_neg );
        max_vec_I( 4, F, &nFC );
            //printf( "\nnFC = %d", nFC );
        // Crea una tabla para acumular los valores
        mem_tabla_result( &TconcrComp, nFC );
        mem_tabla_result( &TaceroComp, nFC );
        // Llena las tablas con ceros
        llenar_ceros( &TconcrComp );
        llenar_ceros( &TaceroComp );
        // Acumula las cantidades correspondientes
        llena_tabla( discr_concr_pos, &TconcrComp,  1.0 );
        llena_tabla( discr_concr_neg, &TconcrComp, -1.0 );
        llena_tabla( discr_acero_pos, &TaceroComp,  1.0 );
        llena_tabla( discr_acero_neg, &TaceroComp, -1.0 );
        // Calcula los centroides
        for( int i = 0; i < nFC; i++ ){
            if( TconcrComp.a[i] > 1e-15 ){
                TconcrComp.cx[i] /= TconcrComp.a[i];
                TconcrComp.cy[i] /= TconcrComp.a[i];
            }
            else{
                TconcrComp.cx[i] = 0.0;
                TconcrComp.cy[i] = 0.0;                
            }         
            if( TaceroComp.a[i] > 1e-15 ){
                TaceroComp.cx[i] /= TaceroComp.a[i];
                TaceroComp.cy[i] /= TaceroComp.a[i];
            }
            else{
                TaceroComp.cx[i] = 0.0;
                TaceroComp.cy[i] = 0.0;
            }
        }
        // Calcula las deformaciones, esfuerzos, fuerzas y momentos
        yMax_zona( &concr_desp->parte_pos, &c );
        yMin_zona( &acero_desp->parte_pos, &a );
            //printf( "\nc = %f", c );
        for( int i = 0; i < nFC; i++ ){
            // Deformaciones
            TconcrComp.def[i] = deformacion( TconcrComp.cy[i], c, a, defMaxConcr, defMaxAcero );
            TaceroComp.def[i] = deformacion( TaceroComp.cy[i], c, a, defMaxConcr, defMaxAcero );
            // Esfuerzos (compresion(-) y tension(+))
            TconcrComp.esf[i] = -esfuerzo( TconcrComp.def[i], concrED );
            TaceroComp.esf[i] = -esfuerzo( TaceroComp.def[i], aceroED );
            // Fuerzas y momentos
            TconcrComp.P[i]  =  TconcrComp.esf[i] * TconcrComp.a[i];
            TconcrComp.Mx[i] =  TconcrComp.P[i] * ( TconcrComp.cy[i] + desp );
            TconcrComp.My[i] = -TconcrComp.P[i] * TconcrComp.cx[i];

            TaceroComp.P[i]  =  TaceroComp.esf[i] * TaceroComp.a[i];
            TaceroComp.Mx[i] =  TaceroComp.P[i] * ( TaceroComp.cy[i] + desp );
            TaceroComp.My[i] = -TaceroComp.P[i] * TaceroComp.cx[i];
        }
        //imp_tabla_result( &TconcrComp, "\n\nCompresion: concreto" );
        //imp_tabla_result( &TaceroComp, "\n\nCompresion: acero" );
    // Tension
        aux1 = max_nFranjas( discr_acero_posT );
        aux2 = max_nFranjas( discr_acero_negT );
        nFT  = aux1;
        if( nFT < aux2 ) nFT = aux2;
            //printf( "\nnFT = %d", nFT );
        mem_tabla_result( &TaceroTens, nFT );
        // Inicia con ceros
        llenar_ceros( &TaceroTens );
        // Acumula las cantidades correspondientes
        llena_tabla( discr_acero_posT, &TaceroTens,  1.0 );
        llena_tabla( discr_acero_negT, &TaceroTens, -1.0 );
        // Calcula los centroides
        for( int i = 0; i < nFT; i++ ){       
            if( TaceroTens.a[i] > 1e-15 ){
                TaceroTens.cx[i] /= TaceroTens.a[i];
                TaceroTens.cy[i] /= TaceroTens.a[i];
            }
            else{
                TaceroTens.cx[i] = 0.0;
                TaceroTens.cy[i] = 0.0;
            }
        }
        // Calcula las deformaciones, esfuerzos, fuerzas y momentos
        for( int i = 0; i < nFT; i++ ){
            // Deformaciones
            TaceroTens.def[i] = deformacion( TaceroTens.cy[i], c, a, defMaxConcr, defMaxAcero );
            // Esfuerzos (compresion(-) y tension(+))
            TaceroTens.esf[i] = esfuerzo( TaceroTens.def[i], aceroED );
            // Fuerzas y momentos
            TaceroTens.P[i]  =  TaceroTens.esf[i] * TaceroTens.a[i];
            TaceroTens.Mx[i] =  TaceroTens.P[i] * ( TaceroTens.cy[i] + desp );
            TaceroTens.My[i] = -TaceroTens.P[i] * TaceroTens.cx[i];
        }
        //imp_tabla_result( &TaceroTens, "\n\nTension: acero" );
    // Resistencia nominal
        *Pn = *Mnx = *Mny = 0.0;
        for( int i = 0; i < nFC; i++ ){
            *Pn  += TconcrComp.P[i];
            *Pn  += TaceroComp.P[i];
            *Mnx += TconcrComp.Mx[i];
            *Mnx += TaceroComp.Mx[i];
            *Mny += TconcrComp.My[i];
            *Mny += TaceroComp.My[i];            
        }
        for( int i = 0; i < nFT; i++ ){
            *Pn  += TaceroTens.P[i];
            *Mnx += TaceroTens.Mx[i];
            *Mny += TaceroTens.My[i];
        }
        /*
        printf( "\nResistencia nominal de la seccion (ejes rotados)" );
        printf( "\n\tPn  = %f", *Pn );
        printf( "\n\tMnx = %f", *Mnx );
        printf( "\n\tMny = %f", *Mny );
        printf( "\n" );
        //*/
    // Libera memoria
    lib_array1D_I( F );
    lib_mem_tabla_result( &TconcrComp );
    lib_mem_tabla_result( &TaceroComp );
    lib_mem_tabla_result( &TaceroTens );
}

int max_nFranjas( discretizacion *discr ){
    int max = 0;
    for( int i = 0; i < discr->nAreas; i++ )
        if( discr->pol[i].nFranjas > max ) max = discr->pol[i].nFranjas;
    return max;
}

void llena_tabla( discretizacion *discr, tabla_result *T, const double f ){
    for( int k = 0; k < discr->nAreas; k++ ){
        for( int i = 0; i < discr->pol[k].nFranjas; i++ ){
            // Areas
            T->a[i] += discr->pol[k].a[i] * f;
            // Para los centroides
            T->cx[i] += discr->pol[k].cx[i] * discr->pol[k].a[i] * f;
            T->cy[i] += discr->pol[k].cy[i] * discr->pol[k].a[i] * f;
        }
    }
}

// Calcula el esfuerzo en base a las graficas de esfuerzo-defor-
// macion dada un deformacion unitaria
double esfuerzo( const double def, graf_EsfDef *G ){
    int hor, vert, pos;
    double esf, x1, x2, y1, y2, m;
    // Cuando no se ha definido correctamente la grafica
    if( G->n <= 1 ){
        esf = 0.0;
    }
    // Cuando la deformacion requerida no esta dentro de la grafica
    else if( def < G->def[0] ){
        esf = 0.0;
    }
    else if( def > G->def[G->n-1] ){
        esf = G->esf[G->n-1];
    }
    // Cuando esta dentro de la grafica
    else{
        pos = 1;
        for( int i = 1; i < G->n; i++ ){
            if( G->def[i-1] <= def && G->def[i] >= def ){
                pos = i;
                break;
            }
        }
        x1 = G->def[pos-1];
        y1 = G->esf[pos-1];
        x2 = G->def[pos];
        y2 = G->esf[pos];
        m = pendiente_recta( x1, y1, x2, y2, &hor, &vert );
        if( vert == 1 ) esf = y1;
        else            esf = y1 + m * ( def - x1 );
        /*
        printf( "\nx1 = %f", x1 );
        printf( "\nx2 = %f", x2 );
        printf( "\ny1 = %f", y1 );
        printf( "\ny2 = %f", y2 );
        printf( "\nm  = %f", m );
        printf( "\nhor  = %d", hor );
        printf( "\nvert = %d", vert );
        printf( "\nG->def[0] = %f", G->def[0] );
        printf( "\ndef = %f", def );
        printf( "\npos = %d", pos );
        printf( "\nesf = %f", esf );
        getchar();
        //*/
    }
    return esf;
}

// Calcula las excentricidades de acuerdo a los momentos y la fuerza axial
void calcExcent( double *ex, double *ey, double P, double Mx, double My ){
     double tol = 0.1;
    // Compresión
    if ( P < -tol ){
        if ( Mx >= 0.0 ){
            *ey = -fabs( Mx / P );
        }
        else{
            *ey =  fabs( Mx / P );
        }

        if ( My >= 0.0 ){
            *ex = -fabs( My / P );
        }
        else{
            *ex =  fabs( My / P );
        }
    }
    // Tensión
    else if ( P > tol ){
        if ( Mx >= 0.0 ){
            *ey =  fabs( Mx / P );
        }
        else{
            *ey = -fabs( Mx / P );
        }
        if ( My >= 0.0 ){
            *ex =  fabs( My / P );
        }
        else{
            *ex = -fabs( My / P );
        }
    }
    // Cuando no hay fuerza axial
    else{
        *ex = 0.0;
        *ey = 0.0;
    }
}

// Devuelve a los ejes originales los momentos
void ejes_ini( const double theta, double *Mnx, double *Mny ){
    double aux1, aux2;
    aux1 = cos(theta)*(*Mnx) - sin(theta)*(*Mny);
    aux2 = sin(theta)*(*Mnx) + cos(theta)*(*Mny);
    *Mnx = aux1;
    *Mny = aux2;
}

void imp_conectividad( miDeque_D *x, miDeque_D *y, const char * const mensaje ){
    printf( "%s", mensaje );
    for( int i = 0; i < x->n; i++ ){
        printf( "\n%f   %f", x->val[i], y->val[i] );
    }
}

// Fija el valor del tamano de la discretizacion
void fija_dh( const double dh, seccion *S ){
    S->dh = dh;
}

// Fija acciones mecanicas externas y excebtricidades buscadas
void fija_cargas( const double P, const double Mx, const double My, seccion *S ){
    S->P  = P;
    S->Mx = Mx;
    S->My = My;
    // Calcula excentricidades que se desean encontrar
    calcExcent( &S->exB, &S->eyB, S->P, S->Mx, S->My );
}

// Calcula la deformacion a una altura dada
double deformacion( const double y, const double c, const double a,
                    const double defMaxConcr, const double defMaxAcero ){
    // c  Ubicacion de la fibra superior de concreto
    // a  Ubicacion de la fibra inferior de acero
    // Cuando el EN esta por debajo de la fibra superior de concreto
    double def = 0.0, delta;
    ///*
    delta = defMaxConcr * fabs( c + a ) / ( defMaxConcr + defMaxAcero );
    if( c > delta ){
        def = fabs(y) * defMaxConcr /  c;
    }
    else{
        def = fabs(y) * defMaxAcero / fabs(a);
    }
    //*/
    //def = fabs(y) * defMaxConcr /  c;
    return def;
}

// Genera la superficie de intaraccion de una seccion
void sup_interaccion( const double deltaT, const double deltaD,
                      const double desp_inf, const double desp_sup,
                       seccion *S, const char * const archivo ){
    double theta, desp;
    FILE *res;
    res = fopen( archivo, "w" );
    if( res == NULL ){
        printf( "\nError: no se pudo crear o abrir el archivo %s", archivo );
        exit( EXIT_FAILURE );
    }
    theta = 0.0;
    while( theta <= 2.0*pi ){
        desp = desp_inf;
        while( desp <= desp_sup ){
            eval_resist_seccion( theta, desp, S );
            fprintf( res, "%lf\t%lf\t%lf\n", S->Mnx, S->Mny, -S->Pn );
            desp += deltaD;
        }
        fprintf( res, "\n" );
        theta += deltaT;
    }
    fclose(res);
}