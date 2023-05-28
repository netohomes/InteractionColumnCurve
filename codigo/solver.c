#include "solver.h"
//=========================================================================
// METODOS DE SOLUCION PARA ENCONTRAR LAS EXCENTRICIDADES
//.........................................................................

void solver( seccion *ss  ){
    // Encuentra la posición del eje neutro para las excentricidades
    // buscadas
    
    //miSolver( ss );

    int n, e1, e2;
    double delta, alpha, gamma, beta, tau, theta_ini, desp_ini, delta_desp;
    double (*f)( double, double, seccion* ), *x, a1, a2, aux;
    x = mem_array1D_D( 2, "{x}" );
    f = &error_tipo1;

    // Parametros para el nelder-mead
    n = 2;       // Numero de variables
    delta = 1.0; // Tamano de paso
    alpha = 1.0; // Parametro de reflexion
    gamma = 1.0; // Parametro de expansion
    beta  = 0.5; // Parametro de contraccion
    tau   = 0.5; // Parametro de contraccion
    theta_ini = 0.0;                     // Angulo inicial
    desp_ini  = despFlexPura( 0.0, ss ); // Desplazamiento inicial
    // Un pequeno desplazamiento
    // Punto inicial para compresion
    if( ss->P < 0.0 ){
        e1 = yMin_zona( &ss->region_concr.parte_pos, &a1 );
        e2 = yMin_zona( &ss->region_acero.parte_pos, &a2 );
        if( e1 == 1 && e2 == 1 ){
            aux = a1;
            if( a1 > a2 ) aux = a2;
        }
        else if( e1 == 1 && e2 == 0 ){
            aux = a1;
        }
        else{
            aux = a2;
        }
        delta_desp = 2.0 * fabs( fabs(aux) - fabs(desp_ini) );
            printf( "\ndesp_FP = %f", desp_ini );
            printf( "\naux     = %f", aux );
            printf( "\ndelta_desp = %f", delta_desp );
            getchar();
        x[0] = theta_ini;
        x[1] = desp_ini - delta_desp;
    }
    // Punto inicial para tension
    else{
        e1 = yMax_zona( &ss->region_concr.parte_pos, &a1 );
        e2 = yMax_zona( &ss->region_acero.parte_pos, &a2 );
        if( e1 == 1 && e2 == 1 ){
            aux = a1;
            if( a1 < a2 ) aux = a2;
        }
        else if( e1 == 1 && e2 == 0 ){
            aux = a1;
        }
        else{
            aux = a2;
        }
        delta_desp = 0.25 * fabs( fabs(aux) - fabs(desp_ini) );
            printf( "\ndesp_FP = %f", desp_ini );
            printf( "\naux     = %f", aux );
            printf( "\ndelta_desp = %f", delta_desp );
            getchar();  
        x[0] = theta_ini;
        x[1] = desp_ini + delta_desp;
        //x[1] = 0.0;
    }

    nelder_mead( n, x, delta, alpha, gamma, beta, tau, ss, f );

    // Imprime datos de solucion
    printf( "\nResistencia nominal de la seccion:" );
    printf( "\n\tPn  = %f", ss->Pn );
    printf( "\n\tMnx = %f", ss->Mnx );
    printf( "\n\tMny = %f", ss->Mny );
    printf( "\n\texB = %f", ss->exB );
    printf( "\n\teyB = %f", ss->eyB );
    printf( "\n\tex  = %f", ss->ex );
    printf( "\n\tey  = %f", ss->ey );
    printf( "\n\ttheta = %f", x[0] );
    printf( "\n\tdesp  = %f", x[1] );
    printf( "\n" );

    // Libera memoria
    lib_array1D_D( x );
}

void miSolver( seccion *ss ){
    // Encamina la búsqueda de acuerdo a las excentricidades que se
    // buscan
    int excBase;
    double porcentaje;
    porcentaje = fabs( ss->exB / ss->eyB );
    // Para escoger una excentricidad dominante
    if( porcentaje > 50.0 ) excBase = 0;
    else                    excBase = 1;
    busca_ex_ey( ss, excBase );
}

void busca_ex_ey( seccion *ss, const int excent ){
    // Para cada giro de calcula el "desp" para obtener la excentri-
    // cidad de acuerdo a la excentricidad dada en los argumentos
    double theta, theta1, theta2, dt, dtIni, tol, err1, err2, desp;
    double errMax, thetaIni;
    int sol, vert, hor, aux = 1;
    double (*busca )( double, seccion*, int* );
    double (*func  )( double, double, seccion* );
    int conta;
    double pi = 3.1415926535897932384;
    // Apuntadores de las funciones que se emplearán
    if( excent == 0 ){
        busca = &busca_ex;
        func  = &funcEy;
    }
    else{
        busca = &busca_ey;
        func  = &funcEx;
    }
        //printf( "\nexcent = %d", excent );
        //getchar( );
    // Variables para el proceso de iteración
    thetaIni = 0.0;
    dtIni = dt = 0.1;  // Incremento en el ángulo
    tol = 0.00001;     // Error mínimo aceptado
    errMax = 1000000.0;// Error máximo en iteraciones
    theta = thetaIni;       // Ángulo inicial
    desp = busca( theta, ss, &sol ); // Desp. inicial
    err1 = func( theta, desp, ss );  // Error inicial
    theta1 = theta;
    conta  = 0;
    printf( "\nGiro  Desp  Err1  Err2  Incr  Sol_valida" );
    printf( "\n%f  %f  %f  %s  %f  %d", theta, desp, err1, "--", dt, sol );
    while( aux  ){
        theta += dt;
        desp = busca( theta, ss, &sol );
        err2 = func( theta, desp, ss );
        theta2 = theta;
        printf( "\n%f  %f  %f  %f  %f  %d", theta, desp, err1, err2, dt, sol );
        if( sol ){
            if( conta != 0 ){
                err1   = err2;
                theta1 = theta2;
            }
            else{
                pendiente_recta( theta1, err1, theta2, err2, &hor, &vert );
                // Cambio de signo (excentricidades grandes)
                if( err1 * err2 < 0.0 ){
                    dt /= 1.5;
                    theta = theta1;
                    if( vert || fabs( err1 ) > errMax ){
                        dt = dtIni;
                        theta  = theta2;
                        theta1 = theta2;
                        err1   = err2;
                    }
                }
                else{
                    err1   = err2;
                    theta1 = theta2;
                }
            }
            conta = 0;
        }
        else{
            dt = dtIni;
            conta++;
        }
        // Para que la fuerza axial sea del mismo tipo
        if( fabs( err1 ) < tol ){
            if( ss->Pn * ss->P > 0.0 ) break;
            else{
                dt = dtIni;
                theta  = theta2;
                err1   = err2;
                theta1 = theta2;
            }
        }
        if( theta - thetaIni > 2.1 * pi ){
            printf("\nSOLUCION NO ENCONTRADA\n" );
            break;
        }
    }
	printf( "\nExcentricidades" );
	printf( "\nex = %f", ss->ex );
	printf( "\ney = %f", ss->ey );
	printf( "\nResistencia nominal" );
	printf( "\nPn  = %f", ss->Pn );
	printf( "\nMnx = %f", ss->Mnx );
	printf( "\nMny = %f", ss->Mny );
	printf( "\nPosicion del EN" );
	printf( "\nGiro = %f", ss->theta );
	printf( "\nDesp = %f", ss->desp );
    printf( "\n\n" );
}

double busca_ey( double theta, seccion *ss, int *sol ){
	// Encuentra la ubicación solo para la excentricidad en "y" a través
	// del método de Newton-Rhapson para un ángulo dado.
	int der1, der2;
	double despFP, dy, y1, y2, y;
	double ( *f )( double, double, seccion * );
	f = &funcEy;
	dy = 0.0001; // OJO: si es muy pequeño no encuentra la solución
	// Punto de flexión pura
	despFP = despFlexPura( theta, ss );
	//cout << "Flexion pura = " << despFP << endl;
	// Verifica en cada dirección
	der1 = newtonRaphson( &y1, f, despFP - dy, -dy, theta, ss ); // A la izquierda
	der2 = newtonRaphson( &y2, f, despFP + dy, dy, theta, ss );  // A la derecha
        /*
        printf( "\ntheta  = %f", theta );
        printf( "\ndespFP = %f", despFP );
        printf( "\ny1 = %f", y1 );
        printf( "\ny2 = %f", y2 );
        printf( "\nder1 = %d", der1 );
        printf( "\nder2 = %d", der2 );
        getchar();
        //*/
	solCorrecta( &y, sol, der1, der2, y1, y2, ss );
        /*
        printf( "\ny = %f", y );
        getchar();
        //*/
	//cout << "y   = " << y << endl;
	//cout << "der = " << der << endl;
	return y;
}

double busca_ex( double theta, seccion *ss, int *sol ){
    // Encuentra la ubicación solo para la excentricidad en "x" a través
    // del método de Newton-Rhapson para un ángulo dado.
    int der;
    double despFP, dx, x;
    double ( *f )( double, double, seccion * );
    f = &funcEx;
    dx = 0.0001;
    // Punto de flexión pura
    despFP = despFlexPura( theta, ss );
    //cout << "Flexion pura = " << despFP << endl;
    // Verifica en cada dirección
    der = newtonRaphson( &x, f, despFP - dx, -dx, theta, ss );     // A la izquierda
    if( !der ) newtonRaphson( &x, f, despFP + dx, dx, theta, ss ); // A la derecha
    *sol = der;
    return x;
}

void solCorrecta( double *y, int *sol, const int der1, const int der2,
                  const double y1, const double y2, seccion *ss ){
    // Escoge la solución correcta obtenida del método de Newton-Raphson
    double Pn1, Pn2;
    if( der1 && der2 ){
            //printf("\nPaso: entrada 1");
        eval_resist_seccion( ss->theta, y1, ss );
        Pn1 = ss->Pn;
        eval_resist_seccion( ss->theta, y2, ss );
        Pn2 = ss->Pn;
        if( ss->P * Pn1 > 0.0 ){
            *y = y1;
            *sol = 1;
        }
        else if( ss->P * Pn2 > 0.0 ){
            *y = y2;
            *sol = 1;
        }
        else{
            *sol = 0;
            *y = 0.0;
        }
    }
    else if( der1 ){
        eval_resist_seccion( ss->theta, y1, ss );
        Pn1 = ss->Pn;
        if( ss->P * Pn1 > 0.0 ){
            *y = y1;
            *sol = 1;
        }
         else{
            *sol = 0;
            *y = 0.0;
        }
    }
    else if( der2 ){
        eval_resist_seccion( ss->theta, y2, ss );
        Pn2 = ss->Pn;
        if( ss->P * Pn2 > 0.0 ){
            *y = y2;
            *sol = 1;
        }
         else{
            *sol = 0;
            *y = 0.0;
        }
    }
    else{
        *sol = 0;
        *y = 0.0;
    }
}

double funcPn( double theta, double desp, seccion *ss ){
    // Determina Pn para un theta fijo y un desp dado
    eval_resist_seccion( theta, desp, ss );
    return ss->Pn;
}

double funcEy( double theta, double desp, seccion *ss ){
    // Determina el error para la exc. en "y"
    eval_resist_seccion( theta, desp, ss );
    return ss->ey - ss->eyB;
}

double funcEx( double theta, double desp, seccion *ss ){
    // Determina el error para la exc. en "x"
    eval_resist_seccion( theta, desp, ss );
    return ss->ex - ss->exB;
}

double funcEtotal( double theta, double desp, seccion *ss ){
    // Calcula una resultante del error de las excentricidades
    double ex, ey;
    eval_resist_seccion( theta, desp, ss );
    ex = ss->ex - ss->exB;
    ey = ss->ey - ss->eyB;
    return sqrt( ex * ex + ey * ey );
}

// Calcula un error total
double error_tipo1( double theta, double desp, seccion *ss ){
    double ex, ey;
    eval_resist_seccion( theta, desp, ss );
    ex = ss->ex - ss->exB;
    ey = ss->ey - ss->eyB;
    return -exp( 1.0 / ( sqrt( ex*ex + ey*ey ) ) );
}

double error_tipo2( double theta, double desp, seccion *ss ){
    double ex, ey;
    eval_resist_seccion( theta, desp, ss );
    ex = ss->ex - ss->exB;
    ey = ss->ey - ss->eyB;
    return sqrt( ex*ex + ey*ey );
}

double despFlexPura( double theta, seccion *ss ){
    // Encuentra el punto de flexión pura con Newton-Raphson
    // para un giro dado
    double x, dx;
    dx = 0.0001;
    // Siempre hay que escribir el nombre la clase
    double ( *f )( double, double, seccion * );
    f = &funcPn;
    newtonRaphson( &x, f, 0.0, dx, theta, ss );
        /*printf( "\ntheta = %f", theta );
        printf( "\nx     = %f", x );
        getchar( );*/
    return x;
}

int newtonRaphson( double *xx, double ( *f ) ( double,
                   double, seccion* ), const double xIni,
                   const double dx, double theta, seccion *ss ){
    // Ejecuta este método para la función dada
    double xAux, fx, fxAux, fpx, tol, x;
    //bool encontrado = true;
    int encontrado = 1;
    int conta, maxIter;
    tol = 0.00001;
    maxIter = 100;
    x  = xIni;
    fx = 1.0;
    conta = 0;
        //printf( "\nxIni  = %f", xIni );
        //printf( "\ntheta = %f", theta );
    while( fabs( fx ) >  tol ){
            //printf( "\nx = %f", x );
        // Para aproximar la derivada
        fx = f( theta, x, ss );
        xAux  = x + dx;
        fxAux = f( theta, xAux, ss );
        fpx   = ( fxAux - fx ) / ( xAux - x );
        //xAnt  = x;// ANTES NO ESTABA COMENTADO (FALTA CHECAR)
        x -= fx / fpx;
        if( fabs( x ) > 10000000.0 || conta > maxIter ){
            //encontrado = false;
            encontrado = 0;
            break;
        }
        //if( fabs( fpx ) < tol ) x = xAnt;  // ANTES NO ESTABA COMENTADO (FALTA CHECAR)
        //cout << x << "   " << fpx << endl;
        //cout << x << "\t" << fx << "\t" << fxAux << endl;
        conta++;
    }
    *xx = x;
    return encontrado;
}

// Checa si un flujo a un archivo se abrio
void checa_flujo( FILE *salida, const char * const archivo ){
    if( salida == NULL ){
        printf( "\nError: no se pudo crear o abrir el archivo %s", archivo );
        exit( EXIT_FAILURE );
    }    
}

// Escribe en un archivo el error en las excentricidades
void escribeError( const double deltaT, const double deltaD, const double limInf,
                   const double limSup, seccion *S, const char * const archivo_ex,
                   const char * const archivo_ey, const char * const archivo_e ){
    double theta, desp, ex, ey, e, max;
    FILE *salida_ex, *salida_ey, *salida_e;
    // Intenta abrir el flujo
    salida_ex = fopen( archivo_ex, "w" ); checa_flujo( salida_ex, archivo_ex );
    salida_ey = fopen( archivo_ey, "w" ); checa_flujo( salida_ey, archivo_ey );
    salida_e  = fopen( archivo_e,  "w" ); checa_flujo( salida_e,  archivo_e );
    theta = 0.0;
    while( theta <= 2.0*pi ){
        desp = limInf;
        yMax_zona( &S->region_concr_desp.parte_pos, &max );
        max = limSup;
        while( desp <= max ){
            // Evalua la resistencia y calcula las excentricidades
            eval_resist_seccion( theta, desp, S );
            // Errores
            ex = S->ex - S->exB;
            ey = S->ey - S->eyB;
            e = -exp( 1.0 / ( sqrt( ex*ex + ey*ey ) ) );
            //e = sqrt( ex*ex + ey*ey );
            // Escribe el error
            fprintf( salida_ex, "%lf\t%lf\t%lf\n", theta, desp, ex );
            fprintf( salida_ey, "%lf\t%lf\t%lf\n", theta, desp, ey );
            //if( e < 100.0 )
                fprintf( salida_e,  "%lf\t%lf\t%lf\n", theta, desp, e );
            // Mueve el eje neutro paralelamente
            desp += deltaD;
        }
        fprintf( salida_ex, "\n" );
        fprintf( salida_ey, "\n" );
        fprintf( salida_e,  "\n" );
        // Gira el eje neutro
        theta += deltaT;
    }
    // Cierra el flujo
    fclose( salida_ex );
    fclose( salida_ey );
    fclose( salida_e );
}

void escribeError_GiroFijo( const double theta, const double deltaD, seccion *S,
                   const char * const archivo_ex, const char * const archivo_ey,
                   const char * const archivo_e ){
    double desp, ex, ey, e, limInf, limSup;
    FILE *salida_ex, *salida_ey, *salida_e;
    // Intenta abrir el flujo
    salida_ex = fopen( archivo_ex, "w" ); checa_flujo( salida_ex, archivo_ex );
    salida_ey = fopen( archivo_ey, "w" ); checa_flujo( salida_ey, archivo_ey );
    salida_e  = fopen( archivo_e,  "w" ); checa_flujo( salida_e,  archivo_e );
    // Limites del desplazamiento del EN
    limInf = -80.0;
    yMax_zona( &S->region_concr_desp.parte_pos, &limSup );
    limSup = 80.0;
    desp = limInf;
    while( desp <= limSup ){
        // Evalua la resistencia y calcula las excentricidades
        eval_resist_seccion( theta, desp, S );
        // Errores
        ex = S->ex - S->exB;
        ey = S->ey - S->eyB;
        e = ex*ex + ey*ey;
        // Escribe el error
        fprintf( salida_ex, "%lf\t%lf\n", desp, ex );
        fprintf( salida_ey, "%lf\t%lf\n", desp, ey );
        fprintf( salida_e,  "%lf\t%lf\n", desp, e );
        // Mueve el eje neutro paralelamente
        desp += deltaD;
    }
    fprintf( salida_ex, "\n" );
    fprintf( salida_ey, "\n" );
    fprintf( salida_e,  "\n" );    
    // Cierra el flujo
    fclose( salida_ex );
    fclose( salida_ey );
    fclose( salida_e );
}

// Ordena de menor a mayor un conjunto de valores en funcion de F
void quickSort( const int ini, const int fin, const int n, double *F, double **X ){
    int p;
    if( fin <= ini ) return;
    p = particion( ini, fin, n, F, X );
    quickSort( ini, p-1, n, F, X );
    quickSort( p+1, fin, n, F, X );
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

// Algoritmo del Nelder-Mead
void nelder_mead( const int n, double *x, const double delta,
                  const double alpha, const double gamma,
                  const double beta, const double tau, seccion *S,
                  double (*f)( double, double, seccion* ) ){
    int conta, cambio;
    double **X, *fX, *xR, fxR, *xE, fxE, *xC, fxC, cte, *x_prom, error;
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
    error = 1.0;
    cambio = 0;
    while( conta < 200 && error > 0.00001 ){
        // Evalua la funcion en los puntos del simplejo y en {x}
        for( int i = 0; i < n + 1; i++ ) fX[i] = f( X[i][0], X[i][1], S );
            //imprime_iter( n+1, fX, X, "\nAntes de ordenar:" );
        // Ordena los puntos de mejor a peor valor de aptitud
        quickSort( 0, n, n, fX, X );
        //if( fabs(fX[0]) > 1e10 ) break;
        if( fabs(fX[0]) > 1e10 && cambio == 0 ){
            f = &error_tipo2;
            for( int i = 0; i < n + 1; i++ ) fX[i] = f( X[i][0], X[i][1], S );
            quickSort( 0, n, n, fX, X );
            cambio = 1;
            printf( "\n--- Afinando solucion -----------------------------" );
        }
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
        fxR = f( xR[0], xR[1], S );
        // Actualiza el simplejo:
        // Cuando la reflexion es mejor que el elite actual
        if( fxR < fX[0] ){
            // Se vuelve a reflejar
            for( int j = 0; j < n; j++ )
                xE[j] = xR[j] + gamma*(xR[j]-x_prom[j]);
            fxE = f( xE[0], xE[1], S );
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
            fxC = f( xC[0], xC[1], S );
            if( fxC < fX[n] ) copia_vec_D( n, X[n], xC );
            else{
                // Nuevo simplejo
                for( int i = 1; i < n + 1; i++ ){
                    for( int j = 0; j < n; j++ )
                        X[i][j] = tau*X[0][j] + (1.0-tau)*X[i][j];
                }
                for( int j = 0; j < n; j++ )
                    X[n][j] = tau*X[0][j] + (1.0-tau)*x[j];
            }
        }
        conta++;
        error = fabs(fX[0]);
    }
    printf( "\nNelder Mead:");
    printf( "\n\tIteraciones = %d", conta );
    printf( "\n\tError final = %f", error );
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

void imprime_iter( const int n, double *f, double **X, char *mensaje ){
    printf( "%s", mensaje );
    for( int i = 0; i < n; i++ ){
        printf( "\n%f", f[i] );
        for( int j = 0; j < n - 1; j++ ) printf( "   %f", X[i][j] );
        printf( "\n" );
    }
}