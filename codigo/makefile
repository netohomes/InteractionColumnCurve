# VARIABLES DE COMPILACIÓN
# Compilador
CC = gcc
# Opciones de compilación
CFLAGS_DEBUG = -std=c99 -Wall -Wextra -pedantic -g
CFLAGS_RELEASE = -std=c99 -Wall -DNDEBUG -O3
LIBS = -lm -fopenmp
# Archivos del proyecto
SOURCES = main.c datos.c geometrias.c materiales.c seccion.c matrix.c memoria.c solver.c
# Nombre del programa
EXE = resist_seccion

all:

debug:
	@echo
	@echo VERSION DEBUG
	@echo
	$(CC) $(SOURCES) $(CFLAGS_DEBUG) $(LIBS) -o $(EXE)

release:
	@echo
	@echo VERSION RELEASE
	@echo
	$(CC) $(SOURCES) $(CFLAGS_RELEASE) $(LIBS) -o $(EXE)	

clean:
	@echo
	@echo BORRANDO EJECUTABLE
	@echo
	rm $(EXE)