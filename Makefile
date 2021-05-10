FLAGS=-O3 -lm
PAR=-fopenmp
EXTRA=

MAIN_SERIAL_C=ballAlg.c
MAIN_PARALLEL_C=ballAlg_omp.c


SERIAL_C=${SERIAL}.c
SERIAL_OUT=${SERIAL}

all: serial mpi 

serial: ballAlg.c vectors.c common.c build_tree_rec.c
	gcc -DSERIAL ${EXTRA} $^ -o ballAlg ${FLAGS} 

mpi: ballAlg-mpi.c common.c vectors.c build_tree_rec.c
	mpicc $^ -o ballAlg-mpi -lm

profile: EXTRA= -pg
profile: all

debug: EXTRA=-g
debug: all

bench: EXTRA= -DSKIP_DUMP
bench: all



query: ballQuery.c
	gcc -O3 -lm ballQuery.c -o ballQuery



.PHONY: clean-serial
clean:
	rm -f ballAlg ballAlg-mpi
