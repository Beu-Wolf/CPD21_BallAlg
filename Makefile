FLAGS=-lm -O3
PAR=-fopenmp
EXTRA=

MAIN_SERIAL_C=ballAlg.c
MAIN_PARALLEL_C=ballAlg_omp.c


SERIAL_C=${SERIAL}.c
SERIAL_OUT=${SERIAL}

SERIAL_REC_OUT=ballAlg_serial_rec
SERIAL_ITR_OUT=ballAlg_serial_itr
PARALL_REC_OUT=ballAlg_parall_rec
PARALL_ITR_OUT=ballAlg_parall_itr


all: serial-rec serial-itr parall-rec parall-itr

serial-rec: ballAlg.c vectors.c common.c build_tree_rec.c
	gcc ${FLAGS} ${EXTRA} $^ -o ${SERIAL_REC_OUT}

serial-itr: ballAlg.c vectors.c common.c build_tree_itr.c
	gcc ${FLAGS} ${EXTRA} $^ -o ${SERIAL_ITR_OUT}

parall-rec: ballAlg_omp.c vectors.c common.c build_tree_rec.c
	gcc ${FLAGS} ${PAR} ${EXTRA} $^ -o ${PARALL_REC_OUT}

parall-itr: ballAlg_omp.c vectors.c common.c build_tree_itr.c
	gcc ${FLAGS} ${PAR} ${EXTRA} $^ -o ${PARALL_ITR_OUT}


profile: EXTRA= -pg
profile: all

debug: EXTRA=-g
debug: all

bench: EXTRA= -DSKIP_DUMP
bench: all



query: ballQuery.c
	gcc -O3 -lm ${QUERY_C} -o ${QUERY_OUT}



.PHONY: clean-serial
clean:
	rm -f ${SERIAL_REC_OUT} ${SERIAL_ITR_OUT} ${PARALL_REC_OUT} ${PARALL_ITR_OUT}
