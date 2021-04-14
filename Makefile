FLAGS=-fopenmp -lm
EXTRA=
OPT=-O3

SERIAL=ballAlg_serial
SERIAL_C=${SERIAL}.c
SERIAL_OUT=${SERIAL}

RECUR=ballAlg_recur
RECUR_C=${RECUR}.c
RECUR_OUT=${RECUR}

ITER=ballAlg_iter
ITER_C=${ITER}.c
ITER_OUT=${ITER}

QUERY=ballQuery
QUERY_C=${QUERY}.c
QUERY_OUT=${QUERY}


all: serial recursive iterative

serial: ${SERIAL_C} gen_points.c
	gcc ${FLAGS} ${EXTRA} ${SERIAL_C} -o ${SERIAL_OUT}

recursive: ${RECUR_C} gen_points.c
	gcc ${FLAGS} ${EXTRA} ${RECUR_C} -o ${RECUR_OUT}

iterative: ${ITER_C} gen_points.c
	gcc ${FLAGS} ${EXTRA} ${ITER_C} -o ${ITER_OUT}


profile: EXTRA=${OPT} -pg
profile: all

debug: EXTRA=-g
debug: all

bench: EXTRA=${OPT} -DSKIP_DUMP
bench: all



query: ballQuery.c
	gcc -O3 -lm ${QUERY_C} -o ${QUERY_OUT}



.PHONY: clean-serial
clean:
	rm -f ${SERIAL_OUT} ${RECUR_OUT} ${ITER_OUT} ${QUERY_OUT}
