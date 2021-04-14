FLAGS=-fopenmp -lm
OPT=-O3

SERIAL_OUT=ballAlg-serial
TASKS_OUT=ballAlg-tasks
OMP_OUT=ballAlg-omp
QUERY_OUT=ballQuery

all: ballAlg-omp.c gen_points.c
	gcc ${FLAGS} ${OPT} ballAlg-omp.c -o ${OMP_OUT}

serial: ballAlg.c gen_points.c
	gcc ${FLAGS} ${OPT} ballAlg.c -o ${SERIAL_OUT}

tasks: ballAlg-onlytasks.c gen_points.c
	gcc ${FLAGS} ${OPT} ballAlg-onlytasks.c -o ${TASKS_OUT}

profile: ballAlg.c gen_points.c
	gcc ${FLAGS} ${OPT} -pg ballAlg.c -o ${SERIAL_OUT}
	gcc ${FLAGS} ${OPT} -pg ballAlg-omp.c -o ${OMP_OUT}

debug: ballAlg.c gen_points.c
	gcc ${FLAGS} -g ballAlg.c -o ${SERIAL_OUT}
	gcc ${FLAGS} -g ballAlg-omp.c -o ${OMP_OUT}

bench: ballAlg-omp.c ballAlg.c
	gcc ${FLAGS} ${OPT} -DSKIP_DUMP ballAlg.c -o ${SERIAL_OUT}
	gcc ${FLAGS} ${OPT} -DSKIP_DUMP ballAlg-omp.c -o ${OMP_OUT}
	gcc ${FLAGS} ${OPT} -DSKIP_DUMP ballAlg-onlytasks.c -o ${TASKS_OUT}


query: ballQuery.c
	gcc -O3 -lm ballQuery.c -o ${QUERY_OUT}


.PHONY: clean-serial
clean:
	rm -f ${SERIAL_OUT} ${QUERY_OUT} ${OMP_OUT}
