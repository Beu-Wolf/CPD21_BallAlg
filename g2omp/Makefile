# choose gcc version
HOST=$(shell sh -c "cat /etc/hostname | sed 's/[0-9]*p.*//'")
GCC=gcc
ifeq ($(HOST),lab)
    GCC=gcc-10
endif

FLAGS=-lm
PAR=-fopenmp
EXTRA= -O3

MAIN_SERIAL=ballAlg
MAIN_PARALLEL=ballAlg-omp

all: serial-rec parall-rec

serial-rec: ballAlg.c vectors.c common.c build_tree_rec.c
	${GCC} -DSERIAL ${EXTRA} $^ -o ${MAIN_SERIAL} ${FLAGS} 

parall-rec: ballAlg-omp.c vectors.c common.c build_tree_rec.c
	${GCC} ${PAR} ${EXTRA} $^ -o ${MAIN_PARALLEL} ${FLAGS}


profile: EXTRA= -pg -O3
profile: all

debug: EXTRA=-g
debug: all

bench: EXTRA= -DSKIP_DUMP -O3
bench: all



query: ballQuery.c
	${GCC} -O3 -lm ${QUERY_C} -o ${QUERY_OUT}



.PHONY: clean-serial
clean:
	rm -f ${MAIN_SERIAL} ${MAIN_PARALLEL}
