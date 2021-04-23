# choose gcc version
HOST=$(shell sh -c "cat /etc/hostname | sed 's/[0-9]*p.*//'")
GCC=gcc
ifeq ($(HOST),lab)
    GCC=gcc-10
endif

FLAGS=-lm
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
	${GCC} -DSERIAL ${EXTRA} $^ -o ${SERIAL_REC_OUT} ${FLAGS} 

serial-itr: ballAlg.c vectors.c common.c build_tree_itr.c
	${GCC} -DSERIAL ${EXTRA} $^ -o ${SERIAL_ITR_OUT} ${FLAGS} 

parall-rec: ballAlg_omp.c vectors.c common.c build_tree_rec.c
	${GCC} ${PAR} ${EXTRA} $^ -o ${PARALL_REC_OUT} ${FLAGS}

parall-itr: ballAlg_omp.c vectors.c common.c build_tree_itr.c
	${GCC} ${PAR} ${EXTRA} $^ -o ${PARALL_ITR_OUT} ${FLAGS}


profile: EXTRA= -pg -O3
profile: all

debug: EXTRA=-g
debug: all

bench: EXTRA= -DSKIP_DUMP -O3
bench: all



query: ballQuery.c
	${GCC} -O3 ballQuery.c -o ballQuery ${FLAGS}



.PHONY: clean-serial
clean:
	rm -f ${SERIAL_REC_OUT} ${SERIAL_ITR_OUT} ${PARALL_REC_OUT} ${PARALL_ITR_OUT}
