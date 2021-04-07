FLAGS=-fopenmp -lm
OPT=-O3


all: ballAlg-omp.c gen_points.c
	gcc ${FLAGS} ${OPT} ballAlg-omp.c -o ballAlg-omp

serial: ballAlg.c gen_points.c
	gcc ${FLAGS} ${OPT} ballAlg.c -o ballAlg-serial

profile: ballAlg.c gen_points.c
	gcc ${FLAGS} ${OPT} -pg ballAlg.c     -o ballAlg
	gcc ${FLAGS} ${OPT} -pg ballAlg-omp.c -o ballAlg-omp

debug: ballAlg.c ballAlg-omp.c gen_points.c
	gcc ${FLAGS} -g ballAlg.c     -o ballAlg-serial
	gcc ${FLAGS} -g ballAlg-omp.c -o ballAlg-omp


query: ballQuery.c
	gcc -O3 -lm ballQuery.c -o ballQuery


.PHONY: clean-serial
clean:
	rm -f ballAlg-serial ballQuery ballAlg-omp
