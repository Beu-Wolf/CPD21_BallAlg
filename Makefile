all: ballAlg-omp.c gen_points.c
	gcc -O3 ballAlg-omp.c -fopenmp -o ballAlg-omp -lm

serial: ballAlg.c gen_points.c
	gcc -O3 ballAlg.c -fopenmp -o ballAlg-serial -lm

profile: ballAlg.c gen_points.c
	gcc -lm -O3 -pg ballAlg.c -fopenmp -o ballAlg

debug: ballAlg.c gen_points.c
	gcc -lstdc++ -lm -g ballAlg.c -fopenmp -o ballAlg

query: ballQuery.c
	gcc -O3 -lm ballQuery.c -o ballQuery

test: test.c vectors.h
	gcc test.c -o test

.PHONY: clean
clean:
	rm -f ballAlg-serial test ballQuery ballAlg-omp
