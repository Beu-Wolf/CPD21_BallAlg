
all: ballAlg.cpp gen_points.c
	gcc -lstdc++ -lm -O3 ballAlg.cpp -fopenmp -o ballAlg

debug: ballAlg.cpp gen_points.c
	gcc -lstdc++ -lm -O3 -g ballAlg.cpp -fopenmp -o ballAlg

test: test.c vectors.h
	gcc test.c -o test

.PHONY: clean
clean:
	rm -f ballAlg test
