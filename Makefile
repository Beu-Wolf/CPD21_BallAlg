all: ballAlg.cpp gen_points.c
	gcc -lstdc++ -lm -O3 ballAlg.cpp -fopenmp -o ballAlg

test: test.c vectors.h
	gcc test.c -o tests

clean:
	rm -f ballAlg tests
