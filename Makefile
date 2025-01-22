build:
	gcc main.c -o main.o -lqiskit_cext -L./lib -I./include

run: build
	LD_LIBRARY_PATH=./lib ./main.o
