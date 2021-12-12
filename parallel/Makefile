MPICC=mpicc
CFLAGS=-O3

default: all

mpi_sort_build:
	@mkdir -p build
	$(MPICC) $(CFLAGS) -c ./src/mpi_sort.c -o ./build/mpi_sort.o
	$(MPICC) $(CFLAGS) -c ./src/quick_select.c -o ./build/quick_select.o
	$(MPICC) $(CFLAGS) -o ./build/mpi_sort.out ./build/mpi_sort.o ./build/quick_select.o

build_quick_select:
	@mkdir -p build
	$(MPICC) $(CFLAGS) -c ./src/quick_select.c -o ./build/quick_select.o
	$(MPICC) $(CFLAGS) -o ./build/quick_select.out ./build/quick_select.o

.PHONY: clean

all: mpi_sort_build

mpi_sort_run:
	mpirun -hostfile hosts ./build/mpi_sort.out ./data/test.bin

run_quick_select: build_quick_select
	./build/quick_select.out

clean:
	rm -rf ./build/*.out
	rm -rf ./build/*.o

