FC= gfortran
FCFLAGS = -O3 -fopenmp

all: spmv33

spmv33: spmv33.f90
	$(FC) $(FCFLAGS) -o spmv33 spmv33.f90

clean:
	rm -f spmv33

