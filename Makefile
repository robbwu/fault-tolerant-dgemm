f95=gfortran
#f95=gcc46

ft_blas.mod ft_blas.o: ft_blas.f95
	$(f95) -c   -O3 ft_blas.f95

test: ft_blas.mod ft_blas.o test.f95
	$(f95) -lblas   -O3 -o test ft_blas.o test.f95
	./test

clean:
	rm *.o *.mod test

