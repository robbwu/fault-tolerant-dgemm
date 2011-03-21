f95=gfortran

ft_blas.mod ft_blas.o: ft_blas.f95
	$(f95) -c  ft_blas.f95

test: ft_blas.mod ft_blas.o
	$(f95) -lblas  -o test ft_blas.o test.f95
	./test


