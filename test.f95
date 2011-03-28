
program test
use ft_blas
implicit none
integer :: i

double precision, dimension(2048, 2048) :: a, b, c


call random_number(a)
call random_number(b)

call ft_dgemm(a, b, c, 3)
! test the result of ft_dgemm against the correct result by matmul function

print *, maxval(c - matmul(a,b))

end program test
