
program test
use ft_blas
implicit none
integer :: i
double precision :: a(6,5) = reshape( (/ (i, i=1,30) /) , (/ 6,5 /))
double precision, dimension(5,3) :: b = reshape( (/ (i*2, i=1,15) /), (/ 5, 3/))
double precision, dimension(6,3) :: d = 0


call ft_dgemm(a, b, d, 2)
! test the result of ft_dgemm against the correct result by matmul function

print *, maxval(d - matmul(a,b))

end program test
