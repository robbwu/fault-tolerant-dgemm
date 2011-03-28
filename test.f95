
program test
use ft_blas
implicit none
integer, parameter :: m = 2048
integer :: i

double precision, dimension(m, m) :: a, b, c,d

call init_random_seed()
call random_number(a)
call random_number(b)
c = 0.0
d = 0.0

call ft_dgemm(a, b, c, 512)
! test the result of ft_dgemm against the correct result by matmul function

!print *, maxval(c - matmul(a,b))
call dgemm('n', 'n', m, m, m, 1.0d+0, a, m, b, m, 0.0d+0, d, m)
print *, maxval(c - d )

contains
subroutine init_random_seed()
implicit none
INTEGER :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed

CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))

CALL SYSTEM_CLOCK(COUNT=clock)

seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)

DEALLOCATE(seed)
END SUBROUTINE


end program test
