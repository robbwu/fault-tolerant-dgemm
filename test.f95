
program test
use ft_blas
implicit none
integer, parameter :: m = 2048, n = 1024
integer :: i
real :: t_start, t_end
double precision :: a(m+1, m), b(m, m+1), c(m+1, m+1), d(m, m), e(m, m), err(4,4)
double precision :: aa(m, m), bb(m, m)

call init_random_seed()
call random_number(a)
call random_number(b)
aa = a(:m, :)
bb = b(:, :m)
c = 0.0
d = 0.0
e = 0.0


call cpu_time(t_start)
call ft_dgemm(a, b, c, 128)
call cpu_time(t_end)
print *, 'ft_dgemm takes ', t_end - t_start, 'seconds'

call cpu_time(t_start)
call dgemm('n', 'n', m, m, m, 1.0d+0, aa, m, bb, m, 0.0d+0, d, m)
call cpu_time(t_end)
print *, 'standard dgemm takes ', t_end-t_start, ' seconds'

call cpu_time(t_start)
call block_dgemm(aa, bb, e, 128)
call cpu_time(t_end)
print *, 'rank-s dgemm takes ', t_end-t_start, ' seconds'


print *, 'max error of ft_dgemm', maxval( abs(c(1:m, 1:m) - d))
print *, 'max error of block_dgemm', maxval(abs(e - d))


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
