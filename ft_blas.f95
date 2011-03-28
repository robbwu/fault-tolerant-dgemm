module ft_blas
implicit none
double precision :: tau

contains
subroutine ft_dgemm(a, b, c, s)
implicit none
double precision, intent(in)    :: a(:,:), b(:, :)
double precision, intent(inout) :: c(:,:)
integer                         :: s,t
integer                         :: m,k,n
integer                         :: i,j
double precision, parameter     :: one = 1.0D+0, zero = 0.0D+0
double precision, allocatable   :: ac(:,:), br(:,:), cf(:,:)
real :: t_start, t_end
!assuming a, b are of shapes (m, k) and (k, n)
m = size(a, 1); k = size(a, 2); n = size(b, 2)

! construct threshold distinguishing from roundoff error
!tau = max(m, n, k) * epsilon(tau) * maxval(a) * maxval(b)
tau = 1.0d-05

allocate(ac(m+1,k), br(k, n+1), cf(m+1, n+1))
ac(1:m,:) = a; br(:, 1:n) = b; cf(1:m, 1:n) = c
ac(m+1,:) = sum(a, dim=1); br(:, n+1) = sum(b, dim=2)
cf(m+1,:) = sum(c, dim=1); cf(:, n+1) = sum(cf, dim=2)

do i=1,k,s
    ! artificial corruption
    cf(i,i) = 2001.d+0
    ! check and recover
    call cpu_time(t_start)
    call check()
    call cpu_time(t_end)
    print *, 'check takes ', t_end-t_start, ' seconds'
    ! perform the rank-s update by calling dgemm
    j = i + s - 1
    t = s
    if (j > k) then 
        j = k
        t = k - i + 1
    endif
    call cpu_time(t_start)
    !cf = cf + matmul(ac(:, i:j), br(i:j, :))
    call dgemm('n', 'n', m+1, n+1, t, one, ac(:, i:j), m+1, br(i:j,:), t, one, cf, m+1)
    call cpu_time(t_end)
    print *, 'matmul takes ', t_end - t_start, ' seconds'
end do

c = cf(1:m, 1:n)
deallocate(ac, br, cf)
contains
subroutine check()
implicit none
    double precision :: dr, dc
    integer ::  ii, jj
    logical :: rowfail, colfail
    rowfail = .false.; colfail = .false.

    ! simple test, assume only one element corrupted
    ! detect and locate the faults
    do ii = 1,m+1
        dr = (cf(ii, n+1) - sum(cf(ii, 1:n)))
        if ( abs(dr) > tau ) then
            rowfail = .true.
            exit
        endif
    end do
    do jj = 1,n+1
        dc = cf(m+1, jj) - sum(cf(1:m, jj)) 
        if ( abs(dc) > tau ) then
            colfail = .true.
            exit
        endif
    end do
    ! attempt to recover or signaling failure
    if (rowfail  .and. colfail ) then
        if (ii <= m .and. jj <= n) then
            if (abs(dr-dc) <= tau) then
                cf(ii, jj) = cf(ii, jj) + dr
            else
                call fail(1)
            endif
        elseif (ii == m+1 .and. jj == n+1) then
            if (abs(dr-dc) <= tau) then
                cf(ii, jj) = cf(ii, jj) - dr
            else
                call fail(2)
            endif
        elseif (ii == m+1 .and. jj <= n) then
            if (abs(dc+dr) <= tau) then
                cf(ii,jj) = cf(ii,jj) + dr
            else
                call fail(3)
            endif
        else
            if (abs(dc+dr) <= tau) then
                cf(ii,jj) = cf(ii,jj) + dc
            else
                call fail(4)
            endif
        endif
    elseif (rowfail .or. colfail ) then
        call fail(5)
    endif


end subroutine check

 !in practical situation this subroutine should abort the current
 !computation and start over again
subroutine fail(i)
    integer, intent(in) :: i
    print *, 'fatal error: cannot recover. need complete recomputation', i
end subroutine fail
end subroutine ft_dgemm

end module ft_blas

