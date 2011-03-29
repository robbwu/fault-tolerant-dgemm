module ft_blas
implicit none
double precision, parameter     :: one = 1.0D+0, zero = 0.0D+0

contains
double precision function infnorm(a)
implicit none
double precision, intent(in) :: a(:, :)
infnorm = maxval( sum(abs(a), dim=2))
end function

! incoming a, b, c should provide room for checksum row and column;
! thus  if a, b, c are of shapes (m, k), (k, n), (m, n) 
! they must allocate extra space to be of shapes (m+1, k), (k, n+1), (m+1, n+1)
subroutine ft_dgemm(a, b, c, s)
implicit none
double precision, intent(inout) :: c(:,:), a(:,:), b(:, :)
integer, intent(in)             :: s
double precision                :: tau
integer                         :: m,k,n
integer                         :: i,j
real                            :: t_start, t_end, cktime
m = size(a, 1)-1; k = size(a, 2); n = size(b, 2)-1
! initialize checksum row/column
a(m+1,:) = sum(a(1:m, :), dim=1); b(:, n+1) = sum(b(:, 1:n), dim=2)
c(m+1,:) = sum(c(1:m, :), dim=1); c(:, n+1) = sum(c(:, 1:n), dim=2)
! threshold to distinguish faults from normal roundoff error
tau = max(m, n, k) * epsilon(tau) * infnorm(a) * infnorm(b)
cktime = 0.0  ! checktime
do i=1,k-mod(k,s),s
    !c(i,i) = 2001.d+0
    call cpu_time(t_start)
    call check()
    call cpu_time(t_end)
    if (i == 1) then
        print *, 'a typical check/recover procedure takes', t_end-t_start, 's'
    endif
    cktime = cktime + t_end - t_start
    call dgemm('n', 'n', m+1, n+1, s, one, a(:, i:i+s-1), m+1, b(i:i+s-1,:), s, one, c, m+1)
end do
j = i + s
if (j <= k) then
    call cpu_time(t_start)
    call check()
    call cpu_time(t_end)
    cktime = cktime + t_end - t_start
    call dgemm('n', 'n', m+1, n+1, k-i, one, a(:, i+1:k), m+1, b(i+1:k, :), k-i, one, c, m+1)
    print *, 'a typical rank-s update dgemm takes ', t_end - t_start, ' seconds'
end if
print *, 'time spend on check/recover', cktime
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
        dr = (c(ii, n+1) - sum(c(ii, 1:n)))
        if ( abs(dr) > tau ) then
            rowfail = .true.
            exit
        endif
    end do
    do jj = 1,n+1
        dc = c(m+1, jj) - sum(c(1:m, jj)) 
        if ( abs(dc) > tau ) then
            colfail = .true.
            exit
        endif
    end do
    ! attempt to recover or signaling failure
    if (rowfail  .and. colfail ) then
        if (ii <= m .and. jj <= n) then
            if (abs(dr-dc) <= tau) then
                c(ii, jj) = c(ii, jj) + dr
            else
                call fail(1)
            endif
        elseif (ii == m+1 .and. jj == n+1) then
            if (abs(dr-dc) <= tau) then
                c(ii, jj) = c(ii, jj) - dr
            else
                call fail(2)
            endif
        elseif (ii == m+1 .and. jj <= n) then
            if (abs(dc+dr) <= tau) then
                c(ii,jj) = c(ii,jj) + dr
            else
                call fail(3)
            endif
        else
            if (abs(dc+dr) <= tau) then
                c(ii,jj) = c(ii,jj) + dc
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

subroutine block_dgemm(a,b,c,s)
implicit none
    double precision, intent(in)    :: a(:,:), b(:, :)
    double precision, intent(inout) :: c(:,:)
    integer, intent(in) :: s
    integer :: m, n, k
    integer :: i, j
    real :: t_start, t_end

    m = size(a, 1); k = size(a, 2); n = size(b, 2)
    do i = 1, k-mod(k,s), s
        call cpu_time(t_start)
        call dgemm('n', 'n', m, n, s, one, a(:, i:i+s-1), m, b(i:i+s-1,:), s, one, c, m)
        call cpu_time(t_end)
        if (i==1) then
            print *, 'typical one rank-k update taks', t_end-t_start, 's'
        end if

    enddo
    j = i + s
    if (j <= k) then
        call dgemm('n', 'n', m, n, j, one, a(:, i+1:k), m, b(i+1:k,  :), j, one, c, m) 
    end if
end subroutine block_dgemm 
end module ft_blas

