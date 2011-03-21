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

!assuming a, b are of shapes (m, k) and (k, n)
m = size(a, 1); k = size(a, 2); n = size(b, 2)

! construct threshold distinguishing from roundoff error
tau = max(m, n, k) * epsilon(tau) * maxval(a) * maxval(b)

allocate(ac(m+1,k), br(k, n+1), cf(m+1, n+1))
ac(1:m,:) = a; br(:, 1:n) = b; cf(1:m, 1:n) = c
ac(m+1,:) = sum(a, dim=1); br(:, n+1) = sum(b, dim=2)
cf(m+1,:) = sum(c, dim=1); cf(:, n+1) = sum(cf, dim=2)

do i=1,k,s
    ! artificial corruption
    cf(i,i) = 2001.d+0
    ! check and recover
    call check(cf)
    ! perform the rank-s update by calling dgemm
    j = i + s - 1
    t = s
    if (j > k) then 
        j = k
        t = k - i + 1
    endif
    call dgemm('n', 'n', m+1, n+1, t, one, ac(:, i:j), m+1, br(i:j,:), t, one, cf, m+1)
    !cf = cf + matmul(ac(:, i:j), br(i:j, :))
end do

c = cf(1:m, 1:n)
deallocate(ac, br, cf)
end subroutine ft_dgemm

subroutine check(cf)
implicit none
    double precision, intent(inout) :: cf(:, :)
    double precision :: dr, dc, drr, dcc
    integer ::  ii, jj, iii, jjj
    integer :: rowfail , colfail 
    integer :: m, n

    rowfail=0;colfail=0;iii=0;jjj=0
    m = size(cf, 1) -1; n = size(cf, 2) -1

    ! I want to ensure that if C didn't pass this test
    ! then there MUST be something wrong that probably
    ! requires a complete recomputation

    ! detect and locate the faults
    do ii = 1,m+1
        drr = (cf(ii, n+1) - sum(cf(ii, 1:n)))
        if ( abs(drr) > tau ) then
            rowfail = rowfail + 1
            iii = ii
            dr = drr
        endif
    end do
    do jj = 1,n+1
        dcc = cf(m+1, jj) - sum(cf(1:m, jj)) 
        if ( abs(dcc) > tau ) then
            colfail = colfail + 1
            jjj = jj
            dc = dcc
        endif
    end do
    ! attempt to recover or signaling failure
    if (rowfail == 1 .and. colfail == 1) then
        if (iii <= m .and. jjj <= n) then
            if (abs(dr-dc) <= tau) then
                cf(iii, jjj) = cf(iii, jjj) + dr
            else
                call fail
            endif
        elseif (iii == m+1 .and. jjj == n+1) then
            if (abs(dr-dc) <= tau) then
                cf(iii, jjj) = cf(iii, jjj) - dr
            else
                call fail
            endif
        elseif (iii == m+1 .and. jjj <= n) then
            if (abs(dc+dr) <= tau) then
                cf(iii,jjj) = cf(iii,jjj) + dr
            else
                call fail
            endif
        else
            if (abs(dc+dr) <= tau) then
                cf(iii,jjj) = cf(iii,jjj) + dc
            else
                call fail
            endif
        endif
    elseif (rowfail > 0 .or. colfail > 0) then
        call fail
    endif


end subroutine check

 !in practical situation this subroutine should abort the current
 !computation and start over again
subroutine fail
    print *, 'fatal error: cannot recover. need complete recomputation'
end subroutine fail

end module ft_blas

