program VPENTAD
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), parameter :: nmax = 32, len = 5
    real(wp), dimension(nmax, 5) :: a, b, c
    real(wp), dimension(len, nmax) :: x, f

    integer(wi) n, ij, seed, im2, im1, ip1, ip2
    real(wp) RAN0, error, sol

! ###################################################################
#define a_a(n) a(n,1)
#define a_b(n) a(n,2)
#define a_c(n) a(n,3)
#define a_d(n) a(n,4)
#define a_e(n) a(n,5)

#define b_a(n) b(n,1)
#define b_b(n) b(n,2)
#define b_c(n) b(n,3)
#define b_d(n) b(n,4)
#define b_e(n) b(n,5)

    seed = 256

! create diagonals
    do n = 1, nmax
!     a(n) = C_3_R/C_11_R/C_4_R
!     b(n) = C_12_R/C_11_R
!     c(n) =-C_51_R/C_22_R
!     d(n) = C_12_R/C_11_R
!     e(n) = C_3_R/C_11_R/C_4_R
        a_a(n) = RAN0(seed)
        a_b(n) = RAN0(seed)
        a_c(n) = RAN0(seed)
        a_d(n) = RAN0(seed)
        a_e(n) = RAN0(seed)
    end do
! check singular case decomposition
!  a(1) = C_0_R; b(1) = C_0_R; c(1) = C_0_R; d(1) = C_0_R; e(1) = C_0_R
!  a(nmax) = C_0_R; b(nmax) = C_0_R; c(nmax) = C_0_R; d(nmax) = C_0_R; e(nmax) = C_0_R
!  c(1) = C_1_R

! padding
    a_a(1) = C_0_R; a_a(2) = C_0_R; 
    a_b(1) = C_0_R
    a_e(nmax - 1) = C_0_R; a_e(nmax) = C_0_R; 
    a_d(nmax) = C_0_R

    do n = 1, nmax
        b_a(n) = a_a(n)
        b_b(n) = a_b(n)
        b_c(n) = a_c(n)
        b_d(n) = a_d(n)
        b_e(n) = a_e(n)
    end do

! create solution
    do n = 1, nmax
        do ij = 1, len
            x(ij, n) = RAN0(seed)
        end do
    end do

! compute forcing term
    do n = 1, nmax
        im2 = mod(n + nmax - 3, nmax) + 1
        im1 = mod(n + nmax - 2, nmax) + 1
        ip1 = mod(n, nmax) + 1
        ip2 = mod(n + 1, nmax) + 1
        write (*, *) im2, im1, n, ip1, ip2
        do ij = 1, len
            f(ij, n) = x(ij, im2)*a_a(n) + x(ij, im1)*a_b(n) + x(ij, n)*a_c(n) + x(ij, ip1)*a_d(n) + x(ij, ip2)*a_e(n)
        end do
    end do

! ###################################################################
! solve system
!  CALL PENTADFS(nmax,      a_a(1), a_b(1), a_c(1), a_d(1), a_e(1))
!  CALL PENTADSS(nmax, len, a_a(1), a_b(1), a_c(1), a_d(1), a_e(1), f)
    call PENTADFS2(nmax, a_a(1), a_b(1), a_c(1), a_d(1), a_e(1))
    call PENTADSS2(nmax, len, a_a(1), a_b(1), a_c(1), a_d(1), a_e(1), f)

! error
    error = C_0_R
    sol = C_0_R
    do n = 1, nmax
        do ij = 1, len
            error = error + (f(ij, n) - x(ij, n))*(f(ij, n) - x(ij, n))
            sol = sol + x(ij, n)*x(ij, n)
        end do
    end do
    write (*, *) 'Solution L2-norm ..:', sqrt(sol)
    write (*, *) 'Relative error ....:', sqrt(error)/sqrt(sol)

    stop
end program VPENTAD
