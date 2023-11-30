#include "dns_error.h"

!########################################################################
! Building blocks to construct FDMs
! Based on Lagrange polynomial for non-uniform grids
! Calculation of RHS for different stencil lengths and bcs (periodic|biased)
!########################################################################
module FDM_PROCS
    use TLAB_CONSTANTS
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP
    implicit none
    private

    public Pi                ! Product function defined over interval given by idx(:), Pi(x-x_j) for all j in idx
    public Pi_p              ! First-order derivative of Pi
    public Pi_pp_3           ! Second-order derivative when idx has only 3 points
    public Lag               ! Lagrange polynomials on idx(:) around i
    public Lag_p             ! First-order derivative of Lag
    public Lag_pp_3          ! Second-order derivative when idx has only 3 points

    public coef_e1n2_biased  ! coefficients for the biased, 2. order approximation to 1. order derivative
    public coef_e1n3_biased  ! coefficients for the biased, 3. order approximation to 1. order derivative

    ! generic cases
    public MatMul_3d            ! Calculate f = B u, assuming B is tridiagonal with center diagonal equal to 1
    public MatMul_3d_add        ! Calculate f = f + B u, assuming B is tridiagonal
    ! special cases, coefficients are constant in the interior points
    public MatMul_3d_antisym    ! Calculate f = B u, assuming B is tridiagonal, antisymmetric with 1. off-diagonal equal to 1
    public MatMul_3d_sym        ! Calculate f = B u, assuming B is tridiagonal, symmetric with 1. off-diagonal equal to 1

    ! generic cases
    public MatMul_5d            ! Calculate f = B u, assuming B is pentadiagonal with center diagonal is 1
    public MatMul_5d_add        ! Calculate f = f + B u, assuming B is pentadiagonal
    ! special cases, coefficients are constant in the interior points
    public MatMul_5d_antisym    ! Calculate f = B u, assuming B is pentadiagonal, antisymmetric with 1. off-diagonal equal to 1
    public MatMul_5d_sym        ! Calculate f = B u, assuming B is pentadiagonal, symmetric with 1. off-diagonal equal to 1

    ! generic cases
    ! tbd when needed
    ! special cases, coefficients are constant in the interior points
    public MatMul_7d_antisym    ! Calculate f = B u, assuming B is heptadiagonal, antisymmetric with 1. off-diagonal equal to 1
    public MatMul_7d_sym        ! Calculate f = B u, assuming B is heptadiagonal, symmetric with 1. off-diagonal equal to 1

    public FDM_Bcs_Neumann      ! Initialize arrays to impose Neumann Bcs
    public FDM_Bcs_Reduce

! ###################################################################
! Compact parameters (1st derivative of 6th-order pentadiagonal); to be removed
! ###################################################################
    real(wp), public :: C1N6M_ALPHA, C1N6M_BETA
    real(wp), public :: C1N6M_ALPHA2, C1N6M_BETA2
    real(wp), public :: C1N6M_A, C1N6M_B, C1N6M_C
    real(wp), public :: C1N6M_AD2, C1N6M_BD4, C1N6M_CD6
    real(wp), public :: C1N6M_BD2, C1N6M_CD3

contains
    !########################################################################
    !########################################################################
    function Pi(x, j, idx) result(f)    ! Product function on interval idx(:) evaluated at x_j
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, idx(:)
        real(wp) f

        integer(wi) k

        f = 1.0_wp
        do k = 1, size(idx)
            f = f*(x(j) - x(idx(k)))
        end do

        return
    end function

    ! -------------------------------------------------------------------
    function Pi_p(x, j, idx) result(f)
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, idx(:)
        real(wp) f

        real(wp) dummy
        integer(wi) k, m

        f = 0.0_wp
        do k = 1, size(idx)
            dummy = 1.0_wp
            do m = 1, size(idx)
                if (m /= k) then
                    dummy = dummy*(x(j) - x(idx(m)))
                end if
            end do
            f = f + dummy
        end do

        return
    end function

    ! -------------------------------------------------------------------
    function Pi_pp_3(x, j, idx) result(f)       ! It assumes idx has only 3 points
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, idx(:)
        real(wp) f

        f = 2.0_wp*(x(j) - x(idx(1)) + x(j) - x(idx(2)) + x(j) - x(idx(3)))

        return
    end function

!########################################################################
!########################################################################
    function Lag(x, j, i, idx) result(f)        ! Lagrange polynomials on idx(:) around i evaluated at x_j
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j, idx(:)
        real(wp) f

        integer(wi) k

        f = 1.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                f = f*(x(j) - x(idx(k)))/(x(i) - x(idx(k)))
            end if
        end do

        return
    end function

    ! -------------------------------------------------------------------
    function Lag_p(x, j, i, idx) result(f)        ! 1. derivative of Lagrange polynomials on idx(:) around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j, idx(:)
        real(wp) f

        integer(wi) k, m
        real(wp) den, dummy

        den = 1.0_wp
        f = 0.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                dummy = 1.0_wp
                do m = 1, size(idx)
                    if (idx(m) /= i .and. m /= k) then
                        dummy = dummy*(x(j) - x(idx(m)))
                    end if
                end do
                f = f + dummy
                den = den*(x(i) - x(idx(k)))
            end if
        end do
        f = f/den

        return
    end function

    ! -------------------------------------------------------------------
    function Lag_pp_3(x, j, i, idx) result(f)    ! It assumes idx has only 3 points
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, i, idx(:)
        real(wp) f

        integer(wi) k

        f = 2.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                f = f/(x(i) - x(idx(k)))
            end if
        end do

        return
    end function

!########################################################################
!########################################################################
    function coef_e1n3_biased(x, i, backwards) result(coef) ! first-order derivative, explicit, 3. order
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i
        logical, intent(in), optional :: backwards
        real(wp) coef(4)

        integer(wi) stencil(4), k

        if (present(backwards)) then
            stencil = [i, i - 1, i - 2, i - 3]
        else
            stencil = [i, i + 1, i + 2, i + 3]
        end if

        do k = 1, size(stencil)
            coef(k) = Lag_p(x, i, stencil(k), stencil(:))
        end do
        ! if uniform, [ -11/6 3 -3/2 1/3 ]/h

        return
    end function

    ! -------------------------------------------------------------------
    function coef_e1n2_biased(x, i, backwards) result(coef) ! first-order derivative, explicit, 2. order
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i
        logical, intent(in), optional :: backwards
        real(wp) coef(3)

        integer(wi) stencil(3), k

        if (present(backwards)) then
            stencil = [i, i - 1, i - 2]
        else
            stencil = [i, i + 1, i + 2]
        end if

        do k = 1, size(stencil)
            coef(k) = Lag_p(x, i, stencil(k), stencil(:))
        end do
        ! if uniform, [ -3/2 2 -1/2 ]/h

        return
    end function

! #######################################################################
! #######################################################################
! Matrix multiplication of n-diagonal matrix with a vector with special boundary conditions.
! The boundary conditions can extend over n/2+2 points
! This allows use to handle systems A y = B x in which A amd B differ by up to 2 diagonals (see notes)

    ! #######################################################################
    ! Calculate f = B u, assuming B is tri-diagonal with center diagonal is 1
    ! Special boundary conditions restricted to 3 points:
    ! r_11 r_12 r_13
    !      r_21 r_22 r_23
    !      r_30 r_31 r_32 r_33
    !                r_41  1.  r_43         <- interior points start here
    !                     ...  ...  ...
    subroutine MatMul_3d(nx, len, r1, r3, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        integer(wi), intent(in) :: nx, len       ! len linear systems or size nx
        real(wp), intent(in) :: r1(nx), r3(nx)   ! RHS diagonals (#=3-1 because center diagonal is 1)
        real(wp), intent(in) :: u(len, nx)       ! function u
        real(wp), intent(out) :: f(len, nx)      ! RHS, f = B u
        integer, optional :: ibc
        real(wp), intent(in), optional :: rhs_b(1:3, 0:3), rhs_t(0:2, 1:4)  ! Special bcs at bottom and top
        real(wp), optional :: bcs_b(len), bcs_t(len)

        ! -------------------------------------------------------------------
        integer(wi) n
        integer ibc_loc

        ! -------------------------------------------------------------------
        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

#define r0_b(j) rhs_b(j,0)
#define r1_b(j) rhs_b(j,1)
#define r2_b(j) rhs_b(j,2)
#define r3_b(j) rhs_b(j,3)

#define r1_t(j) rhs_t(j,1)
#define r2_t(j) rhs_t(j,2)
#define r3_t(j) rhs_t(j,3)
#define r4_t(j) rhs_t(j,4)

        ! -------------------------------------------------------------------
        ! Boundary; the first 3/2+1+1=3 rows might be different
        if (any([BCS_MIN, BCS_BOTH] == ibc_loc)) then
            if (present(bcs_b)) bcs_b(:) = f(:, 1)*r2_b(1) + u(:, 2)*r3_b(1) + u(:, 3)*r1_b(1) ! r1(1) contains extended stencil
            ! f(1) contains the boundary condition
            f(:, 2) = f(:, 1)*r1_b(2) + u(:, 2)*r2_b(2) + u(:, 3)*r3_b(2)
            f(:, 3) = f(:, 1)*r0_b(3) + u(:, 2)*r1_b(3) + u(:, 3)*r2_b(3) + u(:, 4)*r3_b(3)
        else
            f(:, 1) = u(:, 1) + u(:, 2)*r3(1) + u(:, 3)*r1(1)   ! r1(1) contains extended stencil
            f(:, 2) = u(:, 1)*r1(2) + u(:, 2) + u(:, 3)*r3(2)
            f(:, 3) = u(:, 2)*r1(3) + u(:, 3) + u(:, 4)*r3(3)
        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 4, nx - 3
            f(:, n) = u(:, n - 1)*r1(n) + u(:, n) + u(:, n + 1)*r3(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary; the last 3/2+1+1=3 rows might be different
        if (any([BCS_MAX, BCS_BOTH] == ibc_loc)) then
            ! f(nx) contains the boundary condition
            f(:, nx - 2) = u(:, nx - 3)*r1_t(0) + u(:, nx - 2)*r2_t(0) + u(:, nx - 1)*r3_t(0) + f(:, nx)*r4_t(0)
            f(:, nx - 1) = u(:, nx - 2)*r1_t(1) + u(:, nx - 1)*r2_t(1) + f(:, nx)*r3_t(1)
            if (present(bcs_t)) bcs_t(:) = u(:, nx - 2)*r3_t(2) + u(:, nx - 1)*r1_t(2) + f(:, nx)*r2_t(2) ! r3(nx) contains extended stencil
        else
            f(:, nx - 2) = u(:, nx - 3)*r1(nx - 2) + u(:, nx - 2) + u(:, nx - 1)*r3(nx - 2)
            f(:, nx - 1) = u(:, nx - 2)*r1(nx - 1) + u(:, nx - 1) + u(:, nx)*r3(nx - 1)
            f(:, nx) = u(:, nx - 2)*r3(nx) + u(:, nx - 1)*r1(nx) + u(:, nx) ! r3(nx) contains extended stencil
        end if

#undef r0_b
#undef r1_b
#undef r2_b
#undef r3_b

#undef r1_t
#undef r2_t
#undef r3_t
#undef r4_t

        return
    end subroutine MatMul_3d

    ! #######################################################################
    ! Calculate f = f + B u, assuming B is tri-diagonal
    subroutine MatMul_3d_add(nx, len, r1, r2, r3, u, f)
        integer(wi), intent(in) :: nx, len       ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx)
        real(wp), intent(in) :: u(len, nx)       ! function u
        real(wp), intent(inout) :: f(len, nx)    ! RHS, f = B u

        ! -------------------------------------------------------------------
        integer(wi) n

        ! -------------------------------------------------------------------
        ! Boundary
        n = 1
        f(:, n) = f(:, n) + u(:, n)*r2(n) + u(:, n + 1)*r3(n)

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 2, nx - 1
            f(:, n) = f(:, n) + u(:, n - 1)*r1(n) + u(:, n)*r2(n) + u(:, n + 1)*r3(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        n = nx
        f(:, n) = f(:, n) + u(:, n - 1)*r1(n) + u(:, n)*r2(n)

        return
    end subroutine MatMul_3d_add

    ! #######################################################################
    subroutine MatMul_3d_antisym(nx, len, r1, r2, r3, u, f, periodic, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        integer(wi), intent(in) :: nx, len                      ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx)          ! RHS diagonals
        real(wp), intent(in) :: u(len, nx)                      ! function u
        real(wp), intent(inout) :: f(len, nx)                   ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, optional :: ibc
        real(wp), intent(in), optional :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(inout), optional :: bcs_b(len), bcs_t(len)

        ! -------------------------------------------------------------------
        integer(wi) n
        integer ibc_loc

        ! -------------------------------------------------------------------
        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

#define r1_b(j) rhs_b(j,1)
#define r2_b(j) rhs_b(j,2)
#define r3_b(j) rhs_b(j,3)

#define r1_t(j) rhs_t(j,1)
#define r2_t(j) rhs_t(j,2)
#define r3_t(j) rhs_t(j,3)

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) - u(:, nx)
            f(:, 2) = u(:, 3) - u(:, 1)

        else
            if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc_loc)) then
                ! f(1) contains the boundary condition
                if (present(bcs_b)) bcs_b(:) = f(:, 1)*r2_b(1) + u(:, 2)*r3_b(1) + u(:, 3)*r1_b(1) ! r1(1) contains extended stencil
                f(:, 2) = f(:, 1)*r1_b(2) + u(:, 2)*r2_b(2) + u(:, 3)*r3_b(2)

            else
                f(:, 1) = u(:, 1)*r2(1) + u(:, 2)*r3(1) + u(:, 3)*r1(1) ! r1(1) contains extended stencil
                f(:, 2) = u(:, 1)*r1(2) + u(:, 2)*r2(2) + u(:, 3)*r3(2)

            end if

        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 3, nx - 2
            f(:, n) = u(:, n + 1) - u(:, n - 1)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, nx - 1) = u(:, nx) - u(:, nx - 2)
            f(:, nx) = u(:, 1) - u(:, nx - 1)

        else
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc_loc)) then
                ! f(nx) contains the boundary condition
                f(:, nx - 1) = u(:, nx - 2)*r1_t(1) + u(:, nx - 1)*r2_t(1) + f(:, nx)*r3_t(1)
                if (present(bcs_t)) bcs_t(:) = u(:, nx - 2)*r3_t(2) + u(:, nx - 1)*r1_t(2) + f(:, nx)*r2_t(2) ! r3(nx) contains extended stencil

            else
                f(:, nx - 1) = u(:, nx - 2)*r1(nx - 1) + u(:, nx - 1)*r2(nx - 1) + u(:, nx)*r3(nx - 1)
                f(:, nx) = u(:, nx - 2)*r3(nx) + u(:, nx - 1)*r1(nx) + u(:, nx)*r2(nx)  ! r3(nx) contains extended stencil

            end if

        end if

#undef r1_b
#undef r2_b
#undef r3_b

#undef r1_t
#undef r2_t
#undef r3_t

        return
    end subroutine MatMul_3d_antisym

    ! #######################################################################
    ! Calculate f = B u, assuming B is symmetric tri-diagonal with 1. off-diagonal equal to 1
    subroutine MatMul_3d_sym(nx, len, r1, r2, r3, u, f, periodic, ibc)
        integer(wi), intent(in) :: nx, len                   ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx)    ! RHS diagonals
        real(wp), intent(in) :: u(len, nx)                   ! function u
        real(wp), intent(out) :: f(len, nx)                  ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n
        integer ibc_loc

        ! -------------------------------------------------------------------
        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) + u(:, nx) + u(:, 1)*r2(1)

        else
            f(:, 1) = u(:, 1)*r2(1) + u(:, 2)*r3(1) &
                      + u(:, 3)*r1(1)   ! r1(1) contains 2. superdiagonal to allow for longer stencil at boundary

            if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 2, nx - 1
            f(:, n) = u(:, n + 1) + u(:, n - 1) + u(:, n)*r2(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, nx) = u(:, 1) + u(:, nx - 1) + u(:, nx)*r2(nx)

        else
            f(:, nx) = u(:, nx - 2)*r3(nx) & ! r3(nx) contains 2. subdiagonal to allow for longer stencil at boundary
                       + u(:, nx - 1)*r1(nx) + u(:, nx)*r2(nx)

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_3d_sym

    ! #######################################################################
    ! Calculate f = B u, assuming B is penta-diagonal with center diagonal is 1
    subroutine MatMul_5d(nx, len, r1, r2, r4, r5, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        integer(wi), intent(in) :: nx, len       ! len linear systems or size nx
        real(wp), intent(in) :: r1(nx), r2(nx), r4(nx), r5(nx)    ! RHS diagonals (#=3-1 because center diagonal is 1)
        real(wp), intent(in) :: u(len, nx)       ! function u
        real(wp), intent(out) :: f(len, nx)      ! RHS, f = B u
        integer, optional :: ibc
        real(wp), intent(in), optional :: rhs_b(1:4, 0:5), rhs_t(0:3, 1:6)  ! Special bcs at bottom and top
        real(wp), optional :: bcs_b(len), bcs_t(len)

        ! -------------------------------------------------------------------
        integer(wi) n

        integer ibc_loc

        ! -------------------------------------------------------------------
        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

#define r0_b(j) rhs_b(j,0)
#define r1_b(j) rhs_b(j,1)
#define r2_b(j) rhs_b(j,2)
#define r3_b(j) rhs_b(j,3)
#define r4_b(j) rhs_b(j,4)
#define r5_b(j) rhs_b(j,5)

#define r1_t(j) rhs_t(j,1)
#define r2_t(j) rhs_t(j,2)
#define r3_t(j) rhs_t(j,3)
#define r4_t(j) rhs_t(j,4)
#define r5_t(j) rhs_t(j,5)
#define r6_t(j) rhs_t(j,6)

        ! -------------------------------------------------------------------
        ! Boundary; the first 5/2+1+1=4 rows might be different
        if (any([BCS_MIN, BCS_BOTH] == ibc_loc)) then
            ! f(1) contains the boundary condition
            if (present(bcs_b)) bcs_b(:) = f(:, 1)*r3_b(1) + u(:, 2)*r4_b(1) + u(:, 3)*r5_b(1) + u(:, 4)*r1_b(1) ! r1(1) contains extended stencil
            f(:, 2) =                      f(:, 1)*r2_b(2) + u(:, 2)*r3_b(2) + u(:, 3)*r4_b(2) + u(:, 4)*r5_b(2)
            f(:, 3) =                      f(:, 1)*r1_b(3) + u(:, 2)*r2_b(3) + u(:, 3)*r3_b(3) + u(:, 4)*r4_b(3) + u(:, 5)*r5_b(3)
            f(:, 4) =                      f(:, 1)*r0_b(4) + u(:, 2)*r1_b(4) + u(:, 3)*r2_b(4) + u(:, 4)*r3_b(4) + u(:, 5)*r4_b(4) + u(:, 6)*r5_b(4)
        else
            f(:, 1) = u(:, 1) + u(:, 2)*r4(1) + u(:, 3)*r5(1) + u(:, 4)*r1(1)   ! r1(1) contains extended stencil
            f(:, 2) = u(:, 1)*r2(2) + u(:, 2) + u(:, 3)*r4(2) + u(:, 4)*r5(2)
            f(:, 3) = u(:, 1)*r1(3) + u(:, 2)*r2(3) + u(:, 3) + u(:, 4)*r4(3) + u(:, 5)*r5(3)
            f(:, 4) = u(:, 2)*r1(4) + u(:, 3)*r2(4) + u(:, 4) + u(:, 5)*r4(4) + u(:, 6)*r5(4)
        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 5, nx - 4
            f(:, n) = u(:, n - 2)*r1(n) + u(:, n - 1)*r2(n) + u(:, n) + u(:, n + 1)*r4(n) + u(:, n + 2)*r5(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary; the last 5/2+1+1=4 rows might be different
        if (any([BCS_MAX, BCS_BOTH] == ibc_loc)) then
            ! f(nx) contains the boundary condition
            f(:, nx - 3) = u(:, nx - 5)*r1_t(0) + u(:, nx - 4)*r2_t(0) + u(:, nx - 3)*r3_t(0) + u(:, nx - 2)*r4_t(0) + u(:, nx - 1)*r5_t(0) + f(:, nx)*r6_t(0)
            f(:, nx - 2) =                        u(:, nx - 4)*r1_t(1) + u(:, nx - 3)*r2_t(1) + u(:, nx - 2)*r3_t(1) + u(:, nx - 1)*r4_t(1) + f(:, nx)*r5_t(1)
            f(:, nx - 1) =                                               u(:, nx - 3)*r1_t(2) + u(:, nx - 2)*r2_t(2) + u(:, nx - 1)*r3_t(2) + f(:, nx)*r4_t(2)
            if (present(bcs_t)) bcs_t(:) =                               u(:, nx - 3)*r5_t(3) + u(:, nx - 2)*r1_t(3) + u(:, nx - 1)*r2_t(3) + f(:, nx)*r3_t(3) ! r5(nx) contains extended stencil
        else
            f(:, nx - 3) = u(:, nx - 5)*r1(nx - 3) + u(:, nx - 4)*r2(nx - 3) + u(:, nx - 3) + u(:, nx - 2)*r4(nx - 3) + u(:, nx - 1)*r5(nx - 3)
            f(:, nx - 2) = u(:, nx - 4)*r1(nx - 2) + u(:, nx - 3)*r2(nx - 2) + u(:, nx - 2) + u(:, nx - 1)*r4(nx - 2) + u(:, nx)*r5(nx - 2)
            f(:, nx - 1) = u(:, nx - 3)*r1(nx - 1) + u(:, nx - 2)*r2(nx - 1) + u(:, nx - 1) + u(:, nx)*r4(nx - 1)
            f(:, nx) = u(:, nx - 3)*r5(nx) + u(:, nx - 2)*r1(nx) + u(:, nx - 1)*r2(nx) + u(:, nx) ! r5(nx) contains extended stencil
        end if

#undef r0_b
#undef r1_b
#undef r2_b
#undef r3_b
#undef r4_b
#undef r5_b

#undef r1_t
#undef r2_t
#undef r3_t
#undef r4_t
#undef r5_t
#undef r6_t

        return
    end subroutine MatMul_5d

    ! #######################################################################
    ! Calculate f = f + B u, assuming B is pentadiagonal
    subroutine MatMul_5d_add(nx, len, r1, r2, r3, r4, r5, u, f)
        integer(wi), intent(in) :: nx, len       ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx)
        real(wp), intent(in) :: u(len, nx)       ! function u
        real(wp), intent(inout) :: f(len, nx)    ! RHS, f = B u

        ! -------------------------------------------------------------------
        integer(wi) n

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, 1) = f(:, 1) + u(:, 1)*r3(1) + u(:, 2)*r4(1) + u(:, 3)*r5(1)
        f(:, 2) = f(:, 2) + u(:, 1)*r2(2) + u(:, 2)*r3(2) + u(:, 3)*r4(2) + u(:, 4)*r5(2)

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 3, nx - 2
            f(:, n) = f(:, n) + u(:, n - 2)*r1(n) + u(:, n - 1)*r2(n) + u(:, n)*r3(n) + u(:, n + 1)*r4(n) + u(:, n + 2)*r5(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, nx - 1) = f(:, nx - 1) + u(:, nx - 3)*r1(nx - 1) + u(:, nx - 2)*r2(nx - 1) + u(:, nx - 1)*r3(nx - 1) + u(:, nx)*r4(nx - 1)
        f(:, nx) = f(:, nx) + u(:, nx - 2)*r1(nx) + u(:, nx - 1)*r2(nx) + u(:, nx)*r3(nx)

        return
    end subroutine MatMul_5d_add

    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric penta-diagonal with 1. off-diagonal equal to 1
    ! It also assumes equal coefficients in the 2. off-diagonal for the interior points
    subroutine MatMul_5d_antisym(nx, len, r1, r2, r3, r4, r5, u, f, periodic, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        integer(wi), intent(in) :: nx, len          ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx)  ! RHS diagonals
        real(wp), intent(in) :: u(len, nx)          ! function u
        real(wp), intent(inout) :: f(len, nx)       ! RHS, f = B u; f_1 and f_n can contain neumann bcs
        logical, intent(in) :: periodic
        integer, optional :: ibc
        real(wp), optional :: rhs_b(:, :), rhs_t(:, :)
        real(wp), optional :: bcs_b(len), bcs_t(len)

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) r5_loc     ! 2. off-diagonal
        integer ibc_loc

        ! -------------------------------------------------------------------
        r5_loc = r5(4)      ! The first 3 equations, last 3 equations, can be normalized differently

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

#define r1_b(j) rhs_b(j,1)
#define r2_b(j) rhs_b(j,2)
#define r3_b(j) rhs_b(j,3)
#define r4_b(j) rhs_b(j,4)
#define r5_b(j) rhs_b(j,5)

#define r1_t(j) rhs_t(j,1)
#define r2_t(j) rhs_t(j,2)
#define r3_t(j) rhs_t(j,3)
#define r4_t(j) rhs_t(j,4)
#define r5_t(j) rhs_t(j,5)

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) - u(:, nx) + r5_loc*(u(:, 3) - u(:, nx - 1))
            f(:, 2) = u(:, 3) - u(:, 1) + r5_loc*(u(:, 4) - u(:, nx))
            f(:, 3) = u(:, 4) - u(:, 2) + r5_loc*(u(:, 5) - u(:, 1))

        else
            if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc_loc)) then
                ! f(1) contains boundary condition
                if (present(bcs_b)) bcs_b(:) = f(:, 1)*r3_b(1) + u(:, 2)*r4_b(1) + u(:, 3)*r5_b(1) + u(:, 4)*r1_b(1) ! r1(1) with extended stencil
                f(:, 2) = f(:, 1)*r2_b(2) + u(:, 2)*r3_b(2) + u(:, 3)*r4_b(2) + u(:, 4)*r5_b(2)
                f(:, 3) = f(:, 1)*r1_b(3) + u(:, 2)*r2_b(3) + u(:, 3)*r3_b(3) + u(:, 4)*r4_b(3) + u(:, 5)*r5_b(3)

            else
                f(:, 1) = u(:, 1)*r3(1) + u(:, 2)*r4(1) + u(:, 3)*r5(1) + u(:, 4)*r1(1)   ! r1(1) with extended stencil
                f(:, 2) = u(:, 1)*r2(2) + u(:, 2)*r3(2) + u(:, 3)*r4(2) + u(:, 4)*r5(2)
                f(:, 3) = u(:, 1)*r1(3) + u(:, 2)*r2(3) + u(:, 3)*r3(3) + u(:, 4)*r4(3) + u(:, 5)*r5(3)

            end if

        end if

        ! Interior points
        do n = 4, nx - 3
            f(:, n) = u(:, n + 1) - u(:, n - 1) + r5_loc*(u(:, n + 2) - u(:, n - 2))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 2) = u(:, nx - 1) - u(:, nx - 3) + r5_loc*(u(:, nx) - u(:, nx - 4))
            f(:, nx - 1) = u(:, nx) - u(:, nx - 2) + r5_loc*(u(:, 1) - u(:, nx - 3))
            f(:, nx) = u(:, 1) - u(:, nx - 1) + r5_loc*(u(:, 2) - u(:, nx - 2))

        else
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc_loc)) then
                ! f(n) contains boundary condition
                f(:, nx - 2) = u(:, nx - 4)*r1_t(1) + u(:, nx - 3)*r2_t(1) + u(:, nx - 2)*r3_t(1) + u(:, nx - 1)*r4_t(1) + f(:, nx)*r5_t(1)
                f(:, nx - 1) = u(:, nx - 3)*r1_t(2) + u(:, nx - 2)*r2_t(2) + u(:, nx - 1)*r3_t(2) + f(:, nx)*r4_t(2)
                if (present(bcs_b)) bcs_t(:) = u(:, nx - 3)*r5_t(3) + u(:, nx - 2)*r1_t(3) + u(:, nx - 1)*r2_t(3) + f(:, nx)*r3_t(3) ! r5(nx) with extended stencil

            else
                f(:, nx - 2) = u(:, nx - 4)*r1(nx - 2) + u(:, nx - 3)*r2(nx - 2) + u(:, nx - 2)*r3(nx - 2) + u(:, nx - 1)*r4(nx - 2) + u(:, nx)*r5(nx - 2)
                f(:, nx - 1) = u(:, nx - 3)*r1(nx - 1) + u(:, nx - 2)*r2(nx - 1) + u(:, nx - 1)*r3(nx - 1) + u(:, nx)*r4(nx - 1)
                f(:, nx) = u(:, nx - 3)*r5(nx) + u(:, nx - 2)*r1(nx) + u(:, nx - 1)*r2(nx) + u(:, nx)*r3(nx)! r5(nx) with extended stencil
            end if

        end if

#undef r1_b
#undef r2_b
#undef r3_b
#undef r4_b
#undef r5_b

#undef r1_t
#undef r2_t
#undef r3_t
#undef r4_t
#undef r5_t

        return
    end subroutine MatMul_5d_antisym

    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric penta-diagonal with 1. off-diagonal equal to 1
    subroutine MatMul_5d_sym(nx, len, r1, r2, r3, r4, r5, u, f, periodic, ibc)
        integer(wi), intent(in) :: nx, len       ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx)  ! RHS diagonals
        real(wp), intent(in) :: u(len, nx)       ! function u
        real(wp), intent(out) :: f(len, nx)      ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) r3_loc     ! center diagonal
        real(wp) r5_loc     ! 2. off-diagonal
        integer ibc_loc

        ! -------------------------------------------------------------------
        r5_loc = r5(3)      ! The first 2 equations, last 2 equations, are normalized differently
        r3_loc = r3(3)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! Boundary
        if (periodic) then
            f(:, 1) = r3_loc*u(:, 1) + u(:, 2) + u(:, nx) &
                      + r5_loc*(u(:, 3) + u(:, nx - 1))

            f(:, 2) = r3_loc*u(:, 2) + u(:, 3) + u(:, 1) &
                      + r5_loc*(u(:, 4) + u(:, nx))

        else
            f(:, 1) = u(:, 1)*r3(1) + u(:, 2)*r4(1) + u(:, 3)*r5(1) &
                      + u(:, 4)*r1(1)   ! r1(1) contains 3. superdiagonal to allow for longer stencil at boundary

            f(:, 2) = u(:, 1)*r2(2) + u(:, 2)*r3(2) + u(:, 3)*r4(2) + u(:, 4)*r5(2)

            if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

        end if

        ! Interior points
        do n = 3, nx - 2
            f(:, n) = r3_loc*u(:, n) + u(:, n + 1) + u(:, n - 1) &
                      + r5_loc*(u(:, n + 2) + u(:, n - 2))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 1) = r3_loc*u(:, nx - 1) + u(:, nx) + u(:, nx - 2) &
                           + r5_loc*(u(:, 1) + u(:, nx - 3))

            f(:, nx) = r3_loc*u(:, nx) + u(:, 1) + u(:, nx - 1) &
                       + r5_loc*(u(:, 2) + u(:, nx - 2))

        else
            f(:, nx - 1) = u(:, nx - 3)*r1(nx - 1) + u(:, nx - 2)*r2(nx - 1) + u(:, nx - 1)*r3(nx - 1) &
                           + u(:, nx)*r4(nx - 1)

            f(:, nx) = u(:, nx - 3)*r5(nx) & ! r5(nx) contains 3. subdiagonal to allow for longer stencil at boundary
                       + u(:, nx - 2)*r1(nx) + u(:, nx - 1)*r2(nx) + u(:, nx)*r3(nx)

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_5d_sym

    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric hepta-diagonal with 1. off-diagonal equal to 1
    ! It also assumes equal coefficients in the 2. and 3. off-diagonals for the interior points
    subroutine MatMul_7d_antisym(nx, len, r1, r2, r3, r4, r5, r6, r7, u, f, periodic, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        integer(wi), intent(in) :: nx, len          ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx), r6(nx), r7(nx)  ! RHS diagonals
        real(wp), intent(in) :: u(len, nx)          ! function u
        real(wp), intent(inout) :: f(len, nx)       ! RHS, f = B u; f_1 and f_n can contain neumann bcs
        logical, intent(in) :: periodic
        integer, optional :: ibc
        real(wp), optional :: rhs_b(:, :), rhs_t(:, :)
        real(wp), optional :: bcs_b(len), bcs_t(len)

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) r6_loc, r7_loc     ! 2. and 3. off-diagonal
        integer ibc_loc

        ! -------------------------------------------------------------------
        r6_loc = r6(5)      ! The first 4 equations, last 4 equations, can be normalized differently
        r7_loc = r7(5)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

#define r1_b(j) rhs_b(j,1)
#define r2_b(j) rhs_b(j,2)
#define r3_b(j) rhs_b(j,3)
#define r4_b(j) rhs_b(j,4)
#define r5_b(j) rhs_b(j,5)
#define r6_b(j) rhs_b(j,6)
#define r7_b(j) rhs_b(j,7)

#define r1_t(j) rhs_t(j,1)
#define r2_t(j) rhs_t(j,2)
#define r3_t(j) rhs_t(j,3)
#define r4_t(j) rhs_t(j,4)
#define r5_t(j) rhs_t(j,5)
#define r6_t(j) rhs_t(j,6)
#define r7_t(j) rhs_t(j,7)

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) - u(:, nx) + r6_loc*(u(:, 3) - u(:, nx - 1)) + r7_loc*(u(:, 4) - u(:, nx - 2))
            f(:, 2) = u(:, 3) - u(:, 1) + r6_loc*(u(:, 4) - u(:, nx)) + r7_loc*(u(:, 5) - u(:, nx - 1))
            f(:, 3) = u(:, 4) - u(:, 2) + r6_loc*(u(:, 5) - u(:, 1)) + r7_loc*(u(:, 6) - u(:, nx))
            f(:, 4) = u(:, 5) - u(:, 3) + r6_loc*(u(:, 6) - u(:, 2)) + r7_loc*(u(:, 7) - u(:, 1))

        else
            if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc_loc)) then
                ! f(1) contains boundary condition
                if (present(bcs_b)) bcs_b(:) = f(:, 1)*r4_b(1) + u(:, 2)*r5_b(1) + u(:, 3)*r6_b(1) + u(:, 4)*r7_b(1) + u(:, 5)*r1_b(1)   ! r1(1) with extended stencil
                f(:, 2) = f(:, 1)*r3_b(2) + u(:, 2)*r4_b(2) + u(:, 3)*r5_b(2) + u(:, 4)*r6_b(2) + u(:, 5)*r7_b(2)
                f(:, 3) = f(:, 1)*r2_b(3) + u(:, 2)*r3_b(3) + u(:, 3)*r4_b(3) + u(:, 4)*r5_b(3) + u(:, 5)*r6_b(3) + u(:, 6)*r7_b(3)
                f(:, 4) = f(:, 1)*r1_b(4) + u(:, 2)*r2_b(4) + u(:, 3)*r3_b(4) + u(:, 4)*r4_b(4) + u(:, 5)*r5_b(4) + u(:, 6)*r6_b(4) + u(:, 7)*r7_b(4)

            else
                f(:, 1) = u(:, 1)*r4(1) + u(:, 2)*r5(1) + u(:, 3)*r6(1) + u(:, 4)*r7(1) + u(:, 5)*r1(1)   ! r1(1) with extended stencil
                f(:, 2) = u(:, 1)*r3(2) + u(:, 2)*r4(2) + u(:, 3)*r5(2) + u(:, 4)*r6(2) + u(:, 5)*r7(2)
                f(:, 3) = u(:, 1)*r2(3) + u(:, 2)*r3(3) + u(:, 3)*r4(3) + u(:, 4)*r5(3) + u(:, 5)*r6(3) + u(:, 6)*r7(3)
                f(:, 4) = u(:, 1)*r1(4) + u(:, 2)*r2(4) + u(:, 3)*r3(4) + u(:, 4)*r4(4) + u(:, 5)*r5(4) + u(:, 6)*r6(4) + u(:, 7)*r7(4)

            end if

        end if

        ! Interior points
        do n = 5, nx - 4
            f(:, n) = u(:, n + 1) - u(:, n - 1) + r6_loc*(u(:, n + 2) - u(:, n - 2)) + r7_loc*(u(:, n + 3) - u(:, n - 3))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 3) = u(:, nx - 2) - u(:, nx - 4) + r6_loc*(u(:, nx - 1) - u(:, nx - 5)) + r7_loc*(u(:, nx) - u(:, nx - 6))
            f(:, nx - 2) = u(:, nx - 1) - u(:, nx - 3) + r6_loc*(u(:, nx) - u(:, nx - 4)) + r7_loc*(u(:, 1) - u(:, nx - 5))
            f(:, nx - 1) = u(:, nx) - u(:, nx - 2) + r6_loc*(u(:, 1) - u(:, nx - 3)) + r7_loc*(u(:, 2) - u(:, nx - 4))
            f(:, nx) = u(:, 1) - u(:, nx - 1) + r6_loc*(u(:, 2) - u(:, nx - 2)) + r7_loc*(u(:, 3) - u(:, nx - 3))

        else
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc_loc)) then
                ! f(nx) contains boundary condition
                f(:, nx - 3) = u(:, nx - 6)*r1_t(1) + u(:, nx - 5)*r2_t(1) + u(:, nx - 4)*r3_t(1) + u(:, nx - 3)*r4_t(1) + u(:, nx - 2)*r5_t(1) + u(:, nx - 1)*r6_t(1)+ f(:, nx)*r7_t(1)
              f(:, nx - 2) = u(:, nx - 5)*r1_t(2) + u(:, nx - 4)*r2_t(2) + u(:, nx - 3)*r3_t(2) + u(:, nx - 2)*r4_t(2) + u(:, nx - 1)*r5_t(2) + f(:, nx)*r6_t(2)
                f(:, nx - 1) = u(:, nx - 4)*r1_t(3) + u(:, nx - 3)*r2_t(3) + u(:, nx - 2)*r3_t(3) + u(:, nx - 1)*r4_t(3) + f(:, nx)*r5_t(3)
                if (present(bcs_b)) bcs_t(:) = u(:, nx - 4)*r7_t(4) + u(:, nx - 3)*r1_t(4) + u(:, nx - 2)*r2_t(4) + u(:, nx - 1)*r3_t(4) + f(:, nx)*r4_t(4) ! r7(nx) with extended stencil

            else
                f(:, nx - 3) = u(:, nx - 6)*r1(nx - 3) + u(:, nx - 5)*r2(nx - 3) + u(:, nx - 4)*r3(nx - 3) + u(:, nx - 3)*r4(nx - 3) + u(:, nx - 2)*r5(nx - 3) + u(:, nx - 1)*r6(nx - 3)+ u(:, nx)*r7(nx - 3)
                f(:, nx - 2) = u(:, nx - 5)*r1(nx - 2) + u(:, nx - 4)*r2(nx - 2) + u(:, nx - 3)*r3(nx - 2) + u(:, nx - 2)*r4(nx - 2) + u(:, nx - 1)*r5(nx - 2) + u(:, nx)*r6(nx - 2)
                f(:, nx - 1) = u(:, nx - 4)*r1(nx - 1) + u(:, nx - 3)*r2(nx - 1) + u(:, nx - 2)*r3(nx - 1) + u(:, nx - 1)*r4(nx - 1) + u(:, nx)*r5(nx - 1)
                f(:, nx) = u(:, nx - 4)*r7(nx) + u(:, nx - 3)*r1(nx) + u(:, nx - 2)*r2(nx) + u(:, nx - 1)*r3(nx) + u(:, nx)*r4(nx) ! r7(nx) with extended stencil
            end if

        end if

#undef r1_b
#undef r2_b
#undef r3_b
#undef r4_b
#undef r5_b
#undef r6_b
#undef r7_b

#undef r1_t
#undef r2_t
#undef r3_t
#undef r4_t
#undef r5_t
#undef r6_t
#undef r7_t

        return
    end subroutine MatMul_7d_antisym

    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric hepta-diagonal with 1. superdiagonal equal to 1
    subroutine MatMul_7d_sym(nx, len, r1, r2, r3, r4, r5, r6, r7, u, f, periodic, ibc)
        integer(wi), intent(in) :: nx, len       ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx), r6(nx), r7(nx)  ! RHS diagonals
        real(wp), intent(in) :: u(len, nx)       ! function u
        real(wp), intent(out) :: f(len, nx)      ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) r4_loc     ! center diagonal
        real(wp) r6_loc     ! 2. off-diagonal
        real(wp) r7_loc     ! 3. off-diagonal
        integer ibc_loc

        ! -------------------------------------------------------------------
        r7_loc = r7(4)
        r6_loc = r6(4)      ! The first 3 equations, last 3 equations, are normalized differently
        r4_loc = r4(4)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! Boundary
        if (periodic) then
            f(:, 1) = r4_loc*u(:, 1) + u(:, 2) + u(:, nx) &
                      + r6_loc*(u(:, 3) + u(:, nx - 1)) &
                      + r7_loc*(u(:, 4) + u(:, nx - 2))

            f(:, 2) = r4_loc*u(:, 2) + u(:, 3) + u(:, 1) &
                      + r6_loc*(u(:, 4) + u(:, nx)) &
                      + r7_loc*(u(:, 5) + u(:, nx - 1))

            f(:, 3) = r4_loc*u(:, 3) + u(:, 4) + u(:, 2) &
                      + r6_loc*(u(:, 5) + u(:, 1)) &
                      + r7_loc*(u(:, 6) + u(:, nx))
        else
            f(:, 1) = u(:, 1)*r4(1) + u(:, 2)*r5(1) + u(:, 3)*r6(1) + u(:, 4)*r7(1) &
                      + u(:, 5)*r1(1)   ! r1(1) contains 4. superdiagonal to allow for longer stencil at boundary

            f(:, 2) = u(:, 1)*r3(2) + u(:, 2)*r4(2) + u(:, 3)*r5(2) + u(:, 4)*r6(2) + u(:, 5)*r7(2)

            f(:, 3) = u(:, 1)*r2(3) + u(:, 2)*r3(3) + u(:, 3)*r4(3) + u(:, 4)*r5(3) + u(:, 5)*r6(3) + u(:, 6)*r7(3)

            if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

        end if

        ! Interior points
        do n = 4, nx - 3
            f(:, n) = r4_loc*u(:, n) + u(:, n + 1) + u(:, n - 1) &
                      + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                      + r7_loc*(u(:, n + 3) + u(:, n - 3))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 2) = r4_loc*u(:, nx - 2) + u(:, nx - 1) + u(:, nx - 3) &
                           + r6_loc*(u(:, nx) + u(:, nx - 4)) &
                           + r7_loc*(u(:, 1) + u(:, nx - 5))

            f(:, nx - 1) = r4_loc*u(:, nx - 1) + u(:, nx) + u(:, nx - 2) &
                           + r6_loc*(u(:, 1) + u(:, nx - 3)) &
                           + r7_loc*(u(:, 2) + u(:, nx - 4))

            f(:, nx) = r4_loc*u(:, nx) + u(:, 1) + u(:, nx - 1) &
                       + r6_loc*(u(:, 2) + u(:, nx - 2)) &
                       + r7_loc*(u(:, 3) + u(:, nx - 3))
        else
            f(:, nx - 2) = u(:, nx - 5)*r1(nx - 2) + u(:, nx - 4)*r2(nx - 2) + u(:, nx - 3)*r3(nx - 2) + u(:, nx - 2)*r4(nx - 2) + u(:, nx - 1)*r5(nx - 2) &
                           + u(:, nx)*r6(nx - 2)

            f(:, nx - 1) = u(:, nx - 4)*r1(nx - 1) + u(:, nx - 3)*r2(nx - 1) + u(:, nx - 2)*r3(nx - 1) + u(:, nx - 1)*r4(nx - 1) &
                           + u(:, nx)*r5(nx - 1)

            f(:, nx) = u(:, nx - 4)*r7(nx) & ! r7(nx) contains 4. subdiagonal to allow for longer stencil at boundary
                       + u(:, nx - 3)*r1(nx) + u(:, nx - 2)*r2(nx) + u(:, nx - 1)*r3(nx) + u(:, nx)*r4(nx)

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_7d_sym

! #######################################################################
! #######################################################################
    subroutine FDM_Bcs_Neumann(ibc, lhs, rhs, rhs_b, rhs_t)
        integer, intent(in) :: ibc
        real(wp), intent(inout) :: lhs(:, :)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(inout) :: rhs_b(:, :), rhs_t(:, :)

        integer(wi) idl, ndl, idr, ndr, ir, ic, nx
        real(wp) dummy

        ! -------------------------------------------------------------------
        ndl = size(lhs, 2)
        idl = size(lhs, 2)/2 + 1        ! center diagonal in lhs
        ndr = size(rhs, 2)
        idr = size(rhs, 2)/2 + 1        ! center diagonal in rhs
        nx = size(lhs, 1)               ! # grid points

        ! For A_22, we need idl >= idr -1
        if (idl < idr - 1) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. LHS array is too small.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
        ! For b_21, we need idr >= idl
        if (idr < idl) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. RHS array is too small.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

        ! -------------------------------------------------------------------
        if (any([BCS_ND, BCS_NN] == ibc)) then
            rhs_b(1:idr, 1:ndr) = rhs(1:idr, 1:ndr)

            dummy = 1.0_wp/rhs(1, idr)      ! normalize by r11

            ! reduced array B^R_{22}
            rhs_b(1, 1:ndr) = -rhs_b(1, 1:ndr)*dummy
            do ir = 1, idr - 1              ! rows
                do ic = idr + 1, ndr        ! columns
                    rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir) + rhs_b(1 + ir, idr - ir)*rhs_b(1, ic)
                end do
                ! longer stencil at the boundary
                ic = ndr + 1
                rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir) + rhs_b(1 + ir, idr - ir)*rhs_b(1, 1)
            end do

            ! reduced array A^R_{22}
            lhs(1, 1:ndl) = lhs(1, 1:ndl)*dummy
            do ir = 1, idr - 1              ! rows
                do ic = idl + 1, ndl        ! columns
                    lhs(1 + ir, ic - ir) = lhs(1 + ir, ic - ir) - rhs_b(1 + ir, idr - ir)*lhs(1, ic)
                end do
                ! vector a^R_{21} stored in rhs
                ic = idr
                rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir)*lhs(1, idl)
            end do

            ! finalize vector a^R_{21}
            do ir = 1, idl - 1
                ic = idr
                rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir) - lhs(1 + ir, idl - ir)
            end do

            ! store a_11/b_11 in rhs
            rhs_b(1, idr) = lhs(1, idl)

        end if

        if (any([BCS_DN, BCS_NN] == ibc)) then
            rhs_t(1:idr, 1:ndr) = rhs(nx - idr + 1:nx, 1:ndr)

            dummy = 1.0_wp/rhs(nx, idr)     ! normalize by rnn

            ! reduced array B^R_{11}
            rhs_t(idr, 1:ndr) = -rhs_t(idr, 1:ndr)*dummy
            do ir = 1, idr - 1              ! rows
                do ic = 1, idr - 1          ! columns
                    rhs_t(idr - ir, ic + ir) = rhs(nx - ir, ic + ir) + rhs(nx - ir, idr + ir)*rhs_t(idr, ic)
                end do
                ! longer stencil at the boundary
                ic = 0
                rhs_t(idr - ir, ic + ir) = rhs_t(idr - ir, ic + ir) + rhs(nx - ir, idr + ir)*rhs_t(idr, ndr)
            end do

            ! reduced array A^R_{11}
            lhs(nx, 1:ndl) = lhs(nx, 1:ndl)*dummy
            do ir = 1, idr - 1              ! rows
                do ic = 1, idl - 1          ! columns
                    lhs(nx - ir, ic + ir) = lhs(nx - ir, ic + ir) - rhs(nx - ir, idr + ir)*lhs(nx, ic)
                end do
                ! vector a^R_{1n} stored in rhs
                ic = idr
                rhs_t(idr - ir, ic + ir) = rhs_t(idr - ir, ic + ir)*lhs(nx, idl)
            end do

            ! finalize vector a^R_{1n}
            do ir = 1, idl - 1
                ic = idr
                rhs_t(idr - ir, ic + ir) = rhs_t(idr - ir, ic + ir) - lhs(nx - ir, idl + ir)
            end do

            ! store a_nn/b_nn in rhs
            rhs_t(idr, idr) = lhs(nx, idl)

        end if

        return
    end subroutine FDM_Bcs_Neumann

! #######################################################################
! #######################################################################
    subroutine FDM_Bcs_Reduce(ibc, lhs, rhs, rhs_b, rhs_t)
        integer, intent(in) :: ibc
        real(wp), intent(inout) :: lhs(:, :)
        real(wp), intent(in), optional :: rhs(:, :)
        real(wp), intent(inout), optional :: rhs_b(:, 0:), rhs_t(0:, :)

        integer(wi) idl, ndl, idr, ndr, ir, ic, nx, nx_t
        real(wp) dummy

        ! -------------------------------------------------------------------
        ndl = size(lhs, 2)
        idl = size(lhs, 2)/2 + 1        ! center diagonal in lhs
        ndr = size(rhs, 2)
        idr = size(rhs, 2)/2 + 1        ! center diagonal in rhs
        nx = size(lhs, 1)               ! # grid points
        nx_t = idr                      ! # grid points affected by bcs; for clarity

        ! -------------------------------------------------------------------
        if (any([BCS_MIN, BCS_BOTH] == ibc)) then
            dummy = 1.0_wp/lhs(1, idl)      ! normalize by l11

            ! reduced array A^R_{22}
            lhs(1, 1:ndl) = -lhs(1, 1:ndl)*dummy
            do ir = 1, idl - 1              ! rows
                do ic = idl + 1, ndl        ! columns
                    lhs(1 + ir, ic - ir) = lhs(1 + ir, ic - ir) + lhs(1 + ir, idl - ir)*lhs(1, ic)
                end do
                ic = ndl + 1                ! longer stencil at the boundary
                lhs(1 + ir, ic - ir) = lhs(1 + ir, ic - ir) + lhs(1 + ir, idl - ir)*lhs(1, 1)
            end do

            ! reduced array B^R_{22}
            if (present(rhs_b)) then
                if (size(rhs_b, 1) < max(idl, idr + 1) .or. size(rhs_b, 2) < max(ndl, ndr)) then
                    call TLAB_WRITE_ASCII(efile, __FILE__//'. rhs_b array is too small.')
                    call TLAB_STOP(DNS_ERROR_UNDEVELOP)
                end if

                rhs_b(1:max(idl, idr + 1), 1:ndr) = rhs(1:max(idl, idr + 1), 1:ndr)

                rhs_b(1, 1:ndr) = rhs_b(1, 1:ndr)*dummy
                do ir = 1, idl - 1              ! rows
                    do ic = idr, ndr            ! columns; ic = idr corresponds to vector b^R_{21}
                        rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir) - lhs(1 + ir, idl - ir)*rhs_b(1, ic)
                    end do
                    ic = ndr + 1                ! longer stencil at the boundary
                    rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir) - lhs(1 + ir, idl - ir)*rhs_b(1, 1)
                end do
            end if

        end if

        if (any([BCS_MAX, BCS_BOTH] == ibc)) then
            dummy = 1.0_wp/lhs(nx, idl)     ! normalize by lnn

            ! reduced array A^R_{11}
            lhs(nx, 1:ndl) = -lhs(nx, 1:ndl)*dummy
            do ir = 1, idl - 1              ! rows
                ic = 0                      ! longer stencil at the boundary
                lhs(nx - ir, ic + ir) = lhs(nx - ir, ic + ir) + lhs(nx - ir, idl + ir)*lhs(nx, ndl)
                do ic = 1, idl - 1          ! columns
                    lhs(nx - ir, ic + ir) = lhs(nx - ir, ic + ir) + lhs(nx - ir, idl + ir)*lhs(nx, ic)
                end do
            end do

            ! reduced array B^R_{11}
            if (present(rhs_t)) then
                if (size(rhs_t, 1) < max(idl, idr + 1) .or. size(rhs_t, 2) < max(ndl, ndr)) then
                    call TLAB_WRITE_ASCII(efile, __FILE__//'. rhs_t array is too small.')
                    call TLAB_STOP(DNS_ERROR_UNDEVELOP)
                end if

                rhs_t(nx_t - max(idl, idr + 1) + 1:nx_t, 1:ndr) = rhs(nx - max(idl, idr + 1) + 1:nx, 1:ndr)

                rhs_t(nx_t, 1:ndr) = rhs_t(nx_t, 1:ndr)*dummy
                do ir = 1, idl - 1              ! rows
                    ic = 0                      ! columns; ic = 0 corresponds to longer stencil at the boundary
                    rhs_t(nx_t - ir, ic + ir) = rhs_t(nx_t - ir, ic + ir) - lhs(nx - ir, idl + ir)*rhs_t(nx_t, ndr)
                    do ic = 1, idr              ! ic = idr corresponds to vector b^R_{1n}
                        rhs_t(nx_t - ir, ic + ir) = rhs_t(nx_t - ir, ic + ir) - lhs(nx - ir, idl + ir)*rhs_t(nx_t, ic)
                    end do
                end do
            end if

        end if

        return
    end subroutine FDM_Bcs_Reduce

end module FDM_PROCS
