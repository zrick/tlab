#include "dns_error.h"

!########################################################################
! Initialize arrays to calculate boundary-value problems and integrals
! based on the compact schemes of the classes Au' = Bu and Au'' = Bu
!########################################################################

module FDM_Integrate
    use TLAB_CONSTANTS
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP
    use FDM_PROCS
    implicit none
    private

    public FDM_Int1_Initialize      ! Prepare to solve u' +\lambda u = f
    public FDM_Int2_Initialize      ! Prepare to solve u'' - \lamba^2 u = f

contains
!########################################################################
!#
!# Initialize the solver for the BVP
!#
!#     u'_i + \lamba u_i = f_i  N-1 eqns
!#     u_1 or u_N given         1   eqn
!#     Au' = Bu                 N   eqns
!#
!# starting from generic diagonal matrices A (lhs) and B (rhs).
!#
!# The system of N-1 eqns:
!#
!#                    (B + \lambda A)u = Af
!#
!# is established in this routine (see notes).
!#
!# The system is normalized such that the central diagonal in the new rhs is 1
!#
!########################################################################
    subroutine FDM_Int1_Initialize(ibc, lhs, rhs, lambda, lu, rhs_int)
        integer, intent(in) :: ibc              ! Boundary condition, BCS_DN, BCS_ND
        real(wp), intent(in) :: lhs(:, :)       ! diagonals in lhs, or matrix A
        real(wp), intent(in) :: rhs(:, :)       ! diagonals in rhs, or matrix B
        real(wp), intent(in) :: lambda          ! system constant
        real(wp), intent(out) :: lu(:, :)       ! diagonals in new lhs
        real(wp), intent(out) :: rhs_int(:, :)  ! diagonals in new rhs

        ! -------------------------------------------------------------------
        integer(wi) i
        integer(wi) idl, ndl, idr, ndr, ir, ic, nx, nmin, nmax
        real(wp) dummy

        ! -------------------------------------------------------------------
        ndl = size(lhs, 2)
        idl = size(lhs, 2)/2 + 1        ! center diagonal in lhs
        ndr = size(rhs, 2)
        idr = size(rhs, 2)/2 + 1        ! center diagonal in rhs
        nx = size(lhs, 1)               ! # grid points

! ###################################################################
        ! check sizes
        if (size(lu, 2) < ndr) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong array lu size.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
        if (size(rhs_int, 2) < ndl) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong array rhs_int size.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

        ! new rhs diagonals (array A)
        rhs_int(:, 1:ndl) = lhs(:, 1:ndl)

        ! new lhs diagonals (array C = B + h \lambda A)
        lu(:, 1:ndr) = rhs(:, 1:ndr)

        nmin = 1
        nmax = nx
        select case (ibc)
        case (BCS_MIN)
            nmin = nmin + 1
        case (BCS_MAX)
            nmax = nmax - 1
        end select
        lu(nmin:nmax, idr) = lu(nmin:nmax, idr) + lambda*lhs(nmin:nmax, idl)                        ! center diagonal
        do i = 1, idl - 1                                                                           ! off-diagonals
            lu(nmin + i:nx, idr - i) = lu(nmin + i:nx, idr - i) + lambda*lhs(nmin + i:nx, idl - i)      ! skip the top-left corner
            lu(1:nmax - i, idr + i) = lu(1:nmax - i, idr + i) + lambda*lhs(1:nmax - i, idl + i)         ! skip the bottom-right corner
        end do

        ! The first row in C is different from B, but I do not calculate p'_1 or p'_n in the integral operator
        call FDM_Bcs_Reduce_In_Place(ibc, rhs_int(:, 1:ndl), lu(:, 1:ndr))

        select case (ibc)
        case (BCS_MIN)
            rhs_int(2, 1) = 0.0_wp          ! longer stencil at the boundary; because I am using old version of matmul_3d.... to be removed
            do ir = 1, idr - 1              ! change sign in term for nonzero bc
                lu(1 + ir, idr - ir) = -lu(1 + ir, idr - ir)
            end do
        case (BCS_MAX)
            rhs_int(nx - 1, ndl) = 0.0_wp   ! longer stencil at the boundary
            do ir = 1, idr - 1              ! change sign in term for nonzero bc
                lu(nx - ir, idr + ir) = -lu(nx - ir, idr + ir)
            end do
        end select

        ! -------------------------------------------------------------------
        ! normalization such that new central diagonal in rhs is 1
        do ir = 1, nx
            dummy = 1.0_wp/rhs_int(ir, idl)

            rhs_int(ir, 1:ndl) = rhs_int(ir, 1:ndl)*dummy
            lu(ir, 1:ndr) = lu(ir, 1:ndr)*dummy

        end do

        return
    end subroutine FDM_Int1_Initialize

!########################################################################
!#
!# Initialize the solver for the BVP
!#
!#     u''_i - \lamba^2 u_i = f_i  N-2 eqns
!#     u_1 and u_N given           2   eqns
!#     Au'' = Bu                   N   eqns
!#
!# starting from the matrices A (lhs, tridiagonal) and B (rhs, pentadiagonal, with unitary central diagonal).
!#
!# The system of N-2 eqns:
!#
!#                    (B - \lambda^2 A)u = Af = g
!#
!# is established in this routine, giving diagonals a-e and g (see notes).
!#
!# The system is normalized such that the central diagonal in the new rhs is 1
!#
!########################################################################
    subroutine FDM_Int2_Initialize(imax, x, ibc, lhs, rhs, lambda2, lu, f, bvp_rhs)

        integer(wi), intent(in) :: imax             ! original size; here using only 2:imax-1
        real(wp), intent(in) :: x(imax)             ! grid nodes
        integer, intent(in) :: ibc                  ! Boundary condition, BCS_DD, BCS_DN, BCS_ND, BCS_NN
        real(wp), intent(in) :: lhs(imax, 3)        ! matrix A, tridiagonal
        real(wp), intent(in) :: rhs(imax, 4)        ! matrix B, with unitary central diagonal
        real(wp) lambda2                            ! system constatn
        real(wp), intent(out) :: lu(imax, 5)        ! diagonals in new pentadiagonal lhs
        real(wp), intent(out) :: f(imax, 2)         ! forcing terms for the hyperbolic sine
        real(wp), intent(out) :: bvp_rhs(imax, 2)   ! diagonals to calculate new tridiagonalrhs

! -------------------------------------------------------------------
        integer(wi) i
        real(wp) dummy1, dummy2, pprime, coef(5), l2_inv, l2_min, l2_max

! ###################################################################
        ! new prentadiagonal lhs (array C22R)
        lu(:, 1) = rhs(:, 1)
        lu(:, 2) = rhs(:, 2) - lambda2*lhs(:, 1)
        lu(:, 3) = 1.0_wp - lambda2*lhs(:, 2)
        lu(:, 4) = rhs(:, 3) - lambda2*lhs(:, 3)
        lu(:, 5) = rhs(:, 4)

        ! new tridiagonal rhs (array A22R); new central diagonal is 1 and I only need subdiagonal and superdiagonal
        bvp_rhs(:, 1) = lhs(:, 1)
        bvp_rhs(:, 2) = lhs(:, 3)

        ! Boundary corrections
        i = 2
        dummy1 = lhs(2, 1)/lhs(1, 2)
        lu(i, 3:5) = lu(i, 3:5) - dummy1*lu(i - 1, [4, 5, 1])
        l2_min = lhs(i, 2) - dummy1*lhs(1, 3)       ! central diagonal in reduced lhs
        bvp_rhs(i, 1) = 0.0_wp                      ! See MatMul_3d

        i = imax - 1
        dummy2 = lhs(imax - 1, 3)/lhs(imax, 2)
        lu(i, 1:3) = lu(i, 1:3) - dummy2*lu(i + 1, [5, 1, 2])
        l2_max = lhs(i, 2) - dummy2*lhs(imax, 1)    ! central diagonal in reduced lhs
        bvp_rhs(i, 2) = 0.0_wp                      ! See MatMul_3d

        ! Setting the RHS for the hyperbolic sine; the minus sign in included here to save ops
        f(:, :) = 0.0_wp

        f(1, 1) = 1.0_wp            ! b21R; this element is simply the solution at imin of s(-)
        f(2, 1) = -(rhs(2, 2) - dummy1)
        f(3, 1) = -rhs(3, 1)

        f(imax, 2) = 1.0_wp         ! b2nR; this element is simply the solution at imax of s(+)
        f(imax - 1, 2) = -(rhs(imax - 1, 3) - dummy2)
        f(imax - 2, 2) = -rhs(imax - 2, 4)

        ! -------------------------------------------------------------------
        ! Corrections to the BCS_DD to account for Neumann using third-order fdm for derivative at the boundary
        if (any([BCS_ND, BCS_NN] == ibc)) then
            ! Coefficients in FDM p'_1= b_1 p_1 + b_2 p_2 + b_3 p_3 + b_4 p_4 + a_2 p''_2
            coef(:) = 0.0_wp
            ! coef(1:3) = coef_e1n2_biased(x, 1)                ! second-order
            ! coef(1:4) = coef_e1n3_biased(x, 1)                ! third-order
            coef(1:5) = coef_c1n4_biased(x, 1)                ! fourth-order

            ! Data to calculate p_1 in terms of p_2, p_3, p_4 and p'_1 and p''_2
            lu(1, 2:5) = -coef([5, 2, 3, 4])/coef(1)              ! l_12 contains coefficient for p''_2
            lu(1, 3) = lu(1, 3) + lambda2*lu(1, 2)                ! d + lambda^2h^2 e in notes, e only 1 component
            pprime = 1.0_wp/coef(1)

            ! Derived coefficients; see notes
            lu(2, 3:5) = lu(2, 3:5) - f(2, 1)*lu(1, 3:5)          ! in reduced C matrix; the minus sign comes from def of f1
            lu(3, 2:4) = lu(3, 2:4) - f(3, 1)*lu(1, 3:5)

            l2_min = l2_min + lu(1, 2)*f(2, 1)                     ! in reduced A matrix; the plus sign comes from def of f1
            bvp_rhs(3, 1) = bvp_rhs(3, 1) + lu(1, 2)*f(3, 1)

            f(:, 1) = pprime*f(:, 1)                              ! for particular solutions

        end if
        if (any([BCS_DN, BCS_NN] == ibc)) then
            ! Coefficients in FDM p'_n = b_1 p_n + b_2 p_{n-1} + b_3 p_{n-2} +...
            coef(:) = 0.0_wp
            ! coef(1:3) = coef_e1n2_biased(x, imax, backwards=.true.)
            ! coef(1:4) = coef_e1n3_biased(x, imax, backwards=.true.)
            coef(1:5) = coef_c1n4_biased(x, imax, backwards=.true.)

            ! Data to calculate p_n in terms of p_{n-1}, p_{n-2} and p'_n and p''_{n-1}
            lu(imax, [1, 2, 3, 5]) = -coef([4, 3, 2, 5])/coef(1)                    ! l_n5 contains coefficient for p''_{n-1}
            lu(imax, 3) = lu(imax, 3) + lambda2*lu(imax, 5)                         ! d + lambda^2h^2 e in notes, e only 1 component
            pprime = 1.0_wp/coef(1)

            ! Derived coefficients; see notes
            lu(imax - 1, 1:3) = lu(imax - 1, 1:3) - f(imax - 1, 2)*lu(imax, 1:3)      ! in reduced C matrix; the minus sign comes from def of f2
            lu(imax - 2, 2:4) = lu(imax - 2, 2:4) - f(imax - 2, 2)*lu(imax, 1:3)

            l2_max = l2_max + lu(imax, 5)*f(imax - 1, 2)                              ! in reduced A matrix; the plus sign comes from def of f2
            bvp_rhs(imax - 2, 2) = bvp_rhs(imax - 2, 2) + lu(imax, 5)*f(imax - 2, 2)

            f(:, 2) = pprime*f(:, 2)                                                    ! for particular solutions

        end if

        ! -------------------------------------------------------------------
        ! normalization such that new central diagonal in rhs is 1
        i = 2
        l2_inv = 1.0_wp/l2_min

        bvp_rhs(i, :) = bvp_rhs(i, :)*l2_inv
        lu(i, :) = lu(i, :)*l2_inv
        f(i, :) = f(i, :)*l2_inv

        do i = 3, imax - 2
            l2_inv = 1.0_wp/lhs(i, 2)

            bvp_rhs(i, :) = bvp_rhs(i, :)*l2_inv
            lu(i, :) = lu(i, :)*l2_inv
            f(i, :) = f(i, :)*l2_inv

        end do

        i = imax - 1
        l2_inv = 1.0_wp/l2_max

        bvp_rhs(i, :) = bvp_rhs(i, :)*l2_inv
        lu(i, :) = lu(i, :)*l2_inv
        f(i, :) = f(i, :)*l2_inv

        return
    contains
        !########################################################################
        ! 1. derivatie of interpolation polynomial between equations (15) and (16)
        !    p'_1= b_1 p_1 + b_2 p_2 + b_3 p_3 + b_4 p_4 + a_2 p''_2
        !
        ! Notation in Shukla and Zhong (2005), JCP, 204, 404–429 for the interpolation:
        !
        !       +                    I_n: set of points where the function and derivatives are given
        !   +---+---+---+---...
        !   +       +   +            I_m: set of points where only the function is given.
        !########################################################################
        function coef_c1n4_biased(x, i, backwards) result(coef)
            real(wp), intent(in) :: x(:)
            integer(wi), intent(in) :: i
            logical, intent(in), optional :: backwards
            real(wp) coef(5)

            real(wp) a2, b1, b2, b3, b4
            real(wp) dx1, dx3, dx4
            real(wp) D
            integer(wi) set_m(3), i1, i2, i3, i4

            i1 = i
            if (present(backwards)) then
                ! same as fowards, but changing the signs of the increments w.r.t. i
                ! To understand it, e.g., define a new variable k = -j, where k is the
                ! discrete variable moving around i
                i2 = i - 1
                i3 = i - 2
                i4 = i - 3
            else
                i2 = i + 1
                i3 = i + 2
                i4 = i + 3
            end if
            dx1 = x(i2) - x(i1)
            dx3 = x(i2) - x(i3)
            dx4 = x(i2) - x(i4)
            set_m = [i1, i3, i4]

            ! -------------------------------------------------------------------
            a2 = 0.5_wp*(Pi(x, i1, set_m) - dx1*Pi_p(x, i1, set_m))/Pi_p(x, i2, set_m)

            b2 = Pi_p(x, i1, set_m)*(2.0_wp*Pi_p(x, i2, set_m) + dx1*Pi_pp_3(x, i2, set_m)) &
                 - Pi(x, i1, set_m)*Pi_pp_3(x, i2, set_m)
            b2 = 0.5_wp*b2/Pi(x, i2, set_m)/Pi_p(x, i2, set_m)

            ! -------------------------------------------------------------------
            D = Lag(x, i2, i1, set_m) + dx1*Lag_p(x, i2, i1, set_m)
            b1 = Lag(x, i1, i1, set_m)*(Lag(x, i2, i1, set_m) + 2*dx1*Lag_p(x, i2, i1, set_m)) &
                 - dx1*Lag_p(x, i1, i1, set_m)*(Lag(x, i2, i1, set_m) + dx1*Lag_p(x, i2, i1, set_m))
            b1 = -b1/dx1/D

            D = Lag(x, i2, i3, set_m) + dx3*Lag_p(x, i2, i3, set_m)
            b3 = Lag(x, i1, i3, set_m)*(Lag(x, i2, i3, set_m) + 2*dx1*Lag_p(x, i2, i3, set_m)) &
                 - dx1*Lag_p(x, i1, i3, set_m)*(Lag(x, i2, i3, set_m) + dx1*Lag_p(x, i2, i3, set_m))
            b3 = -b3/dx3/D

            D = Lag(x, i2, i4, set_m) + dx4*Lag_p(x, i2, i4, set_m)
            b4 = Lag(x, i1, i4, set_m)*(Lag(x, i2, i4, set_m) + 2*dx1*Lag_p(x, i2, i4, set_m)) &
                 - dx1*Lag_p(x, i1, i4, set_m)*(Lag(x, i2, i4, set_m) + dx1*Lag_p(x, i2, i4, set_m))
            b4 = -b4/dx4/D

            coef = [b1, b2, b3, b4, a2]

            ! if uniform, we should have ( -29/6 54/6 -27/6 2/6 )/h and 3h
            ! print*, [b1, b2, b3, b4]*(x(2)-x(1))
            ! print*, a2/(x(2)-x(1))

            return
        end function

    end subroutine FDM_Int2_Initialize

end module FDM_Integrate
