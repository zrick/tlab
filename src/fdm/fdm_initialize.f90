#include "dns_const.h"
#include "dns_error.h"

subroutine FDM_INITIALIZE(x, g, wrk1d)
    use TLAB_CONSTANTS, only: wp, wi, pi_wp, efile, BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_MIN, BCS_MAX
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
#endif
    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: inb_scal, stagger_on
    use TLAB_VARS, only: visc, schmidt
    use TLAB_PROCS
    use FDM_PROCS
    use FDM_ComX_Direct
    use FDM_Com1_Jacobian
    use FDM_Com2_Jacobian
    use FDM_Integrate

    implicit none

    type(grid_dt), intent(inout) :: g
    real(wp), intent(inout) :: x(g%size, g%inb_grid)
    real(wp), intent(inout) :: wrk1d(g%size, 12)

    target x

! -------------------------------------------------------------------
    integer(wi) i, ib, ip, is, ig, nx, ndl, ndr
    integer(wi) nmin, nmax, nsize, bcs_cases(4)
    real(wp) dummy, coef(5)

    integer, parameter :: i0 = 0, i1 = 1

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'Entering '//__FILE__)
#endif

! ###################################################################
    nx = g%size

    ig = 1 ! Accumulating counter to define pointers inside array x

! ###################################################################
    g%nodes => x(:, ig)

    ig = ig + 1

! ###################################################################
! Jacobians: computational grid is uniform
! ###################################################################
    g%jac => x(:, ig:)

    if (nx == 1) then
        g%jac(:, :) = 1.0_wp
        return
    end if

    ! -------------------------------------------------------------------
    ! first derivative
    g%jac(:, 1) = 1.0_wp

    select case (g%mode_fdm1)
    case (FDM_COM4_JACOBIAN, FDM_COM4_DIRECT)
        call FDM_C1N4_Jacobian(nx, g%jac, wrk1d(:, 1), wrk1d(:, 4), g%nb_diag_1, coef)
        call MatMul_3d_antisym(nx, 1, wrk1d(:, 4), wrk1d(:, 5), wrk1d(:, 6), x, g%jac(:, 1), periodic=.false.)

    case (FDM_COM6_JACOBIAN, FDM_COM6_DIRECT)
        call FDM_C1N6_Jacobian(nx, g%jac, wrk1d(:, 1), wrk1d(:, 4), g%nb_diag_1, coef)
        call MatMul_5d_antisym(nx, 1, wrk1d(:, 4), wrk1d(:, 5), wrk1d(:, 6), wrk1d(:, 7), wrk1d(:, 8), x, g%jac(:, 1), periodic=.false.)

    case (FDM_COM6_JACOBIAN_PENTA)
        call FDM_C1N6_Jacobian_Penta(nx, g%jac, wrk1d(:, 1), wrk1d(:, 6), g%nb_diag_1, coef)
   call MatMul_7d_antisym(nx, 1, wrk1d(:, 6), wrk1d(:, 7), wrk1d(:, 8), wrk1d(:, 9), wrk1d(:, 10), wrk1d(:, 11), wrk1d(:, 12), x, g%jac(:, 1), periodic=.false.)

    end select

    select case (g%nb_diag_1(1))
    case (3)
        call TRIDFS(nx, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
        call TRIDSS(nx, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), g%jac(1, 1))
    case (5)
        call PENTADFS2(nx, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
        call PENTADSS2(nx, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), g%jac(1, 1))
    end select

    ! -------------------------------------------------------------------
    ! second derivative
    g%jac(:, 2) = 1.0_wp

    select case (g%mode_fdm2)
    case (FDM_COM4_JACOBIAN)
        call FDM_C2N4_Jacobian(nx, g%jac(:, 2), wrk1d(:, 1), wrk1d(:, 4), g%nb_diag_2, coef)
        call MatMul_5d_sym(nx, 1, wrk1d(:, 4), wrk1d(:, 5), wrk1d(:, 6), wrk1d(:, 7), wrk1d(:, 8), x, g%jac(:, 2), periodic=.false.)

    case (FDM_COM6_JACOBIAN)
        call FDM_C2N6_Jacobian(nx, g%jac(:, 2), wrk1d(:, 1), wrk1d(:, 4), g%nb_diag_2, coef)
        call MatMul_5d_sym(nx, 1, wrk1d(:, 4), wrk1d(:, 5), wrk1d(:, 6), wrk1d(:, 7), wrk1d(:, 8), x, g%jac(:, 2), periodic=.false.)

    case (FDM_COM6_JACOBIAN_HYPER, FDM_COM6_DIRECT, FDM_COM6_JACOBIAN_PENTA)
        call FDM_C2N6_Hyper_Jacobian(nx, g%jac(:, 2), wrk1d(:, 1), wrk1d(:, 4), g%nb_diag_2, coef)
        call MatMul_7d_sym(nx, 1, wrk1d(:, 4), wrk1d(:, 5), wrk1d(:, 6), wrk1d(:, 7), wrk1d(:, 8), wrk1d(:, 9), wrk1d(:, 10), x, g%jac(:, 2), periodic=.false.)

    end select

    select case (g%nb_diag_2(1))
    case (3)
        call TRIDFS(nx, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
        call TRIDSS(nx, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), g%jac(1, 2))
    case (5)
        call PENTADFS2(nx, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
        call PENTADSS2(nx, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), g%jac(1, 2))
    end select

    ! -------------------------------------------------------------------
    ! Saving operations for the time-stability constraint
    g%jac(:, 3) = 1.0_wp/g%jac(:, 1)
    g%jac(:, 4) = g%jac(:, 3)*g%jac(:, 3)

    ig = ig + 4

! ###################################################################
! first-order derivative
! ###################################################################
    g%lhs1 => x(:, ig:)
    ig = ig + 5
    g%rhs1 => x(:, ig:)
    ig = ig + 7

    select case (g%mode_fdm1)
    case (FDM_COM4_JACOBIAN)
        call FDM_C1N4_Jacobian(nx, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef, g%periodic)

    case (FDM_COM6_JACOBIAN)
        call FDM_C1N6_Jacobian(nx, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef, g%periodic)

    case (FDM_COM6_JACOBIAN_PENTA)
        call FDM_C1N6_Jacobian_Penta(nx, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef, g%periodic)

    case (FDM_COM4_DIRECT)
        call FDM_C1N4_Direct(nx, g%nodes, g%lhs1, g%rhs1, g%nb_diag_1)

    case (FDM_COM6_DIRECT)
        call FDM_C1N6_Direct(nx, g%nodes, g%lhs1, g%rhs1, g%nb_diag_1)

    end select
    ndl = g%nb_diag_1(1)    ! for readability of the source code
    ndr = g%nb_diag_1(2)

    ! -------------------------------------------------------------------
    ! LU decomposition and wave numbers
    g%lu1 => x(:, ig:)

    g%lu1(:, 1:g%nb_diag_1(1)) = g%lhs1(:, 1:g%nb_diag_1(1))
    if (g%periodic) then
        select case (g%nb_diag_1(1))
        case (3)
            call TRIDPFS(nx, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3), g%lu1(1, 4), g%lu1(1, 5))
        case (5)
            call PENTADPFS(nx, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3), g%lu1(1, 4), g%lu1(1, 5), g%lu1(1, 6), g%lu1(1, 7))
        end select

        ig = ig + g%nb_diag_1(1) + 2

        ! -------------------------------------------------------------------
        ! wavenumbers
        do i = 1, nx
            if (i <= nx/2 + 1) then
                wrk1d(i, 1) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
            else
                wrk1d(i, 1) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
            end if
        end do

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        g%mwn1 => x(:, ig)

        if (.not. stagger_on) then

            g%mwn1(:) = 2.0_wp*(coef(3)*sin(wrk1d(:, 1)) + coef(4)*sin(2.0_wp*wrk1d(:, 1)) + coef(5)*sin(3.0_wp*wrk1d(:, 1))) &
                        /(1.0_wp + 2.0_wp*coef(1)*cos(wrk1d(:, 1)) + 2.0_wp*coef(2)*cos(wrk1d(:, 1)))

        else ! staggered case has different modified wavenumbers!

            select case (g%mode_fdm1)

            case DEFAULT
                coef = [9.0_wp/62.0_wp, 0.0_wp, 63.0_wp/62.0_wp, 17.0_wp/62.0_wp, 0.0_wp]

            end select

            g%mwn1(:) = 2.0_wp*(coef(3)*sin(1.0_wp/2.0_wp*wrk1d(:, 1)) + coef(4)/3.0_wp*sin(3.0_wp/2.0_wp*wrk1d(:, 1))) &
                        /(1.0_wp + 2.0_wp*coef(1)*cos(wrk1d(:, 1)))

        end if

        ! final calculations because it is mainly used in the Poisson solver like this
        g%mwn1(:) = (g%mwn1(:)/g%jac(1, 1))**2

        ig = ig + 1

        ! -------------------------------------------------------------------
    else                            ! biased,  different BCs
        bcs_cases(1:4) = [BCS_DD, BCS_ND, BCS_DN, BCS_NN]
        do ib = 1, 4
            ip = (ib - 1)*5

            g%lu1(:, ip + 1:ip + g%nb_diag_1(1)) = g%lhs1(:, 1:g%nb_diag_1(1))

            call FDM_Bcs_Neumann(bcs_cases(ib), g%lu1(:, ip + 1:ip + g%nb_diag_1(1)), g%rhs1(:, 1:g%nb_diag_1(2)), g%rhs1_b, g%rhs1_t)

            nmin = 1; nmax = nx
            if (any([BCS_ND, BCS_NN] == bcs_cases(ib))) nmin = nmin + 1
            if (any([BCS_DN, BCS_NN] == bcs_cases(ib))) nmax = nmax - 1
            nsize = nmax - nmin + 1

            select case (g%nb_diag_1(1))
            case (3)
                call TRIDFS(nsize, g%lu1(nmin:, ip + 1), g%lu1(nmin:, ip + 2), g%lu1(nmin:, ip + 3))
            case (5)
                call PENTADFS2(nsize, g%lu1(nmin:, ip + 1), g%lu1(nmin:, ip + 2), g%lu1(nmin:, ip + 3), g%lu1(nmin:, ip + 4), g%lu1(nmin:, ip + 5))
            end select

            ig = ig + 5

        end do

    end if

! ###################################################################
! first-order integrals (cases lambda = 0.0_wp)
! ###################################################################
    g%lhsi => x(:, ig:)
    ig = ig + 7*2
    g%rhsi => x(:, ig:)
    ig = ig + 5*2

    bcs_cases(1:2) = [BCS_MIN, BCS_MAX]
    do ib = 1, 2
        ip = (ib - 1)*g%nb_diag_1(2)

        call FDM_Int1_Initialize(bcs_cases(ib), g%lhs1(:, 1:ndl), g%rhs1(:, 1:ndr), 0.0_wp, &
                                 g%lhsi(:, ip + 1:ip + ndr), g%rhsi(:, (ib - 1)*ndl + 1:(ib - 1)*ndl + ndl), &
                                 g%rhsi_b((ib - 1)*5 + 1:, :), g%rhsi_t((ib - 1)*5:, :))
        ! LU decomposition
        select case (g%nb_diag_1(2))
        case (3)
            call TRIDFS(g%size - 2, g%lhsi(2:, ip + 1), g%lhsi(2:, ip + 2), g%lhsi(2:, ip + 3))
        case (5)
            call PENTADFS(g%size - 2, g%lhsi(2:, ip + 1), g%lhsi(2:, ip + 2), g%lhsi(2:, ip + 3), &
                          g%lhsi(2:, ip + 4), g%lhsi(2:, ip + 5))
        case (7)
            call HEPTADFS(g%size - 2, g%lhsi(2:, ip + 1), g%lhsi(2:, ip + 2), g%lhsi(2:, ip + 3), &
                          g%lhsi(2:, ip + 4), g%lhsi(2:, ip + 5), g%lhsi(2:, ip + 6), g%lhsi(2:, ip + 7))
        end select

    end do

! ###################################################################
! second-order derivative: LU factorization done in routine TRID*FS
! ###################################################################
    g%lhs2 => x(:, ig:)
    ig = ig + 5
    g%rhs2 => x(:, ig:)
    ig = ig + 7 + 5

    select case (g%mode_fdm2)

    case (FDM_COM4_JACOBIAN)
        call FDM_C2N4_Jacobian(nx, g%jac, g%lhs2, g%rhs2, g%nb_diag_2, coef, g%periodic)
        if (.not. g%uniform) g%need_1der = .true.

    case (FDM_COM6_JACOBIAN)
        call FDM_C2N6_Jacobian(nx, g%jac, g%lhs2, g%rhs2, g%nb_diag_2, coef, g%periodic)
        if (.not. g%uniform) g%need_1der = .true.

    case (FDM_COM6_JACOBIAN_HYPER)
        call FDM_C2N6_Hyper_Jacobian(nx, g%jac, g%lhs2, g%rhs2, g%nb_diag_2, coef, g%periodic)
        if (.not. g%uniform) g%need_1der = .true.

    case (FDM_COM4_DIRECT)
        call FDM_C2N4_Direct(nx, g%nodes, g%lhs2, g%rhs2, g%nb_diag_2)

    case (FDM_COM6_DIRECT)
        call FDM_C2N6_Direct(nx, g%nodes, g%lhs2, g%rhs2, g%nb_diag_2)

    end select

    ! -------------------------------------------------------------------
    ! LU decomposition and wave numbers
    g%lu2 => x(:, ig:)

    g%lu2(:, 1:g%nb_diag_2(1)) = g%lhs2(:, 1:g%nb_diag_2(1))
    if (g%periodic) then
        select case (g%nb_diag_2(1))
        case (3)
            call TRIDPFS(nx, g%lu2(1, 1), g%lu2(1, 2), g%lu2(1, 3), g%lu2(1, 4), g%lu2(1, 5))
        end select
        ig = ig + g%nb_diag_2(1) + 2

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        g%mwn2 => x(:, ig)

        g%mwn2(:) = 2.0_wp*(coef(3)*(1.0_wp - cos(wrk1d(:, 1))) + coef(4)*(1.0_wp - cos(2.0_wp*wrk1d(:, 1))) + coef(5)*(1.0_wp - cos(3.0_wp*wrk1d(:, 1)))) &
                    /(1.0_wp + 2.0_wp*coef(1)*cos(wrk1d(:, 1)) + 2.0_wp*coef(2)*cos(2.0_wp*wrk1d(:, 1)))

        g%mwn2(:) = g%mwn2(:)/(g%jac(1, 1)**2)  ! as used in the Helmholtz solver

        ig = ig + 1

        ! -------------------------------------------------------------------
    else
        select case (g%nb_diag_2(1))
        case (3)
            call TRIDFS(nx, g%lu2(1, 1), g%lu2(1, 2), g%lu2(1, 3))
        end select
        ig = ig + g%nb_diag_2(1)

    end if

! ###################################################################
! LU factorization second-order derivative times the diffusivities
! ###################################################################
    g%lu2d => x(:, ig:)

    ip = 0
    do is = 0, inb_scal ! case 0 for the reynolds number
        if (is == 0) then
            dummy = visc
        else
            dummy = visc/schmidt(is)
        end if

        if (g%nb_diag_2(1) /= 3) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Undeveloped for more than 3 LHS diagonals in 2. order derivatives.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

        if (g%periodic) then                        ! Check routines TRIDPFS and TRIDPSS
            g%lu2d(:, ip + 1) = g%lu2(:, 1)         ! matrix L; 1. subdiagonal
            g%lu2d(:, ip + 2) = g%lu2(:, 2)*dummy   ! matrix L; 1/diagonal
            g%lu2d(:, ip + 3) = g%lu2(:, 3)         ! matrix U is the same
            g%lu2d(:, ip + 4) = g%lu2(:, 4)/dummy   ! matrix L; Additional row/column
            g%lu2d(:, ip + 5) = g%lu2(:, 5)         ! matrix U is the same

            ig = ig + 5
            ip = ip + 5

        else                                        ! Check routines TRIDFS and TRIDSS
            g%lu2d(:, ip + 1) = g%lu2(:, 1)         ! matrix L is the same
            g%lu2d(:, ip + 2) = g%lu2(:, 2)*dummy   ! matrix U; 1/diagonal
            g%lu2d(:, ip + 3) = g%lu2(:, 3)/dummy   ! matrix U; 1. superdiagonal

            ig = ig + 3
            ip = ip + 3

        end if

    end do

! ###################################################################
! LU factorization interpolation, done in routine TRID*FS
! ###################################################################
! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
    if ((stagger_on) .and. g%periodic) then
        g%lu0i => x(:, ig:)

        select case (g%mode_fdm1)
        case DEFAULT
            call FDM_C0INT6P_LHS(nx, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3))
        end select
        call TRIDPFS(nx, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3), g%lu0i(1, 4), g%lu0i(1, 5))
        ig = ig + 5
    end if

! ###################################################################
! LU factorization first interp. derivative, done in routine TRID*FS
! ###################################################################
! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
    if ((stagger_on) .and. g%periodic) then
        g%lu1i => x(:, ig:)

        select case (g%mode_fdm1)
        case DEFAULT
            call FDM_C1INT6P_LHS(nx, g%jac, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3))
        end select
        call TRIDPFS(nx, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3), g%lu1i(1, 4), g%lu1i(1, 5))
        ig = ig + 5
    end if

! ###################################################################
! Density correction in anelastic mode
! ###################################################################
    g%rhoinv => x(:, ig)

    g%anelastic = .false. ! Default; activated in FI_BACKGROUND_INITIALIZE

    ig = ig + 1

! ###################################################################
! Check array sizes
! ###################################################################
    if (ig >= g%inb_grid) then
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Grid size incorrect.')
        call TLAB_STOP(DNS_ERROR_DIMGRID)
    end if

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'Leaving '//__FILE__)
#endif

    return
end subroutine FDM_INITIALIZE
