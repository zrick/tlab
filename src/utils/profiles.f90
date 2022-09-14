#include "dns_const.h"
#include "dns_error.h"

subroutine PROFILES_READBLOCK(bakfile, inifile, block, tag, var)
    use TLAB_TYPES, only: profiles_dt, cp
    use TLAB_CONSTANTS
    use TLAB_PROCS
    implicit none

    character(len=*), intent(in) :: bakfile, inifile, block, tag
    type(profiles_dt), intent(out) :: var

    character(len=512) sRes
    real(cp) derivative

    ! -------------------------------------------------------------------
    call TLAB_WRITE_ASCII(bakfile, '#Profile'//trim(adjustl(tag))//'=<None/Tanh/Erf/Ekman/Parabolic/...>')
    call TLAB_WRITE_ASCII(bakfile, '#'//trim(adjustl(tag))//'=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#YMeanRelative'//trim(adjustl(tag))//'=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Diam'//trim(adjustl(tag))//'=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Thick'//trim(adjustl(tag))//'=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Delta'//trim(adjustl(tag))//'=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#BottomSlope'//trim(adjustl(tag))//'=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#UpperSlope'//trim(adjustl(tag))//'=<value>')

    ! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, inifile, block, 'Profile'//trim(adjustl(tag)), 'none', sRes)
    if (trim(adjustl(sRes)) == 'none') then; var%type = PROFILE_NONE
    else if (trim(adjustl(sRes)) == 'tanh') then; var%type = PROFILE_TANH
    else if (trim(adjustl(sRes)) == 'tanhsymmetric') then; var%type = PROFILE_TANH_SYM
    else if (trim(adjustl(sRes)) == 'tanhantisymmetric') then; var%type = PROFILE_TANH_ANTISYM
    else if (trim(adjustl(sRes)) == 'linear') then; var%type = PROFILE_LINEAR
    else if (trim(adjustl(sRes)) == 'linearcrop') then; var%type = PROFILE_LINEAR_CROP
    else if (trim(adjustl(sRes)) == 'erf') then; var%type = PROFILE_ERF
    else if (trim(adjustl(sRes)) == 'erfsurface') then; var%type = PROFILE_ERF_SURFACE
    else if (trim(adjustl(sRes)) == 'erfantisym') then; var%type = PROFILE_ERF_ANTISYM
    else if (trim(adjustl(sRes)) == 'bickley') then; var%type = PROFILE_BICKLEY
    else if (trim(adjustl(sRes)) == 'gaussian') then; var%type = PROFILE_GAUSSIAN
    else if (trim(adjustl(sRes)) == 'ekman') then; var%type = PROFILE_EKMAN_U
    else if (trim(adjustl(sRes)) == 'ekmanp') then; var%type = PROFILE_EKMAN_U_P
    else if (trim(adjustl(sRes)) == 'parabolic') then; var%type = PROFILE_PARABOLIC
    else if (trim(adjustl(sRes)) == 'mixedlayer') then; var%type = PROFILE_MIXEDLAYER
! the following 2 are used in initialize/flow/pressure_mean; should be cleaned
    else if (trim(adjustl(sRes)) == 'enthalpyerf') then; var%type = -PROFILE_ERF
    else
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong '//trim(adjustl(tag))//' profile.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, block, 'Mean'//trim(adjustl(tag)), 'void', sRes)
    if (trim(adjustl(sRes)) == 'void') then ! Backwards compatibility
        call SCANINIREAL(bakfile, inifile, block, trim(adjustl(tag)), '0.0', var%mean)
    else
        call SCANINIREAL(bakfile, inifile, block, 'Mean'//trim(adjustl(tag)), '0.0', var%mean)
    end if

    call SCANINICHAR(bakfile, inifile, block, 'YMean'//trim(adjustl(tag)), 'void', sRes)
    if (trim(adjustl(sRes)) == 'void') then
        var%relative = .true.
        call SCANINIREAL(bakfile, inifile, block, 'YMeanRelative'//trim(adjustl(tag)), '0.5', var%ymean_rel)    ! Position in relative coordinates
        ! Backwards compatibility
        call SCANINICHAR(bakfile, inifile, block, 'YCoor'//trim(adjustl(tag)), 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            call SCANINIREAL(bakfile, inifile, block, 'YCoor'//trim(adjustl(tag)), '0.5', var%ymean_rel)
            call TLAB_WRITE_ASCII(wfile, 'Update tag YCoor to YMeanRelative.')
        end if
    else
        var%relative = .false.
        call SCANINIREAL(bakfile, inifile, block, 'YMean'//trim(adjustl(tag)), '0.0', var%ymean)         ! Position in absolute coordinates
    end if

    call SCANINIREAL(bakfile, inifile, block, 'Thick'//trim(adjustl(tag)), '0.0', var%thick)
    call SCANINIREAL(bakfile, inifile, block, 'Delta'//trim(adjustl(tag)), '0.0', var%delta)
    ! alternative to provide the variable thick in terms of the maximum derivative
    call SCANINICHAR(bakfile, inifile, block, 'Derivative'//trim(adjustl(tag)), 'void', sRes)
    if (trim(adjustl(sRes)) /= 'void') then
        call SCANINIREAL(bakfile, inifile, block, 'Derivative'//trim(adjustl(tag)), '0.0', derivative)
        call PROFILES_DERTOTHICK(derivative, var)
    end if

    call SCANINIREAL(bakfile, inifile, block, 'BottomSlope'//trim(adjustl(tag)), '0.0', var%lslope)
    call SCANINIREAL(bakfile, inifile, block, 'UpperSlope'//trim(adjustl(tag)), '0.0', var%uslope)
    call SCANINIREAL(bakfile, inifile, block, 'Diam'//trim(adjustl(tag)), '0.0', var%diam)

    call SCANINIREAL(bakfile, inifile, block, 'SurfaceThick'//trim(adjustl(tag)), '1.0', var%parameters(3))
    call SCANINIREAL(bakfile, inifile, block, 'SurfaceDelta'//trim(adjustl(tag)), '0.0', var%parameters(4))
    ! alternative to provide the variable thick in terms of the maximum derivative
    call SCANINICHAR(bakfile, inifile, block, 'SurfaceDerivative'//trim(adjustl(tag)), 'void', sRes)
    if (trim(adjustl(sRes)) /= 'void') then
        call SCANINIREAL(bakfile, inifile, block, 'SurfaceDerivative'//trim(adjustl(tag)), '0.0', derivative)
        call PROFILES_DERTOTHICK(derivative, var)
    end if

    call SCANINIREAL(bakfile, inifile, block, 'ScaleHeight', '0.0', var%parameters(5))

    return
end subroutine PROFILES_READBLOCK

function PROFILES(var, y) result(f)
    use TLAB_TYPES, only: profiles_dt, cp
    use TLAB_CONSTANTS
    implicit none

    type(profiles_dt), intent(in) :: var
    real(cp), intent(in) :: y
    real(cp) f

    ! -------------------------------------------------------------------
    real(cp) yrel, xi, amplify, zamp, cnought

    ! ###################################################################
    yrel = y - var%ymean ! position relative to ycenter
    amplify = 0.0_cp    ! default

    ! -------------------------------------------------------------------
    ! base state varying between two constant levels
    ! -------------------------------------------------------------------
    if (var%thick == 0.0_cp) then
        if (var%type > 0) amplify = 0.5_cp*sign(1.0_cp, yrel)

    else
        xi = yrel/var%thick

        select case (var%type)

        case (PROFILE_LINEAR)
            amplify = -xi

        case (PROFILE_TANH)
            amplify = 0.5_cp*tanh(-0.5_cp*xi)

        case (PROFILE_TANH_SYM)
            amplify = 0.5_cp*(tanh(-0.5_cp*(xi - 0.5_cp*var%diam/var%thick)) &
                              + tanh(0.5_cp*(xi + 0.5_cp*var%diam/var%thick)) - 1.0_cp)

        case (PROFILE_TANH_ANTISYM)
            amplify = 0.25_cp*(tanh(-0.5_cp*(xi - 0.5_cp*var%diam/var%thick)) &
                               - tanh(0.5_cp*(xi + 0.5_cp*var%diam/var%thick)))

        case (PROFILE_ERF, PROFILE_ERF_ANTISYM, PROFILE_ERF_SURFACE)
            amplify = 0.5_cp*erf(-0.5_cp*xi)

        case (PROFILE_PARABOLIC, PROFILE_PARABOLIC_SURFACE)
            amplify = (1.0_cp + 0.5_cp*xi)*(1.0_cp - 0.5_cp*xi)

        case (PROFILE_BICKLEY)
            amplify = 1.0_cp/(cosh(0.5_cp*xi))**2.0_cp

        case (PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_SURFACE)
            amplify = exp(-0.5_cp*xi**2.0_cp)

        case (PROFILE_GAUSSIAN_SYM)
            amplify = exp(-0.5_cp*(xi - 0.5_cp*var%diam/var%thick)**2.0_cp) &
                      + exp(-0.5_cp*(xi + 0.5_cp*var%diam/var%thick)**2.0_cp)

        case (PROFILE_GAUSSIAN_ANTISYM)
            amplify = exp(-0.5_cp*(xi - 0.5_cp*var%diam/var%thick)**2.0_cp) &
                      - exp(-0.5_cp*(xi + 0.5_cp*var%diam/var%thick)**2.0_cp)

        case (PROFILE_EKMAN_U)
            amplify = 1.0_cp - exp(-xi)*cos(xi)

        case (PROFILE_EKMAN_U_P)
            amplify = 1.0_cp - exp(-xi)*cos(xi) ! + perturbation:

            cnought = pi_cp*pi_cp/4.0_cp/4.0_cp       ! Maximum initial Perturbation is at y=pi/2*var%thick
            zamp = sqrt(2.0_cp)*xi*exp(-xi*xi/8.0_cp/cnought)/(var%thick*var%thick*4.0_cp*cnought)**1.5_cp
            amplify = amplify + zamp                  ! Add Perturbations

        case (PROFILE_EKMAN_V)
            amplify = -exp(-xi)*sin(xi)

        end select

    end if

    ! var%mean profile plus two linear-varying layers
    f = var%mean + var%delta*amplify &
        + var%lslope*yrel*0.5_cp*(1.0_cp - sign(1.0_cp, yrel)) &
        + var%uslope*yrel*0.5_cp*(1.0_cp + sign(1.0_cp, yrel))

    ! -------------------------------------------------------------------
    ! special profiles
    ! -------------------------------------------------------------------
    select case (var%type)

    case (PROFILE_LINEAR_CROP)
        if (yrel < 0.0_cp) then
            f = min(var%lslope*yrel, var%lslope*var%thick)
        else
            f = max(var%uslope*yrel, var%uslope*var%thick)
        end if

    case (PROFILE_MIXEDLAYER)
        if (yrel < 0.0_cp) then
            f = min(var%lslope*yrel, var%lslope*var%thick)
        else
            f = max(var%uslope*yrel, var%uslope*var%thick)
        end if
        f = f - 0.25_cp*var%uslope*var%thick*(1.0_cp - sign(1.0_cp, y - var%thick))

    case (PROFILE_ERF_SURFACE)
        xi = y/var%parameters(3)
        f = f + var%parameters(4)*0.5_cp*(1.0_cp + erf(-0.5_cp*xi))

    end select

    return
end function PROFILES

subroutine PROFILES_DERTOTHICK(derivative, var)  ! Obtain thick from the value of the maximum derivative
    use TLAB_TYPES, only: profiles_dt, cp
    use TLAB_CONSTANTS
    use TLAB_PROCS
    implicit none

    real(cp), intent(in) :: derivative
    type(profiles_dt), intent(inout) :: var

    real(cp) thick_ratio    ! for readibility

    select case (var%type)

    case (PROFILE_TANH, PROFILE_TANH_SYM, PROFILE_TANH_ANTISYM)
        thick_ratio = 4.0_cp
        var%thick = -var%delta/derivative/thick_ratio

    case (PROFILE_ERF, PROFILE_ERF_ANTISYM)
        thick_ratio = 2.0_cp*sqrt(pi_cp)
        var%thick = -var%delta/(derivative - var%uslope)/thick_ratio

    case (PROFILE_ERF_SURFACE)
        thick_ratio = 2.0_cp*sqrt(pi_cp)
        var%parameters(3) = -var%parameters(4)/derivative/thick_ratio

    case default
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Undevelop derivative input for this case.')
        call TLAB_STOP(DNS_ERROR_UNDEVELOP)

    end select

    return
end subroutine PROFILES_DERTOTHICK
