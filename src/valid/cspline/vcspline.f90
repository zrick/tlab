#include "types.h"
#include "dns_const.h"

!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2022/03/16 - J. Kostelecky
!#              Created           
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate cubic splines.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################

program CSPLINE

  use TLAB_TYPES, only : grid_dt
  use TLAB_PROCS

  implicit none
 
#include "integers.h"
  
  ! define spline parameters here 
  type(grid_dt)                          :: g, g_int                            ! original and refined grid
  TINTEGER, parameter                    :: inb_grid = 57           
  TINTEGER, parameter                    :: imax     = 11                       ! number of data points
  TINTEGER, parameter                    :: mesh     = 10                       ! mesh refinement factor (mesh=1 for x_int=x)
  TINTEGER, parameter                    :: imax_int = (imax+(mesh-1)*(imax-1)) ! number of spline data points
  TINTEGER                               :: i, j
  TREAL                                  :: xb, xe, lambda
  TREAL                                  :: res_2, res_inf
  logical                                :: periodic
  ! time messurement
  TINTEGER                               :: iter, c1, c2, c3, c11, c22, c33   
  TINTEGER, parameter                    :: start=0, end=50000, step=5000
  TREAL,    dimension((end-start)/step+1):: delta_t
  ! data arrays
  TREAL,    dimension(imax,    inb_grid) :: x
  TREAL,    dimension(imax_int,inb_grid) :: x_int                                                
  TREAL,    dimension(imax             ) :: y
  TREAL,    dimension(imax_int         ) :: y_sp, y_int, delta
  TREAL,    dimension(imax_int         ) :: dydx,ddydx
  ! boundary conditions
  TINTEGER, dimension(2     )            :: bc    ! boundary condition for both endpoints
  TREAL,    dimension(2     )            :: bcval ! boundary values of 1st or 2nd deriv. at endpoints 
  ! working arrays
  TREAL,    dimension(9*imax)            :: wrk  
  TREAL,    DIMENSION(imax,    5)        :: wrk1d
  TREAL,    DIMENSION(imax_int,5)        :: wrk1d_int

! ###################################################################
! set periodicity here
  periodic = .FALSE. ! periodic case not implemented yet

! choose boundary conditions for splines
  bc(1)    = CS_BCS_FIXED_2!CS_BCS_FIXED_1 !CS_BCS_NATURAL
  bc(2)    = bc(1)
  bcval(:) = C_2_R

! initialize original grid with test function
  g%periodic = periodic
  g%size     = imax 
  g%scale    = C_1_R
  lambda     = C_1_R
  g%uniform  = .TRUE.
  g%mode_fdm = FDM_COM6_JACOBIAN 
  do i = 1,imax
    x(i,1) = M_REAL(i-1)/M_REAL(imax-1)*g%scale
    y(i  ) = sin(C_2_R*C_PI_R/g%scale*lambda*x(i,1))
  end do

! set interval for spline approximation
  xb = x(1,1);  xe = x(imax,1)

! create equally spaced array to evaluate spline with boundaries [xb,xe]
  g_int%periodic = g%periodic; g_int%size = imax_int; g_int%scale = g%scale
  g_int%uniform  = g%uniform;  g_int%mode_fdm = g%mode_fdm 
  do i = 1,imax_int
    x_int(i,1) = xb + (xe - xb) * (i - C_1_R) / (imax_int - C_1_R)
    y_int(i  ) = sin(C_2_R*C_PI_R/g_int%scale*lambda*x_int(i,1))
  end do

! initialize grids for fdm calls
  call FDM_INITIALIZE(x,     g,     wrk1d    ) 
  call FDM_INITIALIZE(x_int, g_int, wrk1d_int) 

! cubic spline function
  call CUBIC_SPLINE(bc, bcval, imax, imax_int, g%nodes, y, g_int%nodes, y_sp, wrk)

! compute first and second derivative at boundary points
  if ( (g%mode_fdm == FDM_COM6_JACOBIAN) .and. (.not. periodic)) then
  ! first derivative
    call FDM_C1N6_LHS(imax_int,     i0, i0, g_int%jac, wrk1d_int(1,1),wrk1d_int(1,2),wrk1d_int(1,3))
    call FDM_C1N6_RHS(imax_int, i1, i0, i0, y_sp, dydx)
    call TRIDFS(imax_int,    wrk1d_int(1,1),wrk1d_int(1,2),wrk1d_int(1,3))
    call TRIDSS(imax_int,i1, wrk1d_int(1,1),wrk1d_int(1,2),wrk1d_int(1,3), dydx)
  ! second derivative
    call FDM_C2N6H_LHS(imax_int,     i0, i0, g_int%jac, wrk1d_int(1,1),wrk1d_int(1,2),wrk1d_int(1,3))
    call FDM_C2N6H_RHS(imax_int, i1, i0, i0, y_sp, ddydx)
    call TRIDFS(imax_int,    wrk1d_int(1,1),wrk1d_int(1,2),wrk1d_int(1,3))
    call TRIDSS(imax_int,i1, wrk1d_int(1,1),wrk1d_int(1,2),wrk1d_int(1,3), ddydx)
    write(*,'(1X,A,ES9.2,3X,ES9.2)')'1st derivative values at boundary points = ',  dydx(1),  dydx(imax_int)
    write(*,'(1X,A,ES9.2,3X,ES9.2)')'2nd derivative values at boundary points = ', ddydx(1), ddydx(imax_int)
    call FDM_C1N6_LHS(imax_int,     i0, i0, g_int%jac, wrk1d_int(1,1),wrk1d_int(1,2),wrk1d_int(1,3))
    call FDM_C1N6_RHS(imax_int, i1, i0, i0, dydx, ddydx)
    call TRIDFS(imax_int,    wrk1d_int(1,1),wrk1d_int(1,2),wrk1d_int(1,3))
    call TRIDSS(imax_int,i1, wrk1d_int(1,1),wrk1d_int(1,2),wrk1d_int(1,3), ddydx)
    write(*,'(1X,A,ES9.2,3X,ES9.2)')'2nd derivative values at boundary points = ', ddydx(1), ddydx(imax_int)
  end if 

! ###################################################################
  write(*,*) '================================================================================'
  write(*,*) '================  Validation routines for cubic splines ========================'
  write(*,*) '================================================================================'
! ###################################################################
! ! 0. Validation - scaling of spline functions
!   write(*,*) '0. Validation routine: Messure time consumption'
!   j = 1
!   do iter = start, end, step
!     call system_clock(c1,c2,c3)
!     i = 1
!     do while (i<=iter)
!       i = i + 1
!         call CUBIC_SPLINE(imax, imax_int, x, y, x_int, y_sp, wrk)
!     end do
!     ! time consumption + save
!     call system_clock(c11,c22,c33)
!     write(*,30) i-1, c11 - c1 
!     30 format(1X, I5,' iterations with delta_t (ms) : ',I4) 
!     delta_t(j) = c11 - c1
!     j = j + 1
!   end do

! 1. Validation
  write(*,*) '================================================================================' 
  write(*,*) '1. Validation routine: Check if data points coincident with spline points'
  do i = 1,imax
    if ((abs(y(i) - y_sp(1 + (i-1)*mesh))) <= 1e-10) then
      write(*,40) i, g%nodes(i), y(i), (1 + (i-1)*mesh), abs(y(i) - y_sp(1 + (i-1)*mesh))
      40 format(1X,'Data point ', I3, ' with ', '(', F4.2, ',', F4.2,')', ' is on spline point ', I3,', residua = ' ,ES9.2)
    else
      write(*,50) i, g%nodes(i), y(i), (1 + (i-1)*mesh), g_int%nodes(1 + (i-1)*mesh), y_sp(1 + (i-1)*mesh), abs(y(i) - y_sp(1 + (i-1)*mesh))
      50 format(1X,'Data point ', I3, ' with ', '(', F4.2, ',', F4.2,')', ' is not on spline point ', I3, ' with ', '(', F4.2, ',', F4.2,')',', residua = ' ,ES9.2)
    end if 
  end do

! 2. Validation
  write(*,*) '================================================================================' 
  write(*,*) '2. Validation routine: Error norms of exact solution and spline points'
  do i = 1 , imax_int
    delta(i) = y_int(i) - y_sp(i)
  end do  
  res_2   = norm2(delta)
  res_inf = maxval(abs(delta))
  write(*,'(1X,A,ES9.2)') 'L_2   error norm          :', res_2
  write(*,'(1X,A,ES9.2)') 'L_inf error norm          :', res_inf
  write(*,*) '================================================================================'

! IO - error and function values
  open(20, file='cspline.dat')
  do i = 1,imax_int
    write(20,1000) g_int%nodes(i), y_int(i), y_sp(i), y_int(i) - y_sp(i) 
  end do
  close(20)
  1000 format(6(1x,e16.10))

  ! Validation option in pyhthon with:
  !   from   scipy.interpolate import CubicSpline  
  !   sp = CubicSpline(xorg,yorg,bc_type='natural')

  stop

end program CSPLINE