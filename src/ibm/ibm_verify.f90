#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI 
#include "dns_const_mpi.h"  
#endif

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/05 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   verify that the geometry used satisfies all restrictions
!#    
!# 
!########################################################################
!# ARGUMENTS 
!#
!#
!########################################################################
!# REQUIREMENTS
!#
!# 
!########################################################################

subroutine IBM_VERIFY_GEOMETRY()

  use IBM_VARS 
  use TLAB_PROCS
  use TLAB_VARS,      only : g, icalc_scal
  use TLAB_VARS,      only : imax, jmax, kmax
  use TLAB_CONSTANTS, only : efile
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_VARS,  only : ims_size_i, ims_size_k, ims_err
  use TLAB_MPI_VARS,  only : ims_npro_i, ims_npro_k
#ifdef IBM_DEBUG
  use TLAB_MPI_VARS,  only : ims_pro
#endif
#endif    
   
  implicit none

#include "integers.h"

#ifdef USE_MPI 
  TINTEGER, parameter :: idi = TLAB_MPI_I_PARTIAL 
  TINTEGER, parameter :: idk = TLAB_MPI_K_PARTIAL 
  TREAL               :: dummy
#else
#ifdef IBM_DEBUG
  TINTEGER, parameter :: ims_pro = 0 
#endif
#endif
  TINTEGER            :: nyz, nxz, nxy, dummi    
  TREAL               :: ob_min
  ! ================================================================== !

#ifdef IBM_DEBUG
  if ( ims_pro == 0 ) write(*,*) '================ Verifying the geometry ================='
#endif

  ! npages
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    nyz = ims_size_i(idi)
  else
#endif
  nyz = jmax * kmax 
#ifdef USE_MPI
  end if
#endif

  nxz = imax * kmax     

#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    nxy = ims_size_k(idk)
  else
#endif
  nxy = imax * jmax
#ifdef USE_MPI
  end if
#endif

  ! check if "MaxNumberObj" in dns.ini is set correctly
  dummi = maxval((/nobi_max, nobj_max, nobk_max/))
  if ( dummi > nob_max ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY. MaxNumberObj too small in dns.ini.')
    call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
  end if
  
  ! check if any objects are present
  ob_min = sum(eps)
#ifdef USE_MPI
  dummy = ob_min
  call MPI_ALLREDUCE(dummy, ob_min, i1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
  if ( ob_min == 0 ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY. No objects in flow.')
    call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
  endif

  ! check all objects in each direction
  call IBM_VERIFY(g(1), nyz, isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e) ! x
  call IBM_VERIFY(g(2), nxz, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e) ! y
  call IBM_VERIFY(g(3), nxy, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e) ! z

  ! check if objects on upper boundary are present
  call IBM_VERIFY_UP(eps)

  ! if so, check with scalar options
  if ( icalc_scal == 1 .and. ibm_objup ) then
    call IBM_VERIFY_SCAL()
  end if 
  
  ! final
#ifdef IBM_DEBUG
  if ( ims_pro == 0 ) write(*,*) 'Verification successfull!'
#endif

  return
end subroutine IBM_VERIFY_GEOMETRY

!########################################################################

subroutine IBM_VERIFY(g, nlines, isize_nob, isize_nob_be, nob, nob_b, nob_e)

  use IBM_VARS,       only : nflu
  use TLAB_CONSTANTS, only : efile
  use TLAB_TYPES,     only : grid_dt
  use TLAB_PROCS

  implicit none
  
#include "integers.h"

  type(grid_dt),                     intent(in) :: g
  TINTEGER,                          intent(in) :: nlines, isize_nob, isize_nob_be
  TINTEGER, dimension(isize_nob),    intent(in) :: nob
  TINTEGER, dimension(isize_nob_be), intent(in) :: nob_b, nob_e

  TINTEGER                                      :: ii, ip, ipp, iob
  TINTEGER                                      :: fp_min, fp_l, fp_r, fp_inter
  TINTEGER                                      :: sp_min, sp_ob
  
  ! ================================================================== !
  
  ! min vals solid/fluid
  fp_min = nflu ! 3rd point can be next interface/boundary
  sp_min = i3   ! 3 solid + 2 interface points - 1

  ! index ii (dummy index; for x,y,z: ii == jk,ik,ij)
  do ii = 1, nlines          ! index of ii-plane, loop over plane and check for objects in each line
    if ( nob(ii) /= 0 ) then ! if line contains immersed object(s) --yes-->  spline interpolation
      ip = i0        
      do iob = 1, nob(ii)    ! loop over immersed object(s)
        ! ================================================================== !
        ! check number of fluid points (to borders and between objects)
        if ( iob == 1 ) then                    ! left
          fp_l = nob_b(ip+ii) - i1
          if ( fp_l < fp_min .and. fp_l /= 0 ) then
            call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY. Not enough fluid points between left border and first object. Check also MaxNumberObj in dns.ini.')
            call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
          end if
        end if 
        if ( iob < nob(ii) ) then               ! in between objects
          ipp = ip + nlines                     ! next obj
          fp_inter = nob_b(ipp+ii) - nob_e(ip+ii) - i1
          if ( fp_inter < fp_min ) then
            call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY. Not enough fluid points between objects. Check also MaxNumberObj in dns.ini.')
            call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
          end if
        else if ( iob == nob(ii) ) then         ! right
          fp_r = g%size - nob_e(ip+ii)
          if ( fp_r < fp_min .and. fp_r /= 0 ) then
            call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY. Not enough fluid points between right border and first object. Check also MaxNumberObj in dns.ini.')
            call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
          end if
        end if
        ! ================================================================== !
        ! check number of solid points of object
        sp_ob = nob_e(ip+ii) - nob_b(ip+ii) + i1
        if ( sp_ob < sp_min ) then
          call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY. Not enough solid points in objects. Check also MaxNumberObj in dns.ini.')
          call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
        end if
        ! ================================================================== !        
        ip = ip + nlines
      end do
      ! check for overlapping objects in periodic case
      if ( g%periodic ) then
        if ( nob(ii) > 1 ) then
          fp_l = nob_b(ii)                           ! begin of first object
          fp_r = nob_e((nob(ii) - i1) * nlines + ii) ! end of last object
          if ( (fp_r - fp_l + i1) == g%size ) then
            call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY. Overlapping objects are not allowed. Check also MaxNumberObj in dns.ini.')
            call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
          end if 
        end if
      end if
    end if
  end do

  return
end subroutine IBM_VERIFY

!########################################################################

subroutine IBM_VERIFY_UP(eps)

  use IBM_VARS,       only : ibm_objup, max_height_objlo, max_height_objup
  use TLAB_VARS,      only : isize_field, imax,jmax,kmax
  use TLAB_CONSTANTS, only : efile
  use TLAB_PROCS
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_VARS,  only : ims_err
#endif    

  implicit none
  
#include "integers.h"

  TREAL, dimension(isize_field), intent(in) :: eps

  TINTEGER                                  :: ip_t, ip_b, j, k, nxy
  TREAL                                     :: dummy, top

  ! ================================================================== !

  ! check if objects are present on upper boundary 
  nxy  = imax*jmax
  ip_t = imax*(jmax-1) + 1
  top  = 0
  do k = 1, kmax
    top = top + sum(eps(ip_t:ip_t+imax-1)) 
    ip_t = ip_t + nxy
  end do

#ifdef USE_MPI
  dummy = top
  call MPI_ALLREDUCE(dummy, top, i1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  if ( top > 0 ) then
    ibm_objup = .true. ! objects are present on upper boundary
  elseif ( top == 0 ) then
    ibm_objup = .false.
  else
    call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY. Problem with objects on upper domain boundary.')
    call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
  endif

  ! if so, check height of upper and lower objects
  if ( ibm_objup ) then
    ! upper 
    do j = 2, jmax
      ip_t = imax*(jmax-j) + 1
      dummy = 0
      do k = 1, kmax
        dummy = dummy + sum(eps(ip_t:ip_t+imax-1)) 
        ip_t = ip_t + nxy
      end do  
      if ( dummy == 0 ) exit
      max_height_objup = j   
    end do
#ifdef USE_MPI
    dummy = max_height_objup
    call MPI_ALLREDUCE(dummy, max_height_objup, i1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif
  end if
  ! lower (always assuming objects on lower boundary!)
  do j = 1, jmax
    ip_b = imax*(j-1)+1
    dummy = 0
    do k = 1, kmax
      dummy = dummy + sum(eps(ip_b:ip_b+imax-1)) 
      ip_b = ip_b + nxy
    end do  
    if ( dummy == 0 ) exit
    max_height_objlo = j   
  end do
#ifdef USE_MPI
  dummy = max_height_objlo
  call MPI_ALLREDUCE(dummy, max_height_objlo, i1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif  

  return
end subroutine IBM_VERIFY_UP

!########################################################################

subroutine IBM_VERIFY_SCAL()

  use IBM_VARS,       only : max_height_objlo, max_height_objup, imode_ibm_scal
  use TLAB_VARS,      only : g
  use TLAB_CONSTANTS, only : efile
  use TLAB_PROCS
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_VARS,  only : ims_err
#endif    

  implicit none

  ! ================================================================== !

  if ( imode_ibm_scal == 1 ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY. No objects on upper domain allowed if IBM is turned on for scalars.')
    call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
  end if

  ! no overlapping in vertical direction allowed
  if ( (g(2)%size-max_height_objup) < max_height_objlo ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY. No overlapping of objects in vertical direction if scalars is on.')
    call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
  end if 

  return
end subroutine IBM_VERIFY_SCAL