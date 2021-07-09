#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY / ATHORS
!#
!# 2021/XX/XX - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   called once in dns_main.f90 l.210 to allocate needed memory
!#    
!#    
!#    
!#
!# 
!########################################################################
!# ARGUMENTS 
!#
!# 
!#                            
!#                           
!#
!########################################################################
!# REQUIREMENTS
!#
!# 
!#                            
!#                           
!#
!########################################################################

subroutine ALLOCATE_IBM(allocated)

  use DNS_IBM
  use DNS_CONSTANTS, only: lfile, efile
  use DNS_GLOBAL,    only: imax, jmax, kmax 
  use DNS_GLOBAL,    only: isize_field
  use DNS_GLOBAL,    only: xbars_geo

#ifdef USE_MPI
  use DNS_MPI,       only: ims_size_i, ims_size_j, ims_size_k 
#endif    

  implicit none

#include "integers.h"

#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"
  TINTEGER, parameter       :: idi = DNS_MPI_I_PARTIAL 
  TINTEGER, parameter       :: idj = DNS_MPI_J_PARTIAL 
  TINTEGER, parameter       :: idk = DNS_MPI_K_PARTIAL 
#endif

  logical, intent(inout)    :: allocated       ! flag, just allocate memory space once
  TINTEGER                  :: ierr, inb_ibm
  TINTEGER                  :: nyz, nxz, nxy
  TINTEGER                  :: nob_max
  character(128)            :: str, line

  ! ================================================================== !

  inb_ibm = i1 ! can be also defined in dns_read_local.f90 and module dns_global.f90 (cf. inb_flow ...)

  ! max(nobi_max,nobj_max,nobk_max) from dns.ini file
  nob_max = xbars_geo(1)

  ! npages (cf. dns_mpi_initialize.f90)
#ifdef USE_MPI 
  nyz = ims_size_i(idi) ! local
  nxz = ims_size_j(idj) 
  nxy = ims_size_k(idk) 
#else
  nyz = jmax * kmax     ! global
  nxz = imax * kmax     
  nxy = imax * jmax     
#endif      

  ! ================================================================== !

  ! allocate here all ibm related arrays
  if ( .not. allocated ) then 

    ! eps_aux 3d-array (debugging / generating geometry)
    write(str,*) inb_ibm; line = 'Allocating array IBM eps_aux of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(eps_aux(imax,jmax,kmax), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for eps_aux.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! eps
    write(str,*) inb_ibm; line = 'Allocating array IBM eps of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(eps(isize_field), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for eps.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! ------------------------------------------------------------------ !

    ! epsi
    write(str,*) inb_ibm; line = 'Allocating array IBM epsi of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(epsi(isize_field), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for epsi.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! epsj
    write(str,*) inb_ibm; line = 'Allocating array IBM epsj of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(epsj(isize_field), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for epsj.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! epsk
    write(str,*) inb_ibm; line = 'Allocating array IBM epsk of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(epsk(isize_field), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for epsk.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! ------------------------------------------------------------------ !

    ! nobi
    write(str,*) inb_ibm; line = 'Allocating array IBM nobi of size '//trim(adjustl(str))//'x'
    write(str,*) nyz; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobi(nyz), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobi.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobj
    write(str,*) inb_ibm; line = 'Allocating array IBM nobj of size '//trim(adjustl(str))//'x'
    write(str,*) nxz; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobj(nxz), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobj.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobk
    write(str,*) inb_ibm; line = 'Allocating array IBM nobk of size '//trim(adjustl(str))//'x'
    write(str,*) nxy; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobk(nxy), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobk.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! ------------------------------------------------------------------ !

    ! nobi_b
    write(str,*) inb_ibm; line = 'Allocating array IBM nobi_b of size '//trim(adjustl(str))//'x'
    write(str,*) nyz*nob_max; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobi_b(nyz*nob_max), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobi_b.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobj_b
    write(str,*) inb_ibm; line = 'Allocating array IBM nobj_b of size '//trim(adjustl(str))//'x'
    write(str,*) nxz*nob_max; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobj_b(nxz*nob_max), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobj_b.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobk_b
    write(str,*) inb_ibm; line = 'Allocating array IBM nobk_b of size '//trim(adjustl(str))//'x'
    write(str,*) nxy*nob_max; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobk_b(nxy*nob_max), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobk_b.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! ------------------------------------------------------------------ !

    ! nobi_e
    write(str,*) inb_ibm; line = 'Allocating array IBM nobi_e of size '//trim(adjustl(str))//'x'
    write(str,*) nyz*nob_max; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobi_e(nyz*nob_max), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobi_e.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobj_e
    write(str,*) inb_ibm; line = 'Allocating array IBM nobj_e of size '//trim(adjustl(str))//'x'
    write(str,*) nxz*nob_max; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobj_e(nxz*nob_max), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobj_e.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobk_e
    write(str,*) inb_ibm; line = 'Allocating array IBM nobk_e of size '//trim(adjustl(str))//'x'
    write(str,*) nxy*nob_max; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobk_e(nxy*nob_max), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobk_e.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! ------------------------------------------------------------------ !

    ! set alloc flag: done
    allocated = .true.

  end if

  return
end subroutine ALLOCATE_IBM

!########################################################################