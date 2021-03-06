#include "types.h"
#include "dns_const.h"

PROGRAM VBURGERS

  USE DNS_CONSTANTS
  USE DNS_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(:,:),   ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: a, b, c
  TREAL, DIMENSION(:,:),   ALLOCATABLE :: wrk1d, wrk2d
  TREAL, DIMENSION(:),     ALLOCATABLE :: wrk3d, tmp1

  TINTEGER i, j, k,  bcs(2,2)
  TREAL dummy, error

! ###################################################################
  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL('dns.ini')
#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

  isize_wrk3d = isize_txc_field

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d,inb_wrk1d+1))
  ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))
  ALLOCATE(a(imax,jmax,kmax),b(imax,jmax,kmax),c(imax,jmax,kmax))
  ALLOCATE(tmp1(isize_txc_field),wrk3d(isize_wrk3d))

#include "dns_read_grid.h"

  CALL FI_PROFILES_INITIALIZE(wrk1d)

  bcs = 0

! ###################################################################
! Define forcing term
! ###################################################################
  CALL DNS_READ_FIELDS('field.inp', i1, imax,jmax,kmax, i1,i0, isize_wrk3d, a, wrk3d)

! ###################################################################
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), a,b, c, wrk2d,wrk3d)
  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
!           b(i,j,k) = b(i,j,k) *visc - a(i,j,k) *c(i,j,k)
           b(i,j,k) = b(i,j,k) *visc *ribackground(j)- a(i,j,k) *c(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  CALL OPR_BURGERS_X(i0,i0, imax,jmax,kmax, bcs, g(1), a,a,a, c, tmp1, wrk2d,wrk3d)
  c = c -b

  error = C_0_R
  dummy = C_0_R
  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
           error = error + c(i,j,k)*c(i,j,k)
           dummy = dummy + b(i,j,k)*b(i,j,k)
        ENDDO
     ENDDO
  ENDDO
  WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
!  CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, isize_wrk3d, e, wrk3d)

! ###################################################################
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), a,b, c, wrk2d,wrk3d)
  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
!           b(i,j,k) = b(i,j,k) *visc - a(i,j,k) *c(i,j,k)
           b(i,j,k) = b(i,j,k) *visc *ribackground(j)- a(i,j,k) *c(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  CALL OPR_BURGERS_Y(i0,i0, imax,jmax,kmax, bcs, g(2), a,a,a, c, tmp1, wrk2d,wrk3d)
  c = c -b

  error = C_0_R
  dummy = C_0_R
  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
           error = error + c(i,j,k)*c(i,j,k)
           dummy = dummy + b(i,j,k)*b(i,j,k)
        ENDDO
     ENDDO
  ENDDO
  WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
!  CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, isize_wrk3d, c, wrk3d)

! ###################################################################
  IF ( g(3)%size .GT. 1 ) THEN

  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), a,b, c, wrk2d,wrk3d)
  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
!           b(i,j,k) = b(i,j,k) *visc - a(i,j,k) *c(i,j,k)
           b(i,j,k) = b(i,j,k) *visc *ribackground(j)- a(i,j,k) *c(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  CALL OPR_BURGERS_Z(i0,i0, imax,jmax,kmax, bcs, g(3), a,a,a, c, tmp1, wrk2d,wrk3d)
  c = c -b

  error = C_0_R
  dummy = C_0_R
  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
           error = error + c(i,j,k)*c(i,j,k)
           dummy = dummy + b(i,j,k)*b(i,j,k)
        ENDDO
     ENDDO
  ENDDO
  WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
!  CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, isize_wrk3d, e, wrk3d)

  END IF

  CALL DNS_STOP(0)
END PROGRAM VBURGERS
