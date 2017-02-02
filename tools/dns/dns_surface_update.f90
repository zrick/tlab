#include "types.h"
#include "dns_const.h"
! SUBROUTINE SURFACE UPDATE
! Calculates and provides interactive surface boundary condition 
SUBROUTINE DNS_SURFACE_UPDATE(is,bcs_hb,q,hq,s,hs,tmp1,tmp2,aux,wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax,inb_flow,inb_vars,inb_scal,j1bc 
  USE DNS_GLOBAL, ONLY : isize_field,isize_wrk1d 
  USE DNS_GLOBAL, ONLY : visc,itime,schmidt
  USE DNS_GLOBAL, ONLY : imode_fdm 
  USE DNS_LOCAL,  ONLY : dtime

  IMPLICIT NONE  

#include "integers.h"

  TINTEGER is 
  TREAL, DIMENSION(imax,kmax,inb_vars) :: bcs_hb 
  TREAL, DIMENSION(isize_field,*)      :: q,hq,s,hs 
  TREAL, DIMENSION(isize_field)        :: tmp1,tmp2
  TREAL, DIMENSION(imax,kmax,6),TARGET :: aux 
  TREAL, DIMENSION(isize_wrk1d,*)      :: wrk1d
  TREAL, DIMENSION(*)                  :: wrk2d,wrk3d

  TINTEGER nxy,ip,k 
  TREAL dx(1), dy(1), dz(1)   ! To use uld wrappers for derivatives  

  TREAL, DIMENSION(:,:), POINTER       :: hfx,hfx_anom
  TREAL :: diff,hfx_avg,l_couple 

  TREAL AVG1V2D

  diff = visc/schmidt(is)
  nxy = imax*jmax   

  hfx =>      aux(:,:,1) 
  hfx_anom => aux(:,:,2) 


  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, s(:,is), tmp1, i0,i0, wrk1d,wrk2d,wrk3d) 

  ip=1
  DO k=1,kmax
     hfx(:,k) = diff*tmp1(ip:ip+imax-1); ip=ip+nxy
  ENDDO

  hfx_avg = diff*AVG1V2D(imax,jmax,kmax,1,1,tmp1) 
  hfx_anom = hfx - hfx_avg
  bcs_hb(:,:,inb_flow+is) = bcs_hb(:,:,inb_flow+is) + l_couple*hfx_anom

  WRITE(*,*) 'Heat flux:', hfx_avg, AVG1V2D(imax,1,kmax,1,2,hfx)-hfx_avg**2,AVG1V2D(imax,jmax,kmax,1,2,s)-1

  ! TESTING: solve ds/dt=s => s(t) = exp(t) on the surface boundary
  ! ip=1 
  ! DO k=1,kmax  
  !    bcs_hb(:,k,inb_flow+is) = bcs_hb(:,k,inb_flow+is) + s(ip:ip+imax-1,is)  
  !    ip = ip+nxy
  ! ENDDO

END SUBROUTINE DNS_SURFACE_UPDATE
