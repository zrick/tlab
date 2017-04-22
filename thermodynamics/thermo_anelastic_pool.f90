#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY
!#
!# 2017/03/22 - J.P. Mellado
!#              Extracted from airwater_pool
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate thermodynamic properties from h, p and composition in the
!# incompressible formulation, when the thermodynamic pressure is
!# a given profile
!#
!# s1 is specific static energy
!# s2 is total water specific humidity
!# s3 is liquid water specific humidity, if any
!#
!########################################################################

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_TEMPERATURE(nx,ny,nz, s, e, T)

  USE THERMO_GLOBAL, ONLY : imixture, THERMO_AI

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: T

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL E_LOC

  TREAL Cd, Cdv, Lv0, Cvl

! ###################################################################
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  IF      ( imixture .EQ. 0 ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T(ij) = s(ij,1) - E_LOC
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T(ij) = (s(ij,1) - E_LOC ) / ( Cd + s(ij,2) *Cdv )
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T(ij) = (s(ij,1) - E_LOC + s(ij,3)*Lv0 )  / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )

        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_TEMPERATURE

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_DENSITY(nx,ny,nz, s, e,p, rho)

  USE THERMO_GLOBAL, ONLY : imixture, WGHT_INV, THERMO_AI, MRATIO

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: rho

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL WMEAN_INV, P_LOC, E_LOC, T_LOC
  
  TREAL Rd, Rdv, Cd, Cdv, Lv0, Cvl

! ###################################################################
  Rd = WGHT_INV(2)
  Rdv= WGHT_INV(1) - WGHT_INV(2)
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  IF      ( imixture .EQ. 0 ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = s(ij,1) - E_LOC
           rho(ij) = P_LOC /T_LOC
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC ) / ( Cd + s(ij,2) *Cdv )
           WMEAN_INV = Rd + s(ij,2) *Rdv
           rho(ij) = P_LOC /( WMEAN_INV *T_LOC )
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC + s(ij,3)*Lv0 )  / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )
           WMEAN_INV = Rd + s(ij,2) *Rdv - s(ij,3)*WGHT_INV(1)
           rho(ij) = P_LOC /(WMEAN_INV*T_LOC)
           
        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_DENSITY

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_BUOYANCY(nx,ny,nz, s, e,p,r, b)

  USE THERMO_GLOBAL, ONLY : imixture, WGHT_INV, THERMO_AI, MRATIO
  
  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p,r
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: b

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL WMEAN_INV, P_LOC, E_LOC, R_LOC, R_LOC_INV, T_LOC
  
  TREAL Rd, Rdv, Cd, Cdv, Lv0, Cvl

! ###################################################################
  Rd = WGHT_INV(2)
  Rdv= WGHT_INV(1) - WGHT_INV(2)
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  IF      ( imixture .EQ. 0 ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        R_LOC = r(is)
        R_LOC_INV = C_1_R /R_LOC
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = s(ij,1) - E_LOC
           b(ij) = R_LOC_INV *(R_LOC -P_LOC /T_LOC )
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        R_LOC = r(is)
        R_LOC_INV = C_1_R /R_LOC
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC ) / ( Cd + s(ij,2) *Cdv )
           WMEAN_INV = Rd + s(ij,2) *Rdv
           b(ij) = R_LOC_INV *( R_LOC -P_LOC /(WMEAN_INV*T_LOC) )
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        R_LOC = r(is)
        R_LOC_INV = C_1_R /R_LOC
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC + s(ij,3) *Lv0 ) / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )
           WMEAN_INV = Rd + s(ij,2) *Rdv - s(ij,3)*WGHT_INV(1)
           b(ij) = R_LOC_INV *( R_LOC -P_LOC /(WMEAN_INV*T_LOC) )
           
        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_BUOYANCY

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_QVEQU(nx,ny,nz, s, e,p, T,qvequ)

  USE THERMO_GLOBAL, ONLY : imixture, THERMO_AI, THERMO_PSAT, NPSAT, WGHT_INV, MRATIO

  IMPLICIT NONE

#include "integers.h"

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p 
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: qvequ, T
  
! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is, ipsat
  TREAL psat, E_LOC, P_LOC
  
  TREAL Cd, Cdv, Lv0, Cvl, rd_ov_rv
  
! ###################################################################
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  rd_ov_rv = WGHT_INV(2) /WGHT_INV(1)
  
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO*p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T(ij) = (s(ij,1) - E_LOC + s(ij,3)*Lv0 )  / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )
           
           psat = C_0_R
           DO ipsat = NPSAT,1,-1
              psat = psat*T(ij) + THERMO_PSAT(ipsat)
           ENDDO
           qvequ(ij) = rd_ov_rv *( C_1_R -s(ij,2) ) /( P_LOC/psat -C_1_R )
           
        ENDDO

     ENDDO
  ENDIF
     
  RETURN
END SUBROUTINE THERMO_ANELASTIC_QVEQU

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_RELATIVEHUMIDITY(nx,ny,nz, s, e,p, T,rh)

  USE THERMO_GLOBAL, ONLY : imixture, THERMO_AI, THERMO_PSAT, NPSAT, WGHT_INV, MRATIO

  IMPLICIT NONE

#include "integers.h"

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p 
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: rh, T
  
! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is, ipsat
  TREAL psat, E_LOC, P_LOC
  
  TREAL Rv, Rd, Rdv, Cd, Cdv, Lv0, Cvl
  
! ###################################################################
  Rv = WGHT_INV(1)
  Rd = WGHT_INV(2)
  Rdv= WGHT_INV(1) - WGHT_INV(2)
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO*p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T(ij) = (s(ij,1) - E_LOC + s(ij,3)*Lv0 )  / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )
           
           psat = C_0_R
           DO ipsat = NPSAT,1,-1
              psat = psat*T(ij) + THERMO_PSAT(ipsat)
           ENDDO
           rh(ij) = (s(ij,2) -s(ij,3)) *P_LOC/psat *Rv /( Rd + s(ij,2) *Rdv - s(ij,3) *Rv )
           rh(ij) = rh(ij) *C_100_R
           
        ENDDO

     ENDDO
  ENDIF
     
  RETURN
END SUBROUTINE THERMO_ANELASTIC_RELATIVEHUMIDITY

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_THETA(nx,ny,nz, s, e,p, theta)

  USE THERMO_GLOBAL, ONLY : imixture, WGHT_INV, THERMO_AI, MRATIO, GRATIO

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: theta

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL P_LOC, E_LOC, T_LOC

  TREAL Rd, Rdv, Rl, Cd, Cdv, Cl, Lv0, Cvl, Lv
  
! ###################################################################
  Rd = WGHT_INV(2)               *GRATIO
  Rdv=(WGHT_INV(1) - WGHT_INV(2))*GRATIO
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  IF      ( imixture .EQ. 0 ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = s(ij,1) - E_LOC
           
           Rl = Rd /Cd
          
           theta(ij) = T_LOC / P_LOC**Rl

        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC )  / ( Cd + s(ij,2) *Cdv )
           Cl = C_1_R/( Cd  + s(ij,2) *Cdv )
           Rl =(Rd  + s(ij,2) *Rdv) *Cl
           theta(ij) = T_LOC / P_LOC**Rl

        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC + s(ij,3)*Lv0 ) / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )
           Cl = C_1_R/( Cd  + s(ij,2) *Cdv )
           Rl =(Rd  + s(ij,2) *Rdv) *Cl
           Lv = Lv0 - T_LOC *Cvl
           theta(ij) = T_LOC / P_LOC**Rl *EXP(-Lv *s(ij,3) *Cl / T_LOC )

        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_THETA

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_WEIGHT_INPLACE(nx,ny,nz, weight, a)

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(*),          INTENT(IN)    :: weight
  TREAL, DIMENSION(nx,ny*nz),   INTENT(INOUT) :: a

! -------------------------------------------------------------------
  TINTEGER jk, j

! ###################################################################
  DO jk = 1,ny*nz
     j = MOD(jk-1,ny) +1
     
     a(1:nx,jk) = a(1:nx,jk) *weight(j)
     
  ENDDO
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_WEIGHT_INPLACE

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_WEIGHT_OUTPLACE(nx,ny,nz, weight, a, b)

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(*),          INTENT(IN)  :: weight
  TREAL, DIMENSION(nx,ny*nz),   INTENT(IN)  :: a
  TREAL, DIMENSION(nx,ny*nz),   INTENT(OUT) :: b

! -------------------------------------------------------------------
  TINTEGER jk, j

! ###################################################################
  DO jk = 1,ny*nz
     j = MOD(jk-1,ny) +1
     
     b(1:nx,jk) = a(1:nx,jk) *weight(j)
     
  ENDDO

  RETURN
END SUBROUTINE THERMO_ANELASTIC_WEIGHT_OUTPLACE

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_WEIGHT_ADD(nx,ny,nz, weight, a, b)

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(*),          INTENT(IN)    :: weight
  TREAL, DIMENSION(nx,ny*nz),   INTENT(IN)    :: a
  TREAL, DIMENSION(nx,ny*nz),   INTENT(INOUT) :: b

! -------------------------------------------------------------------
  TINTEGER jk, j

! ###################################################################
  DO jk = 1,ny*nz
     j = MOD(jk-1,ny) +1
     
     b(1:nx,jk) = b(1:nx,jk) +a(1:nx,jk) *weight(j)
     
  ENDDO

  RETURN
END SUBROUTINE THERMO_ANELASTIC_WEIGHT_ADD

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_LWP(nx,ny,nz, g, r, ql, lwp, wrk1d,wrk3d)

  USE DNS_TYPES, ONLY : grid_dt

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TYPE(grid_dt),              INTENT(IN)    :: g
  TREAL, DIMENSION(*),        INTENT(IN)    :: r
  TREAL, DIMENSION(nx*nz,ny), INTENT(IN)    :: ql
  TREAL, DIMENSION(nx,nz),    INTENT(OUT)   :: lwp
  TREAL, DIMENSION(ny),       INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(nx,ny,nz), INTENT(INOUT) :: wrk3d

! -------------------------------------------------------------------
  TINTEGER i,j,k
  TREAL SIMPSON_NU

! ###################################################################
  CALL THERMO_ANELASTIC_WEIGHT_OUTPLACE(nx,ny,nz, r, ql, wrk3d)

  DO k = 1,nz
     DO i = 1,nx
        DO j = 1,ny
           wrk1d(j) = wrk3d(i,j,k)
        ENDDO
        lwp(i,k) = SIMPSON_NU(ny, wrk1d, g%nodes)
     ENDDO
  ENDDO
     
  RETURN
END SUBROUTINE THERMO_ANELASTIC_LWP
