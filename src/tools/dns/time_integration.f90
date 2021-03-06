#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!#
!# Performing the time integration over a given number of steps
!#
!########################################################################
SUBROUTINE TIME_INTEGRATION(q,hq, s,hs, txc, wrk1d,wrk2d,wrk3d, l_q, l_hq, l_txc, l_comm)

  USE DNS_CONSTANTS, ONLY : tag_flow, tag_scal, tag_part, tag_traj, lfile
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, isize_txc_field
  USE DNS_GLOBAL, ONLY : isize_particle
  USE DNS_GLOBAL, ONLY : imode_sim
  USE DNS_GLOBAL, ONLY : icalc_flow, icalc_scal, icalc_part
  USE DNS_GLOBAL, ONLY : visc
  USE DNS_GLOBAL, ONLY : itime, rtime
  USE TIME
  USE DNS_LOCAL
  USE DNS_TOWER
  USE LAGRANGE_GLOBAL, ONLY : itrajectory, l_g
  USE STATISTICS
  USE PARTICLE_TRAJECTORIES
  USE AVG_SCAL_ZT
  USE PLANES
#ifdef LES
  USE LES_GLOBAL, ONLY : iles
#endif

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(isize_field,*) :: q,hq, s,hs
  TREAL, DIMENSION(isize_txc_field,6),    INTENT(INOUT) :: txc
  TREAL, DIMENSION(*)             :: wrk1d, wrk2d, wrk3d

  TREAL, DIMENSION(isize_particle,*) :: l_q, l_hq, l_txc
  TREAL, DIMENSION(*)                :: l_comm

  ! -------------------------------------------------------------------
  CHARACTER*32 fname
  CHARACTER*250 line1
  LOGICAL flag_save

  ! ###################################################################
  WRITE(line1,*) itime; line1 = 'Starting time integration at It'//TRIM(ADJUSTL(line1))//'.'
  CALL IO_WRITE_ASCII(lfile,line1)

  DO
    IF ( itime >= nitera_last   ) EXIT
    IF ( INT(logs_data(1)) /= 0 ) EXIT

    CALL TIME_RUNGEKUTTA(q,hq, s,hs, txc, wrk1d,wrk2d,wrk3d, l_q, l_hq, l_txc, l_comm)

    itime = itime + 1
    rtime = rtime + dtime

    ! -----------------------------------------------------------------------
    IF ( MOD(itime-nitera_first,FilterDomainStep) == 0 ) THEN
      IF ( MOD(itime-nitera_first,nitera_stats)  == 0 ) THEN; flag_save = .TRUE.
      ELSE;                                                   flag_save = .FALSE.
      ENDIF
      CALL DNS_FILTER(flag_save, q,s, txc, wrk1d,wrk2d,wrk3d)
    ENDIF

    ! -----------------------------------------------------------------------
    IF ( flag_viscosity ) THEN          ! Change viscosity if necessary
      visc = visc +visc_rate *dtime
      IF ( rtime .GT. visc_time ) THEN
        visc = visc_stop                ! Fix new value without any roundoff
        flag_viscosity = .FALSE.
      ENDIF
    END IF

    ! -----------------------------------------------------------------------
    CALL TIME_COURANT(q, wrk3d)

    ! ###################################################################
    ! The rest: Logging, postprocessing and saving
    ! ###################################################################
    IF ( MOD(itime-nitera_first,nitera_log) == 0 .OR. INT(logs_data(1)) /= 0 ) THEN ! Log files
      CALL DNS_LOGS(i2)
#ifdef LES
      IF ( iles == 1 ) CALL LES_LOGS(i2)
#endif
    ENDIF

    ! -----------------------------------------------------------------------
    ! Accumulate statistics in spatially evolving cases
    IF ( MOD(itime-nitera_first,nitera_stats_spa) == 0 ) THEN
      CALL AVG_FLOW_ZT_REDUCE(q, hq,txc, mean_flow, wrk2d,wrk3d)
      IF ( icalc_scal == 1 ) THEN
        CALL AVG_SCAL_ZT_REDUCE(q,s, hq,txc, mean_scal, wrk2d,wrk3d)
      ENDIF
    ENDIF

    ! -----------------------------------------------------------------------
    IF ( tower_mode == 1 ) THEN
      CALL DNS_TOWER_ACCUMULATE(q,1,wrk1d)
      CALL DNS_TOWER_ACCUMULATE(s,2,wrk1d)
    ENDIF

    ! -----------------------------------------------------------------------
    IF ( itrajectory /= LAG_TRAJECTORY_NONE ) THEN
      CALL PARTICLE_TRAJECTORIES_ACCUMULATE(q,s, txc, l_g,l_q,l_hq,l_txc,l_comm, wrk2d,wrk3d)
    END IF

    ! -----------------------------------------------------------------------
    IF ( MOD(itime-nitera_first,nitera_stats) == 0 ) THEN ! Calculate statistics
      IF     ( imode_sim == DNS_MODE_TEMPORAL ) THEN
        CALL STATISTICS_TEMPORAL(q,s,hq, txc, wrk1d,wrk2d,wrk3d, l_q,l_txc,l_comm)
      ELSE IF ( imode_sim == DNS_MODE_SPATIAL ) THEN
        CALL STATISTICS_SPATIAL(txc, wrk1d,wrk2d)
      ENDIF
    ENDIF

    ! -----------------------------------------------------------------------
    IF ( MOD(itime-nitera_first,nitera_save) == 0 .OR. &        ! Save restart files
        itime == nitera_last .OR. INT(logs_data(1)) /= 0 ) THEN ! Secure that one restart file is saved

      IF ( icalc_flow == 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
        CALL DNS_WRITE_FIELDS(fname, i2, imax,jmax,kmax, inb_flow, isize_field, q, wrk3d)
      ENDIF

      IF ( icalc_scal == 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
        CALL DNS_WRITE_FIELDS(fname, i1, imax,jmax,kmax, inb_scal, isize_field, s, wrk3d)
      ENDIF

      IF ( tower_mode == 1 ) THEN
        CALL DNS_TOWER_WRITE(wrk3d)
      ENDIF

      IF ( icalc_part == 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(fname))
        CALL IO_WRITE_PARTICLE(fname, l_g, l_q)
        IF ( itrajectory /= LAG_TRAJECTORY_NONE ) THEN
          WRITE(fname,*) itime; fname =TRIM(ADJUSTL(tag_traj))//TRIM(ADJUSTL(fname))
          CALL PARTICLE_TRAJECTORIES_WRITE(fname)
        END IF
      END IF

      IF ( imode_sim == DNS_MODE_SPATIAL .AND. nitera_stats_spa > 0 ) THEN ! Spatial; running averages
        WRITE(fname,*) itime; fname='st'//TRIM(ADJUSTL(fname))
        CALL IO_WRITE_AVG_SPATIAL(fname, mean_flow, mean_scal)
      ENDIF

    ENDIF

    ! -----------------------------------------------------------------------
    IF ( MOD(itime-nitera_first,nitera_pln) == 0 ) THEN
      CALL PLANES_SAVE( q,s, txc(1,1), txc(1,2),txc(1,3),txc(1,4), wrk1d,wrk2d,wrk3d )
    ENDIF

  ENDDO

  RETURN
END SUBROUTINE TIME_INTEGRATION
