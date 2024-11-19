#include "dns_const.h"

program VBURGERS

  use TLAB_CONSTANTS
  use TLAB_VARS
  use TLAB_Workflow
  use TLAB_ARRAYS
  use TLAB_POINTERS_3D, only: tmp1
  use TLab_Memory, only : TLab_Initialize_Memory
#ifdef USE_MPI
  use MPI
  use TLabMPI_PROCS
  use TLabMPI_VARS
#endif
  use IO_FIELDS
  use OPR_PARTIAL
  use OPR_BURGERS
  use OPR_FILTERS
  USE FDM, only : g, FDM_Initialize
  implicit none

#ifdef USE_MPI
  real(wp) error2, dummy2
#else
  integer(wi), parameter :: ims_pro = 0
#endif

  real(wp), dimension(:, :, :), pointer :: a, b, c

  integer(wi) i, j, k, ig, bcs(2, 2)
  real(wp) dummy, error

  integer irun,nrun,stat
  integer clock_0, clock_1,clock_cycle
  integer clock_add0, clock_add1
  CHARACTER(len=64) :: nrun_string 
  real(wp), DIMENSION(:), ALLOCATABLE ::  runtime 
  real(wp) add_time
  
  trans_time = 0.0
  add_time = 0.0 

  call SYSTEM_CLOCK(clock_0,clock_cycle)
  IF ( COMMAND_ARGUMENT_COUNT() .GE. 1 ) THEN 
     call GetArg(1,nrun_String)
     read(nrun_string,*) nrun
  ELSE
     nrun = 1 
  ENDIF
  
  PRINT *,'EXECUTING ',nrun, ' RUNS for Performance Measurement'
  
  ALLOCATE(runtime(nrun))

  
  ! ###################################################################
  call TLAB_START()

  call Tlab_Initialize_Parameters(ifile) 
#ifdef USE_MPI
  call TLabMPI_Initialize()
#endif
  call NavierStokes_Initialize_Parameters(ifile) 
  
  inb_txc = 4

  call TLab_Initialize_Memory(__FILE__)

  a(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 2)
  b(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)
  c(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)

  visc = 1.0_wp/big_wp    ! inviscid

  call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:,1), wrk1d(:,2), wrk1d(:,3))
  call FDM_INITIALIZE(x, g(1), wrk1d(:,1),wrk1d(:,4))
  call FDM_INITIALIZE(y, g(2), wrk1d(:,2),wrk1d(:,4))
  call FDM_INITIALIZE(z, g(3), wrk1d(:,3),wrk1d(:,4))

  call Tlab_Initialize_Background() 

  bcs = 0

  do ig = 1, 3
     call OPR_FILTER_INITIALIZE(g(ig), Dealiasing(ig))
  end do

  ! ###################################################################
  ! Define forcing term
  ! ###################################################################
  call IO_READ_FIELDS('field.inp', IO_SCAL, imax, jmax, kmax, 1, 0, a)

  visc = 1.0_wp/big_wp

  ! Call everything once to get the memory initialized on the APUs
  ! (avoids measuring the first-touch penalty) 
  call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), a, b, c)
  call OPR_BURGERS_X(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(1), a, a, c, tmp1)
  call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), a, b, c)
  call OPR_BURGERS_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(2), a, a, c, tmp1)
  call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), a, b, c)
  call OPR_BURGERS_Z(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(3), a, a, c, tmp1)

  ! reset time for transpositions to avoid measuring first-touch penalty
  ! (else nrun needs to be incremented by one in normalization) 
  trans_time = 0. 
  
  DO irun=1,nrun

  
     call SYSTEM_CLOCK(clock_0) 
     ! ###################################################################
     call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), a, b, c)
     do k = 1, kmax
        do j = 1, jmax
           do i = 1, imax
              b(i, j, k) = b(i, j, k)*visc - a(i, j, k)*c(i, j, k)
           end do
        end do
     end do

     call OPR_BURGERS_X(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(1), a, a, c, tmp1)
     call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), a, b, c)
     CALL SYSTEM_CLOCK(clock_add0) 
     !$omp target teams distribute parallel do collapse(2) default (none) &
     !$omp private(i,j,k) &
     !$omp shared(a,b,c,imax,jmax,kmax,visc) 
     do k = 1, kmax
        do j = 1, jmax
           do i = 1, imax
              b(i, j, k) = b(i, j, k)*visc - a(i, j, k)*c(i, j, k)
           end do
        end do
     end do
     !$omp end target teams distribute parallel do
     CALL SYSTEM_CLOCK(clock_add1)
     add_time = add_time + real(clock_add1 - clock_add0) / clock_cycle 
     call OPR_BURGERS_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(2), a, a, c, tmp1)

     ! ###################################################################
     if (g(3)%size > 1) then
        CALL SYSTEM_CLOCK(clock_add0)
        !$omp target teams distribute parallel do collapse(2) default(none) &
        !$omp private(i,j,k) &
        !$omp shared(a,b,c,imax,jmax,kmax,visc) 
        do k = 1, kmax
           do j = 1, jmax
              do i = 1, imax
                 b(i, j, k) = b(i, j, k)*visc - a(i, j, k)*c(i, j, k)
                 ! b(i, j, k) = b(i, j, k)*visc*ribackground(j) - a(i, j, k)*c(i, j, k)
              end do
           end do
        end do
        !$omp end target teams distribute parallel do
        CALL SYSTEM_CLOCK(clock_add1)
        add_time = add_time + real(clock_add1 - clock_add0) / clock_cycle 
        call OPR_BURGERS_Z(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(3), a, a, c, tmp1)

     end if
     call SYSTEM_CLOCK(clock_1)
     runtime(irun) = real(clock_1-clock_0)/clock_cycle
  ENDDO
  PRINT 100,SUM(runtime)/nrun, MINVAL(runtime),MAXVAL(runtime)
  PRINT 101, 'Transpos',trans_time/nrun, 100*trans_time/SUM(runtime)
  PRINT 101, 'Addition',add_time/nrun, 100*add_time/SUM(runtime) 
100 FORMAT('T MEAN|MIN|MAX [s] : ', F8.4, 1x, F8.4, 1x , F8.4)
101 FORMAT('Time per run in ',A9,'[s]:', F8.4,'s (', F3.0,'%)') 
  call TLAB_STOP(0)
end program VBURGERS
