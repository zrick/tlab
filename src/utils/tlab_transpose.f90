!########################################################################
!# HISTORY
!#
!# 2011/11/01 - C. Ansorge
!#              Created
!#
!########################################################################
!#
!# Transposition using Cache-Blocking and OpenMP
!# transpose of the matrix a and place the transposed matrix in b
!# routine trans below is faster than TRANSPOSE routine from f90
!#
!########################################################################
subroutine TLab_Transpose(a, nra, nca, ma, b, mb)
    use TLab_Constants, only: wp, wi
    use TLab_OpenMP
    use Tlab_VARS, only : trans_time 
    implicit none

    integer(wi), intent(in) :: nra      ! Number of rows in a
    integer(wi), intent(in) :: nca      ! Number of columns in b
    integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
    integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
    real(wp), intent(in)    :: a(ma, *) ! Input array
    real(wp), intent(out)   :: b(mb, *) ! Transposed array

! -------------------------------------------------------------------
    integer(wi) jb, kb
    integer clock_0, clock_1, clock_cycle

#ifdef HLRS_HAWK
    parameter(jb=16, kb=8)
#else
    parameter(jb=64, kb=64)
#endif

    integer(wi) :: srt, end, siz

    integer(wi) k, j, jj, kk
    integer(wi) last_k, last_j

    CALL SYSTEM_CLOCK(clock_0,clock_cycle) 

! -------------------------------------------------------------------
#if defined(USE_MKL)
    call MKL_DOMATCOPY('c', 't', nra, nca, 1.0_wp, a, ma, b, mb)
#elif defined(USE_APU)
    if (  nca < nra .AND. nca < 2e4 ) THEN    ! This 'if' is a workaround for an int-overflow  bug in 
                                              ! OMP implementation of cray in cpe17
       !$omp target teams distribute parallel do collapse(2) default(none) &
       !$omp private(k,j) &
       !$omp shared(a,b,nca,nra)
       do k = 1, nca
          do j = 1, nra 
             b(k, j) = a(j,k)
          end do
       end do
       !$omp end target teams distribute parallel do
    else
       !$omp target teams distribute parallel do default(none) &
       !$omp private(k,j) &
       !$omp shared(a,b,nca,nra)
       do k = 1, nca
          do j = 1, nra 
             b(k, j) = a(j,k)
          end do
       end do
       !$omp end target teams distribute parallel do
    endif
#else 
!$omp parallel default(none) &
!$omp private(k,j,jj,kk,srt,end,siz,last_k,last_j) &
!$omp shared(a,b,nca,nra)

    call TLab_OMP_PARTITION(nca, srt, end, siz)

    kk = 1; jj = 1

    do k = srt, end - kb + 1, kb; 
        do j = 1, nra - jb + 1, jb; 
            do jj = j, j + jb - 1
                do kk = k, k + kb - 1
                    b(kk, jj) = a(jj, kk)
                end do
            end do
        end do
    end do

    last_k = kk
    last_j = jj

    do k = last_k, end
        do j = 1, nra
            b(k, j) = a(j, k)
        end do
    end do

    do k = srt, end
        do j = last_j, nra
            b(k, j) = a(j, k)
        end do
    end do

!$omp end parallel
#endif
    CALL SYSTEM_CLOCK(clock_1)
    trans_time = trans_time + real(clock_1 - clock_0)/ clock_cycle 
    
    return
end subroutine TLab_Transpose

!########################################################################
!########################################################################
subroutine TLab_Transpose_INT1(a, nra, nca, ma, b, mb)
    use TLab_Constants, only: wp, wi
    use TLab_OpenMP
    implicit none

    integer(wi), intent(in) :: nra      ! Number of rows in a
    integer(wi), intent(in) :: nca      ! Number of columns in b
    integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
    integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
    integer(1), intent(in)    :: a(ma, *) ! Input array
    integer(1), intent(out)   :: b(mb, *) ! Transposed array

! -------------------------------------------------------------------
    integer(wi) jb, kb
    parameter(jb=32, kb=32)

    integer(wi) :: srt, end, siz

    integer(wi) k, j, jj, kk
    integer(wi) last_k, last_j

! -------------------------------------------------------------------
!$omp parallel default(none) &
!$omp private(k,j,jj,kk,srt,end,siz,last_k,last_j) &
!$omp shared(a,b,nca,nra)

    call TLab_OMP_PARTITION(nca, srt, end, siz)

    kk = 1; jj = 1

    do k = srt, end - kb + 1, kb; 
        do j = 1, nra - jb + 1, jb; 
            do jj = j, j + jb - 1
                do kk = k, k + kb - 1
                    b(kk, jj) = a(jj, kk)
                end do
            end do
        end do
    end do

    last_k = kk
    last_j = jj

    do k = last_k, end
        do j = 1, nra
            b(k, j) = a(j, k)
        end do
    end do

    do k = srt, end
        do j = last_j, nra
            b(k, j) = a(j, k)
        end do
    end do

!$omp end parallel

    return
end subroutine TLab_Transpose_INT1

!########################################################################
!########################################################################
subroutine TLab_Transpose_COMPLEX(a, nra, nca, ma, b, mb)
    use TLab_Constants, only: wp, wi
    use TLab_OpenMP
    implicit none

    integer(wi), intent(in) :: nra      ! Number of rows in a
    integer(wi), intent(in) :: nca      ! Number of columns in b
    integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
    integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
    complex(wp), intent(in)    :: a(ma, *) ! Input array
    complex(wp), intent(out)   :: b(mb, *) ! Transposed array

! -------------------------------------------------------------------
    integer(wi) jb, kb
#ifdef HLRS_HAWK
    parameter(jb=16, kb=8)
#else
    parameter(jb=64, kb=64)
#endif

    integer(wi) :: srt, end, siz

    integer(wi) k, j, jj, kk
    integer(wi) last_k, last_j

! -------------------------------------------------------------------
!$omp parallel default(none) &
!$omp private(k,j,jj,kk,srt,end,siz,last_k,last_j) &
!$omp shared(a,b,nca,nra)

    call TLab_OMP_PARTITION(nca, srt, end, siz)

    kk = 1; jj = 1

    do k = srt, end - kb + 1, kb; 
        do j = 1, nra - jb + 1, jb; 
            do jj = j, j + jb - 1
                do kk = k, k + kb - 1
                    b(kk, jj) = a(jj, kk)
                end do
            end do
        end do
    end do

    last_k = kk
    last_j = jj

    do k = last_k, end
        do j = 1, nra
            b(k, j) = a(j, k)
        end do
    end do

    do k = srt, end
        do j = last_j, nra
            b(k, j) = a(j, k)
        end do
    end do

!$omp end parallel
    return
end subroutine TLab_Transpose_COMPLEX

