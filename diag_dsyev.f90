!------------------------------------------------------------------------!
! Example of matrix diagonalization using LAPACK routine dsyev.f and the ! 
! Fortran 90 interface diasym.f90.                                       !
!------------------------------------------------------------------------!
! Input from file 'mat.dat' (matrix to be diagonalized):                 !
! line 1      : order of the symmetric matrix [M]                        !
! lines 2-n+1 : rows of the matrix                                       !
!------------------------------------------------------------------------!
! Output in file 'dia.dat':                                              !
! - eigenvalues                                                          !
! - eigenvectors (diagonalizing matrix [D])                              !
! - the original matrix [M] transformed by [D]; [1/D][M][D]              !
!------------------------------------------------------------------------!

!---------------!
 program diatest

!---------------!
 implicit none

 integer :: i,j,n,info
 real(8), allocatable :: m0(:,:),m1(:,:),m2(:,:),eig(:)
 real(8) :: elem,t1,t2

 logical :: gen_mat = .true.
 integer :: mat_len
 integer(8) :: seed(1000000)
 integer :: verbose = 1


 !open(10,file='mat.dat',status='old')
 !
 read(5,*) n
 allocate (m0(n,n))
 allocate (m1(n,n))
 allocate (m2(n,n))
 allocate (eig(n))
 !
 m0 = 0
 !
 !if (.not.gen_mat) then 
 !
 i = 0 ; j=0
 !
 do while(.true.)
   read( 5, *, iostat=info ) i, j, elem
   if ( info /= 0 ) exit   ! exit loop when reach EOF
   m0(i,j) = elem
   m0(j,i) = elem
 enddo
 !
 ! a random matrix is generated
 if (i==0.and.j==0) then 
   !
   do j = 1,n
      do i=j,n
        !if(i>=j)  m0(i,j) = real(rrr(i),8)/real(1099511627776.0,8)
        if(i /= j) then 
          m0(i,j) = real(rrr(j),8)/real(1099511627776.0,8)
          m0(j,i) = m0(i,j)
        endif 
        if(i==j)  m0(i,j) = m0(i,j) + real(100.0,8) * real(j,8)
      enddo
   enddo
   !
 end if 
 !close(10)

 m1(:,:)=m0(:,:)

 t1 = get_real_time()

 call diasym(m1,eig,n)

 t2 = get_real_time()
   
 !open(10,file='dia.dat',status='replace')
 !
 write(6,*)'Eigenvalues:'
 do i=1,n
    write(6,10)i,eig(i)
    10 format(I5,'   ',f15.8)
 enddo
 write(6,*)
 
 if (verbose>=4) then 
  write(6,*)'Eigenvectors:'
  do i=1,n
     write(6,20)i,m1(:,i)
     20 format(i3,'   ',100f14.8)
  enddo
  write(6,*)
 endif
  
 m2=matmul(transpose(m1),m0)
 m0=matmul(m2,m1)

 if (verbose>=4) then 
   write(6,*)'Transformed matrix (check):'
   do i=1,n
      write(6,30)m0(:,i)
      30 format(100(1x,g15.8))
   enddo
   write(6,*)
 endif
 !
 write(6,"('Time = ',f12.1,' s')") t2-t1
 !

 !close(10)

 deallocate(m0); deallocate(m1); deallocate(m2);  deallocate(eig)
 !
 contains 
  !
  integer(8) function  rrr(i)
   integer :: i
   seed(i) = mod((1103515245 * seed(i) + 12345),1099511627776);
   rrr=seed(i)
   return
  end function rrr
  !
  function get_real_time() result(t)
    real(8) :: t
    !
    integer         :: count, count_rate, count_max
    real(8), save :: overflow   =  0
    integer, save   :: last_count = -1
    !
    call system_clock(count,count_rate,count_max)
    !
    ! Try to detect a rollover
    !
    if (count<last_count) then
      overflow = overflow + count_max
    end if
    last_count = count
    !
    ! Convert to seconds
    !
    t = (overflow+count)/count_rate
  end function get_real_time
  !
 end program diatest
!-------------------!

!---------------------------------------------------------!
!Calls the LAPACK diagonalization subroutine DSYEV        !
!input:  a(n,n) = real symmetric matrix to be diagonalized!
!            n  = size of a                               !
!output: a(n,n) = orthonormal eigenvectors of a           !
!        eig(n) = eigenvalues of a in ascending order     !
!---------------------------------------------------------!
 subroutine diasym(a,eig,n)
 implicit none

 integer n,l,inf
 real*8  a(n,n),eig(n),work(n*(3+n/2))

 l=n*(3+n/2)
 call dsyev('V','U',n,a,n,eig,work,l,inf)

 end subroutine diasym
!---------------------!

