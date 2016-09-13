MODULE common_mtx
!=======================================================================
!
! [PURPOSE:] Matrix Functions
!
! [CREATED:] 07/20/2004 Takemasa Miyoshi
! [UPDATED:] 10/16/2004 Takemasa Miyoshi
!
! [PUBLIC:]
!   mtx_eigen  : eigenvalue decomposition
!   mtx_inv    : real symmetric matrix inverse
!   mtx_sqrt   : real symmetric matrix square root
!
! [REFERENCES:]
!    Core subroutines are adapted from netlib.org
!
! [HISTORY:]
!  07/20/2003 Takemasa Miyoshi  Created at University of Maryland, College Park
!
!=======================================================================
  USE common

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mtx_eigen, mtx_inv, mtx_sqrt, mtx_inv_rg

CONTAINS
!=======================================================================
!  Eigenvalue decomposition using subroutine rs
!    INPUT
!      INTEGER :: imode           : mode switch (0: only eiven values)
!      INTEGER :: n               : dimension of matrix
!      REAL(r_size) :: a(n,n)     : input matrix
!    OUTPUT
!      REAL(r_size) :: eival(n)   : eiven values in decending order
!                                   i.e. eival(1) is the largest
!      REAL(r_size) :: eivec(n,n) : eiven vectors
!      INTEGER :: nrank_eff       : number of positive eivenvalues
!=======================================================================
SUBROUTINE mtx_eigen(imode,n,a,eival,eivec,nrank_eff)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: imode ! 0: calculate only eigen values
  INTEGER,INTENT(IN) :: n
  REAL(r_size),INTENT(IN) :: a(1:n,1:n)
  REAL(r_size),INTENT(OUT) :: eival(1:n)
  REAL(r_size),INTENT(OUT) :: eivec(1:n,1:n)
  INTEGER,INTENT(OUT) :: nrank_eff

  REAL(r_dble) :: a8(n,n)
  REAL(r_dble) :: eival8(n)
  REAL(r_dble) :: eivec8(n,n)
  REAL(r_dble) :: wrk1(n)
  REAL(r_dble) :: wrk2(n)
  INTEGER :: ierr,i,j

  a8 = a
  eivec8 = 0.0d0
  CALL rs(n,n,a8,eival8,imode,eivec8,wrk1,wrk2,ierr)
  IF( ierr/=0 ) THEN
    WRITE(6,'(A,I4)') '!!! ERROR (mtx_eigen): rs error code is ',ierr
    STOP 2
  END IF

  nrank_eff = n
  IF( eival8(n) > 0 ) THEN
    DO i=1,n
      IF( eival8(i) < ABS(eival8(n))*SQRT(EPSILON(eival8)) ) THEN
        nrank_eff = nrank_eff - 1
        eival8(i) = 0.0d0
        eivec8(:,i) = 0.0d0
      END IF
    END DO
  ELSE
    WRITE(6,'(A)') '!!! ERROR (mtx_eigen): All Eigenvalues are below 0'
    STOP 2
  END IF

  IF( nrank_eff<n .AND. eival8(1)/=0 ) THEN
    j = 0
    DO i=n,1,-1
      IF( eival8(i) == 0 ) THEN
        eival8(i) = eival8(n-nrank_eff-j)
        eivec(:,i) = eivec8(:,n-nrank_eff-j)
        eival8(n-nrank_eff-j) = 0.0d0
        eivec8(:,n-nrank_eff-j) = 0.0d0
        j = j+1
      END IF
    END DO
  END IF

  DO i=1,n
    eival(i) = eival8(n+1-i)
    eivec(:,i) = eivec8(:,n+1-i)
  END DO

  RETURN
END SUBROUTINE mtx_eigen
!=======================================================================
!  Real symmetric matrix inversion using subroutine dspdi
!    INPUT
!      INTEGER :: n               : dimension of matrix
!      REAL(r_size) :: a(n,n)     : input matrix (real symmetric)
!    OUTPUT
!      REAL(r_size) :: ainv(n,n)  : inverse of a
!=======================================================================
SUBROUTINE mtx_inv(n,a,ainv)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: n
  REAL(r_size),INTENT(IN) :: a(1:n,1:n)
  REAL(r_size),INTENT(OUT) :: ainv(1:n,1:n)

  REAL(r_dble) :: acmp(n*(n+1)/2)
  REAL(r_dble) :: det(2)
  REAL(r_dble) :: work(n)
  INTEGER :: kpvt(n)
  INTEGER :: inert(3)
  INTEGER :: info
  INTEGER :: i,j,k

  IF(n==1) THEN
    ainv(1,1) = 1.0d0 / a(1,1)
  ELSE

!-----------------------------------------------------------------------
!  Packed form of matrix
!-----------------------------------------------------------------------
  k=0
  DO j=1,n
    DO i=1,j
      k = k+1
      acmp(k) = a(i,j)
    END DO
  END DO
!-----------------------------------------------------------------------
!  dspfa
!-----------------------------------------------------------------------
  CALL dspfa(acmp,n,kpvt,info)
  IF(info /= 0) THEN
    WRITE(6,'(A,I4)') '!!! ERROR (mtx_inv): dspfa error code is ',info
    STOP 3
  END IF
!-----------------------------------------------------------------------
!  dspdi
!-----------------------------------------------------------------------
  CALL dspdi(acmp,n,kpvt,det,inert,work,001)
!-----------------------------------------------------------------------
!  unpack matrix
!-----------------------------------------------------------------------
  k=0
  DO j=1,n
    DO i=1,j
      k = k+1
      ainv(i,j) = acmp(k)
    END DO
  END DO

  DO j=1,n
    DO i=j+1,n
      ainv(i,j) = ainv(j,i)
    END DO
  END DO

  END IF

  RETURN
END SUBROUTINE mtx_inv
!=======================================================================
!  Compute square root of real symmetric matrix
!    INPUT
!      INTEGER :: n                : dimension of matrix
!      REAL(r_size) :: a(n,n)      : input matrix (real symmetric)
!    OUTPUT
!      REAL(r_size) :: a_sqrt(n,n) : square root of a
!=======================================================================
SUBROUTINE mtx_sqrt(n,a,a_sqrt)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: n
  REAL(r_size),INTENT(IN) :: a(1:n,1:n)
  REAL(r_size),INTENT(OUT) :: a_sqrt(1:n,1:n)

  REAL(r_size) :: eival(n)   ! holds eivenvalue of a
  REAL(r_size) :: eivec(n,n) ! holds eivenvector of a
  REAL(r_size) :: wk(n,n)
  INTEGER :: i,j,k

  CALL mtx_eigen(1,n,a,eival,eivec,i)

  DO i=1,n
    wk(:,i) = eivec(:,i) * SQRT( eival(i) )
  END DO

!  a_sqrt = matmul(wk,transpose(eivec))
  DO j=1,n
    DO i=1,n
      a_sqrt(i,j) = wk(i,1)*eivec(j,1)
      DO k=2,n
        a_sqrt(i,j) = a_sqrt(i,j) + wk(i,k)*eivec(j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE mtx_sqrt
!=======================================================================
!  Compute inverse of a real matrix (not necessarily symmetric)
!    INPUT
!      INTEGER :: n            : dimension of matrix
!      REAL(r_size) :: aa(n,n) : input matrix (real symmetric)
!    OUTPUT
!      REAL(r_size) :: ff(n,n) : square root of a
!*** COPIED FROM 'A0568.NEW.FORT(MTXINV)' ON 1989.10.1
!    changed to free format by H.Yoshimura 2000.06.27
!    adapted by T.Miyoshi on 2005.10.31
!=======================================================================
SUBROUTINE mtx_inv_rg(n,aa,ff)
!
!##  MATRIX INVERSION
!##  AA IS THE MATRIX TO BE INVERTED
!##  FF IS THE INVERSE OF AA
!
  INTEGER,INTENT(IN) :: n
  REAL(r_size),INTENT(IN) :: aa(n,n)
  REAL(r_size),INTENT(OUT) :: ff(n,n)
!
  REAL(r_size) :: a(n,n)
  REAL(r_size) :: b(n,n)
  REAL(r_size) :: x(n,n)
!
  REAL(r_size) :: c,cc,xx
  INTEGER :: i,j,n1,k,kp,jx,ii,jr,jp
!-------------------------------------------------------
  n1=n-1
!
  do i=1,n
    do j=1,n
      a(i,j)=aa(i,j)
    end do
  end do
!
  do i=1,n
    do j=1,n
      b(i,j)=0.d0
      if( i == j ) b(i,j)=1.d0
    end do
  end do
!
  do j=1,n
    c=abs(a(1,j))
    do i=2,n
      c=max(c,abs(a(i,j)))
    end do
    c=1.d0/c
    do i=1,n
      a(i,j)=a(i,j)*c
    end do
    b(j,j)=b(j,j)*c
  end do
!
  do k=1,n1
    c=abs(a(k,k))
    kp=k+1
    jx=k
    do j=kp,n
      cc=abs(a(k,j))
      if ( cc < c ) cycle
      c=cc
      jx=j
    end do
    do i=k,n
      c=a(i,k)
      a(i,k)=a(i,jx)
      a(i,jx)=c
    end do
    do i=1,n
      c=b(i,k)
      b(i,k)=b(i,jx)
      b(i,jx)=c
    end do
    do j=kp,n
      c=a(k,j)/a(k,k)
      do ii=1,n
        b(ii,j)=b(ii,j)-c*b(ii,k)
      end do
      do i=k,n
        a(i,j)=a(i,j)-c*a(i,k)
      end do
    end do
  end do
!
  do ii=1,n
    x(ii,n)=b(ii,n)/a(n,n)
    do j=1,n1
      jr=n-j
      jp=jr+1
      xx=0.d0
      do i=jp,n
        xx=xx+a(i,jr)*x(ii,i)
      end do
      x(ii,jr)=(b(ii,jr)-xx)/a(jr,jr)
    end do
  end do
!
  do i=1,n
    do j=1,n
      ff(i,j)=x(i,j)
    end do
  end do
!
END SUBROUTINE mtx_inv_rg

END MODULE common_mtx
