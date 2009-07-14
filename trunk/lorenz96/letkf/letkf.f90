PROGRAM letkf
!=======================================================================
! 4D-LETKF with Lorenz-96
!=======================================================================
  USE common
  USE common_letkf
  USE lorenz96
!  USE lorenz96_oro
  USE h_ope

  IMPLICIT NONE

  INTEGER,PARAMETER :: ndays=360
  INTEGER,PARAMETER :: nt=ndays*4
  INTEGER,PARAMETER :: nwindow=1 ! time window for 4D-LETKF
  INTEGER,PARAMETER :: nspinup=180*4 ! time steps for spin-up
  REAL(r_size),PARAMETER :: xlocal=1.0d0 ! localization scale
  REAL(r_size),PARAMETER :: tlocal=3.0d0 ! time localization scale
  REAL(r_size),PARAMETER :: parm_infl=0.04d0 ! inflation parameter
  REAL(r_size) :: xmaxloc
  REAL(r_size) :: obserr=1.0d0
  REAL(r_sngl) :: y4(ny)
  REAL(r_sngl) :: x4(nx)
  REAL(r_size) :: xnature(nx,nt)
  REAL(r_size) :: xa(nx,nbv,nwindow)
  REAL(r_size) :: xf(nx,nbv,nwindow)
  REAL(r_size) :: dxf(nx,nbv,nwindow)
  REAL(r_size) :: xm(nx,nwindow)
  REAL(r_size) :: y(ny,nwindow)
  REAL(r_size) :: d(ny,nwindow)
  REAL(r_size) :: h4d(ny,nx,nwindow)
  REAL(r_size) :: hxf(ny,nbv,nwindow)
  REAL(r_size) :: hxfm(ny,nwindow)
  REAL(r_size) :: hdxf(ny,nbv,nwindow)
  REAL(r_size) :: rdiag_loc(ny*nwindow)
  REAL(r_size) :: d_loc(ny*nwindow)
  REAL(r_size) :: hdxf_loc(ny*nwindow,nbv)
  REAL(r_size) :: trans(nbv,nbv)
  REAL(r_size) :: dist,tdif
  REAL(r_size) :: rmse
  INTEGER :: irmse
  INTEGER :: ktoneday
  INTEGER :: ktcyc
  INTEGER :: i,j,n,nn,it,ios
  INTEGER :: ix
  INTEGER :: ny_loc
  CHARACTER(10) :: initfile='init00.dat'
!-----------------------------------------------------------------------
! model parameters
!-----------------------------------------------------------------------
  dt=0.005d0
  force=8.0d0
  oneday=0.2d0
  ktoneday = INT(oneday/dt)
  ktcyc = ktoneday/4
  xmaxloc = xlocal * 2.0d0 * SQRT(10.0d0/3.0d0)
!-----------------------------------------------------------------------
! nature
!-----------------------------------------------------------------------
  OPEN(10,FILE='nature.dat',FORM='unformatted')
  DO i=1,nt
    READ(10) x4
    xnature(:,i) = REAL(x4,r_size)
  END DO
  CLOSE(10)
!-----------------------------------------------------------------------
! initial conditions 'initXX.dat'
!-----------------------------------------------------------------------
  DO i=1,nbv
    WRITE(initfile(5:6),'(I2.2)') i-1
    OPEN(10,FILE=initfile,FORM='unformatted')
    READ(10) xf(:,i,1)
    CLOSE(10)
  END DO
!-----------------------------------------------------------------------
! main
!-----------------------------------------------------------------------
  !
  ! input files
  !
  OPEN(11,FILE='obs.dat',FORM='unformatted')
  !
  ! output files
  !
  OPEN(90,FILE='guesmean.dat',FORM='unformatted')
  OPEN(91,FILE='analmean.dat',FORM='unformatted')
  OPEN(92,FILE='gues.dat',FORM='unformatted')
  OPEN(93,FILE='anal.dat',FORM='unformatted')
  !>>>
  !>>> LOOP START
  !>>>
  rmse = 0.0d0
  irmse = 0
  it=1
  DO
    !
    ! read obs
    !
    DO i=1,nwindow
      READ(11) y4
      y(:,i) = REAL(y4,r_size)
    END DO
    !
    ! 4d first guess
    !
    IF(nwindow > 1) THEN
      DO i=2,nwindow
        DO j=1,nbv
          CALL tinteg_rk4(ktcyc,xf(:,j,i-1),xf(:,j,i))
        END DO
      END DO
    END IF
    !
    ! ensemble mean -> xm
    !
    DO j=1,nwindow
      DO i=1,nx
        CALL com_mean(nbv,xf(i,:,j),xm(i,j))
      END DO
    END DO
    !
    ! ensemble ptb -> dxf
    !
    DO j=1,nwindow
      DO i=1,nbv
        dxf(:,i,j) = xf(:,i,j) - xm(:,j)
      END DO
    END DO
    !
    ! output first guess
    !
    DO j=1,nwindow
      x4 = xm(:,j)
      WRITE(90) x4
      DO i=1,nbv
        x4 = xf(:,i,j)
        WRITE(92) x4
      END DO
    END DO
    !---------------
    ! analysis step
    !---------------
    !
    ! hxf = H xf
    !
    DO n=1,nwindow
      DO j=1,nbv
        CALL set_h(xf(:,j,n))
        h4d(:,:,n) = h
        hxf(:,j,n) = h4d(:,1,n) * xf(1,j,n)
        DO i=2,nx
          hxf(:,j,n) = hxf(:,j,n) +  h4d(:,i,n) * xf(i,j,n)
        END DO
      END DO
    END DO
    !
    ! hxfm = mean(H xf)
    !
    DO n=1,nwindow
      DO i=1,ny
        CALL com_mean(nbv,hxf(i,:,n),hxfm(i,n))
      END DO
    END DO
    !
    ! d = y - hxfm
    !
    d = y - hxfm
    !
    ! hdxf
    !
    DO n=1,nwindow
      DO i=1,nbv
        hdxf(:,i,n) = hxf(:,i,n) - hxfm(:,n)
      END DO
    END DO
    !
    ! LETKF
    !
    DO nn=1,nwindow
      DO ix=1,nx
        ny_loc = 0
        DO n=1,nwindow
          tdif = REAL(ABS(nn-n),r_size)
          DO i=1,ny
            DO j=1,nx
              IF(MAXVAL(h4d(i,:,n)) == h4d(i,j,n)) EXIT
            END DO
            dist = REAL(MIN(ABS(j-ix),nx-ABS(j-ix)),r_size)
            IF(dist < xmaxloc) THEN
              ny_loc = ny_loc+1
              d_loc(ny_loc) = d(i,n)
              rdiag_loc(ny_loc) = obserr**2 &
                & * exp(0.5 * (dist/xlocal)**2) &! obs localization
                & * exp(0.5 * (tdif/tlocal)**2)  ! time localization
              hdxf_loc(ny_loc,:) = hdxf(i,:,n)
            END IF
          END DO
        END DO
        IF(ny_loc > 0) THEN
          CALL letkf_core(ny*nwindow,ny_loc,hdxf_loc,rdiag_loc,d_loc,parm_infl,trans)
          DO j=1,nbv
            xa(ix,j,nn) = xm(ix,nn)
            DO i=1,nbv
              xa(ix,j,nn) = xa(ix,j,nn) + dxf(ix,i,nn) * trans(i,j)
            END DO
          END DO
        ELSE
          DO j=1,nbv
            xa(ix,j,nn) = xm(ix,nn) + dxf(ix,j,nn)
          END DO
        END IF
      END DO
    END DO
    !
    ! ensemble mean
    !
    DO n=1,nwindow
      DO i=1,nx
        CALL com_mean(nbv,xa(i,:,n),xm(i,n))
      END DO
    END DO
    !
    ! output analysis
    !
    DO n=1,nwindow
      x4 = xm(:,n)
      WRITE(91) x4
      DO i=1,nbv
        x4 = xa(:,i,n)
        WRITE(93) x4
      END DO
    END DO
    !
    ! RMSE
    !
    IF(it > nspinup) THEN
      rmse = rmse + SQRT(SUM((xm(:,nwindow)-xnature(:,it+nwindow-1))**2)/REAL(nx,r_size))
      irmse = irmse + 1
    END IF
    !---------------
    ! forecast step
    !---------------
    DO i=1,nbv
      CALL tinteg_rk4(ktcyc,xa(:,i,nwindow),xf(:,i,1))
    END DO
    it = it+nwindow
    IF(it > nt) EXIT
  END DO
  !<<<
  !<<< LOOP END
  !<<<
  CLOSE(11)
  CLOSE(90)
  CLOSE(91)
  CLOSE(92)
  CLOSE(93)

  rmse = rmse / REAL(irmse,r_size)
  PRINT '(A,F12.5)','RMSE = ',rmse

  STOP
END PROGRAM letkf
