PROGRAM letkf
!=======================================================================
! LETKF with Lorenz-96
!=======================================================================
  USE common
  USE common_letkf
  USE lorenz96
!  USE lorenz96_oro
  USE h_ope

  IMPLICIT NONE

  INTEGER,PARAMETER :: ndays=360
  INTEGER,PARAMETER :: nt=ndays*4
  REAL(r_size),PARAMETER :: xlocal=5.0d0 ! localization scale
  REAL(r_size),PARAMETER :: parm_infl=0.02d0 ! inflation parameter
  REAL(r_size) :: xmaxloc
  REAL(r_size) :: obserr=1.0d0
  REAL(r_sngl) :: y4(ny)
  REAL(r_sngl) :: x4(nx)
  REAL(r_size) :: xa(nx,nbv)
  REAL(r_size) :: xf(nx,nbv)
  REAL(r_size) :: dxf(nx,nbv)
  REAL(r_size) :: xm(nx)
  REAL(r_size) :: y(ny)
  REAL(r_size) :: d(ny)
  REAL(r_size) :: hxf(ny,nbv)
  REAL(r_size) :: hxfm(ny)
  REAL(r_size) :: hdxf(ny,nbv)
  REAL(r_size) :: rdiag_loc(ny)
  REAL(r_size) :: d_loc(ny)
  REAL(r_size) :: hdxf_loc(ny,nbv)
  REAL(r_size) :: trans(nbv,nbv)
  INTEGER :: ktoneday
  INTEGER :: ktcyc
  INTEGER :: i,j,it,ios
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
! initial conditions 'initXX.dat'
!-----------------------------------------------------------------------
  DO i=1,nbv
    WRITE(initfile(5:6),'(I2.2)') i-1
    OPEN(10,FILE=initfile,FORM='unformatted')
    READ(10) xf(:,i)
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
  DO it=1,nt
    !
    ! read obs
    !
    READ(11) y4
    y = REAL(y4,r_size)
    !
    ! ensemble mean -> xm
    !
    DO i=1,nx
      CALL com_mean(nbv,xf(i,:),xm(i))
    END DO
    !
    ! ensemble ptb -> dxf
    !
    DO i=1,nbv
      dxf(:,i) = xf(:,i) - xm
    END DO
    !
    ! output first guess
    !
    x4 = xm
    WRITE(90) x4
    DO i=1,nbv
      x4 = xf(:,i)
      WRITE(92) x4
    END DO
    !---------------
    ! analysis step
    !---------------
    !
    ! hxf = H xf
    !
    DO j=1,nbv
      CALL set_h(xf(:,j))
      hxf(:,j) = h(:,1) * xf(1,j)
      DO i=2,nx
        hxf(:,j) = hxf(:,j) +  h(:,i) * xf(i,j)
      END DO
    END DO
    !
    ! hxfm = mean(H xf)
    !
    DO i=1,ny
      CALL com_mean(nbv,hxf(i,:),hxfm(i))
    END DO
    !
    ! d = y - hxfm
    !
    d = y - hxfm
    !
    ! hdxf
    !
    DO i=1,nbv
      hdxf(:,i) = hxf(:,i) - hxfm
    END DO
    !
    ! LETKF
    !
    DO ix=1,nx
      ny_loc = 0
      DO i=1,ny
        DO j=1,nx
          IF(MAXVAL(h(i,:)) == h(i,j)) EXIT
        END DO
        IF(REAL(ABS(j-ix),r_size) < xmaxloc) THEN
          ny_loc = ny_loc+1
          d_loc(ny_loc) = d(i)
          rdiag_loc(ny_loc) = obserr**2 &
            & * exp(0.5 * (REAL(j-ix,r_size)/xlocal)**2) ! obs localization
          hdxf_loc(ny_loc,:) = hdxf(i,:)
        END IF
      END DO
      IF(ny_loc > 0) THEN
        CALL letkf_core(ny,ny_loc,hdxf_loc,rdiag_loc,d_loc,parm_infl,trans)
        DO j=1,nbv
          xa(ix,j) = xm(ix)
          DO i=1,nbv
            xa(ix,j) = xa(ix,j) + dxf(ix,i) * trans(i,j)
          END DO
        END DO
      ELSE
        DO j=1,nbv
          xa(ix,j) = xm(ix) + dxf(ix,j)
        END DO
      END IF
    END DO
    !
    ! ensemble mean
    !
    DO i=1,nx
      CALL com_mean(nbv,xa(i,:),xm(i))
    END DO
    !
    ! output analysis
    !
    x4 = xm
    WRITE(91) x4
    DO i=1,nbv
      x4 = xa(:,i)
      WRITE(93) x4
    END DO
    !---------------
    ! forecast step
    !---------------
    DO i=1,nbv
      CALL tinteg_rk4(ktcyc,xa(:,i),xf(:,i))
    END DO
  END DO
  !<<<
  !<<< LOOP END
  !<<<
  CLOSE(11)
  CLOSE(90)
  CLOSE(91)
  CLOSE(92)
  CLOSE(93)

  STOP
END PROGRAM letkf
