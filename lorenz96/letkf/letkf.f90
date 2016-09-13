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

  LOGICAL,PARAMETER :: msw_detailout=.FALSE.
  INTEGER,PARAMETER :: ndays=360*5
  INTEGER,PARAMETER :: nt=ndays*4
  INTEGER,PARAMETER :: nwindow=1 ! time window for 4D-LETKF
  INTEGER,PARAMETER :: nspinup=360*3*4 ! time steps for spin-up
  INTEGER,PARAMETER :: msw_local=0 ! localization mode switch
! msw_local : localization mode switch
!  0 : fixed localization
!  1 : adaptive localization
!  2 : combination of both
  REAL(r_size) :: xlocal=3.0d0 ! localization scale
  REAL(r_size),PARAMETER :: tlocal=2.0d0 ! time localization scale
! negative value for no time localization
  REAL(r_size) :: sa=1.0d0 ! adaptive localization parameter
  REAL(r_size) :: sb=1.0d0 ! adaptive localization parameter
  REAL(r_size),PARAMETER :: msw_infl=1.05d0 ! inflation mode switch
! msw_infl : inflation mode switch
!  < 0 : adaptive inflation
!  > 0 : fixed inflation value
  REAL(r_size) :: parm_infl(nx,nt) ! inflation parameter
  REAL(r_size) :: parm
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
  REAL(r_size) :: rloc_loc(ny*nwindow)
  REAL(r_size) :: d_loc(ny*nwindow)
  REAL(r_size) :: hdxf_loc(ny*nwindow,nbv)
  REAL(r_size) :: trans(nbv,nbv)
  REAL(r_size) :: dist,tdif
  REAL(r_size) :: rmse_t(nt),sprd_t(nt),infl_t(nt)
  REAL(r_size) :: rmse_x(nx),sprd_x(nx),infl_x(nx)
  REAL(r_size) :: rmseave,sprdave,inflave
  REAL(r_size) :: obsloc(3),wa,wb
  INTEGER :: irmse
  INTEGER :: ktoneday
  INTEGER :: ktcyc
  INTEGER :: i,j,n,nn,it,ios
  INTEGER :: ix
  INTEGER :: ny_loc
  INTEGER :: nbv2
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
  nbv2 = CEILING(REAL(nbv)/2.0)
  PRINT '(A)'     ,'==========LETKF settings=========='
  PRINT '(A,I8)'  ,' nbv       : ',nbv
  PRINT '(A,I8)'  ,' ny        : ',ny
  PRINT '(A,I8)'  ,' nwindow   : ',nwindow
  PRINT '(A,I8)'  ,' msw_local : ',msw_local
  PRINT '(A,F8.1)',' xlocal    : ',xlocal
  PRINT '(A,F8.1)',' tlocal    : ',tlocal
  PRINT '(A,F8.1)',' sa        : ',sa
  PRINT '(A,F8.1)',' sb        : ',sb
  PRINT '(A,F8.2)',' msw_infl  : ',msw_infl
  PRINT '(A)'     ,'=================================='
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
  irmse = 0
  rmse_t = 0.0d0
  rmse_x = 0.0d0
  sprd_t = 0.0d0
  sprd_x = 0.0d0
  infl_t = 0.0d0
  infl_x = 0.0d0
  parm_infl(:,1) = ABS(msw_infl)
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
    IF(msw_detailout) THEN
      DO j=1,nwindow
        x4 = xm(:,j)
        WRITE(90) x4
        DO i=1,nbv
          x4 = xf(:,i,j)
          WRITE(92) x4
        END DO
      END DO
    END IF
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
        parm = parm_infl(ix,it+nn-1)
        DO n=1,nwindow
          tdif = REAL(ABS(nn-n),r_size)
          IF(tlocal < 0.0d0) tdif = 0.0d0
          DO i=1,ny
            DO j=1,nx
              IF(MAXVAL(h4d(i,:,n)) == h4d(i,j,n)) EXIT
            END DO
            dist = REAL(MIN(ABS(j-ix),nx-ABS(j-ix)),r_size)
            IF(dist < xmaxloc) THEN
              ny_loc = ny_loc+1
              d_loc(ny_loc) = d(i,n)
              rdiag_loc(ny_loc) = obserr**2
              IF(msw_local == 0) THEN ! fixed localization
                rloc_loc(ny_loc) = EXP(-0.5 * (dist/xlocal)**2) &! space
                  & * EXP(-0.5 * (tdif/tlocal)**2) ! time
              ELSE ! adaptive localization
                CALL com_correl(nbv2,hdxf(i,1:nbv2,n),dxf(ix,1:nbv2,nn),obsloc(1))
                CALL com_correl(nbv-nbv2,hdxf(i,nbv2+1:nbv,n),dxf(ix,nbv2+1:nbv,nn),obsloc(2))
                CALL com_correl(nbv,hdxf(i,:,n),dxf(ix,:,nn),obsloc(3))
                wb = ABS(obsloc(3))
                wa = 1.0d0 - (0.5d0*ABS(obsloc(1)-obsloc(2)))
                rloc_loc(ny_loc) = wa**sa * wb**sb
                IF(msw_local == 2) THEN
                  rloc_loc(ny_loc) = rloc_loc(ny_loc) &
                  & * EXP(-0.5 * (dist/xlocal)**2) &! space
                  & * EXP(-0.5 * (tdif/tlocal)**2)  ! time
                END IF
              END IF
              hdxf_loc(ny_loc,:) = hdxf(i,:,n)
              IF(rloc_loc(ny_loc) < 0.0001d0) ny_loc = ny_loc-1
            END IF
          END DO
        END DO
        CALL letkf_core(nbv,ny*nwindow,ny_loc,hdxf_loc,rdiag_loc,rloc_loc,d_loc,parm,trans)
        IF(msw_infl > 0.0d0) parm = msw_infl
        DO j=1,nbv
          xa(ix,j,nn) = xm(ix,nn)
          DO i=1,nbv
            xa(ix,j,nn) = xa(ix,j,nn) + dxf(ix,i,nn) * trans(i,j)
          END DO
        END DO
        parm_infl(ix,it+nn) = parm
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
    IF(msw_detailout) THEN
      DO n=1,nwindow
        x4 = xm(:,n)
        WRITE(91) x4
        DO i=1,nbv
          x4 = xa(:,i,n)
          WRITE(93) x4
        END DO
      END DO
    END IF
    !
    ! RMSE,SPRD
    !
    DO n=1,nwindow
      rmse_t(it+n-1) = SQRT(SUM((xm(:,n)-xnature(:,it+n-1))**2)/REAL(nx,r_size))
      DO i=1,nx
        sprd_t(it+n-1) = sprd_t(it+n-1) + SUM((xa(i,:,n)-xm(i,n))**2)
      END DO
      sprd_t(it+n-1) = SQRT(sprd_t(it+n-1)/REAL(nx*nbv,r_size))
      infl_t(it+n-1) = SUM(parm_infl(:,it+n-1))/REAL(nx,r_size)
      IF(it > nspinup) THEN
        DO i=1,nx
          rmse_x(i) = rmse_x(i) + (xm(i,n)-xnature(i,it+n-1))**2
          sprd_x(i) = sprd_x(i) + SUM((xa(i,:,n)-xm(i,n))**2)/REAL(nbv,r_size)
          infl_x(i) = infl_x(i) + parm_infl(i,it+n-1)
        END DO
        irmse = irmse + 1
      END IF
    END DO
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

  OPEN(10,FILE='infl.dat',FORM='unformatted')
  DO i=1,nt
    x4 = REAL(parm_infl(:,i),r_sngl)
    WRITE(10) x4
  END DO
  CLOSE(10)

  OPEN(10,FILE='rmse_t.dat',FORM='formatted')
  DO i=1,nt
    WRITE(10,'(3F12.4)') REAL(i-1)/4.0,rmse_t(i),sprd_t(i)
  END DO
  CLOSE(10)

  rmse_x = SQRT(rmse_x / REAL(irmse,r_size))
  sprd_x = SQRT(sprd_x / REAL(irmse,r_size))
  OPEN(10,FILE='rmse_x.dat',FORM='formatted')
  DO i=1,nx
    WRITE(10,'(I4,2F12.4)') i,rmse_x(i),sprd_x(i)
  END DO
  CLOSE(10)

  rmseave = SUM(rmse_t(nspinup+1:nt))/REAL(nt-nspinup,r_size)
  sprdave = SUM(sprd_t(nspinup+1:nt))/REAL(nt-nspinup,r_size)
  inflave = SUM(infl_t(nspinup+1:nt))/REAL(nt-nspinup,r_size)
  PRINT '(A,F12.5)','RMSE = ',rmseave
  PRINT '(A,F12.5)','SPRD = ',sprdave
  PRINT '(A,F12.5)','INFL = ',inflave

  STOP
END PROGRAM letkf
