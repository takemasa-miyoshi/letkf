PROGRAM station
  USE common
  USE common_obs
  USE common_roms

  IMPLICIT NONE

  INTEGER,PARAMETER :: msw_obsnet=1 !1:regular, 2:custom
  INTEGER,PARAMETER :: obkind=5
  INTEGER,PARAMETER :: id_hrf   =1
  INTEGER,PARAMETER :: id_prof  =2
  INTEGER,PARAMETER :: id_sst   =3
  INTEGER,PARAMETER :: id_ssh   =4
  INTEGER,PARAMETER :: id_glider=5
  INTEGER :: i,j,k,n,ios
  INTEGER :: nobkind
  LOGICAL :: ex
  INTEGER :: nu,nv,nt,ns,nz
  REAL(r_size) :: a
  REAL(r_size) :: oberr(obkind)
  INTEGER :: obelm(obkind)
  CHARACTER(100) :: cdummy
  !
  ! read table obserr.tbl
  !
  OPEN(11,FILE='obserr.tbl')
  READ(11,'(A)') cdummy
  READ(11,'(A)') cdummy
  nobkind = 0
  DO
    READ(11,'(L1,X,I5,X,ES9.2)',IOSTAT=ios) ex,i,a
    IF(ios /= 0) EXIT
    IF(ex) THEN
      nobkind = nobkind+1
      obelm(nobkind) = i
      oberr(nobkind) = a
    END IF
  END DO
  CLOSE(11)
  !
  ! find n
  !
  DO nu=1,nobkind
    IF(obelm(nu) == id_u_obs) EXIT
  END DO
  DO nv=1,nobkind
    IF(obelm(nv) == id_v_obs) EXIT
  END DO
  DO nt=1,nobkind
    IF(obelm(nt) == id_t_obs) EXIT
  END DO
  DO ns=1,nobkind
    IF(obelm(ns) == id_s_obs) EXIT
  END DO
  DO nz=1,nobkind
    IF(obelm(nz) == id_z_obs) EXIT
  END DO
  !
  ! stdout
  !
  PRINT '(A)','KIND   I   J  K  ELEM     OBER'
  PRINT '(A)','------------------------------'
  SELECT CASE(msw_obsnet)
  CASE(1)
    CALL regular
  CASE(2)
    CALL custom
  END SELECT

  STOP
CONTAINS
!
! regular observing network
!
SUBROUTINE regular
  INTEGER,PARAMETER :: imin=3
  INTEGER,PARAMETER :: imax=nlon-3-2
  INTEGER,PARAMETER :: iint=5
  INTEGER,PARAMETER :: jmin=3
  INTEGER,PARAMETER :: jmax=nlat-2
  INTEGER,PARAMETER :: jint=5
  INTEGER,PARAMETER :: kmin=10
  INTEGER,PARAMETER :: kmax=40
  INTEGER,PARAMETER :: kint=3
  DO j=jmin,jmax,jint
    DO i=imin,imax,iint
      PRINT '(3I4,I3,I6,ES9.2)',id_hrf,i,j,nlev,obelm(nu),oberr(nu)
      PRINT '(3I4,I3,I6,ES9.2)',id_hrf,i,j,nlev,obelm(nv),oberr(nv)
      PRINT '(3I4,I3,I6,ES9.2)',id_ssh,i,j,nlev,obelm(nz),oberr(nz)
      DO k=kmin,kmax,kint
        PRINT '(3I4,I3,I6,ES9.2)',id_prof,i,j,k,obelm(nt),oberr(nt)
        PRINT '(3I4,I3,I6,ES9.2)',id_prof,i,j,k,obelm(ns),oberr(ns)
        IF(k /= nlev) THEN                                            ! UV
          PRINT '(3I4,I3,I6,ES9.2)',id_prof,i,j,k,obelm(nu),oberr(nu) ! UV
          PRINT '(3I4,I3,I6,ES9.2)',id_prof,i,j,k,obelm(nv),oberr(nv) ! UV
        END IF                                                        ! UV
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE regular
!
! custom observing network
!
SUBROUTINE custom
  INTEGER,PARAMETER :: imin=2
  INTEGER,PARAMETER :: imax=nlon-3-1
  INTEGER,PARAMETER :: iint=2
  INTEGER,PARAMETER :: jmin=2
  INTEGER,PARAMETER :: jmax=nlat-1
  INTEGER,PARAMETER :: jint=2
  INTEGER,PARAMETER :: nsounding=0
  INTEGER,PARAMETER :: nglider=0
  INTEGER :: nn
  REAL(r_size),ALLOCATABLE :: rnd(:)
  !
  ! SST
  !
  DO j=jmin,jmax,jint
    DO i=imin,imax,iint
      PRINT '(3I4,I3,I6,ES9.2)',id_sst,i,j,nlev,obelm(nt),oberr(nt)
    END DO
  END DO
  !
  ! SSH
  !
  DO j=jmin,jmax,jint
    DO i=imin,imax,iint
      PRINT '(3I4,I3,I6,ES9.2)',id_ssh,i,j,nlev,obelm(nz),oberr(nz)
    END DO
  END DO
  !
  ! surface flow (only near coast)
  !
  DO j=jmin,jmax
    DO i=imax-30,imax
      PRINT '(3I4,I3,I6,ES9.2)',id_hrf,i,j,nlev,obelm(nu),oberr(nu)
      PRINT '(3I4,I3,I6,ES9.2)',id_hrf,i,j,nlev,obelm(nv),oberr(nv)
    END DO
  END DO
  !
  ! soundings (random locations, deep)
  !
  IF(nsounding /= 0) THEN
    n = nsounding*2
    ALLOCATE(rnd(n))
    CALL com_rand(n,rnd)
    DO nn=1,nsounding
      i = CEILING(rnd(nn*2-1)*(imax-2)+1)
      j = CEILING(rnd(nn*2)*(jmax-2)+1)
      DO k=1,nlev-1
        PRINT '(3I4,I3,I6,ES9.2)',id_prof,i,j,k,obelm(nt),oberr(nt)
        PRINT '(3I4,I3,I6,ES9.2)',id_prof,i,j,k,obelm(ns),oberr(ns)
        PRINT '(3I4,I3,I6,ES9.2)',id_prof,i,j,k,obelm(nu),oberr(nu)
        PRINT '(3I4,I3,I6,ES9.2)',id_prof,i,j,k,obelm(nv),oberr(nv)
      END DO
    END DO
    DEALLOCATE(rnd)
  END IF
  !
  ! gliders (random line-like in x direction near coast, shallow)
  !
  IF(nglider /= 0) THEN
    ALLOCATE(rnd(nglider))
    CALL com_rand(nglider,rnd)
    DO nn=1,nglider
      j = CEILING(rnd(nn)*(jmax-2)+1)
      DO k=20,nlev,2
        DO i=imax-30,imax
          PRINT '(3I4,I3,I6,ES9.2)',id_glider,i,j,k,obelm(nt),oberr(nt)
          PRINT '(3I4,I3,I6,ES9.2)',id_glider,i,j,k,obelm(ns),oberr(ns)
        END DO
      END DO
    END DO
    DEALLOCATE(rnd)
  END IF

  RETURN
END SUBROUTINE custom
END PROGRAM station
