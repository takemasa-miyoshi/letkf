PROGRAM obsdump
  IMPLICIT NONE
  REAL(4) :: wk(6)
  INTEGER :: ios
  CHARACTER(1) :: S
  OPEN(3,FORM='unformatted')
  DO
    READ(3,IOSTAT=ios) wk
    IF(ios /= 0) THEN
      PRINT '(A)','END OF FILE'
      EXIT
    END IF
    PRINT '(I6,2F7.2,F10.2,2ES12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6)
    PRINT '(A)','PRESS "S" TO STOP'
    READ(5,'(A1)') S
    IF(S == 'S') EXIT
  END DO
  CLOSE(3)
  STOP
END PROGRAM obsdump
