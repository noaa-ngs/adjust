C %P%

      SUBROUTINE SYSTIM (ITYPE, TBUF, DIFM, DIFM0)

********************************************************************************
* Get the system time
* System Specific Subroutine - Sun/Solaris unix operating system
* 
* if ITYPE = 0, the system time is returned in the TBUF char string
*               and the saved times are updated
* if ITYPE = 1, the time diffetence in minutes from the the last
*               ITYPE=0 call is returned in TDIF0 and the time
*               difference in minutes from the last call to GETTIM
*               is returned in TDIF.
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*99 SCCSID
      CHARACTER*26 TBUF, CTIME
      INTEGER TMSEC
      INTEGER TIME

      SAVE RMIN0, RLMIN

C      SCCSID='$Id: sysdep_ux.for 44197 2010-07-13 12:26:03Z bruce.tran $	20$Date: 2007/11/20 15:28:37 $ NGS'

*** Get time since Jan 1, 1970 in numerical format

      TMSEC = TIME()
      RMIN = DBLE(TMSEC) / 60.D0

      IF (ITYPE .EQ. 0) THEN

*** Save current minutes

        RMIN0 = RMIN
        RLMIN = RMIN

*** Convert numeric to ASCII string

*** Now put C string into a FORTRAN string

        TBUF = CTIME(TMSEC)

      ELSEIF (ITYPE .EQ. 1) THEN

*** Get time differences

        DIFM0 = RMIN - RMIN0
        DIFM =  RMIN - RLMIN

*** Update RLMIN

        RLMIN = RMIN

      ELSE

*** Illegel ITYPE

        STOP 'ILLEGAL ITYPE IN SUBROUTINE SYSTIM'

      ENDIF

      RETURN
      END


      SUBROUTINE MYFLSH (IUNIT)
********************************************************************************
* FLUSH THE IUNIT BUFFER
* System Specific Subroutine - Sun/Solaris unix operating system
********************************************************************************
      IMPLICIT NONE
      INTEGER IUNIT
C
      CALL FLUSH (IUNIT)
C
      RETURN
      END

