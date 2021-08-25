*********************************************************************
* Modification History -
* 
*********************************************************************

      SUBROUTINE GET_OFIL( LUNIT, FNAME )

*********************************************************************
* OPENS AN OUTPUT FILE ON A UNIX SYSTEM
* TESTS FOR SUCCESSFUL OPEN, EXITS IF OPEN FAILS
*  in - LUNIT : file allocation number
*  in - FNAME : filename
*********************************************************************
      IMPLICIT NONE
      CHARACTER*99 SCCSID
      CHARACTER*80 FNAME
      INTEGER      LUNIT
      INTEGER      IOS

C      SCCSID='$Id: getofil_ux.for 44197 2010-07-13 12:26:03Z bruce.tran $	20$Date: 2007/11/20 15:23:15 $ NGS'

      OPEN (LUNIT,FILE=FNAME,STATUS='NEW',
     *ACCESS='SEQUENTIAL',FORM='FORMATTED',BLANK='ZERO', 
     *IOSTAT=IOS)

      IF (IOS .NE. 0) THEN
        WRITE (*, 10) FNAME
        STOP
      ENDIF
   10 FORMAT (/, 
     + ' Error: Not able to open output file',/,
     + '       ',A80,/,
     + '        Check for existing file of same name',/,
     + '        Program exits',/)
C
      RETURN
      END

