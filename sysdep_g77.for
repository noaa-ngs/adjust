C %P%
C System dependant calls used in ADJUST using PC platform, G77 compiler

      SUBROUTINE SYSTIM( ITYPE,  TBUF, DIFM, DIFM0 )

C********1*********2*********3*********4*********5*********6*********7**
C PURPOSE:      RETURN SYSTEM TIME  FOR PC OPERATING SYSTEM
C
C   in - ITYPE : CONTROL FLAG = 0 = INITIALIZE
C                             = 1 = RELATIVE
C   out- DIFM  : TIME SINCE LAST SYSTIM CALL [ MIN ]
C   out- DIFM0 : TIME SINCE INITIAL SYSTIM CALL [ MIN ]
C   out- TBUF  : ASCII STRING CONTAINING CURRENT TIME
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C IDA         DAY OF MONTH
C IHR         HOUR [ 0 - 23 ]
C IMI         MINUTE [ 0 - 59 ]
C IMO         MONTH [ 1 - 12 ]
C IMS         MILLISECONDS [ 0 - 999 ]
C ISE         SECONDS [ INTEGER ]
C IYR         YEAR
C SEC         SECOND [ INTEGER + FRACTION ]
C TIME        TEMPORARY STORAGE OF CURRENT TIME [MJD]
C TLINE       TEMPORARY STORAGE OF TIME AS AN ASCII STRING
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C CMONTH()    MONTH NAMES
C TIME0       INITIAL TIME [MJD]
C TIMEL       LAST TIME [MJD]
C********1*********2*********3*********4*********5*********6*********7**
      IMPLICIT NONE
C
      CHARACTER*3   CMONTH(12)
      CHARACTER*21  TLINE
      CHARACTER*(*) TBUF

      CHARACTER*99 SCCSID
      CHARACTER*8  DATE
      CHARACTER*9  TIMEC
      CHARACTER*5  ZONE
      INTEGER  ITYPE

      REAL*8 DIFM
      REAL*8 DIFM0
      REAL*8 TIME
      REAL*8 TIME0
      REAL*8 TIMEL
      REAL*8  SEC
      REAL*8  DECSEC
      INTEGER ICY, IYR, IMO, IDA, IHR, IMI, ISEC

C Functions
      REAL*8        DECDAY
      INTEGER*4     MJD
C
      DATA CMONTH / 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
     $              'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /

      SAVE CMONTH, TIME0, TIMEL

      INTRINSIC DATE_AND_TIME

C      SCCSID='$Id: sysdep_g77.for 90675 2016-06-09 18:46:41Z bruce.tran $	20$Date: 2007/11/30 11:30:55 $ NGS'

C
C.........  1.0  GET SYSTEM TIME
C
C DATE  := ccyymmdd
C TIMEC := hhmmss.ss
C ZONE  := Shhmm
C VALUES parameter fails under f77, requires f95

      CALL DATE_AND_TIME( DATE, TIMEC, ZONE )

C Parse to components
      READ (DATE(1:2), '(I2)') ICY
      READ (DATE(3:4), '(I2)') IYR
      READ (DATE(5:6), '(I2)') IMO
      READ (DATE(7:8), '(I2)') IDA
      READ (TIMEC(1:2),'(I2)') IHR
      READ (TIMEC(3:4),'(I2)') IMI
      READ (TIMEC(5:6),'(I2)') ISEC
      READ (TIMEC(8:9),'(F2.0)') DECSEC

      SEC = DBLE(ISEC) + DECSEC/100.0

C
C.........  2.0  STORE TIME AS STRING
C
      WRITE( TLINE, 2000 ) IYR, CMONTH( IMO ), IDA, IHR, IMI, SEC
 2000 FORMAT( I2.2,'/',A3,'/',I2.2, 1X, I2.2,':',I2.2,':',F5.2 )

C DEBUG - creates TLINE directly from DATE_AND_TIME returned variables
C      WRITE( TLINE, 2000 ) DATE(3:4), CMONTH(IMO), DATE(7:8),
C     *     TIMEC(1:2), TIMEC(3:4), TIMEC(5:9)
C 2000 FORMAT( A2,'/',A3,'/',A2,1X, A2,':',A2,':',A5 )

      TBUF= TLINE
C
C.........  3.0  STORE TIME AND COMPARE TO LAST TIME
C
      TIME= MJD( IYR, IMO, IDA ) + DECDAY( IHR, IMI, SEC )
      IF( ITYPE.EQ.0 ) THEN
        TIME0= MJD( IYR, IMO, IDA ) + DECDAY( IHR, IMI, SEC )
        TIMEL= TIME
        DIFM= 0
        DIFM0= 0
      ELSE
        DIFM= ( TIME - TIMEL )*1440.D0
        DIFM0= ( TIME - TIME0 )*1440.D0
      ENDIF
        TIMEL= TIME
C
      RETURN
      END


      SUBROUTINE MYFLSH (IUNIT)
********************************************************************************
* FLUSHES THE IUNIT BUFFER
* 
* http://gcc.gnu.org/onlinedocs/gcc-3.4.6/g77/Flush-Intrinsic.html#Flush-Intrinsic
* Flushes Fortran unit(s) currently open for output. 
* Without the optional argument, all such units are flushed, 
* otherwise just the unit specified by Unit.
* 
* Some non-GNU implementations of Fortran provide this intrinsic as a library 
* procedure that might or might not support the (optional) Unit argument. 
********************************************************************************
      IMPLICIT NONE
      INTEGER IUNIT
      INTRINSIC FLUSH

C  Intrinsic G77 subroutine
C
      CALL FLUSH (IUNIT)
      RETURN
      END



C$$MJD  LAST MODIFICATION:      06SEP89 BY MSS
C
      INTEGER*4 FUNCTION MJD( IYR, IMO, IDA )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:         MJD
C VERSION:      8908.02
C WRITTEN BY:   M. SCHENEWERK
C PURPOSE:      CONVERT DATE TO MODIFIED JULIAN DATE
C
C INPUT PARAMETERS FROM THE ARGUEMENT LIST:
C -----------------------------------------
C IYR           YEAR
C IMO           MONTH
C IDA           DAY
C
C OUTPUT PARAMETERS FROM ARGUEMENT LIST:
C --------------------------------------
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C       THIS MODULE CALLED BY:  GENERAL USE
C
C       THIS MODULE CALLS:      DINT
C
C       INCLUDE FILES USED:
C
C       COMMON BLOCKS USED:     
C
C       REFERENCES:             DUFFETT-SMITH, PETER  1982, "PRACTICAL
C                               ASTRONOMY WITH YOUR CALCULATOR", 2ND
C                               EDITION, CAMBRIDGE UNIVERSITY PRESS,
C                               NEW YORK, P.9
C
C       COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
      IMPLICIT REAL*8 (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER*4     A, B, C, D
      COMMON/UNITS/ LUNIT
C
C........  0.0  EXPLICIT INITIALIZATION
C
      IYRP= IYR
      IF( IYR.LT.1900 ) IYRP= IYR + 1900
      IF(IYR.LT.20) IYRP = IYR + 2000
      IMOP= IMO
      IF( IMO .LT. 3 ) THEN
        IYRP= IYRP - 1
        IMOP= IMOP + 12
      END IF
C
C........  1.0  CALCULATION
C
      A=  IYRP * 0.01D0
      B=  2 - A + DINT( A*0.25D0 )
      C=  365.25D0 * IYRP
      D=  30.6001D0 * (IMOP+1)
      MJD=  B + C + D + IDA - 679006
C
*     WRITE(LUNIT,50) IYRP,IMOP,IDA
*  50 FORMAT(' IYRP,IMOP,IDA= ',i4,2x,i2,2x,i2)    
   
*     WRITE(LUNIT,55) A,B,C,D,MJD  
*  55 FORMAT(' A,B,C,D,MJD= ',i4,2x,i4,2x,i7,2x,i4,2x,i12)    

      RETURN
      END
 
 
C
      REAL*8 FUNCTION DECDAY( IHR, IMI, SEC )
C
*********1*********2*********3*********4*********5*********6*********7**
* CONVERT TIME (HOURS, MINUTES, SECONDS) TO DECIMAL DAY
*   in - IHR : HOUR
*   in - IMI : MINUTE
*   in - SEC : SECOND
*   out- nothing
*   ret- DECDAY
*********1*********2*********3*********4*********5*********6*********7**
      IMPLICIT NONE
      REAL*8   SEC
      INTEGER  IHR
      INTEGER  IMI
C
      DECDAY= ((SEC/60.0D0 + IMI)/60.0D0 + IHR)/24.0D0
C
      RETURN
      END

