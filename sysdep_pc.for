C %P%

      SUBROUTINE SYSTIM( ITYPE,  TBUF, DIFM, DIFM0 )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:         SYSTIM
C VERSION:      9003.29
C WRITTEN BY:   M. SCHENEWERK
C PURPOSE:      RETURN SYSTEM TIME  FOR PC OPERATING SYSTEM
C
C INPUT PARAMETERS:
C -----------------
C ITYPE       CONTROL FLAG = 0 = INITIALIZE
C                          = 1 = RELATIVE
C
C OUTPUT PARAMETERS:
C ------------------
C DIFM        TIME SINCE LAST SYSTIM CALL [ MIN ]
C DIFM0       TIME SINCE INITIAL SYSTIM CALL [ MIN ]
C TBUF        ASCII STRING CONTAINING CURRENT TIME
C
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
C
C
C    THIS MODULE CALLED BY:  ADJUST
C
C    THIS MODULE CALLS:      MOD, TIMEDATE
C
C    INCLUDE FILES USED:
C
C    COMMON BLOCKS USED:
C
C    REFERENCES:
C
C    COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C::MODIFICATION HISTORY
C::9003.29, MSS, CREATION
C::9105.20, SAH, CHANGED TIMEDATE TO IDATE & DOSTIM FOR NDP FORTRAN 3.1.
C::9402.15, SAH, MODIFIED SYSTIM() TO RUN ON AN OS/2 PLATFORM.
C********1*********2*********3*********4*********5*********6*********7**
CE::SYSTIM
C
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
      IMPLICIT INTEGER ( I-N )
C
      CHARACTER*21  TLINE
      CHARACTER*3   CMONTH(12)
      CHARACTER*(*) TBUF
      CHARACTER*99 SCCSID
C
C        code for Digital Fortran
	!DEC IF DEFINED(DF)
	CHARACTER *8 DFDATE
	CHARACTER*10 DFTIME
	CHARACTER*5 DFZONE
	INTEGER DFVALUES(8)
	!DEC ENDIF
CC    INTEGER       IDA, IHR, IMI, IMO, ISE, IYR, IMS
      INTEGER       IDA, IHR, IMI, IMO, ISE, IYR, IHUN
      INTEGER*2     YEAR, MONTH, DAY, HRS, MINS, SECS, HSECS
C
      INTEGER*4     MJD
C
C     INCLUDE 'FSUBLIB.FI'
C
      DATA CMONTH / 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
     $    'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /
C
      SAVE CMONTH, TIME0, TIMEL

C      SCCSID='$Id: sysdep_pc.for 44197 2010-07-13 12:26:03Z bruce.tran $	20$Date: 2007/11/20 15:28:10 $ NGS'

C
C.........  1.0  GET SYSTEM TIME
C
C::NDP1.4E
C::NDP1.4E      CALL TIMEDATE(IYR,IMO,IDA,IHR,IMI,ISE,IMS)
C::NDP1.4E      SEC= ISE + IMS*1.D-3
C::NDP1.4E
C::NDP3.1.0
C          CALL IDATE(IMO,IDA,IYR)
C          CALL DOSTIM(IHR,IMI,ISE,IHUN)
C          SEC = DBLE(ISE) + ( DBLE(IHUN)/100.D0 )
C::NDP3.1.0
C
C::OS/2
* 	!DEC IF DEFINED (WATCOM)
           CALL GETDAT(YEAR,MONTH,DAY)
           CALL GETTIM(HRS,MINS,SECS,HSECS)
           SEC = DBLE(SECS) + ( DBLE(HSECS)/100.D0 )
           IYR = YEAR
           IMO = MONTH
           IDA = DAY
           IHR = HRS
           IMI = MINS
*	!DEC ELSE IF DEFINED(DF)
*		CALL DATE_AND_TIME(DFDATE,DFTIME,DFZONE,DFVALUES)
*		IYR=DFVALUES(1)
* 		IMO=DFVALUES(2)
*		IDA=DFVALUES(3)
*		IHR=DFVALUES(5)
*		IMI=DFVALUES(6)
*		SEC=DBLE(DFVALUES(7)) + (DBLE(DFVALUES(8))/1000.0)
*       !DEC ENDIF


C::OS/2
C
C.........  2.0  STORE TIME AS STRING
C
      WRITE( TLINE, 2000 ) MOD( IYR, 100 ), CMONTH( IMO ), IDA,
     $    IHR, IMI, SEC
 2000 FORMAT( I2.2, '/', A3, '/', I2.2, 1X, I2.2, ':', I2.2, ':', F5.2 )
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

*** FLUSH THE IUNIT BUFFER

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER IUNIT

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C !!!! System Specific Subroutine !!!!!!!!

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C  This is the call in the SUN operating system
C
C      CALL FLUSH (IUNIT)
C
C  This is the call on the HP9000's
C
C     ?
C This is the call using Open WATCOM
      CALL FLUSHUNIT(IUUNIT)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
C::MJD
C
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
 
 
C$$DECDAY       LAST MODIFICATION:      31AUG89 BY MSS
C
      REAL*8 FUNCTION DECDAY( IHR, IMI, SEC )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:         DECDAY
C VERSION:      8908.02
C WRITTEN BY:   M. SCHENEWERK
C PURPOSE:      CONVERT TIME TO DECIMAL DAY
C
C INPUT PARAMETERS FROM THE ARGUEMENT LIST:
C -----------------------------------------
C IHR           HOUR
C IMI           MINUTE
C SEC           SECOND
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
C       THIS MODULE CALLS:
C
C       INCLUDE FILES USED:
C
C       COMMON BLOCKS USED:     
C
C       REFERENCES:
C
C       COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C::DECDAY
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
      DECDAY= ((SEC/60.0D0 + IMI)/60.0D0 + IHR)/24.0D0
C
      RETURN
      END
 
 
