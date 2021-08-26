C %P%
* file tmisc.for v 4.2.1
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
      CHARACTER*99 SCCSID
      COMMON/UNITS/ LUNIT

C      SCCSID="$Id: tmisc.for 81724 2014-12-22 16:17:53Z jarir.saleh $	20$Date: 2007/11/09 11:19:10 $ NGS"

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
 
 
 
CB::LPSEC
C
      SUBROUTINE LPSEC( TH, TL, UTCTAI, GPSUTC )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:         LPSEC
C VERSION:      8912.27
C WRITTEN BY:   M. SCHENEWERK
C PURPOSE:      DETERMINE UTC - TAI, AND GPS - UTC FOR A DESIRED EPOCH
C
C INPUT PARAMETERS FROM THE ARGUEMENT LIST:
C -----------------------------------------
C TH            DATE [JULIAN OR MODIFIED JULIAN DATE]
C TL            UTC TIME [DECIMAL DAYS]
C
C OUTPUT PARAMETERS FROM ARGUEMENT LIST:
C --------------------------------------
C UTCTAI        DIFFERENCE UTC - TAI
C GPSUTC        DIFFERENCE GPS - UTC
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C TMPH          DATE USED [MJD]
C TMPL          TIME USED [DECIMAL DAY]
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
C::MODIFICATION HISTORY
C::8908.31, MSS, DOC. STD. IMPLIMENTED
C::8912.27, MSS, NEW STD. IMPLIMENTED, TIME->TH&TL, MJD OR JD
C********1*********2*********3*********4*********5*********6*********7**
CE::LPSEC
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
      INTEGER*4 TH, TMPH
C
C........  0.0  EXPLICIT INITIALIZATION - VALUES GOOD FROM 01JUL81
C
	GPSUTC    =  0
	UTCTAI    =-19.D0
C
	IF( TH.GT.2400000 ) THEN
	  TMPH= TH - 2400000
	  TMPL= TL - 0.5D0
	ELSE
	  TMPH= TH
	  TMPL= TL
	ENDIF
	TMPH= TMPH + DINT(TMPL)
	TMPL= TMPL - DINT(TMPL)
	IF( TMPL.LT.0 ) THEN
	  TMPH= TMPH - 1
	  TMPL= TMPL + 1.D0
	ENDIF
C
	IF( TMPH.LT.44239 ) THEN
C         WRITE( *,  100 ) TH+TL
 100      FORMAT( T2, 'LPSEC: ', F15.5, ' BEFORE EARLIEST TIME',
     $        ' = 44239 = 01JAN80' )
	  RETURN
	ENDIF
C
C........  1.0  01JAN82
C
	IF( TMPH.LT.44970 ) RETURN
	GPSUTC= GPSUTC + 1.D0
	UTCTAI= UTCTAI - 1.D0
C
C........  2.0  01JUL82
C
	IF( TMPH.LT.45151 ) RETURN
	GPSUTC= GPSUTC + 1.D0
	UTCTAI= UTCTAI - 1.D0
C
C........  3.0  01JUL83
C
	IF( TMPH.LT.45516 ) RETURN
	GPSUTC= GPSUTC + 1.D0
	UTCTAI= UTCTAI - 1.D0
C
C........  4.0  01JUL85
C
	IF( TMPH.LT.46247 ) RETURN
	GPSUTC= GPSUTC + 1.D0
	UTCTAI= UTCTAI - 1.D0
C
C........  5.0  01JAN88
C
	IF( TMPH.LT.47161 ) RETURN
	GPSUTC= GPSUTC + 1.D0
	UTCTAI= UTCTAI - 1.D0
C
C........  6.0  01JAN90
C
	IF( TMPH.LT.47892 ) RETURN
	GPSUTC= GPSUTC + 1.D0
	UTCTAI= UTCTAI - 1.D0
C
	RETURN
	END
 
 
CB::DMY
C
      SUBROUTINE DMY( TH, TL, IYR, IMO, IDA, IHR, IMI, SEC )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:       DMY
C VERSION:    8912.28
C WRITTEN BY: M. SCHENEWERK
C PURPOSE:    CONVERT DATE TO MODIFIED JULIAN DATE
C
C INPUT PARAMETERS FROM THE ARGUEMENT LIST:
C -----------------------------------------
C TH          INTEGER (HIGH) PART OF DATE - (MODIFIED) JULIAN DATE
C TL          FRACTIONAL (LOW) PART OF DATE - DECIMAL DAY
C
C OUTPUT PARAMETERS FROM ARGUEMENT LIST:
C --------------------------------------
C IYR         YEAR
C IMO         MONTH
C IDA         DAY
C IHR         HOUR
C IMI         MINUTE
C SEC         SECOND
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C       THIS MODULE CALLED BY: GENERAL USE
C
C       THIS MODULE CALLS:
C
C       INCLUDE FILES USED:
C
C       COMMON BLOCKS USED:
C
C       REFERENCES:            DUFFETT-SMITH, PETER  1982, "PRACTICAL
C                              ASTRONOMY WITH YOUR CALCULATOR", 2ND
C                              EDITION, CAMBRIDGE UNIVERSITY PRESS,
C                              NEW YORK, P.9
C
C       COMMENTS:               THIS SUBROUTINE IS LIMITED BY THE SIZE
C                               OF THE A REAL*8 VARIABLE.  ACCURACY OF
C                               THE FRACTION SECONDS VALUE IS LIMITED
C                               TO APPROXIMATELY 1 MICROSECOND.
C
C********1*********2*********3*********4*********5*********6*********7**
C::MODIFICATION HISTORY
C::8908.31, MSS, DOC. STD. IMPLIMENTED
C::8912.28, MSS, HANDLES TL<0
C********1*********2*********3*********4*********5*********6*********7**
CE::DMY
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
      INTEGER*4 A, B, C, D, E, G, I, TH
C
C........  0.0  EXPLICIT INITIALIZATION
C
      I= TH
      IF( TH.LT.2400001 ) THEN
	I=  TH + 2400001
      ELSE
	I= TH + 1
      ENDIF
      F= TL
      I= I+DINT(TL)
      F= F-DINT(TL)
      IF( F.LT.0 ) THEN
	I= I - 1
	F= F + 1.D0
      ENDIF
C
C........  1.0  CONVERT DATE
C
      A=  I
      IF( I .GT. 2299160 ) A = (I-1867216.25D0)/36524.25D0
      B=  I + 1 + A - DINT( A*0.25D0 )
      C=  B + 1524
      D=  (C-122.1D0)/365.25D0
      E=  365.25D0*D
      G=  (C-E)/30.6001D0
C
      I=  30.6001D0 * G
      IDA=  DINT(C - E + F - I)
C
      IMO=  G - 1.0D0
      IF( G .GT. 13.5D0 ) IMO = G - 13.0D0
C
      IYR=  D - 6615.0D0
      IF( IMO .GT. 2.5D0 ) IYR = D - 6616.0D0
C
C........  2.0  CONVERT TIME
C
      TEMP=  F * 24.0D0
      IHR=  TEMP
      TEMP=  (TEMP-IHR)*60.0D0
      IMI=  TEMP
      SEC=  (TEMP-IMI)*60.0D0
C
      RETURN
      END
 
C$$IGPSWK
C
      INTEGER FUNCTION IGPSWK( TIMEH, TIMEL, SECWK )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:         IGPSWK
C VERSION:      IGPS03
C WRITTEN BY:   M. SCHENEWERK
C PURPOSE:      CONVERT MODIFIED JULIAN DATE TO GPS WEEK AND SECONDS-
C               OF-WEEK
C
C INPUT PARAMETERS FROM THE ARGUEMENT LIST:
C -----------------------------------------
C TIMEH         INTEGRAL (HIGH) PART OF DATE (MODIFIED JULIAN DATE)
C TIMEL         FRACTIONAL (LOW) PART OF DATE (DECIMAL DAYS)
C
C OUTPUT PARAMETERS FROM ARGUEMENT LIST:
C --------------------------------------
C SECWK         SECONDS INTO THE WEEK
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
C                               NEW YORK, P.12
C
C       COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C::IGPSWK
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
      INTEGER*4 TIMEH
C
      DATA D2S/ 8.64D4 /, W2D/ 7.D0 /, WKOFF/ 349178.D0 /
      DATA P_J2M/ 2400002.D0 /
C
      IGPSWK=  DINT( (TIMEH + P_J2M)/W2D - WKOFF )
      SECWK=  (TIMEH - ((IGPSWK + WKOFF)*W2D - P_J2M) + TIMEL)*D2S
C
      RETURN
      END
 
C$$NORMJD       LAST MODIFICATION:      31AUG89 BY MSS
C
      SUBROUTINE NORMJD( TIMEH, TIMEL )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:         NORMJD
C VERSION:      8908.02
C WRITTEN BY:   M. SCHENEWERK
C PURPOSE:      NORMALIZE (MODIFIED) JULIAN DATE
C
C INPUT PARAMETERS FROM THE ARGUEMENT LIST:
C -----------------------------------------
C TIMEH         INTEGER (HIGH) PORTION OF DATE
C TIMEL         FRACTIONAL (LOW) PORTION OF DATE
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
C       REFERENCES:
C
C       COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C::NORMJD
C
      INTEGER*4 TIMEH
      REAL*8 TIMEL
C
      TIMEH= TIMEH + DINT(TIMEL)
      TIMEL= TIMEL - DINT(TIMEL)
C
      RETURN
      END

CB::TGT
C
      LOGICAL FUNCTION TGT( TH1, TL1, TH2, TL2 )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:         TGT
C VERSION:      TGTX.01
C WRITTEN BY:   M. SCHENEWERK
C PURPOSE:      DETERMINE IF TIME 1 IS GREATER THAN TIME 2
C
C INPUT PARAMETERS:
C -----------------
C TH1      HIGH PART OF TIME 1 [MJD]                            (I*4)
C TL1      LOW PART OF TIME 1 [DECIMAL DAY]                     (D)
C TH2      HIGH PART OF TIME 2 [MJD]                            (I*4)
C TL2      LOW PART OF TIME 2 [DECIMAL DAY]                     (D)
C
C OUTPUT PARAMETERS:
C ------------------
C TGT      IF TIME 1 > TIME 2 = .TRUE., ELSE = .FALSE.          (L)
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C    THIS MODULE CALLED BY:  GENERAL USE IN PAGE AND PREP
C
C    THIS MODULE CALLS:      
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
C    MODIFICATION HISTORY:
C::8910.18, MSS, CREATION
C********1*********2*********3*********4*********5*********6*********7**
CE::TGT
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DOUBLE PRECISION TL1, TL2
      INTEGER*4        TH1, TH2
C
      IF( (TH1.GT.TH2) .OR. (TH1.EQ.TH2 .AND. TL1.GT.TL2) ) THEN
	TGT= .TRUE.
      ELSE
	TGT= .FALSE.
      ENDIF
C
      RETURN
      END

CB::TGE
C
      LOGICAL FUNCTION TGE( TH1, TL1, TH2, TL2 )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:         TGE
C VERSION:      TGEX.01
C WRITTEN BY:   M. SCHENEWERK
C PURPOSE:      DETERMINE IF TIME 1 IS GREATER THAN OR EQUAL TO TIME 2
C
C INPUT PARAMETERS:
C -----------------
C TH1      HIGH PART OF TIME 1 [MJD]                            (I*4)
C TL1      LOW PART OF TIME 1 [DECIMAL DAY]                     (D)
C TH2      HIGH PART OF TIME 2 [MJD]                            (I*4)
C TL2      LOW PART OF TIME 2 [DECIMAL DAY]                     (D)
C
C OUTPUT PARAMETERS:
C ------------------
C TGE      IF TIME 1 >= TIME 2 = .TRUE., ELSE = .FALSE.          (L)
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C    THIS MODULE CALLED BY:  GENERAL USE IN PAGE AND PREP
C
C    THIS MODULE CALLS:      
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
C    MODIFICATION HISTORY:
C::8910.18, MSS, CREATION
C********1*********2*********3*********4*********5*********6*********7**
CE::TGE
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DOUBLE PRECISION TL1, TL2
      INTEGER*4        TH1, TH2
C
      IF( (TH1.GT.TH2) .OR. (TH1.EQ.TH2 .AND. TL1.GE.TL2) ) THEN
	TGE= .TRUE.
      ELSE
	TGE= .FALSE.
      ENDIF
C
      RETURN
      END

CB::TEQ
C
      LOGICAL FUNCTION TEQ( TH1, TL1, TH2, TL2 )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:         TEQ
C VERSION:      TEQX.01
C WRITTEN BY:   M. SCHENEWERK
C PURPOSE:      DETERMINE IF TIME 1 IS EQUAL TO TIME 2
C
C INPUT PARAMETERS:
C -----------------
C TH1      HIGH PART OF TIME 1 [MJD]                            (I*4)
C TL1      LOW PART OF TIME 1 [DECIMAL DAY]                     (D)
C TH2      HIGH PART OF TIME 2 [MJD]                            (I*4)
C TL2      LOW PART OF TIME 2 [DECIMAL DAY]                     (D)
C
C OUTPUT PARAMETERS:
C ------------------
C TEQ      IF TIME 1 = TIME 2 = .TRUE., ELSE = .FALSE.          (L)
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C    THIS MODULE CALLED BY:  GENERAL USE IN PAGE AND PREP
C
C    THIS MODULE CALLS:      
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
C    MODIFICATION HISTORY:
C::8910.18, MSS, CREATION
C********1*********2*********3*********4*********5*********6*********7**
CE::TEQ
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DOUBLE PRECISION TL1, TL2
      INTEGER*4        TH1, TH2
C
      IF( (TH1.EQ.TH2) .AND. (TL1.EQ.TL2) ) THEN
	TEQ= .TRUE.
      ELSE
	TEQ= .FALSE.
      ENDIF
C
      RETURN
      END

CB::TNE
C
      LOGICAL FUNCTION TNE( TH1, TL1, TH2, TL2 )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:         TNE
C VERSION:      TNEX.01
C WRITTEN BY:   M. SCHENEWERK
C PURPOSE:      DETERMINE IF TIME 1 IS EQUAL TO TIME 2
C
C INPUT PARAMETERS:
C -----------------
C TH1      HIGH PART OF TIME 1 [MJD]                            (I*4)
C TL1      LOW PART OF TIME 1 [DECIMAL DAY]                     (D)
C TH2      HIGH PART OF TIME 2 [MJD]                            (I*4)
C TL2      LOW PART OF TIME 2 [DECIMAL DAY]                     (D)
C
C OUTPUT PARAMETERS:
C ------------------
C TNE      IF TIME 1 < TIME 2 OR TIME1 < TIME 2 = .TRUE.,       (L)
C          ELSE = .FALSE.
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C    THIS MODULE CALLED BY:  GENERAL USE IN PAGE AND PREP
C
C    THIS MODULE CALLS:      
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
C    MODIFICATION HISTORY:
C::8910.18, MSS, CREATION
C********1*********2*********3*********4*********5*********6*********7**
CE::TNE
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DOUBLE PRECISION TL1, TL2
      INTEGER*4        TH1, TH2
C
      IF( (TH1.NE.TH2) .OR. (TL1.NE.TL2) ) THEN
	TNE= .TRUE.
      ELSE
	TNE= .FALSE.
      ENDIF
C
      RETURN
      END

CB::TLE
C
      LOGICAL FUNCTION TLE( TH1, TL1, TH2, TL2 )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:         TLE
C VERSION:      TLEX.01
C WRITTEN BY:   M. SCHENEWERK
C PURPOSE:      DETERMINE IF TIME 1 IS LESS THAN OR EQUAL TO TIME 2
C
C INPUT PARAMETERS:
C -----------------
C TH1      HIGH PART OF TIME 1 [MJD]                            (I*4)
C TL1      LOW PART OF TIME 1 [DECIMAL DAY]                     (D)
C TH2      HIGH PART OF TIME 2 [MJD]                            (I*4)
C TL2      LOW PART OF TIME 2 [DECIMAL DAY]                     (D)
C
C OUTPUT PARAMETERS:
C ------------------
C TLE      IF TIME 1 >= TIME 2 = .TRUE., ELSE = .FALSE          (L)
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C    THIS MODULE CALLED BY:  GENERAL USE IN PAGE AND PREP
C
C    THIS MODULE CALLS:      
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
C    MODIFICATION HISTORY:
C::8910.18, MSS, CREATION
C********1*********2*********3*********4*********5*********6*********7**
CE::TLE
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DOUBLE PRECISION TL1, TL2
      INTEGER*4        TH1, TH2
C
      IF( (TH1.LT.TH2) .OR. (TH1.EQ.TH2 .AND. TL1.LE.TL2) ) THEN
	TLE= .TRUE.
      ELSE
	TLE= .FALSE.
      ENDIF
C
      RETURN
      END

CB::TLT
C
      LOGICAL FUNCTION TLT( TH1, TL1, TH2, TL2 )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:         TLT
C VERSION:      TLTX.01
C WRITTEN BY:   M. SCHENEWERK
C PURPOSE:      DETERMINE IF TIME 1 IS LESS THAN TIME 2
C
C INPUT PARAMETERS:
C -----------------
C TH1      HIGH PART OF TIME 1 [MJD]                            (I*4)
C TL1      LOW PART OF TIME 1 [DECIMAL DAY]                     (D)
C TH2      HIGH PART OF TIME 2 [MJD]                            (I*4)
C TL2      LOW PART OF TIME 2 [DECIMAL DAY]                     (D)
C
C OUTPUT PARAMETERS:
C ------------------
C TLT      IF TIME 1 > TIME 2 = .TRUE., ELSE = .FALSE           (L)
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C    THIS MODULE CALLED BY:  GENERAL USE IN PAGE AND PREP
C
C    THIS MODULE CALLS:      
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
C    MODIFICATION HISTORY:
C::8910.18, MSS, CREATION
C********1*********2*********3*********4*********5*********6*********7**
CE::TLT
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DOUBLE PRECISION TL1, TL2
      INTEGER*4        TH1, TH2
C
      IF( (TH1.LT.TH2) .OR. (TH1.EQ.TH2 .AND. TL1.LT.TL2) ) THEN
	TLT= .TRUE.
      ELSE
	TLT= .FALSE.
      ENDIF
C
      RETURN
      END
