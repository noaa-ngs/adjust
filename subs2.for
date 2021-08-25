C %P%
C file $RCSfile: subs2.for,v $

      SUBROUTINE HORDIR (BCARD,IUO,IOBS,B,NX,FATAL,LSN)

*** OBSERVATION EQUATIONS FOR HORIZONTAL DIRECTIONS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N )
      PARAMETER ( LENC = 10 )
      CHARACTER*99 SCCSID
      CHARACTER*80 BCARD
      CHARACTER*4 ASS
      CHARACTER*1 TC
      LOGICAL FATAL,LSN
      LOGICAL GETSSN,GETZ,GETIVF
      LOGICAL ADDCON
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL REJECT,GET82
      DIMENSION IC(LENC),C(LENC)
      DIMENSION B(*),NX(*)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /IZS/    IZ2
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON/UNITS/ LUNIT
      SAVE ITIME, IYR, IMO, IDY

C      SCCSID='$Id: subs2.for 83318 2015-05-28 15:50:02Z jarir.saleh $	20$Date: 2008/10/27 17:18:38 $ NGS'

      READ (BCARD,1) IRT,ISSN,LIST,JSSN,ID,IM,ASS
 1    FORMAT (7X,I2, 1X,I4,    I2,34X,I4, 9X,I3,I2,A4)
      CALL NBLANK (ASS,4,IBLK)
      READ (ASS,3) SS
 3    FORMAT (F4.2)
      IF (IRT .EQ. 20) THEN
        READ (BCARD,5) IYR,IMO,IDY,IHR,IMN,TC
 5      FORMAT (39X,5I2,A1)
        IF (BCARD(40:41) .EQ. '  ') IYR = 84
        IF (BCARD(42:43) .EQ. '  ') IMO = 1
        IF (BCARD(44:45) .EQ. '  ') IDY = 1
        IF (BCARD(46:47) .EQ. '  ') IHR = 0
        IF (BCARD(48:49) .EQ. '  ') IMN = 0
        IF (BCARD(50:50) .EQ. ' ') TC = 'Z'
	CALL GETYR(BCARD,IYR)
*       IF (IYR .LT. 14) THEN
*         IYR = IYR+2000
*       ELSE
*         IYR = IYR+1900
*       ENDIF
        CALL TOMNT (IYR,IMO,IDY,IHR,IMN,TC,ITIME)
      ELSEIF (IRT .EQ. 22) THEN
        READ (BCARD,6) IHR,IMN,TC
 6      FORMAT (45X,2I2,A1)
        IF (BCARD(46:47) .EQ. '  ') IHR = 0
        IF (BCARD(48:49) .EQ. '  ') IMN = 0
        IF (BCARD(50:50) .EQ. ' ') TC = 'Z'
        CALL TOMNT (IYR,IMO,IDY,IHR,IMN,TC,ITIME)
      ENDIF

      IF ( BCARD( 6: 6) .EQ. 'R'  .OR.  BCARD( 6: 6) .EQ. 'O'  .OR.
     &     BCARD( 6: 6) .EQ. 'F' ) THEN
        NREJ = NREJ+1
        REJECT = .TRUE.
      ELSE
        REJECT = .FALSE.
      ENDIF

      IF (GET82(ISSN,I)) RETURN
      IF (GET82(JSSN,I)) RETURN

        IF ( .NOT. GETSSN(ISSN,ISN)) THEN
          CALL LINE (3)
          WRITE (LUNIT,2) BCARD
 2        FORMAT ('0NO *80* RECORD FOR--',A80/)
        ELSEIF ( .NOT. GETSSN(JSSN,JSN)) THEN
          CALL LINE (3)
          WRITE (LUNIT,2) BCARD
        ELSE

*** RETRIEVE THE STD DEV

          CALL STDDEV (BCARD,IRT,SD,REJECT)

*** KIND = 11 HORIZONTAL DIRECTION

          IF ( .NOT. GETZ(ISSN,LIST,IZ)) RETURN
          OBSB = (ID+IM/60.D0+SS/3600.D0)/RAD
          IF (IZ .NE. IZ2) CALL STOREZ (ISN,JSN,B,IZ,OBSB,IGRT)
          KIND = 11
          IOBS = IOBS+1
          NOBS = NOBS+1
          NDIR = NDIR+1
          LSN = .TRUE.
          IF ( .NOT. GETIVF(20,ITIME,IVF)) IVF = 0
          IGRT = 0
          CALL FORMIC (KIND,ISN,JSN,IZ,IC,LENG,IGRT)
          CALL FORMC (KIND,C,B,ISN,JSN,IZ,IGRT)
          CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IZ,IGRT)
          CMO = OBS0 - OBSB
          IF (CMO .GT. PI) THEN
            CMO = CMO - PI - PI
          ELSEIF (CMO .LT. -PI) THEN
            CMO = CMO + PI + PI
          ENDIF
          IF (IMODE .EQ. 0) CMO = 0.D0
          CALL BIGV (KIND,ISN,JSN,IOBS,IVF,CMO,SD,FATAL)

*** UPDATE CONNECTIVITY

          IF ( .NOT. ADDCON(IC,LENG,NX)) THEN
            WRITE (LUNIT,667)
 667        FORMAT ('0INSUFFICIENT STORAGE FOR HOR. DIRECTIONS'/)
             write (lunit,*) 'subs2, 7'
            CALL ABORT2
          ENDIF

          WRITE (IUO) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,IVF,
     &                IZ,IGRT

        ENDIF

      RETURN
      END
      SUBROUTINE IJOBST

*** INITIALIZE JOB STATISTICS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP

      N84   = 0
      N85   = 0
      N86   = 0
      NDIR  = 0
      NANG  = 0
      NGPS  = 0
      NZD   = 0
      NDS   = 0
      NAZ   = 0
      NQQ   = 0
      NREJ  = 0
      NGPSR = 0
      NDOP  = 0

      RETURN
      END

C---------------------------------------------------------------------------------------------
      SUBROUTINE INIT83(T)
 
*** INITIALIZE GPS VECTOR TRANSFORMATIONS TO NAD83
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
 
      COMMON /CONST/  PI, PI2, RAD
      COMMON /XLATE/  RX(30), RY(30), RZ(30), SCALE(30),
     +             DELX(30), DELY(30), DELZ(30), COSRX(30), SINRX(30),
     +             COSRY(30), SINRY(30), COSRZ(30), SINRZ(30)
     
      COMMON/UNITS/ LUNIT     
 
** 07-07-05 this is now passed as a parameter      T  = 1997.0D0
      AM = 4.84813681D0 * 1.0D-9
      AR = 1.0D-12
      
*** GPS GFILE CODE '01' WGS 72 TO NAD83
 
      RX(1) = 0.0D0/3600.D0/RAD
      COSRX(1) = DCOS(RX(1))
      SINRX(1) = DSIN(RX(1))
      RY(1) = 0.0D0/3600.D0/RAD
      COSRY(1) = DCOS(RY(1))
      SINRY(1) = DSIN(RY(1))
      RZ(1) = -0.554D0/3600.D0/RAD
      COSRZ(1) = DCOS(RZ(1))
      SINRZ(1) = DSIN(RZ(1))
      SCALE(1) = +0.2263D0/1.0D6
      DELX(1) = 0.0D0
      DELY(1) = 0.0D0
      DELZ(1) = +4.5D0
 
*** GPS GFILE CODE '02' WGS 84 TO NAD83

      RX(2) = 0.0D0
      COSRX(2) = 1.0D0
      SINRX(2) = 0.0D0
      RY(2) = 0.0D0
      COSRY(2) = 1.0D0
      SINRY(2) = 0.0D0
      RZ(2) = 0.0D0
      COSRZ(2) = 1.0D0
      SINRZ(2) = 0.0D0
      SCALE(2) = 0.0D0
      DELX(2) = 0.0D0
      DELY(2) = 0.0D0
      DELZ(2) = 0.0D0

*** GPS GFILE CODE '03' WGS 72 TO NAD83
 
      RX(3) = 0.0D0/3600.D0/RAD
      COSRX(3) = DCOS(RX(3))
      SINRX(3) = DSIN(RX(3))
      RY(3) = 0.0D0/3600.D0/RAD
      COSRY(3) = DCOS(RY(3))
      SINRY(3) = DSIN(RY(3))
      RZ(3) = -0.554D0/3600.D0/RAD
      COSRZ(3) = DCOS(RZ(3))
      SINRZ(3) = DSIN(RZ(3))
      SCALE(3) = +0.2263D0/1.0D6
      DELX(3) = 0.0D0
      DELY(3) = 0.0D0
      DELZ(3) = +4.5D0

*** GPS GFILE CODE '04' WGS 84 TO NAD83

      RX(4) = 0.D0
      COSRX(4) = 1.D0
      SINRX(4) = 0.D0
      RY(4) = 0.D0
      COSRY(4) = 1.D0
      SINRY(4) = 0.D0
      RZ(4) = 0.D0
      COSRZ(4) = 1.D0
      SINRZ(4) = 0.D0
      SCALE(4) = 0.D0
      DELX(4) = 0.D0
      DELY(4) = 0.D0
      DELZ(4) = 0.D0

*** GPS GFILE CODE '05' ITRF(VLBI)1989 TO NAD83
 
      RX(5) = 0.0275D0/3600.D0/RAD
      COSRX(5) = DCOS(RX(5))
      SINRX(5) = DSIN(RX(5))
      RY(5) = 0.0155D0/3600.D0/RAD
      COSRY(5) = DCOS(RY(5))
      SINRY(5) = DSIN(RY(5))
      RZ(5) = 0.0107D0/3600.D0/RAD
      COSRZ(5) = DCOS(RZ(5))
      SINRZ(5) = DSIN(RZ(5))
      SCALE(5) = 0.0D0/1.0D6
      DELX(5) = 0.9191D0
      DELY(5) = -2.0182D0
      DELZ(5) = -0.4835D0
 
*** GPS GFILE CODE '06' NEOS 1990 (PROVISIONAL) TO NAD83
 
      RX(6) = 0.0229D0/3600.D0/RAD
      COSRX(6) = DCOS(RX(6))
      SINRX(6) = DSIN(RX(6))
      RY(6) = 0.0249D0/3600.D0/RAD
      COSRY(6) = DCOS(RY(6))
      SINRY(6) = DSIN(RY(6))
      RZ(6) = 0.0099D0/3600.D0/RAD
      COSRZ(6) = DCOS(RZ(6))
      SINRZ(6) = DSIN(RZ(6))
      SCALE(6) = 0.0D0/1.0D6
      DELX(6) = 0.8845D0
      DELY(6) = -2.0399D0
      DELZ(6) = -0.4835D0

*** GPS GFILE CODE '12' ITRF93 EPOCH 1995.0 TO NAD83
 
      RX(12) = 0.0264D0/3600.D0/RAD
      COSRX(12) = DCOS(RX(12))
      SINRX(12) = DSIN(RX(12))
      RY(12) = 0.0101D0/3600.D0/RAD
      COSRY(12) = DCOS(RY(12))
      SINRY(12) = DSIN(RY(12))
      RZ(12) = 0.0103D0/3600.D0/RAD
      COSRZ(12) = DCOS(RZ(12))
      SINRZ(12) = DSIN(RZ(12))
      SCALE(12) = 0.0D0/1.0D6
      DELX(12) = 0.9769D0
      DELY(12) = -1.9392D0
      DELZ(12) = -0.5461D0
      
*     write(lunit,886) rx(12),cosrx(12),sinrx(12),ry(12),cosry(12),
*    * sinry(12),rz(12),cosrz(12),sinrz(12),scale(12),delx(12),
*    * dely(12),delz(12)
     
*886   format(' code 12 ',/13('  ',f13.10/))     
 
*** GPS GFILE CODE '15' ITRF94 EPOCH 1996.0 TO NAD83
 
      RX(15) = 0.0275D0/3600.D0/RAD
      COSRX(15) = DCOS(RX(15))
      SINRX(15) = DSIN(RX(15))
      RY(15) = 0.0101D0/3600.D0/RAD
      COSRY(15) = DCOS(RY(15))
      SINRY(15) = DSIN(RY(15))
      RZ(15) = 0.0114D0/3600.D0/RAD
      COSRZ(15) = DCOS(RZ(15))
      SINRZ(15) = DSIN(RZ(15))
      SCALE(15) = 0.0D0/1.0D6
      DELX(15) = 0.9738D0
      DELY(15) = -1.9453D0
      DELZ(15) = -0.5486D0

*     write(lunit,887) rx(15),cosrx(15),sinrx(15),ry(15),cosry(15),
*    * sinry(15),rz(15),cosrz(15),sinrz(15),scale(15),delx(15),
*    * dely(15),delz(15)
     
*887   format(' code 15 ',/13('  ',f13.10/))     


*** GPS GFILE CODE '18' ITRF96 EPOCH 1997.0 TO NAD83

      RX(18) = (125033.0D0 + (258.0D0 * (T - 1997.0D0))) * AR 
      COSRX(18) = DCOS(RX(18))
      SINRX(18) = DSIN(RX(18))
      RY(18) =  (46785.0D0 - (3599.0D0 * (T - 1997.0D0))) * AR
      COSRY(18) = DCOS(RY(18))
      SINRY(18) = DSIN(RY(18))
      RZ(18) = (56529.0D0 - (153.0D0 * (T - 1997.0D0))) * AR
      COSRZ(18) = DCOS(RZ(18))
      SINRZ(18) = DSIN(RZ(18))
      SCALE(18) = 0.0D0/1.0D6
      DELX(18) = 0.9910D0 
      DELY(18) = -1.9072D0 
      DELZ(18) = -0.5129D0
      
*     write(lunit,888) rx(18),cosrx(18),sinrx(18),ry(18),cosry(18),
*    * sinry(18),rz(18),cosrz(18),sinrz(18),scale(18),delx(18),
*    * dely(18),delz(18)
     
*888   format(' code 18 ',/13('  ',f13.10/))     


*** GPS GFILE CODE '19' ITRF97 EPOCH 1997.0 TO NAD83
 
*     RX(19) = 25.915D0*4.84813681D0*1.0D-9
*     COSRX(19) = DCOS(RX(19))
*     SINRX(19) = DSIN(RX(19))
*     RY(19) =  9.426D0*4.84813681D0*1.0D-9
*     COSRY(19) = DCOS(RY(19))
*     SINRY(19) = DSIN(RY(19))
*     RZ(19) = 11.599D0*4.84813681D0*1.0D-9
*     COSRZ(19) = DCOS(RZ(19))
*     SINRZ(19) = DSIN(RZ(19))
*     SCALE(19) = -0.930D0*1.0D-9
*     DELX(19) = 0.9889D0
*     DELY(19) = -1.9074D0
*     DELZ(19) = -0.5030D0

      RX(19) = (25.915D0 + (0.067D0 * (T - 1997.0D0))) * AM 
      COSRX(19) = DCOS(RX(19))
      SINRX(19) = DSIN(RX(19))
      RY(19) =  (9.426D0 - (0.757D0 * (T - 1997.0D0))) * AM 
      COSRY(19) = DCOS(RY(19))
      SINRY(19) = DSIN(RY(19))
      RZ(19) = (11.599D0 - (0.031D0 * (T - 1997.0D0))) * AM 
      COSRZ(19) = DCOS(RZ(19))
      SINRZ(19) = DSIN(RZ(19))
      SCALE(19) = (-0.93D0 - (0.19D0 * (T - 1997.0D0))) * 1.0D-9
      DELX(19) = 0.9889D0 + (0.0007D0 * (T - 1997.0D0))
      DELY(19) = -1.9074D0 - (0.0001D0 * (T - 1997.0D0))
      DELZ(19) = -0.5030D0 - (0.0019D0 * (T- 1997.0D0))
      
      
*     write(lunit,889) rx(19),cosrx(19),sinrx(19),ry(19),cosry(19),
*    * sinry(19),rz(19),cosrz(19),sinrz(19),scale(19),delx(19),
*    * dely(19),delz(19)
     
*889   format(' code 19 ',/13('  ',f13.10/))     


*** GPS GFILE CODE '21' ITRF00 EPOCH 1997.0 TO NAD83
 
*     RX(21) = 25.915D0*4.84813681D0*1.0D-9
*     COSRX(21) = DCOS(RX(21))
*     SINRX(21) = DSIN(RX(21))
*     RY(21) =  9.426D0*4.84813681D0*1.0D-9
*     COSRY(21) = DCOS(RY(21))
*     SINRY(21) = DSIN(RY(21))
*     RZ(21) = 11.599D0*4.84813681D0*1.0D-9
*     COSRZ(21) = DCOS(RZ(21))
*     SINRZ(21) = DSIN(RZ(21))
*     SCALE(21) = +0.620D0*1.0D-9
*     DELX(21) = 0.9956D0
*     DELY(21) = -1.9013D0
*     DELZ(21) = -0.5215D0

      RX(21) = (25.915D0 + (0.067D0 * (T - 1997.0D0))) * AM 
      COSRX(21) = DCOS(RX(21))
      SINRX(21) = DSIN(RX(21))
      RY(21) =  (9.426D0 - (0.757D0 * (T - 1997.0D0))) * AM 
      COSRY(21) = DCOS(RY(21))
      SINRY(21) = DSIN(RY(21))
      RZ(21) = (11.599D0 - (0.051D0 * (T - 1997.0D0))) * AM 
      COSRZ(21) = DCOS(RZ(21))
      SINRZ(21) = DSIN(RZ(21))
      SCALE(21) = (+0.62D0 - (0.18D0 * (T - 1997.0D0))) * 1.0D-9
      DELX(21) = 0.9956D0 + (0.0007D0 * (T - 1997.0D0))
      DELY(21) = -1.9013D0 - (0.0007D0 * (T - 1997.0D0))
      DELZ(21) = -0.5215D0 + (0.0005D0 * (T- 1997.0D0))
      
      
*     write(lunit,889) rx(21),cosrx(21),sinrx(21),ry(21),cosry(21),
*    * sinry(21),rz(21),cosrz(21),sinrz(21),scale(21),delx(21),
*    * dely(21),delz(21)
     
*889   format(' code 21 ',/13('  ',f13.10/))     
*** GPS GFILE CODE '25' ITRF05 EPOCH 1997.0 TO NAD83
 
*     RX(25) = 25.915D0*4.84813681D0*1.0D-9
*     COSRX(25) = DCOS(RX(25))
*     SINRX(25) = DSIN(RX(25))
*     RY(25) =  9.426D0*4.84813681D0*1.0D-9
*     COSRY(25) = DCOS(RY(25))
*     SINRY(25) = DSIN(RY(25))
*     RZ(25) = 11.599D0*4.84813681D0*1.0D-9
*     COSRZ(25) = DCOS(RZ(25))
*     SINRZ(25) = DSIN(RZ(25))
*     SCALE(25) = +1.020D0*1.0D-9
*     DELX(25) = 0.9957D0
*     DELY(25) = -1.9021D0
*     DELZ(25) = -0.5273D0

      RX(25) = (25.915D0 + (0.067D0 * (T - 1997.0D0))) * AM
      COSRX(25) = DCOS(RX(25))
      SINRX(25) = DSIN(RX(25))
      RY(25) =  (9.426D0 - (0.757D0 * (T - 1997.0D0))) * AM
      COSRY(25) = DCOS(RY(25))
      SINRY(25) = DSIN(RY(25))
      RZ(25) = (11.599D0 - (0.051D0 * (T - 1997.0D0))) * AM
      COSRZ(25) = DCOS(RZ(25))
      SINRZ(25) = DSIN(RZ(25))
      SCALE(25) = (+1.02D0 - (0.10D0 * (T - 1997.0D0))) * 1.0D-9
      DELX(25) = 0.9957D0 + (0.0005D0 * (T - 1997.0D0))
      DELY(25) = -1.9021D0 - (0.0006D0 * (T - 1997.0D0))
      DELZ(25) = -0.5273D0 - (0.0013D0 * (T- 1997.0D0))
      
      
*     write(lunit,889) rx(25),cosrx(25),sinrx(25),ry(25),cosry(25),
*    * sinry(25),rz(25),cosrz(25),sinrz(25),scale(25),delx(25),
*    * dely(25),delz(25)
     
*889   format(' code 25 ',/13('  ',f13.10/))     


*** GPS GFILE CODE '26' IGS05 EPOCH 1997.0 TO NAD83
 
*     RX(26) = 25.904D0*4.84813681D0*1.0D-9
*     COSRX(26) = DCOS(RX(26))
*     SINRX(26) = DSIN(RX(26))
*     RY(26) =  9.419D0*4.84813681D0*1.0D-9
*     COSRY(26) = DCOS(RY(26))
*     SINRY(26) = DSIN(RY(26))
*     RZ(26) = 11.598D0*4.84813681D0*1.0D-9
*     COSRZ(26) = DCOS(RZ(26))
*     SINRZ(26) = DSIN(RZ(26))
*     SCALE(26) = -.8353D0*1.0D-9 
*     DELX(26) = 0.9973D0
*     DELY(26) = -1.9023D0
*     DELZ(26) = -0.5249D0

      RX(26) = (25.904D0 + (0.067D0 * (T - 1997.0D0))) * AM
      COSRX(26) = DCOS(RX(26))
      SINRX(26) = DSIN(RX(26))
      RY(26) =  (9.419D0 - (0.757D0 * (T - 1997.0D0))) * AM
      COSRY(26) = DCOS(RY(26))
      SINRY(26) = DSIN(RY(26))
      RZ(26) = (11.598D0 - (0.051D0 * (T - 1997.0D0))) * AM
      COSRZ(26) = DCOS(RZ(26))
      SINRZ(26) = DSIN(RZ(26))
      SCALE(26) = (-.8353D0 - (0.10D0 * (T - 1997.0D0))) * 1.0D-9
      DELX(26) = 0.9973D0 + (0.0005D0 * (T - 1997.0D0))
      DELY(26) = -1.9023D0 - (0.0006D0 * (T - 1997.0D0))
      DELZ(26) = -0.5249D0 - (0.0013D0 * (T- 1997.0D0))
      
      
*     write(lunit,889) rx(26),cosrx(26),sinrx(26),ry(26),cosry(26),
*    * sinry(26),rz(26),cosrz(26),sinrz(26),scale(26),delx(26),
*    * dely(26),delz(26)
     
*889   format(' code 26 ',/13('  ',f13.10/))     

*** GPS GFILE CODE '27' IGS08 EPOCH 2005 (but t0 = 1997.0) TO NAD83(2011)
 
      RX(27)    = (25.91467D0 + (0.06667D0 * (T - 1997.0D0))) * AM
      COSRX(27) = DCOS(RX(27))
      SINRX(27) = DSIN(RX(27))
      RY(27)    = ( 9.42645D0 - (0.75744D0 * (T - 1997.0D0))) * AM
      COSRY(27) = DCOS(RY(27))
      SINRY(27) = DSIN(RY(27))
      RZ(27)    = (11.59935D0 - (0.05133D0 * (T - 1997.0D0))) * AM
      COSRZ(27) = DCOS(RZ(27))
      SINRZ(27) = DSIN(RZ(27))
      SCALE(27) = ( 1.71504D0 - (0.10201D0 * (T - 1997.0D0))) * 1.0D-9
      DELX(27)  =  0.99343D0  + (0.00079D0 * (T - 1997.0D0))
      DELY(27)  = -1.90331D0  - (0.00060D0 * (T - 1997.0D0))
      DELZ(27)  = -0.52655D0  - (0.00134D0 * (T - 1997.0D0))
      
*** GPS GFILE CODE '28' IGb08 EPOCH 2005 (but t0 = 1997.0) TO NAD83(2011)
 
      RX(28)    = (25.91467D0 + (0.06667D0 * (T - 1997.0D0))) * AM
      COSRX(28) = DCOS(RX(28))
      SINRX(28) = DSIN(RX(28))
      RY(28)    = ( 9.42645D0 - (0.75744D0 * (T - 1997.0D0))) * AM
      COSRY(28) = DCOS(RY(28))
      SINRY(28) = DSIN(RY(28))
      RZ(28)    = (11.59935D0 - (0.05133D0 * (T - 1997.0D0))) * AM
      COSRZ(28) = DCOS(RZ(28))
      SINRZ(28) = DSIN(RZ(28))
      SCALE(28) = ( 1.71504D0 - (0.10201D0 * (T - 1997.0D0))) * 1.0D-9
      DELX(28)  =  0.99343D0  + (0.00079D0 * (T - 1997.0D0))
      DELY(28)  = -1.90331D0  - (0.00060D0 * (T - 1997.0D0))
      DELZ(28)  = -0.52655D0  - (0.00134D0 * (T - 1997.0D0))
      
      RETURN
      END
      
      
      SUBROUTINE INTMAP (B)

*** INITIALIZE SINGULARITY MAP

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      CHARACTER*1 AXCHR, TKCHR, STCHR
      COMMON /EXTRMA/ GLAMAX, GLAMIN, GLOMAX, GLOMIN, HTMAX, HTMIN
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /CONST/  PI, PI2, RAD

      DATA AXCHR /'+'/, TKCHR /'+'/, STCHR /'.'/

*** ASSUME 80 COLUMNS BY 48 ROWS MAP SIZE

      DDY = (GLAMAX - GLAMIN )/47.D0
      YMIN = GLAMIN - DDY
      YMAX = GLAMAX + DDY
      DDX = (GLOMAX - GLOMIN )/79.D0
      XMIN = GLOMIN - DDX
      XMAX = GLOMAX + DDX
      CALL SETRNG (XMIN, XMAX, YMIN, YMAX)
      CALL CLRPLT

*** PUT IN AXES

      CALL PUTXAX (YMIN, AXCHR)
      CALL PUTXAX (YMAX, AXCHR)
      CALL PUTYAX (XMIN, AXCHR)
      CALL PUTYAX (XMAX, AXCHR)

*** PUT IN TICK MARKS
*** INCREMENT FOR TICK MARKS IN BOTH LATITUDE AND LONGITUDE IS 1 DEGREE

      X1 = DNINT( XMIN )
      X2 = DNINT( XMAX )
      YC = DNINT ( (X1 + X2) / 2.D0 )
      DINC = 1.D0/RAD

      DO 120 X = X1, X2, DINC
        CALL PUTYTK ( X, YC, DINC, TKCHR)
  120 CONTINUE

*** PUT IN STATION LOCATIONS

      DO 110 ISN = 1, NSTA
        CALL GETGLA (GLA,ISN,B)
        CALL GETGLO (GLO,ISN,B)
        CALL PUTXY (GLO, GLA, STCHR)
  110 CONTINUE

      RETURN
      END
      SUBROUTINE INVIUN (IUNK,I,J,ICODE)

*** GIVEN AN UNKNOWN INDEX NUMBER, GET STATION/PARM NUMBER

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      N0 = 0
      N1 = NAUX
      N2 = N1 + 3 * NGRT
      N3 = N2 + NZ
      N4 = N3 + 3 * NSTA

*** INDEX IS A STATION INDEX  (ICODE = 0)

      IF (IUNK.GT.N3 .AND. IUNK.LE.N4) THEN
        I = (IUNK-N3-1)/3
        J = IUNK-N3-I*3
        I = I+1
        ICODE = 0

*** INDEX IS AUXILIARY PARAMETER INDEX  (ICODE = 1)

      ELSEIF (IUNK.GT.N0 .AND. IUNK.LE.N1) THEN
        I = IUNK
        J = 0
        ICODE = 1

*** INDEX IS AUXILIARY GPS AND DOPPLER ROTATION PARM. INDEX  (ICODE = 2)

      ELSEIF (IUNK.GT.N1 .AND. IUNK.LE.N2) THEN
        I = ( IUNK - N1 - 1 ) / 3
        J = IUNK - N1 - I * 3
        I = I + 1
        ICODE = 2

*** INDEX IS ROTATION PARAMETER INDEX(ICODE = 3)

      ELSEIF (IUNK.GT.N2 .AND. IUNK.LE.N3) THEN
        I = IUNK-N2
        J = 0
        ICODE = 3

*** ILLEGAL INPUT

      ELSE
        WRITE (LUNIT,1) IUNK,NUNK
    1   FORMAT ('0ILLEGAL VALUE IN INVIUN',2I10)
             write (lunit,*) 'subs2, 8'
        CALL ABORT2
      ENDIF

      RETURN
      END
      SUBROUTINE INVRT3 (A)

*** ROUTINE TO INVERT A 3X3 MATRIX IN PLACE

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER ( N3 = 3, N2 = N3 - 1 )
      PARAMETER ( ONE = 1.D0 )
      DIMENSION A(N3,N3)

      DO 3 I = 1, N3
        AII = ONE/A(I,I)
        A(I,I) = AII
        DO 2 J = 1, N3
          IF (I .NE. J) THEN
            AJI = A(J,I)*AII
            A(J,I) = AJI
            DO 1 K = 1, N3
              IF (I .NE. K) THEN
                A(J,K) = A(J,K) - AJI*A(I,K)
                IF (J .EQ. N3) THEN
                  A(I,K) = -AII*A(I,K)
                ENDIF
              ENDIF
    1       CONTINUE
          ENDIF
    2   CONTINUE
    3 CONTINUE
      ANN = -A(N3,N3)
      DO 4 J = 1, N2
        A(N3,J) = ANN*A(N3,J)
    4 CONTINUE

      RETURN
      END
      LOGICAL FUNCTION INVZ (IZ,ISSN,ILIST)

*** GIVEN THE ROTATION UNKNOWN NUMBER, RETURN STA NUM & LIST NUM

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MAXZ =  8000 )
      COMMON /ROTTB1/ IZNS(MAXZ)
      COMMON /ROTTB2/ LOOK(977), IPTRS(MAXZ), NKEY, MAX, IFREE

      IF (IZ.LE.0 .OR. IZ.GE.IFREE) THEN
        INVZ = .FALSE.
      ELSE
        IZZZ = IZNS(IZ)
        ISSN = IZZZ/1000
        ILIST = IZZZ-ISSN*1000
        INVZ = .TRUE.
      ENDIF

      RETURN
      END
      INTEGER FUNCTION INX (I, N)

*** RETURN INDEX FOR DIAGONAL OF FULL TRIANGULAR MATRIX BY ROWS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      NBAR = N*(N+1)/2
      J = N - I + 1
      JBAR = J*(J+1)/2
      INX = NBAR - JBAR + 1

      RETURN
      END
      INTEGER FUNCTION IPOINT (V)

*** LOCATE NEW POSITON IN THE MAX V ARRAY

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRD = 20 )
      COMMON /RESTAT/ VSD20(MXRD), VMAX, VMIN, VSDMAX, VSDMIN, VSUM,
     &                VSDSUM, VSD21, VSD22, VSD23, VSD24, VSD25, VSD26,
     &                VSD27, VSD28, VSD29, VSD210, VSD211, VSD212,
     &                VABS1, VABS2, VABS3, VABS4, VABS5, VABS6, VABS7,
     &                VABS8, VABS9, VABS10, VABS11, VABS12, I20(MXRD),
     &                N0, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10,
     &                N11, N12, NI20
      COMMON /REDUN/  RN1, RN2, RN3, RN4, RN5, RN6, RN7, RN8, RN9,
     &                RN10, RN11, RN12

*** FIND NEW POSITION (CHECK MAX VALUES FIRST)

      DO 1 I = 1,MXRD
      IPOINT = I
      IF (V .GT. VSD20(I)) RETURN
    1 CONTINUE

*** FALL THRU LOOP--NOT A MAXIMAL RESIDUAL

      IPOINT = MXRD + 1

      RETURN
      END
      
      
**v 4.28

      INTEGER FUNCTION IPOINTV (V)

*** LOCATE NEW POSITON IN THE MAX V ARRAY

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRDV = 20 )
      COMMON /RESTA3/ V20(MXRDV), I20V(MXRDV), NI20V  

*** FIND NEW POSITION (CHECK MAX VALUES FIRST)

      DO 1 I = 1,MXRDV
      IPOINTV = I
      IF (V .GT. V20(I)) RETURN
    1 CONTINUE

*** FALL THRU LOOP--NOT A MAXIMAL RESIDUAL

      IPOINTV = MXRDV + 1

      RETURN
      END
      
***************
** v 4.29f
      INTEGER FUNCTION IPOINTU (V)

*** LOCATE NEW POSITON IN THE MAX V ARRAY FOR DU COMPONENT

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRDU = 20 )
      COMMON /RESTA7/ VDU20(MXRDU), I20DU(MXRDU), NI20DU  

*** FIND NEW POSITION (CHECK MAX VALUES FIRST)

      DO 1 I = 1,MXRDU
      IPOINTU = I
      IF (V .GT. VDU20(I)) RETURN
    1 CONTINUE

*** FALL THRU LOOP--NOT A MAXIMAL RESIDUAL

      IPOINTU = MXRDU + 1

      RETURN
      END

*******************************

** v 4.29i

      INTEGER FUNCTION IPOINTL (V)

*** LOCATE NEW POSITON IN THE MAX V ARRAY FOR DL COMPONENT

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRDL = 20 )
      COMMON /RESTA8/ VDL20(MXRDL), I20DL(MXRDL), NI20DL  

*** FIND NEW POSITION (CHECK MAX VALUES FIRST)

      DO 1 I = 1,MXRDL
      IPOINTL = I
      IF (V .GT. VDL20(I)) RETURN
    1 CONTINUE

*** FALL THRU LOOP--NOT A MAXIMAL RESIDUAL

      IPOINTL = MXRDL + 1

      RETURN
      END

*******************************8
      INTEGER FUNCTION ITCODE (TC)

*** FIND HOUR SHIFT TO MOVE TO GMT (ZULU) FOR A NAVY TIME CODE

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*1 TC
      COMMON/UNITS/ LUNIT

*** THE 'IF' CLAUSES ARE USED IN ADDITION TO THE 'ELSEIF' CLAUSES
*** TO PREVENT THE COMPILER FROM COMPLAINING ABOUT THE MAXIMUM
*** ALLOWABLE NUMBER OF NESTED LEVELS.
*** IF THIS ERROR MESSAGE OCCURS, EITHER HAVE THE COMPILER INCREASE
*** ITS LIMIT ON THE NUMBER OF NESTED LEVELS, OR CONVERT MORE OF THE
*** 'ELSEIF' CLAUSES INTO 'IF' CLAUSES WITH 'RETURN' STATEMENTS.

      IF (TC.EQ.'Z') THEN
        ITCODE = 0
        RETURN
      ENDIF
      IF (TC.EQ.'N') THEN
        ITCODE = +1
        RETURN
      ENDIF
      IF (TC.EQ.'O') THEN
        ITCODE = +2
        RETURN
      ENDIF
      IF (TC.EQ.'P') THEN
        ITCODE = +3
      ELSEIF (TC.EQ.'Q') THEN
        ITCODE = +4
      ELSEIF (TC.EQ.'R') THEN
        ITCODE = +5
      ELSEIF (TC.EQ.'S') THEN
        ITCODE = +6
      ELSEIF (TC.EQ.'T') THEN
        ITCODE = +7
      ELSEIF (TC.EQ.'U') THEN
        ITCODE = +8
      ELSEIF (TC.EQ.'V') THEN
        ITCODE = +9
      ELSEIF (TC.EQ.'W') THEN
        ITCODE = +10
      ELSEIF (TC.EQ.'X') THEN
        ITCODE = +11
      ELSEIF (TC.EQ.'Y') THEN
        ITCODE = +12
      ELSEIF (TC.EQ.'A') THEN
        ITCODE = -1
      ELSEIF (TC.EQ.'B') THEN
        ITCODE = -2
      ELSEIF (TC.EQ.'C') THEN
        ITCODE = -3
      ELSEIF (TC.EQ.'D') THEN
        ITCODE = -4
      ELSEIF (TC.EQ.'E') THEN
        ITCODE = -5
      ELSEIF (TC.EQ.'F') THEN
        ITCODE = -6
      ELSEIF (TC.EQ.'G') THEN
        ITCODE = -7
      ELSEIF (TC.EQ.'H') THEN
        ITCODE = -8
      ELSEIF (TC.EQ.'I') THEN
        ITCODE = -9
      ELSEIF (TC.EQ.'K') THEN
        ITCODE = -10
      ELSEIF (TC.EQ.'L') THEN
        ITCODE = -11
      ELSEIF (TC.EQ.'M') THEN
        ITCODE = -12
      ELSE
        WRITE (LUNIT,1) TC
    1   FORMAT ('0ILLEGAL TIME CODE--',A1)
             write (lunit,*) 'subs2, 9'
        CALL ABORT2
      ENDIF

      RETURN
      END
      INTEGER FUNCTION IUNAUX (IAUX)

*** DETERMINE UNKNOWN NUMBER OF AUXILIARY PARAMETER

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

*** PARAMETERS STORED AHEAD OF STATIONS

      IF (IAUX.LE.0 .OR. IAUX.GT.NAUX) THEN
        WRITE (LUNIT,1) IAUX,NAUX
    1   FORMAT ('0ILLEGAL VALUES IN IUNAUX',2I8)
             write (lunit,*) 'subs2, 10'
        CALL ABORT2
      ELSE
        IUNAUX = IAUX
      ENDIF

      RETURN
      END
      INTEGER FUNCTION IUNGRT (IGRT,I)

*** DETERMINE UNKNOWN NUMBER OF AUXILIARY GPS AND DOPPLER ROTATION PARM.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

*** PARAMETERS STORED AHEAD OF STATIONS

      IF ( IGRT.LE.0 .OR. IGRT.GT.NGRT .OR.
     &        I.LT.1 .OR.    I.GT.3 ) THEN
        WRITE (LUNIT,1) IGRT,NGRT,I
    1   FORMAT ('0ILLEGAL VALUES IN IUNGRT',3I8)
             write (lunit,*) 'subs2, 11'
        CALL ABORT2
      ELSE
        IUNGRT = NAUX + ( IGRT - 1 ) * 3 + I
      ENDIF

      RETURN
      END
      INTEGER FUNCTION IUNROT (IZ)

*** DETERMINE UNKNOWN NUMBER FOR ROTATION PARAMETER

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

*** PARAMETER STORED AHEAD OF STATIONS

      IF (IZ.LE.0 .OR. IZ.GT.NZ) THEN
        WRITE (LUNIT,1) IZ,NZ
 1      FORMAT ('0ILLEGAL VALUES IN IUNROT: IZ=',I5,'   NZ=',I5)
             write (lunit,*) 'subs2, 12'
        CALL ABORT2
      ELSE
        IUNROT = NAUX + (NGRT * 3) +IZ
      ENDIF

      RETURN
      END
      INTEGER FUNCTION IUNSHF (ISTA,I)

*** DETERMINE SHIFT INDEX OF A STATION

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

*** SHIFTS STORED THREEWISE--U,V,W IN LOCAL GEOD. HOR.

      IF ( ISTA.LE.0 .OR. ISTA.GT.NSTA .OR.
     &        I.LE.0 .OR.    I.GT.3) THEN
        WRITE (LUNIT,1) ISTA,I,NSTA
    1   FORMAT ('0ILLEGAL VALUES IN IUNSHF',3I8)
             write (lunit,*) 'subs2, 13'
        CALL ABORT2
      ELSE
        IUNSHF = (ISTA-1)*3+I
      ENDIF

      RETURN
      END
      INTEGER FUNCTION IUNSTA (ISTA,I)

*** DETERMINE UNKNOWN NUMBER OF STATION COORDINATE

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

*** STATIONS STORED THREEWISE--U,V,W IN LOCAL GEOD. HOR.
*** STATIONS STORED AFTER AUXILIARY PARAMETERS

      IF ( ISTA.LE.0 .OR. ISTA.GT.NSTA .OR.
     &        I.LE.0 .OR.    I.GT.3) THEN
        WRITE (LUNIT,1) ISTA,I,NSTA
    1   FORMAT ('0ILLEGAL VALUES IN IUNSTA',3I8)
             write (lunit,*) 'subs2, 14'
        CALL ABORT2
      ELSE
        IUNSTA = NAUX + (NGRT * 3) +NZ+(ISTA-1)*3+I
      ENDIF

      RETURN
      END
      SUBROUTINE JOBSTT

*** LIST THE JOB STATISTICS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      LOGICAL L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

*** HEADING

      CALL HEAD
      WRITE (LUNIT,1)
1     FORMAT ('0*** JOB STATISTICS ***', /)

*** PRINT STATS

** v 4.29k
*     WRITE (LUNIT,2) NSTA,N84,N85,N86,NDIR,NANG,NGPS,NDOP,NZD,NDS,
      WRITE (LUNIT,2) NSTA,N85,N86,NDIR,NANG,NGPS,NDOP,NZD,NDS,
     &            NAZ, NCON, NQQ, NREJ
2     FORMAT ('0A.) BLUE-BOOK STATISTICS', /
     &        9X, 'NO. *80* CONTROL RECORDS', I10, /
** v 4.29k
*    &        9X, 'NO. *84* GEOID HT. RECORDS', I8, /
     &        9X, 'NO. *85* DEFLECTION RECORDS', I7, /
     &        9X, 'NO. *86* ELEVATION RECORDS', I8, /
     &        9X, 'NO. DIRECTIONS', I20, /
     &        9X, 'NO. ANGLES', I24, /
     &        9X, 'NO. GPS VECTORS', I19, /
     &        9X, 'NO. DOPPLER OBS.', I18, /
     &        9X, 'NO. ZENITH DISTANCES', I14, /
     &        9X, 'NO. DISTANCES', I21, /
     &        9X, 'NO. AZIMUTHS', I22, /
     &        1X, 'B.) NO. CONSTRAINTS', I23, /
     &        1X, 'C.) NO. ACCURACIES', I24, /
     &        1X, 'D.) NO. REJECTED OBS.', I21)

      IF (L2HLF) THEN
        WRITE (LUNIT,3) N2HLF
    3   FORMAT (' E.) NO. DUAL HEIGHT STA.', I18)
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION LEAP (IYEAR)

*** IS THE YEAR A LEAP YEAR?

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      IF (MOD(IYEAR,400).EQ.0) THEN
        LEAP = .TRUE.
      ELSEIF (MOD(IYEAR,100).EQ.0) THEN
        LEAP = .FALSE.
      ELSEIF (MOD(IYEAR,4).EQ.0) THEN
        LEAP = .TRUE.
      ELSE
        LEAP = .FALSE.
      ENDIF

      RETURN
      END
      SUBROUTINE LINE (I)

*** MAINTAIN A LINE COUNT AND PERFORM PAGE INTERRUPTS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE

      ILINE = ILINE+I
      IF (ILINE.GT.MAXLIN) THEN
        CALL HEAD
        ILINE = ILINE+I
      ENDIF

      RETURN
      END
      SUBROUTINE LINE2 (I)

*** MAINTAIN A LINE COUNT AND PERFORM PAGE INTERRUPTS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE

      ILINE = ILINE+I
      IF (ILINE.GT.MAXLIN) THEN
        CALL HEAD2
        ILINE = ILINE+I
      ENDIF

      RETURN
      END
      SUBROUTINE LINE3 (I)

*** MAINTAIN A LINE COUNT AND PERFORM PAGE INTERRUPTS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE

      ILINE = ILINE+I
      IF (ILINE.GT.MAXLIN) THEN
        CALL HEAD3
        ILINE = ILINE+I
      ENDIF

      RETURN
      END
      SUBROUTINE LINE4 (I)

*** MAINTAIN A LINE COUNT AND PERFORM PAGE INTERRUPTS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE

      ILINE = ILINE+I
      IF (ILINE.GT.MAXLIN) THEN
        CALL HEAD4
        ILINE = ILINE+I
      ENDIF

      RETURN
      END
      SUBROUTINE LINE5 (I)

*** MAINTAIN A LINE COUNT AND PERFORM PAGE INTERRUPTS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE

      ILINE = ILINE+I
      IF (ILINE.GT.MAXLIN) THEN
        CALL HEAD5
        ILINE = ILINE+I
      ENDIF

      RETURN
      END
      SUBROUTINE LINE6 (I)

*** MAINTAIN A LINE COUNT AND PREFORM PAGE INTERUPTS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE

      ILINE = ILINE+I
      IF (ILINE.GT.MAXLIN) THEN
        CALL HEAD6
        ILINE = ILINE+I
      ENDIF

      RETURN
      END
      SUBROUTINE LINE7 (I)

*** MAINTAIN A LINE COUNT AND PREFORM PAGE INTERRUPTS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE

      ILINE = ILINE+I
      IF (ILINE.GT.MAXLIN) THEN
        CALL HEAD7
        ILINE = ILINE+I
      ENDIF

      RETURN
      END
      SUBROUTINE LINE8 (I)

*** MAINTAIN A LINE COUNT AND PERFORM PAGE INTERRUPTS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE

      ILINE = ILINE+I
      IF (ILINE.GT.MAXLIN) THEN
        CALL HEAD8
        ILINE = ILINE+I
      ENDIF

      RETURN
      END
      SUBROUTINE LINE9 (I)

*** MAINTAIN A LINE COUNT AND PERFORM PAGE INTERRUPTS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE

      ILINE = ILINE+I
      IF (ILINE.GT.MAXLIN) THEN
        CALL HEAD9
        ILINE = ILINE+I
      ENDIF

      RETURN
      END
      SUBROUTINE LOADCR (G, GCARD, NR, NC)

*** LOAD CORRELATIONS FROM A RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 GCARD
      DIMENSION G(NR,NC)
      DIMENSION IROW(5), ICOL(5), COR(5)
      COMMON/UNITS/ LUNIT

      READ (GCARD,1) (IROW(I), ICOL(I), COR(I), I = 1, 5)
c     write(*    ,1) (IROW(I), ICOL(I), COR(I), I = 1, 5)
*   1 FORMAT (1X, 5(2I3,F9.7) )
    1 FORMAT (BZ, 1X, 5(2I3,F9.7) )

      DO 10 I = 1, 5
        IR = IROW(I)
        IC = ICOL(I)
        C = COR(I)
        IF (IR .NE. 0  .AND.  IC .NE. 0) THEN
          IF (DABS(C) .GT. 1.D0) THEN
            WRITE (LUNIT,2) IR, IC, C
    2       FORMAT ('0ILLEGAL CORRELATION--',2I8,D20.10)
             write (lunit,*) 'subs2, 15'
            CALL ABORT2
          ELSEIF (IR .LT. 1  .OR.  IR .GT. NR  .OR.
     &            IC .LT. 1  .OR.  IC .GT. NR) THEN
            WRITE (LUNIT,3) IR, IC
    3       FORMAT ('0BAD CORRELATION INDEX--',2I8)
             write (lunit,*) 'subs2, 16'
            CALL ABORT2
          ELSEIF (IR .EQ. IC) THEN
            CALL LINE (1)
            WRITE (LUNIT,4) IR, IC, C
    4       FORMAT (' *** WARNING -- DIAGONAL CORRELATION ',2I8,D20.10)
          ELSE
            G(IR,IC) = C
            G(IC,IR) = C
          ENDIF
        ENDIF
   10 CONTINUE

      RETURN
      END
      SUBROUTINE LOADCV (G, GCARD, NR, NC)

*** LOAD COVARIANCES FROM A RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 GCARD
      DIMENSION G(NR,NC)
      DIMENSION IROW(4), ICOL(4), COV(4)
      COMMON/UNITS/ LUNIT

      READ (GCARD,1) (IROW(I), ICOL(I), COV(I), I = 1, 4)
*   1 FORMAT (1X, 4(2I3,F12.8))
    1 FORMAT (BZ, 1X, 4(2I3,F12.8))

      DO 10 I = 1, 4
        IR = IROW(I)
        IC = ICOL(I)
        C = COV(I)

c  This next line can generate bogus indices into G
c        VAR = G(IR,IR)*G(IC,IC)
c  Move it inside the do-loop Mike Potterfield 10/1/08

        IF(IR .NE. 0 .AND. IC .NE. 0) THEN

c put the assignment statement here
          VAR = G(IR,IR)*G(IC,IC)
c end bug fix 10/1/08

          IF (DABS(C) .GT. VAR) THEN
            WRITE (LUNIT,2) IR, IC, C
    2       FORMAT ('0ILLEGAL COVARIANCE--',2I8,D20.10) 
             write (lunit,*) 'subs2, 17'
            CALL ABORT2
          ELSEIF (IR .LT. 1  .OR.  IR .GT. NR  .OR.
     &            IC .LT. 1  .OR.  IC .GT. NR) THEN
            WRITE (LUNIT,3) IR, IC
    3       FORMAT ('0BAD COVARIANCE INDEX--',2I8)
             write (lunit,*) 'subs2, 18'
            CALL ABORT2
          ELSEIF (IR .EQ. IC) THEN
            CALL LINE (1)
            WRITE (LUNIT,4) IR, IC, C
    4       FORMAT (' *** WARNING -- DIAGONAL COVARIANCE ',2I8,D20.10)
          ELSE
            G(IR,IC) = C
            G(IC,IR) = C
          ENDIF
        ENDIF
   10 CONTINUE
      RETURN
      END

      LOGICAL FUNCTION LOCSSN (ISN,ISSN)

*** SCAN TABLE FOR THE ISN-TH SSN

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      COMMON /SSNTBL/ NSSN, ISSISN(MXSSN)
      SAVE IPOINT
      DATA IPOINT /1/

      IF (ISN.LE.0 .OR. ISN.GT.NSSN) THEN
        LOCSSN = .FALSE.
        RETURN
      ELSE
        DO 1 I = IPOINT, MXSSN
          IF (ISSISN(I).EQ.ISN) THEN
            IPOINT = I
            ISSN = I
            LOCSSN = .TRUE.
            RETURN
          ENDIF
    1   CONTINUE

*** WRAP AROUND

        I1 = IPOINT-1
        DO 2 I = 1,I1
          IF (ISSISN(I).EQ.ISN) THEN
            IPOINT = I
            ISSN = I
            LOCSSN = .TRUE.
            RETURN
          ENDIF
    2   CONTINUE

*** FELL THRU TABLE--ISN NOT LOCATED

        LOCSSN = .FALSE.
        RETURN
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION LOKZ (IZN,ILOOK,IPTR)

*** TABLE LOOKUP & POINTER RETURN

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MAXZ =  8000 )
      COMMON /ROTTB1/ IZNS(MAXZ)
      COMMON /ROTTB2/ LOOK(977), IPTRS(MAXZ), NKEY, MAX, IFREE

*** CHARACTER TYPE HASH FUNCTIONS

      KEY = MOD(IZN,NKEY)+1

*** CHARACTER TYPE HASH FUNCTION

      ILOOK = LOOK(KEY)
      IF (ILOOK.EQ.0) THEN

*** IMMEDIATE FAILURE TO FIND ENTRY
        ILOOK = KEY
        IPTR = 0
        LOKZ = .FALSE.
      ELSE
        IPTR = ILOOK
        ILOOK = 0
 1      IF (IZN.EQ.IZNS(IPTR)) THEN

*** SUCCESSFUL LOCATION OF ENTRY AT IPTR

          LOKZ = .TRUE.
        ELSE
          IF (IPTRS(IPTR).EQ.0) THEN

*** FAILURE TO FIND ENTRY

            LOKZ = .FALSE.
          ELSE

*** CHAIN TO NEXT LOCATION

            IPTR = IPTRS(IPTR)
            GO TO 1
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END
      SUBROUTINE MAPSNG (ISING, ISN, B, CHR)

*** PUT SINGULARITY LOCATION ON MAP

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      CHARACTER*26 ALPHA
      CHARACTER*1 CHR

      DATA ALPHA /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      CALL GETGLA (GLA, ISN, B)
      CALL GETGLO (GLO,ISN,B)
      IC = MOD(ISING, 26)
      IF (IC .EQ. 0) IC = 26
      CHR = ALPHA(IC:IC)
      CALL PUTXY (GLO, GLA, CHR )

      RETURN
      END
      SUBROUTINE NBLANK (A,INUM,IBLK)

*** RETURN PRECISION OF NUMBER AND ZERO FILL

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*(*) A

      LENG = LEN(A)
      L1 = LENG-INUM+1
      IBLK = 0
      DO 1 I = L1,LENG
        IF (A(I:I).EQ.' ') THEN
          IBLK = IBLK+1
          A(I:I) = '0'
        ENDIF
    1 CONTINUE

      RETURN
      END
      SUBROUTINE NEWFOT

*** INITIALIZE THE TABLE TO ZERO LENGTH

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXFOT = 3500 )
      COMMON /FOURTH/ NFOT, NFOTS(MXFOT)

      DO 1 I = 1,MXFOT
        NFOTS(I) = 0
    1 CONTINUE

      RETURN
      END
      

      SUBROUTINE NEWGRP (NVEC, IUNIT, G, NR, NC, FULL, ICM, NICM, KINDS,
** v 4.28**********************
*    &                ISNS, JSNS, LOBS, IOBS, IAUX, B, NVECS, IGRT, I83)
     &                ISNS, JSNS, LOBS, IOBS, IAUX, B, NVECS, IGRT, I83,
     &                ivc, IUO3)
**********************************

*** LOAD AN OBSERVATION GROUP INTO WORK SPACE

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      
**v 4.28
      PARAMETER ( MAXC = 99999)
      
      LOGICAL FULL, LBV, GETSSN
      LOGICAL LEB, LLB, LEG, LLG, LED, LLD
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      CHARACTER*80 GCARD
      CHARACTER*1 ID, CODE
      DIMENSION G(NR,NC)
      
**v 4.28
      CHARACTER*5 SESSID,SESS
      COMMON /SESTAB/ SESS(MAXC)
      DIMENSION ICM(*)
      DIMENSION KINDS(*), ISNS(*), JSNS(*), LOBS(*)
      DIMENSION IC(LENC), B(*)
      COMMON /PECHO/  LEB, LLB, LEG, LLG, LED, LLD
      COMMON /ECHO/   VSD, LBV
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON/UNITS/ LUNIT

*** CLEAR WORK AREA

      DO 1 I = 1, NR
        DO 1 J = 1, NC
          G(I,J) = 0.D0
    1 CONTINUE

*** LOOP OVER NVEC MEMBER RECORDS  ('C' AND/OR  'F')

      DO 10 I = 1, NVEC
        J1 = (I-1)*3 + 1
        J2 = J1 + 1
        J3 = J2 + 1
  100   READ (IUNIT,2,END=666) GCARD
    2   FORMAT (A80)
**v6.2.1
        if (GCARD(1:1) == ' ') goto 100
**fo far for v6.2.1
        READ (GCARD,22) ID
   22   FORMAT (A1)
        IF (ID .EQ. 'C') THEN
**v 4.28

*         READ (GCARD,3) ISSN,JSSN,DX,SDX,DY,SDY,DZ,SDZ,CODE
*   3     FORMAT (1X,2I4,3(F11.4,F5.4),A1)
*       ELSEIF (ID.EQ.'F') THEN
*         READ (GCARD,4) ISSN,JSSN,DX,SDX,DY,SDY,DZ,SDZ,CODE
*   4     FORMAT (1X,2I4,3(F13.4,F5.4),A1)

          READ (GCARD,3) ISSN,JSSN,DX,SDX,DY,SDY,DZ,SDZ,CODE,SESSID
    3     FORMAT (1X,2I4,3(F11.4,F5.4),A1,1X,A5)
        ELSEIF (ID.EQ.'F') THEN
          READ (GCARD,4) ISSN,JSSN,DX,SDX,DY,SDY,DZ,SDZ,CODE,SESSID
    4     FORMAT (1X,2I4,3(F13.4,F5.4),A1,1X,A5)
*****************************************************
        ELSEIF (ID.EQ.'G'.OR.ID.EQ.'H'.OR.ID.EQ.'I') THEN
*****   ELSEIF (ID.EQ.'G'.OR.ID.EQ.'H') THEN
          GO TO 100
        ELSE
          GO TO 666
        ENDIF

c Mike Potterfield 2/13/06 now write the Reject Code for each vector
      WRITE(IUO3) CODE

          IF ( .NOT. GETSSN(ISSN,ISN) ) THEN
            CALL LINE (3)
            WRITE (LUNIT,7) GCARD
    7       FORMAT ('0NO *80* RECORD FOR--',A80/)
             write (lunit,*) 'subs2, 19'
            CALL ABORT2
          ELSEIF ( .NOT. GETSSN(JSSN,JSN) ) THEN
            CALL LINE (3)
            WRITE (LUNIT,7) GCARD
             write (lunit,*) 'subs2, 20'
            CALL ABORT2
          ENDIF

*** CONVERT GPS VECTOR TO DATUM

          IF (I83 .NE. 0) CALL TO83 (DX, DY, DZ, ISN, I83, B)

          ISNS(J1) = ISN
          ISNS(J2) = ISN
          ISNS(J3) = ISN

          JSNS(J1) = JSN
          JSNS(J2) = JSN
          JSNS(J3) = JSN

          KINDS(J1) = 4
          KINDS(J2) = 5
          KINDS(J3) = 6

          G(J1,NC) = DX
          G(J2,NC) = DY
          G(J3,NC) = DZ

          IF ( CODE .NE. 'R'  .AND.  CODE .NE. 'O'  .AND.
     &         CODE .NE. 'F')  THEN
            G(J1,J1) = SDX
            G(J2,J2) = SDY
            G(J3,J3) = SDZ

            IOBS = IOBS + 1
            LOBS(J1) = IOBS
            IOBS = IOBS + 1
            LOBS(J2) = IOBS
            IOBS = IOBS + 1
            LOBS(J3) = IOBS
            NOBS = NOBS + 3
            NGPS = NGPS + 1

*** IF REJECTED VECTOR, THEN USE WEIGHTS OF 100 METERS

          ELSE
            G(J1,J1) = 100.D0
            G(J2,J2) = 100.D0
            G(J3,J3) = 100.D0
            LOBS(J1) = IOBS
            LOBS(J2) = IOBS
            LOBS(J3) = IOBS
            NREJ = NREJ + 3
            NGPSR = NGPSR + 1
          ENDIF
**v 4.28**********************
          SESS(IVC) = SESSID
	  IVC = IVC + 1
        IF(IVC.GT.MAXC) THEN
          WRITE(LUNIT,555) MAXC
 555      FORMAT(/ ' NUMBER OF SESSIONS EXCEED ',I5)
          write (lunit,*) 'subs2, 21'
          CALL ABORT2
         ENDIF

******************************

*** ECHO GPS OBSERVATION

          IF ( .NOT. LLG) THEN
            IF (LGF) THEN
              IF ( CODE .NE. 'R'  .AND.  CODE .NE. 'O'  .AND.
     &             CODE .NE. 'F')  THEN
                CALL LINE (1)
                WRITE (LUNIT,5) IOBS, GCARD
 5              FORMAT (I7,3X,A80)
              ELSE
                CALL LINE (1)
                WRITE (LUNIT,11) GCARD
 11             FORMAT (10X,A80)
              ENDIF
            ENDIF
          ENDIF

          CALL FORMIC (KINDS(J1), ISNS(J1), JSNS(J1), IAUX, IC, L, IGRT)
          CALL PUTICM (IC, L, ICM, NICM, NVECS)

*** THE FOLLOWING FOUR LINES ARE CURRENTLY REDUNDANT WITH THE PRIOR TWO
*** OF COURSE, IF CHANGES THE THE GPS OBSERVATION MODEL ARE MADE,
*** THEN IT MAY BE NECESSARY TO REMOVE COMMENTS ON THESE FOUR LINES

***       CALL FORMIC (KINDS(J2), ISNS(J2), JSNS(J2), IAUX, IC, L, IGRT)
***       CALL PUTICM (IC, L, ICM, NICM, NVECS)
***       CALL FORMIC (KINDS(J3), ISNS(J3), JSNS(J3), IAUX, IC, L, IGRT)
***       CALL PUTICM (IC, L, ICM, NICM, NVECS)

   10 CONTINUE

*** FELL THRU LOOP--NORMAL END OF PROCESSING

      FULL = .TRUE.
      RETURN

*** BAD GFILE STRUCTURE

  666 WRITE (LUNIT,667) GCARD
  667 FORMAT ('0BAD GFILE STRUCTURE--',A80)
          write (lunit,*) 'subs2, 22'
      CALL ABORT2
      RETURN
      END
      SUBROUTINE NEWGRT

*** INITIALIZE THE TABLE TO ZERO LENGTH

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PRMTB2/ NPARMS, ICODES(5), IOLDS(5), INEWS(5)

      NPARMS = 0

      RETURN
      END
      SUBROUTINE NEWI

*** INITIALIZE THE TABLE TO ZERO

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, MAXZ =  8000, MXPRM = 40 )
      PARAMETER ( MXALL = MXPRM + 15 + MAXZ + MXSSN*3)
      COMMON /OBSUM/  NSUM(MXSSN,13), ICMPL(MXALL)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT

      DO 1 I = 1, NSTA
        DO 2 J = 1, 13
          NSUM(I,J) = 0
    2   CONTINUE
    1 CONTINUE

      RETURN
      END
      SUBROUTINE NEWICM (NICMS)

*** INITIALIZE THE TABLE TO ZERO LENGTH

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      NICMS = 0

      RETURN
      END
      SUBROUTINE NEWIVF

*** INITIALIZE THE TABLE TO ZERO LENGTH

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXVF = 40 )
      LOGICAL LFIXS
      COMMON /IVFTB1/ NIVFS, ICODES(MXVF), IOLDS(MXVF), INEWS(MXVF),
     &                LFIXS(MXVF)
      COMMON /NUMVFS/ NVFTOT, NVFREE

      NIVFS = 0
      NVFTOT = 0
      NVFREE = 0

      RETURN
      END
      SUBROUTINE NEWPRM

*** INITIALIZE THE TABLE TO ZERO LENGTH

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXPRM = 40 )
      COMMON /PRMTB1/ NPARMS, ICODES(MXPRM), IOLDS(MXPRM), INEWS(MXPRM)

      NPARMS = 0

      RETURN
      END
      SUBROUTINE NEWSSN

*** INITIALIZE THE TABLE TO ZERO LENGTH

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      COMMON /SSNTBL/ NSSN, ISSISN(MXSSN)
      COMMON /EXTRMA/ GLAMAX, GLAMIN, GLOMAX, GLOMIN, HTMAX, HTMIN

      DO 1 I = 1, MXSSN
        ISSISN(I) = 0
    1 CONTINUE
      NSSN = 0

*** INITIALIZE EXTREMA

      GLAMAX = -1.D22
      GLAMIN =  1.D22
      GLOMAX = -1.D22
      GLOMIN =  1.D22
      HTMAX = -1.D22
      HTMIN =  1.D22

      RETURN
      END
      SUBROUTINE NEWZ

*** INITIALIZE THE ROT. PARM. TABLE

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MAXZ =  8000 )
      COMMON /ROTTB2/ LOOK(977), IPTRS(MAXZ), NKEY, MAX, IFREE

      NKEY = 977
      MAX = MAXZ
      IFREE = 1

      DO 1 I = 1,NKEY
        LOOK(I) = 0
    1 CONTINUE

      RETURN
      END
      SUBROUTINE NEW82

*** INITIALIZE THE *82* RECORD SSN TABLE

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      COMMON /CHILD/ N82,ISS82(MXSSN)

      DO 1 I = 1,MXSSN
        ISS82(I) = 0
    1 CONTINUE
      N82 = 0

      RETURN
      END
      SUBROUTINE NORMAL (IUO, ITER, B, A, G, NX, LAWORK, GOOGE)

*** FORM AND SOLVE NORMAL EQUATIONS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( NVECS = 700, LENC = 10, MXVF = 40, NSING = 200 )
      PARAMETER ( MXSSN = 9999 )
      DIMENSION ICM((NVECS+1)*3+1+3)
      DIMENSION KINDS(NVECS*3), ISNS(NVECS*3), JSNS(NVECS*3)
      DIMENSION LOBS(NVECS*3)
      LOGICAL ADDOBS, OPENN, GETA, SOLVE, ADDCOR
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      DIMENSION IC(LENC), C(LENC)
      DIMENSION C3(3,LENC), EL3(3), COVECF(3,3)
      DIMENSION ISING(NSING), GSING(NSING)
      DIMENSION B(*), A(*), G(*), NX(*), GOOGE(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON /VFTB2/  VFS(MXVF), VTV(MXVF), VFRN(MXVF), VFSS(MXVF)
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT
      SAVE NCON0

*** SINGULARITY TOLERANCE

*     STOL = 1.D-6
      STOL = 1.D-9

*** REWIND OBSERVATION EQUATIONS AND INITIALIZE NORMALS

      REWIND IUO
      IRANK = NUNK
      NONCR = ( 3 - IDIM )*NSTA
      IF (IDIM .EQ. 2  .AND.  N2HLF .GT. 0) NONCR = NONCR - N2HLF
      DENOM = ( ( DBLE( IRANK - NONCR ) )**2 )/2.D0
      INEED = NX(IRANK+1) + IRANK
      AVBW = 100.D0*( INEED - NONCR )/DENOM
      IF ( .NOT. OPENN(A, NX, IRANK, LAWORK) ) THEN
        WRITE (LUNIT,20) LAWORK, INEED, IRANK
   20   FORMAT (/, I9, ' DOUBLE PRECISION WORDS LESS THAN THE ', I9,
     &             ' NEEDED FOR A RANK OF', I9,
     &          /, ' INSUFFICIENT STORAGE -- FATAL!', /)
        IF ( .NOT. L2HLF) THEN
          WRITE (LUNIT,30) IDIM, NSTA, AVBW
   30     FORMAT (/, ' THE AVERAGE BAND WIDTH FOR THE ', I1,
     &               ' DIMENSIONAL ADJUSTMENT OF', I6, ' STATIONS',
     &               ' WOULD BE ', F5.1, '%.')
        ELSE
          WRITE (LUNIT,31) NSTA, AVBW
   31     FORMAT (/, ' THE AVERAGE BAND WIDTH FOR THE 2.5',
     &               ' DIMENSIONAL ADJUSTMENT OF', I6, ' STATIONS',
     &               ' WOULD BE ', F5.1, '%.')
        ENDIF
          write (lunit,*) 'subs2, 23'
        CALL ABORT2
      ELSEIF (ITER .EQ. 0) THEN
        CALL LINE (2)
        IF ( .NOT. L2HLF) THEN
          WRITE (LUNIT,25) IDIM, NSTA, IRANK, AVBW, INEED
   25     FORMAT (/, ' THE AVERAGE BAND WIDTH FOR THE ', I1,
     &               ' DIM. ADJUSTMENT OF', I6, ' STATIONS',
     &               ' AND RANK', I7, ' IS ', F5.1, '%.',
     &               '   D.P. WORDS NEEDED=', I9)
        ELSE
          WRITE (LUNIT,26) NSTA, IRANK, AVBW, INEED
   26     FORMAT (/, ' THE AVERAGE BAND WIDTH FOR THE 2.5',
     &               ' DIM. ADJUSTMENT OF', I6, ' STATIONS',
     &               ' AND RANK', I7, ' IS ', F5.1, '%.',
     &               '   D.P. WORDS NEEDED=', I9)
        ENDIF
      ENDIF

*** LOAD OBSERVATION EQUATIONS AND FORM NORMALS

  100 READ (IUO,END=777) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD,
     &                   IOBS, IVF, IAUX, IGRT
      IF (KIND .LE. 999) THEN
        IF (KIND .LT. 18  .OR.  KIND .GT. 20) THEN
          EL = -CMO
          VAR = SD*SD
          IF (NVFTOT .GT. 0  .AND.  IVF .GT. 0) VAR = VAR*VFS(IVF)
          P = 1.D0/VAR
          IF ( .NOT. ADDOBS(C, IC, LENG, EL, P, A, NX) ) THEN
            WRITE (LUNIT,669)
  669       FORMAT ('0PROFILE ERROR IN NORMAL')
          write (lunit,*) 'subs2, 24'
            CALL ABORT2
          ENDIF
        ELSEIF (KIND .EQ. 18) THEN

*** CORRELATED TYPE (DOPPLER)

          EL3(1) = -CMO
          DO 201 I = 1, LENG
            C3(1,I) = C(I)
  201     CONTINUE
          READ (IUO,END=777) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD,
     &                       IOBS, IVF, IAUX, IGRT
          IF (KIND .NE. 19) STOP 202
          EL3(2) = -CMO
          DO 202 I = 1, LENG
            C3(2,I) = C(I)
  202     CONTINUE
          READ (IUO,END=777) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD,
     &                       IOBS, IVF, IAUX, IGRT
          IF (KIND .NE. 20) STOP 203
          EL3(3) = -CMO
          DO 203 I = 1, LENG
            C3(3,I) = C(I)
  203     CONTINUE
          READ (IUO,END=777) ISIGNL, ( (COVECF(I,J), J = 1,3), I = 1,3 )
          IF (ISIGNL .NE. 2000) STOP 204
          IF (NVFTOT .GT. 0  .AND.  IVF .GT. 0) THEN
            DO 204 I = 1, 3
              DO 205 J = 1, 3
                COVECF(I,J) = COVECF(I,J)*VFS(IVF)
  205         CONTINUE
  204       CONTINUE
          ENDIF
          IF (.NOT. ADDCOR(C3, IC, EL3, COVECF, LENG, 3, A, NX, IFLAG) )
     &    THEN
            IF (IFLAG .EQ. 1) WRITE (LUNIT,661)
  661       FORMAT (' PROGRAMMER ERROR, DOPPLER COV MATRIX COVECF IS',
     &              ' NOT POSITIVE DEFINITE')
            IF (IFLAG .EQ. 2) WRITE (LUNIT,662)
  662       FORMAT (' PROGRAMMER ERROR, PROFILE ERROR IN NORMAL')
          write (lunit,*) 'subs2, 25'
            CALL ABORT2
          ENDIF

        ELSE
          WRITE (LUNIT,664) KIND
  664     FORMAT (/,' PROGRAMMER ERROR IN NORMAL, DOPPLER KIND=', I3,
     &              ' OUT OF ORDER.')
          write (lunit,*) 'subs2, 26'
          CALL ABORT2
        ENDIF

*** DECORRELATED TYPE (GPS)

      ELSE
        NVEC = ISN
        IF ( .NOT. L2HLF) THEN
          LENG = (NVEC+1)*IDIM + 1
        ELSE
          LENG = (NVEC+1)*3 + 1
        ENDIF
        IF (NGRT .GT. 0) LENG = LENG + 3
        NR = 3*NVEC
*       NC = NR + 3 + LENG
        NC = NR + 5 + LENG
	
        READ (IUO,END=666) ICM, NICM, KINDS, ISNS, JSNS, LOBS
        CALL GLOCAT (N1, N2, N3, N4, N5, NVEC, LENG)
        IF (NVFTOT .GT. 0  .AND.  IVF .GT. 0) THEN
          P = 1.D0/VFS(IVF)
        ELSE
          P = 1.D0
        ENDIF
        DO 1 I = 1, NR
          READ (IUO,END=666) (G(J), J = 1, NC)
          EL = -G(N4)
          IF ( .NOT. ADDOBS(G(1+N2), ICM, NICM, EL, P, A, NX) ) THEN
            WRITE (LUNIT,10)
   10       FORMAT ('0ACCUM ERROR IN NORMAL')
          write (lunit,*) 'subs2, 27'
            CALL ABORT2
          ENDIF
    1   CONTINUE
      ENDIF
      GO TO 100

  777 CONTINUE

*** INITIALIZE GOOGE NUMBERS

      DO 2 I = 1, NUNK
        IF ( .NOT. GETA(I, I, VAL, A, NX) ) THEN
          CALL INVIUN (I, ISN, ITYP, ICODE)
          WRITE (LUNIT,3) I, ISN, ITYP
    3     FORMAT ('0FATAL GOOGE ERROR IN NORMALS',3I8)
          write (lunit,*) 'subs2, 28'
          CALL ABORT2
        ELSE
          GOOGE(I) = VAL
        ENDIF
    2 CONTINUE

*** FLUSH STANDARD OUTPUT BUFFER (SYSTEM DEPENDENT)
C***
C***      CALL MYFLSH
C***
*** SOLVE THE NORMALS, DO NOT ABORT ON SINGULARITIES!

      IF ( .NOT. SOLVE(A, NX, STOL, ISING, GSING, LSING) ) THEN
        IF (LSING .LE. 0) THEN
          WRITE (LUNIT,668)
  668     FORMAT ('0HOG STATE ERROR IN NORMAL')
          write (lunit,*) 'subs2, 29'
          CALL ABORT2
        ELSE
          CALL SINGUL (ITER, B, ISING, GSING, LSING, STOL)

*** INCREASE THE NUMBER OF CONSTRAINTS BY LSING

          IF (ITER .EQ. 0) THEN
            NCON = NCON + LSING
            NCON0 = LSING
          ELSE
            IF (LSING .NE. NCON0) THEN
              NCON = NCON - NCON0 + LSING
              NCON0 = LSING
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      RETURN

*** PREMATURE END OF FILE

  666 WRITE (LUNIT,667) NVEC
  667 FORMAT ('0PREMATURE FILE END IN NORMAL -- NVEC =',I8)
          write (lunit,*) 'subs2, 30'
      CALL ABORT2
      RETURN
      END
**********************************
** v 4.30vf - sub obseqw removed
***********************************
      SUBROUTINE OBSSUM (IUO, G, IUO3)

c Mike Potterfield 2/10/06
c The scratch file IUO3 contains flags to identify
c rejected vectors, which are no longer being counted
c in the observational summary

*** PRINTOUT THE OBSERVATIONAL SUMMARY

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, NVECS = 700, LENC = 10 )
      PARAMETER ( MAXZ =  8000, MXPRM = 40 )
      PARAMETER ( MXALL = MXPRM + 15 + MAXZ + MXSSN*3)
      CHARACTER*1 ACON
      CHARACTER*30 NAMES, NAME1
      CHARACTER*1 CREJCT
      CHARACTER*13 TPRJID
      LOGICAL LOCSSN, FATAL
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      DIMENSION ICM((NVECS+1)*3+1+3)
      DIMENSION KINDS(NVECS*3), ISNS(NVECS*3), JSNS(NVECS*3)
      DIMENSION LOBS(NVECS*3)
      DIMENSION G(*), IC(LENC), C(LENC)
      DIMENSION NCONS(MXSSN)
      DIMENSION COVECF(3,3)
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /OBSUM/  NSUM(MXSSN,13), ICMPL(MXALL)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

*** PRINT HEADING

      CALL HEAD6

*** INITIALIZE TABLE

      CALL NEWI

*** INITIALIZE CONSTRAINT ARRAY

      DO 20 I = 1, NSTA
        NCONS(I) = 0
   20 CONTINUE

      REWIND IUO
  100 READ (IUO,END=777) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD,
     &                   IOBS, IVF, IAUX, IGRT

      IF ( (KIND .GE. 4   .AND.  KIND .LE. 17)  .OR.
     &     (KIND .GE. 21  .AND.  KIND .LE. 999) ) THEN
              CALL CONTR (KIND, ISN, JSN, IAUX)
      ELSEIF (KIND .GE. 1  .AND.  KIND .LE. 3) THEN
        NCONS(ISN) = NCONS(ISN) + 1
      ELSEIF (KIND .GE. 1000  .AND.  KIND .LE. 1999) THEN
        NVEC = ISN
        IAUX = JSN
        IF ( .NOT. L2HLF) THEN
          LENG = (NVEC+1)*IDIM + 1
        ELSE
          LENG = (NVEC+1)*3 + 1
        ENDIF
        IF (NGRT .GT. 0) LENG = LENG + 3
        NR = 3*NVEC
*       NC = NR + 3 + LENG
        NC = NR + 5 + LENG

        READ (IUO,END=777) ICM, NICM, KINDS, ISNS, JSNS, LOBS
        READ (IUO3) TPRJID
        DO 10 I2 = 1, NVEC
          READ (IUO3) CREJCT

c Mike Potterfield 2/10/06
c The variable CREJECT is set to non-blank if the vector is rejected

          J1 = (I2-1)*3 + 1
          KIND = KINDS(J1)
          ISN = ISNS(J1)
          JSN = JSNS(J1)
          IF (CREJCT .EQ. ' ') THEN
               CALL CONTR (KIND, ISN, JSN, IAUX)
c          ELSE
c               WRITE(LUNIT,*) 'SKIPPING REJECTED OBSERVATION'
          ENDIF
   10   CONTINUE
        CALL DUMRD (IUO, NR, NC, G)
      ELSEIF (KIND .GE. 18  .AND.  KIND .LE. 20) THEN
        CALL CONTR (KIND, ISN, JSN, IAUX)
        READ (IUO,END=777) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD,
     &                     IOBS, IVF, IAUX, IGRT
        READ (IUO,END=777) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD,
     &                     IOBS, IVF, IAUX, IGRT
        READ (IUO,END=777) KIND, COVECF
      ENDIF
      GO TO 100

  777 FATAL = .FALSE.
      DO 1 I = 1, NSTA
        IF ( .NOT. LOCSSN(I,ISSN) ) THEN
          WRITE (LUNIT,666) I
  666     FORMAT ('0SSN TABLE ERROR IN OBSSUM--',I5)
          write (lunit,*) 'subs2, 31'
          CALL ABORT2
        ENDIF
        NAME1 = NAMES(I)

*** TEST IF STATION IS UNDETER, NO-CHECK, OR CONSTRAINED

        N3DOBS = NSUM(I,7) + NSUM(I,8) + NSUM(I,13)
        IF (IDIM .EQ. 1) THEN
          N1 = NSUM(I, 3) + NSUM(I, 4) + NSUM(I, 5) + NSUM(I, 6) +
     &         N3DOBS + NCONS(I)
          IF (N1 .LE. 1  .AND.  N3DOBS .NE. 1) THEN
            ACON = 'U'
            FATAL = .TRUE.
          ELSEIF ( N1 .EQ. 2  .OR.
     &            (N1 .EQ. 1  .AND.  N3DOBS .EQ. 1) ) THEN
            ACON = 'N'
          ELSEIF ( NCONS(I) .EQ. 3) THEN
            ACON = 'C'
          ELSE
            ACON = ' '
          ENDIF
          IF (ACON .EQ. 'N'  .AND.  N3DOBS .EQ. 2) ACON = ' '
          IUN = IUNSTA(I, 3)
          ICMP = ICMPL(IUN)
        ELSEIF (IDIM .EQ. 2) THEN
          N2 = NSUM(I, 1) + NSUM(I, 2) + NSUM(I, 3) + NSUM(I, 4) +
     &         NSUM(I, 9) + NSUM(I,10) + NSUM(I,11) + NSUM(I,12) +
     &         N3DOBS + NCONS(I)
*         IF ( NSUM(I,9) .GT. 0) N2 = N2 - 1
          IF ( N2 .LE. 1  .AND.  N3DOBS .NE. 1 ) THEN
            ACON = 'U'
            FATAL = .TRUE.
          ELSEIF ( N2 .EQ. 2  .OR.
     &            (N2 .EQ. 1  .AND.  N3DOBS .EQ. 1) ) THEN
            ACON = 'N'
          ELSEIF (NCONS(I) .EQ. 3) THEN
            ACON = 'C'
          ELSE
            ACON = ' '
          ENDIF
          IF (ACON .EQ. 'N'  .AND.  N3DOBS .EQ. 2) ACON = ' '
          IUN = IUNSTA(I, 1)
          ICMP = ICMPL(IUN)
        ELSE
          N3 = NSUM(I, 1) + NSUM(I, 2) + NSUM(I, 3) + NSUM(I, 4) +
     &         NSUM(I, 9) + NSUM(I,10) + NSUM(I,11) + NSUM(I,12) +
     &         NSUM(I, 5) + NSUM(I, 6) + N3DOBS + NCONS(I)
          IF (N3 .LE. 1  .AND.  N3DOBS .NE. 1) THEN
            ACON = 'U'
            FATAL = .TRUE.
          ELSEIF ( N3 .EQ. 2  .OR.
     &            (N3 .EQ. 1  .AND.  N3DOBS .EQ. 1) )THEN
            ACON = 'N'
          ELSEIF (NCONS(I) .EQ. 3) THEN
            ACON = 'C'
          ELSE
            ACON = ' '
          ENDIF
          IF (ACON .EQ. 'N'  .AND.  N3DOBS .EQ. 2) ACON = ' '
          IUN = IUNSTA(I, 1)
          ICMP = ICMPL(IUN)
        ENDIF

        IF ( LOS  .OR.  ( .NOT. LOS  .AND.  ACON .NE. ' ') ) THEN
          CALL LINE6 (1)
          WRITE (LUNIT,5) ISSN, ACON, ICMP, NAME1,
     &                NSUM(I, 9), NSUM(I,10), NSUM(I,11), NSUM(I,12),
     &                NSUM(I, 1), NSUM(I, 2), NSUM(I, 3), NSUM(I, 4),
     &                NSUM(I, 5), NSUM(I, 6), NSUM(I, 7), NSUM(I, 8),
     &                NSUM(I,13)
    5     FORMAT (' ', I5, A1, I4, 2X, A30, I4, I5, 5(I7, I5), I7 )
        ENDIF

    1 CONTINUE

*** NO LONGER TERMINATE IF UNDETERMINED STATIONS EXIST

      IF (FATAL) THEN
        WRITE (LUNIT,667)
  667   FORMAT ('0***** WARNING -- UNDETERMINED(U) STATIONS *****')
          write (lunit,*) 'subs2, 32'
***     CALL ABORT2
      ENDIF

      RETURN
      END
      SUBROUTINE PRTMAP

*** PRINT SINGULARITY MAP

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /EXTRMA/ GLAMAX, GLAMIN, GLOMAX, GLOMIN, HTMAX, HTMIN
      COMMON /CONST/  PI, PI2, RAD
      COMMON/UNITS/ LUNIT

*** PRINT HEADER

      CALL LINE (56)
      WRITE (LUNIT, 160)
  160 FORMAT ('0', 10X, 10('*'), '  SINGULARITY MAP  ', 10('*'), /)

*** PRINT MAP

      CALL GETPLT (0)

*** GET MAP BOUNDARIES IN DEGREES

      YMIN = GLAMIN*RAD
      YMAX = GLAMAX*RAD
      XMIN = GLOMIN*RAD
      XMAX = GLOMAX*RAD
      IF (XMIN .GT. 180) XMIN = XMIN - 360.D0
      IF (XMAX .GT. 180) XMAX = XMAX - 360.D0

*** PRINT MAP INFORMATION

      WRITE (LUNIT, 170) YMIN, YMAX, XMIN, XMAX
  170 FORMAT (/, ' MINIMUM LATITUDE =', F13.7,
     &           '   MAXIMUM LATITUDE =', F13.7,
     &        /, ' MINUMUM LONGITUDE=', F13.7,
     &           '   MAXIMUM LONGITUDE=', F13.7)

      RETURN
      END
      SUBROUTINE PUTALA (ALA,I,B)

*** ROUTINE TO INSERT ASTRO LAT INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      DIMENSION B(*)
      LOGICAL L2HLF, LEHT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTALA')
          write (lunit,*) 'subs2, 33'
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      B(1+N) = ALA

      RETURN
      END
      SUBROUTINE PUTALO (ALO,I,B)

*** ROUTINE TO INSERT ASTRO LON INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      DIMENSION B(*)
      LOGICAL L2HLF, LEHT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTALO')
          write (lunit,*) 'subs2, 34'
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      B(2+N) = ALO

      RETURN
      END
      SUBROUTINE PUTAUX (AUX,I,B)

*** ROUTINE TO INSERT AUX. PARM. INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NAUX) THEN
        WRITE (LUNIT,1) I,NAUX
    1   FORMAT ('0****ILLEGAL IAUX',I5,' FOR NAUX=',I5,' IN PUTAUX')
          write (lunit,*) 'subs2, 35'
        CALL ABORT2
      ENDIF

      B(I) = AUX

      RETURN
      END
      SUBROUTINE PUTECX (ECX,I,B)

*** ROUTINE TO INSERT E.C.F. X INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      DIMENSION B(*)
      LOGICAL L2HLF, LEHT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTECX')
          write (lunit,*) 'subs2, 36'
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      B(5+N) = ECX

      RETURN
      END
      SUBROUTINE PUTECY (ECY,I,B)

*** ROUTINE TO INSERT E.C.F. Y INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      DIMENSION B(*)
      LOGICAL L2HLF, LEHT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTECY')
          write (lunit,*) 'subs2, 37'
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      B(6+N) = ECY

      RETURN
      END
      SUBROUTINE PUTECZ (ECZ,I,B)

*** ROUTINE TO INSERT E.C.F. Z INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      DIMENSION B(*)
      LOGICAL L2HLF, LEHT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTECZ')
          write (lunit,*) 'subs2, 38'
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      B(7+N) = ECZ

      RETURN
      END
C--------------------------------------------------------------------------------------------------
      SUBROUTINE PUTEHT (EHT,I,B)

*** ROUTINE TO INSERT ELLIPSOIDAL HT. INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTGH ')
        write (lunit,*) 'subs2, 39'
        CALL ABORT2
      ENDIF
      N = NAUX+NZ+(I-1)*13
      B(10+N) = EHT

      RETURN
      END
C--------------------------------------------------------------------------------------------------
      SUBROUTINE PUTEHX (EHX,I,B)

*** ROUTINE TO INSERT ELLIP. HT. X INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTEHX')
          write (lunit,*) 'subs2, 40'
        CALL ABORT2
      ENDIF
      N = NAUX+NZ+(I-1)*13
      B(11+N) = EHX

      RETURN
      END
      SUBROUTINE PUTEHY (EHY,I,B)

*** ROUTINE TO INSERT ELLIP. HT. Y INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTEHY')
          write (lunit,*) 'subs2, 41'
        CALL ABORT2
      ENDIF
      N = NAUX+NZ+(I-1)*13
      B(12+N) = EHY

      RETURN
      END
      SUBROUTINE PUTEHZ (EHZ,I,B)

*** ROUTINE TO INSERT ELLIP. HT. Z INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTEHZ')
          write (lunit,*) 'subs2, 42'
        CALL ABORT2
      ENDIF
      N = NAUX+NZ+(I-1)*13
      B(13+N) = EHZ

      RETURN
      END
      SUBROUTINE PUTGH (GH,I,B)

*** ROUTINE TO INSERT GEOID HT. INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      DIMENSION B(*)
      LOGICAL L2HLF, LEHT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTGH ')
          write (lunit,*) 'subs2, 43'
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      B(9+N) = GH

      RETURN
      END
C----------------------------------------------------------------------------------------------------------
      SUBROUTINE PUTGLA (GLA,I,B)

*** ROUTINE TO INSERT GEOD. LAT INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      DIMENSION B(*)
      LOGICAL L2HLF, LEHT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTGLA')
          write (lunit,*) 'subs2, 44'
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      B(3+N) = GLA

      RETURN
      END
C-----------------------------------------------------------------------------------------------------------
      SUBROUTINE PUTGLO (GLO,I,B)

*** ROUTINE TO INSERT GEOD. LON INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      DIMENSION B(*)
      LOGICAL L2HLF, LEHT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTGLO')
          write (lunit,*) 'subs2, 45'
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      B(4+N) = GLO

      RETURN
      END
C-----------------------------------------------------------------------------------------------------------
      LOGICAL FUNCTION PUTGRT (ICODE,IOLD,INEW,IDUP)

*** ADD ENTRY TO TABLE
*** IDUP IS LOCATION OF DUPLICATE ENTRY IN LIST
*** IDUP = 0 INDICATES EXCEEDED MAXIMUM LENGTH OF LIST

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL GETGRT
      COMMON /PRMTB2/ NPARMS, ICODES(5), IOLDS(5), INEWS(5)

      IF (GETGRT(ICODE,IOLD,IPARM)) THEN
        IDUP = IPARM
        PUTGRT = .FALSE.
        RETURN
      ELSEIF (GETGRT(ICODE,INEW,IPARM)) THEN
        IDUP = IPARM
        PUTGRT = .FALSE.
        RETURN
      ELSE
        IDUP = 0
        IF (IPARM.GT.5) THEN
          PUTGRT = .FALSE.
          RETURN
        ELSE
          NPARMS = IPARM
          ICODES(NPARMS) = ICODE
          IOLDS(NPARMS) = IOLD
          INEWS(NPARMS) = INEW
          PUTGRT = .TRUE.
        ENDIF
      ENDIF

      RETURN
      END
      SUBROUTINE PUTICM (IC, L, ICMS, NICMS, NVECS)

*** ADD ENTRIES OF IC INTO TABLE

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      LOGICAL GETICM
      DIMENSION ICMS(*)
      DIMENSION IC(LENC)
      COMMON/UNITS/ LUNIT

      MAXICM = (NVECS+1)*3 + 1 + 3

      DO 1 I = 1, L
        ICM = IC(I)
        IF ( GETICM(ICM, IICM, ICMS, NICMS) ) THEN
          CONTINUE
        ELSE
          IF (IICM .GT. MAXICM) THEN
            WRITE (LUNIT,666)
  666       FORMAT ('0GPS INDEX OVERFLOW IN PUTICM')
          write (lunit,*) 'subs2, 46'
            CALL ABORT2
          ELSE
            NICMS = IICM
            ICMS(NICMS) = ICM
          ENDIF
        ENDIF
    1 CONTINUE

      RETURN
      END
      LOGICAL FUNCTION PUTIVF (ICODE, IOLD, INEW, LFIX, IDUP)

*** ADD ENTRY TO TABLE
*** IDUP IS LOCATION OF DUPLICATE ENTRY IN LIST
*** IDUP = 0 INDICATES EXCEEDED MAXIMUM LENGTH OF LIST

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXVF = 40 )
      LOGICAL GETIVF
      LOGICAL LFIXS, LFIX
      COMMON /IVFTB1/ NIVFS, ICODES(MXVF), IOLDS(MXVF), INEWS(MXVF),
     &                LFIXS(MXVF)
      COMMON /NUMVFS/ NVFTOT, NVFREE

      IF ( GETIVF(ICODE, IOLD, IVF) ) THEN
        IDUP = IVF
        PUTIVF = .FALSE.
        RETURN
      ELSEIF ( GETIVF(ICODE, INEW, IVF) ) THEN
        IDUP = IVF
        PUTIVF = .FALSE.
        RETURN
      ELSE
        IDUP = 0
        IF (IVF.GT.MXVF) THEN
          PUTIVF = .FALSE.
          RETURN
        ELSE
          NIVFS = IVF
          ICODES(NIVFS) = ICODE
          IOLDS(NIVFS) = IOLD
          INEWS(NIVFS) = INEW
          LFIXS(NIVFS) = LFIX
          NVFTOT = NVFTOT + 1
          IF (.NOT.LFIX) NVFREE = NVFREE + 1
          PUTIVF = .TRUE.
        ENDIF
      ENDIF

      RETURN
      END
      SUBROUTINE PUTMSL (GMSL,I,B)

*** ROUTINE TO INSERT MSL INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      DIMENSION B(*)
      LOGICAL L2HLF, LEHT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN PUTMSL')
          write (lunit,*) 'subs2, 47'
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      B(8+N) = GMSL

      RETURN
      END
      LOGICAL FUNCTION PUTPRM (ICODE, IOLD, INEW, IDUP)

*** ADD ENTRY TO TABLE
*** IDUP IS LOCATION OF DUPLICATE ENTRY IN LIST
*** IDUP = 0 INDICATES EXCEEDED MAXIMUM LENGTH OF LIST

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXPRM = 40 )
      LOGICAL GETPRM
      COMMON /PRMTB1/ NPARMS, ICODES(MXPRM), IOLDS(MXPRM), INEWS(MXPRM)

      IF ( GETPRM(ICODE, IOLD, IPARM) ) THEN
        IDUP = IPARM
        PUTPRM = .FALSE.
        RETURN
      ELSEIF ( GETPRM(ICODE, INEW, IPARM) ) THEN
        IDUP = IPARM
        PUTPRM = .FALSE.
        RETURN
      ELSE
        IDUP = 0
        IF (IPARM .GT. MXPRM) THEN
          PUTPRM = .FALSE.
          RETURN
        ELSE
          NPARMS = IPARM
          ICODES(NPARMS) = ICODE
          IOLDS(NPARMS) = IOLD
          INEWS(NPARMS) = INEW
          PUTPRM = .TRUE.
        ENDIF
      ENDIF

      RETURN
      END
      SUBROUTINE PUTRHS (G, NR, NC, N, M)

*** COPY N-TH COLUMN TO M-TH COLUMN

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION G(NR,NC)

      DO 1 I = 1, NR
        G(I,M) = G(I,N)
    1 CONTINUE

      RETURN
      END
      SUBROUTINE PUTROT (ROT,I,B)

*** ROUTINE TO INSERT ROTATION PARM INTO CONTROL PT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NZ) THEN
        WRITE (LUNIT,1) I,NZ
    1   FORMAT ('0****ILLEGAL IZ',I5,' FOR NZ=',I5,' IN PUTROT')
          write (lunit,*) 'subs2, 48'
        CALL ABORT2
      ENDIF

      B(I+NAUX) = ROT

      RETURN
      END
      SUBROUTINE PUTRTG (GRTX,GRTY,GRTZ,I)

*** ROUTINE TO INSERT AUX. GPS AND DOPPLER ROTATION PARAMETERS
*** INTO CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /GRTTB2/ BGRT(5,3)
      COMMON/UNITS/ LUNIT

      IF ( I.LE.0 .OR. I.GT.NGRT ) THEN
        WRITE (LUNIT,1) I,NGRT
    1   FORMAT ('0****ILLEGAL IGRT=',I5,' FOR NGRT=',I5,' IN PUTRTG')
          write (lunit,*) 'subs2, 49'
        CALL ABORT2
      ENDIF

      BGRT(I,1) = GRTX
      BGRT(I,2) = GRTY
      BGRT(I,3) = GRTZ

      RETURN
      END
      LOGICAL FUNCTION PUTSSN (ISSN,IDUP)

*** ADD ENTRY TO TABLE
*** IDUP IS LOCATION OF DUPLICATE ENTRY IN LIST
*** IDUP = 0 INDICATES ILLEGAL ISSN

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      COMMON /SSNTBL/ NSSN, ISSISN(MXSSN)
      LOGICAL GETSSN

      IF (ISSN.LE.0 .OR. ISSN.GT.MXSSN) THEN
        IDUP = 0
        PUTSSN = .FALSE.
      ELSEIF (GETSSN(ISSN,ISN)) THEN
        IDUP = ISN
        PUTSSN = .FALSE.
      ELSE
        IDUP = 0
        NSSN = NSSN+1
        ISSISN(ISSN) = NSSN
        PUTSSN = .TRUE.
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION PUTZ (ISSN,LIST,IDUP)

*** ADD ENTRY TO TABLE
*** IDUP IS LOCATION OF DUPLICATE ENTRY IN LIST
*** IDUP = 0 INDICATES TABLE OVERFLOW

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MAXZ =  8000 )
      COMMON /ROTTB1/ IZNS(MAXZ)
      COMMON /ROTTB2/ LOOK(977), IPTRS(MAXZ), NKEY, MAX, IFREE
      LOGICAL LOKZ

      IZN = ISSN*1000+LIST
      IF (LOKZ(IZN,ILOOK,IPTR)) THEN

*** DUPLICATE ENTRY

        IDUP = IPTR
        PUTZ = .FALSE.
      ELSE
        IF (IFREE.GT.MAX) THEN

*** TABLE OVERFLOW

          IDUP = 0
          PUTZ = .FALSE.
        ELSE

*** PUT ENTRY IN TABLE
          IF (IPTR.EQ.0) THEN
            LOOK(ILOOK) = IFREE
          ELSE
            IPTRS(IPTR) = IFREE
          ENDIF
          IZNS(IFREE) = IZN
          IPTRS(IFREE) = 0
          IFREE = IFREE+1
          PUTZ = .TRUE.
        ENDIF
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION PUT82 (ISSN,IDUP)

*** ADD ENTRY TO THE *82* RECORD TABLE
*** IDUP IS LOCATION OF DUPLICATE ENTRY IN LIST
*** IDUP = 0 INDICATES ILLEGAL ISSN

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      COMMON /CHILD/ N82,ISS82(MXSSN)
      LOGICAL GET82

      IF (ISSN.LE.0 .OR. ISSN.GT.MXSSN) THEN
        IDUP = 0
        PUT82 = .FALSE.
      ELSEIF (GET82(ISSN,ISN)) THEN
        IDUP = ISN
        PUT82 = .FALSE.
      ELSE
        IDUP = 0
        N82  = N82  + 1
        ISS82(ISSN) = N82
        PUT82 = .TRUE.
      ENDIF

      RETURN
      END
      SUBROUTINE RADII (ISN,RMER,RPV,B)

*** COMPUTE RADII OF CURVATURE
*** SEE RAPP, GEOMETRIC GEOD. VOL I, P 19 AND 24

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      DIMENSION B(*)
      COMMON /OPT/ AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &             LMSL, LSS, LUP, LABS, LADJ, LUPI

      CALL GETGLA (GLAT,ISN,B)
      SLAT = DSIN(GLAT)
      SLAT2 = SLAT*SLAT
      W = DSQRT(1.D0-E2*SLAT2)

*** RADIUS OF CURVATURE IN MERIDIAN

      RMER = AX*(1.D0-E2)/(W*W*W)

*** RADIUS OF CURVATURE IN PRIME VERTICAL

      RPV = AX/W

      RETURN
      END
      SUBROUTINE RDOP (IUO, KIND, IOBSX, ISN, JSN, OBSBX, SIGX, B, IAUX,
     &                 IVF, SIGUWT, IC, C, LENG, A, NX, IGRT)

*** LIST ADJUSTED DOPPLER OBSERVATIONS AND RESIDUALS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, LENC = 10, MXVF = 40 )
      DIMENSION A(*), B(*), NX(*)
      DIMENSION IC(LENC), C(LENC)
      DIMENSION C3(3,LENC), COVECF(3,3)
      DIMENSION COVLA(3,3), Q(3,3), COVLBI(3,3), DUMMAT(3,3)
      CHARACTER*30 NAMES,NAME
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      logical lc
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON /VFTB2/  VFS(MXVF), VTV(MXVF), VFRN(MXVF), VFSS(MXVF)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON/UNITS/ LUNIT
      SAVE NAME, ISNX
      save lc
      DATA ISNX /0/

* v 4.28i lc = true if normalized residual exceeds tolerance, print

      IF (KIND .NE. 18) THEN
        WRITE (LUNIT,1) KIND
    1   FORMAT ('0PROGRAMMER ERROR IN RDOP - KIND = ',I5)
          write (lunit,*) 'subs2, 50'
        CALL ABORT2
      ENDIF

      DO 101 I = 1, LENG
        C3(1,I) = C(I)
  101 CONTINUE
      READ (IUO,END=777) KIND, ISN, JSN, IC, C, LENG, CMO, OBSBY, SIGY,
     &                   IOBSY, IVF, IAUX, IGRT
      IF (KIND .NE. 19) STOP 102
      DO 102 I = 1, LENG
        C3(2,I) = C(I)
  102 CONTINUE
      READ (IUO,END=777) KIND, ISN, JSN, IC, C, LENG, CMO, OBSBZ, SIGZ,
     &                   IOBSZ, IVF, IAUX, IGRT
      IF (KIND .NE. 20) STOP 103
      DO 103 I = 1, LENG
        C3(3,I) = C(I)
  103 CONTINUE
      READ (IUO,END=777) ISIGNL, ( (COVECF(I,J), J = 1,3), I = 1,3 )
      IF (ISIGNL .NE. 2000) STOP 104

*** CHECK FOR NEW ISN

      IF (ISNX .NE. ISN) THEN
        ISNX = ISN
	if(lc) then
        WRITE (LUNIT,2)
 2      format(' ')	
        CALL LINE2 (1)
	endif
        NAME = NAMES(ISN)
        IF (IMODE .EQ. 3) THEN
          CALL LINE2 (1)
          WRITE (LUNIT,510) NAME
 510      FORMAT (93X,A30)
        ENDIF
      ENDIF

*** SCALE BY THE VARIANCE FACTOR

      IF (NVFTOT .GT. 0) THEN
        IF (IVF .LT. 0  .OR.  IVF .GT. NVFTOT) THEN
          WRITE (LUNIT,9) IVF,NVFTOT
    9     FORMAT ('0ILLEGAL IVF=',I10,' FOR N=',I10,' IN RDOP')
          write (lunit,*) 'subs2, 51'
          CALL ABORT2
        ELSEIF (IVF .NE. 0) THEN
          SIGFS = DSQRT(VFS(IVF))
          SIGX = SIGX*SIGFS
          SIGY = SIGY*SIGFS
          SIGZ = SIGZ*SIGFS
          DO 104 I = 1, 3
            DO 105 J = 1, 3
              COVECF(I,J) = COVECF(I,J)*SIGFS*SIGFS
  105       CONTINUE
  104     CONTINUE
        ENDIF
      ENDIF

*** COMPUTATION OF COVARIANCE MATRIX FOR ADJUSTED OBS,
*** STD. DEV. OF RESIDUALS, AND REDUNDANCE NUMBER

      IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) THEN

        CALL COMSLA ( A, NX, IC, C3, LENG, COVLA)
        DVARX = COVECF(1,1) - COVLA(1,1)
        SD3X = DSQRT(DVARX)
        IF (DVARX .LT. 1.0D-10) DVARX = 0.D0
        DVARY = COVECF(2,2) - COVLA(2,2)
        IF (DVARY .LT. 1.0D-10) DVARY = 0.D0
        SD3Y = DSQRT(DVARY)
        DVARZ = COVECF(3,3) - COVLA(3,3)
        IF (DVARZ .LT. 1.0D-10) DVARZ = 0.D0
        SD3Z = DSQRT(DVARZ)

*** COMPUTATION OF REDUNDENCY NUMBER

        CALL EQ ( COVECF, COVLBI, 3, 3)
        CALL INVRT3 (COVLBI)
        CALL EQ (COVLBI, DUMMAT, 3, 3)
        CALL GETQ (DUMMAT, COVLA, Q)

        RN = Q(1,1) + Q(2,2) + Q(3,3)
        IF (RN .LT. 1.0D-10) RN = 0.D0
      ENDIF
      IF (IMODE .GE. 0  .AND.  IMODE .LE. 2) THEN
        SD3X = SIGX
        SD3Y = SIGY
        SD3Z = SIGZ
      ENDIF

*** COMPUTATION OF MARGINAL DETECTABLE ERROR(GMD) USING LAMDA=3

      IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) THEN
        IF (RN .LT. 1.0D-8) THEN
          GMDEX = 0.D0
          GMDEY = 0.D0
          GMDEZ = 0.D0
        ELSE
          CALL AB ( COVLBI, Q, DUMMAT, 3, 3, 3)
          GMDEX = 3.D0/DSQRT( DUMMAT(1,1) )
          GMDEY = 3.D0/DSQRT( DUMMAT(2,2) )
          GMDEZ = 3.D0/DSQRT( DUMMAT(3,3) )
        ENDIF
      ENDIF

*** DOPPLER X OBSERVATION

      KIND = 18
      CALL COMPOB (KIND, X, B, OBSBX, ISN, JSN, IAUX, IGRT)
      VX = X - OBSBX
      VSDX = VX/SIGX
      VSD3X = DIVID(VX, SD3X)
      IF (LSS) VSDX = VSDX/SIGUWT
      
****** v 4.28n
      IF (LSS) SD3X = SD3X/SIGUWT
      IF (LSS) VSD3X = VSD3X/SIGUWT
*******************

      IF (LVX) THEN
        IF (DABS(VSD3X) .GE. CRIT) THEN
* v 4.28i
          lc = .true.	
          CALL LINE2 (1)
          IF (IMODE .EQ. 0) THEN
            WRITE (LUNIT,205) IOBSX, X, GMDEX, NAME
 205        FORMAT (' ', I6, ' DOPX', F13.4, F11.3, 20X, A30 )
          ELSEIF (IMODE .EQ. 3) THEN
            WRITE (LUNIT,210) IOBSX, X, OBSBX, VX, SD3X, VSD3X, GMDEX 
 210        FORMAT (' ', I6, ' DOPX', F13.4, F17.4, F19.3, F8.2, F8.2,
     &              F10.4)
          ELSE
            WRITE (LUNIT,215) IOBSX, X, OBSBX, VX, VSDX, NAME
 215        FORMAT (' ', I6, ' DOPX', F13.3, F17.3, F21.3,
     &              F8.1, 1X, A30 )
          ENDIF
	else
	  lc = .false.
        ENDIF
      ENDIF

*** ACCUMULATE STATISTICS

      IF (LSS) VSDX = VSDX*SIGUWT
      
*********v 4.28n
      IF (LSS) SD3X = SD3X*SIGUWT
      IF (LSS) VSD3X = VSD3X*SIGUWT
******************
c IOBSX is being passed to RSTAT so that the observation
c residual summaries will show the correct observation numbers.
c Mike Potterfield 3/15/07

      CALL RSTAT (VX, VSDX, Q(1,1), VSD3X, IOBSX, KIND, IOBSX)

*** DOPPLER Y OBSERVATION

      KIND = 19
      CALL COMPOB (KIND, Y, B, OBSBY, ISN, JSN, IAUX, IGRT)
      VY = Y - OBSBY
      VSDY = VY/SIGY
      VSD3Y = DIVID(VY, SD3Y)
      IF (LSS) VSDY = VSDY/SIGUWT
      
****** v 4.28n
      IF (LSS) SD3Y = SD3Y/SIGUWT
      IF (LSS) VSD3Y = VSD3Y/SIGUWT
*******************
      
      IF (LVX) THEN
        IF (DABS(VSD3Y) .GE. CRIT) THEN
	  lc = .true.
          CALL LINE2 (1)
          IF (IMODE .EQ. 0) THEN
            WRITE (LUNIT,225) IOBSY, Y, GMDEY
 225        FORMAT (' ', I6, ' DOPY', F13.4, F11.3)
          ELSEIF (IMODE .EQ. 3) THEN
            WRITE (LUNIT,230) IOBSY, Y, OBSBY, VY, SD3Y, VSD3Y, GMDEY
 230        FORMAT (' ', I6, ' DOPY', F13.4, F17.4, F19.3, F8.2, F8.2,
     &              F10.4)
          ELSE
            WRITE (LUNIT,235) IOBSY, Y, OBSBY, VY, VSDY
 235        FORMAT (' ', I6, ' DOPY', F13.3, F17.3, F21.3, F8.1 )
          ENDIF
	else
	  lc = .false.
        ENDIF
      ENDIF

*** ACCUMULATE STATISTICS

      IF (LSS) VSDY = VSDY*SIGUWT
      
*********v 4.28n
      IF (LSS) SD3Y = SD3Y*SIGUWT
      IF (LSS) VSD3Y = VSD3Y*SIGUWT
******************
c IOBSX is being passed to RSTAT so that the observation
c residual summaries will show the correct observation numbers.
c Mike Potterfield 3/15/07
      
      CALL RSTAT (VY, VSDY, Q(2,2), VSD3Y, IOBSY, KIND, IOBSY)

*** DOPPLER Z OBSERVATION

      KIND = 20
      CALL COMPOB (KIND, Z, B, OBSBZ, ISN, JSN, IAUX, IGRT)
      VZ = Z - OBSBZ
      VSDZ = VZ/SIGZ
      VSD3Z = DIVID(VZ, SD3Z)
      IF (LSS) VSDZ = VSDZ/SIGUWT
      
****** v 4.28n
      IF (LSS) SD3Z = SD3Z/SIGUWT
      IF (LSS) VSD3Z = VSD3Z/SIGUWT
*******************
      
      IF (LVX) THEN
        IF (DABS(VSD3Z) .GE. CRIT) THEN
	  lc = .true.
          CALL LINE2 (1)
          IF (IMODE .EQ. 0) THEN
            WRITE (LUNIT,245) IOBSZ, Z, GMDEZ, RN
 245        FORMAT (' ', I6, ' DOPZ', F13.4, F11.3, F19.2 )
          ELSEIF (IMODE .EQ. 3) THEN
            WRITE (LUNIT,250)IOBSZ,Z,OBSBZ,VZ,SD3Z,VSD3Z,GMDEZ,RN
 250        FORMAT (' ', I6, ' DOPZ', F13.4, F17.4, F19.3, F8.2, F8.2,
     &              F10.4, F6.2)
          ELSE
            WRITE (LUNIT,255) IOBSZ, Z, OBSBZ, VZ, VSDZ
 255        FORMAT (' ', I6, ' DOPZ', F13.3, F17.3, F21.3,
     *              F8.1 )
          ENDIF
	else
	  lc = .false.
        ENDIF
      ENDIF

*** ACCUMULATE STATISTICS

      IF (LSS) VSDZ = VSDZ*SIGUWT
      
*********v 4.28n
      IF (LSS) SD3Z = SD3Z*SIGUWT
      IF (LSS) VSD3Z = VSD3Z*SIGUWT
******************
c IOBSX is being passed to RSTAT so that the observation
c residual summaries will show the correct observation numbers
c Mike Potterfield 3/15/07
      
      CALL RSTAT (VZ, VSDZ, Q(3,3), VSD3Z, IOBSZ, KIND, IOBSZ)
      IF ( IMODE .EQ. 3  .OR.  IMODE .EQ. 0) THEN
        CALL RSTAT3 (RN)
      ENDIF

      RETURN
  777 WRITE (LUNIT,778)
  778 FORMAT ('0PROGRAMMER ERROR IN RDOP - SHOULD NOT REACH THIS POINT')
      STOP 778
      END
      SUBROUTINE RELACC (ISSN, JSSN, ISN, JSN, A, NX, B, SHIFTS, SIGUWT)

*** COMPUTE RELATIVE ACCURACIES ON SPHERE

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, LENC = 10 )
      LOGICAL PROP
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL LCRIT
      CHARACTER*72 CFLAG
      CHARACTER*30 NAMES, NAME1, NAME2
      CHARACTER*2 CC34
      CHARACTER*1 A2
      DIMENSION A(*), NX(*), B(*), SHIFTS(*)
      DIMENSION IC(LENC), C(LENC)
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /GPSORD/ CC34
      COMMON/UNITS/ LUNIT

*** DISTANCE IS KIND 7

      IAUX = 0
      IGRT = 0
      KIND = 7
**v 4.28***************
      IX1 = 0
      IX2 = 0
**********************
     
      CALL FORMIC (KIND, ISN, JSN, IAUX, IC, LENG, IGRT)
      CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
      CALL COMPOB (KIND, DIST, B, DUMMY, ISN, JSN, IAUX, IGRT)
      IF ( .NOT. PROP(C, IC, LENG, VARDST, A, NX, IFLAG) ) THEN
        IF (IFLAG .EQ. 1) THEN
          WRITE (LUNIT,666)
  666     FORMAT ('0STATE ERROR IN RELACC')
          write (lunit,*) 'subs2, 52'
          CALL ABORT2
        ELSE
          WRITE (LUNIT,667)
  667     FORMAT ('0PROFILE ERROR IN RELACC')
          write (lunit,*) 'subs2, 53'
          CALL ABORT2
        ENDIF
      ENDIF

*** AZIMUTH IS KIND 8

      KIND = 8
      CALL FORMIC (KIND, ISN, JSN, IAUX, IC, LENG, IGRT)
      CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
      CALL COMPOB (KIND, AZM, B, DUMMY, ISN, JSN, IAUX, IGRT)
      CALL DIRDMS (AZM, ID1, IM1, S1)
      IF ( .NOT. PROP(C, IC, LENG, VARAZM, A, NX, IFLAG) ) THEN
        IF (IFLAG .EQ. 1) THEN
          WRITE (LUNIT,666)
          write (lunit,*) 'subs2, 54'
          CALL ABORT2
        ELSE
          WRITE (LUNIT,667)
          write (lunit,*) 'subs2, 55'
          CALL ABORT2
        ENDIF
      ENDIF

*** ZENITH DISTANCE IS KIND 9

      KIND = 9
      CALL FORMIC (KIND, ISN, JSN, IAUX, IC, LENG, IGRT)
      CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
      DO 10 I = 1, LENG
 10   C(I) = -C(I)
      CALL COMPOB (KIND, VERANG, B, DUMMY, ISN, JSN, IAUX, IGRT)
      VERANG = PI2 - VERANG
      CALL GETDMS (VERANG, ID2, IM2, S2, ISIGN)
      IF (ISIGN .GT. 0) THEN
        A2 = 'E'
      ELSE
        A2 = 'D'
      ENDIF
      IF ( .NOT. PROP(C, IC, LENG, VARVER, A, NX, IFLAG) ) THEN
        IF (IFLAG .EQ. 1) THEN
          WRITE (LUNIT,666)
          write (lunit,*) 'subs2, 56'
          CALL ABORT2
        ELSE
          WRITE (LUNIT,667)
          write (lunit,*) 'subs2, 57'
          CALL ABORT2
        ENDIF
      ENDIF

      SD = DSQRT(VARDST)
      SA = DSQRT(VARAZM)*RAD*3600.D0
      SV = DSQRT(VARVER)*RAD*3600.D0
      IF (LSS) THEN
        SD = SD*SIGUWT
        SA = SA*SIGUWT
        SV = SV*SIGUWT
      ENDIF
 
*     IREL1 = IDNINT(DIST/SD)
*     IF ( IREL1 .GT. 99999999 ) IREL1 = 0
      
      REL1 = DIST/SD
      
**v 4.28************************
*     IF(REL1.GE.1.D08) REL1 = 0.D0
*      IF(REL1.GE.1.D08) THEN        
**** v 4.28e
      IF(REL1.GE. 0.0) THEN        
        IX1 = 1 
      ELSE	
        IREL1 = IDNINT(REL1)
      ENDIF	
********************************
*** v 4.29c *****
*     REL1 = DIST/SD
*     IF(REL1.GE.1.D08) REL1 = 0.D0
*     IREL1 = IDNINT(REL1)
***************************
      
      DELX = SHIFTS( IUNSHF(JSN, 1) ) - SHIFTS( IUNSHF(ISN, 1) )
      DELY = SHIFTS( IUNSHF(JSN, 2) ) - SHIFTS( IUNSHF(ISN, 2) )
      SHIFT = DSQRT(DELX*DELX + DELY*DELY)
      IF (SHIFT .LT. 0.0001D0) SHIFT = 0.0001D0
      
*     IREL2 = IDNINT(DIST/SHIFT)
      
      REL2 = DIST/SHIFT
      
**v 4.28 *************************
***v 4.28e
*     IF(REL2.GE.1.D08) REL2 = 0.D0
      IF(REL2.GE.1.D08) THEN          
        IX2 = 1 
      ELSE	
        IREL2 = IDNINT(REL2)
      ENDIF	
************************************	
*** v 4.29c********
*     REL2 = DIST/SHIFT
*     IF(REL2.GE.1.D08) REL2 = 0.D0
*     IREL2 = IDNINT(REL2)
********************

      IF (IMODE .EQ. 0) IREL2 = 0
      NAME1 = NAMES(ISN)
      NAME2 = NAMES(JSN)
      SADST = DIST*SA/3600.D0/RAD
      SVDST = DIST*SV/3600.D0/RAD

*** CALCULATE CRITICAL VALUE -- FLAG IF GREATER THAN 1.96 SD

      IF (CC34 .NE. '  ') THEN
        LCRIT = .TRUE.
        CFLAG( 1:36) = '                                    '
        CFLAG(37:72) = '                                    '
        IF (CC34 .EQ. 'B ') THEN
          CRIT = 0.008D0 + 1.D-6*DIST
          IF (1.96*SD .GT. CRIT ) THEN
            CFLAG( 1:36) = ' M. ***** WARNING - LINE FAILS ORDER'
            CFLAG(37:72) = ' B CRITICAL VALUE TEST *****        '
          ENDIF
        ELSEIF (CC34 .EQ. 'A ') THEN
          CRIT = 0.005D0 + 1.D-7*DIST
          IF (1.96*SD .GT. CRIT ) THEN
            CFLAG( 1:36) = ' M. ***** WARNING - LINE FAILS ORDER'
            CFLAG(37:72) = ' A CRITICAL VALUE TEST *****        '
          ENDIF
        ELSEIF (CC34 .EQ. 'AA') THEN
          CRIT = 0.003D0 + 1.D-8*DIST
          IF (1.96*SD .GT. CRIT ) THEN
            CFLAG( 1:36) = ' M. ***** WARNING - LINE FAILS ORDER'
            CFLAG(37:72) = ' AA CRITICAL VALUE TEST *****       '
          ENDIF
        ELSEIF (CC34 .EQ. '1 ') THEN
          CRIT = 0.010D0 + 1.D-5*DIST
          IF (1.96*SD .GT. CRIT ) THEN
            CFLAG( 1:36) = ' M. ***** WARNING - LINE FAILS ORDER'
            CFLAG(37:72) = ' 1 CRITICAL VALUE TEST *****        '
          ENDIF
        ELSEIF (CC34 .EQ. '21'  .OR.  CC34 .EQ. '2 ') THEN
          CRIT = 0.020D0 + 2.D-5*DIST
          IF (1.96*SD .GT. CRIT ) THEN
            CFLAG( 1:36) = ' M. ***** WARNING - LINE FAILS ORDER'
            CFLAG(37:72) = ' 2-I CRITICAL VALUE TEST *****      '
          ENDIF
        ELSEIF (CC34 .EQ. '22') THEN
          CRIT = 0.030D0 + 5.D-5*DIST
          IF (1.96*SD .GT. CRIT ) THEN
            CFLAG( 1:36) = ' M. ***** WARNING - LINE FAILS ORDER'
            CFLAG(37:72) = ' 2-II CRITICAL VALUE TEST *****     '
          ENDIF
        ELSEIF (CC34 .EQ. '3 ') THEN
          CRIT = 0.050D0 + 1.D-4*DIST
          IF (1.96*SD .GT. CRIT ) THEN
            CFLAG( 1:36) = ' M. ***** WARNING - LINE FAILS ORDER'
            CFLAG(37:72) = ' 3 CRITICAL VALUE TEST *****        '
          ENDIF
        ELSE
          CRIT = 0.D0
          CFLAG( 1:22) = '    ******** ERROR - '''
          CFLAG(23:24) = CC34
          CFLAG(25:72) = ''' IS NOT A VALID CRITICAL VALUE CODE *******'
        ENDIF
      ELSE
        LCRIT = .FALSE.
      ENDIF

*** PRINT RESULTS

      NLINE = 5
      IF (LCRIT) NLINE = NLINE + 1
      CALL LINE4 (NLINE)

      WRITE (LUNIT,1) ISSN, NAME1, JSSN, NAME2
    1 FORMAT ('0FROM: (', I5, ')  ', A30, '    TO: (', I5, ')  ', A30)
    
**v 4.28 *********************************************************
    
*     WRITE (LUNIT,2) DIST, SD, IREL1, SHIFT, IREL2
*   2 FORMAT (9X, 'DISTANCE=', 3X, F13.3, 2X, 'SIGMA=', F9.4, ' M.',
*    &        5X, 'ACC. 1:', I8, 5X, 'HOR. SHIFT=', F7.2, '   1:', I12)
     
      IF(IX1 .EQ. 0 .AND. IX2 .EQ. 0) THEN     
      WRITE (LUNIT,2) DIST, SD, IREL1, SHIFT, IREL2
    2 FORMAT (9X, 'DISTANCE=', 3X, F13.3, 2X, 'SIGMA=', F9.4, ' M.',
*    &        5X, 'ACC. 1:', I8, 5X, 'HOR. SHIFT=', F7.2, '   1:', I12)
     &        5X, 'ACC. 1:',I8,T80, 'HOR. SHIFT=',F7.2, '   1:', I12)
     
     
      ELSEIF(IX1 .EQ. 0 .AND. IX2 .EQ. 1) THEN 
      WRITE (LUNIT,300) DIST, SD, IREL1, SHIFT, REL2
 300  FORMAT (9X, 'DISTANCE=', 3X, F13.3, 2X, 'SIGMA=', F9.4, ' M.',
*    &        5X, 'ACC. 1:', I8, 5X, 'HOR. SHIFT=', F7.2, '   1:', 
*    &        E12.5) 
     &        5X, 'ACC. 1:', I8, T80, 'HOR. SHIFT=', F7.2, '   1:', 
     &        E12.5) 
     
     
      ELSEIF(IX1 .EQ. 1 .AND. IX2 .EQ. 0) THEN 
      WRITE (LUNIT,301) DIST, SD, REL1, SHIFT, IREL2
  301 FORMAT (9X, 'DISTANCE=', 3X, F13.3, 2X, 'SIGMA=', F9.4, ' M.',
*    &        5X, 'ACC. 1:', E12.5, 5X, 'HOR. SHIFT=', F7.2, '   1:',
*    &        I12)
     &        5X, 'ACC. 1:', E8.2, T80, 'HOR. SHIFT=', F7.2, '   1:',
     &        I12)
     
     
      ELSEIF(IX1 .EQ. 1 .AND. IX2 .EQ. 1) THEN 
      WRITE (LUNIT,302) DIST, SD, REL1, SHIFT, REL2
  302 FORMAT (9X, 'DISTANCE=', 3X, F13.3, 2X, 'SIGMA=', F9.4, ' M.',
*    &        5X, 'ACC. 1:', E12.5, 5X, 'HOR. SHIFT=', F7.2, '   1:',
*    &        E12.5)
     &        5X, 'ACC. 1:', E8.2, T80, 'HOR. SHIFT=', F7.2, '   1:',
     &        E12.5)
     
      ENDIF
      
     
**********************************************************************
*** v 4.29c **************
    
*     WRITE (LUNIT,2) DIST, SD, IREL1, SHIFT, IREL2
*   2 FORMAT (9X, 'DISTANCE=', 3X, F13.3, 2X, 'SIGMA=', F9.4, ' M.',
*    &        5X, 'ACC. 1:', I8, 5X, 'HOR. SHIFT=', F7.2, '   1:', I12)
     
*****************************
*v 4.29
*     IF (LCRIT) WRITE (LUNIT,5) CRIT, CFLAG
*   5 FORMAT (9X, 'CRITICAL VALUE=', 18X, F9.4, A72)
**********************************************************************
      WRITE (LUNIT,3) ID1, IM1, S1, SA, SADST
    3 FORMAT (9X, 'AZIMUTH=', 3X, I4, I3, F6.2, 3X, 'SIGMA=', F5.2,
     &        ' SEC./OR', F7.4, ' M.')
      WRITE (LUNIT,4) ID2, IM2, S2, A2, SV, SVDST
    4 FORMAT (9X, 'VERT.ANG.=', 2X, I4, I3, F6.2, A1, 2X, 'SIGMA=',
     &        F5.2, ' SEC./OR', F7.4, ' M.')

      RETURN
      END
      SUBROUTINE RESID (IUO,A,B,NX,G,SIGUWT,IUO3)

*** LIST ADJUSTED OBSERVATIONS AND RESIDUALS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10, MXSSN = 9999 )
      DIMENSION A(*),B(*),NX(*),G(*)
      DIMENSION IC(LENC),C(LENC)
      LOGICAL L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT

*** HEADING

      CALL HEAD2

*** INITIALIZE RESIDUAL STATISTICS

      CALL RSINIT
      ISNX = 0
      JSNX = 0
      KSNX = 0
      
**v 4.28***********
      IVC = 1
*******************

c Mike Potterfield 2/13/06
c Initialize IGPSSL to keep track of GPS solution numbers
      IGPSSL = 0

*** LOOP OVER THE OBSERVATION EQUATIONS

      REWIND IUO
  100 READ (IUO,END=777) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                   IVF,IAUX,IGRT
      IF ( KIND .LE. 999  .AND.
     &    (KIND .LT. 18  .OR.  KIND .GT. 20) ) THEN
        CALL ROBS (KIND,IOBS,ISN,ISNX,JSN,JSNX,KSNX,OBSB,SD,B,IAUX,
     &             IVF,SIGUWT,IC,C,LENG,A,NX,IGRT)
      ELSEIF (KIND .GE. 1000  .AND.  KIND .LT. 1999) THEN

c Mike Potterfield 2/13/06
c Increment the solution number
        IGPSSL = IGPSSL + 1

        NVEC = ISN
        IAUX = JSN
        IF ( .NOT. L2HLF) THEN
          LENG = (NVEC+1)*IDIM+1
        ELSE
          LENG = (NVEC+1)*  3 +1
        ENDIF
        IF (NGRT .GT. 0) LENG = LENG + 3
        NR = 3*NVEC
*       NC = NR+3+LENG
        NC = NR+5+LENG
	
**v 4.28*********************
        CALL RGPS (IUO,IOBS,NVEC,NR,NC,LENG,G,B,IAUX,IVF,SIGUWT,A,NX,
*    &             IGRT)
     &             IGRT,IVC,IUO3,IGPSSL)
******************************
     
      ELSEIF (KIND .GE. 18  .AND.  KIND .LE. 20) THEN
        CALL RDOP (IUO, KIND, IOBS, ISN, JSN, OBSB, SD, B, IAUX,
     &             IVF, SIGUWT, IC, C, LENG, A, NX, IGRT)
      ENDIF
      GO TO 100

*** END OF FILE ENCOUNTERED
*** LIST RESIDUAL STATISTICS

  777 CALL RSOUT (SIGUWT)
      RETURN
      END
      SUBROUTINE RESID2 (IUO, A, B, NX, G, SIGUWT)

*** COMPUTE THE RESIDUALS OF HORIZ OBS GROUPED AROUND INTERSECTION STA.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( NVECS = 700, MXSSN = 9999, MXFOT = 3500)
      PARAMETER ( LENC = 10, MXVF = 40 )
      DIMENSION ICM((NVECS+1)*3+1+3)
      DIMENSION KINDS(NVECS*3), ISNS(NVECS*3), JSNS(NVECS*3)
      DIMENSION LOBS(NVECS*3)
      DIMENSION COVECF(3,3)
      CHARACTER*30 NAMES, NAME1, NAME2, NAME3
      CHARACTER*1 ADIR1, ADIR2, ADIR3
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL PROP
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL L2HLF, LEHT
      DIMENSION B(*), A(*), NX(*), IC(LENC), C(LENC), G(*)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON /VFTB2/  VFS(MXVF), VTV(MXVF), VFRN(MXVF), VFSS(MXVF)
      COMMON /FOURTH/ NFOT, NFOTS(MXFOT)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

      SAVE GLAT, GLON, EHT, RX, RY
      SAVE NAME1, NAME2, NAME3
      SAVE ISN, ISNX, JSN, JSNX

      IF (NFOT .EQ. 0) RETURN

*** HEADING

      CALL HEAD7
      ISNX = 0
      JSNX = 0
      KSNX = 0

*** LOOP OVER THE OBSERVATION EQUATIONS

      DO 777 IFOT = 1, NFOT
        INTT = NFOTS(IFOT)
        REWIND IUO
 100    READ (IUO, END=777) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD,
     &                      IOBS, IVF, IAUX, IGRT
        IF (  ( ( KIND .GT.  0  .AND.  KIND .LT.  18 ) .OR.
     &          ( KIND .GT. 20  .AND.  KIND .LE. 999 ) )  .AND.
     &      INTT .EQ. JSN  ) THEN

*** CHECK FOR NEW ISN

          IF (ISNX .NE. ISN  .AND.  KIND .GT. 0  .AND.  KIND .LT. 7)THEN
            ISNX = ISN
            CALL LINE7 (1)
            CALL GETGLA (GLAT, ISN, B)
            CALL GETGLO (GLON, ISN, B)
            IF (L2HLF  .AND.  I2HLF(ISN) .EQ. 1) THEN
              CALL GETEHT (EHT, ISN, B)
            ELSE
              CALL GETMSL (GMSL, ISN, B)
              CALL GETGH (GHT, ISN, B)
              EHT = GMSL + GHT
            ENDIF
            CALL RADII (ISN, RMER, RPV, B)
            RX = RMER + EHT
            RY = (RPV + EHT)*DCOS(GLAT)
            NAME1 = NAMES(ISN)
          ENDIF
          IF (ISNX .NE. ISN  .AND.  KIND .GE. 7) THEN
            ISNX = ISN
            CALL LINE7 (1)
            NAME1 = NAMES(ISN)
          ENDIF
          IF (JSNX .NE. JSN  .AND.  KIND .GE. 7) THEN
            JSNX = JSN
            CALL LINE7 (1)
            WRITE (LUNIT,1)
    1       FORMAT (' ')
            NAME2 = NAMES(JSN)
          ENDIF
          IF (KIND .EQ. 10  .AND.  KSNX .NE. IAUX) THEN
            KSNX = IAUX
            CALL LINE7 (1)
            NAME3 = NAMES(IAUX)
          ENDIF

*** SCALE BY THE VARIANCE FACTOR

          IF (NVFTOT .GT. 0) THEN
            IF (IVF .LT. 0  .OR.  IVF .GT. NVFTOT) THEN
              WRITE (LUNIT,9) IVF, NVFTOT
    9         FORMAT ('0ILLEGAL IVF=',I10,' FOR N=',I10,' IN RESID2')
          write (lunit,*) 'subs2, 58'
              CALL ABORT2
            ELSEIF (IVF .NE. 0) THEN
              SD = SD*DSQRT( VFS(IVF) )
            ENDIF
          ENDIF

*** COMPUTATION OF S.D. FOR RESIDUALS AND
***  COMPUTATION OF REDUNDANCY NUMBERS

          IF (IMODE .EQ. 3) THEN
            IF ( PROP(C, IC, LENG, VAR, A, NX, IFLAG) ) THEN
              DVAR = SD*SD - VAR
              IF ( (KIND .GE. 8  .AND.  KIND .LE. 12)  .OR.
     &             KIND .EQ. 14 ) THEN
                IF (DVAR .LT. 1.0D-14) DVAR = 0.D0
              ELSE
                IF (DVAR .LT. 1.0D-10) DVAR = 0.D0
              ENDIF
              SD3 = DSQRT(DVAR)
              RN = 1.D0 - VAR/(SD*SD)
              IF (RN .LT. 1.0D-10) RN = 0.D0
            ELSEIF (IFLAG .EQ. 1) THEN
              WRITE (LUNIT,660)
 660          FORMAT ('0SYSTEM NOT INVERTED IN RESID2'/)
          write (lunit,*) 'subs2, 59'
              CALL ABORT2
            ELSE
              WRITE (LUNIT,661)
 661          FORMAT ('0ALL COVARIANCE ELEMENTS NOT WITHIN PROFILE--'
     &                ,'RESID2')
          write (lunit,*) 'subs2, 60'
              CALL ABORT2
            ENDIF
          ENDIF
          IF (IMODE .GT. 0  .AND.  IMODE .LE. 2) SD3 = SD

*** LATITUDE CONSTRAINT

          IF (KIND .EQ. 1) THEN
            CALL GETDMS (GLAT, ID1, IM1, S1, ISIGN)
            IF (ISIGN .GT. 0) THEN
              ADIR1 = 'N'
            ELSE
              ADIR1 = 'S'
            ENDIF
            V = GLAT - OBSB
            VSEC = V*RAD*3600.D0
            VMET = V*RX
            VSD = VMET/SD
            VSD3 = DIVID(VMET, SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL GETDMS (OBSB, ID3, IM3, S3, ISIGN)
            IF (ISIGN .GT. 0) THEN
              ADIR3 = 'N'
            ELSE
              ADIR3 = 'S'
            ENDIF
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,132) IOBS, ID1, IM1, S1, ADIR1, ID3, IM3,S3,
     &                      ADIR3, VSEC, VMET, VSD3, RN, NAME1
 132          FORMAT (' ',I6,' LA',I4,I3,F9.5,A1,I4,I3,F9.5,A1,F9.5,
     &                F8.3,F7.1,F4.1,1X,A30)
            ELSE
              WRITE (LUNIT,13) IOBS, ID1, IM1, S1, ADIR1, ID3, IM3, S3,
     &                     ADIR3, VSEC, VMET, VSD, NAME1
   13         FORMAT (' ',I6,' LA',I4,I3,F9.5,A1,I5,I3,F9.5,A1,
     &                F10.5,F8.3,F8.1,1X,A30)
            ENDIF

*** LONGITUDE CONSTRAINT

          ELSEIF (KIND .EQ. 2) THEN
            CALL GETDMS (GLON,ID2,IM2,S2,ISIGN)
            IF (ISIGN .GT. 0) THEN
              ALONGD=ID2+(IM2/60.D0)+(S2/3600.D0)
              ALONGD=360.D0-ALONGD
              ID2=INT(ALONGD)
              ALONGM=(ALONGD-ID2) * 60.D0
              IM2=INT(ALONGM)
              S2=(ALONGM-IM2) * 60.D0
              ADIR2 = 'W'
            ELSE
              ADIR2 = 'W'
            ENDIF
            V = GLON - OBSB
            VSEC = V*RAD*3600.D0
            VMET = V*RY
            VSD = VMET/SD
            VSD3 = DIVID(VMET,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL GETDMS (OBSB,ID3,IM3,S3,ISIGN)
            IF (ISIGN .GT. 0) THEN
              ALONGD=ID3+(IM3/60.D0)+(S3/3600.D0)
              ALONGD=360.D0-ALONGD
              ID3=INT(ALONGD)
              ALONGM=(ALONGD-ID3) * 60.D0
              IM3=INT(ALONGM)
              S3=(ALONGM-IM3) * 60.D0
              ADIR3 = 'W'
            ELSE
              ADIR3 = 'W'
            ENDIF
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,142) IOBS,ID2,IM2,S2,ADIR2,ID3,IM3,S3,ADIR3,
     &                      VSEC,VMET,VSD3,RN,NAME1
 142          FORMAT (' ',I6,' LO',I4,I3,F9.5,A1,I4,I3,F9.5,A1,F9.5,
     &                F8.3,F7.1,F4.1,1X,A30)
            ELSE
              WRITE (LUNIT,14) IOBS,ID2,IM2,S2,ADIR2,ID3,IM3,S3,ADIR3,
     &                     VSEC,VMET,VSD,NAME1
   14         FORMAT (' ',I6,' LO',I4,I3,F9.5,A1,I5,I3,F9.5,A1,
     &                F10.5,F8.3,F8.1,1X,A30)
            ENDIF

*** HEIGHT CONSTRAINT

          ELSEIF (KIND .EQ. 3) THEN
            V = EHT - OBSB
            VSD = V/SD
            VSD3 = DIVID(V,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,152) IOBS,EHT,OBSB,V,VSD3,RN,NAME1
 152          FORMAT (' ',I6,' EH',F14.3,F17.3,F20.3,F7.1,F4.1,1X,A30)
            ELSE
              WRITE (LUNIT,15) IOBS,EHT,OBSB,V,VSD,NAME1
   15         FORMAT (' ',I6,' EH',F14.3,F18.3,F21.3,F8.1,1X,A30)
            ENDIF

*** MARK TO MARK DISTANCE

          ELSEIF (KIND .EQ. 7) THEN
            CALL COMPOB (KIND,S,B,OBSB,ISN,JSN,IAUX,IGRT)
            V = S - OBSB
            VSD = V/SD
            VSD3 = DIVID(V,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,172) IOBS,S,OBSB,V,VSD3,RN,NAME1,NAME2
 172          FORMAT (' ',I6,' DIST',F13.4,F17.4,F19.3,F7.1,F4.1,
     &                2(1X,A30) )
            ELSE
              WRITE (LUNIT,17) IOBS,S,OBSB,V,VSD,NAME1,NAME2
  17          FORMAT (' ',I6,' DIST',F12.3,F18.3,F21.3,
     &                F8.1,2(1X,A30) )
            ENDIF

*** ASTRONOMIC AZIMUTH

          ELSEIF (KIND .EQ. 8) THEN
            CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
            CALL DIRDMS (OBS0,ID0,IM0,S0)
            V = OBS0 - OBSB
            IF (V .GT. PI) THEN
              V = V - PI - PI
            ELSEIF (V .LT. -PI) THEN
              V = V + PI + PI
            ENDIF
            VSEC = V*RAD*3600.D0
            VSD = V/SD
            VSD3 = DIVID(V,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL DIRDMS (OBSB,IDB,IMB,SB)
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,192) IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,VSD3,
     &                      RN,NAME1,NAME2
 192          FORMAT (' ',I6,' AZ',I4,I3,F6.2,'N',I7,I3,F6.2,'N',
     &                F9.2,F18.1,F4.1,2(1X,A30) )
            ELSE
              WRITE (LUNIT,19) IOBS,ID0,IM0,S0,IDB,IMB,SB,
     &                     VSEC,VSD,NAME1,NAME2
 19           FORMAT (' ',I6,' AZ',I4,I3,F6.2,'N',I8,I3,F5.1,
     &                'N',4X,F10.5,8X,F8.1,2(1X,A30) )
            ENDIF

*** ZENITH DISTANCE

          ELSEIF (KIND .EQ. 9) THEN
            CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
            CALL VERDMS (OBS0,ID0,IM0,S0,ISIGN)
            V = OBS0 - OBSB
            VSEC = V*RAD*3600.D0
            VSD = V/SD
            VSD3 = DIVID(V,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL VERDMS (OBSB,IDB,IMB,SB,ISIGN)
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,212) IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,VSD3,
     &                      RN,NAME1,NAME2
 212          FORMAT (' ',I6,' ZD',I4,I3,F6.2,' ',I7,I3,F6.2,
     &                F10.2,F18.1,F4.1,2(1X,A30) )
            ELSE
              WRITE (LUNIT,21) IOBS,ID0,IM0,S0,IDB,IMB,SB,
     &                     VSEC,VSD,NAME1,NAME2
  21          FORMAT (' ',I6,' ZD',I4,I3,F6.2,I9,I3,F6.2,
     &                4X,F10.5,8X,F8.1,2(1X,A30) )
            ENDIF

*** HORIZONTAL ANGLE

          ELSEIF (KIND .EQ. 10) THEN
            CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
            CALL DIRDMS (OBS0,ID0,IM0,S0)
            V = OBS0 - OBSB
            IF (V .GT. PI) THEN
              V = V - PI - PI
            ELSEIF (V .LT. -PI) THEN
              V = V + PI + PI
            ENDIF
            VSEC = V*RAD*3600.D0
            VSD = V/SD
            VSD3 = DIVID(V,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL DIRDMS (OBSB,IDB,IMB,SB)
            CALL LINE7 (2)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,232) IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,VSD3,
     &                      RN,NAME1,NAME2,NAME3
 232          FORMAT (' ',I6,' HA',I4,I3,F6.2,' ',I7,I3,F6.2,
     &                F10.2,F18.1,F4.1,2(1X,A30)/T103,A30)
            ELSE
              WRITE (LUNIT,23) IOBS,ID0,IM0,S0,IDB,IMB,SB,
     &                     VSEC,VSD,NAME1,NAME2,NAME3
 23           FORMAT (' ',I6,' HA',I4,I3,F6.2,' ',I8,I3,F6.2,
     &                4X,F10.5,8X,F8.1,2(1X,A30)/T102,A30)
            ENDIF

*** HORIZONTAL DIRECTION

          ELSEIF (KIND .EQ. 11) THEN
            CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
            CALL DIRDMS (OBS0,ID0,IM0,S0)
            V = OBS0 - OBSB
            IF (V .GT. PI) THEN
              V = V - PI - PI
            ELSEIF (V .LT. -PI) THEN
              V = V + PI + PI
            ENDIF
            VSEC = V*RAD*3600.D0
            VSD = V/SD
            VSD3 = DIVID(V,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL DIRDMS (OBSB,IDB,IMB,SB)
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,252) IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,VSD3,
     &                      RN,NAME1,NAME2
 252          FORMAT (' ',I6,' HD',I4,I3,F6.2,I8,I3,F6.2,
     &                F10.2,F18.1,F4.1,2(1X,A30) )
            ELSE
              WRITE (LUNIT,25) IOBS,ID0,IM0,S0,IDB,IMB,SB,
     &                     VSEC,VSD,NAME1,NAME2
 25           FORMAT (' ',I6,' HD',I4,I3,F6.2,' ',I8,I3,F6.2,
     &                4X,F10.5,8X,F8.1,2(1X,A30) )
            ENDIF

*** CONSTRAINED ASTRONOMIC AZIMUTH

          ELSEIF (KIND .EQ. 12) THEN
            I8 = 8
            CALL COMPOB (I8,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
            CALL DIRDMS (OBS0,ID0,IM0,S0)
            V = OBS0 - OBSB
            IF (V .GT. PI) THEN
              V = V - PI - PI
            ELSEIF (V .LT. -PI) THEN
              V = V + PI + PI
            ENDIF
            VSEC = V*RAD*3600.D0
            VSD = V/SD
            VSD3 = DIVID(V,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL DIRDMS (OBSB,IDB,IMB,SB)
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,272) IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,VSD3,
     &                      RN,NAME1,NAME2
 272          FORMAT (' ',I6,' CA',I4,I3,F6.2,'N',I7,I3,F6.2,
     &                'N',2X,F10.5,7X,F8.1,F4.1,2(1X,A30) )
            ELSE
              WRITE (LUNIT,27) IOBS,ID0,IM0,S0,IDB,IMB,SB,
     &                     VSEC,VSD,NAME1,NAME2
 27           FORMAT (' ',I6,' CA',I4,I3,F6.2,'N',I8,I3,F6.2,
     &                'N',3X,F10.5,8X,F8.1,2(1X,A30) )
            ENDIF

*** CONSTRAINED MARK TO MARK DISTANCE

          ELSEIF (KIND .EQ. 13) THEN
            I7 = 7
            CALL COMPOB (I7,S,B,OBSB,ISN,JSN,IAUX,IGRT)
            V = S - OBSB
            VSD = V/SD
            VSD3 = DIVID(V,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,292) IOBS,S,OBSB,V,VSD3,RN,NAME1,NAME2
 292          FORMAT (' ',I6,' CD',F15.4,F17.4,F19.3,F7.1,F4.1,
     &                 2(1X,A30) )
            ELSE
              WRITE (LUNIT,29) IOBS,S,OBSB,V,VSD,NAME1,NAME2
  29          FORMAT (' ',I6,' CD ',F14.4,F18.4,F20.3,
     &                F8.1,2(1X,A30) )
            ENDIF

*** CONSTRAINED ZENITH DISTANCE

          ELSEIF (KIND .EQ. 14) THEN
            I9 = 9
            CALL COMPOB (I9,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
            CALL VERDMS (OBS0,ID0,IM0,S0,ISIGN)
            V = OBS0 - OBSB
            VSEC = V*RAD*3600.D0
            VSD = V/SD
            VSD3 = DIVID(V,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL VERDMS (OBSB,IDB,IMB,SB,ISIGN)
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,312) IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,VSD3,
     &                      RN,NAME1,NAME2
 312          FORMAT (' ',I6,' CZ',I4,I3,F6.2,' ',I7,I3,F6.2,
     &                3X,F10.5,7X,F8.1,F4.1,2(1X,A30) )
            ELSE
              WRITE (LUNIT,31) IOBS,ID0,IM0,S0,IDB,IMB,SB,
     &                     VSEC,VSD,NAME1,NAME2
  31          FORMAT (' ',I6,' CZ',I4,I3,F6.2,I9,I3,F6.2,
     &                4X,F10.5,8X,F8.1,2(1X,A30) )
            ENDIF

*** CONSTRAINED GEOID HEIGHT DIFFERENCE

          ELSEIF (KIND .EQ. 15) THEN
            CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
            V = OBS0 - OBSB
            VSD = V/SD
            VSD3 = DIVID(V,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,412) IOBS,OBS0,OBSB,V,VSD3,RN,NAME1,NAME2
 412          FORMAT (' ',I6,' DN',F15.4,F17.4,F19.3,F7.1,F4.1,
     &                 2(1X,A30) )
            ELSE
              WRITE (LUNIT,41) IOBS,OBS0,OBSB,V,VSD,NAME1,NAME2
 41           FORMAT (' ',I6,' DN',F14.3,F18.3,F21.3,F8.1,
     &                2(1X,A30) )
            ENDIF

*** CONSTRAINED ORTHOMETRIC HEIGHT DIFFERENCE

          ELSEIF (KIND .EQ. 16) THEN
            CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
            V = OBS0 - OBSB
            VSD = V/SD
            VSD3 = DIVID(V,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,452) IOBS,OBS0,OBSB,V,VSD3,RN,NAME1,NAME2
 452          FORMAT (' ',I6,' DO',F15.4,F17.4,F19.3,F7.1,F4.1,
     &                 2(1X,A30) )
            ELSE
              WRITE (LUNIT,45) IOBS,OBS0,OBSB,V,VSD,NAME1,NAME2
 45           FORMAT (' ',I6,' DO',F14.3,F18.3,F21.3,F8.1,
     &                2(1X,A30) )
            ENDIF

*** CONSTRAINED ELLIPSOIDAL HEIGHT DIFFERENCE

          ELSEIF (KIND .EQ. 17) THEN
            CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
            V = OBS0 - OBSB
            VSD = V/SD
            VSD3 = DIVID(V,SD3)
            IF (LSS) VSD = VSD/SIGUWT
            CALL LINE7 (1)
            IF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,492) IOBS,OBS0,OBSB,V,VSD3,RN,NAME1,NAME2
 492          FORMAT (' ',I6,' DE',F15.4,F17.4,F19.3,F7.1,F4.1,
     &                 2(1X,A30) )
            ELSE
              WRITE (LUNIT,49) IOBS,OBS0,OBSB,V,VSD,NAME1,NAME2
 49           FORMAT (' ',I6,' DE',F14.3,F18.3,F21.3,F8.1,
     &                2(1X,A30) )
            ENDIF

*** ILLEGAL KIND

          ELSE
            WRITE (LUNIT,667) KIND
  667       FORMAT ('0ILLEGAL KIND IN RESID2--',I5)
          write (lunit,*) 'subs2, 61'
            CALL ABORT2
          ENDIF

        ELSEIF (KIND .EQ. 18  .OR.  KIND .EQ. 19  .OR.
     &          KIND .EQ. 20) THEN
          READ (IUO,END=777) KIND, ISN, JSN, IC, C, LENG, CMO, OBSBY,
     &                       SIGY, IOBSY, IVF, IAUX, IGRT
          READ (IUO,END=777) KIND, ISN, JSN, IC, C, LENG, CMO, OBSBZ,
     &                       SIGZ, IOBSZ, IVF, IAUX, IGRT
          READ (IUO,END=777) ISIGNL, ( (COVECF(I,J), J = 1,3), I = 1,3 )

        ELSEIF (KIND .GE. 1000) THEN
          NVEC = ISN
          IAUX = JSN
          IF ( .NOT. L2HLF) THEN
            LENG = (NVEC+1)*IDIM + 1
          ELSE
            LENG = (NVEC+1)*3 + 1
          ENDIF
          IF (NGRT .GT. 0) LENG = LENG + 3
          NR = 3*NVEC
*         NC = NR + 3 + LENG
          NC = NR + 5 + LENG
	  
          READ (IUO,END=666) ICM,NICM,KINDS,ISNS,JSNS,LOBS
          CALL DUMRD (IUO,NR,NC,G)

        ENDIF
        GO TO 100
 777  CONTINUE

      RETURN

*** PREMATURE END OF FILE

  666 WRITE (LUNIT,669) NVEC
  669 FORMAT ('0PREMATURE FILE END IN RESID2 -- NVEC=',I5)
          write (lunit,*) 'subs2, 62'
      CALL ABORT2
      RETURN
      END


      SUBROUTINE RGPS (IUO, IOBS, NVEC, NR, NC, LENG, G, B, IAUX, IVF,
**v 4.28*******************
*    &                 SIGUWT, A, NX, IGRT)
     &                 SIGUWT, A, NX, IGRT,IVC, IUO3,IGPSSL)
*******************************

c Mike Potterfield 2/13/06
c The IUO3 parameter is the scratch file with rejects and project ids
c The IGPSSL parameter is the running count of the current solution

*** LIST ADJUSTED OBSERVATIONS AND RESIDUALS

* * * MODIFIED 1/2/91 TO PRINT ALL RESIDUALS REGARDLESS OF THE SIZE
* 11-15-01  print all residuals for a vector whenever residual for one
* component exceeds tolerance in afile (PP rec cc 6-8 = XX.x mm)

** v 4.29 
* 7-16-02   print all residuals for a vector whenever residual for one
* component exceeds tolerance in afile (PP rec cc 24-26 = XX.x mm)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
c Mike Potterfield 2/13/06
c The next two variables will identify the project
c and rejected vectors in the adjusted obs output
      CHARACTER*1 CREJCT
      CHARACTER*13 TPRJID
      PARAMETER ( NVECS = 700, MXSSN = 9999 )

c Mike Potterfield 3/22/06, modified 3/15/07
c New variable IGPSOS used for component count
      SAVE IGPSOS
      
**v 4.28
      PARAMETER ( MAXC = 99999 )
      
      DIMENSION ICM((NVECS+1)*3+1+3)
      DIMENSION KINDS(NVECS*3), ISNS(NVECS*3), JSNS(NVECS*3)
      DIMENSION LOBS(NVECS*3)
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL FATAL, NOBIGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      LOGICAL LNEWX


c  Mike Potterfield new arrays below added to compute the 
c  std devs of residuals in N, E, U 1/28/06 - 2/2/06
c  first an array to recreate the XYZ cov mat
      DIMENSION XYZCOV(3,3)
c next an array for cov of adjusted XYZ
      DIMENSION XYZADJ(3,3)
c next an array for cov of XYZ resids
      DIMENSION XYZRSD(3,3)
c next an array to hold the covmat of NEU residuals
      DIMENSION TNEUCV(3,3)
c next the array for rotations
      DIMENSION RM(3,3)
c and finally a work vector for ABAT
      DIMENSION WRKSTD(3)


**v 4.28h
      logical lc,lp
********
** v 4.29f
      logical locssn
***************
      CHARACTER*30 NAMES, NAME1, NAME2
      CHARACTER*30  NAME1X, NAME2X
      CHARACTER*30  NAME1Y, NAME2Y
      CHARACTER*30  NAME1Z, NAME2Z
      
**v 4.28
      CHARACTER*5 SESSID,SESS
      COMMON /SESTAB/ SESS(MAXC)
***************
      DIMENSION G(NR,NC), B(*), A(*), NX(*)
      DIMENSION GMDES(NVECS*3)
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
* v 4.29     
      COMMON /OPRIN2/ CRIT2
     
      COMMON /GRTTB1/ GLAT0, GLON0, R0(3,3)
      COMMON/UNITS/ LUNIT

c Mike Potterfield 2/2/06 add the common block below to access NUNK
c in accessing the covariance matrix of adjusted XYZ components for
c the purpose of computing std devs of residuals in dNEU
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT

      SAVE ISNX, JSNX, ISN, JSN
      SAVE X1, Y1, Z1, X2, Y2, Z2, RNSUM
      SAVE NAME1, NAME2
      
* v 4.28h
      save lp
** v 4.29f
      save isssn,jsssn
**************
** v 4.29f
      save iobsl
**********
      
      DATA ISNX /0/
      DATA JSNX /0/

*** DO NOT ABORT/PRINT BIG RESIDUALS

      NOBIGV = .TRUE.
      
**v 4.28h print all components if one component exceeds tolerance
* lc = true,  a component residual exceeds tolerance crit2, print all components 
* lp = true, to print blank line before first DX
*      components.

      lc = .false.
      
* v 4.29      
*     crit2 = crit * 0.001
      crit3 = crit2 * 0.001
      
      READ (IUO,END=666) ICM, NICM, KINDS, ISNS, JSNS, LOBS
      
c Mike Potterfield 2/13/06
c Read the Project ID from scratch file IUO3
      READ(IUO3) TPRJID

*** LOAD WORK SPACE (G)

      DO 1 I = 1, NR
        READ (IUO,END=666) (G(I, J), J = 1, NC)
    1 CONTINUE

*** COMPUTE AND DECORRELATE OBS. EQ. AND MISCLOSURE
*** IF MODE 3, COMPUTE STD. DEV. OF RESIDUALS AND REDUNDANCY NUMBERS
*** IF MODE 0, COMPUTE MARGINALLY DETECT. ERROR AND REDUNDANCY NUMBERS
*** IF MODE 3, COMPUTE MDE, PATCH WORK

      IF (IMODE .EQ. 3  .OR.  IMODE .EQ. 0) THEN
        CALL GPSSDV (G, NR, NC, NVEC, LENG, B, ICM, NICM, KINDS, ISNS,
     &               JSNS, LOBS, IAUX, IVF, FATAL, NOBIGV, A, NX, N2,
     &               N3, N4, GMDES, IGRT)
      ELSE
        CALL GOBSEQ (G, NR, NC, NVEC, LENG, B, ICM, NICM, KINDS, ISNS,
     &               JSNS, LOBS, IAUX, IVF, FATAL, NOBIGV, N4, IGRT)
      ENDIF

*** UNLOAD WORK SPACE

      RNSUM = 0.D0
      DO 500 I = 1, NR
        OBSB = G(I,NC)
        IF (IMODE .EQ. 3) THEN
          VSD = G(I,N2)
          SD3 = G(I,N4)
** v 4.29f	  
*         IF (SD3 .EQ. 0.D0) SD3 = 1.0D + 20
          IF (SD3 .EQ. 0.D0) SD3 = 1.D0 + 20
******	  
	  
          RN = G(I,N3)
        ELSE
          VSD = G(I,N4)
          VSD3 = G(I,N4)
          RN = 0.D0
        ENDIF

        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) THEN
          RN = G(I,N3)
          IF (RN .EQ. 0.D0) THEN
            GMDE = 0.D0
          ELSE
            GMDE = 3.D0*DSQRT( GMDES(I) )
          ENDIF
        ENDIF

        LNEWX = .FALSE.
        KIND = KINDS(I)
        ISN = ISNS(I)
        JSN = JSNS(I)
        IOBS = LOBS(I)
        IF (ISN .NE. ISNX) THEN
          ISNX = ISN
          LNEWX = .TRUE.
          IF (L2HLF) THEN
            CALL GETEHX (X1, ISN, B)
            CALL GETEHY (Y1, ISN, B)
            CALL GETEHZ (Z1, ISN, B)
          ELSE
            CALL GETECX (X1, ISN, B)
            CALL GETECY (Y1, ISN, B)
            CALL GETECZ (Z1, ISN, B)
          ENDIF
          NAME1 = NAMES(ISN)
** v 4.29f
          if(.not.locssn(isn,isssn)) then
	     write(lunit,966) isn
 966         format('0SSN TABLE ERROR IN RGPS-- ',I8)
             write (lunit,*) 'subs2, 1'
             call abort2
	  endif
*************	     	     
	  	  
        ENDIF
        IF (JSN .NE. JSNX) THEN
          JSNX = JSN
          LNEWX = .TRUE.
          IF (L2HLF) THEN
            CALL GETEHX (X2, JSN, B)
            CALL GETEHY (Y2, JSN, B)
            CALL GETEHZ (Z2, JSN, B)
          ELSE
            CALL GETECX (X2, JSN, B)
            CALL GETECY (Y2, JSN, B)
            CALL GETECZ (Z2, JSN, B)
          ENDIF
          NAME2 = NAMES(JSN)
** v 4.29f
          if(.not.locssn(jsn,jsssn)) then
	     write(lunit,966) jsn
             write (lunit,*) 'subs2, 2'
             call abort2
	  endif
*************	     	     
	  
        ENDIF
        RNSUM = RNSUM + RN

*** GPS DEL X

        IF (KIND .EQ. 4) THEN
          OBS0X = X2 - X1
          IF (IAUX .GT. 0) THEN
            CALL GETAUX (VAL, IAUX, B)
            OBS0X = OBS0X - OBS0X*VAL
          ENDIF
          IF (IGRT .GT. 0) THEN
            DELY = Y2 - Y1
            IF (IAUX .GT. 0) DELY = DELY - DELY * VAL
            DELZ = Z2 - Z1
            IF (IAUX .GT. 0) DELZ = DELZ - DELZ * VAL
            CALL GETRTG (VALX, VALY, VALZ, IGRT)
            OBS0X = OBS0X
     &             + ( R0(1,3)*VALX + R0(2,3)*VALY + R0(3,3)*VALZ )*DELY
     &             - ( R0(1,2)*VALX + R0(2,2)*VALY + R0(3,2)*VALZ )*DELZ
          ENDIF
          IF (IMODE .EQ. 0) THEN
            OBSBX = OBS0X
          ELSE
            OBSBX = OBSB
          ENDIF

c  Initialize IGPSOS here so that the summaries have the correct
c  observation numbers.
c  The following algorithm gives the obsno of the first vector
c  component - 1, where NCON is the number of constraints and
c  IVC is the vector number.
c  Mike Potterfield 3/15/07

          IGPSOS = NCON + ((IVC-1)*3)

          VX = OBS0X - OBSBX
          IF (IMODE .EQ. 3) VSD3 = DIVID(VX, SD3)
          IF (SD3 .EQ. 1.0D+20) SD3 = 0.D0
c New additional argument to RSTAT Mike Potterfield 3/15/07
          CALL RSTAT (VX, VSD, RN, VSD3, IOBS, KIND,1+IGPSOS)
          IF (LSS) VSD = VSD/SIGUWT
	  
*** v 4.28n
          IF (LSS) SD3 = SD3/SIGUWT
          IF (LSS) VSD3 = VSD3/SIGUWT
**************
	  	  
          IF (IMODE .EQ. 3) THEN
            VSDX = SD3
          ELSE
            VSDX = VSD
          ENDIF
*v 4.28h
	      if(dabs(vx) .ge. crit3) lc = .true.  
              IF (IMODE .EQ. 0) THEN
*v 4.28h	      
                 iobsx = iobs
		 gmdex = gmde
		 name1x = name1
		 name2x = name2
              ELSEIF (IMODE .EQ. 3) THEN
	         iobsx = iobs
		 sd3x = sd3
		 vsd3x = vsd3
		 gmdex = gmde
		 name1x = name1
              ELSE
	         iobsx = iobs
		 vsdx = vsd
		 name1x = name1
		 name2x = name2 
              ENDIF

*** DEL Y GPS

        ELSEIF (KIND .EQ. 5) THEN
          OBS0Y = Y2 - Y1
          IF (IAUX .GT. 0) THEN
            CALL GETAUX (VAL, IAUX, B)
            OBS0Y = OBS0Y - OBS0Y*VAL
          ENDIF
          IF (IGRT .GT. 0) THEN
            DELX = X2 - X1
            IF (IAUX .GT. 0) DELX = DELX - DELX * VAL
            DELZ = Z2 - Z1
            IF (IAUX .GT. 0) DELZ = DELZ - DELZ * VAL
            CALL GETRTG (VALX, VALY, VALZ, IGRT)
            OBS0Y = OBS0Y
     &             - ( R0(1,3)*VALX + R0(2,3)*VALY + R0(3,3)*VALZ )*DELX
     &             + ( R0(1,1)*VALX + R0(2,1)*VALY + R0(3,1)*VALZ )*DELZ
          ENDIF
          IF (IMODE .EQ. 0) THEN
            OBSBY = OBS0Y
          ELSE
            OBSBY = OBSB
          ENDIF
          VY = OBS0Y - OBSBY
          IF (IMODE .EQ. 3) VSD3 = DIVID(VY, SD3)
          IF (SD3 .EQ. 1.0D+20) SD3 = 0.D0
c New additional argument to RSTAT Mike Potterfield 3/15/07
          CALL RSTAT (VY, VSD, RN, VSD3, IOBS, KIND,2+IGPSOS)
          IF (LSS) VSD = VSD/SIGUWT
	  
*** v 4.28n
          IF (LSS) SD3 = SD3/SIGUWT
          IF (LSS) VSD3 = VSD3/SIGUWT
**************
	  
          IF (IMODE .EQ. 3) THEN
            VSDY = SD3
          ELSE
            VSDY = VSD
          ENDIF
            if(dabs(vy) .ge. crit3) lc = .true.
              IF (IMODE .EQ. 0) THEN
	        iobsy = iobs
		gmdey = gmde
              ELSEIF (IMODE .EQ. 3) THEN
	        iobsy = iobs
		sd3y = sd3
		vsd3y = vsd3
		gmdey = gmde
		name2y = name2
              ELSE
                IF (LNEWX) THEN
		  iobsy = iobs
		  vsdy = vsd
		  name1y = name1
		  name2y = name2
                ELSE
		  iobsy = iobs
		  vsdy = vsd
                ENDIF
              ENDIF

*** DEL Z GPS

        ELSEIF (KIND .EQ. 6) THEN

**v 4.28*****************************
          SESSID = SESS(IVC)
	  IVC = IVC + 1
***************************************
	  	
          OBS0Z = Z2 - Z1
          IF (IAUX .GT. 0) THEN
            CALL GETAUX (VAL, IAUX, B)
            OBS0Z = OBS0Z - OBS0Z*VAL
          ENDIF
          IF (IGRT .GT. 0) THEN
            DELX = X2 - X1
            IF (IAUX .GT. 0) DELX = DELX - DELX * VAL
            DELY = Y2 - Y1
            IF (IAUX .GT. 0) DELY = DELY - DELY * VAL
              CALL GETRTG (VALX, VALY, VALZ, IGRT)
            OBS0Z = OBS0Z
     &             + ( R0(1,2)*VALX + R0(2,2)*VALY + R0(3,2)*VALZ )*DELX
     &             - ( R0(1,1)*VALX + R0(2,1)*VALY + R0(3,1)*VALZ )*DELY
          ENDIF
          IF (IMODE .EQ. 0) THEN
            OBSBZ = OBS0Z
          ELSE
            OBSBZ = OBSB
          ENDIF
          VZ = OBS0Z - OBSBZ
          IF (IMODE .EQ. 3) VSD3 = DIVID(VZ, SD3)
          IF (SD3 .EQ. 1.0D+20) SD3 = 0.D0
c New additional argument to RSTAT Mike Potterfield 3/15/07
          CALL RSTAT (VZ, VSD, RN, VSD3, IOBS, KIND,3+IGPSOS)
          IF (I .EQ. NR) THEN
            IF (IMODE .EQ. 3  .OR.  IMODE .EQ. 0) THEN
              CALL RSTAT2 (RNSUM)
            ENDIF
          ENDIF
          IF (LSS) VSD = VSD/SIGUWT
	  
	  
*** v 4.28n
          IF (LSS) SD3 = SD3/SIGUWT
          IF (LSS) VSD3 = VSD3/SIGUWT
**************
	  
          IF (IMODE .EQ. 3) THEN
            VSDZ = SD3
          ELSE
            VSDZ = VSD
          ENDIF
              if(dabs(vz) .ge. crit3) lc = .true.
              IF (I .EQ. NR) THEN
                IF (IMODE .EQ. 0) THEN
                   iobsz = iobs
		   gmdez = gmde
                ELSEIF (IMODE .EQ. 3) THEN
		   iobsz = iobs
		   sd3z = sd3
		   vsd3z = vsd3
		   gmdez = gmde
                ELSE
                  IF (LNEWX) THEN
		     iobsz = iobs
		     vsdz = vsd
		     name1z = name1
		     name2z = name2
                  ELSE
		     iobsz = iobs
		     vsdz = vsd
                  ENDIF
                ENDIF
              ELSE
                IF (IMODE .EQ. 0) THEN
		  iobsz = iobs
		  gmdez = gmde
                ELSEIF (IMODE .EQ. 3) THEN
		  iobsz = iobs
		  sd3z = sd3
		  vsd3z = vsd3
		  gmdez = gmde
                ELSE
                  IF (LNEWX) THEN
		     iobsz = iobs
		     vsdz = vsd
		     name1z = name1
		     name2z = name2
                  ELSE
		     iobsz = iobs
		     vsdz = vsd
                  ENDIF
                ENDIF
              ENDIF

***  N,E,U RESIDUALS AFTER Z COMPONENT

                CALL TOLGH (OBS0X, OBS0Y, OBS0Z, OBS0N, OBS0E, OBS0U,
     &                      ISN, JSN, B)
                CALL TOLGH (OBSBX, OBSBY, OBSBZ, OBSBN, OBSBE, OBSBU,
     &                      ISN, JSN, B)
                CALL TOLGH (VX, VY, VZ, VN, VE, VU, ISN, JSN, B)
***v 4.27
c Mike Potterfield 2/2/06 add XYZCOV and RM to COVLGH parameters
                CALL COVLGH (SN,SE,SU,ISN,JSN,B,I,G,NR,NC,LENG,
     &               XYZCOV,RM)
****           		

c Mike Potterfield 2/13/06
c Now read the reject code
      READ(IUO3) CREJCT

c Now compute the vector number
c I is always the z (third) component, so
      IVCNUM = I/3

c Mike Potterfield 2/2/06
c Now load in the covariance matrix for the adjusted observation
c I is the Z component for the current vector
c The A work array from GPSSDV only stores the upper half 
      NB1 = NX(NUNK + 1)
      DO 119 J = 1,3
         XYZADJ(J,J) = A(INX(J+I-3,NR)+NB1)
         XYZRSD(J,J) = XYZCOV(J,J) - XYZADJ(J,J)
         DO 118 K = J,3
            IF (J .NE. K) THEN
                XYZADJ(J,K) = A(NB1+(INX(J+I-3,NR)-(J+I-3))+K+I-3)
                XYZADJ(K,J) = XYZADJ(J,K)
                XYZRSD(J,K) = XYZCOV(J,K) - XYZADJ(J,K)
                XYZRSD(K,J) = XYZCOV(K,J) - XYZADJ(K,J)
            ENDIF
  118    CONTINUE
  119 CONTINUE

c Mike Potterfield 2/2/06
c      WRITE(LUNIT,*)'Now the covariance matrix of XYZ adjusted:'
c      WRITE(LUNIT,116) XYZADJ(1,1), XYZADJ(1,2), XYZADJ(1,3)
c      WRITE(LUNIT,116) XYZADJ(2,1), XYZADJ(2,2), XYZADJ(2,3)
c      WRITE(LUNIT,116) XYZADJ(3,1), XYZADJ(3,2), XYZADJ(3,3)
c      WRITE(LUNIT,*)'Now the covariance matrix of XYZ residuals:'
c      WRITE(LUNIT,116) XYZRSD(1,1), XYZRSD(1,2), XYZRSD(1,3)
c      WRITE(LUNIT,116) XYZRSD(2,1), XYZRSD(2,2), XYZRSD(2,3)
c      WRITE(LUNIT,116) XYZRSD(3,1), XYZRSD(3,2), XYZRSD(3,3)

c Now rotate XYZRSD into TNEUCV
      CALL ABAT(RM,XYZRSD,TNEUCV,WRKSTD,3,3)

c      WRITE(LUNIT,*)'Now the covariance matrix of NEU residuals:'
c      WRITE(LUNIT,116) TNEUCV(1,1), TNEUCV(1,2), TNEUCV(1,3)
c      WRITE(LUNIT,116) TNEUCV(2,1), TNEUCV(2,2), TNEUCV(2,3)
c      WRITE(LUNIT,116) TNEUCV(3,1), TNEUCV(3,2), TNEUCV(3,3)
		
                CALL TOLGH (VSDX, VSDY, VSDZ, VSDN, VSDE, VSDU,
     &                      ISN, JSN, B)
     
* v 4.28h
*               OBS0H = DSQRT(OBS0N**2 + OBS0E**2)
*               OBSBH = DSQRT(OBSBN**2 + OBSBE**2)
***********************************		
** v 4.29g		
**               VL = DSQRT(VN**2 + VE**2)
**              SH = DSQRT(SN**2 + SE**2)
*               VL = OBS0H - OBSBH       
*               SH = DSQRT(SN**2 + SE**2)
*************************************		

** v 4.29k		
                OBS0L = DSQRT(OBS0N**2 + OBS0E**2)
                OBSBL = DSQRT(OBSBN**2 + OBSBE**2)
                VL = DSQRT(VN**2 + VE**2)
                SL = DSQRT(SN**2 + SE**2)
*************************************		

	
	        if(vn .ge. crit3) lc = .true.	
	        if(ve .ge. crit3) lc = .true.	
	        if(vu .ge. crit3) lc = .true.	
	        if(vl .ge. crit3) lc = .true.	
		
		   	
*** ACCULULATE N,E,U STATS, BYPASS DOWNWEIGHTED OUTLIERS AND NO CHECK

              IF ( DABS(VX+VY+VZ) .LE. 1.D-6 ) THEN
                CONTINUE
              ELSEIF ( IMODE .EQ. 3              .AND.
     &                 VSDZ .GT. 2.D0) THEN
                CONTINUE
              ELSEIF ( ( IMODE .EQ. 1  .OR.  IMODE .EQ. 2 )  .AND.
     &                   DABS(VSDZ) .LE. 0.02) THEN
                CONTINUE
              ELSE
                CALL RSTAT4 (VN, VE, VU)

**v 4.28
** v 4.29h
c Fix the observation component numbers Mike Potterfield 10/26/08
c The original code used IOBS as the observation component number,
c but this has been modified so that rejected observation
c components are numbered sequentially.  The variable IGPSOS
c corrects the observation component number to match the rest
c of the Adjustment output.
c	           CALL RSTAT5(90,VN,IOBS)
	           CALL RSTAT5(90,VN,1+IGPSOS)
c	           CALL RSTAT5(91,VE,IOBS)
	           CALL RSTAT5(91,VE,2+IGPSOS)
c	           CALL RSTAT5(92,VU,IOBS)
	           CALL RSTAT5(92,VU,3+IGPSOS)
************		   
**v 4.29f  
c                   CALL RSTAT5(93,VL,IOBS)
                   CALL RSTAT5(93,VL,3+IGPSOS)
c End of changes by Mike Potterfield 10/26/08
*******************		   

** v 4.29f
c New additional argument to acumdu Mike Potterfield 1/15/07
                if(iobs.ne.iobsl) then
                call acum2du(iobs,vu,3+IGPSOS)
		endif
		
c New additional argument to acumdl Mike Potterfield 1/15/07
** v 4.29i
                if(iobs.ne.iobsl) then
                call acum2dl(iobs,vl,3+IGPSOS)
		endif
*********		
		
              ENDIF

* v 4.28h print all components of a vector if one exceeds tolerance

*** GPS DEL X

        if(.not.lp) then
         CALL LINE2 (1)
         WRITE (LUNIT,2) 
    2   FORMAT (' ')
        lp = .true.
        endif
	

      if(lc) then
      
        if(lvg) then
	  if(lgv) then
             call line2 (8)
          else  
             call line2 (4)         
          endif
         endif	  

	 if(lvg) then
              IF (IMODE .EQ. 0) THEN
** v 4.29f	      
*               WRITE (LUNIT,161) IOBSX, OBS0X, GMDEX,NAME1X,NAME2X
* 161         FORMAT (' ',I6,' DX',F15.4,F12.4,18X,
*    &                  2(1X,A30) )
                WRITE (LUNIT,161) IOBSX, OBS0X, GMDEX,isssn,NAME1X,
     &                jsssn,name2x		
  161         FORMAT (' ',I6,' DX',F15.4,F12.4,17X,
** v 4.29i
*    &                  2(2x,i5,2x,A30) )
     &                  2(1x,'(',i5,')',1x,A30) )
**************
     
              ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	      
*               WRITE (LUNIT,162) IOBSX,OBS0X,OBSBX,VX,SD3X,VSD3X,
*    &                       GMDEX,NAME1
*  162           FORMAT (' ',I6,' DX',F15.4,F17.4,F20.4,F9.4,F10.4,
*     &                  F10.4,8X,A30)
c Mike Potterfield 2/14/06
c The component numbering is revised to give unique numbers to all
c components, including rejected vectors
c                WRITE (LUNIT,162) IOBSX,OBS0X,OBSBX,VX,SD3X,VSD3X,

c Mike Potterfield 3/22/06 compute the offset here
c It is now being computed above 3/15/07
      WRITE (LUNIT,162)1+IGPSOS,OBS0X,OBSBX,VX,SD3X,
     &      VSD3X,GMDEX,isssn,NAME1
  162           FORMAT (' ',I6,' DX',F15.4,F17.4,F20.4,F9.4,F10.4,
** v 4.29i
*    &                  F10.4,t99,i5,2x,A30)
     &                  F10.4,t99,'(',i5,')',1x,A30)
******************      
              ELSE
** v 4.29f	      
*                WRITE (LUNIT,16) IOBSX,OBS0X,OBSBX,VX,VSDX,NAME1X,
*    &                  NAME2X	      
*  16           FORMAT (' ',I6,' DX',F15.4,F18.4,F21.4,F8.1,2(1X,A30) )
   
                 WRITE (LUNIT,16) IOBSX,OBS0X,OBSBX,VX,VSDX,isssn,
     &                  name1x,jsssn,NAME2X	      
   16           FORMAT (' ',I6,' DX',F15.4,F18.4,F21.4,F7.1,1x,
** v 4.29i
*    &                 2(i5,2x,a30))
     &                 2('(',i5,')',1x,a30))
***********************
   
              ENDIF

*** DEL Y GPS

              IF (IMODE .EQ. 0) THEN
                WRITE (LUNIT,171) IOBSY, OBS0Y, GMDEY
  171           FORMAT (' ',I6,' DY',F15.4,F12.4)
              ELSEIF (IMODE .EQ. 3) THEN
	      
** v 4.29f	      
*               WRITE (LUNIT,172) IOBSY,OBS0Y,OBSBY,VY,SD3Y,VSD3Y,
*    &                        GMDEY,NAME2Y
* 172           FORMAT (' ',I6,' DY',F15.4,F17.4,F20.4,F9.4,F10.4,
*    &                    F10.4,8X,A30)
c Mike Potterfield 2/14/06
c The component numbering is revised to give unique numbers to all
c components, including rejected vectors
c               WRITE (LUNIT,172) IOBSY,OBS0Y,OBSBY,VY,SD3Y,VSD3Y,
c Mike Potterfield component number updated by IGPSOS
      WRITE (LUNIT,172)2+IGPSOS,OBS0Y,OBSBY,VY,SD3Y,
     &        VSD3Y,GMDEY,jsssn,NAME2Y
  172           FORMAT (' ',I6,' DY',F15.4,F17.4,F20.4,F9.4,F10.4,
** v 4.29i
*    &                    F10.4,t99,i5,2x,A30)
     &                    F10.4,t99,'(',i5,')',1x,A30)
************     
              ELSE
                IF (LNEWX) THEN
** v 4.29f		
*                 WRITE (LUNIT,173) IOBSY,OBS0Y,OBSBY,VY,VSDY,
*    &                          NAME1Y,NAME2Y
* 173             FORMAT (' ',I6,' DY',F15.4,F18.4,F21.4,F7.1,2(1X,A30))
                  WRITE (LUNIT,173) IOBSY,OBS0Y,OBSBY,VY,VSDY,
     &                          isssn,NAME1Y,jsssn,NAME2Y
  173             FORMAT (' ',I6,' DY',F15.4,F18.4,F21.4,F7.1,
** v 4.29i
*    &                   2(2x,i5,2x,a30)) 
     &                   2(1x,'(',i5,')',1x,a30)) 
***************************
  
                ELSE
                  WRITE (LUNIT,17) IOBSY, OBS0Y, OBSBY, VY, VSDY
   17             FORMAT (' ',I6,' DY',F15.4,F18.4,F21.4,F7.1)
                ENDIF
              ENDIF

*** DEL Z GPS

              IF (I .EQ. NR) THEN
                IF (IMODE .EQ. 0) THEN
                   WRITE (LUNIT,181) IOBSZ, OBS0Z, GMDEZ, RNSUM,
     &             SESSID		  
  181             FORMAT (' ',I6,' DZ',F15.4,F12.4,F18.2,T58,A5)
                ELSEIF (IMODE .EQ. 3) THEN
c Mike Potterfield 2/14/06
c The component numbering is revised to give unique numbers to all
c components, including rejected vectors
c                   WRITE (LUNIT,182) IOBSZ,OBS0Z,OBSBZ,VZ,SD3Z,
      WRITE (LUNIT,182) 3+IGPSOS,OBS0Z,OBSBZ,VZ,SD3Z,
     &                         VSD3Z,GMDEZ,RNSUM,SESSID
** v 4.29f  
* 182             FORMAT (' ',I6,' DZ',F15.4,F17.4,F20.4,F9.4,F10.4,
*    &                     F10.4,F6.2,2X,A5)
  182             FORMAT (' ',I6,' DZ',F15.4,F17.4,F20.4,F9.4,F10.4,
     &                     F10.4,F6.2,2X,A5)
*******************************     
     
                ELSE
                  IF (LNEWX) THEN
** v 4.29f		  
*                    WRITE (LUNIT,186) IOBSZ,OBS0Z,OBSBZ,VZ,
*    &                       VSDZ,NAME1Z, NAME2Z,SESSID
* 186               FORMAT (' ',I6,' DZ',F15.4,F18.4,F21.4,F7.1,
*    &                     2(1X,A30), 2X, A5)		  
                     WRITE (LUNIT,186) IOBSZ,OBS0Z,OBSBZ,VZ,
     &                       VSDZ,isssn,NAME1Z,jsssn, NAME2Z,SESSID
  186               FORMAT (' ',I6,' DZ',F15.4,F18.4,F21.4,F7.1,
** v 4.29i
*    &                     2(1X,i5,2x,A30), 2X, A5)		  
     &                     2(1X,'(',i5,')',1x,A30), 2X, A5)		  
********************     
                  ELSE
                     WRITE (LUNIT,18) IOBSZ,OBS0Z,OBSBZ,VZ,VSDZ,
     &                       SESSID		    
   18               FORMAT (' ',I6,' DZ',F15.4,F18.4,F21.4,F7.1,2X,A5)
                  ENDIF
                ENDIF
              ELSE
                IF (IMODE .EQ. 0) THEN
                   WRITE (LUNIT,183) IOBSZ, OBS0Z, GMDEZ,SESSID
  183             FORMAT (' ',I6,' DZ',F15.4,F12.4,T58,A5)
                ELSEIF (IMODE .EQ. 3) THEN
c Mike Potterfield 2/14/06
c The component numbering is revised to give unique numbers to all
c components, including rejected vectors
c                  WRITE (LUNIT,184) IOBSZ,OBS0Z,OBSBZ,VZ,SD3Z,
c Mike Potterfield 3/22/06 component number updated by offset
      WRITE (LUNIT,184)3+IGPSOS,OBS0Z,OBSBZ,VZ,SD3Z,
     &                             VSD3Z,GMDEZ, SESSID		  
  184           FORMAT(' ',I6,' DZ',F15.4,F17.4,F20.4,F9.4,F10.4,F10.4,
     &                 8X,A5)
                ELSE
                  IF (LNEWX) THEN
** v 4.29f		  
*                    WRITE (LUNIT,187) IOBSZ,OBS0Z,OBSBZ,VZ,VSDZ,
*    &                           NAME1Z,NAME2Z,SESSID
* 187               FORMAT (' ',I6,' DZ',F15.4,F18.4,F21.4,F7.1,
*    &                      2(1X,A30), 4X,A5 )
                     WRITE (LUNIT,187) IOBSZ,OBS0Z,OBSBZ,VZ,VSDZ,
     &                           isssn,NAME1Z,jsssn,NAME2Z,SESSID
  187               FORMAT (' ',I6,' DZ',F15.4,F18.4,F21.4,F7.1,
** v 4.29i
*    &                      2(1X,i5,2x,A30), 4X,A5 )
     &                      2(1X,'(',i5,')',1x,A30), 4X,A5 )
*****************************     
                  ELSE
                     WRITE (LUNIT,185) IOBSZ, OBS0Z, OBSBZ, VZ,
     &                     VSDZ, SESSID		    
** v 4.29f
* 185               FORMAT (' ',I6,' DZ',F15.4,F18.4,F21.4,F7.1,
*    &                     4X,A5)
  185               FORMAT (' ',I6,' DZ',F15.4,F18.4,F21.4,F7.1,
     &                     2X,A5)
********************     
     
                  ENDIF
                ENDIF
          endif
	   
	  if(.not.lgv) then 
	    write(lunit,3)
    3       format(' ')	    
	  else  

*** PRINT N,E,U RESIDUALS AFTER Z COMPONENT

                CALL TOLGH (OBS0X, OBS0Y, OBS0Z, OBS0N, OBS0E, OBS0U,
     &                      ISN, JSN, B)
                CALL TOLGH (OBSBX, OBSBY, OBSBZ, OBSBN, OBSBE, OBSBU,
     &                      ISN, JSN, B)
                CALL TOLGH (VX, VY, VZ, VN, VE, VU, ISN, JSN, B)
***v 4.27
c Mike Potterfield 2/2/06 XYZCOV and RM are added to the parameters
                CALL COVLGH (SN,SE,SU,ISN,JSN,B,I,G,NR,NC,LENG,
     &                   XYZCOV,RM)
****           		
		
                CALL TOLGH (VSDX, VSDY, VSDZ, VSDN, VSDE, VSDU,
     &                      ISN, JSN, B)
     
* v 4.28h
*               OBS0H = DSQRT(OBS0N**2 + OBS0E**2)
*               OBSBH = DSQRT(OBSBN**2 + OBSBE**2)
** v 4.29g		
**               VL = DSQRT(VN**2 + VE**2)
**              SH = DSQRT(SN**2 + SE**2)
*               VL = OBS0H - OBSBH 
*               SH = DSQRT(SN**2 + SE**2)
**********************		
** v 4.29k		

                OBS0L = DSQRT(OBS0N**2 + OBS0E**2)
                OBSBL = DSQRT(OBSBN**2 + OBSBE**2)
                VL = DSQRT(VN**2 + VE**2)
                SL = DSQRT(SN**2 + SE**2)
**************************
		   	
                IF (IMODE .EQ. 0) THEN
                  WRITE (LUNIT,191) ' DN', OBS0N
                  WRITE (LUNIT,191) ' DE', OBS0E
* v 4.28h
		  
*                  write( lunit,191) ' DL', OBS0H	  
**********
* v 4.28k
                  write( lunit,191) ' DL', OBS0L	  
*****************		  
		  
		  
                  WRITE (LUNIT,191) ' DU', OBS0U
  191             FORMAT (' ', 6X, A3, F15.4)
                  write(lunit,3)
                ELSEIF (IMODE .EQ. 3) THEN
**4.27
c Mike Potterfield 2/2/06
c In the new and revised code below, the variables
c VSDN, VSDE, VSDU, and VSDL have been resurrected
c to fulfill their originally intended roles to 
c express the standard deviations of the residuals
c in delta North, East, Length and Up.
                  IF (TNEUCV(1,1) .GT. 0.) THEN
                      VSDN = DSQRT(TNEUCV(1,1))
                      ELSE 
                          VSDN = 0
                      ENDIF
                  IF (TNEUCV(2,2) .GT. 0.) THEN
                      VSDE = DSQRT(TNEUCV(2,2))
                      ELSE 
                           VSDE = 0
                      ENDIF
                  IF (TNEUCV(3,3) .GT. 0.) THEN
                      VSDU = DSQRT(TNEUCV(3,3))
                      ELSE 
                          VSDU = 0
                      ENDIF
                  VSDL = TNEUCV(1,1) + TNEUCV(2,2)
                  IF (VSDL .GT. 0.) THEN
                       VSDL = DSQRT(VSDL)
                       ELSE 
                             VSDL = 0.
                       ENDIF
		
                 WRITE (LUNIT,192) ' DN', OBS0N, OBSBN, VN, VSDN
c                  WRITE (LUNIT,192) ' DN', OBS0N, OBSBN, VN, SN   
**v4.27		  
c                 WRITE (LUNIT,192) ' DE', OBS0E, OBSBE, VE, VSDE
c Output the vector and solution numbers on the DE line
      WRITE(LUNIT,392)' DE', OBS0E, OBSBE, VE, VSDE,IVCNUM,IGPSSL
c                 WRITE (LUNIT,192) ' DE', OBS0E, OBSBE, VE, SE  
** v 4.29k		  
		  
*                 WRITE (LUNIT,192) ' DL', OBS0H, OBSBH, VL, SH  
c                  WRITE (LUNIT,292) ' DL',  VL, VSDL
c Output the Project ID on the DL line  
                  WRITE (LUNIT,393) ' DL',  VL, VSDL, TPRJID  
***********
**v4.27
c If the vector is rejected, identify it as such, otherwise
c no output
      if (CREJCT .EQ. ' ') THEN
                 WRITE (LUNIT,192) ' DU', OBS0U, OBSBU, VU, VSDU
      ELSE
        WRITE (LUNIT, 394) ' DU', OBS0U, OBSBU, VU, VSDU,CREJCT
      ENDIF
c                 WRITE (LUNIT,192) ' DU', OBS0U, OBSBU, VU, SU  
		  
c Mike Potterfield 2/2/06 
c This concludes the revisions required to output correct values
c for standard deviations in residuals for dNELU

  192             FORMAT (' ',6X, A3, F15.4, F17.4, F20.4, F9.4)
  392 FORMAT (' ',6X, A3, F15.4, F17.4, F20.4, F9.4,14X,'Vector',
     &         I4,' Solution',I7)

** v 4.29k
c  292             FORMAT (' ',6X, A3, 32x , F20.4, F9.4)
  393             FORMAT (' ',6X, A3, 32x , F20.4, F9.4,14X,
     & 'Project ID = ',13A)

  394   FORMAT(' ',6X, A3, F15.4, F17.4, F20.4, F9.4,14X,
     &         'Rejected vector code = ',A1)
  


                  write(lunit,3)
                ELSE
**v 4.27		
*                 WRITE (LUNIT,193) ' DN', OBS0N, OBSBN, VN, VSDN
                  WRITE (LUNIT,193) ' DN', OBS0N, OBSBN, VN, SN  

**v4.27
*                 WRITE (LUNIT,193) ' DE', OBS0E, OBSBE, VE, VSDE
                  WRITE (LUNIT,193) ' DE', OBS0E, OBSBE, VE, SE  
* v 4.28m
** v 4.29k
*                 WRITE (LUNIT,193) ' DL', OBS0H, OBSBH, VL, SH  
                  WRITE (LUNIT,293) ' DL',  VL, SL  
*******		  
		  
**v4.27
*                 WRITE (LUNIT,193) ' DU', OBS0U, OBSBU, VU, VSDU
                  WRITE (LUNIT,193) ' DU', OBS0U, OBSBU, VU, SU  

  193             FORMAT (' ', 6X, A3, F15.4, F18.4, F21.4, F7.1)
** v 4.29k
  293             FORMAT (' ', 6X, A3, 33x, F21.4, F7.1)
  
                  write(lunit,3)
  
                ENDIF
             endif     
           endif 
	   lc = .false.
	   
** v 4.29f
            iobsl = iobs
**************	    	   

         endif	      
        endif
  500 CONTINUE

      RETURN

*** PREMATURE END OF FILE

  666 WRITE (LUNIT,667) NVEC
  667 FORMAT ('0PREMATURE FILE END IN RGPS -- NVEC=',I5)
          write (lunit,*) 'subs2, 63'
      CALL ABORT2
      RETURN
      END
      SUBROUTINE ROBS (KIND,IOBS,ISN,ISNX,JSN,JSNX,KSNX,OBSB,SD,B,IAUX,
     *                 IVF,SIGUWT,IC,C,LENG,A,NX,IGRT)

*** LIST ADJUSTED OBSERVATIONS AND RESIDUALS FOR ALL OBS BUT GPS AND DOP

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, LENC = 10, MXVF = 40 )
      CHARACTER*30 NAMES,NAME1,NAME2,NAME3
      CHARACTER*1 ADIR1,ADIR2,ADIR3
      LOGICAL PROP
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL L2HLF, LEHT
      
* v 4.28i
      logical lc
      
** v 4.29f
      logical locssn
***********

      DIMENSION B(*),A(*),NX(*),IC(LENC),C(LENC)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON /VFTB2/  VFS(MXVF), VTV(MXVF), VFRN(MXVF), VFSS(MXVF)
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT
      SAVE GLAT,GLON,EHT,RX,RY
      SAVE NAME1,NAME2,NAME3
* v 4.28i
      save lc      
* lc = true, residual exceeds tolerance, print residual
*** CHECK FOR NEW ISN
** v 4.29f
      save isssn,jsssn,ksssn
**********      

      IF ( ISNX .NE. ISN  .AND.  KIND .GT. 0  .AND.  KIND .LT. 7 ) THEN
        ISNX = ISN
	if(lc) then
        CALL LINE2 (1)
        WRITE (LUNIT,1)
    1   FORMAT (' ')
        endif
        CALL GETGLA (GLAT,ISN,B)
        CALL GETGLO (GLON,ISN,B)
        IF (L2HLF  .AND.  I2HLF(ISN) .EQ. 1) THEN
          CALL GETEHT (EHT,ISN,B)
        ELSE
          CALL GETMSL (GMSL,ISN,B)
          CALL GETGH (GHT,ISN,B)
          EHT = GMSL+GHT
        ENDIF
        CALL RADII (ISN,RMER,RPV,B)
        RX = RMER+EHT
        RY = (RPV+EHT)*DCOS(GLAT)
        NAME1 = NAMES(ISN)
** v 4.29f
        if(.not.locssn(isn,isssn)) then
	  write(lunit,966) isn
 966      format('0SSN TABLE ERROR IN ROBS-- ',i8)
             write (lunit,*) 'subs2, 3'
          call abort2
	endif
********************	  	  	

        IF ( KIND .GT. 3  .AND.  IMODE .EQ. 3 ) THEN
          CALL LINE2 (1)
** v 4.29f	  
*         WRITE (LUNIT,510) NAME1
*510      FORMAT (93X,A30)
          WRITE (LUNIT,510) isssn,NAME1
** v 4.29i	  
*510      FORMAT (t99,i5,2x,A30)
 510      FORMAT (t99,'(',i5,')',1x,A30)
************************** 
 
        ENDIF
      ENDIF
      IF ( ISNX .NE. ISN  .AND.  KIND .GE. 7  ) THEN
        ISNX = ISN
	if(lc) then
        WRITE (LUNIT,2)
 2      format(' ')
        CALL LINE2 (1)
	endif
        NAME1 = NAMES(ISN)
	
** v 4.29f
        if(.not.locssn(isn,isssn)) then
	  write(lunit,966) isn
             write (lunit,*) 'subs2, 4'
          call abort2
	endif
********************	  	  	
	
        IF ( IMODE .EQ. 3 ) THEN
          CALL LINE2 (1)
** v 4.29f	  
*         WRITE (LUNIT,510) NAME1
          WRITE (LUNIT,510) isssn,NAME1
*****************	  
	  
        ENDIF
      ENDIF
      IF ( JSNX .NE. JSN  .AND.  KIND .GE. 7 ) THEN
        JSNX = JSN
        NAME2 = NAMES(JSN)
	
** v 4.29f
        if(.not.locssn(jsn,jsssn)) then
	  write(lunit,966) isn
             write (lunit,*) 'subs2, 5'
          call abort2
	endif
********************	  	  	
	
      ENDIF
      IF ( KIND .EQ. 10  .AND.  KSNX .NE. IAUX ) THEN
        KSNX = IAUX
        NAME3 = NAMES(IAUX)
	
** v 4.29f
        if(.not.locssn(iaux,ksssn)) then
	  write(lunit,966)iaux 
             write (lunit,*) 'subs2, 6'
          call abort2
	endif
********************	  	  	
	
      ENDIF

*** SCALE BY THE VARIANCE FACTOR

      IF ( NVFTOT .GT. 0 ) THEN
        IF (IVF .LT. 0  .OR.  IVF .GT. NVFTOT ) THEN
          WRITE (LUNIT,9) IVF,NVFTOT
    9     FORMAT ('0ILLEGAL IVF=',I10,' FOR N=',I10,' IN ROBS')
          write (lunit,*) 'subs2, 64'
          CALL ABORT2
        ELSEIF ( IVF .NE. 0 ) THEN
          SD = SD*DSQRT(VFS(IVF))
        ENDIF
      ENDIF

*** COMPUTATION OF S.D. FOR RESIDUALS AND
***  COMPUTATION OF REDUNDANCY NUMBERS

      IF ( IMODE .EQ. 0  .OR.  IMODE .EQ. 3 ) THEN
        IF ( PROP(C,IC,LENG,VAR,A,NX,IFLAG) ) THEN
          DVAR = SD*SD-VAR
          IF ( (KIND .GE. 8  .AND.  KIND .LE. 12)  .OR.
     &         KIND .EQ. 14 ) THEN
            IF (DVAR .LT. 1.0D-14) DVAR = 0.D0
          ELSE
            IF (DVAR .LT. 1.0D-10) DVAR = 0.D0
          ENDIF
          SD3 = DSQRT(DVAR)
          RN = 1.D0-VAR/(SD*SD)
          IF (RN .LT. 1.0D-10) RN = 0.D0
        ELSEIF (IFLAG .EQ. 1) THEN
          WRITE (LUNIT,660)
 660      FORMAT ('0SYSTEM NOT INVERTED IN ROBS'/)
          write (lunit,*) 'subs2, 65'
          CALL ABORT2
        ELSE
          WRITE (LUNIT,661)
 661      FORMAT ('0ALL COVARIANCE ELEMENTS NOT WITHIN PROFILE--ROBS')
          write (lunit,*) 'subs2, 66'
          CALL ABORT2
        ENDIF
      ENDIF
      IF (IMODE .GE. 0  .AND.  IMODE .LE. 2) SD3 = SD

***COMPUTE THE MARGINAL DETECTABLE ERROR(GMDE) USING 3*SIGMA

      IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) THEN
        IF (RN .LT. 1.0D-8) THEN
          GMDE = 0.D0
        ELSE
          GMDE = (3.D0*SD)/DSQRT(RN)
        ENDIF
      ENDIF

*** AUXILIARY PARAMETER CONSTRAINTS

      IF (KIND .EQ. 0) THEN
        IAUX = ISN
        CALL GETAUX (AUX,IAUX,B)
        V = AUX-OBSB
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************

        IF (LVC) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
* v 4.28i
            lc = .true.	  
            CALL LINE2 (2)
            IF (IMODE .EQ. 0) THEN
              WRITE (LUNIT,11) IOBS,AUX,GMDE,RN,IAUX
 11           FORMAT ('0',I6,' AP',1PD17.2,1PD12.2,0PF16.2,
     *                 ' AUX PARM #',I4)
            ELSEIF (IMODE .EQ. 3) THEN
              WRITE (LUNIT,12) IOBS,AUX,OBSB,V,SD3,VSD3,GMDE,RN,IAUX
 12           FORMAT ('0',I6,' AP',1PD17.2,1PD17.2,1PD17.2,0PF8.2,F8.2,
     &                F10.4,F6.2,'  AUX PARM #',I4)
            ELSE
              WRITE (LUNIT,10) IOBS,AUX,OBSB,V,VSD,IAUX
   10         FORMAT ('0',I6,' AP',1PD17.2,1PD18.2,1PD18.2,0PF8.1,
     *                ' AUX PARM #',I4)
            ENDIF
	    
* v 4.28i
          else
	    lc = .false.	    
	    
          ENDIF
        ENDIF

*** LATITUDE CONSTRAINT

      ELSEIF (KIND .EQ. 1) THEN
        CALL GETDMS (GLAT,ID1,IM1,S1,ISIGN)
        IF (ISIGN .GT. 0) THEN
          ADIR1 = 'N'
        ELSE
          ADIR1 = 'S'
        ENDIF
        CALL GETDMS (OBSB,ID3,IM3,S3,ISIGN)
        IF (ISIGN .GT. 0) THEN
          ADIR3 = 'N'
        ELSE
          ADIR3 = 'S'
        ENDIF
        V = GLAT-OBSB
        VSEC = V*RAD*3600.D0
        VMET = V*RX
        VSD = VMET/SD
        VSD3 = DIVID(VMET,SD3)
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        IF (LVC) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,131) IOBS,ID1,IM1,S1,ADIR1,GMDE,RN,NAME1
*131          FORMAT (' ',I6,' LA',I4,I3,F9.5,A1,F9.3,F19.2,1X,A30)
              WRITE (LUNIT,131) IOBS,ID1,IM1,S1,ADIR1,GMDE,RN,
     &	            isssn,name1 
** v 4.29i
*131          FORMAT (' ',I6,' LA',I4,I3,F9.5,A1,F9.3,F19.2,1x,i5,
*    &               2x,a30)
 131          FORMAT (' ',I6,' LA',I4,I3,F9.5,A1,F9.3,F19.2,1x,'(',i5,
     &               ')',1x,a30)
*************** 
 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,132) IOBS,ID1,IM1,S1,ADIR1,ID3,IM3,S3,ADIR3,
*    &                      VSEC,VMET,SD3,VSD3,GMDE,RN,NAME1
*132          FORMAT (' ',I6,' LA',I4,I3,F9.5,A1,I4,I3,F9.5,A1,F9.5,
*    &                F8.3,F8.2,F8.2,F10.4,F6.2,2X,A30)
     
              WRITE (LUNIT,132) IOBS,ID1,IM1,S1,ADIR1,ID3,IM3,S3,ADIR3,
     &                      VSEC,VMET,SD3,VSD3,GMDE,RN,isssn,NAME1
 132          FORMAT (' ',I6,' LA',I4,I3,F9.5,A1,I4,I3,F9.5,A1,F9.5,
** v 4.29i
*    &                F8.3,F8.2,F8.2,F10.4,F6.2,t99,i5,2x,A30)
     &                F8.3,F8.2,F8.2,F10.4,F6.2,t99,'(',i5,')',1x,A30)
**********************************    
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,13) IOBS,ID1,IM1,S1,ADIR1,ID3,IM3,S3,ADIR3,
*    *                     VSEC,VMET,VSD,NAME1
*  13         FORMAT (' ',I6,' LA',I4,I3,F9.5,A1,I5,I3,F9.5,A1,
*    *                F10.5,F8.3,F8.1,1X,A30)
              WRITE (LUNIT,13) IOBS,ID1,IM1,S1,ADIR1,ID3,IM3,S3,ADIR3,
     *                     VSEC,VMET,VSD,isssn,NAME1
   13         FORMAT (' ',I6,' LA',I4,I3,F9.5,A1,I5,I3,F9.5,A1,
** v 4.29i
*    *                F10.5,F8.3,F8.1,1x,i5,2x,A30)
     *                F10.5,F8.3,F8.1,1x,'(',i5,')',1x,A30)
******************
     
            ENDIF
	  else
	    lc = .false.  
          ENDIF
        ENDIF

*** LONGITUDE CONSTRAINT

      ELSEIF (KIND .EQ. 2) THEN
        CALL GETDMS (GLON,ID2,IM2,S2,ISIGN)
        IF (ISIGN .GT. 0) THEN
              ALONGD=ID2+(IM2/60.D0)+(S2/3600.D0)
              ALONGD=360.D0-ALONGD
              ID2=INT(ALONGD)
              ALONGM=(ALONGD-ID2) * 60.D0
              IM2=INT(ALONGM)
              S2=(ALONGM-IM2) * 60.D0
              ADIR2 = 'W'
        ELSE
          ADIR2 = 'W'
        ENDIF
        CALL GETDMS (OBSB,ID3,IM3,S3,ISIGN)
        IF (ISIGN .GT. 0) THEN
              ALONGD=ID3+(IM3/60.D0)+(S3/3600.D0)
              ALONGD=360.D0-ALONGD
              ID3=INT(ALONGD)
              ALONGM=(ALONGD-ID3) * 60.D0
              IM3=INT(ALONGM)
              S3=(ALONGM-IM3) * 60.D0
              ADIR3 = 'W'
        ELSE
          ADIR3 = 'W'
        ENDIF
        V = GLON-OBSB
        VSEC = V*RAD*3600.D0
        VMET = V*RY
        VSD = VMET/SD
        VSD3 = DIVID(VMET,SD3)
        IF (LSS) VSD = VSD/SIGUWT
       
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        IF (LVC) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,141) IOBS,ID2,IM2,S2,ADIR2,GMDE,RN,NAME1
*141          FORMAT (' ',I6,' LO',I4,I3,F9.5,A1,F9.3,F19.2,1X,A30)
              WRITE (LUNIT,141) IOBS,ID2,IM2,S2,ADIR2,GMDE,RN,
     &                          isssn,name1	      
 141          FORMAT (' ',I6,' LO',I4,I3,F9.5,A1,F9.3,F19.2,1x,
** v 4.29i
*    &                i5,2x,a30)
     &                '(',i5,')',1x,a30)
**************************** 
 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,142) IOBS,ID2,IM2,S2,ADIR2,ID3,IM3,S3,ADIR3,
*    &                      VSEC,VMET,SD3,VSD3,GMDE,RN,NAME1
*142          FORMAT (' ',I6,' LO',I4,I3,F9.5,A1,I4,I3,F9.5,A1,F9.5,
*    &                F8.3,F8.2,F8.2,F10.4,F6.2,2X,A30)
              WRITE (LUNIT,142) IOBS,ID2,IM2,S2,ADIR2,ID3,IM3,S3,ADIR3,
     &                      VSEC,VMET,SD3,VSD3,GMDE,RN,isssn,NAME1
 142          FORMAT (' ',I6,' LO',I4,I3,F9.5,A1,I4,I3,F9.5,A1,F9.5,
** v 4.29i
*    &                F8.3,F8.2,F8.2,F10.4,F6.2,t99,i5,2x,A30)
     &                F8.3,F8.2,F8.2,F10.4,F6.2,t99,'(',i5,')',1x,A30)
***********************************     
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,14) IOBS,ID2,IM2,S2,ADIR2,ID3,IM3,S3,ADIR3,
*    *                     VSEC,VMET,VSD,NAME1
*  14         FORMAT (' ',I6,' LO',I4,I3,F9.5,A1,I5,I3,F9.5,A1,
*    *                F10.5,F8.3,F8.1,1X,A30)
              WRITE (LUNIT,14) IOBS,ID2,IM2,S2,ADIR2,ID3,IM3,S3,ADIR3,
     *                     VSEC,VMET,VSD,isssn,NAME1
   14         FORMAT (' ',I6,' LO',I4,I3,F9.5,A1,I5,I3,F9.5,A1,
** v 4.29i
*    *                F10.5,F8.3,F8.1,1x,i5,2x,A30)
     *                F10.5,F8.3,F8.1,1x,'(',i5,')',1x,A30)
***************************
     
            ENDIF
	  else
	    lc = .false.  
          ENDIF
        ENDIF

*** HEIGHT CONSTRAINT

      ELSEIF (KIND .EQ. 3) THEN
        V = EHT-OBSB
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        IF (LVC) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,151) IOBS,EHT,GMDE,RN,NAME1
*151          FORMAT (' ',I6,' EH',F14.3,F12.3,F19.2,1X,A30)
              WRITE (LUNIT,151) IOBS,EHT,GMDE,RN,isssn,NAME1
** v 4.29i	      
*151          FORMAT (' ',I6,' EH',F14.3,F12.3,F19.2,1x,i5,2x,A30)
 151          FORMAT (' ',I6,' EH',F14.3,F12.3,F19.2,1x,'(',i5,')',1x,
     *                A30)
*************************
 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,152) IOBS,EHT,OBSB,V,SD3,VSD3,GMDE,RN,NAME1          
*152          FORMAT (' ',I6,' EH',F14.3,F17.3,F20.3,F8.2,F8.2,F10.4,
*    &                F6.2,2X,A30)
              WRITE (LUNIT,152) IOBS,EHT,OBSB,V,SD3,VSD3,GMDE,RN,         
     &                          isssn,name1	      
 152          FORMAT (' ',I6,' EH',F14.3,F17.3,F20.3,F8.2,F8.2,F10.4,
** v 4.29i
*    &                F6.2,t99,i5,2x,A30)
     &                F6.2,t99,'(',i5,')',1x,A30)
****************
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,15) IOBS,EHT,OBSB,V,VSD,NAME1
*  15         FORMAT (' ',I6,' EH',F14.3,F18.3,F21.3,F8.1,1X,A30)
              WRITE (LUNIT,15) IOBS,EHT,OBSB,V,VSD,isssn,NAME1
   15         FORMAT (' ',I6,' EH',F14.3,F18.3,F21.3,F8.1,1x,
** v 4.29i
*    &                i5,2x,a30)
     &                '(',i5,')',1x,a30)
*************************************
   
            ENDIF
	  else
	    lc = .false.
          ENDIF
        ENDIF

*** MARK TO MARK DISTANCE

      ELSEIF (KIND .EQ. 7) THEN
        CALL COMPOB (KIND,S,B,OBSB,ISN,JSN,IAUX,IGRT)
        V = S-OBSB
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        IF (LVS) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,171) IOBS,S,GMDE,RN,NAME1,NAME2
*171          FORMAT (' ',I6,' DIST',F13.4,F11.3,F19.2,2(1X,A30))
              WRITE (LUNIT,171) IOBS,S,GMDE,RN,isssn,NAME1,jsssn,NAME2
** v 4.29i	      
*171          FORMAT (' ',I6,' DIST',F13.4,F11.3,F19.2,2(1X,i5,2x,A30))
 171          FORMAT (' ',I6,' DIST',F13.4,F11.3,F19.2,2(1X,'(',i5,')',
     *                1x,A30))
****************************
 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,172) IOBS,S,OBSB,V,SD3,VSD3,GMDE,RN,NAME2
*172          FORMAT (' ',I6,' DIST',F13.4,F17.4,F19.3,F8.2,F8.2,F10.4,
*    &                F6.2,4X,A30)
              WRITE (LUNIT,172) IOBS,S,OBSB,V,SD3,VSD3,GMDE,RN,
     &	                        jsssn,name2
 172          FORMAT (' ',I6,' DIST',F13.4,F17.4,F19.3,F8.2,F8.2,F10.4,
** v 4.20i
*    &                F6.2,t99,i5,2x,A30)
     &                F6.2,t99,'(',i5,')',1x,A30)
*************************
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,17) IOBS,S,OBSB,V,VSD,NAME1,NAME2
* 17          FORMAT (' ',I6,' DIST',F12.3,F18.3,F21.3,
*    *                F8.1,2(1X,A30))
              WRITE (LUNIT,17) IOBS,S,OBSB,V,VSD,
     &                         isssn,name1,jsssn,name2	      
  17          FORMAT (' ',I6,' DIST',F12.3,F18.3,F21.3,
** v 4.29i
*    *                F8.1,2(1X,i5,2x,A30))
     *                F8.1,2(1X,'(',i5,')',1x,A30))
*****************************     
     
            ENDIF
	  else
	    lc = .false.
          ENDIF
        ENDIF

*** ASTRONOMIC AZIMUTH

      ELSEIF (KIND .EQ. 8) THEN
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        CALL DIRDMS (OBS0,ID0,IM0,S0)
        V = OBS0-OBSB
        IF (V .GT. PI) THEN
          V = V-PI-PI
        ELSEIF (V .LT. -PI) THEN
          V = V+PI+PI
        ENDIF
        VSEC = V*RAD*3600.D0
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) GMDE = GMDE*RAD*3600.D0
        IF (IMODE .EQ. 3) SD3 = SD3*RAD*3600.D0
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        CALL DIRDMS (OBSB,IDB,IMB,SB)
        IF (LVR) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,191) IOBS,ID0,IM0,S0,GMDE,RN,NAME1,NAME2
*191          FORMAT (' ',I6,' AZ',I4,I3,F6.2,'N',F11.2,F20.2,2(1X,A30))
              WRITE (LUNIT,191) IOBS,ID0,IM0,S0,GMDE,RN,
     &                           isssn,name1,jsssn,name2	      
 191          FORMAT (' ',I6,' AZ',I4,I3,F6.2,'N',F11.2,F20.2,
** v 4.29i
*    &                2(1x,i5,2x,a30))
     &                2(1x,'(',i5,')',1x,a30))
****************************
 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,192)IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,SD3,VSD3,
*    &                      GMDE,RN,NAME2
*192          FORMAT (' ',I6,' AZ',I4,I3,F6.2,'N',I7,I3,F6.2,'N',
*    &                F9.2,F19.2,F8.2,F8.2,F8.2,4X,A30)
              WRITE (LUNIT,192)IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,SD3,VSD3,
     &                      GMDE,RN,jsssn,NAME2
 192          FORMAT (' ',I6,' AZ',I4,I3,F6.2,'N',I7,I3,F6.2,'N',
** v 4.29i
*    &                F9.2,F19.2,F8.2,F8.2,F8.2,t99,i5,2x,A30)
     &                F9.2,F19.2,F8.2,F8.2,F8.2,t99,'(',i5,')',1x,A30)
**************************     
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,19) IOBS,ID0,IM0,S0,IDB,IMB,SB,
*    &                     VSEC,VSD,NAME1,NAME2
*19           FORMAT (' ',I6,' AZ',I4,I3,F6.2,'N',I8,I3,F5.1,
*    &                'N',4X,F10.5,8X,F8.1,2(1X,A30))
              WRITE (LUNIT,19) IOBS,ID0,IM0,S0,IDB,IMB,SB,
     &                     VSEC,VSD,isssn,NAME1,jsssn,NAME2
 19           FORMAT (' ',I6,' AZ',I4,I3,F6.2,'N',I8,I3,F5.1,
** v 4.29i
*    &                'N',4X,F10.5,8X,F8.1,2(1X,i5,2x,A30))
     &                'N',4X,F10.5,8X,F8.1,2(1X,'(',i5,')',1x,A30))
*****************************
     
            ENDIF
	  else
	    lc = .false.
          ENDIF
        ENDIF

*** ZENITH DISTANCE

      ELSEIF (KIND .EQ. 9) THEN
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        CALL VERDMS (OBS0,ID0,IM0,S0,ISIGN)
        V = OBS0-OBSB
        VSEC = V*RAD*3600.D0
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) GMDE = GMDE*RAD*3600.D0
        IF (IMODE .EQ. 3) SD3 = SD3*RAD*3600.D0
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        CALL VERDMS (OBSB,IDB,IMB,SB,ISIGN)
        IF (LVZ) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,211) IOBS,ID0,IM0,S0,GMDE,RN,NAME1,NAME2
*211          FORMAT (' ',I6,' ZD',I4,I3,F6.2,' ',F11.2,F20.2,2(1X,A30))
              WRITE (LUNIT,211) IOBS,ID0,IM0,S0,GMDE,RN,
     &                           isssn,name1,jsssn,name2	      
 211          FORMAT (' ',I6,' ZD',I4,I3,F6.2,' ',F11.2,F20.2,
** v 4.29i
*    &                  2(1x,i5,2x,a30))
     &                  2(1x,'(',i5,')',1x,a30))
*************************
 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,212)IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,SD3,VSD3,
*    &                      GMDE,RN,NAME2
*212          FORMAT (' ',I6,' ZD',I4,I3,F6.2,' ',I7,I3,F6.2,
*    &                F10.2,F19.2,F8.2,F8.2,F8.2,4X,A30)
              WRITE (LUNIT,212)IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,SD3,VSD3,
     &                      GMDE,RN,jsssn,NAME2
 212          FORMAT (' ',I6,' ZD',I4,I3,F6.2,' ',I7,I3,F6.2,
** v 4.29i
*    &                F10.2,F19.2,F8.2,F8.2,F8.2,t99,i5,2x,A30)
     &                F10.2,F19.2,F8.2,F8.2,F8.2,t99,'(',i5,')',1x,A30)
*******************************
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,21) IOBS,ID0,IM0,S0,IDB,IMB,SB,
*    &                     VSEC,VSD,NAME1,NAME2
* 21          FORMAT (' ',I6,' ZD',I4,I3,F6.2,I9,I3,F6.2,
*    &                4X,F10.5,8X,F8.1,2(1X,A30))
              WRITE (LUNIT,21) IOBS,ID0,IM0,S0,IDB,IMB,SB,
     &                     VSEC,VSD,isssn,NAME1,jsssn,NAME2
  21          FORMAT (' ',I6,' ZD',I4,I3,F6.2,I9,I3,F6.2,
** v 4.29i
*    &                4X,F10.5,8X,F8.1,2(1X,i5,2x,A30))
     &                4X,F10.5,8X,F8.1,2(1X,'(',i5,')',1x,A30))
******************************
            ENDIF
	  else
	    lc = .false.
          ENDIF
        ENDIF

*** HORIZONTAL ANGLE

      ELSEIF (KIND .EQ. 10) THEN
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        CALL DIRDMS (OBS0,ID0,IM0,S0)
        V = OBS0-OBSB
        IF (V .GT. PI) THEN
          V = V-PI-PI
        ELSEIF (V .LT. -PI) THEN
          V = V+PI+PI
        ENDIF
        VSEC = V*RAD*3600.D0
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) GMDE = GMDE*RAD*3600.D0
        IF (IMODE .EQ. 3) SD3 = SD3*RAD*3600.D0
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        CALL DIRDMS (OBSB,IDB,IMB,SB)
        IF (LVA) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (2)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,231)IOBS,ID0,IM0,S0,GMDE,RN,NAME1,NAME2,NAME3 
*231          FORMAT (' ',I6,' HA',I4,I3,F6.2,' ',F11.2,F20.2,2(1X,A30)/
*    &                T86,A30)
              WRITE (LUNIT,231)IOBS,ID0,IM0,S0,GMDE,RN, 
     &                         isssn,name1,jsssn,name2,ksssn,name3	      
 231          FORMAT (' ',I6,' HA',I4,I3,F6.2,' ',F11.2,F20.2,2(1X,A30)/
** v 4.29i
*    &                2(1x,i5,2x,a30)/ T86,i5,2x,A30)
     &                2(1x,'(',i5,')',1x,a30)/ T86,'(',i5,')',1x,A30)
*****************************************
     
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,232)IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,SD3,VSD3,
*    &                      GMDE,RN,NAME2,NAME3
*232          FORMAT (' ',I6,' HA',I4,I3,F6.2,' ',I7,I3,F6.2,
*    &                F10.2,F19.2,F8.2,F8.2,F8.2,4X,A30/T96,A30)
              WRITE (LUNIT,232)IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,SD3,VSD3,
     &                      GMDE,RN,jsssn,NAME2,ksssn,NAME3
 232          FORMAT (' ',I6,' HA',I4,I3,F6.2,' ',I7,I3,F6.2,
** v 4.29i
*    &          F10.2,F19.2,F8.2,F8.2,F8.2,t99,i5,2x,A30/T99,i5,2x,A30)
     &          F10.2,F19.2,F8.2,F8.2,F8.2,t99,'(',i5,')',1x,A30/T99,
     *          '(',i5,')',1x,A30)
*******************************     
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,23) IOBS,ID0,IM0,S0,IDB,IMB,SB,
*    &                     VSEC,VSD,NAME1,NAME2,NAME3
*23           FORMAT (' ',I6,' HA',I4,I3,F6.2,' ',I8,I3,F6.2,
*    &                4X,F10.5,8X,F8.1,2(1X,A30)/T102,A30)
              WRITE (LUNIT,23) IOBS,ID0,IM0,S0,IDB,IMB,SB,
     &               VSEC,VSD,isssn,NAME1,jsssn,NAME2,ksssn,NAME3
 23           FORMAT (' ',I6,' HA',I4,I3,F6.2,' ',I8,I3,F6.2,
** v 4.29i
*    &           4X,F10.5,8X,F8.1,2(1X,i5,2x,A30)/T102,i5,2x,A30)
     &           4X,F10.5,8X,F8.1,2(1X,'(',i5,')',1x,A30)/T102,
     *            '(',i5,')',1x,A30)
****************
            ENDIF
	  else
	    lc = .false.
          ENDIF
        ENDIF

*** HORIZONTAL DIRECTION

      ELSEIF (KIND .EQ. 11) THEN
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        CALL DIRDMS (OBS0,ID0,IM0,S0)
        V = OBS0-OBSB
        IF (V .GT. PI) THEN
          V = V-PI-PI
        ELSEIF (V .LT. -PI) THEN
          V = V+PI+PI
        ENDIF
        VSEC = V*RAD*3600.D0
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) GMDE = GMDE*RAD*3600.D0
        IF (IMODE .EQ. 3) SD3 = SD3*RAD*3600.D0
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        CALL DIRDMS (OBSB,IDB,IMB,SB)
        IF (LVD) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,251) IOBS,ID0,IM0,S0,GMDE,RN,NAME1,NAME2
*251          FORMAT (' ',I6,' HD',I4,I3,F6.2,F12.2,F20.2,2(1X,A30))
              WRITE (LUNIT,251) IOBS,ID0,IM0,S0,GMDE,RN,
     &                         isssn,name1,jsssn,name2	      
 251          FORMAT (' ',I6,' HD',I4,I3,F6.2,F12.2,F20.2,
** v 4.29i
*    &               2(1x,i5,2x,a30)) 
     &               2(1x,'(',i5,')',1x,a30)) 
********************* 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,252)IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,SD3,VSD3,
*    &                      GMDE,RN,NAME2
*252          FORMAT (' ',I6,' HD',I4,I3,F6.2,I8,I3,F6.2,
*    &                F10.2,F19.2,F8.2,F8.2,F8.2,4X,A30)
              WRITE (LUNIT,252)IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,SD3,VSD3,
     &                      GMDE,RN,jsssn,NAME2
 252          FORMAT (' ',I6,' HD',I4,I3,F6.2,I8,I3,F6.2,
** v 4.29i
*    &                F10.2,F19.2,F8.2,F8.2,F8.2,t99,i5,2x,A30)
     &                F10.2,F19.2,F8.2,F8.2,F8.2,t99,'(',i5,')',1x,A30)
*******************************
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,25) IOBS,ID0,IM0,S0,IDB,IMB,SB,
*    &                     VSEC,VSD,NAME1,NAME2
*25           FORMAT (' ',I6,' HD',I4,I3,F6.2,' ',I8,I3,F6.2,
*    &                4X,F10.5,8X,F8.1,2(1X,A30))
              WRITE (LUNIT,25) IOBS,ID0,IM0,S0,IDB,IMB,SB,
     &                     VSEC,VSD,isssn,NAME1,jsssn,NAME2
 25           FORMAT (' ',I6,' HD',I4,I3,F6.2,' ',I8,I3,F6.2,
** v 4.29i
*    &                4X,F10.5,8X,F8.1,2(1X,i5,2x,A30))
     &                4X,F10.5,8X,F8.1,2(1X,'(',i5,')',1x,A30))
*********************
            ENDIF
	  else
	    lc = .false.
          ENDIF
        ENDIF

*** CONSTRAINED ASTRONOMIC AZIMUTH

      ELSEIF (KIND .EQ. 12) THEN
        I8 = 8
        CALL COMPOB (I8,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        CALL DIRDMS (OBS0,ID0,IM0,S0)
        V = OBS0-OBSB
        IF (V .GT. PI) THEN
          V = V-PI-PI
        ELSEIF (V .LT. -PI) THEN
          V = V+PI+PI
        ENDIF
        VSEC = V*RAD*3600.D0
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) GMDE = GMDE*RAD*3600.D0
        IF (IMODE .EQ. 3) SD3 = SD3*RAD*3600.D0
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        CALL DIRDMS (OBSB,IDB,IMB,SB)
        IF (LVR) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,271) IOBS,ID0,IM0,S0,GMDE,RN,NAME1,NAME2
*271          FORMAT (' ',I6,' CA',I4,I3,F6.2,'N',F11.2,F20.2,2(1X,A30))
              WRITE (LUNIT,271) IOBS,ID0,IM0,S0,GMDE,RN,
     &                         isssn,name1,jsssn,name2	      
 271          FORMAT (' ',I6,' CA',I4,I3,F6.2,'N',F11.2,F20.2,
** v 4.29i
*    &                2(1x,i5,2x,a30))
     &                2(1x,'(',i5,')',1x,a30))
***************
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE(LUNIT,272)IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,SD3,VSD3,
*    &                      GMDE,RN,NAME2
*272          FORMAT (' ',I6,' CA',I4,I3,F6.2,'N',I7,I3,F6.2,
*    &                'N',2X,F7.2,10X,F9.2,F8.2,F8.2,F8.2,4X,A30)
              WRITE(LUNIT,272)IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,SD3,VSD3,
     &                      GMDE,RN,jsssn,NAME2
 272          FORMAT (' ',I6,' CA',I4,I3,F6.2,'N',I7,I3,F6.2,
** v 4.29i
*    &            'N',2X,F7.2,10X,F9.2,F8.2,F8.2,F8.2,t99,i5,2x,A30)
     &            'N',2X,F7.2,10X,F9.2,F8.2,F8.2,F8.2,t99,'(',i5,')',
     *             1x,A30)
*****************************
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,27) IOBS,ID0,IM0,S0,IDB,IMB,SB,
*    &                     VSEC,VSD,NAME1,NAME2
*27           FORMAT (' ',I6,' CA',I4,I3,F6.2,'N',I8,I3,F6.2,
*    &                'N',3X,F7.2,11X,F8.1,2(1X,A30))
              WRITE (LUNIT,27) IOBS,ID0,IM0,S0,IDB,IMB,SB,
     &                     VSEC,VSD,isssn,NAME1,jsssn,NAME2
 27           FORMAT (' ',I6,' CA',I4,I3,F6.2,'N',I8,I3,F6.2,
** v 4.29i
*    &                'N',3X,F7.2,11X,F8.1,2(1X,i5,2x,A30))
     &                'N',3X,F7.2,11X,F8.1,2(1X,'(',i5,')',1x,A30))
*******************
            ENDIF
	  else
	    lc = .false.
          ENDIF
        ENDIF

*** CONSTRAINED MARK TO MARK DISTANCE

      ELSEIF (KIND .EQ. 13) THEN
        I7 = 7
        CALL COMPOB (I7,S,B,OBSB,ISN,JSN,IAUX,IGRT)
        V = S-OBSB
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        IF (LVS) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,291) IOBS,S,GMDE,RN,NAME1,NAME2
*291          FORMAT (' ',I6,' CD',F15.4,F11.3,F19.2,2(1X,A30))
              WRITE (LUNIT,291) IOBS,S,GMDE,RN,
     &                          isssn,name1,jsssn,name2	      
 291          FORMAT (' ',I6,' CD',F15.4,F11.3,F19.2,
** v 4.29i
*    &               2(1x,i5,2x,a30))
     &               2(1x,'(',i5,')',1x,a30))
********************************* 
 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,292) IOBS,S,OBSB,V,SD3,VSD3,GMDE,RN,NAME2
*292          FORMAT (' ',I6,' CD',F15.4,F17.4,F19.3,F8.2,F8.2,
*    &                F10.4,F6.2,4X,A30)
              WRITE (LUNIT,292) IOBS,S,OBSB,V,SD3,VSD3,GMDE,RN,
     &                           jsssn,name2	      
 292          FORMAT (' ',I6,' CD',F15.4,F17.4,F19.3,F8.2,F8.2,
** v 4.29i
*    &                F10.4,F6.2,t99,i5,2x,A30)
     &                F10.4,F6.2,t99,'(',i5,')',1x,A30)
******************************   
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,29) IOBS,S,OBSB,V,VSD,NAME1,NAME2
* 29          FORMAT (' ',I6,' CD ',F14.4,F18.4,F20.3,
*    *                F8.1,2(1X,A30))
              WRITE (LUNIT,29) IOBS,S,OBSB,V,VSD,
     &	                      isssn,name1,jsssn,name2   
  29          FORMAT (' ',I6,' CD ',F14.4,F18.4,F20.3,
** v 4.29i
*    *                F8.1,2(1X,i5,2x,A30))
     *                F8.1,2(1X,'(',i5,')',1x,A30))
*****************
            ENDIF
	  else
	    lc = .false.
          ENDIF
        ENDIF

*** CONSTRAINED ZENITH DISTANCE

      ELSEIF (KIND .EQ. 14) THEN
        I9 = 9
        CALL COMPOB (I9,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        CALL VERDMS (OBS0,ID0,IM0,S0,ISIGN)
        V = OBS0-OBSB
        VSEC = V*RAD*3600.D0
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) GMDE = GMDE*RAD*3600.D0
        IF (IMODE .EQ. 3) SD3 = SD3*RAD*3600.D0
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        CALL VERDMS (OBSB,IDB,IMB,SB,ISIGN)
        IF (LVZ) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,311) IOBS,ID0,IM0,S0,GMDE,RN,NAME1,NAME2
*311          FORMAT (' ',I6,' CZ',I4,I3,F6.2,' ',F11.2,F20.2,2(1X,A30))
              WRITE (LUNIT,311) IOBS,ID0,IM0,S0,GMDE,RN,
     &              isssn,name1,jsssn,name2     
 311          FORMAT (' ',I6,' CZ',I4,I3,F6.2,' ',F11.2,F20.2,
** v 4.29i
*    &               2(1x,i5,2x,a30)) 
     &               2(1x,'(',i5,')',1x,a30)) 
******************************** 
 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE(LUNIT,312)IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,SD3,VSD3,
*    &                      GMDE,RN,NAME2
*312          FORMAT (' ',I6,' CZ',I4,I3,F6.2,' ',I7,I3,F6.2,
*    &                3X,F7.2,10X,F9.2,F8.2,F8.2,F8.2,4X,A30)
              WRITE(LUNIT,312)IOBS,ID0,IM0,S0,IDB,IMB,SB,VSEC,SD3,VSD3,
     &                      GMDE,RN,jsssn,NAME2
 312          FORMAT (' ',I6,' CZ',I4,I3,F6.2,' ',I7,I3,F6.2,
** v 4.29i
*    &             3X,F7.2,10X,F9.2,F8.2,F8.2,F8.2,t99,i5,2x,A30)
     &             3X,F7.2,10X,F9.2,F8.2,F8.2,F8.2,t99,'(',i5,')',1x,
     *             A30)
*****************************************     
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,31) IOBS,ID0,IM0,S0,IDB,IMB,SB,
*    &                     VSEC,VSD,NAME1,NAME2
* 31          FORMAT (' ',I6,' CZ',I4,I3,F6.2,I9,I3,F6.2,
*    &                4X,F7.2,11X,F8.1,2(1X,A30))
              WRITE (LUNIT,31) IOBS,ID0,IM0,S0,IDB,IMB,SB,
     &                     VSEC,VSD,isssn,NAME1,jsssn,NAME2
  31          FORMAT (' ',I6,' CZ',I4,I3,F6.2,I9,I3,F6.2,
** v 4.29i
*    &                4X,F7.2,11X,F8.1,2(1X,i5,2x,A30))
     &                4X,F7.2,11X,F8.1,2(1X,'(',i5,')',1x,A30))
**********************************     
     
            ENDIF
	  else
	    lc = .false.
          ENDIF
        ENDIF

*** CONSTRAINED GEOID HEIGHT DIFFERENCE

      ELSEIF (KIND .EQ. 15) THEN
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        V = OBS0-OBSB
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        IF (LVC) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,411) IOBS,OBS0,GMDE,RN,NAME1,NAME2
*411          FORMAT (' ',I6,' DN',F15.4,F12.4,F18.2,2(1X,A30))
              WRITE (LUNIT,411) IOBS,OBS0,GMDE,RN,
     &                          isssn,name1,jsssn,name2	      
 411          FORMAT (' ',I6,' DN',F15.4,F12.4,F18.2,
** v 4.29i
*    &                2(1x,i5,2x,a30))
     &                2(1x,'(',i5,')',1x,a30))
*************************************** 
 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,412) IOBS,OBS0,OBSB,V,SD3,VSD3,GMDE,RN,NAME2
*412          FORMAT (' ',I6,' DN',F15.4,F17.4,F19.3,F8.2,F8.2,
*    &                 F10.4,F6.2,4X,A30)
              WRITE (LUNIT,412) IOBS,OBS0,OBSB,V,SD3,VSD3,GMDE,RN,
     &                          jsssn,name2	      
 412          FORMAT (' ',I6,' DN',F15.4,F17.4,F19.3,F8.2,F8.2,
** v 4.29i
*    &                 F10.4,F6.2,t99,i5,2x,A30)
     &                 F10.4,F6.2,t99,'(',i5,')',1x,A30)
******************************
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,41) IOBS,OBS0,OBSB,V,VSD,NAME1,NAME2
*41           FORMAT (' ',I6,' DN',F14.3,F18.3,F21.3,F8.1,
*    &                2(1X,A30))
              WRITE (LUNIT,41) IOBS,OBS0,OBSB,V,VSD,
     &                         isssn,name1,jsssn,name2	      
 41           FORMAT (' ',I6,' DN',F14.3,F18.3,F21.3,F8.1,
** v 4.29i
*    &                2(1X,i5,2x,A30))
     &                2(1X,'(',i5,')',1x,A30))
**********************************
     
            ENDIF
	  else
	    lc = .false.
          ENDIF
        ENDIF

*** CONSTRAINED ORTHOMETRIC HEIGHT DIFFERENCE

      ELSEIF (KIND .EQ. 16) THEN
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        V = OBS0-OBSB
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        IF (LVC) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,451) IOBS,OBS0,GMDE,RN,NAME1,NAME2
*451          FORMAT (' ',I6,' DO',F15.4,F12.4,F18.2,2(1X,A30))
              WRITE (LUNIT,451) IOBS,OBS0,GMDE,RN,
     &                          isssn,name1,jsssn,name2	      
 451          FORMAT (' ',I6,' DO',F15.4,F12.4,F18.2,
** v 4.29i
*    &                 2(1x,i5,2x,a30)) 
     &                 2(1x,'(',i5,')',1x,a30)) 
******************************* 
 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,452) IOBS,OBS0,OBSB,V,SD3,VSD3,GMDE,RN,NAME2
*452          FORMAT (' ',I6,' DO',F15.4,F17.4,F19.3,F8.2,F8.2,
*    &                 F10.4,F6.2,4X,A30)
              WRITE (LUNIT,452) IOBS,OBS0,OBSB,V,SD3,VSD3,GMDE,RN,
     &	                             jsssn,name2 
 452          FORMAT (' ',I6,' DO',F15.4,F17.4,F19.3,F8.2,F8.2,
** v 4.29i
*    &                 F10.4,F6.2,t99,i5,2x,A30)
     &                 F10.4,F6.2,t99,'(',i5,')',1x,A30)
************************
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,45) IOBS,OBS0,OBSB,V,VSD,NAME1,NAME2
*45           FORMAT (' ',I6,' DO',F14.3,F18.3,F21.3,F8.1,
*    &                2(1X,A30))
              WRITE (LUNIT,45) IOBS,OBS0,OBSB,V,VSD,
     &                         isssn,name1,jsssn,name2	      
 45           FORMAT (' ',I6,' DO',F14.3,F18.3,F21.3,F8.1,
** v 4.29i
*    &                2(1X,i5,2x,A30))
     &                2(1X,'(',i5,')',1x,A30))
***********************************     
     
            ENDIF
	  else
	    lc = .false.
          ENDIF
        ENDIF

*** CONSTRAINED ELLIPSOIDAL HEIGHT DIFFERENCE

      ELSEIF (KIND .EQ. 17) THEN
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        V = OBS0-OBSB
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
        IF (LVC) THEN
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,491) IOBS,OBS0,GMDE,RN,NAME1,NAME2
*491          FORMAT (' ',I6,' DE',F15.4,F12.4,F18.2,2(1X,A30))
              WRITE (LUNIT,491) IOBS,OBS0,GMDE,RN,
     &                          isssn,name1,jsssn,name2	      
 491          FORMAT (' ',I6,' DE',F15.4,F12.4,F18.2,
** v 4.29i
*    &                2(1x,i5,2x,a30)) 
     &                2(1x,'(',i5,')',1x,a30)) 
***********************	    
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,492) IOBS,OBS0,OBSB,V,SD3,VSD3,GMDE,RN,NAME2 
*492          FORMAT (' ',I6,' DE',F15.4,F17.4,F19.3,F8.2,F8.2,
*    &                 F10.4,F6.2,4X,A30)
              WRITE (LUNIT,492) IOBS,OBS0,OBSB,V,SD3,VSD3,GMDE,RN, 
     &                           jsssn,name2	      
 492          FORMAT (' ',I6,' DE',F15.4,F17.4,F19.3,F8.2,F8.2,
** v 4.29i
*    &                 F10.4,F6.2,t99,i5,2x,A30)
     &                 F10.4,F6.2,t99,'(',i5,')',1x,A30)
******************************
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,49) IOBS,OBS0,OBSB,V,VSD,NAME1,NAME2
*49           FORMAT (' ',I6,' DE',F14.3,F18.3,F21.3,F8.1,
*    &                2(1X,A30))
              WRITE (LUNIT,49) IOBS,OBS0,OBSB,V,VSD,
     &                         isssn,name1,jsssn,name2	      
 49           FORMAT (' ',I6,' DE',F14.3,F18.3,F21.3,F8.1,
** v 4.29i
*    &                2(1X,i5,2x,A30))
     &                2(1X,'(',i5,')',1x,A30))
*********************************
     
            ENDIF
	  else
	    lc = .false.
          ENDIF
        ENDIF

*** CONSTRAINED NORTH COORD DIFFERENCE

      ELSEIF (KIND .EQ. 21) THEN
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        V = OBS0-OBSB
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
          IF (DABS(VSD3) .GE. CRIT) THEN
	    lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
*             WRITE (LUNIT,771) IOBS,OBS0,GMDE,RN,NAME1,NAME2
*771          FORMAT (' ',I6,' CP-N',F13.4,F11.3,F19.2,2(1X,A30))
              WRITE (LUNIT,771) IOBS,OBS0,GMDE,RN,
     &                          isssn,name1,jsssn,name2	      
 771          FORMAT (' ',I6,' CP-N',F13.4,F11.3,F19.2,
** v 4.29i
*    &                2(1x,i5,2x,a30)) 
     &                2(1x,'(',i5,')',1x,a30)) 
********************************** 
 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,772) IOBS,OBS0,OBSB,V,SD3,VSD3,GMDE,RN,NAME2
*772          FORMAT (' ',I6,' CP-N',F13.4,F17.4,F19.3,F8.2,F8.2,F10.4,
*    &                F6.2,4X,A30)
              WRITE (LUNIT,772) IOBS,OBS0,OBSB,V,SD3,VSD3,GMDE,RN,
     &                          jsssn,name2	      
 772          FORMAT (' ',I6,' CP-N',F13.4,F17.4,F19.3,F8.2,F8.2,F10.4,
** v 4.29i
*    &                F6.2,t99,i5,2x,A30)
     &                F6.2,t99,'(',i5,')',1x,A30)
*********************
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,77) IOBS,OBS0,OBSB,V,VSD,NAME1,NAME2
* 77          FORMAT (' ',I6,' CP-N',F12.3,F18.3,F21.3,
*    &                F8.1,2(1X,A30))
              WRITE (LUNIT,77) IOBS,OBS0,OBSB,V,VSD,
     &	                       isssn,name1,jsssn,name2      
  77          FORMAT (' ',I6,' CP-N',F12.3,F18.3,F21.3,
** v 4.29i
*    &                F8.1,2(1X,i5,2x,A30))
     &                F8.1,2(1X,'(',i5,')',1x,A30))
*****************************
     
            ENDIF
	  else
	    lc = .false.
          ENDIF

*** CONSTRAINED EAST COORD DIFFERENCE

      ELSEIF (KIND .EQ. 22) THEN
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        V = OBS0-OBSB
        VSD = V/SD
        VSD3 = DIVID(V,SD3)
        IF (LSS) VSD = VSD/SIGUWT
	
***** ver 4.28n
        IF (LSS) SD3 = SD3/SIGUWT
        IF (LSS) VSD3 = VSD3/SIGUWT
************
	
          IF (DABS(VSD3) .GE. CRIT) THEN
	  lc = .true.
            CALL LINE2 (1)
            IF (IMODE .EQ. 0) THEN
** v 4.29f	    
              WRITE (LUNIT,871) IOBS,OBS0,GMDE,RN,
     &                             isssn,name1,jsssn,name2	      
 871          FORMAT (' ',I6,' CP-E',F13.4,F11.3,F19.2,
** v 4.29i
*    &                 2(1x,i5,2x,a30)) 
     &                 2(1x,'(',i5,')',1x,a30)) 
***************************** 
 
            ELSEIF (IMODE .EQ. 3) THEN
** v 4.29f	    
*             WRITE (LUNIT,872) IOBS,OBS0,OBSB,V,SD3,VSD3,GMDE,RN,NAME2
*872          FORMAT (' ',I6,' CP-E',F13.4,F17.4,F19.3,F8.2,F8.2,F10.4,
*    &                F6.2,4X,A30)
              WRITE (LUNIT,872) IOBS,OBS0,OBSB,V,SD3,VSD3,GMDE,RN,
     &                jsssn,name2	      
 872          FORMAT (' ',I6,' CP-E',F13.4,F17.4,F19.3,F8.2,F8.2,F10.4,
** v 4.29i
*    &                F6.2,t99,i5,2x,A30)
     &                F6.2,t99,'(',i5,')',1x,A30)
*****************************
     
            ELSE
** v 4.29f	    
*             WRITE (LUNIT,87) IOBS,OBS0,OBSB,V,VSD,NAME1,NAME2
* 87          FORMAT (' ',I6,' CP-E',F12.3,F18.3,F21.3,
*    *                F8.1,2(1X,A30))
              WRITE (LUNIT,87) IOBS,OBS0,OBSB,V,VSD,
     &	                        isssn,name1,jsssn,name2
  87          FORMAT (' ',I6,' CP-E',F12.3,F18.3,F21.3,
** v 4.29i
*    &                F8.1,2(1X,i5,2x,A30))
     &                F8.1,2(1X,'(',i5,')',1x,A30))
********************************
     
            ENDIF
	  else
	    lc = .false.
          ENDIF

*** ILLEGAL KIND

      ELSE
        WRITE (LUNIT,667) KIND
  667   FORMAT ('0ILLEGAL KIND IN ROBS--',I5)
          write (lunit,*) 'subs2, 67'
        CALL ABORT2
      ENDIF

*** ACCUMULATE STATISTICS

      IF (LSS) VSD = VSD*SIGUWT
      
***** ver 4.28n
        IF (LSS) SD3 = SD3*SIGUWT
        IF (LSS) VSD3 = VSD3*SIGUWT
************
c New additional argument to RSTAT Mike Potterfield 3/15/07

      IF ( KIND .LE.  3  .OR.  KIND .EQ.  7  .OR.
     &     KIND .EQ. 13  .OR.  KIND .GE. 15 ) THEN
        CALL RSTAT (V,VSD,RN,VSD3,IOBS,KIND,IOBS)
      ELSE
        CALL RSTAT (VSEC,VSD,RN,VSD3,IOBS,KIND,IOBS)
      ENDIF

      RETURN
      END
      SUBROUTINE RSINIT

*** INITIALIZE RESIDUAL STATISTICS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRD = 20 )
**v 4.28
      PARAMETER ( MXRDV = 20 )
**v 4.29f
      PARAMETER ( MXRDU = 20 )
**v 4.29i
      PARAMETER ( MXRDL = 20 )
      
      COMMON /RESTAT/ VSD20(MXRD), VMAX, VMIN, VSDMAX, VSDMIN, VSUM,
     &                VSDSUM, VSD21, VSD22, VSD23, VSD24, VSD25, VSD26,
     &                VSD27, VSD28, VSD29, VSD210, VSD211, VSD212,
     &                VABS1, VABS2, VABS3, VABS4, VABS5, VABS6, VABS7,
     &                VABS8, VABS9, VABS10, VABS11, VABS12, I20(MXRD),
     &                N0, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10,
     &                N11, N12, NI20
**4.24
      COMMON /RESTA5/ VDXMAX, VDYMAX, VDZMAX, VDNMAX, VDEMAX, VDUMAX,
     &                VDXMIN, VDYMIN, VDZMIN, VDNMIN, VDEMIN, VDUMIN,
     &                IDXMAX, IDYMAX, IDZMAX, IDNMAX, IDEMAX, IDUMAX,
     &                IDXMIN, IDYMIN, IDZMIN, IDNMIN, IDEMIN, IDUMIN
****
**v 4.29f
*     COMMON /RESTA6/ VDLMAX, IDLMAX
************
** v 4.29g
      COMMON /RESTA6/ VDLMAX, VDLMIN, IDLMAX, IDLMIN
************

      COMMON /REDUN/  RN1, RN2, RN3, RN4, RN5, RN6, RN7, RN8, RN9,
     &                RN10, RN11, RN12
      COMMON /RESTA4/ VN2, VE2, VU2, VNABS, VEABS, VUABS, NNN, NNE, NNU
**v 4.28
      COMMON /RESTA3/ V20(MXRDV),I20V(MXRDV),NI20V  
**v 4.29f
      COMMON /RESTA7/ VDU20(MXRDU),I20DU(MXRDU),NI20DU  
**v 4.29i
      COMMON /RESTA8/ VDL20(MXRDL),I20DL(MXRDL),NI20DL  


      DO 1 I = 1,MXRD
        VSD20(I) = 0.D0
        I20(I) = 0
    1 CONTINUE
      NI20 = 0

    
**v 4.28
      DO 2 I = 1,MXRDV
        V20(I) = 0.D0
        I20V(I) = 0
    2 CONTINUE
      NI20V = 0
***********
**v 4.29f
      DO 3 I = 1,MXRDU
        VDU20(I) = 0.D0
        I20DU(I) = 0
    3 CONTINUE
      NI20DU = 0
***********
**v 4.29i
      DO 4 I = 1,MXRDL
        VDL20(I) = 0.D0
        I20DL(I) = 0
    4 CONTINUE
      NI20DL = 0
************

      N0  = 0
      N1  = 0
      N2  = 0
      N3  = 0
      N4  = 0
      N5  = 0
      N6  = 0
      N7  = 0
      N8  = 0
      N9  = 0
      N10 = 0
      N11 = 0
      N12 = 0

      RN1  = 0.D0
      RN2  = 0.D0
      RN3  = 0.D0
      RN4  = 0.D0
      RN5  = 0.D0
      RN6  = 0.D0
      RN7  = 0.D0
      RN8  = 0.D0
      RN9  = 0.D0
      RN10 = 0.D0
      RN11 = 0.D0
      RN12 = 0.D0

      VMAX   = -1.D50
      VMIN   =  1.D50
      VSDMAX = -1.D50
      VSDMIN =  1.D50
      VSUM   =  0.D0
      VSDSUM =  0.D0
**4.24
      VDXMAX   = -1.D50
      VDYMAX   = -1.D50
      VDZMAX   = -1.D50
      VDNMAX   = -1.D50
      VDEMAX   = -1.D50
      VDUMAX   = -1.D50
** v 4.29f
      VDLMAX   = -1.D50
**********
** v 4.29g
      VDLMIN   =  1.D50
**********
      
      VDXMIN   =  1.D50
      VDYMIN   =  1.D50
      VDZMIN   =  1.D50
      VDNMIN   =  1.D50
      VDEMIN   =  1.D50
      VDUMIN   =  1.D50

      IDXMAX   = 0
      IDYMAX   = 0
      IDZMAX   = 0
      IDNMAX   = 0
      IDEMAX   = 0
      IDUMAX   = 0
** v 4.29f
      IDLMAX   = 0
**************      
** v 4.29g
      IDLMIN   = 0
**************      
      IDXMIN   = 0
      IDYMIN   = 0
      IDZMIN   = 0
      IDNMIN   = 0
      IDEMIN   = 0
      IDUMIN   = 0
*****
      VSD21  = 0.D0
      VSD22  = 0.D0
      VSD23  = 0.D0
      VSD24  = 0.D0
      VSD25  = 0.D0
      VSD26  = 0.D0
      VSD27  = 0.D0
      VSD28  = 0.D0
      VSD29  = 0.D0
      VSD210 = 0.D0
      VSD211 = 0.D0
      VSD212 = 0.D0

      VABS1  = 0.D0
      VABS2  = 0.D0
      VABS3  = 0.D0
      VABS4  = 0.D0
      VABS5  = 0.D0
      VABS6  = 0.D0
      VABS7  = 0.D0
      VABS8  = 0.D0
      VABS9  = 0.D0
      VABS10 = 0.D0
      VABS11 = 0.D0
      VABS12 = 0.D0

      VN2   = 0.D0
      VE2   = 0.D0
      VU2   = 0.D0
      VNABS = 0.D0
      VEABS = 0.D0
      VUABS = 0.D0
      NNN   = 0
      NNE   = 0
      NNU   = 0

      RETURN
      END
      SUBROUTINE RSOUT (SIGUWT)

*** LIST RESIDUAL STATISTICS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRD = 20 )
      
**v 4.28
      PARAMETER ( MXRDV = 20 )
      
**v 4.29f
      PARAMETER ( MXRDU = 20 )
      
**v 4.29i
      PARAMETER ( MXRDL = 20 )
      
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
**4.24
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
******
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /RESTAT/ VSD20(MXRD), VMAX, VMIN, VSDMAX, VSDMIN, VSUM,
     &                VSDSUM, VSD21, VSD22, VSD23, VSD24, VSD25, VSD26,
     &                VSD27, VSD28, VSD29, VSD210, VSD211, VSD212,
     &                VABS1, VABS2, VABS3, VABS4, VABS5, VABS6, VABS7,
     &                VABS8, VABS9, VABS10, VABS11, VABS12, I20(MXRD),
     &                N0, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10,
     &                N11, N12, NI20
**4.24**
      COMMON /RESTA5/ VDXMAX, VDYMAX, VDZMAX, VDNMAX, VDEMAX, VDUMAX,
     &                VDXMIN, VDYMIN, VDZMIN, VDNMIN, VDEMIN, VDUMIN,
     &                IDXMAX, IDYMAX, IDZMAX, IDNMAX, IDEMAX, IDUMAX,
     &                IDXMIN, IDYMIN, IDZMIN, IDNMIN, IDEMIN, IDUMIN
********
** v 4.29f
*     COMMON /RESTA6/ VDLMAX, IDLMAX
************
** v 4.29g
      COMMON /RESTA6/ VDLMAX, VDLMIN, IDLMAX, IDLMIN
************

      COMMON /REDUN/  RN1, RN2, RN3, RN4, RN5, RN6, RN7, RN8, RN9,
     &                RN10, RN11, RN12
      COMMON /CONST/  PI, PI2, RAD
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON /RESTA4/ VN2, VE2, VU2, VNABS, VEABS, VUABS, NNN, NNE, NNU
      COMMON/UNITS/ LUNIT
**v 4.28
      COMMON /RESTA3/ V20(MXRDV),I20V(MXRDV),NI20V  

**v 4.29f
      COMMON /RESTA7/ VDU20(MXRDU),I20DU(MXRDU),NI20DU  
      
**v 4.29i
      COMMON /RESTA8/ VDL20(MXRDL),I20DL(MXRDL),NI20DL  
      
*** CORRECTION TO GPS REDUNDANCY NUMBERS IF REJECTIONS EXIST

      IF (NGPSR .GT. 0) THEN
        RN1 = RN1 - NGPSR
        RN2 = RN2 - NGPSR
        RN3 = RN3 - NGPSR
      ENDIF

*** COMPLETE COMPUTATION OF STATS

      N = N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8 + N9 + N10 + N11 + N12
      RN = RN1 + RN2 + RN3 + RN4 + RN5 + RN6 + RN7 + RN8 + RN9 +
     &     RN10 + RN11 + RN12
      VSD2 = VSD21 + VSD22 + VSD23 + VSD24 + VSD25 + VSD26 + VSD27 +
     &       VSD28 + VSD29 + VSD210 + VSD211 + VSD212
      VMEAN = DIVIDE(VSUM, N)
      VSDMN = DIVIDE(VSDSUM, N)

      RMSV   = DSQRT( DIVIDE( VSD2,   N)  )
      RMSV1  = DSQRT( DIVIDE( VSD21, N1)  )
      RMSV2  = DSQRT( DIVIDE( VSD22, N2)  )
      RMSV3  = DSQRT( DIVIDE( VSD23, N3)  )
      RMSV4  = DSQRT( DIVIDE( VSD24, N4)  )
      RMSV5  = DSQRT( DIVIDE( VSD25, N5)  )
      RMSV6  = DSQRT( DIVIDE( VSD26, N6)  )
      RMSV7  = DSQRT( DIVIDE( VSD27, N7)  )
      RMSV8  = DSQRT( DIVIDE( VSD28, N8)  )
      RMSV9  = DSQRT( DIVIDE( VSD29, N9)  )
      RMSV10 = DSQRT( DIVIDE(VSD210, N10) )
      RMSV11 = DSQRT( DIVIDE(VSD211, N11) )
      RMSV12 = DSQRT( DIVIDE(VSD212, N12) )
      RMSVN  = DSQRT( DIVIDE(   VN2, NNN) )
      RMSVE  = DSQRT( DIVIDE(   VE2, NNE) )
      RMSVU  = DSQRT( DIVIDE(   VU2, NNU) )

      ABSV1  = DIVIDE( VABS1, N1)
      ABSV2  = DIVIDE( VABS2, N2)
      ABSV3  = DIVIDE( VABS3, N3)
      ABSV4  = DIVIDE( VABS4, N4)
      ABSV5  = DIVIDE( VABS5, N5)
      ABSV6  = DIVIDE( VABS6, N6)
      ABSV7  = DIVIDE( VABS7, N7)
      ABSV8  = DIVIDE( VABS8, N8)
      ABSV9  = DIVIDE( VABS9, N9)
      ABSV10 = DIVIDE(VABS10, N10)
      ABSV11 = DIVIDE(VABS11, N11)
      ABSV12 = DIVIDE(VABS12, N12)
      ABSVN  = DIVIDE( VNABS, NNN)
      ABSVE  = DIVIDE( VEABS, NNE)
      ABSVU  = DIVIDE( VUABS, NNU)

*** HEADING

      CALL HEAD
      CALL LINE (3)
      WRITE (LUNIT,1)
    1 FORMAT ('0RESIDUAL STATISTICS', /)

*** MAXIMUM RESIDUALS

**v 4.28

      IF (IMODE .NE. 0) THEN
        CALL LINE (2)
        IF (IMODE .NE. 3) THEN
          WRITE (LUNIT,52) NI20V
   52     FORMAT (' OBSERVATION NUMBERS OF', I3, ' GREATEST',
     &            ' RESIDUALS (V)', /)
        ELSE
          WRITE (LUNIT,51) NI20V
   51     FORMAT (' OBSERVATION NUMBERS OF', I3, ' GREATEST',
     &            ' RESIDUALS (V)', /)
        ENDIF
        CALL LINE (2)
        WRITE (LUNIT,53) (I20V(I), I = 1, NI20V)
   53   FORMAT (1X, 20I6, /)
      ENDIF

*********************************

      IF (IMODE .NE. 0) THEN
        CALL LINE (2)
        IF (IMODE .NE. 3) THEN
          WRITE (LUNIT,2) NI20
    2     FORMAT (' OBSERVATION NUMBERS OF', I3, ' GREATEST',
     &            ' QUASI-NORMALIZED RESIDUALS (V/SDV)', /)
        ELSE
          WRITE (LUNIT,21) NI20
   21     FORMAT (' OBSERVATION NUMBERS OF', I3, ' GREATEST',
     &            ' STANDARDIZED RESIDUALS (V/SDV)', /)
        ENDIF
        CALL LINE (2)
        WRITE (LUNIT,3) (I20(I), I = 1, NI20)
    3   FORMAT (1X, 20I6, /)


** v 4.29f
      if(NI20DU .gt. 0)  then

      IF (IMODE .NE. 0) THEN
        CALL LINE (2)
        IF (IMODE .NE. 3) THEN
          WRITE (LUNIT,152) NI20DU
  152     FORMAT (' OBSERVATION NUMBERS OF', I3, ' GREATEST',
     &            ' GPS DU COMPONENT RESIDUALS (V)', /)
        ELSE
          WRITE (LUNIT,151) NI20DU
  151     FORMAT (' OBSERVATION NUMBERS OF', I3, ' GREATEST',
     &            ' GPS DU COMPONENT RESIDUALS (V)', /)
        ENDIF
        CALL LINE (2)
        WRITE (LUNIT,153) (I20DU(I), I = 1, NI20DU)
  153   FORMAT (1X, 20I6, /)
      ENDIF
      endif

********************************

** v 4.29i
      if(NI20DL .gt. 0)  then

      IF (IMODE .NE. 0) THEN
        CALL LINE (2)
        IF (IMODE .NE. 3) THEN
          WRITE (LUNIT,162) NI20DU
  162     FORMAT (' OBSERVATION NUMBERS OF', I3, ' GREATEST',
     &            ' GPS DL COMPONENT RESIDUALS (V)', /)
        ELSE
          WRITE (LUNIT,161) NI20DU
  161     FORMAT (' OBSERVATION NUMBERS OF', I3, ' GREATEST',
     &            ' GPS DL COMPONENT RESIDUALS (V)', /)
        ENDIF
        CALL LINE (2)
        WRITE (LUNIT,163) (I20DL(I), I = 1, NI20DL)
  163   FORMAT (1X, 20I6, /)
      ENDIF
      endif

***************************
*** EXTREMA

        IF (LSS) THEN
          VSDMAX = VSDMAX/SIGUWT
          VSDMIN = VSDMIN/SIGUWT
          VSDMN = VSDMN/SIGUWT
        ENDIF
        CALL LINE (6)
        WRITE (LUNIT,4) N, N0, VMAX, VSDMAX, VMIN, VSDMIN, VMEAN, VSDMN 
    4   FORMAT ('  TOTAL=', I6, 8X, 'NO-CHECK=', I5, /,
     &          '  MAX V=', 1PD9.1, 4X, 'MAX V/SDV=', 0PF9.3, /,
     &          '  MIN V=', 1PD9.1, 4X, 'MIN V/SDV=', 0PF9.3, /,
     &          ' MEAN V=', 1PD9.1, 3X, 'MEAN V/SDV=', 0PF9.3, /)

**4.24***
        IF( .NOT. LGPS) THEN
        CALL LINE (12)
*       WRITE (LUNIT,14) VDXMAX,IDXMAX,VDXMIN,IDXMIN, 
*    &                   VDYMAX,IDYMAX,VDYMIN,IDYMIN, 
*    &                   VDZMAX,IDZMAX,VDZMIN,IDZMIN, 
*    &                   VDNMAX,IDNMAX,VDNMIN,IDNMIN, 
*    &                   VDEMAX,IDEMAX,VDEMIN,IDEMIN, 
*    &                   VDUMAX,IDUMAX,VDUMIN,IDUMIN  

*  14   FORMAT ('           MAX V     OBS #      MIN V     OBS #',/,       
*    &       '  DX ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/,      
*    &       '  DY ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/,      
*    &       '  DZ ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/,      
*    &       '  DN ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/,      
*    &       '  DE ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/,      
*    &       '  DU ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/)      

*********
** v 4.29f/g
        WRITE (LUNIT,14) VDXMAX,IDXMAX,VDXMIN,IDXMIN, 
     &                   VDYMAX,IDYMAX,VDYMIN,IDYMIN, 
     &                   VDZMAX,IDZMAX,VDZMIN,IDZMIN, 
     &                   VDNMAX,IDNMAX,VDNMIN,IDNMIN, 
     &                   VDEMAX,IDEMAX,VDEMIN,IDEMIN, 
     &                   VDLMAX,IDLMAX,VDLMIN,IDLMIN,
     &                   VDUMAX,IDUMAX,VDUMIN,IDUMIN  

   14   FORMAT ('           MAX V     OBS #      MIN V     OBS #',/,       
     &       '  DX ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/,      
     &       '  DY ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/,      
     &       '  DZ ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/,      
     &       '  DN ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/,      
     &       '  DE ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/,      
     &       '  DL ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/,
     &       '  DU ',2X,1PD9.1,2X,I8,2X,1PD9.1,2X,I8,/)      
     
********************

        ENDIF
*******************
*** STATS

        IF (IMODE .NE. 3) THEN
          CALL LINE (11)
          WRITE (LUNIT,5) N1,  VSD21,  RMSV1,  ABSV1,
     &                N2,  VSD22,  RMSV2,  ABSV2,
     &                N3,  VSD23,  RMSV3,  ABSV3,
     &                N10, VSD210, RMSV10, ABSV10,
     &                N11, VSD211, RMSV11, ABSV11,
     &                N12, VSD212, RMSV12, ABSV12,
     &                N9,  VSD29,  RMSV9,  ABSV9,
     &                N8,  VSD28,  RMSV8,  ABSV8,
     &                N7,  VSD27,  RMSV7,  ABSV7,
     &                N5,  VSD25,  RMSV5,  ABSV5,
     &                N6,  VSD26,  RMSV6,  ABSV6,
     &                N4,  VSD24,  RMSV4,  ABSV4,
     &                N,   VSD2,   RMSV
    5     FORMAT (15X, 'N', 6X, 'VTPV', 6X, 'RMS', 3X, 'MEAN ABS', /,
     &            31X, 'VTPV', 3X, 'RESIDUAL', /,
     &            ' DELTA X  ', I6, F10.1, F9.2, F11.3, ' (METERS)', /,
     &            ' DELTA Y  ', I6, F10.1, F9.2, F11.3, ' (METERS)', /,
     &            ' DELTA Z  ', I6, F10.1, F9.2, F11.3, ' (METERS)', /,
     &            ' DOPPLER X', I6, F10.1, F9.2, F11.3, ' (METERS)', /,
     &            ' DOPPLER Y', I6, F10.1, F9.2, F11.3, ' (METERS)', /,
     &            ' DOPPLER Z', I6, F10.1, F9.2, F11.3, ' (METERS)', /,
     &            ' DIRECTION', I6, F10.1, F9.2, F11.3, ' (SECONDS)', /,
     &            ' H ANGLE  ', I6, F10.1, F9.2, F11.3, ' (SECONDS)', /,
     &            ' ZEN DIST ', I6, F10.1, F9.2, F11.3, ' (SECONDS)', /,
     &            ' DISTANCE ', I6, F10.1, F9.2, F11.3, ' (METERS)', /,
     &            ' AZIMUTH  ', I6, F10.1, F9.2, F11.3, ' (SECONDS)', /,
     &            ' OTHER    ', I6, F10.1, F9.2, F11.3, /,
     &            ' TOTAL    ', I6, F10.1, F9.2, /)

*** GPS STATS IN LOCAL HORIZON

          IF ( NNN+NNE+NNU .GT. 0) THEN
            CALL LINE (6)
            WRITE (LUNIT,25) NNN, RMSVN, ABSVN,
     &                   NNE, RMSVE, ABSVE,
     &                   NNU, RMSVU, ABSVU
   25       FORMAT (15X, 'N', 10X, 'RMS', 9X, 'MEAN ABS', /,
     &              10X, 'CONTRIB.', 6X, 'RESIDUAL', 6X, 'RESIDUAL', /,
     &             ' NORTH    ', I6, F15.3, 4X, F11.3, ' (METERS)', /,
     &             ' EAST     ', I6, F15.3, 4X, F11.3, ' (METERS)', /,
     &             ' UP       ', I6, F15.3, 4X, F11.3, ' (METERS)', /)
          ENDIF

        ELSE

          VRN1  = DIVID(VSD21,  RN1)
          VRN2  = DIVID(VSD22,  RN2)
          VRN3  = DIVID(VSD23,  RN3)
          VRN9  = DIVID(VSD29,  RN9)
          VRN8  = DIVID(VSD28,  RN8)
          VRN7  = DIVID(VSD27,  RN7)
          VRN5  = DIVID(VSD25,  RN5)
          VRN6  = DIVID(VSD26,  RN6)
          VRN4  = DIVID(VSD24,  RN4)
          VRN10 = DIVID(VSD210, RN10)
          VRN11 = DIVID(VSD211, RN11)
          VRN12 = DIVID(VSD212, RN12)
          VRN   = DIVID(VSD2,   RN)

          CALL LINE (11)
          WRITE (LUNIT,15) N1,  VSD21,  RMSV1,  RN1,  VRN1,  ABSV1,
     &                 N2,  VSD22,  RMSV2,  RN2,  VRN2,  ABSV2,
     &                 N3,  VSD23,  RMSV3,  RN3,  VRN3,  ABSV3,
     &                 N10, VSD210, RMSV10, RN10, VRN10, ABSV10,
     &                 N11, VSD211, RMSV11, RN11, VRN11, ABSV11,
     &                 N12, VSD212, RMSV12, RN12, VRN12, ABSV12,
     &                 N9,  VSD29,  RMSV9,  RN9,  VRN9,  ABSV9,
     &                 N8,  VSD28,  RMSV8,  RN8,  VRN8,  ABSV8,
     &                 N7,  VSD27,  RMSV7,  RN7,  VRN7,  ABSV7,
     &                 N5,  VSD25,  RMSV5,  RN5,  VRN5,  ABSV5,
     &                 N6,  VSD26,  RMSV6,  RN6,  VRN6,  ABSV6,
     &                 N4,  VSD24,  RMSV4,  RN4,  VRN4,  ABSV4,
     &                 N,   VSD2,   RMSV,   RN,   VRN
   15     FORMAT (15X, 'N', 6X, 'VTPV', 6X, 'RMS', 7X, 'RN', 8X,
     &            'VTPV/RN', 6X, 'MEAN ABS', /,
     &            31X, 'VTPV', 30X, 'RESIDUAL', /,
     &    ' DELTA X  ', I6, F10.1, F9.2,F9.2,F15.2,F11.3,' (METERS)',/,
     &    ' DELTA Y  ', I6, F10.1, F9.2,F9.2,F15.2,F11.3,' (METERS)',/,
     &    ' DELTA Z  ', I6, F10.1, F9.2,F9.2,F15.2,F11.3,' (METERS)',/,
     &    ' DOPPLER X', I6, F10.1, F9.2,F9.2,F15.2,F11.3,' (METERS)',/,
     &    ' DOPPLER Y', I6, F10.1, F9.2,F9.2,F15.2,F11.3,' (METERS)',/,
     &    ' DOPPLER Z', I6, F10.1, F9.2,F9.2,F15.2,F11.3,' (METERS)',/,
     &    ' DIRECTION', I6, F10.1, F9.2,F9.2,F15.2,F11.3,' (SECONDS)',/,
     &    ' H ANGLE  ', I6, F10.1, F9.2,F9.2,F15.2,F11.3,' (SECONDS)',/,
     &    ' ZEN DIST ', I6, F10.1, F9.2,F9.2,F15.2,F11.3,' (SECONDS)',/,
     &    ' DISTANCE ', I6, F10.1, F9.2,F9.2,F15.2,F11.3,' (METERS)',/,
     &    ' AZIMUTH  ', I6, F10.1, F9.2,F9.2,F15.2,F11.3,' (SECONDS)',/,
     &    ' OTHER    ', I6, F10.1, F9.2,F9.2,F15.2,F11.3, /,
     &    ' TOTAL    ', I6, F10.1, F9.2,F9.2,F15.2, /)

*** GPS STATS IN LOCAL HORIZON

          IF ( NNN+NNE+NNU .GT. 0) THEN
            CALL LINE (6)
            WRITE (LUNIT,26) NNN, RMSVN, ABSVN,
     &                   NNE, RMSVE, ABSVE,
     &                   NNU, RMSVU, ABSVU
   26       FORMAT (15X, 'N', 10X, 'RMS', 36X, 'MEAN ABS', /,
     &              10X, 'CONTRIB.', 6X, 'RESIDUAL', 33X, 'RESIDUAL', /,
     &             ' NORTH    ', I6, F15.3, 28X, F11.3, ' (METERS)', /,
     &             ' EAST     ', I6, F15.3, 28X, F11.3, ' (METERS)', /,
     &             ' UP       ', I6, F15.3, 28X, F11.3, ' (METERS)', /)
          ENDIF
        ENDIF

*** COMPUTE VARIANCES

        IDOF   = NOBS + NCON - NUNK
        SUMPVV = VSD2
        VARUWT = DIVIDE(SUMPVV, IDOF)
        SIGUWT = DSQRT(VARUWT)

*** PRINT VARIANCE OF UNIT WEIGHT

        WRITE (LUNIT,6) IDOF, SUMPVV, SIGUWT, VARUWT
    6   FORMAT ('0DEGREES OF FREEDOM      =', I17, /,
     &          ' VARIANCE SUM            =', F19.1, /,
** v 4.29k
*     &          ' STD.DEV.OF UNIT WEIGHT  =', F20.2, /,
     &          ' STD.DEV.OF UNIT WEIGHT  =', F21.3, /,
     &          ' VARIANCE OF UNIT WEIGHT =', F20.2)
        IF (SIGUWT .LT. 0.01D0  .AND.  LSS) THEN
          LSS = .FALSE.
          CALL LINE (2)
          WRITE (LUNIT,23)
   23     FORMAT ('0***** WILL NOT SCALE BY THE SIGMA OF UNIT WEIGHT',
     &            ' DUE TO ITS LOW VALUE *****')
        ENDIF
      ELSE

*** EXTREMA

        CALL LINE (2)
        WRITE (LUNIT,13) N, N0
   13   FORMAT ('     TOTAL=', I7, T31, 'NO-CHECK=', I5)

*** STATS

        CALL LINE (11)
        WRITE (LUNIT,11) N1,  RN1,
     &               N2,  RN2,
     &               N3,  RN3,
     &               N10, RN10,
     &               N11, RN11,
     &               N12, RN12,
     &               N9,  RN9,
     &               N8,  RN8,
     &               N7,  RN7,
     &               N5,  RN5,
     &               N6,  RN6,
     &               N4,  RN4,
     &               N,   RN
   11   FORMAT (T16, 'N', T26, 'RN', /,
     &          ' DELTA X  ', I6, F10.2, /,
     &          ' DELTA Y  ', I6, F10.2, /,
     &          ' DELTA Z  ', I6, F10.2, /,
     &          ' DOPPLER X', I6, F10.2, /,
     &          ' DOPPLER Y', I6, F10.2, /,
     &          ' DOPPLER Z', I6, F10.2, /,
     &          ' DIRECTION', I6, F10.2, /,
     &          ' H ANGLE  ', I6, F10.2, /,
     &          ' ZEN DIST ', I6, F10.2, /,
     &          ' DISTANCE ', I6, F10.2, /,
     &          ' AZIMUTH  ', I6, F10.2, /,
     &          ' OTHER    ', I6, F10.2, /,
     &          ' TOTAL    ', I6, F10.2, /)

*** COMPUTE THE D.OF F. AND PRINTOUT RESULT

        IDOF = NOBS + NCON - NUNK
        WRITE (LUNIT,12) IDOF
   12   FORMAT ('0DEGREES OF FREEDOM=', I15)
      ENDIF

      RETURN
      END
c New additional argument to RSTAT Mike Potterfield 3/15/07
c This new argument gives the correct observation numbers for
c the residual summaries.
      SUBROUTINE RSTAT (V,VSD,RN,VSD3,IOBS,KIND,IOBSNO)

*** ACCUMULATE RESIDUAL STATISTICS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRD = 20 )
**v4.28
      PARAMETER ( MXRDV = 20 )
      
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /RESTAT/ VSD20(MXRD), VMAX, VMIN, VSDMAX, VSDMIN, VSUM,
     &                VSDSUM, VSD21, VSD22, VSD23, VSD24, VSD25, VSD26,
     &                VSD27, VSD28, VSD29, VSD210, VSD211, VSD212,
     &                VABS1, VABS2, VABS3, VABS4, VABS5, VABS6, VABS7,
     &                VABS8, VABS9, VABS10, VABS11, VABS12, I20(MXRD),
     &                N0, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10,
     &                N11, N12, NI20
      COMMON /REDUN/  RN1, RN2, RN3, RN4, RN5, RN6, RN7, RN8, RN9,
     &                RN10, RN11, RN12
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
**v 4.28
      COMMON /RESTA3/ V20(MXRDV), I20V(MXRDV), NI20V
      COMMON /UNITS/ LUNIT
******************

**v 4.28
         SAVE IOBSL
*************
 
*** ACCUMULATE NO-CHECK

      IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) THEN
        IF ( RN .LT. 1.0D-5  .AND.
     &      ( (KIND .GE.  0  .AND.  KIND .LT.   4)  .OR.
     &        (KIND .GT.  6  .AND.  KIND .LT.  18)  .OR.
     *        (KIND .GT. 20  .AND.  KIND .LE. 999) )  ) N0 = N0+1
      ELSE
        IF (V .EQ. 0.D0) N0 = N0+1
      ENDIF

*** ACCUMULATE EXTREMA

** v 4.28

      IF(IOBS.NE.IOBSL) THEN
      
c      CALL RSTAT5 (KIND, V, IOBS) 
c Mike Potterfield 10/26/08
c In order to make the observation component number match the rest
c of the adjustment output, the third parameter in the call to 
c RSTAT5 has been modified to be IOBSNO instead of IOBS.
      CALL RSTAT5 (KIND, V, IOBSNO) 
c End of changes 10/26/08

**************

      IF (V .GT. VMAX) VMAX = V
      IF (V .LT. VMIN) VMIN = V
      IF (IMODE .NE. 3) THEN
        IF (VSD .GT. VSDMAX) VSDMAX = VSD
        IF (VSD .LT. VSDMIN) VSDMIN = VSD
      ELSE
        IF (VSD3 .GT. VSDMAX) VSDMAX = VSD3
        IF (VSD3 .LT. VSDMIN) VSDMIN = VSD3
      ENDIF
      ENDIF
*******************************

*** ACCUMULATE SUM

      VABS = DABS(V)
      VSD2 = VSD*VSD
      VSUM = VSUM+V
      VSDSUM = VSDSUM+VSD

*** ACCUMULATE BY KIND

      IF (KIND .EQ. 4) THEN
        N1 = N1+1
        VABS1 = VABS1+VABS
        VSD21 = VSD21+VSD2
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) RN1 = RN1+RN
      ELSEIF (KIND .EQ. 5) THEN
        N2 = N2+1
        VABS2 = VABS2+VABS
        VSD22 = VSD22+VSD2
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) RN2 = RN2+RN
      ELSEIF (KIND .EQ. 6) THEN
        N3 = N3+1
        VABS3 = VABS3+VABS
        VSD23 = VSD23+VSD2
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) RN3 = RN3+RN
      ELSEIF (KIND .EQ. 7) THEN
        N5 = N5+1
        VABS5 = VABS5+VABS
        VSD25 = VSD25+VSD2
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) RN5 = RN5+RN
      ELSEIF (KIND .EQ. 8) THEN
        N6 = N6+1
        VABS6 = VABS6+VABS
        VSD26 = VSD26+VSD2
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) RN6 = RN6+RN
      ELSEIF (KIND .EQ. 9) THEN
        N7 = N7+1
        VABS7 = VABS7+VABS
        VSD27 = VSD27+VSD2
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) RN7 = RN7+RN
      ELSEIF (KIND .EQ. 10) THEN
        N8 = N8+1
        VABS8 = VABS8+VABS
        VSD28 = VSD28+VSD2
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) RN8 = RN8+RN
      ELSEIF (KIND .EQ. 11) THEN
        N9 = N9+1
        VABS9 = VABS9+VABS
        VSD29 = VSD29+VSD2
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) RN9 = RN9+RN
      ELSEIF (KIND .EQ. 18) THEN
        N10 = N10+1
        VABS10 = VABS10+VABS
        VSD210 = VSD210+VSD2
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) RN10 = RN10+RN
      ELSEIF (KIND .EQ. 19) THEN
        N11 = N11+1
        VABS11 = VABS11+VABS
        VSD211 = VSD211+VSD2
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) RN11 = RN11+RN
      ELSEIF (KIND .EQ. 20) THEN
        N12 = N12+1
        VABS12 = VABS12+VABS
        VSD212 = VSD212+VSD2
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) RN12 = RN12+RN
      ELSE
        N4 = N4+1
        VABS4 = VABS4+VABS
        VSD24 = VSD24+VSD2
        IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) RN4 = RN4+RN
      ENDIF

*** ACCUMULATE MXRD (20) LARGEST

**v 4.28
      IF(IOBS.NE.IOBSL) THEN

      IF (IMODE .NE. 3) THEN
        CALL ACUM20 (IOBS,VSD,IOBSNO)
      ELSE
        CALL ACUM20 (IOBS,VSD3,IOBSNO)
      ENDIF
      
**v 4.28
*** ACCUMULATE MXRDV (20) LARGEST

      IF (IMODE .NE. 3) THEN
        CALL ACUM20V (IOBS,V,IOBSNO)
      ELSE
        CALL ACUM20V (IOBS,V,IOBSNO)
      ENDIF
      ENDIF
*************************************
      IOBSL = IOBS

      RETURN
      END
      SUBROUTINE RSTAT2 (RNSUM)

*** ACCUMULATE GPS NO-CHECK FOR MODE = 0 & MODE=3

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRD = 20 )
      COMMON /RESTAT/ VSD20(MXRD), VMAX, VMIN, VSDMAX, VSDMIN, VSUM,
     &                VSDSUM, VSD21, VSD22, VSD23, VSD24, VSD25, VSD26,
     &                VSD27, VSD28, VSD29, VSD210, VSD211, VSD212,
     &                VABS1, VABS2, VABS3, VABS4, VABS5, VABS6, VABS7,
     &                VABS8, VABS9, VABS10, VABS11, VABS12, I20(MXRD),
     &                N0, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10,
     &                N11, N12, NI20

      IF (RNSUM .LT. 3.0D-5) N0 = N0 + 3

      RETURN
      END
      SUBROUTINE RSTAT3 (RNSUM)

*** ACCUMULATE DOPPLER NO-CHECK FOR MODE = 0 & MODE=3

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRD = 20 )
      COMMON /RESTAT/ VSD20(MXRD), VMAX, VMIN, VSDMAX, VSDMIN, VSUM,
     &                VSDSUM, VSD21, VSD22, VSD23, VSD24, VSD25, VSD26,
     &                VSD27, VSD28, VSD29, VSD210, VSD211, VSD212,
     &                VABS1, VABS2, VABS3, VABS4, VABS5, VABS6, VABS7,
     &                VABS8, VABS9, VABS10, VABS11, VABS12, I20(MXRD),
     &                N0, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10,
     &                N11, N12, NI20

      IF ( RNSUM .LT. 3.0D-5 ) N0 = N0 + 3

      RETURN
      END
      SUBROUTINE RSTAT4 (VN, VE, VU)

*** ACCUMULATE LGH GPS STATS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /RESTA4/ VN2, VE2, VU2, VNABS, VEABS, VUABS, NNN, NNE, NNU

      NNN = NNN + 1
      VN2 = VN2 + VN*VN
      VNABS = VNABS + DABS(VN)

      NNE = NNE + 1
      VE2 = VE2 + VE*VE
      VEABS = VEABS + DABS(VE)

      NNU = NNU + 1
      VU2 = VU2 + VU*VU
      VUABS = VUABS + DABS(VU)

      RETURN
      END

**4.28**************

      SUBROUTINE RSTAT5 (KIND,VXX,IXX)

*** ACCUMULATE MAX DX, DY, DZ, DN, DE, DU RESIDUALS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      COMMON /RESTA5/ VDXMAX, VDYMAX, VDZMAX, VDNMAX, VDEMAX, VDUMAX,
     &                VDXMIN, VDYMIN, VDZMIN, VDNMIN, VDEMIN, VDUMIN,
     &                IDXMAX, IDYMAX, IDZMAX, IDNMAX, IDEMAX, IDUMAX,
     &                IDXMIN, IDYMIN, IDZMIN, IDNMIN, IDEMIN, IDUMIN
     
** v 4.29f
*     COMMON /RESTA6/ VDLMAX, IDLMAX
*********
** v 4.29g
      COMMON /RESTA6/ VDLMAX, VDLMIN, IDLMAX, IDLMIN
*********
     
      COMMON /UNITS/ LUNIT

      IF(KIND .EQ. 4) THEN
         IF(VDXMAX .LT. VXX) THEN
            VDXMAX = VXX
            IDXMAX = IXX
         ENDIF
         IF(VDXMIN .GT. VXX) THEN
            VDXMIN = VXX
            IDXMIN = IXX
         ENDIF
      ELSEIF(KIND .EQ. 5) THEN
         IF(VDYMAX .LT. VXX) THEN
            VDYMAX = VXX
            IDYMAX = IXX
         ENDIF
         IF(VDYMIN .GT. VXX) THEN
            VDYMIN = VXX
            IDYMIN = IXX
         ENDIF
      ELSEIF(KIND .EQ. 6) THEN
         IF(VDZMAX .LT. VXX) THEN
            VDZMAX = VXX
            IDZMAX = IXX
         ENDIF
         IF(VDZMIN .GT. VXX) THEN
            VDZMIN = VXX
            IDZMIN = IXX
         ENDIF
** v 4.29h	 
      ELSEIF(KIND .EQ. 90) THEN
         IF(VDNMAX .LT. VXX) THEN
            VDNMAX = VXX
            IDNMAX = IXX
         ENDIF
         IF(VDNMIN .GT. VXX) THEN
            VDNMIN = VXX
            IDNMIN = IXX
         ENDIF
      ELSEIF(KIND .EQ. 91) THEN
         IF(VDEMAX .LT. VXX) THEN
            VDEMAX = VXX
            IDEMAX = IXX
         ENDIF
         IF(VDEMIN .GT. VXX) THEN
            VDEMIN = VXX
            IDEMIN = IXX
         ENDIF
      ELSEIF(KIND .EQ. 92) THEN
         IF(VDUMAX .LT. VXX) THEN
            VDUMAX = VXX
            IDUMAX = IXX
         ENDIF
         IF(VDUMIN .GT. VXX) THEN
            VDUMIN = VXX
            IDUMIN = IXX
         ENDIF
**v 4.29f	 
      ELSEIF(KIND .EQ. 93) THEN
         IF(VDLMAX .LT. VXX) THEN
            VDLMAX = VXX
            IDLMAX = IXX
         ENDIF
**************************
** v 4.29g	 
         IF(VDLMIN .GT. VXX) THEN
            VDLMIN = VXX
            IDLMIN = IXX
         ENDIF
******************	 
      ENDIF

      RETURN
      END

************
      SUBROUTINE SECADJ (IUNIT,IUO,IOBS,B,NX,FATAL)

*** SECOND TRIP THRU ADJUSTMENT FILE

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL FATAL
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      CHARACTER*80 ACARD
      CHARACTER*2 ID
      DIMENSION B(*),NX(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON/UNITS/ LUNIT

**v6.3

      parameter  (MXSSN = 9999)
      integer*4  CC_counter
      LOGICAL    overwrite_const_coord_in_newbb
      character  CC_records*80
      COMMON/MM6/overwrite_const_coord_in_newbb
      COMMON/MM6_array/CC_records(MXSSN)                      

**So far for v6.3

      IF (LCS) THEN
        CALL HEAD
        CALL LINE (4)
        WRITE (LUNIT,1)
    1   FORMAT (' ******** CONSTRAINTS *************'/
     &          '0 OBS #'/)
      ENDIF
      CALL DIMCON (IOBS,IUO,B)

**V6.3

C  Initialize array CC_records

      CC_counter = 0
      do iii = 1,MXSSN
        CC_records(iii) = ' '
      enddo

**sofar for v6.3

  100 READ (IUNIT,2,END=777) ACARD
    2 FORMAT (A80)
      READ (ACARD,5) ID
    5 FORMAT (A2)

      IF (ID .EQ. 'CC') THEN

**v6.3

c       if (overwrite_const_coord_in_newbb) then
          CC_counter             = CC_counter + 1
          CC_records(CC_counter) = ACARD
c       endif

**so far for v6.3

        CALL SECCC (ACARD,IUO,IOBS,B,FATAL)
      ELSEIF (ID .EQ. 'QQ') THEN
        CALL SECQQ (ACARD,NX)
      ELSEIF (ID .EQ. 'SS') THEN
        CALL SECSS (ACARD,IUO,IOBS,B)
      ELSEIF (ID .EQ. 'CA') THEN
        IF (IDIM .NE. 1) CALL SECCA (ACARD,IUO,IOBS,B,NX,FATAL)
      ELSEIF (ID .EQ. 'CD') THEN
        CALL SECCD (ACARD,IUO,IOBS,B,NX,FATAL)
      ELSEIF (ID .EQ. 'CH') THEN
        CALL SECCH (ACARD,IUO,IOBS,B,NX,FATAL)
      ELSEIF (ID .EQ. 'CP') THEN
        CALL SECCP (ACARD,IUO,IOBS,B,NX,FATAL)
      ELSEIF (ID .EQ. 'CZ') THEN
        CALL SECCZ (ACARD,IUO,IOBS,B,NX,FATAL)
      ENDIF
      GO TO 100

*** END OF PROCESSING -- END OF FILE ENCOUNTERED

 777  IF (LCS) THEN
        CALL LINE (2)
        WRITE (LUNIT,3)
    3   FORMAT ('0************ END OF CONSTRAINTS *************')
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------------------------------
      SUBROUTINE SECBB (IUNIT,IUO,IOBS,B,NX,FATAL)

*** FORM OBS EQ. FOR NON-GPS, NON-DOPPLER BB RECORDS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 BCARD
      CHARACTER*2 IRT
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL FATAL,LSN
      DIMENSION B(*),NX(*)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON/UNITS/ LUNIT

*** PRINT HEADING, THEN PROCESS BBOOK

      IF (LBB) THEN
        CALL HEAD
        CALL LINE (4)
        WRITE (LUNIT,3)
  3     FORMAT (' ************ BLUE BOOK ************'/
     &          '0 OBS #'/)
      ENDIF

*** LOOP OVER RECORDS OF BLUE BOOK

  100 READ (IUNIT,2,END=777) BCARD
  2   FORMAT (A80)
      READ (BCARD,5) IRT
  5   FORMAT (7X,A2)

      LSN = .FALSE.

      IF (IRT .EQ. '20'  .OR.  IRT .EQ. '22') THEN
          IF ( .NOT. LDIR) CALL HORDIR (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
        ELSEIF (IRT .EQ. '30'  .OR.  IRT .EQ. '32') THEN
          IF ( .NOT. LANG) CALL HORANG (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
        ELSEIF (IRT .EQ. '40'  .OR.  IRT .EQ. '42') THEN
          IF ( .NOT. LZEN) CALL VERTAN (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
        ELSEIF (IRT .EQ. '52'  .OR.  IRT .EQ. '54') THEN
          IF ( .NOT. LDIS) CALL DISTOB (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
        ELSEIF (IRT .EQ. '60') THEN
          IF ( .NOT. LAZI) CALL ASTRAZ (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
        ELSEIF (IRT .EQ. '12') THEN
          CALL FIR12(BCARD)
      ENDIF

*** ECHO THE BLUE BOOK OBSERVATIONS

      IF (LBB) CALL ECHOBB (BCARD,IOBS,LSN)
      GO TO 100

*** END OF PROCESSING -- END OF FILE ENCOUNTERED ***

  777 IF (LBB) THEN
        CALL LINE (2)
        WRITE (LUNIT,4)
  4     FORMAT (' ******** END OF BLUE BOOK ********')
      ENDIF

      RETURN
      END
      SUBROUTINE SECCA (ACARD,IUO,IOBS,B,NX,FATAL)

*** WRITE CONSTRAINTS FOR AZIMUTH RECORDS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      CHARACTER*80 ACARD
      CHARACTER*4 ASS
      LOGICAL GETSSN,FATAL
      LOGICAL ADDCON
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      DIMENSION B(*),NX(*)
      DIMENSION IC(LENC),C(LENC)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) ISSN,JSSN,ID,IM,ASS,SD
 1    FORMAT (2X,2I4, 3X,I3,I2,A4,F5.2)
      CALL NBLANK (ASS,4,IBLK)
      READ (ASS,2) SS
 2    FORMAT (F4.2)

*** STD DEV DEFAULT IS 0.01 ARC SECONDS

      IF (ACARD(23:27) .EQ. '     ') SD = 0.01D0
      SD = SD/(3600.D0*RAD)
      NAZ = NAZ+1

      IF ( .NOT. GETSSN(ISSN,ISN)) THEN
        CALL LINE (3)
        WRITE(LUNIT,4) ACARD
 4      FORMAT ('0NO *80* RECORD FOR -- ',A80/)
      ELSEIF ( .NOT. GETSSN(JSSN,JSN)) THEN
        CALL LINE (3)
        WRITE (LUNIT,4) ACARD
      ELSEIF (ACARD(14:22) .NE. '         ') THEN
        KIND = 12
        NCON = NCON+1
        IOBS = IOBS+1
        IAUX = 0
        IVF = 0
        IGRT = 0
        OBSB = (ID+IM/60.D0+SS/3600.D0)/RAD
        IF (OBSB .LT. 0.D0) OBSB = OBSB + PI + PI
        IF (OBSB .GE. PI+PI) OBSB = OBSB - PI - PI
        CALL FORMIC (KIND,ISN,JSN,IDUMMY,IC,LENG,IGRT)
        CALL FORMC (KIND,C,B,ISN,JSN,IDUMMY,IGRT)
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IDUMMY,IGRT)
        CMO = OBS0 - OBSB
        IF (CMO .GT. PI) THEN
          CMO = CMO - PI - PI
        ELSEIF (CMO .LT. -PI) THEN
          CMO = CMO + PI + PI
        ENDIF
        IF (IMODE .EQ. 0) CMO = 0.D0
        VSD = CMO/SD
        IF (DABS(VSD) .GT. VP) THEN
          CALL LINE (1)
          WRITE (LUNIT,19) IOBS,VSD
   19     FORMAT (1X,'   OBS# =',I5,F70.1,' *** LARGE MISCLOSURE')
          IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
          IF ( .NOT. LCS) THEN
            CALL LINE (1)
            WRITE (LUNIT,12) ACARD
   12       FORMAT (10X, A80)
          ENDIF
        ENDIF

*** UPDATE CONNECTIVITY

        IF ( .NOT. ADDCON(IC,LENG,NX)) THEN
          WRITE (LUNIT,667)
 667      FORMAT ('0INSUFFICIENT STORAGE FOR CONSTRAINED ASTRO.
     &            AZIMUTH'/)
          write (lunit,*) 'subs2, 68'
          CALL ABORT2
        ENDIF
        WRITE (IUO) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,IVF,
     &              IAUX,IGRT

      ELSEIF (ACARD(14:22) .EQ. '         ') THEN
        KIND = 12
        NCON = NCON+1
        IOBS = IOBS+1
        IAUX = 0
        IVF = 0
        IGRT = 0
        CALL FORMIC (KIND,ISN,JSN,IDUMMY,IC,LENG,IGRT)
        CALL FORMC (KIND,C,B,ISN,JSN,IDUMMY,IGRT)
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IDUMMY,IGRT)
        OBSB = OBS0
        CMO = 0.D0

*** UPDATE CONNECTIVITY

        IF ( .NOT. ADDCON(IC,LENG,NX)) THEN
          WRITE (LUNIT,667)
          write (lunit,*) 'subs2, 69'
          CALL ABORT2
        ENDIF
        WRITE (IUO) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,IVF,
     &              IAUX,IGRT

      ENDIF

*** ECHO CONSTRAINT

      IF (LCS) THEN
        CALL LINE (1)
        WRITE (LUNIT,7) IOBS,ACARD
 7      FORMAT (I7,3X,A80)
      ENDIF

      RETURN
      END
      
C---------------------------------------------------------------------------------------------
      SUBROUTINE SECCC (ACARD,IUO,IOBS,B,FATAL)

*** WRITE CONSTRAINTS FOR FIXED COORDINATE RECORD

**v 4.28 read ellip ht in CC rec cc 81-87 xxxx(.)xxx

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, LENC = 10 )
      CHARACTER*80 ACARD
      CHARACTER*1 ADLA,ADLO,ACODE
      CHARACTER*7 AHT,ASLA,ASLO
      LOGICAL GETSSN,FATAL
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL ELFLAG,DFFLAG
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL L2HLF, LEHT
      DIMENSION B(*)
      DIMENSION IC(LENC),C(LENC)
      COMMON /FLAGS/  ELFLAG(MXSSN), DFFLAG(MXSSN)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) ISSN,SDLA,SDLO,SDHT,
     &               IDLA,IMLA,ASLA,ADLA,
** v 4.29j
     &               IDLO,IMLO,ASLO,ADLO,AHT,ACODE
*    &               IDLO,IMLO,ASLO,ADLO,AHT,EHT  
     
    1 FORMAT (10X,I4,    3F6.2,
     &        12X,2I2,A7,A1,
     &             I3,I2,A7,A1,A7,A1)
*    &             I3,I2,A7,A1,A7,T81,F7.3)

C  Convert coordinate standard deviations from cm back to mm

      SDLA = SDLA*0.01d0
      SDLO = SDLO*0.01d0
      SDHT = SDHT*0.01d0

C  Done with unit conversion from cm to mm (the way it was before the modification)
     
      IAUX = 0
      IVF = 0
      IGRT = 0

      CALL NBLANK (ASLA,5,IBLK)
      READ (ASLA,10) ISLA
      CALL NBLANK (ASLO,5,IBLK)
      READ (ASLO,10) ISLO
   10 FORMAT (I7)
      CALL NBLANK (AHT,3,IBLK)
      READ (AHT,5) HT
    5 FORMAT (F7.3)

      IF ( .NOT. GETSSN(ISSN,ISN)) THEN
        CALL LINE (3)
        WRITE (LUNIT,2) ACARD
    2   FORMAT ('0NO *80* RECORD FOR --',A80/)
      ELSE

*** STD. DEV DEFAULT WEIGHT IS 0.1 MM

        IF (ACARD(15:20) .EQ. '      ') SDLA = 0.0001D0
        IF (ACARD(21:26) .EQ. '      ') SDLO = 0.0001D0
        IF (ACARD(27:32) .EQ. '      ') SDHT = 0.0001D0

*** KIND = 1 -- NORTH COORDINATE SHIFT CONSTRAINT (LATITUDE)(DIM=2&3)

        IF (ACARD(45:55) .NE. '          '  .AND.  IDIM .NE. 1) THEN
          KIND = 1
          NCON = NCON+1
          IOBS = IOBS+1
          IF (ADLA .EQ. 'S') THEN
            SIGN = -1.D0
          ELSE
            SIGN = 1.D0
          ENDIF
          CALL GETRAD (IDLA,IMLA,ISLA,SIGN,GLATB)
          CALL FORMIC (KIND,ISN,IDUMMY,IDUMM2,IC,LENG,IGRT)
          CALL FORMC (KIND,C,B,IDUMM1,IDUMM2,IDUMM3,IGRT)
          CALL COMPOB (KIND,OBS0,B,GLATB,ISN,IDUMMY,IDUMM2,IGRT)
          OBSB = GLATB
          CMO = OBS0
          IF (IMODE .EQ. 0) CMO = 0.D0
          VSD = CMO/SDLA
          IF (DABS(VSD) .GT. VP) THEN
            CALL LINE (1)
            WRITE (LUNIT,11) IOBS,VSD
   11       FORMAT (5X,I5,F18.1,' *** LARGE MISCLOSURE ')
            IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
            IF ( .NOT. LCS) THEN
              CALL LINE (1)
              WRITE (LUNIT,12) ACARD
   12         FORMAT (10X, A80)
            ENDIF
          ENDIF
          WRITE (IUO) KIND,ISN,IDUMMY,IC,C,LENG,CMO,OBSB,SDLA,IOBS,IVF,
     &                IAUX,IGRT
        ENDIF

*** KIND = 2 -- EAST COORDINATE SHIFT CONSTRAINT (LONGITUDE)(DIM=2&3)

        IF (ACARD(57:68) .NE. '            '  .AND.  IDIM .NE. 1) THEN
          KIND = 2
          NCON = NCON+1
          IOBS = IOBS+1
          IF (ADLO .EQ. ' '  .OR.  ADLO .EQ. 'W') THEN
            SIGN = -1.D0
          ELSE
            SIGN = 1.D0
          ENDIF
          CALL GETRAD (IDLO,IMLO,ISLO,SIGN,GLONB)
          CALL FORMIC (KIND,ISN,IDUMMY,IDUMM2,IC,LENG,IGRT)
          CALL FORMC (KIND,C,B,IDUMM1,IDUMM2,IDUMM3,IGRT)
          CALL COMPOB (KIND,OBS0,B,GLONB,ISN,IDUMMY,IDUMM2,IGRT)
          OBSB = GLONB
          CMO = OBS0
          IF (IMODE .EQ. 0) CMO = 0.D0
          VSD = CMO/SDLO
          IF (DABS(VSD) .GT. VP) THEN
            CALL LINE (1)
            WRITE (LUNIT,11) IOBS,VSD
            IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
            IF ( .NOT. LCS) THEN
              CALL LINE (1)
              WRITE (LUNIT,12) ACARD
            ENDIF
          ENDIF
          WRITE (IUO) KIND,ISN,IDUMMY,IC,C,LENG,CMO,OBSB,SDLO,IOBS,IVF,
     &                IAUX,IGRT
        ENDIF

*** KIND = 3 -- VERTICAL COORDINATE SHIFT CONSTRAINT (UP)
***             (DIM = 1 & 3)  OR  (DIM=2.5)
*** v 4.28
*** orthometric height

        IF ( ACARD(70:76) .NE. '       '  .AND.
     &       ( IDIM .NE. 2  .OR.
     &         ( IDIM .EQ. 2  .AND.  L2HLF  .AND.  I2HLF(ISN) .EQ. 1 )
     &       )  ) THEN
          KIND = 3
          NCON = NCON+1
          IOBS = IOBS+1
**deleted out for v 4.28 and put back in for v 4.29j	  
          IF (ACODE.EQ.'E'.OR.ACODE.EQ.'e') THEN
            EHB = HT
          ELSEIF (ELFLAG(ISN)) THEN
            CALL GETGH (GHT0,ISN,B)
            EHB = HT+GHT0
          ELSE
            CALL GETMSL (GMSL0,ISN,B)
            EHB = GMSL0+HT
          ENDIF
***BIGADJUST code****************
*
*         IF (ELFLAG(ISN)) THEN
*           CALL GETGH (GHT0,ISN,B)
*           EHB = HT+GHT0
*         ELSE
*           CALL GETMSL (GMSL0,ISN,B)
*           EHB = GMSL0+HT
*         ENDIF
*
***************************************
          CALL FORMIC (KIND,ISN,IDUMMY,IDUMM2,IC,LENG,IGRT)
          CALL FORMC (KIND,C,B,IDUMM1,IDUMM2,IDUMM3,IGRT)
          CALL COMPOB (KIND,OBS0,B,ADUMMY,ISN,IDUMMY,IDUMM2,IGRT)
          OBSB = EHB
          CMO = OBS0-OBSB
          IF (IMODE .EQ. 0) CMO = 0.D0
          VSD = CMO/SDHT
          IF (DABS(VSD) .GT. VP) THEN
            CALL LINE (1)
            WRITE (LUNIT,11) IOBS,VSD
            IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
            IF ( .NOT. LCS) THEN
              CALL LINE (1)
              WRITE (LUNIT,12) ACARD
            ENDIF
          ENDIF
          WRITE (IUO) KIND,ISN,IDUMMY,IC,C,LENG,CMO,OBSB,SDHT,IOBS,IVF,
     &                IAUX,IGRT
        ENDIF

**v 4.28
*** ellipsoidal height

*       IF ( ACARD(81:87) .NE. '       '  .AND.
*    &       ( IDIM .NE. 2  .OR.
*    &         ( IDIM .EQ. 2  .AND.  L2HLF  .AND.  I2HLF(ISN) .EQ. 1 )
*    &       )  ) THEN
*         KIND = 3
*         NCON = NCON+1
*         IOBS = IOBS+1
*         EHB = EHT
*         CALL FORMIC (KIND,ISN,IDUMMY,IDUMM2,IC,LENG,IGRT)
*         CALL FORMC (KIND,C,B,IDUMM1,IDUMM2,IDUMM3,IGRT)
*         CALL COMPOB (KIND,OBS0,B,ADUMMY,ISN,IDUMMY,IDUMM2,IGRT)
*         OBSB = EHB
*         CMO = OBS0-OBSB
*         IF (IMODE .EQ. 0) CMO = 0.D0
*         VSD = CMO/SDHT
*         IF (DABS(VSD) .GT. VP) THEN
*           CALL LINE (1)
*           WRITE (LUNIT,11) IOBS,VSD
*           IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
*           IF ( .NOT. LCS) THEN
*             CALL LINE (1)
*             WRITE (LUNIT,12) ACARD
*           ENDIF
*         ENDIF
*         WRITE (IUO) KIND,ISN,IDUMMY,IC,C,LENG,CMO,OBSB,SDHT,IOBS,IVF,
*    &                IAUX,IGRT
*       ENDIF
*	
****************************


*** ECHO OBSERVATION NUMBERS

      IF (LCS) THEN
        CALL LINE (1)
        WRITE (LUNIT,3) IOBS,ACARD
      ENDIF
    3 FORMAT (1X,I6,3X,A80)
      ENDIF
      
      RETURN
      END
C---------------------------------------------------------------------------------------------------
      SUBROUTINE SECCD (ACARD,IUO,IOBS,B,NX,FATAL)

*** WRITE CONSTRAINTS FOR DISTANCE RECORDS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      CHARACTER*80 ACARD
      LOGICAL GETSSN,FATAL
      LOGICAL ADDCON
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      DIMENSION B(*),NX(*)
      DIMENSION IC(LENC),C(LENC)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) ISSN,JSSN,OBSB,SD
 1    FORMAT (2X,2I4,F12.4,F5.4)

*** STD DEV DEFAULT IS 0.1 MM

      IF (ACARD(23:27) .EQ. '     ') SD = 0.0001D0
      NDS = NDS+1

      IF ( .NOT. GETSSN(ISSN,ISN)) THEN
        CALL LINE (3)
        WRITE (LUNIT,4) ACARD
 4      FORMAT ('0NO *80* RECORD FOR -- ',A80/)
      ELSEIF ( .NOT. GETSSN(JSSN,JSN)) THEN
        CALL LINE (3)
        WRITE (LUNIT,4) ACARD
      ELSEIF (ACARD(11:22) .NE. '            ') THEN
        KIND = 13
        NCON = NCON+1
        IOBS = IOBS+1
        IAUX = 0
        IVF = 0
        IGRT = 0
        CALL FORMIC (KIND,ISN,JSN,IAUX,IC,LENG,IGRT)
        CALL FORMC (KIND,C,B,ISN,JSN,IAUX,IGRT)
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        CMO = OBS0-OBSB
        IF (IMODE .EQ. 0) CMO = 0.D0
        VSD = CMO/SD
        IF (DABS(VSD) .GT. VP) THEN
          CALL LINE (1)
          WRITE (LUNIT,11) IOBS,VSD
   11     FORMAT (1X,'   OBS# =',I5,F70.1,' *** LARGE MISCLOSURE')
          IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
          IF ( .NOT. LCS) THEN
            CALL LINE (1)
            WRITE (LUNIT,12) ACARD
   12       FORMAT (10X, A80)
          ENDIF
        ENDIF

*** UPDATE CONNECTIVITY

        IF ( .NOT. ADDCON(IC,LENG,NX)) THEN
          WRITE (LUNIT,667)
 667      FORMAT ('0INSUFFICIENT STORAGE FOR CONSTRAINED DISTANCES'/)
          write (lunit,*) 'subs2, 70'
          CALL ABORT2
        ENDIF
        WRITE (IUO) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,IVF,
     &              IAUX,IGRT

      ELSEIF (ACARD(11:22) .EQ. '            ') THEN
        KIND = 13
        NCON = NCON+1
        IOBS = IOBS+1
        IAUX = 0
        IVF = 0
        IGRT = 0
        CALL FORMIC (KIND,ISN,JSN,IAUX,IC,LENG,IGRT)
        CALL FORMC (KIND,C,B,ISN,JSN,IAUX,IGRT)
        CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
        OBSB = OBS0
        CMO = 0.D0

*** UPDATE CONNECTIVITY

        IF ( .NOT. ADDCON(IC,LENG,NX)) THEN
          WRITE (LUNIT,667)
          write (lunit,*) 'subs2, 71'
          CALL ABORT2
        ENDIF
        WRITE (IUO) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,IVF,
     &              IAUX,IGRT

      ENDIF

*** ECHO CONSTRAINT

      IF (LCS) THEN
        CALL LINE (1)
        WRITE (LUNIT,7) IOBS,ACARD
 7      FORMAT (I7,3X,A80)
      ENDIF

      RETURN
      END
      SUBROUTINE SECCH (ACARD, IUO, IOBS, B, NX, FATAL)

*** WRITE CONSTRAINTS FOR HEIGHT DIFFERENCES

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, LENC = 10 )
      CHARACTER*80 ACARD
      LOGICAL GETSSN, FATAL
      LOGICAL ADDCON
      LOGICAL ELFLAG, DFFLAG
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      DIMENSION B(*), NX(*)
      DIMENSION IC(LENC), C(LENC)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /FLAGS/  ELFLAG(MXSSN), DFFLAG(MXSSN)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) ISSN, JSSN, OBSB, SD
    1 FORMAT (2X, 2I4, F12.4, F5.4)

      IF ( .NOT. GETSSN(ISSN, ISN) ) THEN
        CALL LINE (3)
        WRITE (LUNIT,4) ACARD
    4   FORMAT ('0NO *80* RECORD FOR -- ', A80, /)
      ELSEIF ( .NOT. GETSSN(JSSN, JSN) ) THEN
        CALL LINE (3)
        WRITE (LUNIT,4) ACARD
      ELSEIF ( IDIM .EQ. 2  .AND.  L2HLF  .AND.
     &         ( I2HLF(ISN) .NE. 1  .OR.  I2HLF(JSN) .NE. 1 )
     &       ) THEN

*** IF THIS IS A 2 DIM ADJUSTMENT AND ONE OF THE STATIONS IS NOT A
*** DUAL HEIGHT STATION THEN RETURN

        RETURN

      ELSE

*** STD DEV DEFAULT IS 0.1 MM

        IF (ACARD(23:27) .EQ. '     ') SD = 0.0001D0

        IF (ACARD(11:22) .NE. '            ') THEN
          IF ( ELFLAG(ISN) ) THEN
            IF ( ELFLAG(JSN) ) THEN
              KIND = 16
            ELSE
              KIND = 17
            ENDIF
          ELSE
            IF ( ELFLAG(JSN) ) THEN
              KIND = 17
            ELSE
              KIND = 15
            ENDIF
          ENDIF
          NCON = NCON + 1
          IOBS = IOBS + 1
          IAUX = 0
          IVF = 0
          IGRT = 0
          CALL FORMIC (KIND, ISN, JSN, IAUX, IC, LENG, IGRT)
          CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
          CALL COMPOB (KIND, OBS0, B, OBSB, ISN, JSN, IAUX, IGRT)
          CMO = OBS0 - OBSB
          IF (IMODE .EQ. 0) CMO = 0.D0
          VSD = CMO/SD
          IF (DABS(VSD) .GT. VP) THEN
            CALL LINE (1)
            WRITE (LUNIT,11) IOBS, VSD
   11       FORMAT (1X, '   OBS# =', I5, F70.1, ' *** LARGE MISCLOSURE')
            IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
            IF ( .NOT. LCS) THEN
              CALL LINE (1)
              WRITE (LUNIT,12) ACARD
   12         FORMAT (10X, A80)
            ENDIF
          ENDIF

*** UPDATE CONNECTIVITY

          IF ( .NOT. ADDCON(IC, LENG, NX) ) THEN
            WRITE (LUNIT,667)
 667        FORMAT ('0INSUFFICIENT STORAGE FOR CONSTRAINED
     &              HEIGHT DIFFERENCES', /)
          write (lunit,*) 'subs2, 72'
            CALL ABORT2
          ENDIF
          WRITE (IUO) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD, IOBS,
     &                IVF, IAUX, IGRT

        ELSEIF (ACARD(11:22) .EQ. '            ') THEN
          IF ( ELFLAG(ISN) ) THEN
            IF ( ELFLAG(JSN) ) THEN
              KIND = 16
            ELSE
              KIND = 17
            ENDIF
          ELSE
            IF ( ELFLAG(JSN) ) THEN
              KIND = 17
            ELSE
              KIND = 15
            ENDIF
          ENDIF
          NCON = NCON + 1
          IOBS = IOBS + 1
          IAUX = 0
          IVF = 0
          IGRT = 0
          CALL FORMIC (KIND, ISN, JSN, IAUX, IC, LENG, IGRT)
          CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
          CALL COMPOB (KIND, OBS0, B, OBSB, ISN, JSN, IAUX, IGRT)
          OBSB = OBS0
          CMO = 0.D0

*** UPDATE CONNECTIVITY

          IF ( .NOT. ADDCON(IC, LENG, NX) ) THEN
            WRITE (LUNIT,667)
          write (lunit,*) 'subs2, 73'
            CALL ABORT2
          ENDIF
          WRITE (IUO) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD, IOBS,
     &                IVF, IAUX, IGRT

        ENDIF

*** ECHO CONSTRAINT

        IF (LCS) THEN
          CALL LINE (1)
          WRITE (LUNIT,7) IOBS, ACARD
 7        FORMAT (I7, 3X, A80)
        ENDIF

      ENDIF

      RETURN
      END

      SUBROUTINE SECCP (ACARD, IUO, IOBS, B, NX, FATAL)

*** WRITE CONSTRAINTS FOR HORIZONTAL POSITION DIFFERENCES

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, LENC = 10 )
      CHARACTER*80 ACARD
      LOGICAL GETSSN, FATAL
      LOGICAL ADDCON
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      DIMENSION B(*), NX(*)
      DIMENSION IC(LENC), C(LENC)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) ISSN, JSSN, DN, SDN, DE, SDE
    1 FORMAT (2X, 2I4, 4F10.4)

      IF ( .NOT. GETSSN(ISSN, ISN) ) THEN
        CALL LINE (3)
        WRITE (LUNIT,4) ACARD
    4   FORMAT ('0NO *80* RECORD FOR -- ', A80, /)
      ELSEIF ( .NOT. GETSSN(JSSN, JSN) ) THEN
        CALL LINE (3)
        WRITE (LUNIT,4) ACARD
      ELSEIF (IDIM .EQ. 1) THEN
        RETURN

      ELSE

*** STD DEV DEFAULT IS 0.1 MM

        IF (ACARD(21:30) .EQ. '          ') SDN = 0.0001D0
        IF (ACARD(41:50) .EQ. '          ') SDE = 0.0001D0

*** DEFAULT VALUES FOR DN AND DE ARE ZERO

        IF (ACARD(11:20) .EQ. '          ') DN = 0.D0
        IF (ACARD(31:40) .EQ. '          ') DE = 0.D0

*** DELTA NORTH

        KIND = 21
        OBSB = DN
        SD = SDN
        NCON = NCON + 1
        IOBS = IOBS + 1
        IAUX = 0
        IVF = 0
        IGRT = 0
        CALL FORMIC (KIND, ISN, JSN, IAUX, IC, LENG, IGRT)
        CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
        CALL COMPOB (KIND, OBS0, B, OBSB, ISN, JSN, IAUX, IGRT)
        CMO = OBS0 - OBSB
        IF (IMODE .EQ. 0) CMO = 0.D0
        VSD = CMO/SD
        IF (DABS(VSD) .GT. VP) THEN
          CALL LINE (1)
          WRITE (LUNIT,11) IOBS, VSD
   11     FORMAT (1X, '   OBS# =', I5, F70.1, ' *** LARGE MISCLOSURE')
          IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
          IF ( .NOT. LCS) THEN
            CALL LINE (1)
            WRITE (LUNIT,12) ACARD
   12       FORMAT (10X, A80)
          ENDIF
        ENDIF

*** UPDATE CONNECTIVITY

        IF ( .NOT. ADDCON(IC, LENG, NX) ) THEN
          WRITE (LUNIT,667)
  667     FORMAT ('0INSUFFICIENT STORAGE FOR CONSTRAINED
     &              HORIZONTAL POSITION DIFFERENCES', /)
          write (lunit,*) 'subs2, 74'
          CALL ABORT2
        ENDIF
        WRITE (IUO) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD, IOBS,
     &              IVF, IAUX, IGRT

*** DELTA EAST

        KIND = 22
        OBSB = DE
        SD = SDE
        NCON = NCON + 1
        IOBS = IOBS + 1
        IAUX = 0
        IVF = 0
        IGRT = 0
        CALL FORMIC (KIND, ISN, JSN, IAUX, IC, LENG, IGRT)
        CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
        CALL COMPOB (KIND, OBS0, B, OBSB, ISN, JSN, IAUX, IGRT)
        CMO = OBS0 - OBSB
        IF (IMODE .EQ. 0) CMO = 0.D0
        VSD = CMO/SD
        IF (DABS(VSD) .GT. VP) THEN
          CALL LINE (1)
          WRITE (LUNIT,11) IOBS, VSD
          IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
          IF ( .NOT. LCS) THEN
            CALL LINE (1)
            WRITE (LUNIT,12) ACARD
          ENDIF
        ENDIF

*** UPDATE CONNECTIVITY

        IF ( .NOT. ADDCON(IC, LENG, NX) ) THEN
          WRITE (LUNIT,667)
          write (lunit,*) 'subs2, 75'
          CALL ABORT2
        ENDIF
        WRITE (IUO) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD, IOBS,
     &              IVF, IAUX, IGRT


*** ECHO CONSTRAINT

        IF (LCS) THEN
          CALL LINE (1)
          WRITE (LUNIT,7) IOBS, ACARD
 7        FORMAT (I7, 3X, A80)
        ENDIF

      ENDIF

      RETURN
      END
      SUBROUTINE SECCZ (ACARD, IUO, IOBS, B, NX, FATAL)

*** WRITE CONSTRAINTS FOR ZENITH DISTANCES

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, LENC = 10 )
      CHARACTER*80 ACARD
      CHARACTER*4 ASS
      LOGICAL GETSSN, FATAL
      LOGICAL ADDCON
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      DIMENSION B(*), NX(*)
      DIMENSION IC(LENC), C(LENC)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) ISSN, JSSN, ID, IM, ASS, SD
    1 FORMAT (2X, 2I4, 3X, I3, I2, A4, F5.2)
      CALL NBLANK (ASS, 4, IBLK)
      READ (ASS,2) SS
    2 FORMAT (F4.2)

      IF ( .NOT. GETSSN(ISSN, ISN) ) THEN
        CALL LINE (3)
        WRITE (LUNIT,4) ACARD
    4   FORMAT ('0NO *80* RECORD FOR -- ', A80, /)
      ELSEIF ( .NOT. GETSSN(JSSN, JSN) ) THEN
        CALL LINE (3)
        WRITE (LUNIT,4) ACARD
      ELSEIF ( IDIM .EQ. 2  .AND.  L2HLF  .AND.
     &         ( I2HLF(ISN) .EQ. 0  .OR.  I2HLF(JSN) .EQ. 0 )
     &       ) THEN

*** IF THIS IS A 2 DIM ADJUSTMENT AND ONE OF THE STATIONS IS NOT A
*** DUAL HEIGHT STATION THEN RETURN

        RETURN

      ELSE

*** STD DEV DEFAULT IS 0.01 ARC SECONDS

        IF (ACARD(23:27) .EQ. '     ') SD = 0.01D0
        SD = SD/(3600.D0*RAD)
        NZD = NZD + 1
        IF (ACARD(14:22) .NE. '         ') THEN
          KIND = 14
          NCON = NCON + 1
          IOBS = IOBS + 1
          IAUX = 0
          IVF = 0
          IGRT = 0
          OBSB = (ID + IM/60.D0 + SS/3600.D0)/RAD
          CALL FORMIC (KIND, ISN, JSN, IAUX, IC, LENG, IGRT)
          CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
          CALL COMPOB (KIND, OBS0, B, OBSB, ISN, JSN, IAUX, IGRT)
          CMO = OBS0 - OBSB
          IF (IMODE .EQ. 0) CMO = 0.D0
          VSD = CMO/SD
          IF (DABS(VSD) .GT. VP) THEN
            CALL LINE (1)
            WRITE (LUNIT,11) IOBS, VSD
   11       FORMAT (1X, '   OBS# =', I5, F70.1, ' *** LARGE MISCLOSURE')
            IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
            IF ( .NOT. LCS) THEN
              CALL LINE (1)
              WRITE (LUNIT,12) ACARD
   12         FORMAT (10X, A80)
            ENDIF
          ENDIF

*** UPDATE CONNECTIVITY

          IF ( .NOT. ADDCON(IC, LENG, NX) ) THEN
            WRITE (LUNIT,667)
  667       FORMAT ('0INSUFFICIENT STORAGE FOR CONSTRAINED ZEN.
     &              DISTANCES', /)
          write (lunit,*) 'subs2, 76'
            CALL ABORT2
          ENDIF
          WRITE (IUO) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD, IOBS,
     &                IVF, IAUX, IGRT

        ELSEIF (ACARD(14:22) .EQ. '         ') THEN
          KIND = 14
          NCON = NCON + 1
          IOBS = IOBS + 1
          IAUX = 0
          IVF = 0
          IGRT = 0
          CALL FORMIC (KIND, ISN, JSN, IAUX, IC, LENG, IGRT)
          CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
          CALL COMPOB (KIND, OBS0, B, OBSB, ISN, JSN, IAUX, IGRT)
          OBSB = OBS0
          CMO = 0.D0

*** UPDATE CONNECTIVITY

          IF ( .NOT. ADDCON(IC, LENG, NX) ) THEN
            WRITE (LUNIT,667)
          write (lunit,*) 'subs2, 77'
            CALL ABORT2
          ENDIF
          WRITE (IUO) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD, IOBS,
     &                IVF, IAUX, IGRT

        ENDIF

*** ECHO CONSTRAINT

        IF (LCS) THEN
          CALL LINE (1)
          WRITE (LUNIT,7) IOBS, ACARD
    7     FORMAT (I7, 3X, A80)
        ENDIF

      ENDIF

      RETURN
      END
      SUBROUTINE SECDOP (IUNIT, IUO, IOBS, B, NX, FATAL)

*** FORM OBS EQ. DOPPLER RECORDS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      DIMENSION IC(LENC), C(LENC)
      DIMENSION B(*), NX(*)
      DIMENSION COVECF(3,3)
      CHARACTER*89 DCARD
      CHARACTER*2 ID
      CHARACTER*1 TC, AREJ
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL GETPRM, GETIVF, GETGRT, GETSSN, ADDCON
      LOGICAL FATAL, LSN
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON/UNITS/ LUNIT

*** PRINT HEADING, THEN PROCESS DFILE

      IF (LDF) THEN
        CALL HEAD
        CALL LINE (4)
        WRITE (LUNIT,3)
  3     FORMAT (' ************ DOPPLER OBSERVATIONS ************'/
     &          '0 OBS #'/)
      ENDIF

*** ENTER PROCESSING LOOP--EXIT ON END OF FILE

  100 READ (IUNIT,1,END=777) DCARD
    1 FORMAT (A89)
      READ (DCARD,5) ID
    5 FORMAT (A2)
      LSN = .FALSE.

      IF (ID .EQ. 'DP') THEN
        READ (DCARD,10) I83, ISSN, X, Y, Z, SIGN, SIGE, SIGU,
     &                  CORNE, CORNU, COREU, AREJ, IYR, IMO, IDY
   10   FORMAT (BZ, 2X, I1, I5, 3F11.3, 3F5.3, 3F8.6, A1, I4, I2, I2)
        TC = 'Z'
        IF (IMO .EQ. 0) IMO = 7
        IF (IDY .EQ. 0) IDY = 1
        IHR = 12
        IMN = 0
        ICODE = 99
        JSN = 0
        NDOP = NDOP+1

        IF (AREJ .EQ. 'R'  .OR.  AREJ .EQ. 'O'  .OR.
     &      AREJ .EQ. 'F') THEN
          NREJ = NREJ + 1
        ELSE
          IF ( .NOT. GETSSN(ISSN,ISN) ) THEN
            CALL LINE (3)
            WRITE (LUNIT,2) DCARD
 2          FORMAT ('0NO *80* RECORD FOR--', A89/)
          ELSEIF (I83 .NE. 1) THEN
            WRITE (LUNIT, 130) I83, DCARD
  130       FORMAT (//, ' ERROR - DOPPLER DATA REFERENCE FRAME CODE',
     &              ' IS', I2, ' FOR THE FOLLOWING RECORD',
     &               /, ' - MUST BE A ''1'' FOR NAD83.', /, A89)
          write (lunit,*) 'subs2, 78'
            CALL ABORT2
          ELSE

            LSN = .TRUE.
            CALL TOMNT (IYR, IMO, IDY, IHR, IMN, TC, ITIME)
            IF ( .NOT. GETPRM(ICODE, ITIME, IAUX) ) IAUX = 0
            IF ( .NOT. GETIVF(ICODE, ITIME, IVF) ) IVF = 0
            IF ( .NOT. GETGRT(ICODE, ITIME, IGRT) ) IGRT = 0
            CALL DOPCOV (ISSN, X, Y, Z, SIGN, SIGE, SIGU,
     &                   CORNE, CORNU, COREU, COVECF)


*** KIND = 18 DOPPLER X

            KIND = 18
            IOBS = IOBS+1
            NOBS = NOBS+1
            CALL FORMIC (KIND, ISN, JSN, IAUX, IC, LENG, IGRT)
            CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
            CALL COMPOB (KIND, OBS0, B, X, ISN, JSN, IAUX, IGRT)
            CMO = OBS0 - X
            SIGX = DSQRT( COVECF(1,1) )
            IF (IMODE .EQ. 0) CMO = 0.D0
            CALL BIGV (KIND, ISN, JSN, IOBS, IVF, CMO, SIGX, FATAL)
            WRITE (IUO) KIND, ISN, JSN, IC, C, LENG, CMO, X, SIGX, IOBS,
     &                  IVF, IAUX, IGRT

*** UPDATE CONNECTIVITY (ONLY NEED TO DO ONCE PER DOPPLER RECORD)

            IF ( .NOT. ADDCON(IC, LENG, NX) ) THEN
              WRITE (LUNIT,666)
 666          FORMAT ('0INSUFFICIENT STORAGE FOR DOPPLER OBS.'/)
          write (lunit,*) 'subs2, 79'
              CALL ABORT2
            ENDIF

*** KIND = 19 DOPPLER Y

            KIND = 19
            IOBS = IOBS+1
            NOBS = NOBS+1
            CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
            CALL COMPOB (KIND, OBS0, B, Y, ISN, JSN, IAUX, IGRT)
            CMO = OBS0 - Y
            SIGY = DSQRT( COVECF(2,2) )
            IF (IMODE .EQ. 0) CMO = 0.D0
            CALL BIGV (KIND, ISN, JSN, IOBS, IVF, CMO, SIGY, FATAL)
            WRITE (IUO) KIND, ISN, JSN, IC, C, LENG, CMO, Y, SIGY, IOBS,
     &                  IVF, IAUX, IGRT

*** KIND = 20 DOPPLER Z

            KIND = 20
            IOBS = IOBS+1
            NOBS = NOBS+1
            CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
            CALL COMPOB (KIND, OBS0, B, Z, ISN, JSN, IAUX, IGRT)
            CMO = OBS0 - Z
            SIGZ = DSQRT( COVECF(3,3) )
            IF (IMODE .EQ. 0) CMO = 0.D0
            CALL BIGV (KIND, ISN, JSN, IOBS, IVF, CMO, SIGZ, FATAL)
            WRITE (IUO) KIND, ISN, JSN, IC, C, LENG, CMO, Z, SIGZ, IOBS,
     &                  IVF, IAUX, IGRT

*** WRITE THE COVARIANCE FUNCTION

            ISIGNL = 2000
            WRITE (IUO) ISIGNL, ( (COVECF(I,J), J = 1,3), I = 1,3 )

          ENDIF
        ENDIF
      ENDIF

*** ECHO THE DFILE OBSERVATIONS

      IF (LDF) CALL ECHODF (DCARD, IOBS, LSN, AREJ)

      GO TO 100

*** END OF PROCESSING -- END OF FILE ENCOUNTERED ***

  777 IF (LDF) THEN
        CALL LINE (2)
        WRITE (LUNIT,4)
    4   FORMAT ('0*********** END OF DOPPLER OBSERVATIONS *********')
      ENDIF

      RETURN
      END
      SUBROUTINE SECGPS (IUNIT, IUO, IUO3, IOBS, B, NX, G, FATAL)

c Mike Potterfield 1/10/06
c This subroutine has a new argument IOU3, a scratch file
c identifying rejected vectors, so that OBSSUM won't count
c them in the observational summary

*** FORM OBS. EQ FROM G FORMAT

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MAXB = 49999 )
      PARAMETER ( NVECS = 700, MXSSN = 9999 )
      DIMENSION ICM((NVECS+1)*3+1+3)
      DIMENSION KINDS(NVECS*3), ISNS(NVECS*3), JSNS(NVECS*3)
      DIMENSION LOBS(NVECS*3)
      LOGICAL FULL, FATAL, GETPRM, GETIVF, GETGRT
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LEB, LLB, LEG, LLG, LED, LLD
      LOGICAL L2HLF, LEHT
      LOGICAL CCFLAG,CVFLAG
      CHARACTER*120 GCARD
      CHARACTER*13 PROJID
      CHARACTER*1 ID, TC
      DIMENSION B(*), NX(*), G(*)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /GPS/    MAXVEC
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /PECHO/  LEB, LLB, LEG, LLG, LED, LLD
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT
      COMMON/BREC/ NCCNT(MAXB)

*** INITIALIZE GPS VECTOR TRANSFORMATIONS TO NAD83

***07-07-05      CALL INIT83
c Keep the f90 compiler from complaining
      TTIME = 1997.000
      CALL INIT83(TTIME)
*************

*** TEST FOR MAXIMUM NUMBER OF VECTORS ALLOWED IN GROUPS

      IF (MAXVEC .GT. NVECS) THEN
        WRITE (LUNIT,222) MAXVEC, NVECS
 222    FORMAT (' THE NUMBER OF VECTORS IN A GROUP =',I3,
     &          ' HAS EXCEEDED THE MAXIMUM ALLOCATED =',I3,/,
     &      ' *** FATAL -- INCREASE NVECS IN PARAMETER STATEMENTS ***'/)
          write (lunit,*) 'subs2, 80'
        CALL ABORT2
      ENDIF

      IF (LGF) THEN
        CALL HEAD
        CALL LINE (2)
        WRITE (LUNIT,3)
    3   FORMAT (' ******** GPS OBSERVATIONS *************'/)
      ENDIF

      FULL = .FALSE.
      CALL NEWICM (NICM)

*** ENTER PROCESSING LOOP--EXIT ON END OF FILE

      NOB = 0
** v 4.28*****************
      IVC = 1
**************************
*** 8-20-02	
        CCFLAG = .FALSE.
        CVFLAG = .FALSE.
*************************

  100 READ (IUNIT,1,END=777) GCARD
    1 FORMAT (A120)
**v6.2.1
      if (GCARD(1:1) == ' ') goto 100
**so far for v6.2.1
      READ (GCARD,5) ID
    5 FORMAT (A1)
      IF ( .NOT. LEG) THEN
        IF ( .NOT. LLG) THEN
          IF (LGF) THEN
            CALL LINE (1)
            WRITE (LUNIT,2) GCARD
    2       FORMAT (10X,A80)
          ENDIF
        ENDIF
      ENDIF

*** GROUP HEADER RECORD

      IF (ID .EQ. 'B') THEN
        NOB = NOB + 1
*** 8-20-02	
*       CCFLAG = .FALSE.
*       CVFLAG = .FALSE.
**********************8-26-02
*       IF (FULL) CALL OBSEQW (IUO, FULL, G, B, NX, NVEC, NR, NC, LENG,
*    &                         NICM, ICM, KINDS, ISNS, JSNS, LOBS,
*    &                         IAUX, IVF, FATAL, IOBS, IGRT, CVFLAG)
        IF (FULL) THEN
          CALL OBSEQW (IUO, FULL, G, B, NX, NVEC, NR, NC, LENG,
     &                         NICM, ICM, KINDS, ISNS, JSNS, LOBS,
     &                         IAUX, IVF, FATAL, IOBS, IGRT, CVFLAG,
     &                         PROJID)

c Mike Potterfield 2/10/06
c The scratch file IUO3 is being used to identify rejected vectors

          CCFLAG = .FALSE.
          CVFLAG = .FALSE.
        ENDIF 
**********************************************************************

        READ (GCARD,10) IYR1, IMO1, IDY1, IHR1, IMN1,
     &                  IYR2, IMO2, IDY2, IHR2, IMN2
   10   FORMAT (BZ,1X,2(I4,4I2))
        READ (GCARD,110) NVEC, I83, PROJID
  110   FORMAT (25X,I2,24X,I2,T91, A13)	
        TC = 'Z'
        ICODE = 25

*** ASSIGN C REC COUNT TO NVEC IF B REC CC 26-27 IS BLANK

        IF (NVEC .LT. 1) NVEC = NCCNT(NOB)

        IF (NVEC .LT. 1  .OR.  NVEC .GT. MAXVEC) THEN
          WRITE (LUNIT,11) GCARD, NVEC, MAXVEC
   11     FORMAT (5X,A80,5X,'NVEC =',I5,' EXCEEDS MAXVEC=',I5)
          write (lunit,*) 'subs2, 81'
          CALL ABORT2
        ENDIF
        CALL TOMNT (IYR1, IMO1, IDY1, IHR1, IMN1, TC, IOLD)
        CALL TOMNT (IYR2, IMO2, IDY2, IHR2, IMN2, TC, INEW)
        ITIME = (IOLD + INEW)/2

***07-07-05  
        T = IYR1 + IMO1/12.
        CALL INIT83(T)
*********

        IF ( .NOT. GETPRM(ICODE, ITIME, IAUX) ) IAUX = 0
        IF ( .NOT. GETIVF(ICODE, ITIME, IVF) ) IVF = 0
        IF ( .NOT. GETGRT(ICODE, ITIME, IGRT) ) IGRT = 0

        IF (I83 .ge. 7 .and. i83 .le. 11) then
            I83 = 5
        elseif (i83 .ge. 12 .and. i83 .le. 14) then
            i83 = 12
        elseif (i83 .ge. 15 .and. i83 .le. 17) then
            i83 = 15
        elseif (i83 .eq. 18) then
            i83 = 18	    
        elseif (i83 .eq. 19) then
            i83 = 19	    
	  elseif (i83 .eq. 20) then
	    i83 = 19    
        elseif (i83 .ge. 21 .and. i83 .le. 24) then
	    i83 = 21
c   this is new code 11/6/07
        elseif (i83 .eq. 26) then
**         Do not do anything
        elseif (i83 .eq. 27) then
**         Do not do anything
        elseif (i83 .eq. 28) then
**         Do not do anything
        endif

 	
        NR = 3*NVEC
        IF ( .NOT. L2HLF) THEN
          LENG = (NVEC+1)*IDIM + 1
        ELSE
          LENG = (NVEC+1)*3 + 1
        ENDIF
        IF (NGRT .GT. 0) LENG = LENG + 3
*       NC = NR + 3 + LENG
        NC = NR + 5 + LENG

c Mike Potterfield 2/13/06
c For each GPS solution, write PROJID and then for each vector
c write the reject code to IUO3
      WRITE(IUO3) PROJID
	
**v 4.28***************	
        CALL NEWGRP (NVEC, IUNIT, G, NR, NC, FULL, ICM, NICM, KINDS,
*    &               ISNS, JSNS, LOBS, IOBS, IAUX, B, NVECS, IGRT, I83)
     &               ISNS, JSNS, LOBS, IOBS, IAUX, B, NVECS, IGRT, I83,
     *               ivc, IUO3) 
****************************

*** EXTRA MEMBER RECORDS

      ELSEIF (ID.EQ.'C'.OR.ID.EQ.'F') THEN
        WRITE (LUNIT,20) GCARD
   20   FORMAT ('0TOO MANY VECTOR RECORDS'/
     &          ' BAD GFILE STRUCTURE--',A80)
          write (lunit,*) 'subs2, 82'
        CALL ABORT2

*** CORRELATION AND COVARIANCE RECORDS

      ELSEIF (ID .EQ. 'D') THEN
*       IF (CCFLAG.AND.CVFLAG) THEN
        IF (CVFLAG) THEN
          WRITE(LUNIT,30) GCARD
   30     FORMAT('0BAD GFILE STRUCTURE'/
     *           ' CANT MIX D AND E RECORDS -- ',A80)
          write (lunit,*) 'subs2, 83'
          CALL ABORT2
        ELSE
          CALL LOADCR (G, GCARD, NR, NC)
**** 8-20-02	  
*         CCFLAG = .TRUE.
*         CVFLAG = .FALSE.
          CCFLAG = .TRUE.
********************************

        ENDIF
      ELSEIF (ID .EQ. 'E') THEN
** 8-20-02
*       IF (.NOT.CVFLAG.AND.CCFLAG) THEN
        IF (CCFLAG) THEN
******************	
          WRITE(LUNIT,30) GCARD
          write (lunit,*) 'subs2, 84'
          CALL ABORT2
        ELSE
          CALL LOADCV (G, GCARD, NR, NC)
**  8-20-02	  
*         CCFLAG = .TRUE.
*         CVFLAG = .TRUE.
          CVFLAG = .TRUE.
***************************	  
        ENDIF
      ENDIF

      GO TO 100

*** END OF PROCESSING--END OF FILE ENCOUNTERED

  777 IF (FULL) CALL OBSEQW (IUO, FULL, G, B, NX, NVEC, NR, NC, LENG,
     &                       NICM, ICM, KINDS, ISNS, JSNS, LOBS, IAUX,
     &                       IVF, FATAL, IOBS, IGRT, CVFLAG,
     &                       PROJID)

c Mike Potterfield 1/10/06
c The scratch file IUO3 is being used to identify rejected vectors
c for the observational summary in OBSSUM

      IF (LGF) THEN
        CALL LINE (2)
        WRITE (LUNIT,4)
    4   FORMAT ('0*********** END OF GPS OBSERVATIONS *********')
      ENDIF

      RETURN
      END
      SUBROUTINE SECOND (IUO, B, G, NX, LNWORK, IUO3)

c Mike Potterfield 2/10/06
c The new scratch file IUO3 is being used to tract rejected GPS
c vectors.

*** READ DATA AND WRITE FIRST OBSERVATION EQUATIONS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL FATAL
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      DIMENSION B(*), G(*), NX(*)
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7  ADJFIL, NAFILE
      CHARACTER*26 TBUF
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      COMMON/UNITS/ LUNIT

      IOBS = 0
      FATAL = .FALSE.

*** INITIALIZE CONNECTION MATRIX

      N = NUNK
      CALL OPENG (N, NX, LNWORK)

*** ADJUSTMENT FILE CONSTRAINTS

      IF (AFILE .EQ. 'NOAFILE') GOTO 1
      IUNIT = 11
      OPEN (IUNIT,ERR=1,STATUS='OLD',FILE=AFILE,IOSTAT=IOS)
      CALL SECADJ (IUNIT, IUO, IOBS, B, NX, FATAL)
      CLOSE (IUNIT)
    1 CONTINUE

*** BLUE-BOOK OBSERVATIONS

      ITYPE = 1
      CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM1)

      IUNIT = 12
      OPEN (IUNIT,STATUS='OLD',FILE=BBOOK)
      CALL SECBB (IUNIT, IUO, IOBS, B, NX, FATAL)
      CLOSE (IUNIT)

*** GFILE OBSERVATIONS

      IF ( .NOT. LGPS) THEN
        IUNIT = 13
        OPEN (IUNIT,ERR=2,STATUS='OLD',FILE=GFILE,IOSTAT=IOS)

c Mike Potterfield 2/10/06
c The scratch file IUO3 is being passed to SECGPS, and thence
c to NEWGRP, to keep track of rejected vectors

        CALL SECGPS (IUNIT, IUO, IUO3, IOBS, B, NX, G, FATAL)
        CLOSE (IUNIT)
      ENDIF
    2 CONTINUE

*** DFILE OBSERVATIONS

      IF ( .NOT. LDOP) THEN
        IUNIT = 14
        OPEN (IUNIT,ERR=4,STATUS='OLD',FILE=DFILE,IOSTAT=IOS)
        CALL SECDOP (IUNIT, IUO, IOBS, B, NX, FATAL)
        CLOSE (IUNIT)
      ENDIF
    4 CONTINUE

      CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM2)
      DIFM = DIFM2 - DIFM1
      CALL LINE (2)
      WRITE (LUNIT,55) DIFM
   55 FORMAT (/, ' MINUTES TO READ BBOOK, GFILE, AND DFILE', F7.1)

*** ABORT DUE TO LARGE MISCLOSURES

      IF (FATAL) THEN
        CALL LINE (3)
        WRITE (LUNIT,3) VM
    3   FORMAT ('0TERMINATED DUE TO MISCLOSURES (C-O)/SD EXCEEDING ',
     &          F7.1/)
          write (lunit,*) 'subs2, 85'
        CALL ABORT2
      ENDIF

*** REORDER THE UNKNOWNS

      ITYPE = 1
      CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM1)

      CALL REORDR (NX)

      CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM2)
      DIFM = DIFM2 - DIFM1
      CALL LINE (2)
      WRITE (LUNIT,65) DIFM
   65 FORMAT (/, ' MINUTES TO REORDER', F7.1)

*** GET THE COMPONENT LIST FOR THE UNKNOWNS

      CALL CMPNT (NX)

      RETURN
      END
      SUBROUTINE SECQQ (ACARD,NX)

*** ADD CONNECTIONS FOR AN ACCURACY RECORD

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      LOGICAL GETSSN
      LOGICAL ADDCON
      CHARACTER*80 ACARD
      DIMENSION NX(*)
      DIMENSION IC(LENC)
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) ISSN,JSSN
    1 FORMAT (10X,I4,36X,I4)

      IF (GETSSN(ISSN,ISN)) THEN
        IF (GETSSN(JSSN,JSN)) THEN
          NQQ = NQQ+1
          CALL CONEC (ISN,JSN,IC,LENG)
          IF ( .NOT. ADDCON(IC,LENG,NX)) THEN
            WRITE (LUNIT,666)
 666        FORMAT ('0INSUFFICIENT STORAGE FOR ACCURACY COMPUTATION')
          write (lunit,*) 'subs2, 86'
            CALL ABORT2
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END
      SUBROUTINE SECSS (ACARD, IUO, IOBS, B)

*** WRITE CONSTRAINT OBS. EQS. FOR AUXILIARY PARAMETERS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      LOGICAL GETPRM
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      CHARACTER*80 ACARD
      CHARACTER*1 TC1, TC2
      DIMENSION IC(LENC), C(LENC)
      DIMENSION B(*)
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON/UNITS/ LUNIT

*** EXTRACT INFORMATION FROM RECORD

      READ (ACARD,1) ICODE, IYR1, IMO1, IDY1, IHR1, IMN1, TC1,
     &                      IYR2, IMO2, IDY2, IHR2, IMN2, TC2, IVAL, ISD
    1 FORMAT (2X, I2, 2(I4, 4I2, A1), 2I5)
      IF (ICODE .EQ. 42) ICODE = 40
      IF (ICODE .EQ. 54) ICODE = 52
      IF (ICODE .GE. 26  .AND.  ICODE .LE. 29) ICODE = 25
      IF (ICODE .EQ. 25) THEN
        IF (LGPS) RETURN
      ELSEIF (ICODE .EQ. 40) THEN
        IF (LZEN) RETURN
      ELSEIF (ICODE .EQ. 52) THEN
        IF (LDIS) RETURN
      ELSEIF (ICODE .EQ. 99) THEN
        IF (LDOP) RETURN
      ELSE
        RETURN
      ENDIF

      IF (ACARD( 5: 8) .EQ. '    ') IYR1 = 1801
      IF (ACARD( 9:10) .EQ. '  ')   IMO1 = 1
      IF (ACARD(11:12) .EQ. '  ')   IDY1 = 1
      IF (ACARD(17:17) .EQ. ' ')    TC1 = 'Z'
      IF (ACARD(18:21) .EQ. '    ') IYR2 = 2099
      IF (ACARD(22:23) .EQ. '  ')   IMO2 = 12
      IF (ACARD(24:25) .EQ. '  ')   IDY2 = 31
      IF (ACARD(26:27) .EQ. '  ')   IHR2 = 23
      IF (ACARD(28:29) .EQ. '  ')   IMN2 = 59
      IF (ACARD(30:30) .EQ. ' ')    TC2 = 'Z'
      IF (ACARD(36:40) .EQ. '     ') ISD = 100

      CALL TOMNT (IYR1, IMO1, IDY1, IHR1, IMN1, TC1, IOLD)
      CALL TOMNT (IYR2, IMO2, IDY2, IHR2, IMN2, TC2, INEW)
      ITIME = (IOLD+INEW)/2

*** WRITE CONSTRAINT OBS. EQ.

      IF ( .NOT. GETPRM(ICODE, ITIME, IAUX) ) THEN
        CALL LINE (1)
        WRITE (LUNIT,2) ACARD
    2   FORMAT (1X, A80, ' *** RECORD NOT IN PARAMETER TABLE')
      ELSEIF (ACARD(31:40) .NE. '          ') THEN
        KIND = 0
        NCON = NCON+1
        IVF = 0
        IGRT = 0
        IOBS = IOBS+1
        OBSB = DBLE(IVAL)*1.0D-8
        SD = DBLE(ISD)*1.0D-8
        CALL FORMIC (KIND, IAUX, IDUMM1, IDUMM2, IC, LENG, IGRT)
        CALL FORMC (KIND, C, B, IAUX, IDUMM1, IDUMM2, IGRT)
        CALL COMPOB (KIND, OBS0, B, ADUMMY, IAUX, IDUMM1, IDUMM2, IGRT)
        CMO = OBS0-OBSB
        IF (IMODE .EQ. 0) CMO = 0.D0
        WRITE (IUO) KIND, IAUX, IDUMM1, IC, C, LENG, CMO, OBSB, SD,
     &              IOBS, IVF, IAUX, IGRT

*** ECHO OBSERVATION NUMBERS

        IF (LCS) THEN
          CALL LINE (1)
          WRITE (LUNIT,3) IOBS, ACARD
    3     FORMAT (1X, I6, 3X, A80)
        ENDIF
      ENDIF

      RETURN
      END
      SUBROUTINE SINGUL (ITER, B, ISING, GSING, LSING, STOL)

*** FATAL TERMINATE DUE TO SINGULARITY

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      PARAMETER ( MXSSN = 9999, MAXZ =  8000, MXPRM = 40, NSING = 200 )
      PARAMETER ( MXALL = MXPRM + 15 + MAXZ + MXSSN*3)
      LOGICAL GETSSN, INVZ, LOCSSN
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      CHARACTER*30 NAME, NAMES
      CHARACTER*1 CHR
      DIMENSION ISING(NSING), GSING(NSING)
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OBSUM/  NSUM(MXSSN,13), ICMPL(MXALL)
      COMMON/UNITS/ LUNIT

      CALL LINE (2)
      WRITE (LUNIT,1) LSING, STOL
    1 FORMAT ('0THE FOLLOWING', I5,' UNKNOWNS FALL BELOW TOLERANCE OF ',
     *        1PD8.2)
      IF ( .NOT. LOS) THEN
        CALL LINE (5)
        WRITE (LUNIT, 201)
  201   FORMAT (/, ' *** STATION CONNECTIONS FOR SINGULAR STATIONS ***',
     &          /, 8X, 'SSN',
     &             7X, 'DIR', 9X, 'ANG', 9X, 'AZI',
     &             9X, 'DIS', 9X, 'ZD', 10X, 'GPS', 7X, 'DOP',
     &          /, 8X, 'CMP',
     &             5X, 'FRM   TO', 4X, 'FRM   TO', 4X, 'FRM   TO',
     &             4X, 'FRM   TO', 4X, 'FRM   TO', 4X, 'FRM   TO', /)
      ENDIF

*** INITIALIZE SINGULARITY MAP

      IF (ITER .EQ. 0) THEN
        CALL INTMAP (B)
      ENDIF

      IF (LSING .GT. NSING) THEN
        WRITE (LUNIT,10)
   10   FORMAT (' TOO MANY SINGULARITIES -- CANNOT CONTINUE')
          write (lunit,*) 'subs2, 87'
        CALL ABORT2
      ENDIF

      CHR = ' '
      DO 4 I = 1, LSING
        CALL INVIUN (ISING(I), ISN, IUTYPE, IUCODE)

*** SINGULAR COORDINATES

        IF (IUCODE .EQ. 0) THEN
          NAME = NAMES(ISN)
          IF ( .NOT. LOCSSN(ISN,ISSN) ) THEN
            WRITE (LUNIT,34) ISN, ISSN
   34       FORMAT ('0ILLEGAL STATION SERIAL NUMBER IN SINGUL', 2I6)
          ELSEIF (IUTYPE .EQ. 1) THEN
            CALL LINE (1)
            IF (ITER .EQ. 0) CALL MAPSNG (I, ISN, B, CHR)
            WRITE (LUNIT,11) I, CHR, ISSN, NAME, GSING(I)
   11       FORMAT (1X, I3, A1, I6, 2X, A30, ' NORTH/SOUTH SHIFT GOOGE',
     &              ' NUMBER IS', 1PD10.2)
            IF ( .NOT. LOS) THEN
              IUN = IUNSTA(ISN, 1)
              ICMP = ICMPL(IUN)
              CALL LINE (1)
              WRITE (LUNIT, 205) ICMP,
     &           NSUM(ISN, 9), NSUM(ISN,10), NSUM(ISN,11), NSUM(ISN,12),
     &           NSUM(ISN, 1), NSUM(ISN, 2), NSUM(ISN, 3), NSUM(ISN, 4),
     &           NSUM(ISN, 5), NSUM(ISN, 6), NSUM(ISN, 7), NSUM(ISN, 8),
     &           NSUM(ISN,13)
  205         FORMAT (' ', 4X, I6, 4X, I4, I5, 5(I7, I5), I7 )
            ENDIF
          ELSEIF (IUTYPE .EQ. 2) THEN
            CALL LINE (1)
            IF (ITER .EQ. 0) CALL MAPSNG (I, ISN, B, CHR)
            WRITE (LUNIT,12) I, CHR, ISSN, NAME, GSING(I)
   12       FORMAT (1X, I3, A1, I6, 2X, A30, ' EAST/WEST   SHIFT GOOGE',
     &              ' NUMBER IS', 1PD10.2)
            IF ( .NOT. LOS) THEN
              IUN = IUNSTA(ISN, 2)
              ICMP = ICMPL(IUN)
              CALL LINE (1)
              WRITE (LUNIT, 205) ICMP,
     &           NSUM(ISN, 9), NSUM(ISN,10), NSUM(ISN,11), NSUM(ISN,12),
     &           NSUM(ISN, 1), NSUM(ISN, 2), NSUM(ISN, 3), NSUM(ISN, 4),
     &           NSUM(ISN, 5), NSUM(ISN, 6), NSUM(ISN, 7), NSUM(ISN, 8),
     &           NSUM(ISN,13)
            ENDIF
          ELSEIF (IUTYPE .EQ. 3) THEN
            CALL LINE (1)
            IF (ITER .EQ. 0) CALL MAPSNG (I, ISN, B, CHR)
            WRITE (LUNIT,13) I, CHR, ISSN, NAME, GSING(I)
   13       FORMAT (1X, I3, A1, I6, 2X, A30, ' UP/DOWN     SHIFT GOOGE',
     &              ' NUMBER IS', 1PD10.2)
            IF ( .NOT. LOS) THEN
              IUN = IUNSTA(ISN, 3)
              ICMP = ICMPL(IUN)
              CALL LINE (1)
              WRITE (LUNIT, 205) ICMP,
     &           NSUM(ISN, 9), NSUM(ISN,10), NSUM(ISN,11), NSUM(ISN,12),
     &           NSUM(ISN, 1), NSUM(ISN, 2), NSUM(ISN, 3), NSUM(ISN, 4),
     &           NSUM(ISN, 5), NSUM(ISN, 6), NSUM(ISN, 7), NSUM(ISN, 8),
     &           NSUM(ISN,13)
            ENDIF
          ELSE
            WRITE (LUNIT,3) IUTYPE
    3       FORMAT ('0ILLEGAL TYPE CODE', I5)
          write (lunit,*) 'subs2, 88'
            CALL ABORT2
          ENDIF

*** SINGULAR AUXILIARY PARAMETERS

        ELSEIF (IUCODE .EQ. 1) THEN
          IUN = IUNAUX(ISN)
          ICMP = ICMPL(IUN)
          CALL LINE (2)
          WRITE (LUNIT,21) I, ISN, GSING(I), ICMP
   21     FORMAT (1X, I3, ' ', 2X, 'THE', I6,
     &            '-TH AUX. PARM. GOOGE NUMBER IS', 1PD10.2,
     &            /, 16X, 'ITS COMPONENT # IS', I3)

*** SINGULAR AUXILIARY GPS AND DOPPLER ROTATION PARAMETERS

        ELSEIF (IUCODE .EQ. 2) THEN
          CALL LINE (2)
          IF (IUTYPE .EQ. 1) THEN
            IUN = IUNGRT(ISN, 1)
            ICMP = ICMPL(IUN)
            WRITE (LUNIT,41) I, ISN, GSING(I), ICMP
   41       FORMAT (1X, I3, ' ', 2X, 'THE X COMPONENT OF THE', I6,
     &              '-TH AUX. GPS ROT. PARM. GOOGE NUMBER IS', 1PD10.2,
     &            /, 16X, 'ITS COMPONENT # IS', I3)
          ELSEIF (IUTYPE .EQ. 2) THEN
            IUN = IUNGRT(ISN, 2)
            ICMP = ICMPL(IUN)
            WRITE (LUNIT,42) I, ISN, GSING(I)
   42       FORMAT (1X, I3, ' ', 2X, 'THE Y COMPONENT OF THE', I6,
     &              '-TH AUX. GPS ROT. PARM. GOOGE NUMBER IS', 1PD10.2,
     &            /, 16X, 'ITS COMPONENT # IS', I3)
          ELSEIF (IUTYPE .EQ. 3) THEN
            IUN = IUNGRT(ISN, 3)
            ICMP = ICMPL(IUN)
            WRITE (LUNIT,43) I, ISN, GSING(I)
   43       FORMAT (1X, I3, ' ', 2X, 'THE Z COMPONENT OF THE', I6,
     &              '-TH AUX. GPS ROT. PARM. GOOGE NUMBER IS', 1PD10.2,
     &            /, 16X, 'ITS COMPONENT # IS', I3)
          ENDIF

*** SINGULAR ROTATION PARAMETERS

        ELSEIF (IUCODE .EQ. 3) THEN
          IZ = ISN
          IF ( .NOT. INVZ(IZ, ISSN, ILIST) ) THEN
            WRITE (LUNIT,32) I, IZ
   32       FORMAT ('0ILLEGAL ROTATION UNKNOWN IN SINGUL', 2I6)
          write (lunit,*) 'subs2, 89'
            CALL ABORT2
          ELSEIF ( .NOT. GETSSN(ISSN, ISN) ) THEN
            WRITE (LUNIT,33) IZ, ISSN, ILIST
   33       FORMAT ('0ILLEGAL STATION SERIAL NUMBER IN SINGUL', 3I6)
          ELSE
            IUN = IUNROT(IZ)
            ICMP = ICMPL(IUN)
            NAME = NAMES(ISN)
            CALL LINE (1)
            IF (ITER .EQ. 0) CALL MAPSNG (I, ISN, B, CHR)
            WRITE (LUNIT,31) I, CHR, IZ, GSING(I), ILIST, ISSN, NAME
   31       FORMAT (1X, I3, A1, 2X, 'THE', I6, '-TH ROT. UNK. GOOGE',
     &              ' NUMBER IS', 1PD10.2, '  FOR LIST #',
     &              0P, I6, ' AT SSN', I5, ' [', A30, ']')
            IF ( .NOT. LOS) THEN
              CALL LINE (1)
              WRITE (LUNIT, 205) ICMP,
     &           NSUM(ISN, 9), NSUM(ISN,10), NSUM(ISN,11), NSUM(ISN,12),
     &           NSUM(ISN, 1), NSUM(ISN, 2), NSUM(ISN, 3), NSUM(ISN, 4),
     &           NSUM(ISN, 5), NSUM(ISN, 6), NSUM(ISN, 7), NSUM(ISN, 8),
     &           NSUM(ISN,13)
            ENDIF
          ENDIF

        ELSE
          WRITE (LUNIT,2) I, ISING(I), ISN, IUTYPE, IUCODE
    2     FORMAT ('0ILLEGAL IUCODE IN SINGUL', 5I6)
          write (lunit,*) 'subs2, 90'
          CALL ABORT2
        ENDIF
    4 CONTINUE

      WRITE (LUNIT,5)
    5 FORMAT ('0THE ADJUSTMENT IS SINGULAR !! ', 30('*') )

*** PRINT THE SINGULARITY MAP

      IF (ITER .EQ. 0) THEN
        CALL PRTMAP
      ENDIF

      IF (LABS) then
          write (lunit,*) 'subs2, 91'
        CALL ABORT2
      endif

      RETURN
      END
      SUBROUTINE STDDEV (BCARD,IRT,SD,REJECT)

*** DETERMINE THE STD DEV OF A OBSERVATION

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL REJECT
      CHARACTER*3 AM,AMM
      CHARACTER*4 ASS
      CHARACTER*80 BCARD
      COMMON /CONST/  PI, PI2, RAD

*** NOT FULLY IMPLEMENTED YET
*** STANDARD DEVIATIONS FOR DIRECTIONS, AZIMUTHS AND ZEN. DIST.

      IF (IRT .EQ. 20  .OR.  IRT .EQ. 22  .OR.  IRT .EQ. 40  .OR.
     &    IRT .EQ. 42  .OR.  IRT .EQ. 60) THEN
        IF (REJECT) THEN
          SD = PI2
        ELSEIF (BCARD(77:80) .NE. '    ') THEN
          READ (BCARD,1) ASS
  1       FORMAT (76X,A4)
          CALL NBLANK (ASS,4,IBLK)
          READ (ASS,2) SD
  2       FORMAT (F4.2)
          SD = SD/(3600.D0*RAD)
        ELSE
          SD = 5.D0/(3600.D0*RAD)
        ENDIF

*** STD DEV FOR REDUCED DISTANCES

      ELSEIF (IRT .EQ. 52) THEN
        IF (REJECT) THEN
          SD = 100.D0
        ELSEIF (BCARD(77:80) .NE. '   ') THEN
          READ (BCARD,3) ASS
  3       FORMAT (76X,A4)
          CALL NBLANK (ASS,4,IBLK)
          READ (ASS,4) SD
  4       FORMAT (F4.1)
          SD = SD/1000.D0
        ELSE
          SD = 0.05D0
        ENDIF

*** STD DEV FOR REDUCED LONG DISTANCES

      ELSEIF (IRT .EQ. 54) THEN
        IF (REJECT) THEN
          SD = 100.D0
        ELSEIF (BCARD(78:80) .NE. '   ') THEN
          READ (BCARD,5) AM
  5       FORMAT (77X,A3)
          CALL NBLANK (AM,3,IBLK)
          READ (AM,6) SD
  6       FORMAT (F3.2)
        ELSE
          SD = 0.5D0
        ENDIF

*** STD DEV FOR HORIZONTAL ANGLES

      ELSEIF (IRT .EQ. 30  .OR.  IRT .EQ. 32) THEN
        READ (BCARD,7) AMM
 7      FORMAT (68X,A3)
        CALL NBLANK (AMM,3,IBLK)
        IF (REJECT) THEN
          SD = PI2
        ELSEIF (IBLK .EQ. 3) THEN
          SD = 60.D0/(3600.D0*RAD)
        ELSEIF (IBLK .EQ. 2) THEN
          SD = 10.D0/(3600.D0*RAD)
        ELSE
          SD = 5.D0/(3600.D0*RAD)
        ENDIF

      ENDIF
      RETURN
      END
      SUBROUTINE STOREZ (ISN,JSN,B,IZ,OBSB,IGRT)

*** STORE THE ROTATION PARM

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /IZS/    IZ2

      I8 = 8
      CALL COMPOB (I8,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
      ROT = OBS0 - OBSB
 1    IF (ROT .GE. PI+PI) THEN
        ROT = ROT - PI - PI
        GO TO 1
      ENDIF
 2    IF (ROT .LT. 0.D0) THEN
        ROT = ROT + PI + PI
        GO TO 2
      ENDIF
      CALL PUTROT (ROT,IZ,B)
      IZ2 = IZ

      RETURN
      END
      SUBROUTINE TINYG (ISN, IUTYPE, IUCODE, GI)

*** VERY SMALL GOOGE NUMBER

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, MAXZ =  8000, MXPRM = 40 )
      PARAMETER ( MXALL = MXPRM + 15 + MAXZ + MXSSN*3)
      LOGICAL GETSSN, INVZ, LOCSSN
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      CHARACTER*30 NAME, NAMES
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /OBSUM/  NSUM(MXSSN,13), ICMPL(MXALL)
      COMMON/UNITS/ LUNIT

*** COORDINATES PARAMETERS WITH SMALL GOOGE NUMBERS

      IF (IUCODE .EQ. 0) THEN
        NAME = NAMES(ISN)
        IF ( .NOT. LOCSSN(ISN, ISSN) ) THEN
          WRITE (LUNIT,34) ISN, ISSN
   34     FORMAT ('0ILLEGAL STATION SERIAL NUMBER IN TINYG', 2I6)
        ELSEIF (IUTYPE .EQ. 1) THEN
          CALL LINE (1)
          WRITE (LUNIT,11) ISSN, NAME, GI
   11     FORMAT (6X, I6, 2X, A30, ' NORTH/SOUTH SHIFT GOOGE',
     &            ' NUMBER IS', 1PD10.2)
          IF ( .NOT. LOS) THEN
            IUN = IUNSTA(ISN, 1)
            ICMP = ICMPL(IUN)
            CALL LINE (1)
            WRITE (LUNIT, 205) ICMP,
     &           NSUM(ISN, 9), NSUM(ISN,10), NSUM(ISN,11), NSUM(ISN,12),
     &           NSUM(ISN, 1), NSUM(ISN, 2), NSUM(ISN, 3), NSUM(ISN, 4),
     &           NSUM(ISN, 5), NSUM(ISN, 6), NSUM(ISN, 7), NSUM(ISN, 8),
     &           NSUM(ISN,13)
  205       FORMAT (' ', 5X, I6, 4X, I4, I5, 5(I7, I5), I7 )
          ENDIF
        ELSEIF (IUTYPE .EQ. 2) THEN
          CALL LINE (1)
          WRITE (LUNIT,12) ISSN, NAME, GI
   12     FORMAT (6X, I6, 2X, A30, ' EAST/WEST   SHIFT GOOGE',
     &            ' NUMBER IS', 1PD10.2)
          IF ( .NOT. LOS) THEN
            IUN = IUNSTA(ISN, 2)
            ICMP = ICMPL(IUN)
            CALL LINE (1)
            WRITE (LUNIT, 205) ICMP,
     &           NSUM(ISN, 9), NSUM(ISN,10), NSUM(ISN,11), NSUM(ISN,12),
     &           NSUM(ISN, 1), NSUM(ISN, 2), NSUM(ISN, 3), NSUM(ISN, 4),
     &           NSUM(ISN, 5), NSUM(ISN, 6), NSUM(ISN, 7), NSUM(ISN, 8),
     &           NSUM(ISN,13)
          ENDIF
        ELSEIF (IUTYPE .EQ. 3) THEN
          CALL LINE (1)
          WRITE (LUNIT,13) ISSN, NAME, GI
   13     FORMAT (6X, I6, 2X, A30, ' UP/DOWN     SHIFT GOOGE',
     &            ' NUMBER IS', 1PD10.2)
          IF ( .NOT. LOS) THEN
            IUN = IUNSTA(ISN, 3)
            ICMP = ICMPL(IUN)
            CALL LINE (1)
            WRITE (LUNIT, 205) ICMP,
     &           NSUM(ISN, 9), NSUM(ISN,10), NSUM(ISN,11), NSUM(ISN,12),
     &           NSUM(ISN, 1), NSUM(ISN, 2), NSUM(ISN, 3), NSUM(ISN, 4),
     &           NSUM(ISN, 5), NSUM(ISN, 6), NSUM(ISN, 7), NSUM(ISN, 8),
     &           NSUM(ISN,13)
          ENDIF
        ELSE
          WRITE (LUNIT,3) IUTYPE, IUCODE
    3     FORMAT ('0ILLEGAL TYPE CODE', I5, ' FOR IUCODE', I5)
          write (lunit,*) 'subs2, 92'
          CALL ABORT2
        ENDIF

*** AUXILIARY PARAMETERS WITH SMALL GOOGE NUMBERS

      ELSEIF (IUCODE .EQ. 1) THEN
        IUN = IUNAUX(ISN)
        ICMP = ICMPL(IUN)
        CALL LINE (2)
        WRITE (LUNIT,21) ISN, GI, ICMP
   21   FORMAT (8X, 'THE', I6, '-TH AUX. PARM. GOOGE NUMBER IS',1PD10.2,
     &          /, 17X, 'ITS COMPONENT # IS', I3)

*** AUXILIARY GPS ROTATION PARAMETERS WITH SMALL GOOGE NUMBERS

      ELSEIF (IUCODE .EQ. 2) THEN
        CALL LINE (2)
        IF (IUTYPE .EQ. 1) THEN
          IUN = IUNGRT(ISN, 1)
          ICMP = ICMPL(IUN)
          WRITE (LUNIT,41) ISN, GI, ICMP
   41     FORMAT (8X, 'THE X COMPONENT OF THE', I6,
     &            '-TH AUX. GPS ROT. PARM. GOOGE NUMBER IS', 1PD10.2,
     &            /, 17X, 'ITS COMPONENT # IS', I3)
        ELSEIF (IUTYPE .EQ. 2) THEN
          IUN = IUNGRT(ISN, 2)
          ICMP = ICMPL(IUN)
          WRITE (LUNIT,42) ISN, GI, ICMP
   42     FORMAT (8X, 'THE Y COMPONENT OF THE', I6,
     &            '-TH AUX. GPS ROT. PARM. GOOGE NUMBER IS', 1PD10.2,
     &            /, 17X, 'ITS COMPONENT # IS', I3)
        ELSEIF (IUTYPE .EQ. 3) THEN
          IUN = IUNGRT(ISN, 3)
          ICMP = ICMPL(IUN)
          WRITE (LUNIT,43) ISN, GI, ICMP
   43     FORMAT (8X, 'THE Z COMPONENT OF THE', I6,
     &            '-TH AUX. GPS ROT. PARM. GOOGE NUMBER IS', 1PD10.2,
     &            /, 17X, 'ITS COMPONENT # IS', I3)
        ELSE
          WRITE (LUNIT,3) IUTYPE, IUCODE
          write (lunit,*) 'subs2, 93'
          CALL ABORT2
        ENDIF

*** ROTATION PARAMETERS WITH SMALL GOOGE NUMBERS

      ELSEIF (IUCODE .EQ. 3) THEN
        IZ = ISN
        IF ( .NOT. INVZ(IZ, ISSN, ILIST) ) THEN
          WRITE (LUNIT,32) IZ
   32     FORMAT ('0ILLEGAL ROTATION UNKNOWN IN TINYG', 2I6)
          write (lunit,*) 'subs2, 94'
          CALL ABORT2
        ELSEIF ( .NOT. GETSSN(ISSN, ISN) ) THEN
          WRITE (LUNIT,33) IZ, ISSN, ILIST
   33     FORMAT ('0ILLEGAL STATION SERIAL NUMBER IN TINYG', 3I6)
        ELSE
          NAME = NAMES(ISN)
          CALL LINE (1)
          WRITE (LUNIT,31) IZ, GI, ILIST, NAME
   31     FORMAT (8X, 'THE', I6, '-TH ROT. UNK. GOOGE NUMBER',
     &            ' IS', 1PD10.2, '  FOR LIST #', 0P, I6, ' AT ', A30)
          IF ( .NOT. LOS) THEN
            IUN = IUNROT(IZ)
            ICMP = ICMPL(IUN)
            CALL LINE (1)
            WRITE (LUNIT, 205) ICMP,
     &           NSUM(ISN, 9), NSUM(ISN,10), NSUM(ISN,11), NSUM(ISN,12),
     &           NSUM(ISN, 1), NSUM(ISN, 2), NSUM(ISN, 3), NSUM(ISN, 4),
     &           NSUM(ISN, 5), NSUM(ISN, 6), NSUM(ISN, 7), NSUM(ISN, 8),
     &           NSUM(ISN,13)
          ENDIF
        ENDIF

      ELSE
        WRITE (LUNIT,2) GI, ISN, IUTYPE, IUCODE
    2   FORMAT ('0ILLEGAL IUCODE IN TINYG', 4I6)
          write (lunit,*) 'subs2, 95'
        CALL ABORT2
      ENDIF

      RETURN
      END
      SUBROUTINE TO83 (DX, DY, DZ, ISN, I83, B)

*** CONVERT GPS VECTOR TO NAD83

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      LOGICAL L2HLF, LEHT
      DIMENSION B(*)
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT

*** GET 83 AND CONVERT TO SATELLITE SYSTEM

      IF (L2HLF) THEN
        CALL GETEHX (X, ISN, B)
        CALL GETEHY (Y, ISN, B)
        CALL GETEHZ (Z, ISN, B)
      ELSE
        CALL GETECX (X, ISN, B)
        CALL GETECY (Y, ISN, B)
        CALL GETECZ (Z, ISN, B)
      ENDIF
      CALL FROM83 (I83, X, Y, Z, XI, YI, ZI)

*** ADD VECTOR IN SATELLITE SYSTEM

      XJ = XI + DX
      YJ = YI + DY
      ZJ = ZI + DZ

*** CONVERT ENDPOINTS BACK TO 83

      CALL GET83 (I83, XI, YI, ZI, X1, Y1, Z1)
      CALL GET83 (I83, XJ, YJ, ZJ, X2, Y2, Z2)

*** GET VECTOR IN 83

      DX = X2 - X1
      DY = Y2 - Y1
      DZ = Z2 - Z1

      RETURN
      END
      SUBROUTINE TOGEO2 (X, Y, Z, GLAT, GLON)

*** CONVERT X, Y, Z INTO GEODETIC LAT, LON
*** REF: EQ A.4B, P. 132, APPENDIX A, OSU #370
*** REF: GEOM GEOD NOTES GS 658, RAPP

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (MAXINT = 10, TOL = 1.D-13)
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON/UNITS/ LUNIT

      AE2 = AX*E2

*** COMPUTE INITIAL ESTIMATE OF REDUCED LATITUDE  (EHT=0)

      P = DSQRT( X*X + Y*Y )
      ICOUNT = 0
      TGLA = Z/P/(1.D0 - E2)

*** ITERATE TO CONVERGENCE, OR TO MAX # ITERATIONS

    1 IF (ICOUNT .LE. MAXINT) THEN
        TGLAX = TGLA
        TGLA = Z/(  P - ( AE2/DSQRT( 1.D0 + (1.D0 - E2)*TGLA*TGLA ) )  )
        ICOUNT = ICOUNT + 1
        IF ( DABS(TGLA - TGLAX) .GT. TOL ) GO TO 1

*** CONVERGENCE ACHIEVED

        GLAT = DATAN(TGLA)
        GLON = DATAN2(Y, X)

*** TOO MANY ITERATIONS

      ELSE
        WRITE(LUNIT,2)
    2   FORMAT (' FAILURE TO CONVERGE IN TOGEO2 -- FATAL')
        STOP
      ENDIF

      RETURN
      END
      SUBROUTINE TOLGH( X, Y, Z, N, E, U, ISN, JSN, B)

*** CONVERT RESIDUALS TO LGH

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION B(*)
      DOUBLE PRECISION X, Y, Z, N, E, U
      SAVE GLA1, GLO1, GLA2, GLO2
      SAVE SB, CB, SL, CL
      SAVE ISNX, JSNX
      DATA ISNX /0/, JSNX /0/

      IF (ISN .NE. ISNX  .OR.  JSN .NE. JSNX) THEN
        IF (ISN .NE. ISNX) THEN
          CALL GETGLA (GLA1, ISN, B)
          CALL GETGLO (GLO1, ISN, B)
          ISNX = ISN
        ENDIF
        IF (JSN .NE. JSNX) THEN
          CALL GETGLA (GLA2, JSN, B)
          CALL GETGLO (GLO2, JSN, B)
          JSNX = JSN
        ENDIF

*** COMPUTE MEAN ORIENTATION

        GLA = (GLA1 + GLA2)/2.D0
        GLO = (GLO1 + GLO2)/2.D0

        SB = DSIN(GLA)
        CB = DCOS(GLA)
        SL = DSIN(GLO)
        CL = DCOS(GLO)
      ENDIF

      N = -SB*CL*X - SB*SL*Y + CB*Z
      E = -   SL*X +    CL*Y
      U =  CB*CL*X + CB*SL*Y + SB*Z

      RETURN
      END
      SUBROUTINE TOMNT (IYEAR,IMO,IDY,IHR,IMN,TC,MINS)

*** CONVERT TIME TO MINUTES

*** CAUTION:  ONLY STRICTLY CORRECT FOR RECENT TIMES
***           FOR DEFINITIVE CONVERSION,USE MODIFIED JULIAN DAY ROUTINE

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL LEAP
      CHARACTER*1 TC
      DIMENSION IDAYS(12),LDAYS(12)

***   THESE ARRAYS SHOULD BE PARAMETERS AS THEY ARE FIXED VALUES

      DATA IDAYS /0,31,59,90,120,151,181,212,243,273,304,334/
      DATA LDAYS /0,31,60,91,121,152,182,213,244,274,305,335/

      IF (LEAP(IYEAR)) THEN
        IDAY = LDAYS(IMO)+IDY
      ELSE
        IDAY = IDAYS(IMO)+IDY
      ENDIF
      ID = (IYEAR-1900)*365.25+IDAY
      IHRS = IHR+ITCODE(TC)
      IH = ID*24+IHRS
      MINS = IH*60+IMN

      RETURN
      END
C--------------------------------------------------------------------------------------
      SUBROUTINE TOXYZ (GLAT,GLON,EHT,X,Y,Z)

*** COMPUTE X,Y,Z
*** REF P.17 GEOMETRIC GEODESY NOTES VOL 1, OSU, RAPP

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /OPT/  AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &              LMSL, LSS, LUP, LABS, LADJ, LUPI

      SLAT = DSIN(GLAT)
      CLAT = DCOS(GLAT)
      W = DSQRT(1.D0-E2*SLAT*SLAT)
      EN = AX/W

      X = (EN+EHT)*CLAT*DCOS(GLON)
      Y = (EN+EHT)*CLAT*DSIN(GLON)
      Z = (EN*(1.D0-E2)+EHT)*SLAT

      RETURN
      END
C------------------------------------------------------------------------------------
      SUBROUTINE UP80 (BCARD, B)

*** UPDATE THE *80* CONTROL PT RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      CHARACTER*80 BCARD
      CHARACTER*1 ADIR1, ADIR2
      CHARACTER*2 AD1, AM1, AM2
      CHARACTER*3 AD2
      CHARACTER*6 AMSL
      CHARACTER*6 AEHT
      CHARACTER*7 AS1, AS2
      LOGICAL ELFLAG, DFFLAG, GETSSN
      LOGICAL L2HLF, LEHT
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /FLAGS/  ELFLAG(MXSSN), DFFLAG(MXSSN)
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

**v6.3

      logical    overwrite_const_coord_in_newbb
      character  CC_records*80
      character  lat_lon_char*25
      logical    un_constrained_Hz
      common/MM6/overwrite_const_coord_in_newbb
      common/MM6_array/CC_records(MXSSN)

**so far for v6.3

      READ (BCARD,1) ISSN
    1 FORMAT (10X, I4)

*** UPDATE LAT. AND LONG.

      IF ( GETSSN(ISSN, ISN) ) THEN
        IF (IDIM .NE. 1) THEN
**v6.3
          if (.not. overwrite_const_coord_in_newbb) then
**so far for v6.3
            CALL GETGLA (GLAT, ISN, B)
            CALL GETGLO (GLON, ISN, B)
            CALL GETDMS (GLAT, ID1, IM1, S1, ISIGN)
            IF (ISIGN .GT. 0) THEN
              ADIR1 = 'N'
            ELSE
              ADIR1 = 'S'
            ENDIF
            CALL GETDMS (GLON, ID2, IM2, S2, ISIGN)
            IF (ISIGN .GT. 0) THEN
              ALONGD=ID2+(IM2/60.D0)+(S2/3600.D0)
              ALONGD=360.D0-ALONGD
              ID2=INT(ALONGD)
              ALONGM=(ALONGD-ID2) * 60.D0
              IM2=INT(ALONGM)
              S2=(ALONGM-IM2) * 60.D0
              ADIR2 = 'W'
            ELSE
              ADIR2 = 'W'
            ENDIF
            IS1 = IDNINT(S1*100000.D0)
            IS2 = IDNINT(S2*100000.D0)

            WRITE (AD1,4) ID1
            WRITE (AM1,4) IM1
            WRITE (AS1,2) IS1
            WRITE (AD2,3) ID2
            WRITE (AM2,4) IM2
            WRITE (AS2,2) IS2
    4       FORMAT (I2.2)
    2       FORMAT (I7.7)
    3       FORMAT (I3.3)

            BCARD(45:46) = AD1
            BCARD(47:48) = AM1
            BCARD(49:55) = AS1
            BCARD(56:56) = ADIR1
            BCARD(57:59) = AD2
            BCARD(60:61) = AM2
            BCARD(62:68) = AS2
            BCARD(69:69) = ADIR2
**v6.3

          else
            do ii=1,MXSSN
              read (CC_records(ii)(11:14),'(i4)') iissn
              if (ISSN == iissn) then               
                read(CC_records(ii)(45:69),'(a25)') lat_lon_char
                un_constrained_Hz = (lat_lon_char == ' ')
                if (.not. un_constrained_Hz) then
                  BCARD(45:69) = CC_records(ii)(45:69)
                  exit
                endif
              endif
            enddo
          endif

**so far for v6.3

        ENDIF

*** UPDATE HEIGHT

        IF (ELFLAG(ISN)  .AND.  IDIM .NE. 2) THEN
          CALL GETMSL (GMSL, ISN, B)
          MSL = IDNINT(GMSL*100.D0)
          WRITE (AMSL,5) MSL
    5     FORMAT (I6.3)
          BCARD(70:75) = AMSL
        ENDIF

**v 4.28 *******************************
*** IF 2.5DIM ADJUSTMENT AND ELLIPSOIDAL HEIGHTS ARE READ FROM THE
*** EXTENDED *80* RECORDS THEN UPDATE ELLIPSOIDAL HEIGHTS - BIGADJUST

*       IF (IDIM .EQ. 2  .AND.  L2HLF  .AND.  LEHT) THEN
*         IF (I2HLF(ISN) .EQ. 1) THEN
*           CALL GETEHT (EHT, ISN, B)
*           IEHT = IDNINT(EHT*100.D0)
*           WRITE (AEHT,7) IEHT
*   7       FORMAT (I6.3)
*           BCARD(83:88) = AEHT
*         ENDIF
*       ENDIF

      ELSE
        WRITE (LUNIT,666) BCARD
  666   FORMAT ('0ILLEGAL ISN IN UP80--', A80)
          write (lunit,*) 'subs2, 95'
        CALL ABORT2
      ENDIF

      RETURN
      END
C-----------------------------------------------------------------------------------------------
      SUBROUTINE UP84 (BCARD, B)

*** UPDATE *84* GEOID HEIGHT RECORDS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      CHARACTER*80 BCARD
      CHARACTER*6 AGH
      LOGICAL ELFLAG, DFFLAG, GETSSN
      DIMENSION B(*)
      COMMON /FLAGS/  ELFLAG(MXSSN), DFFLAG(MXSSN)
      COMMON/UNITS/ LUNIT

      READ (BCARD,1) ISSN
    1 FORMAT (10X, I4)
      IF ( GETSSN(ISSN, ISN) ) THEN
        IF ( .NOT. ELFLAG(ISN) ) THEN
          CALL GETGH (GGH, ISN, B)
          IGH = IDNINT(GGH*1.0D2)
          WRITE (AGH,5) IGH
    5     FORMAT (I6.2)
          BCARD(70:75) = AGH
        ENDIF
      ELSE
        WRITE (LUNIT,666) BCARD
  666   FORMAT ('0ILLEGAL ISN IN UP84--', A80)
          write (lunit,*) 'subs2, 96'
        CALL ABORT2
      ENDIF

      RETURN
      END
C---------------------------------------------------------------------------------------------
      SUBROUTINE UP86 (BCARD, B)

*** UPDATE *86* ELEVATION RECORDS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      CHARACTER*80 BCARD
      CHARACTER*7 AGH, AMSL, AEHT
      LOGICAL L2HLF, LEHT
      LOGICAL ELFLAG, DFFLAG, GETSSN
      DIMENSION B(*)
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /FLAGS/  ELFLAG(MXSSN), DFFLAG(MXSSN)
      COMMON/UNITS/ LUNIT

**V6.3

      logical    overwrite_const_coord_in_newbb
      character  HT_char*7
      logical    un_constrained_HT
      character  CC_records*80
      common/MM6/overwrite_const_coord_in_newbb
      common/MM6_array/CC_records(MXSSN)

**so far for v6.3

      READ (BCARD,1) ISSN
    1 FORMAT (10X, I4)
      IF (GETSSN(ISSN, ISN) ) THEN
**v6.3
        if (.not. overwrite_const_coord_in_newbb) then
**so far for v6.3
          CALL GETMSL (GMSL0, ISN, B)
          CALL GETGH (GHT0, ISN, B)
          IF (L2HLF  .AND.  I2HLF(ISN) .EQ. 1) THEN
            CALL GETEHT (EHT, ISN, B)
            IF ( ELFLAG(ISN) ) THEN
              GMSL = EHT - GHT0
              GHT = GHT0
            ELSE
              GHT = EHT - GMSL0
              GMSL = GMSL0
            ENDIF
          ELSE
            GMSL = GMSL0
            GHT = GHT0
            EHT = GMSL + GHT
          ENDIF
  
          MSL = IDNINT(GMSL*1000.D0)
          WRITE (AMSL,3) MSL
  3       FORMAT (I7.3)
          BCARD(17:23) = AMSL
          IGH = IDNINT(GHT*1000.D0)
          WRITE (AGH,3) IGH
          BCARD(36:42) = AGH
          IEHT = IDNINT(EHT*1000.D0)
          WRITE (AEHT,3) IEHT
          BCARD(46:52) = AEHT
**v6.3

        else
          do ii=1,MXSSN
            read (CC_records(ii)(11:14),'(i4)') iissn
            if (ISSN == iissn) then               
              read(CC_records(ii)(70:76),'(a7)') HT_char
              un_constrained_HT = (HT_char == ' ')
              if(.not. un_constrained_HT) then          !The station is constrained includin the height
                BCARD(46:52) = CC_records(ii)(70:76)
                exit
              endif
            endif
          enddo
        endif

**so far for v6.3

      ELSE
        WRITE (LUNIT,666) BCARD
  666   FORMAT ('0ILLEGAL ISN IN UP86--', A80)
        write (lunit,*) 'subs2, 97'
        CALL ABORT2
      ENDIF

      RETURN
      END
C-----------------------------------------------------------------------------------------------
      SUBROUTINE UPAFIL (B, ITER)

*** UPDATE AFILE PARAMETER RECORDS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      CHARACTER*80 ACARD
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7  ADJFIL, NAFILE
      CHARACTER*2 CC12
      DIMENSION B(*)
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

*** OPEN OLD AFILE

      IOLD = 14
      OPEN (IOLD,ERR=666,STATUS='OLD',FILE=AFILE)

*** OPEN NEW AFILE

      INEW = 15
      OPEN (INEW,ERR=667,STATUS='UNKNOWN',FILE=NAFILE)

*** READ THE AFILE RECORDS

 100  READ (IOLD,1,END=668) ACARD
 1    FORMAT (A80)
      READ (ACARD,4) CC12
 4    FORMAT (A2)

      IF (CC12 .EQ. 'RR') THEN
        CALL UPRR (ACARD)
      ELSEIF (CC12 .EQ. 'SS') THEN
        CALL UPSS (ACARD, B)
      ELSEIF (CC12 .EQ. 'VV') THEN
        CALL UPVV (ACARD)
      ELSEIF (CC12 .EQ. 'II') THEN
        CALL UPII (ACARD, ITER)
      ELSEIF (CC12 .EQ. 'DD') THEN
        CALL UPDD (ACARD)
      ENDIF

      WRITE (INEW,2) ACARD
 2    FORMAT (A80)
      GO TO 100

 668  CLOSE (IOLD)
      CLOSE (INEW)
      CALL LINE (2)
      WRITE (LUNIT,5) NAFILE
  5   FORMAT ('0UPDATED CONTROL POINT RECORDS IN FILE -- ', A7)
      RETURN

*** NO OLD AFILE FOUND

 666  WRITE (LUNIT,699)
 699  FORMAT ('0NO OLD AFILE FOUND'/)
          write (lunit,*) 'subs2, 98'
      CALL ABORT2
      RETURN

*** NOT ABLE TO OPEN NEW AFILE

 667  WRITE (LUNIT,698)
 698  FORMAT ('0NOT ABLE TO OPEN NEW AFILE FILE'/)
          write (lunit,*) 'subs2, 99'
      CALL ABORT2

      RETURN
      END
      SUBROUTINE UPDD (ACARD)

*** UPDATE DIMENSIONALITY RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      CHARACTER*80 ACARD
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      LOGICAL L2HLF, LEHT
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT

*** UPDATE ONLY IF A 2.5 DIMENSIONAL ADJUSTMENT

      IF (IDIM .NE. 2) RETURN
      IF (.NOT. L2HLF) RETURN

*** IF 2.5-DIM ADJUSTMENT MAKE SURE THAT ELLIPIPSOID HEIGHTS WILL
*** BE READ FROM 80 RECS

** v 4.28
*     IF (ACARD(6:6) .NE. 'E' ) THEN
*       ACARD(6:6) = 'E'
*       LEHT = .TRUE.
*     ENDIF
**************

      RETURN
      END
      SUBROUTINE UPII (ACARD, ITER)

*** UPDATE ITERATION RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI

      READ (ACARD,1) ITMAXX
    1 FORMAT (2X, I2)
      IF (ACARD(3:4) .NE. '  ') THEN
         ITMAX = IABS(ITMAXX) - (ITER + 1)
      ELSE
         ITMAX = 5 - (ITER + 1)
      ENDIF
      WRITE (ACARD(3:4), '(I2)') ITMAX
      RETURN
      END
      SUBROUTINE UPRR (ACARD)

*** UPDATE AUXILIARY GPS AND DOPPLER ROTATION PARAMETER RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      CHARACTER*2 ACODE
      CHARACTER*1 TC1, TC2
      LOGICAL GETGRT
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /CONST/  PI, PI2, RAD
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) ACODE, IYR1, IMO1, IDY1, IHR1, IMN1, TC1,
     &                      IYR2, IMO2, IDY2, IHR2, IMN2, TC2
    1 FORMAT (2X, A2, 2(I4, 4I2, A1) )

      IF (ACODE .EQ. '  ') THEN
        ICODE = 25
      ELSE
        READ (ACODE,5) ICODE
    5   FORMAT (I2)
      ENDIF

      IF (ICODE .EQ. 25) THEN
        IF (LGPS) RETURN
      ELSEIF (ICODE .EQ. 99) THEN
        IF (LDOP) RETURN
      ELSE
        CALL LINE (3)
        WRITE (LUNIT,2) ICODE
    2   FORMAT ('0ILLEGAL OBSERVATIONAL TYPE', I4/)
        RETURN
      ENDIF


*** SET DEFAULTS

      IF (ACARD( 5: 8) .EQ. '    ') IYR1 = 1801
      IF (ACARD( 9:10) .EQ. '  ')   IMO1 = 1
      IF (ACARD(11:12) .EQ. '  ')   IDY1 = 1
      IF (ACARD(17:17) .EQ. ' ')    TC1 = 'Z'
      IF (ACARD(18:21) .EQ. '    ') IYR2 = 2099
      IF (ACARD(22:23) .EQ. '  ')   IMO2 = 12
      IF (ACARD(24:25) .EQ. '  ')   IDY2 = 31
      IF (ACARD(26:27) .EQ. '  ')   IHR2 = 23
      IF (ACARD(28:29) .EQ. '  ')   IMN2 = 59
      IF (ACARD(30:30) .EQ. ' ')    TC2 = 'Z'

      CALL TOMNT (IYR1, IMO1, IDY1, IHR1, IMN1, TC1, IOLD)
      CALL TOMNT (IYR2, IMO2, IDY2, IHR2, IMN2, TC2, INEW)
      ITIME = (IOLD + INEW)/2

      IF ( .NOT. GETGRT(ICODE, ITIME, IGRT) ) THEN
        CALL LINE (1)
        WRITE (LUNIT,4) ACARD
    4   FORMAT (1X, A80, ' *** RECORD NOT IN PARAMETER TABLE')
      ELSE
        CALL GETRTG (VALX, VALY, VALZ, IGRT)

***   CONVERT GPS AND DOPPLER ROT. VALUES TO 1.D-5 ARC SEC. FROM RADIANS

        SECRAD = PI / ( 60.D0 * 60.D0 * 180.D0 * 1.D5 )
        IVALX = IDNINT( VALX / SECRAD )
        IVALY = IDNINT( VALY / SECRAD )
        IVALZ = IDNINT( VALZ / SECRAD )

*** CHECK FOR SIZE OF NUMBERS
*** ( MAXIMUM INTEGER*4 IS 2147483647 )

        IMAX = 2147483647
        IF (IVALX .GT. IMAX) THEN
          IVALX = IMAX
          WRITE (LUNIT,6) VALX, ACARD
    6     FORMAT (' ***** WARNING, GPS ROTATION PARAMETER VALUE =',
     &            F14.9, /, ' RESET -- TOO LARGE FOR RECORD', /, ' ',
     &            A80)
        ELSEIF (IVALX .LT. -999999999) THEN
          IVALX = -999999999
          WRITE (LUNIT,6) VALX, ACARD
        ENDIF
        IF (IVALY .GT. IMAX) THEN
          IVALY = IMAX
          WRITE (LUNIT,6) VALY, ACARD
        ELSEIF (IVALY .LT. -999999999) THEN
          IVALY = -999999999
          WRITE (LUNIT,6) VALY, ACARD
        ENDIF
        IF (IVALZ .GT. IMAX) THEN
          IVALZ = IMAX
          WRITE (LUNIT,6) VALZ, ACARD
        ELSEIF (IVALZ .LT. -999999999) THEN
          IVALZ = -999999999
          WRITE (LUNIT,6) VALZ, ACARD
        ENDIF

        WRITE (ACARD(31:60),3) IVALX, IVALY, IVALZ
    3   FORMAT (3I10)
      ENDIF

      RETURN
      END
      SUBROUTINE UPSS (ACARD, B)

*** UPDATE AUXILIARY PARAMETER RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      CHARACTER*15 PNAME
      LOGICAL GETPRM
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      CHARACTER*1 TC1, TC2
      DIMENSION B(*)
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) ICODE, IYR1, IMO1, IDY1, IHR1, IMN1, TC1,
     &                      IYR2, IMO2, IDY2, IHR2, IMN2, TC2,
     &               PNAME
    1 FORMAT (2X, I2, 2(I4, 4I2, A1), 10X, A15)
      IF (ICODE .EQ. 42) ICODE = 40
      IF (ICODE .EQ. 54) ICODE = 52
      IF (ICODE .GE. 26  .AND.  ICODE .LE. 29) ICODE = 25

      IF (ICODE .EQ. 25) THEN
        IF (LGPS) RETURN
      ELSEIF (ICODE .EQ. 40) THEN
        IF (LZEN) RETURN
      ELSEIF (ICODE .EQ. 52) THEN
        IF (LDIS) RETURN
      ELSEIF (ICODE .EQ. 99) THEN
        IF (LDOP) RETURN
      ELSE
        CALL LINE (3)
        WRITE (LUNIT,2) ICODE
 2      FORMAT ('0ILLEGAL OBSERVATIONAL TYPE', I4/)
        RETURN
      ENDIF

*** SET DEFAULTS

      IF (ACARD( 5: 8) .EQ. '    ') IYR1 = 1801
      IF (ACARD( 9:10) .EQ. '  ')   IMO1 = 1
      IF (ACARD(11:12) .EQ. '  ')   IDY1 = 1
      IF (ACARD(17:17) .EQ. ' ')    TC1 = 'Z'
      IF (ACARD(18:21) .EQ. '    ') IYR2 = 2099
      IF (ACARD(22:23) .EQ. '  ')   IMO2 = 12
      IF (ACARD(24:25) .EQ. '  ')   IDY2 = 31
      IF (ACARD(26:27) .EQ. '  ')   IHR2 = 23
      IF (ACARD(28:29) .EQ. '  ')   IMN2 = 59
      IF (ACARD(30:30) .EQ. ' ')    TC2 = 'Z'

      CALL TOMNT (IYR1, IMO1, IDY1, IHR1, IMN1, TC1, IOLD)
      CALL TOMNT (IYR2, IMO2, IDY2, IHR2, IMN2, TC2, INEW)
      ITIME = (IOLD + INEW)/2

      IF ( .NOT. GETPRM(ICODE, ITIME, IAUX) ) THEN
        CALL LINE (1)
        WRITE (LUNIT,4) ACARD
    4   FORMAT (1X, A80, ' *** RECORD NOT IN PARAMETER TABLE')
      ELSE
        CALL GETAUX (VAL, IAUX, B)
        IVAL = IDNINT( VAL/1.0D-8 )

*** CHECK FOR SIZE OF NUMBERS

        IF (IVAL .GT. 99999) THEN
          IVAL = 99999
          WRITE (LUNIT,6) VAL, ACARD
    6     FORMAT (' ***** WARNING, AUXILIARY PARAMETER VALUE =', F14.9,
     &            /, ' RESET -- TOO LARGE FOR RECORD', /, ' ',
     &            A80)
        ELSEIF (IVAL .LT. -9999) THEN
          IVAL = -9999
          WRITE (LUNIT,6) ACARD
        ENDIF

        WRITE (ACARD(31:35),3) IVAL
    3   FORMAT (I5.5)
      ENDIF

      RETURN
      END
      SUBROUTINE UPVV (ACARD)

*** UPDATE VARIANCE FACTOR RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXVF = 40 )
      CHARACTER*80 ACARD
      LOGICAL GETIVF
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      CHARACTER*1 TC1, TC2
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /VFTB2/  VFS(MXVF), VTV(MXVF), VFRN(MXVF), VFSS(MXVF)
      COMMON/UNITS/ LUNIT

**** 10-23-03
      if(acard(3:4) .eq. 'HU') return
      if(acard(3:4) .eq. 'hu') return
**********

      READ (ACARD,1) ICODE, IYR1, IMO1, IDY1, IHR1, IMN1, TC1,
     &                      IYR2, IMO2, IDY2, IHR2, IMN2, TC2
    1 FORMAT (2X, I2, 2(I4, 4I2, A1) )
      IF (ICODE .EQ. 22) ICODE = 20
      IF (ICODE .EQ. 32) ICODE = 30
      IF (ICODE .EQ. 42) ICODE = 40
      IF (ICODE .EQ. 54) ICODE = 52
      IF (ICODE .GE. 26  .AND.  ICODE .LE. 29) ICODE = 25

*** TEST FOR BYPASS PROCESSING

      IF (ICODE .EQ. 25) THEN
        IF (LGPS) RETURN
      ELSEIF (ICODE .EQ. 20) THEN
        IF (LDIR) RETURN
      ELSEIF (ICODE .EQ. 30) THEN
        IF (LANG) RETURN
      ELSEIF (ICODE .EQ. 40) THEN
        IF (LZEN) RETURN
      ELSEIF (ICODE .EQ. 52) THEN
        IF (LDIS) RETURN
      ELSEIF (ICODE .EQ. 60) THEN
        IF (LAZI) RETURN
      ELSEIF (ICODE .EQ. 99) THEN
        IF (LDOP) RETURN
      ELSE
        CALL LINE (3)
        WRITE (LUNIT,2)
    2   FORMAT ('0ILLEGAL OBSERVATIONAL TYPE'/)
        RETURN
      ENDIF

*** SET DEFAULTS

      IF (ACARD( 5: 8) .EQ. '    ') IYR1 = 1801
      IF (ACARD( 9:10) .EQ. '  ')   IMO1 = 1
      IF (ACARD(11:12) .EQ. '  ')   IDY1 = 1
      IF (ACARD(17:17) .EQ. ' ')    TC1 = 'Z'
      IF (ACARD(18:21) .EQ. '    ') IYR2 = 2099
      IF (ACARD(22:23) .EQ. '  ')   IMO2 = 12
      IF (ACARD(24:25) .EQ. '  ')   IDY2 = 31
      IF (ACARD(26:27) .EQ. '  ')   IHR2 = 23
      IF (ACARD(28:29) .EQ. '  ')   IMN2 = 59
      IF (ACARD(30:30) .EQ. ' ')    TC2 = 'Z'

      CALL TOMNT (IYR1, IMO1, IDY1, IHR1, IMN1, TC1, IOLD)
      CALL TOMNT (IYR2, IMO2, IDY2, IHR2, IMN2, TC2, INEW)
      ITIME = (IOLD + INEW)/2

      IF ( .NOT. GETIVF(ICODE, ITIME, IVF) ) THEN
        CALL LINE (1)
        WRITE (LUNIT,4) ACARD
    4   FORMAT (1X, A80, ' *** RECORD NOT IN VARIANCE FACTOR TABLE')
      ELSE
        VAL = VFS(IVF)
        IVAL = IDNINT(VAL/0.01D0)

*** CHECK FOR SIZE OF NUMBERS

        IF (IVAL .GT. 99999) THEN
          WRITE (LUNIT,6) IVAL, ACARD
    6     FORMAT (' ***** WARNING, VARIANCE FACTOR VALUE (', I10, ')',/,
     &            6X, ' IS NOT CORRECT IN NEW AFILE --',
     &            ' VALUE TOO LARGE FOR RECORD', /, ' ',
     &            A80)
          IVAL = 99999
        ELSEIF (IVAL .LT. -9999) THEN
          IVAL = -9999
          WRITE (LUNIT,6) ACARD
        ENDIF

        WRITE (ACARD(31:35),3) IVAL
    3   FORMAT (I5)

      ENDIF

      RETURN
      END
      
      LOGICAL FUNCTION VALCH1(RECORD,I1,I2,I3)
      
*** ROUTINE TO CHECK RECORD(I1:I2) WITH CHARACTER STRING REFERENED BY I3
*** VALCH1 = TRUE IF FOUND IN STRING(S)
*** I3 = 1 - 0123456789
*** I3 = 9 - 18,19,20

      COMMON/UNITS/ LUNIT
      
      CHARACTER*80 RECORD
      CHARACTER*10 DIGITS
      INTEGER I1,I2,I3
      
      VALCH1 = .TRUE.
      DIGITS = '0123456789'
      
      IF(I3.EQ.1) THEN
         DO 20 I = I1,I2
	   IF(INDEX(DIGITS,RECORD(I:I)).EQ.0) THEN
	      VALCH1 = .FALSE.
	      RETURN
	   ENDIF
 20      CONTINUE
      ELSEIF(I3.EQ.9) THEN
         IF((RECORD(I1:I2).NE.'18').AND.
     *     (RECORD(I1:I2).NE.'19').AND.
     *     (RECORD(I1:I2).NE.'20')) THEN
             VALCH1 = .FALSE.
	     RETURN
	 ENDIF 
      ENDIF
      END
      
      SUBROUTINE VERDMS (VAL,ID,IM,S,ISIGN)

*** CONVERT ZENITH DISTANCE RADIANS TO DEG, MIN, SEC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /CONST/ PI, PI2, RAD

    1 IF (VAL.GT.PI) THEN
        VAL = VAL-PI-PI
        GO TO 1
      ENDIF

    2 IF (VAL.LT.-PI) THEN
        VAL = VAL+PI+PI
        GO TO 2
      ENDIF

      IF (VAL.LT.0.D0) THEN
        ISIGN = -1
      ELSE
        ISIGN = +1
      ENDIF

      S = DABS(VAL*RAD)
      ID = IDINT(S)
      S = (S-ID)*60.D0
      IM = IDINT(S)
      S = (S-IM)*60.D0

      RETURN
      END
      SUBROUTINE VERTAN (BCARD, IUO, IOBS, B, NX, FATAL, LSN)

*** OBSERVATION EQUATIONS FOR ZENITH DISTANCES

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      PARAMETER ( LENC = 10 )
      CHARACTER*80 BCARD
      CHARACTER*1 TC
      CHARACTER*3 ASS
      LOGICAL FATAL, LSN
      LOGICAL GETSSN, GETPRM, GETIVF
      LOGICAL ADDCON
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      LOGICAL REJECT,GET82
      DIMENSION B(*), NX(*)
      DIMENSION IC(LENC), C(LENC)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT
      SAVE IYR, IMO, IDY

      READ (BCARD,1) IRT, ISSN, JSSN, ID, IM, ASS
    1 FORMAT (7X, I2, 1X, I4, 36X, I4, 9X, I3, I2, A3)
      CALL NBLANK (ASS, 3, IBLK)
      READ (ASS,3) SS
    3 FORMAT (F3.1)

      IF (IRT .EQ. 40) THEN
        READ (BCARD,5) IYR, IMO, IDY, IHR, IMN, TC
    5   FORMAT (39X, 5I2, A1)
        IF (BCARD(40:41) .EQ. '  ') IYR = 84
        IF (BCARD(42:43) .EQ. '  ') IMO = 1
        IF (BCARD(44:45) .EQ. '  ') IDY = 1
        IF (BCARD(46:47) .EQ. '  ') IHR = 0
        IF (BCARD(48:49) .EQ. '  ') IMN = 0
        IF (BCARD(50:50) .EQ. ' ') TC = 'Z'
	CALL GETYR(BCARD,IYR)
*       IYR = IYR + 1900
      ELSE
        READ (BCARD,6) IHR, IMN, TC
    6   FORMAT (45X, 2I2, A1)
        IF (BCARD(46:47) .EQ. '  ') IHR = 0
        IF (BCARD(48:49) .EQ. '  ') IMN = 0
        IF (BCARD(50:50) .EQ. ' ')  TC = 'Z'
      ENDIF
      IF ( BCARD( 6: 6) .EQ. 'R'  .OR.  BCARD( 6: 6) .EQ. 'O'  .OR.
     &     BCARD( 6: 6) .EQ. 'F' ) THEN
        NREJ = NREJ + 1
        REJECT = .TRUE.
      ELSE
        REJECT = .FALSE.
      ENDIF

      IF (GET82(ISSN,I)) RETURN
      IF (GET82(JSSN,I)) RETURN

        IF ( .NOT. GETSSN(ISSN, ISN) ) THEN
          CALL LINE (3)
          WRITE (LUNIT,2) BCARD
    2     FORMAT ('0NO *80* RECORD FOR--', A80, /)
        ELSEIF ( .NOT. GETSSN(JSSN, JSN) ) THEN
          CALL LINE (3)
          WRITE (LUNIT,2) BCARD
        ELSEIF ( IDIM .EQ. 2  .AND.  L2HLF  .AND.

*** IF THIS IS A 2 DIM ADJUSTMENT AND ONE OF THE STATIONS IS NOT A
*** DUAL HEIGHT STATION THEN RETURN

     &           ( I2HLF(ISN) .EQ. 0  .OR. I2HLF(JSN) .EQ. 0) ) THEN
          RETURN
        ELSE

*** RETRIEVE THE STD DEV

          CALL STDDEV (BCARD, IRT, SD, REJECT)

*** KIND = 9 ZENITH DISTANCE

          KIND = 9
          IOBS = IOBS + 1
          NOBS = NOBS + 1
          NZD = NZD + 1
          LSN = .TRUE.
          OBSB = (ID + IM/60.D0 + SS/3600.D0)/RAD
          CALL TOMNT (IYR, IMO, IDY, IHR, IMN, TC, ITIME)
          IF ( .NOT. GETPRM(40, ITIME, IAUX) ) IAUX = 0
          IF ( .NOT. GETIVF(40, ITIME, IVF) ) IVF = 0
          IGRT = 0
          CALL FORMIC (KIND, ISN, JSN, IAUX, IC, LENG, IGRT)
          CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
          CALL COMPOB (KIND, OBS0, B, OBSB, ISN, JSN, IAUX, IGRT)
          CMO = OBS0 - OBSB
          IF (IMODE .EQ. 0) CMO = 0.D0
          CALL BIGV (KIND, ISN, JSN, IOBS, IVF, CMO, SD, FATAL)

*** UPDATE CONNECTIVITY

          IF ( .NOT. ADDCON(IC, LENG, NX) ) THEN
            WRITE (LUNIT,667)
  667       FORMAT ('0INSUFFICIENT STORAGE FOR ZENITH DIST.', /)
          write (lunit,*) 'subs2, 100'
            CALL ABORT2
          ENDIF

          WRITE (IUO) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD, IOBS,
     &                IVF, IAUX, IGRT

        ENDIF

      RETURN
      END
      LOGICAL FUNCTION VFCVRG (IUO,IUO2,B,G,A,NX,GOOGE)

*** REFORM OBS. EQS., COMPUTE INVERSE AND V.F.'S, AND TEST CONVERGE

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10, MXVF = 40 )
      LOGICAL GETA
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL FATAL,FATAL2,FIXVF,PROP,INVERT
      DIMENSION B(*),G(*),A(*),NX(*),GOOGE(*)
      DIMENSION IC(LENC),C(LENC)
      DIMENSION COVECF(3,3), C3(3,LENC)
      DIMENSION COVLBI(3,3), DUMMAT(3,3), COVLA(3,3), Q(3,3)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /CONST/  PI, PI2, RAD
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON /VFTB2/  VFS(MXVF), VTV(MXVF), VFRN(MXVF), VFSS(MXVF)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      VFCVRG = .TRUE.
      VFSTOL = 1.D-5
      VFCTOL = 1.D-2

***  COMPLETE THE COMPUTATION OF THE GOOGE NUMBERS

      DO 44 I = 1,NUNK
        IF ( .NOT. GETA(I,I,VAL,A,NX)) THEN
          CALL INVIUN (I,ISN,ITYP,ICODE)
          WRITE (LUNIT,55) I,ISN,ITYP,ICODE
   55     FORMAT ('0 FATAL ERROR IN GOOGE COMPUTATION IN VFCVRG',4I8)
          write (lunit,*) 'subs2, 101'
          CALL ABORT2
        ENDIF

***  ALLOW FOR GOOGE OF ZERO IF SOLVING A SINGULAR SYSTEM

        VAL2 = VAL*VAL
        GOOGE(I) = DIVID( VAL2, GOOGE(I) )
   44 CONTINUE

*** INVERT WITHIN PROFILE

      IF (IMODE .NE. 1) THEN
        IF ( .NOT. INVERT(A,NX)) THEN
          WRITE (LUNIT,666)
  666     FORMAT ('0STATE ERROR IN INVERSION OF VFCVRG')
          write (lunit,*) 'subs2, 102'
          CALL ABORT2
        ENDIF
      ENDIF

*** INITIALIZE VARIANCE FACTORS

      DO 1 IVF = 1, NVFTOT
        IF ( .NOT. FIXVF(IVF)) THEN
          VTV(IVF) = 0.D0
          VFRN(IVF) = 0.D0
        ENDIF
    1 CONTINUE

*** LOOP OVER THE OBSERVATIONS

      FATAL2 = .FALSE.
      FATAL = .FALSE.
      REWIND IUO
      REWIND IUO2

  100 READ (IUO,END=777) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                   IVF,IAUX,IGRT
      IF (KIND .LE. 999) THEN
        IF (KIND .LT. 18  .OR.  KIND .GT. 20) THEN
          CALL FORMC (KIND,C,B,ISN,JSN,IAUX,IGRT)
          CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
          IF (KIND .EQ. 1  .OR.  KIND .EQ. 2) THEN
            CMO = OBS0
          ELSE
            CMO = OBS0 - OBSB
            IF ( KIND .EQ.  8  .OR.  KIND .EQ. 10  .OR.
     &           KIND .EQ. 11  .OR.  KIND .EQ. 12) THEN
              IF (CMO .GT. PI) THEN
                CMO = CMO - PI - PI
              ELSEIF (CMO .LT. -PI) THEN
                CMO = CMO + PI + PI
              ENDIF
            ENDIF
          ENDIF

          IF (IVF .NE. 0) THEN
            IF ( .NOT. FIXVF(IVF)) THEN
              IF (PROP(C,IC,LENG,VARL0,A,NX,IFLAG)) THEN
                VARLB = SD*SD*VFS(IVF)
                RN = 1.D0-VARL0/VARLB
                IF (RN .LT. 1.0D-6) RN = 0.D0
                VFRN(IVF) = VFRN(IVF)+RN
                VTV(IVF) = VTV(IVF)+CMO*CMO/VARLB
              ELSEIF (IFLAG .EQ. 1) THEN
                WRITE (LUNIT,668)
  668           FORMAT ('0STATE ERROR IN VFCVRG')
          write (lunit,*) 'subs2, 103'
                CALL ABORT2
              ELSE
                WRITE (LUNIT,669)
  669           FORMAT ('0PROFILE ERROR IN VFCVRG')
          write (lunit,*) 'subs2, 104'
                CALL ABORT2
              ENDIF
            ENDIF
          ENDIF

          CALL BIGL (KIND,ISN,JSN,IOBS,IVF,CMO,SD,FATAL)
          WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                 IVF,IAUX,IGRT

        ELSEIF (KIND .EQ. 18) THEN

*** CORRELATED TYPE (DOPPLER)

          OBSBX = OBSB
          SDX = SD
          IOBSX = IOBS
          CALL FORMC (KIND,C,B,ISN,JSN,IAUX,IGRT)
          DO 101 I = 1, LENG
            C3(1,I) = C(I)
  101     CONTINUE
          CALL COMPOB (KIND,OBS0X,B,OBSBX,ISN,JSN,IAUX,IGRT)
          CMOX = OBS0X - OBSBX
          CALL BIGL (KIND,ISN,JSN,IOBSX,IVF,CMOX,SDX,FATAL)
          WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMOX,OBSBX,SDX,IOBSX,
     &                 IVF,IAUX,IGRT
          READ (IUO,END=777) KIND,ISN,JSN,IC,C,LENG,CMO,OBSBY,SDY,IOBSY,
     &                       IVF,IAUX,IGRT
          CALL FORMC (KIND,C,B,ISN,JSN,IAUX,IGRT)
          DO 102 I = 1, LENG
            C3(2,I) = C(I)
  102     CONTINUE
          CALL COMPOB (KIND,OBS0Y,B,OBSBY,ISN,JSN,IAUX,IGRT)
          CMOY = OBS0Y - OBSBY
          CALL BIGL (KIND,ISN,JSN,IOBSY,IVF,CMOY,SDY,FATAL)
          WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMOY,OBSBY,SDY,IOBSY,
     &                 IVF,IAUX,IGRT
          READ (IUO,END=777) KIND,ISN,JSN,IC,C,LENG,CMO,OBSBZ,SDZ,IOBSZ,
     &                       IVF,IAUX,IGRT
          CALL FORMC (KIND,C,B,ISN,JSN,IAUX,IGRT)
          DO 103 I = 1, LENG
            C3(3,I) = C(I)
  103     CONTINUE
          CALL COMPOB (KIND,OBS0Z,B,OBSBZ,ISN,JSN,IAUX,IGRT)
          CMOZ = OBS0Z - OBSBZ
          CALL BIGL (KIND,ISN,JSN,IOBS,IVF,CMOZ,SDZ,FATAL)
          WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMOZ,OBSBZ,SDZ,IOBSZ,
     &                 IVF,IAUX,IGRT
          READ (IUO,END=777) ISIGNL, ( (COVECF(I,J), J = 1,3), I = 1,3 )
          WRITE (IUO2) ISIGNL, ( (COVECF(I,J), J = 1,3), I = 1,3 )

*** CALCULATION OF THE VARIANCE FACTOR

          IF (IVF .NE. 0) THEN
            IF ( .NOT. FIXVF(IVF)) THEN
              CALL COMSLA ( A, NX, IC, C3, LENG, COVLA)

*** COMPUTATION OF REDUNDENCY NUMBER

              CALL EQ ( COVECF, COVLBI, 3, 3)
              CALL INVRT3 (COVLBI)
              CALL EQ (COVLBI, DUMMAT, 3, 3)
              CALL GETQ (DUMMAT, COVLA, Q)

              RN = Q(1,1) + Q(2,2) + Q(3,3)
              IF (RN .LT. 1.0D-6) RN = 0.D0
              VFRN(IVF) = VFRN(IVF) + RN
              VTV(IVF) = VTV(IVF) + CMOX*CMOX/( SDX*SDX*VFS(IVF) )
     &                            + CMOY*CMOY/( SDY*SDY*VFS(IVF) )
     &                            + CMOZ*CMOZ/( SDZ*SDZ*VFS(IVF) )
            ENDIF
          ENDIF
        ENDIF

*** GPS OBSERVATIONS

      ELSE
        NVEC = ISN
        IAUX = JSN
        NR = 3*NVEC
        WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &               IVF,IAUX,IGRT
        CALL FORMG2 (IUO,IUO2,NVEC,NR,G,B,A,NX,FATAL,IVF,IAUX,IGRT)
      ENDIF
      GO TO 100

*** ABORT DUE TO LARGE MISCLOSURES

  777 IF (FATAL) THEN
        CALL LINE (3)
        WRITE (LUNIT,3) VM
    3   FORMAT ('0TERMINATED DUE TO MISCLOSURES (C-O)/SD EXCEEDING ' ,
     &          F11.1/)
          write (lunit,*) 'subs2, 106'
        CALL ABORT2
      ENDIF

      REWIND IUO
      REWIND IUO2

*** EXCHANGE PRIMARY/SECONDARY OBS EQ FILE INDICATOR

      ITEMP = IUO
      IUO = IUO2
      IUO2 = ITEMP

*** HEADING

      CALL LINE (3)
      WRITE (LUNIT,20)
   20 FORMAT ('0NUM    VAR. FACTOR    DEG. OF FREE.',
     &         T44,'V.F. RATIO(COMP/INIT)'/)

*** COMPUTE VARIANCE FACTOR AND TEST CONVERGENCE

      DO 2 IVF = 1,NVFTOT
        IF (FIXVF(IVF)) THEN
          CALL LINE (1)
          WRITE (LUNIT,10) IVF,VFS(IVF)
   10     FORMAT (1X,I5,F14.3,'   *** FIXED ***')
        ELSEIF (VFRN(IVF) .LT. VFSTOL) THEN
          CALL LINE (1)
          WRITE (LUNIT,11) IVF,VFRN(IVF)
   11     FORMAT (1X,I5,'  FDF = ',1PD8.2,' **** IS SINGULAR ********')
          FATAL2 = .TRUE.
        ELSE
          VF = VTV(IVF)/VFRN(IVF)
          IF (DABS(VF-1.D0) .GT. VFCTOL) VFCVRG = .FALSE.
          VFS(IVF) = VF*VFS(IVF)
          VFRAT = VFS(IVF)/VFSS(IVF)
          CALL LINE (1)
          WRITE (LUNIT,12) IVF,VFS(IVF),VFRN(IVF),VFRAT
   12     FORMAT (1X,I5,2F14.3,11X,F15.3)
        ENDIF
   2   CONTINUE

*** ABORT IF SINGULAR VARIANCE FACTORS

      IF (FATAL2) THEN
        CALL LINE (3)
        WRITE (LUNIT,4)
    4   FORMAT ('0TERMINATED DUE TO VARIANCE FACTOR SINGULARITIES  ')
          write (lunit,*) 'subs2, 107'
        CALL ABORT2
      ENDIF

      RETURN
      END
      SUBROUTINE VFGPS (G, NR, NC, NVEC, LENG, B, ICM, NICM, KINDS,ISNS,
     &            JSNS, LOBS, IAUX, IVF, FATAL, A, NX, N2, N3, N4, IGRT)

*** COMPUTE AND DECORRELATE OBS. EQ. AND MISCLOSURE

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( NVECS = 700, LENC = 10, MXVF = 40, MXSSN = 9999 )
      DIMENSION ICM((NVECS+1)*3+1+3)
      DIMENSION KINDS(NVECS*3), ISNS(NVECS*3), JSNS(NVECS*3)
      DIMENSION LOBS(NVECS*3)
      LOGICAL FATAL, PROP, FIXVF
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      DIMENSION CM((NVECS+1)*3+1+3)
      DIMENSION B(*), G(NR,NC), A(*), NX(*)
      DIMENSION IC(LENC), C(LENC)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /LASTSN/ ISNX, JSNX
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /VFTB2/  VFS(MXVF), VTV(MXVF), VFRN(MXVF), VFSS(MXVF)
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

*** LENGTH OF COEFF. ARRAYS

      IF ( .NOT. L2HLF) THEN
        LENGC = 2*IDIM
      ELSE
        LENGC = 2*3
      ENDIF
      IF (IAUX .GT. 0) LENGC = LENGC + 1
      IF (IGRT .GT. 0) LENGC = LENGC + 3

*** SET LAST STATION ENCOUNTERED TO NULL

      ISNX = 0
      JSNX = 0

*** DETERMINE WORK ARRAY COLUMN POINTERS

      CALL GLOCAT (N1, N2, N3, N4, N5, NVEC, LENG)

*** (STEP 1.)  COMPUTE COEFFICIENTS

      DO 10 I = 1, NR
        KIND = KINDS(I)
        ISN = ISNS(I)
        JSN = JSNS(I)
        CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
        CALL GETGIC (ISN, JSN, IAUX, ICM, IC, N2, NICM, IGRT)

        DO 2 J = 1, LENG
          K = J + N2
          G(I,K) = 0.D0
    2   CONTINUE

        DO 20 J = 1, LENGC
          G(I,IC(J) ) = C(J)
   20   CONTINUE
   10 CONTINUE

*** (STEP 2.)  PROPAGATE FOR STD. DEV. OF RESIDUALS


*** (STEP 3.)  DECORRELATE COEFFICIENTS

      DO 30 I = 1, LENG
        CALL GETRHS (G, NR, NC, N2, N2+I)
        CALL COMRHS (G, NR)
        CALL PUTRHS (G, NR, NC, N2, N2+I)
   30 CONTINUE

*** (STEP 4.)  PROPAGATE FOR REDUNDANCY NUMBER (IF VARIANCE FACTOR)

      IF (IVF .GT. 0) THEN
        IF ( .NOT. FIXVF(IVF) ) THEN
          DO 80 IROW = 1, NR

*** BYPASS REJECTION WHEN COMPUTING A REDUNDANCY NUMBER

            IF (G(IROW,1) .GE. 100.D0) GO TO 80
            DO 70 I = 1, NICM
              CM(I) = G(IROW,N2+I)
   70       CONTINUE

            IF ( .NOT. PROP(CM, ICM, NICM, W, A, NX, IFLAG) ) THEN
              IF (IFLAG .EQ. 1) THEN
                WRITE (LUNIT,668)
  668           FORMAT ('0STATE ERROR IN PROP OF VFGPS')
          write (lunit,*) 'subs2, 108'
                CALL ABORT2
              ELSE
                WRITE (LUNIT,669)
  669           FORMAT ('0PROFILE ERROR IN PROP OF VFGPS')
          write (lunit,*) 'subs2, 109'
                CALL ABORT2
              ENDIF
            ENDIF

            W = W/VFS(IVF)
            RN = 1.D0 - W
            IF (RN .LT. 1.D-6) RN = 0.D0
            VFRN(IVF) = VFRN(IVF) + RN
   80     CONTINUE
        ENDIF
      ENDIF

*** (STEP 5.)  COMPUTE OBSERVATIONS

      DO 90 I = 1, NR
        KIND = KINDS(I)
        ISN = ISNS(I)
        JSN = JSNS(I)
        CALL COMPOB (KIND, OBS0, B, ADUMMY, ISN, JSN, IAUX, IGRT)
        CMO = OBS0 - G(I,NC)
        G(I,N2) = CMO
   90 CONTINUE

*** (STEP 6.)  DECORRELATE MISCLOSURE

      CALL COMRHS (G, NR)
      CALL PUTRHS (G, NR, NC, N2, N4)
      SD = 1.D0
      DO 16 I = 1, NR
        CMO = G(I,N2)
        IF (IVF .GT. 0) THEN
          IF ( .NOT. FIXVF(IVF) ) VTV(IVF) = VTV(IVF) + CMO*CMO/VFS(IVF)
        ENDIF
        CALL BIGV (KINDS(I), ISNS(I), JSNS(I), LOBS(I), IVF, CMO, SD,
     &             FATAL)
   16 CONTINUE

      RETURN
      END
      SUBROUTINE VSHIFT (IP,IOBS,V, IOBSNO)

*** SHIFT ARRAY CONTENTS AND LOAD ARRAYS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRD = 20 )
      COMMON /RESTAT/ VSD20(MXRD), VMAX, VMIN, VSDMAX, VSDMIN, VSUM,
     &                VSDSUM, VSD21, VSD22, VSD23, VSD24, VSD25, VSD26,
     &                VSD27, VSD28, VSD29, VSD210, VSD211, VSD212,
     &                VABS1, VABS2, VABS3, VABS4, VABS5, VABS6, VABS7,
     &                VABS8, VABS9, VABS10, VABS11, VABS12, I20(MXRD),
     &                N0, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10,
     &                N11, N12, NI20
      COMMON /REDUN/  RN1, RN2, RN3, RN4, RN5, RN6, RN7, RN8, RN9,
     &                RN10, RN11, RN12

*** PROTECT AGAINST ILLEGAL VALUES

      IF (IP .LE. 0  .OR.  IP .GT. MXRD) RETURN

*** SHIFT ARRAY CONTENTS

      IF (IP .LT. MXRD) THEN
        DO 1 I = MXRD-1,IP,-1
          VSD20(I+1) = VSD20(I)
          I20(I+1) = I20(I)
    1   CONTINUE
      ENDIF

*** LOAD THE SHIFTED ARRAY

      VSD20(IP) = V
      I20(IP) = IOBSNO
      IF (NI20 .LT. MXRD) NI20 = NI20+1

      RETURN
      END
      
**v 4.28

      SUBROUTINE VSHIFTV (IP,IOBS,V,IOBSNO)

*** SHIFT ARRAY CONTENTS AND LOAD ARRAYS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRDV = 20 )
      COMMON /RESTA3/ V20(MXRDV), I20V(MXRDV), NI20V  

*** PROTECT AGAINST ILLEGAL VALUES

      IF (IP .LE. 0  .OR.  IP .GT. MXRDV) RETURN

*** SHIFT ARRAY CONTENTS

      IF (IP .LT. MXRDV) THEN
        DO 1 I = MXRDV-1,IP,-1
          V20(I+1) = V20(I)
          I20V(I+1) = I20V(I)
    1   CONTINUE
      ENDIF

*** LOAD THE SHIFTED ARRAY

      V20(IP) = V
      I20V(IP) = IOBSNO
      IF (NI20V .LT. MXRDV) NI20V = NI20V+1

      RETURN
      END
      
      
************
** v 4.29f

      SUBROUTINE VSHIFTU (IP,IOBS,V,IOBSNO)

*** SHIFT ARRAY CONTENTS AND LOAD ARRAYS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRDU = 20 )
      COMMON /RESTA7/ VDU20(MXRDU), I20DU(MXRDU), NI20DU  

*** PROTECT AGAINST ILLEGAL VALUES

      IF (IP .LE. 0  .OR.  IP .GT. MXRDU) RETURN

*** SHIFT ARRAY CONTENTS

      IF (IP .LT. MXRDU) THEN
        DO 1 I = MXRDU-1,IP,-1
          VDU20(I+1) = VDU20(I)
          I20DU(I+1) = I20DU(I)
    1   CONTINUE
      ENDIF

*** LOAD THE SHIFTED ARRAY

      VDU20(IP) = V
      I20DU(IP) = IOBSNO
      IF (NI20DU .LT. MXRDU) NI20DU = NI20DU+1

      RETURN
      END

***********************

** v 4.29i

      SUBROUTINE VSHIFTL (IP,IOBS,V,IOBSNO)

*** SHIFT ARRAY CONTENTS AND LOAD ARRAYS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRDL = 20 )
      COMMON /RESTA8/ VDL20(MXRDL), I20DL(MXRDL), NI20DL  

*** PROTECT AGAINST ILLEGAL VALUES

      IF (IP .LE. 0  .OR.  IP .GT. MXRDL) RETURN

*** SHIFT ARRAY CONTENTS

      IF (IP .LT. MXRDL) THEN
        DO 1 I = MXRDL-1,IP,-1
          VDL20(I+1) = VDL20(I)
          I20DL(I+1) = I20DL(I)
    1   CONTINUE
      ENDIF

*** LOAD THE SHIFTED ARRAY

      VDL20(IP) = V
      I20DL(IP) = IOBSNO
      IF (NI20DL .LT. MXRDL) NI20DL = NI20DL+1

      RETURN
      END

***********************
      SUBROUTINE WR3X3 (A, NAME, IUNIT)

***   THIS SUBROUTINE PRINTS A 3 BY 3 MATRIX.  THE PARAMETER NAME IS THE
***   TITLE THAT WILL BE PRINTED ABOVE THE MATRIX.  IT MUST CONSIST OF
***   AN ALPHANUMERIC STRING.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION A(3,3)
      CHARACTER*(*) NAME
      COMMON/UNITS/ LUNIT

      WRITE (IUNIT,10) NAME
 10   FORMAT (//, ' ',  5X, A)

      IF (IUNIT .GT. 0) THEN
        DO 20 I = 1, 3
          WRITE (IUNIT,32) (A(I,J), J=1, 3)
   32     FORMAT ( 3(1X, G16.9) )
 20     CONTINUE
      ELSE
        DO 30 I = 1, 3
          WRITE (LUNIT,32) (A(I,J), J=1, 3)
 30     CONTINUE
      ENDIF

      RETURN
      END

