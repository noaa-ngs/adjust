C %P%

*** file newsubs.for

      SUBROUTINE AFNL (ACARD)
********************************************************************************
* SET PROCESSING MODE FOR NETWORK AND LOCAL ACCURACIES
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*99 SCCSID
      CHARACTER*80 ACARD
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7 ADJFIL, NAFILE
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      LOGICAL NLMODE,NLSCL
      COMMON /NLOPT/ NLMODE,NLSCL,NLUN1,NLUN2,NLUN3,NLUN4,NLUN5

C      SCCSID='$Id: newsubs.for 81724 2014-12-22 16:17:53Z jarir.saleh $	20$Date: 2008/08/25 16:02:42 $ NGS'
      NLMODE = .TRUE.
      IF (ACARD(3:3) .EQ. 'Y') THEN
        NLSCL=.TRUE.
      ELSEIF (ACARD(3:3) .EQ. 'N' .OR. ACARD(3:3) .EQ. ' ') THEN
        NLSCL = .FALSE.
      ENDIF
      RETURN
      END

      SUBROUTINE CIRC95(SIGMX,SIGMN,H95)
********************************************************************************
* COMPUTE THE RADIUS OF A CIRCLE THAT CONTAINS 95% PROBABILITY
* IN A TWO DIMENSIONAL NORMAL DISTRIBUTION
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )

*** 95% values (cf. leenhouts, navigation, 32(1), spring 1985, pp.16-28)
      parameter(q0=1.960790d0)
      parameter(q1=0.004071d0)
      parameter(q2=0.114276d0)
      parameter(q3=0.371625d0)
      parameter(q4=1.960d0)

      c1=sigmn/sigmx
      c2=c1*c1
      c3=c2*c1
      cayp=q0+q1*c1+q2*c2+q3*c3
      h95=cayp*sigmx
      RETURN
      END
C************************************************************************************
      SUBROUTINE HEAD10
********************************************************************************
* PRINT A HEADING AND SUBHEADING #10
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE
      LOGICAL NLMODE,NLSCL
      COMMON /NLOPT/ NLMODE,NLSCL,NLUN1,NLUN2,NLUN3,NLUN4,NLUN5
      COMMON/UNITS/ LUNIT

      CALL HEAD
      IF(NLSCL) THEN
        WRITE(LUNIT,6000)
 6000   FORMAT(//T1,'NETWORK & LOCAL ACCURACIES (95% CONFIDENCE LEVEL)'/
     *,'COMPUTED USING A POSTERIORI STANDARD DEVIATION OF UNIT WEIGHT')
        WRITE(LUNIT,*) '   '
        WRITE(LUNIT,*) '   '
      ELSE
        WRITE(LUNIT,6001)
 6001   FORMAT(//T1,'NETWORK & LOCAL ACCURACIES (95% CONFIDENCE LEVEL)'/
     *,'COMPUTED USING A PRIORI STANDARD DEVIATION OF UNIT WEIGHT')
        WRITE(LUNIT,*) '   '
        WRITE(LUNIT,*) '   '
      ENDIF
c     WRITE(LUNIT,6005)
 6005 FORMAT(/T10, 'STATION',T52,'TO STATION',T103,'ACCURACIES (CM.)')
c     WRITE(LUNIT,6006)
 6006 FORMAT('  ISN  SSN   NAME', T42, 'TYPE     ISN  SSN  NAME',
     *   T100,'    HORIZ  ELLIP')
      ILINE = 10
      RETURN
      END
C************************************************************************************
      LOGICAL FUNCTION ISNADD(ISN,ISNS,LISN,NSTA)
********************************************************************************
* ADD STATION NUMBER ISN TO THE LIST ISNS OF STATIONS IN A SESSION
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION ISNS(*)
      IF (LISN.EQ.0) THEN
        ISNS(1)=ISN
        LISN=1
        ISNADD=.TRUE.

      ELSEIF(LISN.LE.NSTA) THEN
        DO 10 I=1,LISN
          IF (ISNS(I).EQ.ISN) THEN
            ISNADD=.TRUE.
            RETURN
          ENDIF
   10   CONTINUE

*** REACHED THE END OF LIST WITHOUT FINDING ISN. 
*** TRY TO ADD IT AT END
        IF (LISN.LT.NSTA) THEN
          LISN=LISN+1
          ISNS(LISN)=ISN
          ISNADD=.TRUE.
          RETURN
        ELSE
          ISNADD=.FALSE.
          RETURN
        ENDIF
      ENDIF
      END

      SUBROUTINE NL21 (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
********************************************************************************
* READ A BLUE BOOK RECORD AND GENERATE CONNECTIONS FOR
* BLUEBOOK RECORD WITH OBSERVED SSN IN CC51-54
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
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

      READ (BCARD,1) ISSN,JSSN
 1    FORMAT (T11,I4,T51,I4)

      IF ( BCARD( 6: 6) .EQ. 'R'  .OR.  BCARD( 6: 6) .EQ. 'O'  .OR.
     &     BCARD( 6: 6) .EQ. 'F' ) THEN
***11-30-00
*       NREJ = NREJ+1
        REJECT = .TRUE.
      ELSE
        REJECT = .FALSE.
      ENDIF

      IF (GET82(ISSN,I)) RETURN
      IF (GET82(JSSN,I)) RETURN


      IF ( .NOT. GETSSN(ISSN,ISN)) THEN
        CALL LINE (3)
        WRITE (LUNIT,2) BCARD
 2      FORMAT ('0NO *80* RECORD FOR--',A80/)
      ELSEIF ( .NOT. GETSSN(JSSN,JSN)) THEN
        CALL LINE (3)
        WRITE (LUNIT,2) BCARD
      ELSE

*** ADD THE CONNECTION
        IC(1)=ISN
        IC(2)=JSN
        LENG=2
        IF ( .NOT. ADDCON(IC,LENG,NX)) THEN
          WRITE (LUNIT,667)
 667      FORMAT ('0INSUFFICIENT STORAGE FOR HOR. DIRECTIONS'/)
          CALL ABORT2
        ENDIF
      ENDIF
      RETURN
      END

      SUBROUTINE NL22 (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
********************************************************************************
* READ A BLUE BOOK RECORD AND GENERATE CONNECTIONS FOR
* BLUEBOOK RECORD WITH OBSERVED SSN IN CC46-49 
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
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

      READ (BCARD,1) ISSN,JSSN
 1    FORMAT (T11,I4,T46,I4)

      IF ( BCARD( 6: 6) .EQ. 'R'  .OR.  BCARD( 6: 6) .EQ. 'O'  .OR.
     &     BCARD( 6: 6) .EQ. 'F' ) THEN
***11-30-00     
*       NREJ = NREJ+1
        REJECT = .TRUE.
      ELSE
        REJECT = .FALSE.
      ENDIF

      IF (GET82(ISSN,I)) RETURN
      IF (GET82(JSSN,I)) RETURN

      IF ( .NOT. GETSSN(ISSN,ISN)) THEN
        CALL LINE (3)
        WRITE (LUNIT,2) BCARD
 2      FORMAT ('0NO *80* RECORD FOR--',A80/)
      ELSEIF ( .NOT. GETSSN(JSSN,JSN)) THEN
        CALL LINE (3)
        WRITE (LUNIT,2) BCARD
      ELSE

*** ADD THE CONNECTION
        IC(1)=ISN
        IC(2)=JSN
        LENG=2
        IF ( .NOT. ADDCON(IC,LENG,NX)) THEN
          WRITE (LUNIT,667)
 667      FORMAT ('0INSUFFICIENT STORAGE FOR HOR. DIRECTIONS'/)
          CALL ABORT2
        ENDIF
      ENDIF
      RETURN
      END

      SUBROUTINE NL3 (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
********************************************************************************
* READ A BLUE BOOK RECORD AND GENERATE CONNECTIONS FOR
*    BLUEBOOK RECORD WITH TARGET SSN'S IN CC51-54 AND 72-75 
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
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

      READ (BCARD,1) ISSN,JSSN,KSSN
 1    FORMAT (T11,I4,T51,I4,T72,I4)

      IF ( BCARD( 6: 6) .EQ. 'R'  .OR.  BCARD( 6: 6) .EQ. 'O'  .OR.
     &     BCARD( 6: 6) .EQ. 'F' ) THEN
***11-30-00
*       NREJ = NREJ+1
        REJECT = .TRUE.
      ELSE
        REJECT = .FALSE.
      ENDIF

      IF (GET82(ISSN,I)) RETURN
      IF (GET82(JSSN,I)) RETURN

      IF ( .NOT. GETSSN(ISSN,ISN)) THEN
        CALL LINE (3)
        WRITE (LUNIT,2) BCARD
 2      FORMAT ('0NO *80* RECORD FOR--',A80/)
      ELSEIF ( .NOT. GETSSN(JSSN,JSN)) THEN
        CALL LINE (3)
        WRITE (LUNIT,2) BCARD
      ELSEIF ( .NOT. GETSSN(KSSN,KSN)) THEN
        CALL LINE (3)
        WRITE (LUNIT,2) BCARD
      ELSE

*** ADD THE CONNECTION
        IC(1)=ISN
        IC(2)=JSN
        IC(3)=KSN
        LENG=3
        IF ( .NOT. ADDCON(IC,LENG,NX)) THEN
          WRITE (LUNIT,667)
 667      FORMAT ('0INSUFFICIENT STORAGE FOR BLUEBOOK CONNNECTIONS'/)
          CALL ABORT2
        ENDIF
      ENDIF
      RETURN
      END
C*******************************************************************************
      SUBROUTINE NLACUR (A, B, NX, SHIFTS, GOOGE, SIGUWT)

* COMPUTE NETWORK AND LOCAL ACCUARCIES
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )

*** 95% values (cf. leenhouts, navigation, 32(1), spring 1985, pp.16-28)
      parameter(q0=1.960790d0)
      parameter(q1=0.004071d0)
      parameter(q2=0.114276d0)
      parameter(q3=0.371625d0)
      parameter(q4=1.960d0)

      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LOCSSN, GETA
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      LOGICAL ELFLAG, DFFLAG
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*30 NAMES, NAME1,NAME2
      CHARACTER*7 ADJFIL, NAFILE
      CHARACTER*1 ADIR1, ADIR2, ADX, ADY, AEHT
      CHARACTER   dummy*10,dummy1*10,dashes*80

      real*8      MedH(MXSSN),MedV(MXSSN),MedD(MXSSN)  !for replacing the mean of the accuracies by their median

** 5-29-03

      character   accode*1,PIDs*6,ST*2,SSN_char*4,SSN_char1*4

      DIMENSION A(*), B(*), NX(*), SHIFTS(*), GOOGE(*)
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /CONST/  PI, PI2, RAD
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /FLAGS/  ELFLAG(MXSSN), DFFLAG(MXSSN)
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      LOGICAL NLMODE,NLSCL
      COMMON /NLOPT/ NLMODE,NLSCL,NLUN1,NLUN2,NLUN3,NLUN4,NLUN5
      COMMON/UNITS/ LUNIT
      COMMON/PIDs / PIDs (MXSSN),ST(MXSSN)

*** HEADING

      CALL HEAD10
      REWIND NLUN1       !contains SSNs of all stations connected to a certain station
      REWIND NLUN2      
      REWIND NLUN3       
      REWIND NLUN4       
      REWIND NLUN5       

C  Was a scratch file. therefore it was not found at this point

      OPEN(NLUN2,file='NLUN2',FORM='FORMATTED')
      OPEN(NLUN3,file='NLUN3',FORM='FORMATTED')
      OPEN(NLUN4,file='NLUN4',FORM='FORMATTED')
      OPEN(NLUN5,file='NLUN5',FORM='FORMATTED')

C  Write a header for the 95% network and local accuracies in ADJUST.OUT 

      write(NLUN4,*) 'NETWORK ACCURACIES (CENTIMETERS)'
      write(NLUN4,*) '    '                                       
      WRITE(NLUN4,6006)
 6006 FORMAT('  ISN  SSN   STATION DESIGNATION',23x,'HORIZ    UP ')
      write(NLUN5,*)
     &  'LOCAL ACCURACIES (CENTIMETERS) AND DISTANCE (KILOMETERS)'
      write(NLUN5,*) '    '                                       
      WRITE(NLUN5,6007)
 6007 FORMAT('          ISN  SSN   STATION DESIGNATION')

C  Start loop over all stations in the network

      DO ISN=1,NSTA
        IF (.NOT. LOCSSN(ISN,ISSN)) THEN
          WRITE(LUNIT,666) ADJFIL, ISN
 666      FORMAT ('0SSN TABLE ERROR IN ',A8,' --',I8)
          CALL ABORT2
        ENDIF

        NAME1 = NAMES(ISN)

C  Extract network accuracy. 

        CALL pellip (isn,sigmx,sigmn,theta,a,nx,b,vardn,varde,coven)
        SD_dn  = sqrt(vardn)
        SD_de  = sqrt(varde)

C  v6.0, added at 05/20/2012

        corren = coven/SD_dn/SD_de
        if (NLSCL) then
          SD_dn = SD_dn*SIGUWT
          SD_de = SD_de*SIGUWT
        endif

C  so far for v6.0 (05/20/2012)

        CALL CIRC95(SIGMX,SIGMN,H95)

C  Grab vertical accuracy of the base node of the baseline

        IZ = IUNSTA(ISN,3)
        IF ( .NOT. GETA(IZ, IZ, CZ, A, NX) ) THEN
          WRITE (LUNIT,667) ADJFIL, ISN, IZ
  667     FORMAT ('0GET SIGMA ERROR IN ', A8, ' --', 2I8)
          CALL ABORT2
        ENDIF
        VARU = CZ
        SD_U = sqrt(VARU)

C v6.0,  added at 05/20/2012

        if (NLSCL) then
          SD_U = SD_U*SIGUWT
        endif

C  so far for v6.0, end addition at 05/20/2012

        v95 =q4*dsqrt(varu)
        IF(NLSCL) THEN
          H95=H95*SIGUWT
          V95=V95*SIGUWT
        ENDIF
        H95CM=H95*100.
        V95CM=V95*100.

c       WRITE(LUNIT,6010) ISN,ISSN,NAME1,H95CM,V95CM
c6010   FORMAT(2I5, 1X,A30,'NETWORK',T100,2F8.2)

        WRITE(NLUN4,6010) ISN,ISSN,NAME1,H95CM,V95CM
 6010   FORMAT(2I5, 1X,A30,10x,2F8.2,1x,'CM')

C  Write a line about the from station in the local accuracy file

c       WRITE(NLUN5,6011) ISN,ISSN,NAME1
c6011   FORMAT('FROM    ',2I5, 1X,A30,5x,'HORIZ    UP       DISTANCE')
*
********WRITE RECORDS FOR NEW BLUE BOOK
** 5-29-03
*      WRITE(NLUN2,6012) ISSN,H95CM,V95CM
* 6012 FORMAT(6X,'*91*',I4,6X,F10.2,F10.2,40X)
        if(nlscl) accode = 'Y'
        if(.not.nlscl) accode = 'N'   
*     WRITE(NLUN2,6012) ISSN,H95CM,V95CM,accode
*6012 FORMAT(6X,'*91*',I4,6X,F10.2,F10.2,4X,a1)

C  Added at 05/20/2012

        write (SSN_char,'(i4)') ISSN
        do i=1,4
          if (SSN_char(i:i) == ' ') SSN_char(i:i) = '0'
        enddo

C  Convert corren to character to replace the leading zero with a space

        write (dummy,'(f10.8)') corren
        if (dummy(1:1) == '0') dummy(1:1) = ' '

C  Write the *91* record

        WRITE (NLUN2,6012) PIDs(ISN),SSN_char,SD_dn*100.d0,SD_de*100.d0,
     &        dummy,SD_U*100.d0,accode
 6012   FORMAT(a6,'*91*',a4,6X,F10.2,F10.2,a10,f10.2,4X,a1)

C  End addition at 05/20/2012

*********COMPUTE LOCAL ACCURACIES FOR EACH CONNECTION

        NLOCAL = 0
        SLOCH  = 0.0
        SLOCV  = 0.0
        do ii=1,MXSSN
          MedH(ii) = 0.d0
          MedV(ii) = 0.d0
        enddo

  100   READ(NLUN1) I,JSN  
        IF(I.NE.ISN) THEN
          WRITE(LUNIT,6100) I,ISN
 6100     FORMAT(' TABLE ERROR IN NLLIST', 2I6)
          CALL ABORT2
        ENDIF
        IF (JSN.GT.0) THEN
          IF (.NOT. LOCSSN(JSN,JSSN)) THEN
            WRITE(LUNIT,669) ADJFIL, JSN
 669        FORMAT ('SSN TABLE ERROR IN ',A8,' --',I8)
            CALL ABORT2
          ENDIF
          NAME2 =NAMES(JSN)
          CALL rellip(isn,jsn,sigmx,sigmn,theta,a,nx,b,vardn1,varde1,
     &                 coven1)
          SD_dn1  = sqrt(vardn1)
          SD_de1  = sqrt(varde1)

C  Added at 05/20/2012

          corren1 = coven1/SD_dn1/SD_de1
          if (NLSCL) then
            SD_dn1 = SD_dn1*SIGUWT
            SD_de1 = SD_de1*SIGUWT
          endif

C  End addition at 05/20/2012

          CALL CIRC95(SIGMX,SIGMN,H95)

C  Vertical local accuracy of the rover node of the baseline

          IZ1=IUNSTA(JSN,3)
          IF ( .NOT. GETA(IZ1, IZ1, CZ1, A, NX) ) THEN
            WRITE (LUNIT,667) ADJFIL, ISN, IZ1
            CALL ABORT2
          ENDIF
          VARU1=CZ1
          SD_U1 = sqrt(VARU1)

C  The covariance between the vertical of the base and the vertical of the rover         

          IF ( .NOT. GETA(IZ, IZ1, COVV , A, NX) ) THEN
            WRITE (LUNIT,667) ADJFIL, ISN, IZ,IZ1
            CALL ABORT2
          ENDIF
          VARU12 = CZ + CZ1 - 2.d0*COVV
c         write (*,'(f8.3,2x,f8.3,2x,f8.3)') VARU,COVV,VARU1
          SD_U1  = sqrt(VARU12)

C  Added at 05/20/2012

          if (NLSCL) then
            SD_U1 = SD_U1*SIGUWT
          endif

C  End addition at 05/20/2012

          v95=q4*dsqrt(varu12)
          IF(NLSCL) THEN
            H95=H95*SIGUWT
            V95=V95*SIGUWT
          ENDIF
          H95CM=H95*100.
          V95CM=V95*100.
c         WRITE(LUNIT,6020) JSN,JSSN,NAME2,H95CM,V95CM
c6020     FORMAT(41X,'LOCAL',2I6,2X,A30,T100,2F8.2)

C  Compute the distance (baseline length) between the two stations related to these local errors

          CALL GETECX(XPOS1,ISN,B)
          CALL GETECY(YPOS1,ISN,B)
          CALL GETECZ(ZPOS1,ISN,B)
          CALL GETECX(XPOS2,JSN,B)
          CALL GETECY(YPOS2,JSN,B)
          CALL GETECZ(ZPOS2,JSN,B)
          DX = XPOS2 - XPOS1
          DY = YPOS2 - YPOS1
          DZ = ZPOS2 - ZPOS1
          DIST = 0.001d0*sqrt(DX**2 + DY**2 + DZ**2)   !Distance in KM

C  Write a line about the from station in the local accuracy file

          if (JSSN > ISSN) then
            NLOCAL=NLOCAL+1
          else
            goto 100
          endif
          if (NLOCAL == 1) WRITE(NLUN5,6011) ISN,ISSN,NAME1
 6011     FORMAT('FROM    ',2I5, 1X,A30,5x,'HORIZ    UP       DISTANCE')
*
C  Write local accuracies to scratch file

**v6.2.1
          if (JSSN > ISSN) then
**so far for v6.2.1
            WRITE(NLUN5,6020) JSN,JSSN,NAME2,H95CM,V95CM,DIST
 6020       FORMAT('  TO    ',2I5,1X,A30,1x,2F8.2,1x,'CM',1x,f7.2,1x,
     &           'KM')
**v6.2.1
          endif
**so far for v6.2.1

********WRITE RECORDS FOR NEW BLUE BOOK
** 5-29-03
*        WRITE(NLUN2,6022) ISSN,JSSN,H95CM,V95CM
* 6022   FORMAT(6X,'*92*',I4,1X,I4,1X,F10.2,F10.2,40X)
*        WRITE(NLUN2,6022) ISSN,JSSN,H95CM,V95CM,accode
*6022    FORMAT(6X,'*92*',I4,1X,I4,1X,F10.2,F10.2,4X,a1)

C  Added at 05/20/2012

         if (ISSN < JSSN) then
           write (SSN_char,'(i4)') ISSN
           do i=1,4
             if (SSN_char(i:i) == ' ') SSN_char(i:i) = '0'
           enddo
           write (SSN_char1,'(i4)') JSSN
           do i=1,4
             if (SSN_char1(i:i) == ' ') SSN_char1(i:i) = '0'
           enddo

C  Convert corren to character to replace the leading zero with a space

           write (dummy1,'(f10.8)') corren1
           if (dummy1(1:1) == '0') dummy1(1:1) = ' '

           WRITE(NLUN3,6022) PIDs(ISN),SSN_char,SSN_char1,SD_dn1*100.d0,
     &                SD_de1*100.d0,dummy1,SD_U1*100.d0,accode,PIDs(JSN)
 6022      FORMAT(a6,'*92*',a4,2X,a4,2X,F10.2,F10.2,a10,F10.2,4X,a1,3x,
     &            a6)
         endif

C  End addition at 05/20/2012
*****************
          if (JSSN < ISSN) goto 100          
c         SLOCH=SLOCH+H95CM
c         SLOCV=SLOCV+V95CM
          MedH(NLOCAL) = H95CM
          MedV(NLOCAL) = V95CM
          MedD(NLOCAL) = DIST
          GO TO 100
        ELSEIF(JSN.NE.-ISN) THEN
          WRITE(LUNIT,6202) ISN,JSN
 6202     FORMAT(' TABLE2 ERROR IS NLLIST',2I6)
          CALL ABORT2
        ENDIF
C        write (*,*) "NLOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOCAL",NLOCAL
        IF (NLOCAL.GT.0) THEN
c         AVELH=SLOCH/FLOAT(NLOCAL)
c         AVELV=SLOCV/FLOAT(NLOCAL)
          call MDIAN1(MedH,NLOCAL,Hmedian)
          call MDIAN1(MedV,NLOCAL,Vmedian)
          call MDIAN1(MedD,NLOCAL,Dmedian)
c         WRITE (LUNIT,6030) AVELH,AVELV
c         WRITE (LUNIT,6030) Hmedian,Vmedian 
          dashes = "-------------------------------------------"
          WRITE (NLUN5,6029) dashes
          WRITE (NLUN5,6030) Hmedian,Vmedian,Dmedian
          WRITE (NLUN5,*) '   '
c       else
c         WRITE (NLUN5,*) '   '
        ENDIF
c6030   FORMAT(T42,'LOCAL AVERAGE',T100,2F8.2)
 6029   FORMAT(a80)                                                     
 6030   FORMAT(T1,'MEDIAN OF LOCAL ACCURACY',T51,2F8.2,1x,'CM',f8.2,1x,
     &         'KM')

c       WRITE(LUNIT,6031)
 6031   FORMAT()
 
      enddo     

      CLOSE(NLUN1)
      CLOSE(NLUN2)
      CLOSE(NLUN3)
      CLOSE(NLUN4)
      CLOSE(NLUN5)

      RETURN
      END
C----------------------------------------------------------------------------------------------------------
      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=500000)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=l-1
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C----------------------------------------------------------------------------------------
      SUBROUTINE sort2_one_real8_one_char80 (n,arr,brr)
      INTEGER n,M,NSTACK
      REAL*8     arr(n)
      character  brr(n)*80
      PARAMETER  (M=7,NSTACK=50000)
      INTEGER    i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8     a,temp1
      character  b*80,temp2*80    

      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=l-1
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp1=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp1
        temp2=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp2
        if(arr(l).gt.arr(ir))then
          temp1=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp1
          temp2=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp2
        endif
        if(arr(l+1).gt.arr(ir))then
          temp1=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp1
          temp2=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp2
        endif
        if(arr(l).gt.arr(l+1))then
          temp1=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp1
          temp2=brr(l)
          brr(l)=brr(l+1)
          brr(l+1)=temp2
        endif
        i=l+1
        j=ir
        a=arr(l+1)
        b=brr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp1=arr(i)
        arr(i)=arr(j)
        arr(j)=temp1
        temp2=brr(i)
        brr(i)=brr(j)
        brr(j)=temp2
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        brr(l+1)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C----------------------------------------------------------------------------------------
      subroutine MDIAN1(X,N,XMED)

      implicit real*8 (a-h,o-z)

      dimension   X(N)

      call sort (N,X)
      N2 = N/2
      if (2*N2 == N) then
        XMED = 0.5d0*(X(N2) + X(N2+1))
      else
        XMED = X(N2+1)
      endif

      return
      end
C-----------------------------------------------------------------------------------------
      SUBROUTINE NLBB (IUNIT,IUO,IOBS,B,NX,FATAL)
********************************************************************************
* FORM OBS EQ. FOR NON-GPS, NON-DOPPLER BB RECORDS
********************************************************************************
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

*** LOOP OVER RECORDS OF BLUE BOOK

  100 READ (IUNIT,2,END=777) BCARD
  2   FORMAT (A80)
      READ (BCARD,5) IRT
  5   FORMAT (7X,A2)

      LSN = .FALSE.

      IF (IRT .EQ. '20'  .OR.  IRT .EQ. '22') THEN
          IF ( .NOT. LDIR) CALL NL21 (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
        ELSEIF (IRT .EQ. '30'  .OR.  IRT .EQ. '32') THEN
          IF ( .NOT. LANG) CALL NL3 (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
        ELSEIF (IRT .EQ. '40'  .OR.  IRT .EQ. '42') THEN
          IF ( .NOT. LZEN) CALL NL21 (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
        ELSEIF (IRT .EQ. '52'  .OR.  IRT .EQ. '54') THEN
          IF ( .NOT. LDIS) CALL NL22 (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
        ELSEIF (IRT .EQ. '60') THEN
          IF ( .NOT. LAZI) CALL NL21 (BCARD,IUO,IOBS,B,NX,FATAL,LSN)
      ENDIF
      GO TO 100
*
  777 CONTINUE
      RETURN
      END

      SUBROUTINE NLGPS (IUNIT, IUO, IOBS, B, NX, G, FATAL)
********************************************************************************
* DETERMINE STATION CONNECTIONS FROM G FORMAT
* FOR COMPUTATION OF NETWORK AND LOCAL ACCURACIES
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MAXB = 19999 )
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
      LOGICAL ADDCON
      CHARACTER*80 GCARD
      CHARACTER*1 ID, TC
      DIMENSION B(*), NX(*), G(*)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /GPS/    MAXVEC
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /PECHO/  LEB, LLB, LEG, LLG, LED, LLD
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT
      COMMON/BREC / NCCNT(MAXB)

*** TEST FOR MAXIMUM NUMBER OF VECTORS ALLOWED IN GROUPS

      IF (MAXVEC .GT. NVECS) THEN
        WRITE (LUNIT,222) MAXVEC, NVECS
 222    FORMAT (' THE NUMBER OF VECTORS IN A GROUP =',I3,
     &          ' HAS EXCEEDED THE MAXIMUM ALLOCATED =',I3,/,
     &      ' *** FATAL -- INCREASE NVECS IN PARAMETER STATEMENTS ***'/)
        CALL ABORT2
      ENDIF

*** ENTER PROCESSING LOOP--EXIT ON END OF FILE

      NOB = 0

  100 READ (IUNIT,1,END=777) GCARD
    1 FORMAT (A80)
**v6.2.1
      if (GCARD(1:1) == ' ') goto 100
**so far for v6.2.1
      READ (GCARD,5) ID
    5 FORMAT (A1)

*** GROUP HEADER RECORD
      IF (ID.EQ.'A') GO TO 100
      IF (ID .EQ. 'B') THEN
        NOB = NOB + 1
        READ (GCARD,10) IYR1, IMO1, IDY1, IHR1, IMN1,
     &                  IYR2, IMO2, IDY2, IHR2, IMN2
   10   FORMAT (BZ,1X,2(I4,4I2))
        READ (GCARD,110) NVEC, I83
  110   FORMAT (25X,I2,24X,I2)
 
*** ASSIGN C REC COUNT TO NVEC IF B REC CC 26-27 IS BLANK

        IF ( NVEC .LT. 1) NVEC = NCCNT(NOB)
        IF( NVEC .LT. 1  .OR.  NVEC .GT. MAXVEC) THEN
          WRITE (LUNIT,11) GCARD, NVEC, MAXVEC
   11     FORMAT (5X,A80,5X,'NVEC =',I5,' EXCEEDS MAXVEC=',I5)
          CALL ABORT2
        ENDIF

        CALL NLGRP (NVEC, IUNIT, G, ISNS,LISN)

*** EXTRA MEMBER RECORDS

      ELSEIF (ID.EQ.'C'.OR.ID.EQ.'F') THEN
        WRITE (LUNIT,20) GCARD
   20   FORMAT ('0TOO MANY VECTOR RECORDS'/
     &          ' BAD GFILE STRUCTURE--',A80)
        CALL ABORT2

*** CORRELATION AND COVARIANCE RECORDS

      ELSEIF (ID .EQ. 'D') THEN
      ELSEIF (ID .EQ. 'E') THEN
      ENDIF

*** ADD THE CONNECTION
      IF ( .NOT. ADDCON(ISNS,LISN,NX)) THEN
        WRITE (LUNIT,667)
 667    FORMAT ('0INSUFFICIENT STORAGE FOR BLUEBOOK CONNNECTIONS'/)
        CALL ABORT2
      ENDIF

      GO TO 100

*** END OF PROCESSING--END OF FILE ENCOUNTERED

  777 CONTINUE
      RETURN
      END

      SUBROUTINE NLGRP (NVEC, IUNIT, G, ISNS,LISN)
********************************************************************************
* LOAD AN OBSERVATION GROUP INTO WORK SPACE
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      LOGICAL FULL, LBV, GETSSN
      LOGICAL LEB, LLB, LEG, LLG, LED, LLD
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      CHARACTER*80 GCARD
      CHARACTER*1 ID, CODE
      DIMENSION ISNS(*)
      DIMENSION G(*)
      COMMON /PECHO/  LEB, LLB, LEG, LLG, LED, LLD
      COMMON /ECHO/   VSD, LBV
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON/UNITS/ LUNIT
      LOGICAL ISNADD

********  INITIAL LIST OF STATIONS IN THIS SESSION
      DO 5 I=1,NSTA
    5   ISNS(I)=0    
      LISN=0


*** LOOP OVER NVEC MEMBER RECORDS  ('C' AND/OR  'F')

      DO 10 I = 1, NVEC
  100   READ (IUNIT,2,END=666) ID,GCARD
    2   FORMAT (A1,T1,A80)
**v6.2.1
        if (GCARD(1:1) == ' ') goto 100
**so far for v6.2.1
        IF (ID .EQ. 'C') THEN
          READ (GCARD,3) ISSN,JSSN,DX,SDX,DY,SDY,DZ,SDZ,CODE
    3     FORMAT (1X,2I4,3(F11.4,F5.4),A1)
        ELSEIF (ID.EQ.'F') THEN
          READ (GCARD,4) ISSN,JSSN,DX,SDX,DY,SDY,DZ,SDZ,CODE
    4     FORMAT (1X,2I4,3(F13.4,F5.4),A1)
        ELSEIF (ID.EQ.'G'.OR.ID.EQ.'H'.OR.ID .EQ.'D'.OR.ID.EQ.'E'
     *    .OR.ID.EQ.'I') THEN
          GO TO 100
        ELSE
          GO TO 666
        ENDIF

        IF ( .NOT. GETSSN(ISSN,ISN) ) THEN
          CALL LINE (3)
          WRITE (LUNIT,7) GCARD
    7     FORMAT ('0NO *80* RECORD FOR--     ',A80/)
          write (lunit,*) 'It is ISSN'
          CALL ABORT2
        ELSEIF ( .NOT. GETSSN(JSSN,JSN) ) THEN
          CALL LINE (3)
          WRITE (LUNIT,7) GCARD
          write (lunit,*) 'It is JSSN'
          CALL ABORT2
        ENDIF
*
        IF (.NOT. ISNADD(ISN,ISNS,LISN,NSTA))THEN
          WRITE(LUNIT,664) ISN
  664 FORMAT('0ERROR ADDING CONNECTION IN NLGRP',I5)
          CALL ABORT2
        ENDIF
        IF (.NOT. ISNADD(JSN,ISNS,LISN,NSTA))THEN
          WRITE(LUNIT,664) JSN
          CALL ABORT2
        ENDIF
   10 CONTINUE

*** FELL THRU LOOP--NORMAL END OF PROCESSING
      RETURN

*** BAD GFILE STRUCTURE

  666 WRITE (LUNIT,667) GCARD
  667 FORMAT ('0BAD GFILE STRUCTURE--',A80)
      CALL ABORT2
      RETURN
      END


      SUBROUTINE NLLIST( B, G, NX, LNWORK)
********************************************************************************
* READ DATA AND WRITE FIRST OBSERVATION EQUATIONS
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL FATAL
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      DIMENSION B(*), G(*), NX(*)
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7 ADJFIL, NAFILE
      CHARACTER*26 TBUF
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      LOGICAL NLMODE,NLSCL
      COMMON /NLOPT/ NLMODE,NLSCL,NLUN1,NLUN2,NLUN3,NLUN4,NLUN5
      INTEGER GFIRST,GNEXT
      COMMON/UNITS/ LUNIT

*** INITIALIZE CONNECTION MATRIX

      N = NSTA
      CALL OPENG (NSTA, NX, LNWORK)

*** BLUE-BOOK OBSERVATIONS

      IUNIT = 12
      OPEN (IUNIT,STATUS='OLD',FILE=BBOOK)
      CALL NLBB (IUNIT, IUO, IOBS, B, NX, FATAL)
      CLOSE (IUNIT)

*** GFILE OBSERVATIONS

      IF ( .NOT. LGPS) THEN
        IUNIT = 13
        OPEN (IUNIT,ERR=2,STATUS='OLD',FILE=GFILE,IOSTAT=IOS)
        CALL NLGPS (IUNIT, IUO, IOBS, B, NX, G, FATAL)
        CLOSE (IUNIT)
      ENDIF
    2 CONTINUE

*** DFILE OBSERVATIONS

      IF ( .NOT. LDOP) THEN
        IUNIT = 14
        OPEN (IUNIT,ERR=4,STATUS='OLD',FILE=DFILE,IOSTAT=IOS)
*        CALL NLDOP (IUNIT, IUO, IOBS, B, NX, FATAL)
        CLOSE (IUNIT)
      ENDIF
    4 CONTINUE
*
*  WRITE THE NEIGHBOR LIST OUT TO A SCRATCH FILE
      OPEN (NLUN1,STATUS='SCRATCH',FORM='UNFORMATTED')

      DO 200 ISTA=1,NSTA
        IPTR=GFIRST(ISTA,NX)
*** IF IPTR IS NEGATIVE THEN END OF LIST IS REACHED
   50   IF (IPTR.LT.0) THEN 
          GO TO 100
        ELSE
          JSTA=GNEXT(IPTR,NX)
          WRITE(NLUN1) ISTA,JSTA
*          WRITE(LUNIT,6010) ISTA,JSTA
          GO TO 50
        ENDIF
  100 CONTINUE
      WRITE(NLUN1) ISTA,IPTR
*     WRITE(LUNIT,6010) ISTA,IPTR
 6010 FORMAT(2I6)
  200 CONTINUE
      REWIND NLUN1
*      CALL DUMPT(NX,LUNIT)
      RETURN
      END

C---------------------------------------------------------------------------------
      subroutine pellip (isn,sigmx,sigmn,theta,a,nx,b,vardn,varde,coven)
********************************************************************************
* compute point error ellipse
********************************************************************************
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      logical prop,propcv
      dimension a(*),nx(*),b(*)
      dimension icn(2),cn(2)
      dimension ice(2),ce(2)
      common/units/ lunit

*** dummy variables -- placeholders for formic(), and formc()

      iaux=0
      igrt=0

***  north is kind 1

      kind=1
      call formic(kind,isn,isn,iaux,icn,leng,igrt)
      call formc (kind,cn,b,isn,isn,iaux,igrt)
      if(.not.prop(cn,icn,leng,vardn,a,nx,iflag)) then
        if(iflag.eq.1) then
          write(lunit,*) '0state error 1 in rellip'
c         write(*    ,*) '0state error 1 in rellip'
          call abort2
        else
          write(lunit,*) '0profile error 1 in rellip'
c         write(*    ,*) '0profile error 1 in rellip'
          call abort2
        endif
      endif

***  east is kind 2

      kind=2
      call formic(kind,isn,isn,iaux,ice,leng2,igrt)
      call formc(kind,ce,b,isn,isn,iaux,igrt)
      if(.not.prop(ce,ice,leng2,varde,a,nx,iflag)) then
        if(iflag.eq.1) then
          write(lunit,*) '0state error 2 in rellip'
c         write(*    ,*) '0state error 2 in rellip'
          call abort2
        else
          write(lunit,*) '0profile error 2 in rellip'
c         write(*    ,*) '0profile error 2 in rellip'
          call abort2
        endif
      endif

*** compute covariance

      if(.not.propcv(ce,ice,leng,cn,icn,leng2,coven,a,nx,iflag)) then
        if(iflag.eq.1) then
          write(lunit,*) '0state error 3 in rellip'
c         write(*    ,*) '0state error 3 in rellip'
          call abort2
        else
          write(lunit,*) '0profile error 3 in rellip'
c         write(*    ,*) '0profile error 3 in rellip'
          call abort2
        endif
      endif

*** compute error ellipse parameters
*** theta counter-clockwise from local horizon x axis (radians)
*** note: vardn is on y axis, varde is on x axis
*** cf. mikhail

      t1=(varde+vardn)/2.D0
      diff=varde-vardn
      t2=dsqrt(diff*diff/4.D0+coven*coven)
      varmx=t1+t2
      varmn=t1-t2
      sigmx=dsqrt(varmx)
      sigmn=dsqrt(varmn)
      theta2=datan2((coven+coven),(varde-vardn))
      theta=theta2/2.D0
      
      return
      end

C*******************************************************************************
      subroutine rellip(isn,jsn,sigmx,sigmn,theta,a,nx,b,vardn,
     &       varde,coven)
********************************************************************************
* compute relative error ellipse
********************************************************************************
      implicit double precision(a-h,o-z)
      implicit integer(i-n)

      logical prop,propcv
      dimension a(*),nx(*),b(*)
      dimension icn(2),cn(2)
      dimension ice(2),ce(2)
      common/units/ lunit


*** dummy variables -- placeholders for formic(), and formc()

      iaux=0
      igrt=0

*** delta north is kind 21

      kind=21
      call formic(kind,isn,jsn,iaux,icn,leng,igrt)
      call formc(kind,cn,b,isn,jsn,iaux,igrt)
      if(.not.prop(cn,icn,leng,vardn,a,nx,iflag)) then
        if(iflag.eq.1) then
          write(lunit,*) '0state error 1 in rellip'
c         write(*    ,*) '0state error 1 in rellip'
          call abort2
        else
          write(lunit,*) '0profile error 1 in rellip'
c         write(*    ,*) '0profile error 1 in rellip'
          call abort2
        endif
      endif

*** delta east is kind 22

      kind=22
      call formic(kind,isn,jsn,iaux,ice,leng2,igrt)
      call formc(kind,ce,b,isn,jsn,iaux,igrt)
      if(.not.prop(ce,ice,leng2,varde,a,nx,iflag)) then
        if(iflag.eq.1) then
          write(lunit,*) '0state error 2 in rellip'
c         write(*    ,*) '0state error 2 in rellip'
          call abort2
        else
          write(lunit,*) '0profile error 2 in rellip'
c         write(*    ,*) '0profile error 2 in rellip'
          call abort2
        endif
      endif

*** compute covariance

      if(.not.propcv(ce,ice,leng,cn,icn,leng2,coven,a,nx,iflag)) then
        if(iflag.eq.1) then
          write(lunit,*) '0state error 3 in rellip'
c         write(*    ,*) '0state error 3 in rellip'
          call abort2
        else
          write(lunit,*) '0profile error 3 in rellip'
c         write(*    ,*) '0profile error 3 in rellip'
          call abort2
        endif
      endif

*** compute error ellipse parameters
*** theta counter-clockwise from local horizon x axis (radians)
*** note: vardn is on y axis, varde is on x axis
*** cf. mikhail

      t1=(varde+vardn)/2.d0
      diff=varde-vardn
      t2=dsqrt(diff*diff/4.d0+coven*coven)
      varmx=t1+t2
      varmn=t1-t2
      sigmx=dsqrt(varmx)
      sigmn=dsqrt(varmn)
      theta2=datan2((coven+coven),(varde-vardn))
      theta=theta2/2.d0
      
      return
      end

C***********************************************************************************
      SUBROUTINE UPNL(INEW)

      implicit integer (i-n)
      parameter (MXSSN = 9999)

      CHARACTER*80 BCARD
      LOGICAL      NLMODE,NLSCL
      
      real*8       distance                            
      real*8       dist(MXSSN)
      character    cards(MXSSN)*80

      COMMON /NLOPT/ NLMODE,NLSCL,NLUN1,NLUN2,NLUN3,NLUN4,NLUN5
      COMMON/UNITS/ LUNIT

** v 4.28

      OPEN(NLUN2,file='NLUN2',FORM='FORMATTED')
      OPEN(NLUN3,file='NLUN3',FORM='FORMATTED')
      OPEN(NLUN4,file='NLUN4',FORM='FORMATTED')
      OPEN(NLUN5,file='NLUN5',FORM='FORMATTED')

      ncount = 0

      REWIND NLUN2
  100 READ(NLUN2,'(a80)',ERR=666,END=700) BCARD
      WRITE(INEW,'(a80)') BCARD
      GO TO 100

  700 REWIND NLUN3
  101 READ(NLUN3,'(a80)',ERR=666,END=701) BCARD
      WRITE(INEW,'(a80)') BCARD
      GO TO 101

  701 REWIND NLUN4
  102 READ(NLUN4,'(a80)',ERR=666,END=702) BCARD
      WRITE(LUNIT,'(a80)') BCARD
      GO TO 102

  702 WRITE(LUNIT,'(///)')
      REWIND NLUN5
c 103 READ(NLUN5,'(a80)',ERR=666,END=703) BCARD
c     WRITE(LUNIT,'(a80)') BCARD
c     GO TO 103

  103 READ(NLUN5,'(a80)',ERR=666,END=703) BCARD
      if (BCARD(1:6)==' LOCAL'.or.BCARD(1:4)=='FROM'.or.
     &    BCARD(1:6)=='      ') then
        WRITE(LUNIT,'(a80)') BCARD
      elseif (BCARD(1:4) == '  TO') then
        ncount = ncount + 1
        cards(ncount) = BCARD
        read (BCARD(70:77),'(f8.2)') distance      
        dist(ncount) = distance
      elseif (BCARD(1:6) == 'MEDIAN') then
        call sort2_one_real8_one_char80 (ncount,dist,cards)
        do ii=1,ncount
          write (LUNIT,'(a80)') cards(ii)
        enddo
c       WRITE(LUNIT,'(a80)') BCARD
        ncount = 0
      endif
      GO TO 103

  666 WRITE (LUNIT,667)
  667 FORMAT ('0NO FILE OF NETWORK AND LOCAL ACCURACIES FOUND', /)
      CALL ABORT2
      RETURN

***********3-3-04
  703 close(nlun2)
      close(nlun3)
      close(nlun4)
      close(nlun5)

C  Delete scratch files

c#ifdef NGS_UNIX_ENV
c      call system ('rm NLUN2')
c      call system ('rm NLUN3')
c      call system ('rm NLUN4')
c      call system ('rm NLUN5')
c#else
      call system ('del NLUN2')
      call system ('del NLUN3')
      call system ('del NLUN4')
      call system ('del NLUN5')
c#endif

      return
      END
C*********************************************************************************
