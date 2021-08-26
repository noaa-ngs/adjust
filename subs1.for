* %P%
* file subs1.for  ver $RCSfile: subs1.for,v $
  
      SUBROUTINE AB (A, B, R, L, M, N)

***  FORM THE MATRIX PRODUCT R=AB
***  THE MATRICES A AND B ARE RETURNED UNCHANGED

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION A(L,M), B(M,N), R(L,N)
      CHARACTER*99 SCCSID

C      SCCSID='$Id: subs1.for 107411 2018-10-30 14:54:00Z jarir.saleh $	20$Date: 2008/08/25 16:04:51 $ NGS'

      DO 5 I = 1, L
        DO 4 J = 1, N
          R(I,J) = 0.D0
          DO 3 K = 1, M
            R(I,J) = R(I,J) + A(I,K)*B(K,J)
    3     CONTINUE
    4   CONTINUE
    5 CONTINUE

      RETURN
      END
      SUBROUTINE ABAT (A, B, R, W, M, N)

*** FORM THE MATRIX PRODUCT R=ABA'
*** THE MATRICES A AND B ARE RETURNED UNCHANGED
*** W IS A SCRATCH VECTOR

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION A(M,N), B(N,N), R(M,M), W(N)

      DO 15 I = 1, M

        DO 5 K = 1, N
          W(K) = 0.D0
          DO 4 L = 1, N
            W(K) = W(K) + A(I,L)*B(L,K)
    4     CONTINUE
    5   CONTINUE

        DO 14 J = I, M
          R(I,J) = 0.D0
          DO 10 L = 1, N
            R(I,J) = R(I,J) + W(L)*A(J,L)
   10     CONTINUE
          R(J,I) = R(I,J)
   14   CONTINUE

   15 CONTINUE

      RETURN
      END
      SUBROUTINE ABORT2

*** PRINT MESSAGE OF FATAL TERMINATION

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*26 TBUF
      COMMON/UNITS/ LUNIT

      WRITE (LUNIT,1)
    1 FORMAT ('0')
      WRITE (LUNIT,2)
    2 FORMAT (1X, 130('*') )
      WRITE (LUNIT,3)
    3 FORMAT (' ***** FATAL TERMINATION -- FATAL TERMINATION !!'/
     *        ' ***** THIS DUMP IS INTENTIONAL !'/
     *        ' ***** REFER TO PRIOR ERROR MESSAGES')
      WRITE (LUNIT,2)
      ITYPE = 0
      CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
      WRITE (LUNIT,4) TBUF(1:24)
    4 FORMAT ('0SYSTEM TIME IS ', A24)

      STOP 666
      END
      SUBROUTINE ACCUR (A, NX, B, SHIFTS, SIGUWT)

*** COMPUTE AND LIST ACCURACIES
*** NOTE: ORTHOMETRIC OR GEOID HEIGHT UPDATED IF 2.5DIM ADJUSTMENT

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      CHARACTER*80 ACARD
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7  ADJFIL, NAFILE
      CHARACTER*2 CC12, CC34
      LOGICAL GETSSN
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      LOGICAL ELFLAG, DFFLAG
      DIMENSION A(*), NX(*), B(*), SHIFTS(*)
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /FLAGS/  ELFLAG(MXSSN), DFFLAG(MXSSN)
      COMMON /GPSORD/ CC34
      COMMON/UNITS/ LUNIT

*** DO NOT COMPUTE ACCURACIES UNLESS INVERSE COMPUTED

      IF (IMODE .NE. 1) THEN
        IUNIT = 11
        OPEN (IUNIT,ERR=200,STATUS='OLD',FILE=AFILE,IOSTAT=IOS)
        CALL HEAD4

***** UPDATE HEIGHT IF 2.5DIM ADJUSTMENT

        IF (L2HLF) THEN
          CALL LINE4 (4)
          WRITE (LUNIT,10)
  10      FORMAT (/, ' WARNING - FOR RELATIVE LENGTH ACCURACY',
     &               ' COMPUTATIONS, DUAL HEIGHT STATIONS',
     &            /, ' ARE TREATED AS 2-DIM STATIONS',
     &               ' AND USE THE ORIGINAL STATION HEIGHT', /)
        ENDIF

***** LOOP OVER ADJUSTMENT FILE FOR ACCURACY CARDS (QQ)

  100   READ (IUNIT,1,END=777) ACARD
    1   FORMAT (A80)
        READ (ACARD,5) CC12
    5   FORMAT (A2)

        IF (CC12 .EQ. 'QQ') THEN
          READ (ACARD,2) ISSN, JSSN
    2     FORMAT (10X, I4, 36X, I4)
          CC34(1:2) = ACARD(3:4)
          IF (CC34(1:1) .EQ. ' '  .AND.  CC34(2:2) .NE. ' ') THEN
            CC34(1:1) = CC34(2:2)
            CC34(2:2) = ' '
          ENDIF
          IF (GETSSN( ISSN, ISN)  .AND.  GETSSN(JSSN, JSN) ) THEN
            CALL RELACC (ISSN, JSSN, ISN, JSN, A, NX, B, SHIFTS, SIGUWT)
          ENDIF
        ENDIF
        GO TO 100

  777   CLOSE (IUNIT)
        CALL LINE4 (2)
        WRITE (LUNIT,3)
    3   FORMAT ('0************ END OF ACCURACY PROCESSING **********')
      ENDIF
  200 CONTINUE

      RETURN
      END

c The new argument IOBSNO is being passed to this subroutine so that
c it can be passed to the VSHIFT subroutines.  This observation number
c reflects the numbering of vector components in sequence, without
c skipping rejected vectors.  Mike Potterfield 3/15/07.

      SUBROUTINE ACUM20 (IOBS, VSD,IOBSNO)

*** ACCUMULATE A LARGE RESIDUAL  (NORMALIZED)

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
      COMMON /REDUN/  RN1, RN2, RN3, RN4, RN5, RN6, RN7, RN8, RN9,
     &                RN10, RN11, RN12

*** DEAL ONLY WITH ABSOLUTE VALUES .GT. SMALLEST

      V = DABS(VSD)
      IF (V .LE. VSD20(MXRD) ) RETURN

*** POINT TO NEW ARRAY LOCATION

      IP = IPOINT(V)

*** SHIFT ARRAY CONTENTS AND LOAD ARRAYS

      CALL VSHIFT (IP, IOBS, V, IOBSNO)

      RETURN
      END
      
      
      
**v 4.28
c The new argument IOBSNO is being passed to this subroutine so that
c it can be passed to the VSHIFT subroutines.  This observation number
c reflects the numbering of vector components in sequence, without
c skipping rejected vectors.  Mike Potterfield 3/15/07.

      SUBROUTINE ACUM20V (IOBS, V, IOBSNO)

*** ACCUMULATE A LARGE RESIDUAL

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRDV = 20 )
      COMMON /RESTA3/ V20(MXRDV), I20V(MXRDV),NI20V                

      COMMON /UNITS/ LUNIT
      
*** DEAL ONLY WITH ABSOLUTE VALUES .GT. SMALLEST

      V2 = DABS(V)
 
      IF (V2 .LE. V20(MXRDV) ) RETURN

*** POINT TO NEW ARRAY LOCATION

      IP = IPOINTV(V2)

*** SHIFT ARRAY CONTENTS AND LOAD ARRAYS

      CALL VSHIFTV (IP, IOBS, V2, IOBSNO)

      RETURN
      END

******************

** v 4.29f
c The new argument IOBSNO is being passed to this subroutine so that
c it can be passed to the VSHIFT subroutines.  This observation number
c reflects the numbering of vector components in sequence, without
c skipping rejected vectors.  Mike Potterfield 3/15/07.

      SUBROUTINE ACUM2DU (IOBS, V, IOBSNO)

*** ACCUMULATE A LARGE RESIDUAL FOR DU COMPONENT

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRDU = 20 )
      COMMON /RESTA7/ VDU20(MXRDU), I20DU(MXRDU),NI20DU                

      COMMON /UNITS/ LUNIT
      
*** DEAL ONLY WITH ABSOLUTE VALUES .GT. SMALLEST

      V2 = DABS(V)
 
      IF (V2 .LE. VDU20(MXRDU) ) RETURN

*** POINT TO NEW ARRAY LOCATION

      IP = IPOINTU(V2)

*** SHIFT ARRAY CONTENTS AND LOAD ARRAYS

      CALL VSHIFTU (IP, IOBS, V2, IOBSNO)

      RETURN
      END
*********************

** v 4.29i
c The new argument IOBSNO is being passed to this subroutine so that
c it can be passed to the VSHIFT subroutines.  This observation number
c reflects the numbering of vector components in sequence, without
c skipping rejected vectors.  Mike Potterfield 3/15/07.

      SUBROUTINE ACUM2DL (IOBS, V, IOBSNO)

*** ACCUMULATE A LARGE RESIDUAL FOR DL COMPONENT

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXRDL = 20 )
      COMMON /RESTA8/ VDL20(MXRDL), I20DL(MXRDL),NI20DL                

      COMMON /UNITS/ LUNIT
      
*** DEAL ONLY WITH ABSOLUTE VALUES .GT. SMALLEST

      V2 = DABS(V)
 
      IF (V2 .LE. VDL20(MXRDL) ) RETURN

*** POINT TO NEW ARRAY LOCATION

      IP = IPOINTL(V2)

*** SHIFT ARRAY CONTENTS AND LOAD ARRAYS

      CALL VSHIFTL (IP, IOBS, V2, IOBSNO)

      RETURN
      END

***********************
      SUBROUTINE ADJAUX (A, NX, GOOGE, B, SIGUWT)

*** LIST ADJUSTED AUXILIARY PARAMETERS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL GETA
      DIMENSION A(*), NX(*), GOOGE(*), B(*)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

*** HEADING

      CALL HEAD5
      NLINE = 2
      IF (IMODE .NE. 1) NLINE = NLINE + 1
      IF (LPG) NLINE = NLINE + 1

*** LOOP OVER PARAMETERS

      DO 2 IAUX = 1, NAUX
        CALL GETAUX (AUX, IAUX, B)
        CALL LINE5 (NLINE)
        WRITE (LUNIT,3) IAUX, AUX
    3   FORMAT ('0#', I5, ' AUXILIARY PARAMETER IS ', 1PD10.2)

*** SIGMA

        IF (IMODE .NE. 1) THEN
          IA = IUNAUX(IAUX)
          IF ( .NOT. GETA(IA, IA, SA, A, NX) ) THEN
            WRITE (LUNIT,667) IAUX, IA
  667       FORMAT ('0GET SIGMA ERROR IN ADJAUX--', 2I8)
            CALL ABORT2
          ENDIF
          SA = DSQRT(SA)
          IF (LSS) THEN
            SA = SA*SIGUWT
            AUXSA = AUX/SA
            WRITE (LUNIT,5) SA, AUXSA
    5       FORMAT (17X, 'SCALED SIGMA =', 1PD10.2, 7X, 'VAL/SD=',
     &              0PF7.1)
          ELSE
            AUXSA = AUX/SA
            WRITE (LUNIT,9) SA, AUXSA
    9       FORMAT (15X, 'UNSCALED SIGMA =', 1PD10.2, 7X, 'VAL/SD=',
     &              0PF7.1)
          ENDIF
        ENDIF

*** GOOGE

        IF (LPG) THEN
          IA = IUNAUX(IAUX)
          WRITE (LUNIT,6) GOOGE(IA)
    6     FORMAT (23X, 'GOOGE = ', 1PD10.2)
        ENDIF

    2 CONTINUE

      RETURN
      END
      SUBROUTINE ADJGRT (A, NX, GOOGE, SIGUWT)

*** LIST ADJUSTED GPS AND DOPPLER ROTATION AUXILIARY PARAMETERS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL GETA
      DIMENSION A(*), NX(*), GOOGE(*)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /CONST/  PI, PI2, RAD
      COMMON/UNITS/ LUNIT

*** HEADING

      CALL HEAD8
      NLINE = 2
      IF (IMODE .NE. 1) NLINE = NLINE + 1
      IF (LPG) NLINE = NLINE + 1

*** LOOP OVER PARAMETERS - CHANGE FROM RADIANS TO ARC SECONDS

      RADSEC = ( 60.D0 * 60.D0 * 180.D0 ) / PI
      DO 2 IGRT = 1, NGRT
        CALL GETRTG (GRTX, GRTY, GRTZ, IGRT)
        GRTX = GRTX * RADSEC
        GRTY = GRTY * RADSEC
        GRTZ = GRTZ * RADSEC
        CALL LINE8 (NLINE)
        WRITE (LUNIT,3) IGRT, GRTX, GRTY, GRTZ
    3   FORMAT ('0#', I3, 9X, 'VALUE =', 1PD10.2, 15X, 1PD10.2, 15X,
     &          1PD10.2, '   (ARC SECONDS)')

*** SIGMA

        IF (IMODE .NE. 1) THEN
          IA = IUNGRT(IGRT, 1)
          IF ( .NOT. GETA(IA, IA, SAX, A, NX) ) THEN
            WRITE (LUNIT,667) IGRT, IA
  667       FORMAT ('0GET SIGMA ERROR IN ADJGRT--', 2I8)
            CALL ABORT2
          ENDIF
          IA = IUNGRT(IGRT, 2)
          IF ( .NOT. GETA(IA, IA, SAY, A, NX) ) THEN
            WRITE (LUNIT,667) IGRT, IA
            CALL ABORT2
          ENDIF
          IA = IUNGRT(IGRT, 3)
          IF ( .NOT. GETA(IA, IA, SAZ, A, NX) ) THEN
            WRITE (LUNIT,667) IGRT, IA
            CALL ABORT2
          ENDIF
          SAX = DSQRT(SAX) * RADSEC
          SAY = DSQRT(SAY) * RADSEC
          SAZ = DSQRT(SAZ) * RADSEC
          IF (LSS) THEN
            SAX = SAX*SIGUWT
            SAY = SAY*SIGUWT
            SAZ = SAZ*SIGUWT
            GRTSAX = GRTX/SAX
            GRTSAY = GRTY/SAY
            GRTSAZ = GRTZ/SAZ
            WRITE (LUNIT,5) SAX, SAY, SAZ, GRTSAX, GRTSAY, GRTSAZ
    5       FORMAT (7X, 'SCALED SIGMA =', 1PD10.2, 15X, 1PD10.2, 15X,
     &              1PD10.2, /, 13X, 'VAL/SD =', 0PF7.1, 18X, 0PF7.1,
     &              18X, 0PF7.1)
          ELSE
            GRTSAX = GRTX/SAX
            GRTSAY = GRTY/SAY
            GRTSAZ = GRTZ/SAZ
            WRITE (LUNIT,9) SAX, SAY, SAZ, GRTSAX, GRTSAY, GRTSAZ
    9       FORMAT (5X, 'UNSCALED SIGMA =', 1PD10.2, 15X, 1PD10.2,
     &              15X, 1PD10.2,
     &              /, 13X, 'VAL/SD =', 0PF7.1, 18X, 0PF7.1, 18X,0PF7.1)
          ENDIF
        ENDIF

*** GOOGE

        IF (LPG) THEN
          IAX = IUNGRT(IGRT, 1)
          IAY = IUNGRT(IGRT, 2)
          IAZ = IUNGRT(IGRT, 3)
          WRITE (LUNIT,6) GOOGE(IAX), GOOGE(IAY), GOOGE(IAZ)
    6     FORMAT (14X, 'GOOGE =', 1PD10.2, 15X, 1PD10.2, 15X, 1PD10.2)
        ENDIF

    2 CONTINUE

      RETURN
      END
C--------+---------+---------+---------+---------+---------+---------+------------------------
      SUBROUTINE ADJPOS (A, B, NX, SHIFTS, GOOGE, SIGUWT)

*** LIST ADJUSTED POSITIONS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LOCSSN, GETA
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      LOGICAL ELFLAG, DFFLAG
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*30 NAMES, NAME
      CHARACTER*7 ADJFIL, NAFILE
      CHARACTER*1 ADIR1, ADIR2, ADX, ADY, AEHT
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
      COMMON/UNITS/ LUNIT

**v6.3

      CHARACTER   la_s1*7,lo_s1*7
      CHARACTER   la_s2*1,lo_s2*1
      character   lat_lon_char*25,HT_char*7
      logical     un_constrained_HT,un_constrained_Hz
      LOGICAL     overwrite_const_coord_in_newbb
      character   CC_records*80,hemi*1
      character   CHAR_DX*8,CHAR_DY*8,CHAR_DZ*8
      character   CHAR_ADX*1,CHAR_ADY*1
      character   CHAR_IAZ*3
      character   CHAR_DH*7,CHAR_TOT*7 
      COMMON/MM6/overwrite_const_coord_in_newbb
      COMMON/MM6_array/CC_records(MXSSN)

**so far for v6.3

*** HEADING

      IF (LAP) THEN
        CALL HEAD3
        NLINE = 2
        IF (LPS) NLINE = NLINE + 1
        IF (IMODE .GE. 2) NLINE = NLINE + 1
        IF (LPG) NLINE = NLINE + 1
      ENDIF

*** OPEN ADJUSTED POSITION SUMMARY FILE

      IF (LADJ) THEN
        LUS = 23
        OPEN (LUS, FILE=ADJFIL, STATUS='UNKNOWN')
        WRITE (LUS, 102)
  102   FORMAT (' ADJUSTED POSITION SUMMARY', ///,
     &          '  SSN   LATITUDE    LONGITUDE     GMSL     GHT',
     &          '      EHT      DLAT     DLON     DHT    AZ   DH',
     &          '   DTOT  SIG LA   SIG LO   SIG HT',
     &          '   GOOGE X  GOOGE Y  GOOGE Z   NAME')
      ENDIF

*** GET POSITIONS

      DO 2 ISN = 1, NSTA
        CALL GETGLA (GLAT, ISN, B)
        CALL GETGLO (GLON, ISN, B)
        CALL GETMSL (GMSL0, ISN, B)
        CALL GETECX (XPOS,ISN,B)
        CALL GETECY (YPOS,ISN,B)
        CALL GETECZ (ZPOS,ISN,B)
        CALL GETGH  (GHT0, ISN, B)
        IF (L2HLF  .AND.  I2HLF(ISN) .EQ. 1) THEN
c         write (LUNIT,*) 'GEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEZ'
          CALL GETEHT (EHT, ISN, B)
          IF ( ELFLAG(ISN) ) THEN
c           write (LUNIT,*) 'TEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEZ'
            GMSL = EHT - GHT0
            GHT = GHT0
          ELSE
            GHT = EHT - GMSL
            GMSL = GMSL0
          ENDIF
          AEHT = '*'
        ELSE
          GMSL = GMSL0
          GHT = GHT0
          EHT = GMSL + GHT
          AEHT = ' '
        ENDIF
        NAME = NAMES(ISN)
        CALL GETDMS (GLAT, ID1, IM1, S1, ISIGN)
        IF (ISIGN .GT. 0) THEN
          ADIR1 = 'N'
        ELSE
          ADIR1 = 'S'
        ENDIF
        DLAT = GLAT*RAD*ISIGN
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
        DLON = GLON*RAD*ISIGN
        IF ( .NOT. LOCSSN(ISN, ISSN) ) THEN
          WRITE (LUNIT,666) ADJFIL, ISN
  666     FORMAT ('0SSN TABLE ERROR IN ', A8, ' --', I8)
          CALL ABORT2
        ENDIF

        IF (LAP) THEN
          CALL LINE3 (NLINE)

**v6.3 Replace adjusted coordinates by constrained ones in the ADJUST output file

c         if (overwrite_const_coord_in_newbb) then
c           do ii=1,MXSSN
c             if (CC_records(ii)(1:2) == 'CC') then
c               read (CC_records(ii)(11:14),'(i4)') iissn
c               if (ISSN == iissn) then
c                 read(CC_records(ii)(45:69),'(a25)') lat_lon_char
c                 read(CC_records(ii)(70:76),'(a7)') HT_char
c                 un_constrained_Hz = (lat_lon_char == ' ')
c                 un_constrained_HT = (HT_char      == ' ')

c                 if (.not. un_constrained_Hz) then
c                   if(.not. un_constrained_HT) then          !The station is constrained includin the height
c                     read(CC_records(ii)(45:76),301) 
c    &                ID1,IM1,la_s1,la_s2,ID2,IM2,lo_s1,lo_s2,iEHT
c                     do iii=1,7
c                       if(la_s1(iii:iii)==' ') la_s1(iii:iii)='0'
c                       if(lo_s1(iii:iii)==' ') lo_s1(iii:iii)='0'
c                     enddo
c                     read (la_s1,'(f7.5)') s1
c                     read (lo_s1,'(f7.5)') s2

C  Latitude

c                     if (la_s2 == 'N') then
c                       dlat = ID1+(IM1/60.D0)+(S1/3600.D0)
c                     else
c                       dlat = -(ID1+(IM1/60.D0)+(S1/3600.D0))
c                     endif
c                     glat  = dlat/RAD
c                     ADIR1 = la_s2

C  Longitude

c                     if (lo_s2 == 'W') then
c                       dlon = ID2+(IM2/60.D0)+(S2/3600.D0)
c                     else
c                       ALONGD=ID2+(IM2/60.D0)+(S2/3600.D0)
c                       ALONGD=360.D0-ALONGD
c                       ID2=INT(ALONGD)
c                       ALONGM=(ALONGD-ID2) * 60.D0
c                       IM2=INT(ALONGM)
c                       S2=(ALONGM-IM2) * 60.D0
c                       dlon = ID2+(IM2/60.D0)+(S2/3600.D0)
c                       lo_s2 = 'W'
c                     endif
c                     glon = (360.d0 - dlon)/RAD
c                     ADIR2 = lo_s2

c  Height

c                     EHT = 0.001d0*iEHT

c  Cartesian coordinates

c                     call TOXYZ (GLAT,GLON,EHT,XPOS,YPOS,ZPOS)

c                     exit
c                   else         !The station's Hz position is constrained but not the height 
c                     read(CC_records(ii)(45:76),301) 
c    &                ID1,IM1,la_s1,la_s2,ID2,IM2,lo_s1,lo_s2
c                     do iii=1,7
c                       if(la_s1(iii:iii)==' ') la_s1(iii:iii)='0'
c                       if(lo_s1(iii:iii)==' ') lo_s1(iii:iii)='0'
c                     enddo
c                     read (la_s1,'(f7.5)') s1
c                     read (lo_s1,'(f7.5)') s2

c  Latitude

c                     if (la_s2 == 'N') then
c                       dlat = ID1+(IM1/60.D0)+(S1/3600.D0)
c                     else
c                       dlat = -(ID1+(IM1/60.D0)+(S1/3600.D0))
c                     endif
c                     glat   = dlat/RAD
c                     ADIR1  = la_s2

c  Longitude

c                     if (lo_s2 == 'W') then
c                       dlon = ID2+(IM2/60.D0)+(S2/3600.D0)
c                     else
c                       ALONGD=ID2+(IM2/60.D0)+(S2/3600.D0)
c                       ALONGD=360.D0-ALONGD
c                       ID2=INT(ALONGD)
c                       ALONGM=(ALONGD-ID2) * 60.D0
c                       IM2=INT(ALONGM)
c                       S2=(ALONGM-IM2) * 60.D0
c                       dlon = ID2+(IM2/60.D0)+(S2/3600.D0)
c                       lo_s2 = 'W'
c                     endif
c                     glon = (360.d0 - dlon)/RAD
c                     ADIR2 = lo_s2

c  Height

c                     !EHT = EHT      Take the adjusted ellipsoidal height of the station

c  Cartesian coordinates

c                     call TOXYZ (GLAT,GLON,EHT,XPOS,YPOS,ZPOS)

c                     exit
c                   endif
c                 else
c                   read(CC_records(ii)(70:76),'(i7)') iEHT
c                   EHT = 0.001d0*iEHT
c                   !glat and glon are the adjusted ones
c                   call TOXYZ (GLAT,GLON,EHT,XPOS,YPOS,ZPOS)
c                   exit
c                 endif
c               endif
c             else
c               cycle
c             endif
c           enddo
 301        format (2i2,a7,a1,i3,i2,a7,a1,i7)
c         endif

**so far for v6.3

          WRITE (LUNIT,3) ISN, ISSN, NAME, ID1, IM1, S1, ADIR1,
     *               ID2, IM2, S2, ADIR2, GMSL, GHT, EHT, AEHT
    3     FORMAT ('0', I5, I5, 1X, A30,    2I3, F9.5, A1,
     *             1X, I4, I3, F9.5, A1, 2X, F8.3, F9.3, F10.3, A1)
        ENDIF

*** PRINT SHIFTS IF ENABLED

        IF (LPS) THEN
          IF (IMODE .NE. 0) THEN
            I1 = IUNSHF(ISN, 1)
            I2 = IUNSHF(ISN, 2)
            I3 = IUNSHF(ISN, 3)
            DX = SHIFTS(I1)
            DY = SHIFTS(I2)
            DZ = SHIFTS(I3)
            DXS = DX
            DYS = DY
            DZS = DZ
            DX2 = DX*DX
            DY2 = DY*DY
            DH = DSQRT(DX2 + DY2)
            IF (DABS(DX) .LE. 1.D-3  .OR.  DABS(DY) .LE. 1.D-3) THEN
              IAZ = 0
            ELSE
              IAZ = IDNINT( DATAN2(DY, DX)*RAD + 0.5D0 )
              IF (IAZ .LT. 0) IAZ = IAZ + 360
            ENDIF
            DTOT = DSQRT(DX2 + DY2 + DZ*DZ)
            IF (DX .LT. 0.D0) THEN
              DX = -DX
              ADX = 'S'
            ELSE
              ADX = 'N'
            ENDIF
            IF (DY .LT. 0.D0) THEN
              DY = -DY
              ADY = 'W'
            ELSE
              ADY = 'E'
            ENDIF

**v6.3 Replace actual shifts by asterisks for constrained stations at which the adjusted positions 
**     were replaced (this is weighted constraints) by their constrained coordinates.

            icounter = 0
c           if (overwrite_const_coord_in_newbb) then
c             do ii=1,MXSSN
c               if (CC_records(ii)(1:2) == 'CC') then
c                 read (CC_records(ii)(11:14),'(i4)') iissn
c                 if (ISSN == iissn) then
c                   read(CC_records(ii)(45:69),'(a25)') lat_lon_char
c                   read(CC_records(ii)(70:76),'(a7)') HT_char
c                   un_constrained_Hz = (lat_lon_char == ' ')
c                   un_constrained_HT = (HT_char      == ' ')
c                   if (.not. un_constrained_Hz) then
c                     if(.not. un_constrained_HT) then          !The station is constrained includin the height
c                       icounter = icounter + 1
c                       CHAR_DX  = '********' 
c                       CHAR_ADX = '*'
c                       CHAR_DY  = '********' 
c                       CHAR_ADY = '*'
c                       CHAR_DZ  = '********' 
c                       CHAR_IAZ = '***'
c                       CHAR_DH  = '*******' 
c                       CHAR_TOT = '*******' 
c                       IF (LAP) then
c                         WRITE(LUNIT,41) CHAR_DX,CHAR_ADX,CHAR_DY,
c    &                                    CHAR_ADY,CHAR_DZ,CHAR_IAZ, 
c    &                                    CHAR_DH, CHAR_TOT
c41                       FORMAT(32X,'SHIFTS (M.)',4X,a8,A1,9X,a8,A1,
c    &                           22X,a8,' AZ=',a3,' H=',a7,' T=',a7)
c                       endif
c                       exit
c                     else
c                       icounter = icounter + 1
c                       CHAR_DX  = '********' 
c                       CHAR_ADX = '*'
c                       CHAR_DY  = '********' 
c                       CHAR_ADY = '*'
c                       CHAR_IAZ = '***'
c                       CHAR_DH  = '*******' 
c                       DTOT     = DZ                            
c                       IF (LAP) then
c                         WRITE(LUNIT,43) CHAR_DX,CHAR_ADX,CHAR_DY,
c    &                                    CHAR_ADY,DZ,CHAR_IAZ, 
c    &                                    CHAR_DH, DTOT
c43                       FORMAT(32X,'SHIFTS (M.)',4X,a8,A1,9X,a8,A1,
c    &                        22X,f8.3,' AZ=',a3,' H=',a7,' T=',f7.3)
c                       endif
c                       exit
c                     endif
c                   else
c                     icounter = icounter + 1
c                     CHAR_DZ  = '********' 
c                     DTOT     = DSQRT(DX2 + DY2)
c--------+---------+---------+---------+---------+---------+---------+------------------------
c                     IF (LAP) then
c                       WRITE(LUNIT,42) DX,ADX,DY,ADY,CHAR_DZ,IAZ,DH,
c    &                                  DTOT
c42                     FORMAT(32X,'SHIFTS (M.)',4X,f8.3,A1,9X,f8.3,A1,
c    &                         22X,a8,' AZ=',I3,' H=',f7.3,' T=',f7.3)
c                     endif
c                     exit
c                   endif
c                 else
c                   cycle
c                 endif
c               else
c                 cycle
c               endif
c             enddo
c           endif  
            if (icounter == 0) then
              IF (LAP) WRITE(LUNIT,4) DX,ADX,DY,ADY,DZ,IAZ,DH,DTOT
 4            FORMAT(32X,'SHIFTS (M.)',4X,F8.3,A1,9X,F8.3,A1,22X,
** v 4.29g
**   *              F8.3, ' AZ=', I3, ' HOR=', F5.1, ' TOT=', F5.1)
     *              F8.3, ' AZ=', I3, ' H=', F7.3, ' T=', F7.3)
            endif
********************     
          ENDIF
        ENDIF

** so far for v6.3

*** PRINT INVERSE ELEMENTS IF ENABLED

        IF (IMODE .NE. 1) THEN
          IX = IUNSTA(ISN, 1)
          IY = IUNSTA(ISN, 2)
          IZ = IUNSTA(ISN, 3)
          IF ( .NOT. GETA(IX, IX, SX, A, NX) ) THEN
            WRITE (LUNIT,667) ADJFIL, ISN, IX
  667       FORMAT ('0GET SIGMA ERROR IN ', A8, ' --', 2I8)
            CALL ABORT2
          ENDIF
          IF ( .NOT. GETA(IY, IY, SY, A, NX) ) THEN
            WRITE (LUNIT,667) ISN, IY
            CALL ABORT2
          ENDIF
          IF ( .NOT. GETA(IZ, IZ, SZ, A, NX) ) THEN
            WRITE (LUNIT,667) ISN, IZ
            CALL ABORT2
          ENDIF
          SX = DSQRT(SX)
          SY = DSQRT(SY)
          SZ = DSQRT(SZ)
          IF (LSS) THEN
            SX = SX*SIGUWT
            SY = SY*SIGUWT
            SZ = SZ*SIGUWT
            IF (LAP) WRITE (LUNIT,5) SX, SY, SZ
    5       FORMAT (25X, 'SCALED SIGMAS (M.)', 4X, F8.3, 10X, F8.3, 23X,
     &              F8.3)
          ELSE
            IF (LAP) WRITE (LUNIT,9) SX, SY, SZ
    9       FORMAT (23X, 'UNSCALED SIGMAS (M.)', 4X, F8.3, 10X, F8.3,
     *              T97, F8.3)
          ENDIF
        ENDIF

*** PRINT GOOGE NUMBERS IF ENABLED

        IF (LPG) THEN
          IX = IUNSTA(ISN, 1)
          IY = IUNSTA(ISN, 2)
          IZ = IUNSTA(ISN, 3)
          GX = GOOGE(IX)
          GY = GOOGE(IY)
          GZ = GOOGE(IZ)
          IF (LAP) WRITE (LUNIT,6) GX, GY, GZ
    6     FORMAT (32X, 'GOOGES', 11X, 1PD8.1, 10X, 1PD8.1, 23X, 1PD8.1)
        ENDIF

*** PRINT ADJUSTED POSITION OUTPUT FILE, ENSURE THAT NO * ARE
*** PRINTED, REGARDLESS OF THE OPERATING SYSTEM

        IF (LADJ) THEN
          IF ( DXS .GE. 1.D4  .OR.  DXS .LE. -1.D3  .OR.
     &         DABS(DXS) .LT. 0.001D0 ) DXS = 0.D0
          IF ( DYS .GE. 1.D4  .OR.  DYS .LE. -1.D3  .OR.
     &         DABS(DYS) .LT. 0.001D0 ) DYS = 0.D0
          IF ( DZS .GE. 1.D4  .OR.  DZS .LE. -1.D3  .OR.
     &         DABS(DZS) .LT. 0.001D0 ) DZS = 0.D0
          IF ( DH .GE. 1.D3  .OR.  DH .LE. -1.D2  .OR.
     &         DABS(DH) .LT. 0.1D0 ) DH = 0.D0
          IF ( DTOT .GE. 1.D3  .OR.  DTOT .LE. -1.D2  .OR.
     &         DABS(DTOT) .LT. 0.1D0 ) DTOT = 0.D0
          IF ( SX .GE. 1.D4  .OR.  SX .LE. -1.D3  .OR.
     &         DABS(SX) .LT. 0.001D0 ) SX = 0.D0
          IF ( SY .GE. 1.D4  .OR.  SY .LE. -1.D3  .OR.
     &         DABS(SY) .LT. 0.001D0 ) SY = 0.D0
          IF ( SZ .GE. 1.D4  .OR.  SZ .LE. -1.D3  .OR.
     &         DABS(SZ) .LT. 0.001D0 ) SZ = 0.D0
          IF ( GX .GE. 1.D6  .OR.  GX .LE. -1.D5  .OR.
     &         DABS(GX) .LT. 0.1D0 ) GX = 0.D0
          IF ( GY .GE. 1.D6  .OR.  GY .LE. -1.D5  .OR.
     &         DABS(GY) .LT. 0.1D0 ) GY = 0.D0
          IF ( GZ .GE. 1.D6  .OR.  GZ .LE. -1.D5  .OR.
     &         DABS(GZ) .LT. 0.1D0 ) GZ = 0.D0

          WRITE (LUS, 100) ISSN, DLAT, DLON,
     &                    GMSL, GHT, EHT, AEHT,
     &                    DXS, DYS, DZS,
     &                    IAZ, DH, DTOT,
     &                    SX, SY, SZ,
     &                    GX, GY, GZ, NAME
  100     FORMAT (I5, ':', F11.7, ':', F12.7, ':',
     &            F8.3, ':', F8.3, ':', F8.3, A1, ':',
     &            F8.3, ':', F8.3, ':', F8.3, ':',
     &            I3, ':', F5.1, ':', F5.1, ':',
     &            F8.3, ':', F8.3, ':', F8.3, ':',
     &            1PD8.1, ':', 1PD8.1, ':', 1PD8.1, ':', 2X, A30)
        ENDIF

*** PRINT OUT THE EARTH CENTERED COORDINATES

        WRITE(LUNIT,33) XPOS,YPOS,ZPOS
 33     FORMAT(1X,T17,'X = ',F15.4,10X,'Y = ',F15.4,
     *         10X,'Z = ',F15.4)
        CALL LINE3 (1)
    2 CONTINUE

      IF (LADJ) CLOSE (LUS)

      RETURN
      END
C-------------------------------------------------------------------------------------------------
**v6.3
      SUBROUTINE constrained_stations_shifts (A, B, NX)

*** LIST N-E-U SHIFTS OF CONSTRAINED STATIONS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      LOGICAL LOCSSN, GETA
      CHARACTER*30 NAMES, NAME
      DIMENSION A(*), B(*), NX(*)
      double precision R(3,3)
      double precision DXYZ(3),DNEU(3)
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV

      CHARACTER   la_s1*7,lo_s1*7
      CHARACTER   la_s2*1,lo_s2*1
      character   CC_records*80
      character   lat_lon_char*25,HT_char*7
      logical     un_constrained_HT,un_constrained_Hz

      COMMON /OPRINT/ CRIT,LBB,LGF,LCS,LVD,LVA,LVZ,LVS,LVR,LVG,
     &                LVC,LIS,LPS,LPG,LDR,LOS,LAP,LDF,LVX,LGV
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /STRUCT/NSTA,NAUX,NUNK,IDIM,NSTAS,NOBS,NCON,NZ,NGRT
      COMMON /CONST/  PI, PI2, RAD
      COMMON/UNITS/ LUNIT
      COMMON/MM6_array/CC_records(MXSSN)

*** HEADING

c     CALL HEAD3
      NLINE = 2
      IF (LPS) NLINE = NLINE + 1
      IF (IMODE .GE. 2) NLINE = NLINE + 1
      IF (LPG) NLINE = NLINE + 1

C  write header

      write (LUNIT,'(///)')
      write (LUNIT,660) 
660   format
     &('COMPARISON OF ADJUSTED AND CONSTRAINED COORDINATES AND HEIGHTS')
      write (LUNIT,661) 
661   format('ADJUSTED MINUS CONSTRAINED COORDINATE SHIFTS (METERS)')
      write (LUNIT,6611) 
6611  format('CURRENTLY NO VALUES ARE PRINTED FOR VERTICAL CONSTRAINTS'/
     &)
      write(LUNIT,662)'ISN','SSN','STATION DESIGNATION','DELTA NORTH',
     &'DELTA EAST','DELTA HORIZ','DELTA UP'  
662   format (1x,a3,3x,a3,10x,a19,18x,a11,4x,a10,3x,a11,3x,a8)

*** GET POSITIONS

      DO ISN = 1, NSTA
        CALL GETECX (XPOS,ISN,B)
        CALL GETECY (YPOS,ISN,B)
        CALL GETECZ (ZPOS,ISN,B)
        CALL GETGLA (GLAT,ISN,B)
        CALL GETGLO (GLON,ISN,B)
        CALL GETEHT (EHT ,ISN,B)

        NAME = NAMES(ISN)

        IF ( .NOT. LOCSSN(ISN, ISSN) ) THEN
          WRITE (LUNIT,666) ADJFIL, ISN
  666     FORMAT ('0SSN TABLE ERROR IN ', A8, ' --', I8)
          CALL ABORT2
        ENDIF

c       CALL LINE3 (NLINE)

        do ii=1,MXSSN
          if (CC_records(ii)(1:2) == 'CC') then
            read (CC_records(ii)(11:14),'(i4)') iissn
            if (ISSN == iissn) then
              read(CC_records(ii)(45:69),'(a25)') lat_lon_char
              read(CC_records(ii)(70:76),'(a7)') HT_char

C  This modification was done on 09/17/2014
C  The following two lines exclude vertically constrained stations from being listed in the ADJUST.out
C  file in the "ADJUSTED MINUS CONSTRAINED COORDINATE SHIFTS" table.
C  This is the last thing I did for Michael Dennis before he leaves in 09/2014.
C  The reason this was done is that this table was messed up for vertical constrains.
C  while it looks perfectly fine for all other constraints.

              if (CC_records(ii)(77:77) /= 'E') then
                write (LUNIT,6671) ISN,ISSN,NAME,"Vertical constraints"
6671            format (i4,2x,i4,8x,a30,8x,a20)                           
                cycle                 
              endif

C  End of modification on 09/17/2014

              un_constrained_Hz = (lat_lon_char == ' ')
              un_constrained_HT = (HT_char      == ' ')

              if (.not. un_constrained_Hz) then
                if(.not. un_constrained_HT) then          !The station is constrained includin the height
                  read(CC_records(ii)(45:76),301) 
     &            ID1,IM1,la_s1,la_s2,ID2,IM2,lo_s1,lo_s2,iEHT
                  do iii=1,7
                    if (la_s1(iii:iii) == ' ') la_s1(iii:iii) = '0'
                    if (lo_s1(iii:iii) == ' ') lo_s1(iii:iii) = '0'
                  enddo
                  read (la_s1,'(f7.5)') s1
                  read (lo_s1,'(f7.5)') s2

C  Latitude

                  if (la_s2 == 'N') then
                    dlat = ID1+(IM1/60.D0)+(S1/3600.D0)
                  else
                    dlat = -(ID1+(IM1/60.D0)+(S1/3600.D0))
                  endif
                  glat  = dlat/RAD

C  Longitude

                  if (lo_s2 == 'W') then
                    dlon = ID2+(IM2/60.D0)+(S2/3600.D0)
                  else
                    ALONGD=ID2+(IM2/60.D0)+(S2/3600.D0)
                    ALONGD=360.D0-ALONGD
                    ID2=INT(ALONGD)
                    ALONGM=(ALONGD-ID2) * 60.D0
                    IM2=INT(ALONGM)
                    S2=(ALONGM-IM2) * 60.D0
                    dlon = ID2+(IM2/60.D0)+(S2/3600.D0)
                    lo_s2 = 'W'
                  endif
                  glon = (360.d0 - dlon)/RAD

C  Height

                  EHT = 0.001d0*iEHT

C  Cartesian coordinates

                  call TOXYZ (GLAT,GLON,EHT,XPOS_C,YPOS_C,ZPOS_C)
                  DX   = XPOS - XPOS_C
                  DY   = YPOS - YPOS_C
                  DZ   = ZPOS - ZPOS_C
                  DXYZ(1) = DX
                  DXYZ(2) = DY
                  DXYZ(3) = DZ
                  call BUILDR (R,ISN,B)
          
                  DN   = R(1,1)*DX + R(1,2)*DY + R(1,3)*DZ
                  DE   = R(2,1)*DX + R(2,2)*DY 
                  DU   = R(3,1)*DX + R(3,2)*DY + R(3,3)*DZ
                  DH   = sqrt(DN**2 + DE**2)

                  write (LUNIT,667) ISN,ISSN,NAME,DN,DE,DH,DU
667               format (i4,2x,i4,8x,a30,8x,f11.4,2x,f11.4,2x,f11.4,2x,
     &                f11.4)
                  exit
                else
                  read(CC_records(ii)(45:76),301) 
     &            ID1,IM1,la_s1,la_s2,ID2,IM2,lo_s1,lo_s2
                  do iii=1,7
                    if (la_s1(iii:iii) == ' ') la_s1(iii:iii) = '0'
                    if (lo_s1(iii:iii) == ' ') lo_s1(iii:iii) = '0'
                  enddo
                  read (la_s1,'(f7.5)') s1
                  read (lo_s1,'(f7.5)') s2

C  Latitude

                  if (la_s2 == 'N') then
                    dlat = ID1+(IM1/60.D0)+(S1/3600.D0)
                  else
                    dlat = -(ID1+(IM1/60.D0)+(S1/3600.D0))
                  endif
                  glat  = dlat/RAD

C  Longitude

                  if (lo_s2 == 'W') then
                    dlon = ID2+(IM2/60.D0)+(S2/3600.D0)
                  else
                    ALONGD=ID2+(IM2/60.D0)+(S2/3600.D0)
                    ALONGD=360.D0-ALONGD
                    ID2=INT(ALONGD)
                    ALONGM=(ALONGD-ID2) * 60.D0
                    IM2=INT(ALONGM)
                    S2=(ALONGM-IM2) * 60.D0
                    dlon = ID2+(IM2/60.D0)+(S2/3600.D0)
                    lo_s2 = 'W'
                  endif
                  glon = (360.d0 - dlon)/RAD

C  Height

                  !EHT = EHT       The adjusted height

C  Cartesian coordinates

                  call TOXYZ (GLAT,GLON,EHT,XPOS_C,YPOS_C,ZPOS_C)
                  DX   = XPOS - XPOS_C
                  DY   = YPOS - YPOS_C
                  DZ   = ZPOS - ZPOS_C
                  DXYZ(1) = DX
                  DXYZ(2) = DY
                  DXYZ(3) = DZ
                  call BUILDR (R,ISN,B)
          
                  DN   = R(1,1)*DX + R(1,2)*DY + R(1,3)*DZ
                  DE   = R(2,1)*DX + R(2,2)*DY 
                  DU   = R(3,1)*DX + R(3,2)*DY + R(3,3)*DZ
                  DH   = sqrt(DN**2 + DE**2)

                  write (LUNIT,667) ISN,ISSN,NAME,DN,DE,DH
                  exit
                endif
              else
                read(CC_records(ii)(70:76),'(i7)') iEHT
                !glat and glon are adjusted  
                EHT = 0.001d0*iEHT

C  Cartesian coordinates

                call TOXYZ (GLAT,GLON,EHT,XPOS_C,YPOS_C,ZPOS_C)
                DX   = XPOS - XPOS_C
                DY   = YPOS - YPOS_C
                DZ   = ZPOS - ZPOS_C
                DXYZ(1) = DX
                DXYZ(2) = DY
                DXYZ(3) = DZ
                call BUILDR (R,ISN,B)
          
c               DN   = R(1,1)*DX + R(1,2)*DY + R(1,3)*DZ
c               DE   = R(2,1)*DX + R(2,2)*DY 
                DU   = R(3,1)*DX + R(3,2)*DY + R(3,3)*DZ
c               DH   = sqrt(DN**2 + DE**2)

c  only the vertical is constrained, so we do not need to display the Hz shifts   

                write (LUNIT,668) ISN,ISSN,NAME,DU
668             format (i4,2x,i4,8x,a30,47x,f11.4)
                exit
              endif
            endif
          else
            cycle
          endif
        enddo
 301    format (2i2,a7,a1,i3,i2,a7,a1,i7)
      enddo

      RETURN
      END

**v6.3 so far for v6.3
C-------------------------------------------------------------------------------------------------
      SUBROUTINE AFAA (ACARD)

*** ELLIPSOID PARAMETER RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) AXX, E2X
    1 FORMAT (2X, F10.3, F18.18)
      IF (ACARD(3:12) .NE. '          ') AX = AXX
      IF (ACARD(13:30) .NE. '                  ') E2 = E2X

      IF (AX.LT. 6300000.D0 .OR. AX.GT. 6400000.D0) WRITE (LUNIT,2) AX
 2    FORMAT (' ***********NOTE: SEMIMAJOR AXIS =', F12.3)
      IF (E2 .LT. 0.006D0  .OR.  E2 .GT. 0.007D0) WRITE (LUNIT,3) E2
 3    FORMAT (' ***********NOTE: ECCENTRICITY =', F20.18)

      RETURN
      END
      SUBROUTINE AFBB (ACARD)

*** BYPASS RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (NAUX .GT. 0  .OR.  NGRT .GT. 0) THEN
        WRITE (LUNIT,666)
 666   FORMAT ('0THE BB RECORD MUST PROCEDE THE DD,SS,RR AND VV ',
     &         ' RECORDS')
        CALL ABORT2
      ENDIF

      IF (ACARD(3:3) .NE. ' ') LDIR = .TRUE.
      IF (ACARD(4:4) .NE. ' ') LANG = .TRUE.
      IF (ACARD(5:5) .NE. ' ') LZEN = .TRUE.
      IF (ACARD(6:6) .NE. ' ') LDIS = .TRUE.
      IF (ACARD(7:7) .NE. ' ') LAZI = .TRUE.
      IF (ACARD(8:8) .NE. ' ') LGPS = .TRUE.
      IF (ACARD(9:9) .NE. ' ') LDOP = .TRUE.

      RETURN
      END
      SUBROUTINE AFDD (ACARD)

*** DIMENSIONALITY RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      CHARACTER*80 ACARD
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      LOGICAL L2HLF, LEHT
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

*** GET DIMENSION OF ADJUSTMENT

      IF (ACARD(3:3) .EQ. '1') THEN
        IDIM = 1
        LDIR = .TRUE.
        LANG = .TRUE.
        LAZI = .TRUE.
      ELSEIF (ACARD(3:3) .EQ. '2') THEN
        IDIM = 2
      ELSEIF (ACARD(3:3) .EQ. '3') THEN
        IDIM = 3
      ELSE
        CALL LINE (3)
        WRITE (LUNIT,1) ACARD
    1   FORMAT ('0ILLEGAL DIMENSION --', A80, /)
      ENDIF

*** IF 2-DIM ADJUSTMENT, DETERMINE IF 2.5-DIM ADJUSTMENT

      IF (IDIM .EQ. 2  .AND.  ACARD(5:5) .EQ. '5') THEN
        L2HLF = .TRUE.
      ELSE
        L2HLF = .FALSE.
      ENDIF

**v 4.28*****used in BIGADJUST**
*** IF 2.5-DIM ADJUSTMENT, CHECK IF ELLIP HTS TO BE READ FROM 80 RECS
*
*     IF (L2HLF .AND. ACARD(6:6) .EQ. 'E' ) THEN
*        LEHT = .TRUE.
*     ELSE
*        LEHT = .FALSE.
*     ENDIF
********************************************************

*** IF 2-DIM ADJUSTMENT AND NOT 2.5 DIM, MAKE LZEN TRUE
    
      IF (IDIM .EQ. 2  .AND.  .NOT.L2HLF) THEN
        LZEN = .TRUE.
      ENDIF

      RETURN
      END
      SUBROUTINE AFEE (ACARD)

*** DEFAULT MEAN SEA LEVEL ELEVATION RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI

      READ (ACARD,1) DMSLX
    1 FORMAT (2X, F7.3)
      IF (ACARD(3:9) .NE. '       ') DMSL = DMSLX

      RETURN
      END
      SUBROUTINE AFGG (ACARD)

*** DEFAULT GEOID HEIGHT RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI

      READ (ACARD,1) DGHX
    1 FORMAT (2X, F7.3)
      IF (ACARD(3:9) .NE. '       ') DGH = DGHX

      RETURN
      END
      SUBROUTINE AFHD (ACARD)

***DEFAULT HEIGHT ADJUSTMENT RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /OPT/   AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &               LMSL, LSS, LUP, LABS, LADJ, LUPI

      IF (ACARD(3:3) .EQ. 'E') THEN
        LMSL = .TRUE.
      ELSEIF (ACARD(3:3) .EQ. 'G') THEN
        LMSL = .FALSE.
      ENDIF

      RETURN
      END
      SUBROUTINE AFII (ACARD)

*** ITERATION RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7 ADJFIL, NAFILE
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
***10-24-03
      common /opt2/ itmax2
*********      
      READ (ACARD,1) ITMAXX, VMAXX, VPRTX, CTOLX
    1 FORMAT (2X, I2, F4.0, F3.0, F6.3)

      IF (ACARD(3:4) .NE. '  ') ITMAX = IABS(ITMAXX)
      IF (ACARD(5:8) .NE. '    ') VM = DABS(VMAXX)
      IF (ACARD(9:11) .NE. '   ') VP = DABS(VPRTX)
      IF (VM .LT. VP) VM = VP
      IF (VM .GT. 9998.D0) VM = 99999.0D0
      IF (ACARD(12:17) .NE. '      ') CTOL = DABS(CTOLX)
      IF (ACARD(18:18) .EQ. 'Y') LDR = .TRUE.
      IF (ACARD(19:19) .EQ. 'Y') LUPI = .TRUE.

      IF (VM .GE. 9999.D0) VM = 1.D20
      IF (VP .GE. 999.D0) VP = 1.D20

      IF ( ACARD(20:26) .EQ. '       ' ) THEN
        NAFILE = 'NEWAF'
      ELSE
        READ (ACARD,2) NAFILE
    2   FORMAT (20X, A7)
      ENDIF
***10-24-03
       itmax2=itmax
******

      RETURN
      END
      SUBROUTINE AFILE2

*** ROUTINE TO PROCESS THE A FILE FOR SECOND TIME

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7 ADJFIL, NAFILE
      CHARACTER*2 CC12
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE

*** OPEN ADJUSTMENT FILE (IUNIT=11)

      IF (AFILE .EQ. 'NOAFILE') GOTO 200
      IUNIT = 11
      OPEN (IUNIT,ERR=200,STATUS='OLD',FILE=AFILE,IOSTAT=IOS)

*** READ A FILE -- PROCESS RECORDS WITH SSN FIELDS

   10 READ (IUNIT,11,END=100) ACARD
   11 FORMAT (A80)
      READ (ACARD,12) CC12
   12 FORMAT (A2)

      IF (CC12 .EQ. 'HC') THEN
        CALL ASHC (ACARD)
      ENDIF
      GO TO 10

*** END OF A FILE PROCESSING

  100 CLOSE (IUNIT)

*** A FILE NOT PRESENT

  200 RETURN
      END
C*************************************************************************************
      SUBROUTINE AFMM (ACARD)

*** ADJUSTMENT MODE RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7 ADJFIL, NAFILE
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL overwrite_const_coord_in_newbb    ! The last is for V6.3
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      COMMON /UNITS/ LUNIT
      COMMON /MM6/overwrite_const_coord_in_newbb    !V6.3

      IF (ACARD(3:3) .EQ. '0') THEN
        IMODE = 0
      ELSEIF (ACARD(3:3) .EQ. '1') THEN
        IMODE = 1
      ELSEIF (ACARD(3:3) .EQ. '2') THEN
        IMODE = 2
      ELSEIF (ACARD(3:3) .EQ. '3') THEN
        IMODE = 3
      ENDIF

      IF (ACARD(4:4) .EQ. 'Y') THEN
        LSS = .TRUE.
      ELSE
        LSS = .FALSE.
      ENDIF

      IF (ACARD(5:5) .EQ. 'Y') THEN
        LUP = .TRUE.
      ELSE
        LUP = .FALSE.
      ENDIF

**v6.3

c     IF(ACARD(6:6)=='Y'.or.ACARD(6:6)=='y'.or.ACARD(6:6)==' ') THEN
      IF(ACARD(6:6)=='Y'.or.ACARD(6:6)=='y') THEN
c     IF(ACARD(6:6)/='N'.and.ACARD(6:6)/='n') THEN
        overwrite_const_coord_in_newbb = .TRUE.
c       overwrite_const_coord_in_newbb = .FALSE.              !This is temporary. Michael is leaving for university studies
      ELSE                                                    !Therefore, this was disabled for now. Later, spend some time and fix it.
        overwrite_const_coord_in_newbb = .FALSE.
      ENDIF

**so far for  v6.3

**v 4.28
*     IF (ACARD(6:12) .EQ. '       ') THEN
*       BBNAM = 'NEWBB'
*     ELSE
*       READ (ACARD,2) BBNAM
*   2   FORMAT (5X, A7)
*     ENDIF
*     write(lunit,99) bbnam
*99    format(' afmm- ',a7)     

      IF (ACARD(3:3) .EQ. '0') THEN
        LSS = .FALSE.
        LUP = .FALSE.
      ENDIF

      IF (ACARD(13:13) .EQ. 'N') THEN
        LABS = .FALSE.
      ELSE
        LABS = .TRUE.
      ENDIF

      IF (ACARD(14:14) .EQ. 'Y') THEN
        LADJ = .TRUE.
        IF (ACARD(15:21) .NE. '       ') THEN
          READ (ACARD,3) ADJFIL
    3     FORMAT ( 14X, A7) 
        ELSE
          ADJFIL = 'ADJPOS'
        ENDIF
      ELSE
        LADJ = .FALSE.
      ENDIF

      RETURN
      END
C**************************************************************************************
      SUBROUTINE AFPP (ACARD)

*** PRINT OUTPUT RECORD
* v 4.29 change PP rec format

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      CHARACTER*3 ACRIT
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LEB, LLB, LEG, LLG, LED, LLD
      COMMON /PECHO/  LEB, LLB, LEG, LLG, LED, LLD
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
     
* v 4.29     
      COMMON /OPRIN2/ CRIT2

      READ (ACARD,1) ACRIT
 1    FORMAT (8X, A3)
      CALL NBLANK (ACRIT, 3, IBLK)
      READ (ACRIT,2) CRIT
 2    FORMAT (F3.1)
 
      READ (ACARD,3) ACRIT
 3    FORMAT (5X, A3)
      CALL NBLANK (ACRIT, 3, IBLK)
      READ (ACRIT,2) CRIT2
 

      IF (ACARD( 3: 3) .NE. ' '  .AND.  ACARD( 3: 3) .NE. '0'  .AND.
     &    ACARD( 3: 3) .NE. '1'  .AND.  ACARD( 3: 3) .NE. '2')
     &  LBB = .FALSE.
      IF (ACARD( 3: 3) .EQ. '1') LEB = .TRUE.
      IF (ACARD( 3: 3) .EQ. '2') LLB = .TRUE.
      IF (ACARD( 4: 4) .NE. ' '  .AND.  ACARD( 4: 4) .NE. '0'  .AND.
     &    ACARD( 4: 4) .NE. '1'  .AND.  ACARD( 4: 4) .NE. '2')
     &  LGF = .FALSE.
      IF (ACARD( 4: 4) .EQ. '1') LEG = .TRUE.
      IF (ACARD( 4: 4) .EQ. '2') LLG = .TRUE.
      IF (ACARD( 5: 5) .NE. ' ') LCS = .FALSE.
      IF (ACARD( 12: 12) .NE. ' ') LVD = .FALSE.
      IF (ACARD(13:13) .NE. ' ') LVA = .FALSE.
      IF (ACARD(14:14) .NE. ' ') LVZ = .FALSE.
      IF (ACARD(15:15) .NE. ' ') LVS = .FALSE.
      IF (ACARD(16:16) .NE. ' ') LVR = .FALSE.
      IF (ACARD(17:17) .NE. ' ') LVG = .FALSE.
      IF (ACARD(18:18) .NE. ' ') LVC = .FALSE.
      IF (ACARD(19:19) .NE. ' ') LIS = .FALSE.
      IF (ACARD(20:20) .NE. ' ') LPS = .FALSE.
      IF (ACARD(21:21) .NE. ' ') LPG = .FALSE.
      IF (ACARD(22:22) .NE. ' ') LOS = .FALSE.
      IF (ACARD(23:23) .NE. ' ') LAP = .FALSE.
      IF (ACARD(24:24) .NE. ' '  .AND.  ACARD(24:24) .NE. '0'  .AND.
     &    ACARD(24:24) .NE. '1'  .AND.  ACARD(24:24) .NE. '2')
     &  LDF = .FALSE.
      IF (ACARD(24:24) .EQ. '1') LED = .TRUE.
      IF (ACARD(24:24) .EQ. '2') LLD = .TRUE.
      IF (ACARD(25:25) .NE. ' ') LVX = .FALSE.
      IF (ACARD(26:26) .NE. ' ') LGV = .FALSE.

      IF (ACARD( 6: 8) .EQ. '   ') CRIT2 = 0.D0
      IF (ACARD( 9:11) .EQ. '   ') CRIT = 0.D0


      RETURN
      END
      SUBROUTINE AFRR (ACARD, FATAL)

*** AUXILIARY GPS AND DOPPLER ROTATION PARAMETER RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      CHARACTER*2 ACODE
      CHARACTER*1 TC1, TC2
      LOGICAL PUTGRT, FATAL
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT        
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /CONST/  PI, PI2, RAD
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) ACODE, IYR1, IMO1, IDY1, IHR1, IMN1, TC1,
     &                      IYR2, IMO2, IDY2, IHR2, IMN2, TC2,
     &                      IVALX, IVALY, IVALZ
    1 FORMAT (2X, A2, 2(I4, 4I2, A1), 3I10)

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

***   CONVERT GPS AND DOPPLER ROT. VALUES FROM 1.D-5 ARC SEC. TO RADIANS

      SECRAD = PI / ( 60.D0 * 60.D0 * 180.D0 * 1.D5 )
      IF (ACARD(31:40) .EQ. '          ') THEN
        VALX = 0.D0
      ELSE
        VALX = IVALX * SECRAD
      ENDIF

      IF (ACARD(41:50) .EQ. '          ') THEN
        VALY = 0.D0
      ELSE
        VALY = IVALY * SECRAD
      ENDIF

      IF (ACARD(51:60) .EQ. '          ') THEN
        VALZ = 0.D0
      ELSE
        VALZ = IVALZ * SECRAD
      ENDIF

      IF ( .NOT. PUTGRT(ICODE, IOLD, INEW, IDUP) ) THEN
        IF (IDUP .EQ. 0) THEN
          CALL LINE (3)
          WRITE (LUNIT,3)
    3     FORMAT ('0GPS AND DOPPLER ROTATION PARM TABLE OVERFLOWED',
     &            '--RECORD IGNORED'/)
        ELSE
          CALL LINE (3)
          FATAL = .TRUE.
          WRITE (LUNIT,4) IDUP
    4     FORMAT ('0DUPLICATES THE ', I2, '-TH GPS AND DOPPLER',
     &            ' ROTATION PARMETER RECORD'/)
        ENDIF
      ELSE
        NGRT = NGRT + 1
        CALL PUTRTG (VALX, VALY, VALZ, NGRT)
      ENDIF

      RETURN
      END
      SUBROUTINE AFSS (ACARD, B, FATAL)

*** AUXILIARY PARAMETER RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      LOGICAL PUTPRM, FATAL
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      CHARACTER*1 TC1, TC2
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) ICODE, IYR1, IMO1, IDY1, IHR1, IMN1, TC1,
     &                      IYR2, IMO2, IDY2, IHR2, IMN2, TC2, IVAL
    1 FORMAT (2X, I2, 2(I4, 4I2, A1), I5)
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

      IF ( .NOT. PUTPRM(ICODE, IOLD, INEW, IDUP) ) THEN
        IF (IDUP .EQ. 0) THEN
          CALL LINE (3)
          WRITE (LUNIT,3)
    3     FORMAT ('0AUXILIARY PARM TABLE OVERFLOWED--RECORD IGNORED'/)
        ELSE
          CALL LINE (3)
          FATAL = .TRUE.
          WRITE (LUNIT,4) IDUP
    4     FORMAT ('0DUPLICATES THE ', I2, '-TH AUX. PARAMETER RECORD'/)
        ENDIF
      ELSE
        VAL = DBLE(IVAL)*1.0D-8
        NAUX = NAUX + 1
        CALL PUTAUX (VAL, NAUX, B)
      ENDIF

      RETURN
      END
      SUBROUTINE AFVV (ACARD)

*** VARIANCE FACTOR RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXVF = 40 )
      CHARACTER*80 ACARD
      LOGICAL PUTIVF, FATAL, LFIX
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      CHARACTER*1 TC1, TC2, CON
** v 4.30vf
      CHARACTER*2 CCODE
**********
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON /VFTB2/  VFS(MXVF), VTV(MXVF), VFRN(MXVF), VFSS(MXVF)
      COMMON/UNITS/ LUNIT
** v 4.32vf
       LOGICAL IVCGPS,VSlogic
      COMMON/VCREC/SIGH,SIGU, VFHOR,VFUP, IVCGPS
      CHARACTER*13 VSPROJ
      PARAMETER (MAXPRJ=2500)
      COMMON/VSREC/SIGHS(MAXPRJ), SIGUS(MAXPRJ), ISETHU, VSPROJ(MAXPRJ),
     &             VSlogic
         
      READ (ACARD,101) CCODE
  101 FORMAT(2X,A2)
      IF (CCODE.EQ.'HU') THEN
C   COMPUTE VARIANCE FACTORS FOR GPS HORIZONTAL AND UP COMPONENTS OF OBSERVATIONS
        IVCGPS=.TRUE.
        RETURN
      ENDIF
**********

      READ (ACARD,1) ICODE, IYR1, IMO1, IDY1, IHR1, IMN1, TC1,
     &                      IYR2, IMO2, IDY2, IHR2, IMN2, TC2, IVAL, CON
    1 FORMAT (2X, I2, 2(I4, 4I2, A1), I5, A1)
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
      IF (CON.EQ.'Y') THEN
        LFIX = .TRUE.
      ELSE
        LFIX = .FALSE.
      ENDIF
      IF (IVAL .LE. 0) IVAL = 100

      CALL TOMNT (IYR1, IMO1, IDY1, IHR1, IMN1, TC1, IOLD)
      CALL TOMNT (IYR2, IMO2, IDY2, IHR2, IMN2, TC2, INEW)

      IF ( .NOT. PUTIVF(ICODE, IOLD, INEW, LFIX, IDUP) ) THEN
        IF (IDUP .EQ. 0) THEN
          CALL LINE (3)
          WRITE (LUNIT,3)
    3     FORMAT ('0VARIANCE FACTOR TABLE OVERFLOWED--RECORD IGNORED'/)
        ELSE
          CALL LINE (3)
          WRITE (LUNIT,4) IDUP
    4     FORMAT ('0DUPLICATES THE ', I2, '-TH VARIANCE FACTOR RECORD'/)
          FATAL = .TRUE.
        ENDIF
      ELSE
        VAL = DBLE(IVAL)*0.01D0
        VFS(NVFTOT) = VAL
        VFSS(NVFTOT) = VAL
      ENDIF

      RETURN
      END
      SUBROUTINE ALLXYZ (B, GLAAVE, GLOAVE)

*** COMPUTE X,Y,Z FOR ALL THE STATIONS AND AVERAGE OF NETWORK LOCAT'S

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      DIMENSION B(*)
      LOGICAL L2HLF, LEHT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT

      GLASUM = 0.D0
      GLOSUM = 0.D0

      DO 1 I = 1, NSTA
        CALL GETGLA (GLAT, I, B)
        CALL GETGLO (GLON, I, B)
        CALL GETMSL (GMSL, I, B)
        CALL GETGH (GHT, I, B)
        EHT = GMSL + GHT

        CALL TOXYZ (GLAT, GLON, EHT, X, Y, Z)

        CALL PUTECX (X, I, B)
        CALL PUTECY (Y, I, B)
        CALL PUTECZ (Z, I, B)

        GLASUM = GLASUM + GLAT
        GLOSUM = GLOSUM + GLON

    1 CONTINUE

*** CALCULATE AVERGAE OF LATITUDES, LONGITUDES, AND ELLIP. HEIGHTS

      GLAAVE = GLASUM / NSTA
      GLOAVE = GLOSUM / NSTA

*** IF 2.5DIM ADJUSTMENT, THEN THE X,Y,Z COMPONENTS MUST BE OBTAINED
*** FOR THE ELLIPSOIDAL HEIGHTS

      IF (L2HLF) THEN
        DO 2 I = 1, NSTA
            CALL GETGLA (GLAT, I, B)
            CALL GETGLO (GLON, I, B)
            CALL GETEHT (EHT, I, B)

            CALL TOXYZ (GLAT, GLON, EHT, X, Y, Z)

            CALL PUTEHX (X, I, B)
            CALL PUTEHY (Y, I, B)
            CALL PUTEHZ (Z, I, B)
    2   CONTINUE
      ENDIF

      RETURN
      END
      SUBROUTINE ALOCAT (ID1, ID2, ID3, II4, ID5, LAWORK, LNWORK)

*** ALLOCATE STORAGE

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      LOGICAL L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /GPS/    MAXVEC
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

*** COMPUTE STRUCTURE PARAMETERS

      NSTAS = 3*NSTA
      NUNK = NAUX + (3*NGRT) + NZ + (3*NSTA)
      NOBS = 0
      NCON = 0
      NMAX = NSTA

      IF ( .NOT. L2HLF) THEN
        LENG = (MAXVEC + 1)*IDIM + 1 + 3
      ELSE
        LENG = (MAXVEC + 1)*3 + 1 + 3
      ENDIF
      MAXROW = 3*MAXVEC
*     MAXCOL = MAXROW + LENG + 3
      MAXCOL = MAXROW + LENG + 5
      LENGPS = MAXROW*MAXCOL

*** IF 2.5DIM ADJUSTMENT, THEN LEAVE ROOM FOR EHT, EHX, EHY, EHZ

      IF ( .NOT. L2HLF) THEN
        ID1 = 9*NSTA + NAUX + NZ
      ELSE
        ID1 = 13*NSTA + NAUX + NZ
      ENDIF
      ID2 = ID1 + NUNK
      ID3 = ID2 + IUNSHF(NMAX, 3)
      ID4 = ID3 + LENGPS
      ID5 = ID4 + (3*NUNK + 3)/2

      II4 = ID4 + ID4

      LAWORK = LAWORK - ID5
      LNWORK = LNWORK - II4

*** SEE IF STORAGE EXCEEDED

      IF (LAWORK .LE. 2*NUNK) THEN
        WRITE (LUNIT,1) LAWORK, NUNK
    1   FORMAT ('0', I6, ' DOUBLE PRECISION WORDS IS TOO SMALL FOR',
     &          I5, ' PARAMETERS')
        CALL ABORT2
      ENDIF

      RETURN
      END
      SUBROUTINE ASHC (ACARD)

*** CONTROL POINT HEIGHT ADJUSTMENT RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      LOGICAL GETSSN
      LOGICAL ELFLAG, DFFLAG
      CHARACTER*80 ACARD
      CHARACTER*1 HID
      COMMON /FLAGS/  ELFLAG(MXSSN), DFFLAG(MXSSN)
      COMMON/UNITS/ LUNIT

      READ (ACARD,1) HID, ISSN
    1 FORMAT (2X, A1, I4)

      IF ( .NOT. GETSSN(ISSN, ISN) ) THEN
        CALL LINE (3)
        WRITE (LUNIT,2) ACARD
    2   FORMAT ('0NO *80* RECORD FOR --', A80/)
      ELSE
        IF (HID .EQ. 'G') THEN
          ELFLAG(ISN) = .FALSE.
        ELSE
          ELFLAG(ISN) = .TRUE.
        ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE ASTRAZ (BCARD, IUO, IOBS, B, NX, FATAL, LSN)

*** OBSERVATION EQUATIONS FOR ASTRONOMIC AZIMUTHS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      CHARACTER*80 BCARD
      CHARACTER*3 ASS
      CHARACTER*4 PVC
      CHARACTER*1 TC, AC
      LOGICAL FATAL, LSN
      LOGICAL GETSSN, GETIVF
      LOGICAL ADDCON
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL REJECT,GET82
      DIMENSION IC(LENC), C(LENC)
      DIMENSION B(*), NX(*)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON/UNITS/ LUNIT

      READ (BCARD,1) IRT, ISSN, PVC, AC, JSSN, ID, IM, ASS
 1    FORMAT (7X, I2, 1X, I4, A4, A1, 31X, I4,  9X, I3, I2, A3)
      CALL NBLANK (ASS, 3, IBLK)
      READ (ASS,3) SS
 3    FORMAT (F3.1)
      CALL NBLANK (PVC, 3, IBLK)
      READ (PVC,4) ETA
 4    FORMAT (F4.1)
      IF (BCARD(19:19) .EQ. 'W') ETA = -ETA

      READ (BCARD,5) IYR, IMO, IDY, IHR, IMN, TC
 5    FORMAT (39X, 5I2, A1)
        IF (BCARD(40:41) .EQ. '  ') IYR = 84
        IF (BCARD(42:43) .EQ. '  ') IMO = 1
        IF (BCARD(44:45) .EQ. '  ') IDY = 1
        IF (BCARD(46:47) .EQ. '  ') IHR = 0
        IF (BCARD(48:49) .EQ. '  ') IMN = 0
        IF (BCARD(50:50) .EQ. ' ') TC = 'Z'
      CALL GETYR(BCARD,IYR)
*     IF (IYR .LT. 14) THEN
*       IYR = IYR + 2000
*     ELSE
*       IYR = IYR + 1900
*     ENDIF
      CALL TOMNT (IYR, IMO, IDY, IHR, IMN, TC, ITIME)
      IF ( BCARD( 6: 6) .EQ. 'R'  .OR.  BCARD( 6: 6) .EQ. 'O'  .OR.
     &     BCARD( 6: 6) .EQ. 'F' ) THEN
        NREJ = NREJ + 1
        REJECT = .TRUE.
      ELSE
        REJECT = .FALSE.
      ENDIF

      IF (GET82(ISSN,I)) return
      IF (GET82(JSSN,I)) RETURN

*** CHECK IF AZIMUTH ORIGIN IS NORTH OR SOUTH

        IF (BCARD(72:72) .EQ. 'S') THEN
          IF (ID .GE. 0  .AND.  ID .LE. 179) THEN
            ID = ID + 180
          ELSEIF (ID .GE. 180  .AND.  ID .LE. 359) THEN
            ID = ID-180
          ELSE
            WRITE (LUNIT,666) ID
 666        FORMAT ('0ILLEGAL AZIMUTH DEGREES =', I4)
            CALL ABORT2
          ENDIF
        ENDIF

        IF ( .NOT. GETSSN(ISSN, ISN) ) THEN
          CALL LINE (3)
          WRITE (LUNIT,2) BCARD
 2        FORMAT ('0NO *80* RECORD FOR--', A80/)
        ELSEIF ( .NOT. GETSSN(JSSN, JSN) ) THEN
          CALL LINE (3)
          WRITE (LUNIT,2) BCARD
        ELSE

*** RETRIEVE THE STD DEV

          CALL STDDEV (BCARD, IRT, SD, REJECT)

*** KIND = 8  ASTRONOMIC AZIMUTH

          KIND = 8
          IOBS = IOBS + 1
          NOBS = NOBS + 1
          NAZ = NAZ + 1
          LSN = .TRUE.
          IAUX = 0
          OBSB = (ID + IM/60.D0 + SS/3600.D0)/RAD

*** CHANGE TO ASTRO AZ IF OBS IS LAPLACE AZIMUTH

          IF (AC .NE. 'A'  .AND.  BCARD(15:19) .NE. '     ') THEN
            CALL GETGLA (GLA, ISN, B)
            OBSB = OBSB + DTAN(GLA)*( ETA/(3600.D0*RAD) )
          ENDIF
          IF (OBSB .LT. 0.D0) OBSB = OBSB + PI + PI
          IF (OBSB .GE. PI+PI) OBSB = OBSB - PI - PI
          IF ( .NOT. GETIVF(60, ITIME, IVF) ) IVF = 0
          IGRT = 0
          CALL FORMIC (KIND, ISN, JSN, IDUMMY, IC, LENG, IGRT)
          CALL FORMC (KIND, C, B, ISN, JSN, IDUMMY, IGRT)
          CALL COMPOB (KIND, OBS0, B, OBSB, ISN, JSN, IDUMMY, IGRT)
          CMO = OBS0 - OBSB
          IF (CMO .GT. PI) THEN
            CMO = CMO - PI - PI
          ELSEIF (CMO .LT. -PI) THEN
            CMO = CMO + PI + PI
          ENDIF
          IF (IMODE .EQ. 0) CMO = 0.D0
          CALL BIGV (KIND, ISN, JSN, IOBS, IVF, CMO, SD, FATAL)

*** UPDATE CONNECTIVITY

          IF ( .NOT. ADDCON(IC, LENG, NX) ) THEN
            WRITE (LUNIT,667)
 667        FORMAT ('0INSUFFICIENT STORAGE FOR LAPLACE AZIMUTHS'/)
            CALL ABORT2
          ENDIF

          WRITE (IUO) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD, IOBS,
     &                IVF, IAUX, IGRT

        ENDIF

      RETURN
      END
      SUBROUTINE AVERT (A, SCR, N)

*** INVERT TRIANGULAR CHOLESKY FACTOR IN PLACE

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION A(*), SCR(*)

*** PUT CHOLESKY RECIPROCALS ON DIAGONAL

      DO 1 I = 1, N
        IDEX = INX(I, N)
        A(IDEX) = 1.D0/A(IDEX)
    1 CONTINUE

*** COMPUTE GAUSSIAN FACTOR

      DO 2 I = 1, N
        DIAG = A( INX(I, N) )
        DO 2 J = I, N
          INDEX = INX(I, N) + J - I
          A(INDEX) = DIAG*A(INDEX)
    2 CONTINUE

*** PROCEED BY ROWS -- BOTTOM TO TOP

      DO 5 I = N-1, 1, -1
        I1 = I + 1
        IDEX = INX(I, N)
        DO 3 J = I1, N
          SCR(J) = -A(INX(I, N) + J - I)
    3   CONTINUE
        DO 4 J = I1, N
          INDEX = INX(I, N) + J - I
          A(INDEX) = 0.D0
          DO 4 K = I1, N
            IF (K .LE. J) THEN
              JNDEX = INX(K, N) + J - K
            ELSE
              JNDEX = INX(J, N) + K - J
            ENDIF
            A(INDEX) = A(INDEX) + SCR(K)*A(JNDEX)
    4   CONTINUE
        DO 5 J = I1, N
          A(IDEX) = A(IDEX) + SCR(J)*A(INX(I, N) + J - I)
    5   CONTINUE

      RETURN
      END
      SUBROUTINE BIGL (KIND,ISN,JSN,IOBS,IVF,CMO,SD,FATAL)

*** LIST LARGE RESIDUALS (MISCLOSURES)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXVF = 40 )
      LOGICAL FATAL
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON /VFTB2/  VFS(MXVF), VTV(MXVF), VFRN(MXVF), VFSS(MXVF)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON/UNITS/ LUNIT

      VSD = CMO/SD
      IF (NVFTOT .GT. 0) THEN
        IF (IVF .LT. 0  .OR.  IVF .GT. NVFTOT) THEN
          WRITE (LUNIT,3) IVF,NVFTOT
 3        FORMAT ('0ILLEGAL IVF=',I10,'  FOR N=',I10,' IN BIGL')
          CALL ABORT2
        ELSEIF (IVF .NE. 0) THEN
          VSD = VSD/DSQRT(VFS(IVF))
        ENDIF
      ENDIF
      IF (DABS(VSD) .GT. VP) THEN
        CALL LINE (1)
        WRITE (LUNIT,1) IOBS,VSD
 1      FORMAT (' ','   OBS#=',I5,F70.1,' *** LARGE MISCLOSURE')
        IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
      ENDIF

      RETURN
      END
      SUBROUTINE BIGV (KIND,ISN,JSN,IOBS,IVF,CMO,SD,FATAL)

*** LIST LARGE RESIDUALS (MISCLOSURES)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXVF = 40 )
      LOGICAL FATAL, LBV, LOCSSN
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /ECHO/   VSD, LBV
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON /VFTB2/  VFS(MXVF), VTV(MXVF), VFRN(MXVF), VFSS(MXVF)
      COMMON/UNITS/ LUNIT

      LBV = .FALSE.
      VSD = CMO/SD
      IF (NVFTOT .GT. 0) THEN
        IF (IVF .LT. 0  .OR.  IVF .GT. NVFTOT) THEN
          WRITE (LUNIT,3) IVF,NVFTOT
 3        FORMAT ('0ILLEGAL IVF=',I10,'  FOR N=',I10,' FOR OBS#',
     &             I10,'  IN BIGV')
          CALL ABORT2
        ELSEIF (IVF .NE. 0) THEN
          VSD = VSD/DSQRT(VFS(IVF))
        ENDIF
      ENDIF
      IF (KIND .GE. 7  .AND.  (KIND .LE. 17  .OR.  KIND .GE. 21) ) THEN
        IF (DABS(VSD) .GT. VP) LBV = .TRUE.
        IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
      ELSEIF (KIND .LE. 6) THEN
        IF (DABS(VSD) .GT. VP) THEN
          CALL LINE (1)
          LBV = .TRUE.
          IF ( .NOT. LOCSSN(ISN,ISSN)) THEN
            WRITE (LUNIT,666) ISN
 666        FORMAT ('0SSN TABLE ERROR IN BIGV--',I5)
            CALL ABORT2
          ENDIF
          IF ( .NOT. LOCSSN(JSN,JSSN)) THEN
            WRITE (LUNIT,666) JSN
            CALL ABORT2
          ENDIF
          WRITE (LUNIT,1) IOBS,ISSN,JSSN,VSD
 1        FORMAT (' ','   OBSERVATION #',I5,' BETWEEN STATIONS',I5,
     &            ' AND',I5,F30.1,' *** LARGE MISCLOSURE')
          IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
        ENDIF
      ELSEIF (KIND .GE. 18  .AND.  KIND .LE. 20) THEN
        IF (DABS(VSD) .GT. VP) THEN
          CALL LINE (1)
          LBV = .TRUE.
          IF ( .NOT. LOCSSN(ISN,ISSN)) THEN
            WRITE (LUNIT,666) ISN
            CALL ABORT2
          ENDIF
          WRITE (LUNIT,4) IOBS,ISSN,VSD
 4        FORMAT (' ','   OBSERVATION #',I5,' AT STATIONS',I5,
     &            F30.1,' *** LARGE MISCLOSURE')
          IF (DABS(VSD) .GT. VM) FATAL = .TRUE.
        ENDIF
      ELSE
        STOP 'PROGRAMMER ERROR - ILLEGAL KIND IN BIGV'
      ENDIF

      RETURN
      END
C---------------------------------------------------------------------------------
      SUBROUTINE BUILDR (R,ISN,B)

*** RETURN THE GEODETIC HORIZON MATRIX

*** RAPP -- GEOMETRIC GEODESY, VOL II, P.115, EQ (6)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*),R(3,3)

      CALL GETGLA (GLAT,ISN,B)
      CALL GETGLO (GLON,ISN,B)
      SLAT = DSIN(GLAT)
      CLAT = DCOS(GLAT)
      SLON = DSIN(GLON)
      CLON = DCOS(GLON)

      R(1,1) = -SLAT*CLON
      R(1,2) = -SLAT*SLON
      R(1,3) = CLAT

      R(2,1) = -SLON
      R(2,2) = CLON
      R(2,3) = 0.D0

      R(3,1) = CLAT*CLON
      R(3,2) = CLAT*SLON
      R(3,3) = SLAT

      RETURN
      END
C-----------------------------------------------------------------------------------
      SUBROUTINE CMPNT (NX)

***   GET COMPONENT LIST FOR UNKNOWNS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, MAXZ =  8000, MXPRM = 40 )
      PARAMETER ( MXALL = MXPRM + 15 + MAXZ + MXSSN*3)
      INTEGER BWIDTH, NWIDTH
      DIMENSION ICMPT(MXALL), NTCMP(MXALL)
      DIMENSION NX(*)
      LOGICAL L2HLF, LEHT

      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /OBSUM/  NSUM(MXSSN,13), ICMPL(MXALL)
      COMMON/UNITS/ LUNIT

*** FOR DEBUGGING, PRINT THE NX ARRAY

*     CALL PRNX (NX, NUNK)

*** ZERO COMPONTENT TABLE

      DO 50 I = 1, MXALL
         ICMPT(I) = 0
         NTCMP(I) = 0
   50 CONTINUE

*** NUMBER THE COMPONENTS

      N1 = NUNK + 1
      NCMP = 0
      BWIDTH = 0
      DO 100 I = 1, NUNK-1
         IF (I .LE. BWIDTH) THEN
            NWIDTH = NX(I+1) - NX(I) + I - 1
            IF (NWIDTH .GT. BWIDTH) BWIDTH = NWIDTH
         ELSE
            NCMP = NCMP + 1
            BWIDTH = NX(I+1) - NX(I) + I - 1
         ENDIF
         ICMPL( NX(N1 + I) ) = NCMP
         ICMPT(NCMP) = ICMPT(NCMP) + 1
  100 CONTINUE
      IF (BWIDTH .LT. NUNK) NCMP = NCMP + 1
      ICMPL( NX(N1 + NUNK) ) = NCMP
      ICMPT(NCMP) = ICMPT(NCMP) + 1

*** IF A COMPONENT INCLUDES ONLY ONE UNKNOWN, ZERO THE ELEMENT
*** IN THE UNKNOWN COMPONENT LIST

      DO 200 I = 1, NUNK
         IF ( ICMPT( ICMPL(I) ) .EQ. 1) THEN
            ICMPL(I) = 0
         ENDIF
  200 CONTINUE

*** GET THE NUMBER OF NON-TRIVIAL COMPONENTS AND
*** SET UP THE RENUMBERING ARRAY FOR THE NON-TRIVIAL COMPONENTS

      NNTCMP = 0
      DO 300 I = 1, NCMP
         IF (ICMPT(I) .GT. 1) THEN
             NNTCMP = NNTCMP + 1
             NTCMP(I) = NNTCMP
         ENDIF
  300 CONTINUE

*** RENUMBER THE NON-TRIVIAL COMPONENTS

      DO 400 I = 1, NUNK
         IF ( ICMPL(I) .NE. 0) THEN
            ICMPL(I) = NTCMP( ICMPL(I) )
         ENDIF
  400 CONTINUE

      CALL LINE (2)
      WRITE (LUNIT,500) NNTCMP
  500 FORMAT (/,' THE NUMBER OF NON-TRIVIAL COMPONENTS IS', I3, '.')

      RETURN
      END
      SUBROUTINE PRNX (NX, NUNK)

***   PRINT NX ARRAY

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION NX(*)
      COMMON/UNITS/ LUNIT

      N1 = NUNK + 1
      N2 = N1 + NUNK
      N3 = N2 + NUNK

      WRITE (LUNIT,65) NUNK
   65 FORMAT(' NUMBER OF UNKNOWNS=', I4, /, ' --- PROFILE ---')
      WRITE (LUNIT,66) (NX(I), I = 1, N1)
   66 FORMAT (8I10)
      WRITE (LUNIT,67)
   67 FORMAT (' --- INTERNAL TO EXTERNAL TABLE ---')
      WRITE (LUNIT,66) (NX(I), I = N1 + 1, N2)
      WRITE (LUNIT,68)
   68 FORMAT (' --- EXTERNAL TO INTERNAL TABLE ---')
      WRITE (LUNIT,66) (NX(I), I = N2 + 1, N3)

      RETURN
      END
      SUBROUTINE COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)

*** COMPUTE THE OBSERVATION FROM THE PARAMETERS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      LOGICAL L2HLF, LEHT
      DIMENSION B(*)
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /CONST/  PI, PI2, RAD
      COMMON /GRTTB1/ GLAT0, GLON0, R0(3,3)
      COMMON/UNITS/ LUNIT

*** KIND 0 IS AN AUXILIARY PARAMETER CONSTRAINT

      IF (KIND .EQ. 0) THEN
        IAUX = ISN
        CALL GETAUX (AUX,IAUX,B)
        OBS0 = AUX

*** KIND 1 IS LATITUDE CONSTRAINT (OBSB = GEOD. LAT)

      ELSEIF (KIND .EQ. 1) THEN
        CALL GETGLA (GLA0,ISN,B)
        CALL GETMSL (GMSL0,ISN,B)
        CALL GETGH (GHT0,ISN,B)
        EHT0 = GMSL0+GHT0
        CALL RADII (ISN,RMER,RPV,B)
        CMO = (GLA0-OBSB)*(RMER+EHT0)
        OBS0 = CMO

*** KIND 2 IS LONGITUDE CONSTRAINT (OBSB = GEOD. LONG.)

      ELSEIF (KIND .EQ. 2) THEN
        CALL GETGLA (GLA0,ISN,B)
        CALL GETGLO (GLO0,ISN,B)
        CALL GETMSL (GMSL0,ISN,B)
        CALL GETGH (GHT0,ISN,B)
        EHT0 = GMSL0+GHT0
        CALL RADII (ISN,RMER,RPV,B)
        CMO = (GLO0-OBSB)*(RPV+EHT0)*DCOS(GLA0)
        OBS0 = CMO

*** KIND 3 IS HEIGHT CONSTRAINT

      ELSEIF (KIND .EQ. 3) THEN
        IF (N2HLF .GT. 0) THEN
          CALL GETEHT (EHT0,ISN,B)
        ELSE
          CALL GETMSL (GMSL0,ISN,B)
          CALL GETGH (GHT0,ISN,B)
          EHT0 = GMSL0+GHT0
        ENDIF
        OBS0 = EHT0

*** KIND 4 IS DIFFERENTIAL X (GPS)

      ELSEIF (KIND .EQ. 4) THEN
        IF (L2HLF) THEN
          CALL GETEHX (XI,ISN,B)
          CALL GETEHX (XJ,JSN,B)
        ELSE
          CALL GETECX (XI,ISN,B)
          CALL GETECX (XJ,JSN,B)
        ENDIF
        DELX = XJ-XI
        OBS0 = DELX
        IF (IAUX .GT. 0) THEN
          CALL GETAUX (VAL,IAUX,B)
          OBS0 = OBS0-OBS0*VAL
        ENDIF
        IF (IGRT .GT. 0) THEN
          IF (L2HLF) THEN
            CALL GETEHY (YI,ISN,B)
            CALL GETEHY (YJ,JSN,B)
            CALL GETEHZ (ZI,ISN,B)
            CALL GETEHZ (ZJ,JSN,B)
          ELSE
            CALL GETECY (YI,ISN,B)
            CALL GETECY (YJ,JSN,B)
            CALL GETECZ (ZI,ISN,B)
            CALL GETECZ (ZJ,JSN,B)
          ENDIF
          DELY = YJ-YI
          IF (IAUX .GT. 0) DELY = DELY - DELY * VAL
          DELZ = ZJ-ZI
          IF (IAUX .GT. 0) DELZ = DELZ - DELZ * VAL
          CALL GETRTG (VALX,VALY,VALZ,IGRT)
          OBS0 = OBS0
     &           + ( R0(1,3)*VALX + R0(2,3)*VALY + R0(3,3)*VALZ )*DELY
     &           - ( R0(1,2)*VALX + R0(2,2)*VALY + R0(3,2)*VALZ )*DELZ
        ENDIF

*** KIND 5 IS DIFFERENTIAL Y (GPS)

      ELSEIF (KIND .EQ. 5) THEN
        IF (L2HLF) THEN
          CALL GETEHY (YI,ISN,B)
          CALL GETEHY (YJ,JSN,B)
        ELSE
          CALL GETECY (YI,ISN,B)
          CALL GETECY (YJ,JSN,B)
        ENDIF
        OBS0 = YJ-YI
        IF (IAUX .GT. 0) THEN
          CALL GETAUX (VAL,IAUX,B)
          OBS0 = OBS0-OBS0*VAL
        ENDIF
        IF (IGRT .GT. 0) THEN
          IF (L2HLF) THEN
            CALL GETEHX (XI,ISN,B)
            CALL GETEHX (XJ,JSN,B)
            CALL GETEHZ (ZI,ISN,B)
            CALL GETEHZ (ZJ,JSN,B)
          ELSE
            CALL GETECX (XI,ISN,B)
            CALL GETECX (XJ,JSN,B)
            CALL GETECZ (ZI,ISN,B)
            CALL GETECZ (ZJ,JSN,B)
          ENDIF
          DELX = XJ-XI
          IF (IAUX .GT. 0) DELX = DELX - DELX * VAL
          DELZ = ZJ-ZI
          IF (IAUX .GT. 0) DELZ = DELZ - DELZ * VAL
          CALL GETRTG (VALX,VALY,VALZ,IGRT)
          OBS0 = OBS0
     &           - ( R0(1,3)*VALX + R0(2,3)*VALY + R0(3,3)*VALZ )*DELX
     &           + ( R0(1,1)*VALX + R0(2,1)*VALY + R0(3,1)*VALZ )*DELZ
        ENDIF

*** KIND 6 IS DIFFERENTIAL Z (GPS)

      ELSEIF (KIND .EQ. 6) THEN
        IF (L2HLF) THEN
          CALL GETEHZ (ZI,ISN,B)
          CALL GETEHZ (ZJ,JSN,B)
        ELSE
          CALL GETECZ (ZI,ISN,B)
          CALL GETECZ (ZJ,JSN,B)
        ENDIF
        OBS0 = ZJ-ZI
        IF (IAUX .GT. 0) THEN
          CALL GETAUX (VAL,IAUX,B)
          OBS0 = OBS0-OBS0*VAL
        ENDIF
        IF (IGRT .GT. 0) THEN
          IF (L2HLF) THEN
            CALL GETEHX (XI,ISN,B)
            CALL GETEHX (XJ,JSN,B)
            CALL GETEHY (YI,ISN,B)
            CALL GETEHY (YJ,JSN,B)
          ELSE
            CALL GETECX (XI,ISN,B)
            CALL GETECX (XJ,JSN,B)
            CALL GETECY (YI,ISN,B)
            CALL GETECY (YJ,JSN,B)
          ENDIF
          DELX = XJ-XI
          IF (IAUX .GT. 0) DELX = DELX - DELX * VAL
          DELY = YJ-YI
          IF (IAUX .GT. 0) DELY = DELY - DELY * VAL
          CALL GETRTG (VALX,VALY,VALZ,IGRT)
          OBS0 = OBS0
     &           + ( R0(1,2)*VALX + R0(2,2)*VALY + R0(3,2)*VALZ )*DELX
     &           - ( R0(1,1)*VALX + R0(2,1)*VALY + R0(3,1)*VALZ )*DELY
        ENDIF

*** KIND 7 IS MARK TO MARK DISTANCE

      ELSEIF (KIND .EQ. 7) THEN
        CALL GETLAH (ISN,JSN,B,P1,Q1,R1,T1,P2,Q2,R2,T2,S)
        OBS0 = S
        IF (IAUX .GT. 0) THEN
          CALL GETAUX (VAL,IAUX,B)
          OBS0 = OBS0-OBS0*VAL
        ENDIF

*** KIND 8 IS AN ASTRONOMIC AZIMUTH

      ELSEIF (KIND .EQ. 8) THEN
        CALL GETLAH (ISN,JSN,B,P1,Q1,R1,T1,P2,Q2,R2,T2,S)
        OBS0 = DATAN2(Q1,P1)
        IF (OBS0 .LT. 0.D0) OBS0 = OBS0 + PI + PI

*** KIND 9 IS A ZENITH DISTANCE

      ELSEIF (KIND .EQ. 9) THEN
        CALL GETLAH (ISN,JSN,B,P1,Q1,R1,T1,P2,Q2,R2,T2,S)
        OBS0 = DATAN2(T1,R1)
        OBS0 = PI2 - OBS0
        IF (IAUX .GT. 0) THEN
          CALL GETAUX (VAL,IAUX,B)
          OBS0 = OBS0 - VAL*R1
        ENDIF

*** KIND 10 IS A HORIZONTAL ANGLE

      ELSEIF (KIND .EQ. 10) THEN
        KSN = IAUX
        CALL GETLAH (ISN,JSN,B,P1,Q1,R1,T1,P2,Q2,R2,T2,S)
        DIRL = DATAN2(Q1,P1)
        CALL GETLAH (ISN,KSN,B,P1,Q1,R1,T1,P2,Q2,R2,T2,S)
        DIRR = DATAN2(Q1,P1)
        OBS0 = DIRR-DIRL
 18     IF (OBS0 .GE. PI+PI) THEN
          OBS0 = OBS0 - PI - PI
          GO TO 18
        ENDIF
 19     IF (OBS0 .LT. 0.D0) THEN
          OBS0 = OBS0 + PI + PI
          GO TO 19
        ENDIF

*** KIND 11 IS A HORIZONTAL DIRECTION

      ELSEIF (KIND .EQ. 11) THEN
        CALL GETLAH (ISN,JSN,B,P1,Q1,R1,T1,P2,Q2,R2,T2,S)
        CALL GETROT (ROT,IAUX,B)
        OBS0 = DATAN2(Q1,P1)-ROT
 21     IF (OBS0 .GE. PI+PI) THEN
          OBS0 = OBS0-PI-PI
          GO TO 21
        ENDIF
 23     IF (OBS0 .LT. 0.D0) THEN
          OBS0 = OBS0+PI+PI
          GO TO 23
        ENDIF

*** KIND 12 IS A CONSTRAINED AZIMUTH

      ELSEIF (KIND .EQ. 12) THEN
        CALL GETLAH (ISN,JSN,B,P1,Q1,R1,T1,P2,Q2,R2,T2,S)
        OBS0 = DATAN2(Q1,P1)
        IF (OBS0 .LT. 0.D0) OBS0 = OBS0+PI+PI

*** KIND 13 IS CONSTRAINED MARK TO MARK DISTANCE

      ELSEIF (KIND .EQ. 13) THEN
        CALL GETLAH (ISN,JSN,B,P1,Q1,R1,T1,P2,Q2,R2,T2,S)
        OBS0 = S

*** KIND 14 IS A CONSTRAINED ZENITH DISTANCE

      ELSEIF (KIND .EQ. 14) THEN
        CALL GETLAH (ISN,JSN,B,P1,Q1,R1,T1,P2,Q2,R2,T2,S)
        OBS0 = DATAN2(T1,R1)
        OBS0 = PI2-OBS0

*** KIND 15 IS A CONSTRAINED GEOID HEIGHT DIFFERENCE

      ELSEIF (KIND .EQ. 15) THEN
        CALL GETGH (DGH1,ISN,B)
        CALL GETGH (DGH2,JSN,B)
        OBS0 = DGH2-DGH1

*** KIND 16 IS A CONSTRAINED ORTHOMETRIC HEIGHT DIFFERENCE

      ELSEIF (KIND .EQ. 16) THEN
        CALL GETMSL (GMSL1,ISN,B)
        CALL GETMSL (GMSL2,JSN,B)
        OBS0 = GMSL2-GMSL1

*** KIND 17 IS A CONSTRAINED ELLIPSOIDAL HEIGHT DIFFERENCE

      ELSEIF (KIND .EQ. 17) THEN
        CALL GETMSL (GMSL1,ISN,B)
        CALL GETGH (DGH1,ISN,B)
        CALL GETMSL (GMSL2,JSN,B)
        CALL GETGH (DGH2,JSN,B)
        OBS0 = (GMSL2+DGH2)-(GMSL1+DGH1)

*** KIND 18 IS E.C.F. X (DOPPLER)

      ELSEIF (KIND .EQ. 18) THEN
        IF (L2HLF) THEN
          CALL GETEHX (XI,ISN,B)
        ELSE
          CALL GETECX (XI,ISN,B)
        ENDIF
        DELX = +XI
        OBS0 = DELX
        IF (IAUX .GT. 0) THEN
          CALL GETAUX (VAL,IAUX,B)
          OBS0 = OBS0-OBS0*VAL
        ENDIF
        IF (IGRT .GT. 0) THEN
          IF (L2HLF) THEN
            CALL GETEHY (YI,ISN,B)
            CALL GETEHZ (ZI,ISN,B)
          ELSE
            CALL GETECY (YI,ISN,B)
            CALL GETECZ (ZI,ISN,B)
          ENDIF
          DELY = +YI
          IF (IAUX .GT. 0) DELY = DELY - DELY * VAL
          DELZ = +ZI
          IF (IAUX .GT. 0) DELZ = DELZ - DELZ * VAL
          CALL GETRTG (VALX,VALY,VALZ,IGRT)
          OBS0 = OBS0
     &           + ( R0(1,3)*VALX + R0(2,3)*VALY + R0(3,3)*VALZ )*DELY
     &           - ( R0(1,2)*VALX + R0(2,2)*VALY + R0(3,2)*VALZ )*DELZ
        ENDIF

*** KIND 19 IS E.C.F. Y (DOPPLER)

      ELSEIF (KIND .EQ. 19) THEN
        IF (L2HLF) THEN
          CALL GETEHY (YI,ISN,B)
        ELSE
          CALL GETECY (YI,ISN,B)
        ENDIF
        OBS0 = +YI
        IF (IAUX .GT. 0) THEN
          CALL GETAUX (VAL,IAUX,B)
          OBS0 = OBS0-OBS0*VAL
        ENDIF
        IF (IGRT .GT. 0) THEN
          IF (L2HLF) THEN
            CALL GETEHX (XI,ISN,B)
            CALL GETEHZ (ZI,ISN,B)
          ELSE
            CALL GETECX (XI,ISN,B)
            CALL GETECZ (ZI,ISN,B)
          ENDIF
          DELX = +XI
          IF (IAUX .GT. 0) DELX = DELX - DELX * VAL
          DELZ = +ZI
          IF (IAUX .GT. 0) DELZ = DELZ - DELZ * VAL
          CALL GETRTG (VALX,VALY,VALZ,IGRT)
          OBS0 = OBS0
     &           - ( R0(1,3)*VALX + R0(2,3)*VALY + R0(3,3)*VALZ )*DELX
     &           + ( R0(1,1)*VALX + R0(2,1)*VALY + R0(3,1)*VALZ )*DELZ
        ENDIF

*** KIND 20 IS E.C.F. Z (DOPPLER)

      ELSEIF (KIND .EQ. 20) THEN
        IF (L2HLF) THEN
          CALL GETEHZ (ZI,ISN,B)
        ELSE
          CALL GETECZ (ZI,ISN,B)
        ENDIF
        OBS0 = +ZI
        IF (IAUX .GT. 0) THEN
          CALL GETAUX (VAL,IAUX,B)
          OBS0 = OBS0-OBS0*VAL
        ENDIF
        IF (IGRT .GT. 0) THEN
          IF (L2HLF) THEN
            CALL GETEHX (XI,ISN,B)
            CALL GETEHY (YI,ISN,B)
          ELSE
            CALL GETECX (XI,ISN,B)
            CALL GETECY (YI,ISN,B)
          ENDIF
          DELX = +XI
          IF (IAUX .GT. 0) DELX = DELX - DELX * VAL
          DELY = +YI
          IF (IAUX .GT. 0) DELY = DELY - DELY * VAL
          CALL GETRTG (VALX,VALY,VALZ,IGRT)
          OBS0 = OBS0
     &           + ( R0(1,2)*VALX + R0(2,2)*VALY + R0(3,2)*VALZ )*DELX
     &           - ( R0(1,1)*VALX + R0(2,1)*VALY + R0(3,1)*VALZ )*DELY
        ENDIF

*** KIND 21 IS A CONSTRAINED NORTH COORD DIFFERENCE

      ELSEIF(KIND.EQ.21) THEN
        CALL GETLAH(ISN,JSN,B,P1,Q1,R1,T1,P2,Q2,R2,T2,S)
        OBS0 = P1

*** KIND 22 IS A CONSTRAINED EAST COORD DIFFERENCE

      ELSEIF(KIND.EQ.22) THEN
        CALL GETLAH(ISN,JSN,B,P1,Q1,R1,T1,P2,Q2,R2,T2,S)
        OBS0 = Q1

*** ILLEGAL KIND

      ELSE
        WRITE (LUNIT,1) KIND
    1   FORMAT ('0ILLEGAL KIND IN COMPOB =',I5)
        CALL ABORT2
      ENDIF

      RETURN
      END
      SUBROUTINE COMRHS (G, N)

*** SINGLE SUBSTITUTION FOR DECORRELATION FROM A CHOL. FACTOR

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION G(N,*)

      N1 = N + 1
      DO 1 I = 1, N
        I1 = I - 1
        TEMP = 0.D0
        IF (I1 .GT. 0) THEN
          DO 2 K = 1, I1
            TEMP = TEMP + G(K,I)*G(K,N1)
    2     CONTINUE
        ENDIF
        G(I,N1) = (G(I,N1)-TEMP)/G(I,I)
    1 CONTINUE

      RETURN
      END
      SUBROUTINE COMSLA ( A, NX, IC, C3, LENG, COVLA)

*** COMPUTATION OF COVARIANCE MATRIX FOR ADJUSTED DOPPLER OBSERVATIONS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      DIMENSION A(*), NX(*)
      DIMENSION IC(LENC), CI(LENC), CJ(LENC)
      DIMENSION C3(3,LENC)
      DIMENSION COVLA(3,3)
      LOGICAL PROP, PROPCV
      COMMON/UNITS/ LUNIT

      DO 110 I = 1, 3
        DO 105 K = 1, LENG
          CI(K) = C3(I,K)
  105   CONTINUE
        DO 120 J = I, 3
          IF (J .EQ. I) THEN
            IF (PROP(CI,IC,LENG,VAR,A,NX,IFLAG)) THEN
              COVLA(I,I) = VAR
            ELSEIF (IFLAG .EQ. 1) THEN
              WRITE (LUNIT,660)
 660          FORMAT ('0SYSTEM NOT INVERTED IN COMSLA'/)
              CALL ABORT2
            ELSE
              WRITE (LUNIT,661)
 661          FORMAT ('0ALL COVARIANCE ELEMENTS NOT WITHIN PROFILE',
     &                '--COMSLA')
              CALL ABORT2
            ENDIF
          ELSE
            DO 135 K = 1, LENG
              CJ(K) = C3(J,K)
  135       CONTINUE
            IF (PROPCV(CI,IC,LENG,CJ,IC,LENG,COV,A,NX,IFLAG)) THEN
              COVLA(I,J) = COV
              COVLA(J,I) = COV
            ELSEIF (IFLAG .EQ. 1) THEN
              WRITE (LUNIT,660)
              CALL ABORT2
            ELSE
              WRITE (LUNIT,661)

              CALL ABORT2
            ENDIF
          ENDIF
  120   CONTINUE
  110 CONTINUE

      RETURN
      END
      SUBROUTINE CONEC (I,J,IC,LENG)

*** FILL CONNECTION MATRIX FOR TWO STATIONS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      DIMENSION IC(LENC)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT

      IF (IDIM .EQ. 1) THEN
        LENG = 2
        IC(1) = IUNSTA(I,3)
        IC(2) = IUNSTA(J,3)
      ELSEIF (IDIM .EQ. 2) THEN
        LENG = 4
        IC(1) = IUNSTA(I,1)
        IC(2) = IUNSTA(I,2)
        IC(3) = IUNSTA(J,1)
        IC(4) = IUNSTA(J,2)
      ELSE
        LENG = 6
        IC(1) = IUNSTA(I,1)
        IC(2) = IUNSTA(I,2)
        IC(3) = IUNSTA(I,3)
        IC(4) = IUNSTA(J,1)
        IC(5) = IUNSTA(J,2)
        IC(6) = IUNSTA(J,3)
      ENDIF

      RETURN
      END
      SUBROUTINE CONTR (KIND,I,J,K)

*** KEEP COUNT OF EACH OBSERVATIONAL TYPE

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, MAXZ =  8000, MXPRM = 40 )
      PARAMETER ( MXALL = MXPRM + 15 + MAXZ + MXSSN*3)
      COMMON /OBSUM/  NSUM(MXSSN,13), ICMPL(MXALL)

*** COUNT ALL THE AZIMUTHS FROM(*,1) AND TO(*,2)

      IF (KIND .EQ. 8  .OR.  KIND .EQ. 12) THEN
        NSUM(I,1) = NSUM(I,1)+1
        NSUM(J,2) = NSUM(J,2)+1

*** COUNT ALL THE DISTANCES FROM(*,3) AND TO(*,4)

      ELSEIF (KIND .EQ. 7  .OR.  KIND .EQ. 13) THEN
        NSUM(I,3) = NSUM(I,3)+1
        NSUM(J,4) = NSUM(J,4)+1

*** COUNT ALL THE ZENITH DISTANCES FROM(*,5) AND TO(*,6)

      ELSEIF (KIND .EQ. 9  .OR.  KIND .EQ. 14) THEN
        NSUM(I,5) = NSUM(I,5)+1
        NSUM(J,6) = NSUM(J,6)+1

*** COUNT ALL THE GPS FROM(*,7) AND TO(*,8)

      ELSEIF (KIND .EQ. 4) THEN
        NSUM(I,7) = NSUM(I,7)+1
        NSUM(J,8) = NSUM(J,8)+1

*** COUNT THE DIRECTIONS FROM(*,9) AND TO(*,10)

      ELSEIF (KIND .EQ. 11) THEN
        NSUM(I,9) = NSUM(I,9)+1
        NSUM(J,10) = NSUM(J,10)+1

*** COUNT THE ANGLES FROM(*,11) AND TO(*,12)

      ELSEIF (KIND .EQ. 10) THEN
        NSUM(I,11) = NSUM(I,11)+1
        NSUM(J,12) = NSUM(J,12)+1
        NSUM(K,12) = NSUM(K,12)+1

*** COUNT ALL THE DOPPLER AT(*,13)

      ELSEIF (KIND .EQ. 18) THEN
        NSUM(I,13) = NSUM(I,13)+1
      ENDIF

      RETURN
      END
C-----------------------------------------------------------------------------------------
      LOGICAL FUNCTION CONVRG (A, NX, B, SHIFTS, SIGUWT, ITER)

*** UPDATE THE UNKNOWNS AND TEST FOR CONVERGENCE

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      CHARACTER*30 NAME1,NAMES
      LOGICAL GETA
      LOGICAL ELFLAG,DFFLAG
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL L2HLF, LEHT
      DIMENSION A(*),NX(*),B(*),SHIFTS(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /FLAGS/  ELFLAG(MXSSN), DFFLAG(MXSSN)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT
      COMMON /CONST/  PI, PI2, RAD
      SAVE ISNX,GLAT,GLON,EHT,RX,RY
      SAVE SHFMX,ISHFMX,IISN

**v6.3

      LOGICAL    overwrite_const_coord_in_newbb
      character  CC_records*80,hemi*1
      COMMON/MM6/overwrite_const_coord_in_newbb
      COMMON/MM6_array/CC_records(MXSSN)

**So far for v6.3

      SSQ = 0.D0
      N1 = NUNK+1
      SHFMX = 0.D0
      IISN = 0
      ISHFMX = 0

*** LOOP OVER ALL UNKNOWNS

      ISNX = 0
      DO 2 I = 1,NUNK
        IF ( .NOT. GETA(I,N1,VAL,A,NX)) THEN
          WRITE (LUNIT,1)
    1     FORMAT ('0PROGRAMMER ERROR IN CONVRG')
          CALL ABORT2
        ELSE
          CALL INVIUN (I,ISN,IUTYPE,IUCODE)
          IF (IUCODE .EQ. 0) THEN

*** GET COORDINATES AND RADII

            IF (ISN .NE. ISNX) THEN
              ISNX = ISN
              CALL GETGLA (GLAT,ISN,B)
              CALL GETGLO (GLON,ISN,B)
              IF (L2HLF  .AND.  I2HLF(ISN) .EQ. 1 ) THEN
                CALL GETEHT (EHT,ISN,B)
              ELSE
                CALL GETMSL (GMSL,ISN,B)
                CALL GETGH (GHT,ISN,B)
                EHT = GMSL+GHT
              ENDIF
              CALL RADII (ISN,RMER,RPV,B)
              RX = RMER+EHT
              RY = (RPV+EHT)*DCOS(GLAT)
            ENDIF

*** UPDATE LATITUDE

            IF ( IUTYPE .EQ. 1  .AND.
     &          (IDIM .EQ. 2  .OR.  IDIM .EQ. 3) ) THEN

              IF (LPS) THEN
                J = IUNSHF(ISN,IUTYPE)
                SHIFTS(J) = SHIFTS(J)+VAL
              ENDIF
              SSQ = SSQ+VAL*VAL
              GLAT = GLAT+VAL/RX
              CALL PUTGLA (GLAT,ISN,B)
              IF ( .NOT. DFFLAG(ISN)) CALL PUTALA (GLAT,ISN,B)

*** UPDATE LONGITUDE

            ELSEIF ( IUTYPE .EQ. 2  .AND.
     &              (IDIM .EQ. 2  .OR.  IDIM .EQ. 3) ) THEN

              IF (LPS) THEN
                J = IUNSHF(ISN,IUTYPE)
                SHIFTS(J) = SHIFTS(J)+VAL
              ENDIF
              SSQ = SSQ+VAL*VAL
              GLON = GLON+VAL/RY
              CALL PUTGLO (GLON,ISN,B)
              IF ( .NOT. DFFLAG(ISN)) CALL PUTALO (GLON,ISN,B)

*** UPDATE HEIGHT

            ELSEIF ( IUTYPE .EQ. 3  .AND.
     &              (IDIM .EQ. 1  .OR.  IDIM .EQ. 3) ) THEN

              IF (LPS) THEN
                J = IUNSHF(ISN,IUTYPE)
                SHIFTS(J) = SHIFTS(J)+VAL
              ENDIF
              SSQ = SSQ+VAL*VAL
              IF ( .NOT. ELFLAG(ISN)) THEN
                GHT = GHT+VAL
                CALL PUTGH (GHT,ISN,B)
              ELSE
                GMSL = GMSL+VAL
                CALL PUTMSL (GMSL,ISN,B)
              ENDIF

*** UPDATE HEIGHT FOR 2.5DIM ADJUSTMENTS

            ELSEIF (IUTYPE .EQ. 3  .AND.
     &              L2HLF  .AND.  I2HLF(ISN) .EQ. 1 ) THEN
              IF (LPS) THEN
                J = IUNSHF(ISN,IUTYPE)
                SHIFTS(J) = SHIFTS(J)+VAL
              ENDIF
              SSQ = SSQ+VAL*VAL
              EHT = EHT + VAL
              CALL PUTEHT (EHT,ISN,B)
            ENDIF

*** LOCATE THE MAXIMUM SHIFT

            IF (LPS) THEN
              IF (DABS(VAL) .GT. SHFMX) THEN
                SHFMX = DABS(VAL)
                ISHFMX = I
                IISN = ISN
              ENDIF
            ENDIF

*** AUXILIARY PARAMETER TYPES

          ELSEIF (IUCODE .EQ. 1) THEN
            IAUX = ISN
            CALL GETAUX (AUX,IAUX,B)
            AUX = AUX+VAL
            CALL PUTAUX (AUX,IAUX,B)

*** AUXILIARY GPS AND DOPPLER ROTATION PARAMETER TYPES

          ELSEIF (IUCODE .EQ. 2) THEN
            IGRT = ISN
            CALL GETRTG (VALX,VALY,VALZ,IGRT)
            IF (IUTYPE .EQ. 1) THEN
              VALX = VALX+VAL
            ELSEIF (IUTYPE .EQ. 2) THEN
              VALY = VALY+VAL
            ELSEIF (IUTYPE .EQ. 3) THEN
              VALZ = VALZ+VAL
            ENDIF
            CALL PUTRTG (VALX,VALY,VALZ,IGRT)

*** ROTATION PARAMETER TYPES

          ELSEIF (IUCODE .EQ. 3) THEN
            IZ = ISN
            CALL GETROT (AUX,IZ,B)
            AUX = AUX+VAL
            CALL PUTROT (AUX,IZ,B)

          ELSE
            WRITE (LUNIT,666) IUCODE
  666       FORMAT ('0ILLEGAL PARAMETER CODE IN CONVRG--',I5)
            CALL ABORT2
          ENDIF
        ENDIF
    2 CONTINUE

*** COMPUTE E.C.F. X,Y,Z

      CALL ALLXYZ (B,GLAAVE,GLOAVE)

*** TEST FOR CONVERGENCE

      IF (N2HLF .GT. 0) THEN
        VAL = DSQRT( SSQ/( NSTA*IDIM ) )
      ELSE
        VAL = DSQRT( SSQ/( NSTA*IDIM + N2HLF ) )
      ENDIF
      IF (VAL .LE. CTOL) THEN
        CONVRG = .TRUE.
      ELSE
        CONVRG = .FALSE.
      ENDIF

*** PRINT RMS SHIFT AND V'PV

      IF ( .NOT. GETA(N1,N1,VTPV,A,NX)) THEN
        WRITE (LUNIT,667) N1
  667   FORMAT ('0BASEMENT WINDOW ERROR--',I5)
        CALL ABORT2
      ENDIF
      IDOF = NOBS+NCON-NUNK
      IF (VTPV .LT. 0.0) THEN
        WRITE (LUNIT,5)
   5    FORMAT ('0*** WARNING: NEGATIVE VTPV ***')
        VARUWT = DIVIDE(DABS(VTPV),IDOF)
      ELSE
        VARUWT = DIVIDE(VTPV,IDOF)
      ENDIF
      SIGUWT = DSQRT(VARUWT)
      CALL LINE (2)
      WRITE (LUNIT,3) ITER,VAL,VTPV,IDOF,VARUWT
    3 FORMAT ('0ITERATION #',I2,' THE RMS CORRECTION IS ',F10.3,
     *        ' METERS --- VTPV=',F15.3,'  DF=',I6,
     *        '  VARIANCE=',F15.2)
      IF (LPS) THEN
        CALL LINE (2)
        CALL INVIUN (ISHFMX,IISN,ITP,ICODE)
        NAME1 = NAMES(IISN)
        IF (ITP .EQ. 1) THEN
          WRITE (LUNIT,6) NAME1,SHFMX
   6      FORMAT ('0 MAXIMUM SHIFT - STATION: ',A30,10X,
     *            'LATITUDE SHIFT= ',F10.3,' METERS')
        ELSEIF (ITP .EQ. 2) THEN
          WRITE (LUNIT,7) NAME1,SHFMX
   7      FORMAT ('0 MAXIMUM SHIFT - STATION: ',A30,10X,
     *            'LONGITUDE SHIFT= ',F10.3,' METERS')
        ELSEIF (ITP .EQ. 3) THEN
          WRITE (LUNIT,8) NAME1,SHFMX
   8      FORMAT ('0 MAXIMUM SHIFT - STATION: ',A30,10X,
     *            'VERTICAL SHIFT= ',F10.3,' METERS')
        ENDIF
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------------------------------------
***v 4.27

c      SUBROUTINE COVLGH(SN,SE,SU,ISN,JSN,B,I,G,NR,NC,LENG)
c 2/2/06 Mike Potterfield add SXYZ and RM to the argument list
c so that the latter two arrays are returned to SUBROUTINE
c RGPS for the purpose of computing the covariance matrix for
c the dNEU residuals
      SUBROUTINE COVLGH(SN,SE,SU,ISN,JSN,B,I,G,NR,NC,LENG,SXYZ,RM)
      
*** CONVERT COVARIANCE MATRIX TO LGH

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      DIMENSION SXYZ(3,3),SNEU(3,3)
      DIMENSION RM(3,3),W(3)
      DIMENSION G(NR,NC)
      SAVE GLA1, GLO1, GLA2, GLO2
      SAVE SB, CB, SL, CL
      SAVE ISNX, JSNX
      DATA ISNX /0/, JSNX /0/
      
      IF( ISN .NE. ISNX .OR. JSN .NE. JSNX) THEN
        IF ( ISN .NE. ISNX) THEN
          CALL GETGLA (GLA1, ISN, B)
          CALL GETGLO (GLO1, ISN, B)
          ISNX = ISN
        ENDIF
        IF ( JSN .NE. JSNX) THEN
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
      
*** GET COVARIANCE MATRIX
c  I is always the Z component
      NVEC = NR/3
      CALL GLOCAT(N1,N2,N3,N4,N5,NVEC,LENG)
      
      SX=G(I-2,N5)
      SY=G(I-1,N5)
      SZ=G(I,N5)
      RHOXY=G(I-2,N5+1)
      RHOYZ=G(I-1,N5+1)
      RHOXZ=G(I,N5+1)
      SXYZ(1,1)=SX*SX    
      SXYZ(1,2)=SX*SY*RHOXY    
      SXYZ(1,3)=SX*SZ*RHOXZ    
      SXYZ(2,1)=SXYZ(1,2)    
      SXYZ(2,2)=SY*SY    
      SXYZ(2,3)=SY*SZ*RHOYZ    
      SXYZ(3,1)=SXYZ(1,3)    
      SXYZ(3,2)=SXYZ(2,3)    
      SXYZ(3,3)=SZ*SZ     
      
*** PUT TOGETHER ROTATION MATRIX --- CMP

      RM(1,1)=-SB*CL
      RM(2,1)=-SL
      RM(3,1)=CB*CL
      RM(1,2)=-SB*SL
      RM(2,2)=CL 
      RM(3,2)=CB*SL
      RM(1,3)=CB
      RM(2,3)=0 
      RM(3,3)=SB
      
*** USE THE SUBROUTINE TO DO THE TRANSFORMATION

      N=3
      CALL ABAT(RM,SXYZ,SNEU,W,N,N)
      
*** CALCULATE THE RESIDUALS

      SN=DSQRT(SNEU(1,1))
      SE=DSQRT(SNEU(2,2))
      SU=DSQRT(SNEU(3,3))
      
      RETURN
      END
      
***********************
      
      SUBROUTINE DIMCON (IOBS,IUO,B)

*** ADD CONSTRAINT OBS EQS TO SATISFY DIMENSIONALITY

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, LENC = 10 )
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      DIMENSION IC(LENC),C(LENC)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT

      SD = 0.001D0
      IAUX = 0
      IVF = 0
      IGRT = 0

      IF (IDIM .EQ. 1  .OR.  IDIM .EQ. 2) THEN
        DO 100 ISN = 1,NSTA
          IF (IDIM .EQ. 1) THEN
            KIND = 1
            NCON = NCON+1
            IOBS = IOBS+1
            CALL GETGLA (GLATB,ISN,B)
            CALL FORMIC (KIND,ISN,IDUMMY,IDUMM2,IC,LENG,IGRT)
            CALL FORMC (KIND,C,B,IDUMM1,IDUMM2,IDUMM3,IGRT)
            CALL COMPOB (KIND,OBS0,B,GLATB,ISN,IDUMMY,IDUMM2,IGRT)
            OBSB = GLATB
            CMO = OBS0
            IF (IMODE .EQ. 0) CMO = 0.D0
            WRITE (IUO) KIND,ISN,IDUMMY,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                  IVF,IAUX,IGRT

            KIND = 2
            NCON = NCON+1
            IOBS = IOBS+1
            CALL GETGLO (GLONB,ISN,B)
            CALL FORMIC (KIND,ISN,IDUMMY,IDUMM2,IC,LENG,IGRT)
            CALL FORMC (KIND,C,B,IDUMM1,IDUMM2,IDUMM3,IGRT)
            CALL COMPOB (KIND,OBS0,B,GLONB,ISN,IDUMMY,IDUMM2,IGRT)
            OBSB = GLONB
            CMO = OBS0
            IF (IMODE .EQ. 0) CMO = 0.D0
            WRITE (IUO) KIND,ISN,IDUMMY,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                  IVF,IAUX,IGRT
          ELSE
            IF ( .NOT. L2HLF) THEN
              KIND = 3
              NCON = NCON+1
              IOBS = IOBS+1
              CALL GETGH (GHT0,ISN,B)
              CALL GETMSL (GMSL0,ISN,B)
              EHB = GMSL0+GHT0
              CALL FORMIC (KIND,ISN,IDUMMY,IDUMM2,IC,LENG,IGRT)
              CALL FORMC (KIND,C,B,IDUMM1,IDUMM2,IDUMM3,IGRT)
              CALL COMPOB (KIND,OBS0,B,ADUMMY,ISN,IDUMMY,IDUMM2,IGRT)
              OBSB = EHB
              CMO = OBS0-OBSB
              IF (IMODE .EQ. 0) CMO = 0.D0
              WRITE (IUO) KIND,ISN,IDUMMY,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                    IVF,IAUX,IGRT
            ELSEIF (L2HLF  .AND.  I2HLF(ISN) .EQ. 0) THEN
              KIND = 3
              NCON = NCON+1
              IOBS = IOBS+1
              CALL GETGH (GHT0,ISN,B)
              CALL GETMSL (GMSL0,ISN,B)
              EHB = GMSL0+GHT0
              CALL FORMIC (KIND,ISN,IDUMMY,IDUMM2,IC,LENG,IGRT)
              CALL FORMC (KIND,C,B,IDUMM1,IDUMM2,IDUMM3,IGRT)
              CALL COMPOB (KIND,OBS0,B,ADUMMY,ISN,IDUMMY,IDUMM2,IGRT)
              OBSB = EHB
              CMO = OBS0-OBSB
              IF (IMODE .EQ. 0) CMO = 0.D0
              WRITE (IUO) KIND,ISN,IDUMMY,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                    IVF,IAUX,IGRT
            ENDIF
          ENDIF
  100   CONTINUE
      ENDIF

      RETURN
      END
      SUBROUTINE DIRDMS (VAL,ID,IM,S)

*** CONVERT DIRECTION,ANGLE,AZIMUTH RADIANS TO DEG, MIN, SEC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /CONST/ PI, PI2, RAD

    1 IF (VAL.GT.PI+PI) THEN
        VAL = VAL-PI-PI
        GO TO 1
      ENDIF

    2 IF (VAL.LT.0.D0) THEN
        VAL = VAL+PI+PI
        GO TO 2
      ENDIF

      S = DABS(VAL*RAD)
      ID = IDINT(S)
      S = (S-ID)*60.D0
      IM = IDINT(S)
      S = (S-IM)*60.D0

      RETURN
      END
      SUBROUTINE DISTOB (BCARD, IUO, IOBS, B, NX, FATAL, LSN)

*** OBSERVATION EQUATIONS FOR DISTANCES

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      CHARACTER*80 BCARD
      CHARACTER*10 ADIS54
      CHARACTER*9 ADIS52
      CHARACTER*1 TC
      LOGICAL FATAL, LSN
      LOGICAL GETSSN, GETPRM, GETIVF
      LOGICAL ADDCON
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL REJECT,GET82
      DIMENSION B(*), NX(*)
      DIMENSION IC(LENC), C(LENC)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON/UNITS/ LUNIT

      READ (BCARD,1) IRT, ISSN, IYR, IMO, IDY, IHR, IMN, TC, JSSN
 1    FORMAT (7X, I2, 1X, I4, 20X, 5I2, A1, I4)
      IF (IRT .EQ. 52) THEN
        READ (BCARD,3) ADIS52
 3      FORMAT (63X, A9)
        CALL NBLANK (ADIS52, 4, IBLK)
        READ (ADIS52,5) OBSB
 5      FORMAT (F9.4)
      ELSE
        READ (BCARD,4) ADIS54
 4      FORMAT (63X, A10)
        CALL NBLANK (ADIS54, 3, IBLK)
        READ (ADIS54,6) OBSB
 6      FORMAT (F10.3)
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

        IF (BCARD(35:36) .EQ. '  ') IYR = 84
        IF (BCARD(37:38) .EQ. '  ') IMO = 1
        IF (BCARD(39:40) .EQ. '  ') IDY = 1
        IF (BCARD(41:42) .EQ. '  ') IHR = 0
        IF (BCARD(43:44) .EQ. '  ') IMN = 0
        IF (BCARD(45:45) .EQ. ' ') TC = 'Z'
        CALL GETYR(BCARD,IYR)
*       IYR = IYR+1900
        IF ( .NOT. GETSSN(ISSN, ISN) ) THEN
          CALL LINE (3)
          WRITE (LUNIT,2) BCARD
 2        FORMAT ('0NO *80* RECORD FOR--', A80, /)
        ELSEIF ( .NOT. GETSSN(JSSN, JSN) ) THEN
          CALL LINE (3)
          WRITE (LUNIT,2) BCARD
        ELSE

***   RETRIEVE A STD DEV

          CALL STDDEV (BCARD, IRT, SD, REJECT)

*** KIND = 7 MARK TO MARK DISTANCE

          KIND = 7
          IOBS = IOBS+1
          NOBS = NOBS+1
          NDS = NDS+1
          LSN = .TRUE.
          CALL TOMNT (IYR, IMO, IDY, IHR, IMN, TC, ITIME)
          IF ( .NOT. GETPRM(52, ITIME, IAUX) ) IAUX = 0
          IF ( .NOT. GETIVF(52, ITIME, IVF) ) IVF = 0
          IGRT = 0
          CALL FORMIC (KIND, ISN, JSN, IAUX, IC, LENG, IGRT)
          CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
          CALL COMPOB (KIND, OBS0, B, OBSB, ISN, JSN, IAUX, IGRT)
          CMO = OBS0-OBSB
          IF (IMODE .EQ. 0) CMO = 0.D0
          CALL BIGV (KIND, ISN, JSN, IOBS, IVF, CMO, SD, FATAL)

*** UPDATE CONNECTIVITY

          IF ( .NOT. ADDCON(IC, LENG, NX) ) THEN
            WRITE (LUNIT,666)
 666        FORMAT ('0INSUFFICIENT STORAGE FOR REDUCED DIST.'/)
            CALL ABORT2
          ENDIF

          WRITE (IUO) KIND, ISN, JSN, IC, C, LENG, CMO, OBSB, SD, IOBS,
     &                IVF, IAUX, IGRT
        ENDIF

      RETURN
      END
      DOUBLE PRECISION FUNCTION DIVID (X,Y)

*** DIVIDE X BY Y -- ALLOW FOR Y = 0

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      IF (Y.EQ.0.D0) THEN
        DIVID = 0.D0
      ELSE
        DIVID = X/Y
      ENDIF

      RETURN
      END
      DOUBLE PRECISION FUNCTION DIVIDE (X,N)

*** DIVIDE X BY N -- ALLOW FOR N = 0

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      IF (N.EQ.0) THEN
        DIVIDE = 0.D0
      ELSE
        DIVIDE = X/N
      ENDIF

      RETURN
      END
      SUBROUTINE DOPCOV (ISSN, X, Y, Z, SIGN, SIGE, SIGU,
     &                   CORNE, CORNU, COREU, COVECF)

*** COMPUTE GEOCENTRIC CARTESIAN COV. MATRIX GIVEN LOCAL HORIZON VALUES

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION G(3,3), COVLHS(3,3), COVECF(3,3), WORK(3)
      DATA TOL /1.D-10/

*** OBTAIN LATITUDE AND LONGITUDE, CONSTRUCT TRANSFORMATION
*** (ROTATION) MATRIX

      CALL TOGEO2 (X, Y, Z, GLA, GLO)

      SLA = DSIN(GLA)
      CLA = DCOS(GLA)
      SLO = DSIN(GLO)
      CLO = DCOS(GLO)

      G(1,1) = -SLA*CLO
      G(1,2) = -SLO
      G(1,3) =  CLA*CLO

      G(2,1) = -SLA*SLO
      G(2,2) =  CLO
      G(2,3) =  CLA*SLO

      G(3,1) =  CLA
      G(3,2) =  0.D0
      G(3,3) =  SLA

*** CONSTRUCT THE LOCAL HORIZON COV MATRIX  (FULL, SQUARE)

      COVLHS(1,1) = SIGN*SIGN
      COVLHS(2,2) = SIGE*SIGE
      COVLHS(3,3) = SIGU*SIGU

      COVLHS(1,2) = CORNE*SIGN*SIGE
      COVLHS(1,3) = CORNU*SIGN*SIGU
      COVLHS(2,3) = COREU*SIGE*SIGU

      COVLHS(2,1) = COVLHS(1,2)
      COVLHS(3,1) = COVLHS(1,3)
      COVLHS(3,2) = COVLHS(2,3)

*** LINEAR ERROR PROPAGATION

      CALL ABAT (G, COVLHS, COVECF, WORK, 3, 3)

      RETURN
      END
      SUBROUTINE DUALHT (A,B,NX,SHIFTS,GOOGE,SIGUWT)

*** LIST DUAL HEIGHT DIFFERENCES

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LOCSSN
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      LOGICAL ELFLAG,DFFLAG
      CHARACTER*30 NAMES,NAME
      DIMENSION A(*),B(*),NX(*),SHIFTS(*),GOOGE(*)
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /CONST/  PI, PI2, RAD
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /FLAGS/  ELFLAG(MXSSN), DFFLAG(MXSSN)
      COMMON/UNITS/ LUNIT

*** HEADING

      CALL HEAD9
      NLINE = 2

*** GET POSITIONS

      DO 2 ISN = 1,NSTA
        IF (I2HLF(ISN) .EQ. 1) THEN
          CALL GETMSL (GMSL0,ISN,B)
          CALL GETGH (GHT0,ISN,B)
          EHT0 = GMSL0+GHT0
          CALL GETEHT (EHT1,ISN,B)
          IF (ELFLAG(ISN)) THEN
            GMSL1 = EHT1 - GHT0
            GHT1 = GHT0
          ELSE
            GHT1 = EHT1 - GMSL0
            GMSL1 = GMSL0
          ENDIF
          DIF = EHT0 - EHT1
          NAME = NAMES(ISN)
          IF ( .NOT. LOCSSN(ISN,ISSN)) THEN
            WRITE (LUNIT,666) ISN
  666       FORMAT ('0SSN TABLE ERROR IN DUALHT--',I8)
            CALL ABORT2
          ENDIF

          CALL LINE9 (NLINE)
          WRITE (LUNIT,3) ISSN, NAME, GMSL0, GHT0, EHT0,
     &                                 GMSL1, GHT1, EHT1, DIF
    3     FORMAT ('0', I5, 1X, A30, F9.3, F10.3, F11.3,
     &                          2X, F9.3, F10.3, F11.3, F11.3)

        ENDIF
    2 CONTINUE

      RETURN
      END
      SUBROUTINE DUMRD (IUO,NR,NC,G)

*** LOAD WORK SPACE

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION G(NR,NC)
      COMMON/UNITS/ LUNIT

      DO 1 I = 1, NR
        READ (IUO,END=666) ( G(I,J), J = 1, NC )
   1  CONTINUE

      RETURN

*** PREMATURE END OF FILE

  666 WRITE (LUNIT,667) NR
  667 FORMAT ('0PREMATURE FILE END IN DUMRD -- NR=',I5)
      CALL ABORT2
      RETURN
      END
      SUBROUTINE ECHOBB (BCARD,IOBS,LSN)

*** ECHO THE BLUE BOOK RECORDS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 BCARD
      LOGICAL LSN, LBV
      LOGICAL LEB, LLB, LEG, LLG, LED, LLD
      COMMON /PECHO/  LEB, LLB, LEG, LLG, LED, LLD
      COMMON /ECHO/   VSD, LBV
      COMMON/UNITS/ LUNIT

*** ECHO ALL OF THE BLUE BOOK

      IF ( .NOT. LEB) THEN
        IF ( .NOT. LLB) THEN
          IF (LSN) THEN
            IF ( BCARD(6:6) .NE. 'R'  .AND.  BCARD(6:6) .NE. 'O'  .AND.
     &           BCARD(6:6) .NE. 'F') THEN
              IF (LBV) THEN
                CALL LINE (2)
                WRITE (LUNIT,1) IOBS,VSD,BCARD
 1              FORMAT (1X,I6,3X,F10.1,' *MISCLOSURE',/,10X,A80)
              ELSE
                CALL LINE (1)
                WRITE (LUNIT,3) IOBS,BCARD
 3              FORMAT (1X,I6,3X,A80)
              ENDIF
            ELSE
              CALL LINE (1)
              WRITE (LUNIT,2) BCARD
            ENDIF
          ELSE
            CALL LINE (1)
            WRITE (LUNIT,2) BCARD
 2          FORMAT (10X,A80)
          ENDIF
        ENDIF
      ENDIF

*** ECHO OBSERVATIONS ONLY

      IF (LEB) THEN
        IF (LSN) THEN
          IF ( BCARD(6:6) .NE. 'R'  .AND.  BCARD(6:6) .NE. 'O'  .AND.
     &         BCARD(6:6) .NE. 'F') THEN
            IF (LBV) THEN
              CALL LINE (2)
              WRITE (LUNIT,1) IOBS,VSD,BCARD
            ELSE
              CALL LINE (1)
              WRITE (LUNIT,3) IOBS,BCARD
            ENDIF
          ENDIF
        ENDIF
      ENDIF

*** ECHO LARGE MISCLOSURES ONLY

      IF (LLB) THEN
        IF (LSN) THEN
          IF ( BCARD(6:6) .NE. 'R'  .AND.  BCARD(6:6) .NE. 'O'  .AND.
     &         BCARD(6:6) .NE. 'F') THEN
            IF (LBV) THEN
              CALL LINE (2)
              WRITE (LUNIT,1) IOBS,VSD,BCARD
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END
      SUBROUTINE ECHODF (DCARD, IOBS, LSN, AREJ)

*** ECHO THE DFILE RECORDS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*89 DCARD
      CHARACTER*1 AREJ
      LOGICAL LSN, LBV
      LOGICAL LEB, LLB, LEG, LLG, LED, LLD
      COMMON /PECHO/  LEB, LLB, LEG, LLG, LED, LLD
      COMMON /ECHO/   VSD, LBV
      COMMON/UNITS/ LUNIT

*** ECHO ALL OF THE DFILE

      IF ( .NOT. LED) THEN
        IF ( .NOT. LLD) THEN
          IF (LSN) THEN
            IF (AREJ .NE. 'R'  .AND . AREJ .NE. 'O'  .AND.
     &          AREJ .NE. 'F') THEN
              CALL LINE (1)
              WRITE (LUNIT,3) IOBS, DCARD
 3            FORMAT (1X, I6, 3X, A89)
            ELSE
              CALL LINE (1)
              WRITE (LUNIT,2) DCARD
            ENDIF
          ELSE
            CALL LINE (1)
            WRITE (LUNIT,2) DCARD
 2          FORMAT (10X, A89)
          ENDIF
        ENDIF
      ENDIF

*** ECHO OBSERVATIONS ONLY

      IF (LED) THEN
        IF (LSN) THEN
          IF (AREJ .NE. 'R'  .AND.  AREJ .NE. 'O'  .AND.
     &        AREJ .NE. 'F') THEN
            CALL LINE (1)
            WRITE (LUNIT,3) IOBS, DCARD
          ENDIF
        ENDIF
      ENDIF

*** ECHO LARGE MISCLOSURES ONLY

      IF (LLD) THEN
        IF (LSN) THEN
          IF (AREJ .NE. 'R'  .AND.  AREJ .NE. 'O'  .AND.
     &        AREJ .NE. 'F') THEN
            IF (LBV) THEN
              CALL LINE (1)
              WRITE (LUNIT,3) IOBS, DCARD
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END
      SUBROUTINE EFINIT

*** INITIALIZE THE ELEVATION ADJUSTMENT FLAGS AND RECONCILE FLAGS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7 ADJFIL, NAFILE
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL ELFLAG,DFFLAG
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /FLAGS/  ELFLAG(MXSSN), DFFLAG(MXSSN)
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP

*** INITIALIZE HEIGHT ADJ. FLAGS (TRUE = ADJUST MSL)

      DO 1 I = 1, MXSSN
        ELFLAG(I) = LMSL
    1 CONTINUE

*** SET IMODE = 3 FOR VAR. FACTORS

      IF (NVFTOT .GT. 0) IMODE = 3
      IF (IMODE .EQ. 0) LPS = .FALSE.

*** IF NO GPS OR DOPPLER FILES, SET FLAGS

      IF (GFILE .EQ. 'NOGFILE') LGPS = .TRUE.
      IF (DFILE .EQ. 'NODFILE') LDOP = .TRUE.

      RETURN
      END
      SUBROUTINE EQ (A, B, N, M)

*** SET ONE MATRIX EQUAL TO ANOTHER

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION A(N,M), B(N,M)

      DO 10 I = 1, N
        DO 5 J = 1, M
          B(I,J) = A(I,J)
    5   CONTINUE
   10 CONTINUE

      RETURN
      END
      
      SUBROUTINE FIR12 (BCARD)
      
*** FIRST PASS OF 12 RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 BCARD
      LOGICAL VALCH1
      COMMON/UNITS/ LUNIT
      COMMON /YEARS/ ICBEG,IYBEG,ICEND,IYEND
      
      IF(.NOT.VALCH1(BCARD,11,12,9)) THEN
         WRITE(LUNIT,10) BCARD
 10      FORMAT(/' ',A80/' *** ERROR - B-FILE *12* CENTURY FIELD',
     *   ' OBS BEGAB CC 11-12 IS INVALID'/)
         CALL ABORT2
      ELSEIF(.NOT.VALCH1(BCARD,17,18,9)) THEN
         WRITE(LUNIT,15) BCARD
 15      FORMAT(/' ',A80/' *** ERROR - B-FILE *12* CENTURY FIELD',
     *   ' OBS END CC 17-18 IS INVALID'/)
         CALL ABORT2
      ELSEIF(.NOT.VALCH1(BCARD,13,14,1)) THEN
         WRITE(LUNIT,20) BCARD
 20      FORMAT(/' ',A80/' *** ERROR - B-FILE *12* YEAR FIELD',
     *   ' OBS BEGAN CC 13-14 IS NON-INTEGER'/)
         CALL ABORT2
      ELSEIF(.NOT.VALCH1(BCARD,19,20,1)) THEN
         WRITE(LUNIT,25) BCARD
 25      FORMAT(/' ',A80/' *** ERROR - B-FILE *12* YEAR FIELD',
     *   ' OBS END CC 19-20 IS NON-INTEGER'/)
         CALL ABORT2
         ELSE
         READ (BCARD,30) ICBEG,IYBEG,ICEND,IYEND
 30      FORMAT (10X,I2,I2,2X,I2,I2)
 
         IF(((ICEND*100)+IYEND).LT.((ICBEG*100)+IYBEG)) THEN
           WRITE (LUNIT,35) BCARD
 35        FORMAT (/' ',A80/' ***ERROR - B-FILE *12* END YEAR',
     *     ' CC 17-20 IS LESS THAN BEGAN YEAR CC 11-14'/)
           CALL ABORT2
         ENDIF
         IF(((ICEND*100)+IYEND)-((ICBEG*100)+IYBEG).GE. 100) THEN
           WRITE (LUNIT,40) BCARD
 40        FORMAT (/' ',A80/' ***ERROR - B-FILE *12* END YEAR',
     *     ' CC 17-20 EXCEEDS BEGAN YEAR CC 11-14 BY MORE THAN 99'/)
           CALL ABORT2
         ENDIF
         ENDIF
      
         RETURN
         END
          
      SUBROUTINE FIR20 (BCARD)

*** FIRST PASS OF DIRECTION SET RECORDS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 BCARD
      LOGICAL DOIT, DUNIT
      LOGICAL GET82
      COMMON /ZUNKS/  ISSN, LIST, DOIT, DUNIT

      DOIT = .FALSE.
      DUNIT = .FALSE.

      READ (BCARD,1) ISSN, LIST, JSSN
 1    FORMAT (10X, I4, I2, 34X, I4)

      IF (GET82(ISSN,I)) RETURN
      IF (GET82(JSSN,I)) RETURN

      DOIT = .TRUE.

      RETURN
      END
      SUBROUTINE FIR22 (BCARD, FATAL)

*** FIRST PASS OF DIRECTION RECORDS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 BCARD
      LOGICAL DOIT, DUNIT, FATAL, PUTZ
      LOGICAL GET82
      COMMON /ZUNKS/  ISSN, LIST, DOIT, DUNIT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      READ(BCARD,1) ISSN, JSSN
    1 FORMAT(10X,I4,36X,I4)
      IF(GET82(ISSN,I)) RETURN
      IF(GET82(JSSN,I)) RETURN

      IF (DUNIT) RETURN

      IF ( .NOT. DOIT) THEN
        DOIT = .TRUE.
      ELSE
        DUNIT = .TRUE.
        IF ( .NOT. PUTZ(ISSN, LIST, IDUP) ) THEN
          IF (IDUP .EQ. 0) THEN
            CALL LINE (3)
            WRITE (LUNIT,666) BCARD, ISSN
 666        FORMAT ('0', A80, I5.3, ' TOO MANY LISTS -- FATAL', /)
            FATAL = .TRUE.
          ELSE
            CALL LINE (3)
            WRITE (LUNIT,667) BCARD, ISSN, IDUP
 667        FORMAT ('0', A80, I5.3, ' DUPLICATES THE', I5, 'ENTRY', /)
            FATAL = .TRUE.
          ENDIF
        ELSE
          NZ = NZ + 1
        ENDIF
      ENDIF

      RETURN
      END

C---------------------------------------------------------------------------------------
      SUBROUTINE FIR80 (BCARD, B, FATAL)

*** FIRST ENCOUNTER OF CONTROL POINT RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, MXFOT = 3500 )
      LOGICAL FATAL, PUTSSN
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      CHARACTER*80 BCARD
      CHARACTER*30 NAME, NAMES
      CHARACTER*7 ASLA, ASLO
      CHARACTER*6 AGMSL,PIDs,staPID
      CHARACTER*1 ADLA, ADLO, OT
      CHARACTER*2 state,ST          
      DIMENSION B(*)
      COMMON /OPT/ AX,E2,DMSL,DGH,VM,VP,CTOL,ITMAX,IMODE,
     &             LMSL,LSS,LUP,LABS,LADJ,LUPI
      COMMON /NAMTAB/ NAMES(MXSSN)
      COMMON /STRUCT/ NSTA,NAUX,NUNK,IDIM,NSTAS,NOBS,NCON,NZ,NGRT
      COMMON /FOURTH/ NFOT, NFOTS(MXFOT)
      COMMON /EXTRMA/ GLAMAX, GLAMIN, GLOMAX, GLOMIN, HTMAX, HTMIN
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT
      COMMON/PIDs / PIDs (MXSSN),ST(MXSSN)

      READ(BCARD,1,ERR=666,IOSTAT=IOS) staPID,ISSN,NAME,IDLA,IMLA,ASLA,
     &                                 ADLA,IDLO,IMLO,ASLO,ADLO,
     &                                 AGMSL,OT,State
    1 FORMAT(a6,4X,I4,A30,2I2,A7,A1,I3,I2,A7,A1,A6,3X,A1,a2)
      CALL NBLANK (AGMSL, 2, IBLK)
      READ (AGMSL,6,ERR=666,IOSTAT=IOS) GMSL
    6 FORMAT (F6.2)
      IF (BCARD(70:75) .EQ. '      ') GMSL = DMSL

      CALL NBLANK (ASLA, 5, IBLK)
      CALL NBLANK (ASLO, 5, IBLK)
      READ (ASLA,5,ERR=666,IOSTAT=IOS) ISLA
      READ (ASLO,5,ERR=666,IOSTAT=IOS) ISLO
    5 FORMAT (I7)

      IF ( .NOT. PUTSSN(ISSN, IDUP) ) THEN
        IF (IDUP .EQ. 0) THEN
          CALL LINE (3)
          WRITE (LUNIT,2) ISSN
    2     FORMAT ('0', I5, ' IS ILLEGAL SSN -- FATAL', /)
          FATAL = .TRUE.
        ELSE
          CALL LINE (3)
          WRITE (LUNIT,3) ISSN, IDUP
    3     FORMAT ('0', I5, ' DUPLICATES THE ', I5,
     &            '-TH ENTRY -- FATAL', /)
          FATAL = .TRUE.
        ENDIF

      ELSE
        NSTA       = NSTA + 1
        PIDs(NSTA) = staPID
        ST  (NSTA) = state
        IF (OT .EQ. '4') THEN
          NFOT = NFOT + 1
          IF (NFOT .GT. MXFOT) THEN
            WRITE (LUNIT,4) MXFOT
    4       FORMAT ('0OVER', I4, ' LANDMARKS')
            CALL ABORT2
          ENDIF
          NFOTS(NFOT) = NSTA
        ENDIF

*** IF 2.5DIM ADJUSTMENT, GET ELLIP HEIGHTS EITHER FROM EXTENDED
*** *80* RECORDS OR INITIALIZE WITH ORTHOMETRIC HEIGHTS - BIGADJUST

        IF (L2HLF) THEN
*         IF (LEHT) THEN
*           READ (BCARD,20,ERR=30) EHT
*  20       FORMAT (82X, F6.2)
*         ELSE
*           EHT = GMSL
*         ENDIF
*         GOTO 40

*** TAKE CARE OF THE POSSIBILITY THAT *80* RECORD DOES NOT HAVE EHT

   30     EHT = GMSL
*  40     CONTINUE

        ENDIF

*** LATITUDES ARE POSITIVE NORTH IN RADIANS
*** LONGITUDES ARE POSITIVE EAST IN RADIANS

        IF (ADLA .EQ. 'S') THEN
          GLASGN = -1.D0
        ELSE
          GLASGN = 1.D0
        ENDIF
        IF (ADLO .EQ. 'W') THEN
          GLOSGN = -1.D0
        ELSE
          GLOSGN = 1.D0
        ENDIF
        CALL GETRAD (IDLA, IMLA, ISLA, GLASGN, GLA)
        CALL GETRAD (IDLO, IMLO, ISLO, GLOSGN, GLO)

*** LOAD CONTROL POINT VALUES AND DEFAULTS

        CALL PUTALA (GLA, NSTA, B)
        CALL PUTALO (GLO, NSTA, B)
        CALL PUTGLA (GLA, NSTA, B)
        CALL PUTGLO (GLO, NSTA, B)
        CALL PUTMSL (GMSL, NSTA, B)
        CALL PUTGH (DGH, NSTA, B)
        IF (L2HLF) CALL PUTEHT (EHT, NSTA, B)
        I2HLF(NSTA) = 0
        NAMES(NSTA) = NAME

*** UPDATE MAXIMUM AND MINIMUM VALUES

        IF (GLA .GT. GLAMAX ) GLAMAX = GLA
        IF (GLA .LT. GLAMIN ) GLAMIN = GLA
        IF (GLO .GT. GLOMAX ) GLOMAX = GLO
        IF (GLO .LT. GLOMIN ) GLOMIN = GLO
        IF (GMSL .GT. HTMAX ) HTMAX = GMSL
        IF (GMSL .LT. HTMIN ) HTMIN = GMSL
      ENDIF

  990 RETURN
  666 WRITE (LUNIT,667) IOS, BCARD
  667 FORMAT (//, ' FORTRAN ERROR #', I5, ' IN SUBROUTINE FIR80 WHEN',          
     &           ' READING THE FOLLOWING RECORD', A80)
      GOTO 990
      END
C----------------------------------------------------------------------------------
      SUBROUTINE FIR82 (BCARD, FATAL)

*** FIRST ENCOUNTER OF *82* RECORD

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL PUT82,FATAL
      CHARACTER*80 BCARD
      COMMON/UNITS/ LUNIT

      READ(BCARD,1) ISSN
    1 FORMAT(10X,I4)
      IF(.NOT.PUT82(ISSN,IDUP)) THEN
        IF(IDUP.EQ.0) THEN
          CALL LINE (3)
          WRITE(LUNIT,2) ISSN
    2     FORMAT('0',I5,' IS ILLEGAL SSS -- FATAL',/)
          FATAL = .TRUE.
        ELSE
          CALL LINE (3)
          WRITE(LUNIT,3) ISSN, IDUP
    3     FORMAT('0',I5,' DUPLICATES THE ',I5,
     &           '-TH ENTRY -- FATAL'/)
          FATAL = .TRUE.
        ENDIF
      ENDIF

      RETURN
      END
      SUBROUTINE FIR84 (BCARD, B)

*** FIRST ENCOUNTER OF GEOID HEIGHT RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      CHARACTER*80 BCARD
      LOGICAL GETSSN
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      DIMENSION B(*)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

      READ (BCARD,1,ERR=666,IOSTAT=IOS) ISSN, GH
    1 FORMAT (10X, I4, 55X, F6.2)
      IF (BCARD(70:75) .EQ. '     ') GH = DGH

      IF ( GETSSN(ISSN, ISN) ) THEN
        CALL PUTGH (GH, ISN, B)
        N84 = N84 + 1
      ELSE
        CALL LINE (3)
        WRITE (LUNIT,2) BCARD
    2   FORMAT ('0NO *80* RECORD FOR --', A80, /)
        GOTO 990
      ENDIF

*** 2.5DIM ADJUSTMENT AND THE ELLIPSOID HEIGHTS ARE NOT FROM
*** THE EXTENDED *80* RECORD BUT ARE INITIALIZED TO BE THE
*** SUM OF THE ORTHOMETRIC AND GEOID HEIGHTS 

      IF (L2HLF  .AND.  ( .NOT. LEHT) ) THEN
        CALL GETEHT (EHT, ISN, B)
        EHT = EHT - GH
        CALL PUTEHT (EHT, ISN, B)
      ENDIF

  990 RETURN
  666 WRITE (LUNIT,667) IOS, BCARD
  667 FORMAT (//, ' FORTRAN ERROR #', I5, ' IN SUBROUTINE FIR84 WHEN',         
     &            ' READING THE FOLLOWING RECORD', A80)
      GOTO 990
      END
      SUBROUTINE FIR85 (BCARD, B)

*** FIRST ENCOUNTER OF DEFLECTION RECORD

*****************************************************************
*  CODE ASSUMES DATUM IS NAD 1983 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
******************************************************************

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      LOGICAL GETSSN
      CHARACTER*80 BCARD
      CHARACTER*1 ADXI, ADETA
      CHARACTER*5 AXI, AETA
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL ELFLAG, DFFLAG
      DIMENSION B(*)
      COMMON /FLAGS/  ELFLAG(MXSSN), DFFLAG(MXSSN)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /CONST/  PI, PI2, RAD
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON/UNITS/ LUNIT

      READ (BCARD,1,ERR=666,IOSTAT=IOS) ISSN, AXI, ADXI, AETA, ADETA
    1 FORMAT (10X, I4, 48X, A5, A1, 3X, A5, A1)

*** POS. MER. DEFL. -> ASTRO NORTH OF GEOD.
*** POS. PRIME V.   -> ASTRO EAST OF GEOD.
*** (NOTE: PRIME VERTICAL OPPOSITE SENSE OF DEFLECTION CHART)

      CALL NBLANK (AXI, 2, IBLK)
      READ (AXI,5,ERR=666,IOSTAT=IOS) XI
      CALL NBLANK (AETA, 2, IBLK)
      READ (AETA,5,ERR=666,IOSTAT=IOS) ETA
    5 FORMAT (F5.2)

      IF (ADXI .EQ. 'S') XI = -XI
      IF (ADETA .EQ. 'W') ETA = -ETA
      IF (GETSSN(ISSN, ISN) ) THEN
        XI = XI/(3600.D0*RAD)
        ETA = ETA/(3600.D0*RAD)
        CALL GETGLA (GLAT, ISN, B)
        CALL GETGLO (GLON, ISN, B)

        ALAT = GLAT + XI
        ALON = GLON + ETA/DCOS(GLAT)

        CALL PUTALA (ALAT, ISN, B)
        CALL PUTALO (ALON, ISN, B)
        DFFLAG(ISN) = .TRUE.
        N85 = N85 + 1
      ELSE
        CALL LINE (3)
        WRITE (LUNIT,2) BCARD
    2   FORMAT ('0NO *80* RECORD FOR --', A80, /)
      ENDIF

  990 RETURN
  666 WRITE (LUNIT,667) IOS, BCARD
  667 FORMAT (//, ' FORTRAN ERROR #', I5, ' IN SUBROUTINE FIR85 WHEN',
     &           ' READING THE FOLLOWING RECORD', A80)
      GOTO 990
      END
C--------------------------------------------------------------------------------------------------
      SUBROUTINE FIR86 (BCARD, B)

*** FIRST ENCOUNTER OF HEIGHT RECORD

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      LOGICAL GETSSN
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      CHARACTER*80 BCARD
      DIMENSION B(*)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

      READ (BCARD,1,ERR=666,IOSTAT=IOS) ISSN, GMSL, GH, EHT
    1 FORMAT (BZ, 10X, I4, 2X, F7.3, 12X, F7.3, 3X, F7.3)

*** DEFAULTS

      IF (BCARD(17:23) .EQ. '      ') GMSL = DMSL
      IF (BCARD(36:42) .EQ. '      ')   GH = DGH
      IF (BCARD(46:52) .EQ. '      ')  EHT = GMSL + GH

*** LOAD CONTROL POINT VALUES OR DEFAULTS

      IF ( GETSSN(ISSN, ISN) ) THEN
        CALL PUTMSL (GMSL, ISN, B)
        CALL PUTGH (GH, ISN, B)
        IF (L2HLF) CALL PUTEHT (EHT, ISN, B)
        N86 = N86 + 1
      ELSE
        CALL LINE (3)
        WRITE (LUNIT,2) BCARD
    2   FORMAT ('0NO *80* RECORD FOR --', A80, /)
        GOTO 990
      ENDIF

  990 RETURN

  666 WRITE (LUNIT,667) IOS, BCARD
  667 FORMAT (//, ' FORTRAN ERROR #',  I5, ' IN SUBROUTINE FIR86 WHEN',
     &           ' READING THE FOLLOWING RECORD', A80)
      GOTO 990
      END
C------------------------------------------------------------------------------------
      SUBROUTINE FIRST (B)

*** ROUTINE FOR FIRST PASS OF BLUE-BOOK, GFILE, AND DFILE
*** DETERMINE NUM OF STATIONS, PARAMETERS, AND PROBLEM STRUCTURE

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( NVECS = 700, MXFOT = 3500, MXSSN = 9999 )
      PARAMETER ( MAXB =49999 )

      LOGICAL FATAL, GETSSN
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL L2HLF, LEHT
      CHARACTER*89 DCARD
      CHARACTER*80 BCARD, GCARD
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7 ADJFIL, NAFILE
      CHARACTER*2 BBID, CC2,ST
      CHARACTER*1 CC1
      CHARACTER*6 PIDs
      DIMENSION B(*)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /FOURTH/ NFOT, NFOTS(MXFOT)
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      COMMON /GPS/    MAXVEC
      COMMON /GRTTB1/ GLAT0, GLON0, R0(3,3)
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT
      COMMON/BREC / NCCNT (MAXB)
      COMMON/PIDs / PIDs (MXSSN),ST(MXSSN)

*** PROCESS BLUE BOOK
**v6.0
      do i=1,MXSSN
        PIDs(i) = '      '
      enddo
**so far for v6.0
      IUNIT = 12
      OPEN (IUNIT,ERR=666,STATUS='OLD',FILE=BBOOK)
      CALL NEWSSN
      NSTA = 0
      FATAL = .FALSE.
      CALL NEWZ
      NZ = 0
      CALL NEW82
      CALL NEWFOT
      NFOT = 0

    1 READ (IUNIT,11,END=2) BCARD
      READ (BCARD,12) BBID
      IF (BBID .EQ. '82') CALL FIR82 (BCARD, FATAL)
      GO TO 1
    2 REWIND (IUNIT)

   10 READ (IUNIT,11,END=100) BCARD
   11 FORMAT (A80)
      READ (BCARD,12) BBID
   12 FORMAT (7X, A2)

      IF (BBID .EQ. '80') THEN
        CALL FIR80 (BCARD, B, FATAL)
      ELSEIF (BBID .EQ. '84') THEN
** v 4.29k
*       CALL FIR84 (BCARD, B)
            WRITE (LUNIT, 13) BCARD
   13       FORMAT (/, ' WARNING -  *84* RECORD IGNORED - IT IS ',
     &             'NO LONGER USED IN THE BBOOK FILE' /, A80)
*****
      ELSEIF (BBID .EQ. '85') THEN
        CALL FIR85 (BCARD, B)
      ELSEIF (BBID .EQ. '86') THEN
        CALL FIR86 (BCARD, B)
      ELSEIF (BBID .EQ. '20') THEN
        IF ( .NOT. LDIR) CALL FIR20 (BCARD)
      ELSEIF (BBID .EQ. '22') THEN
        IF ( .NOT. LDIR) CALL FIR22 (BCARD, FATAL)

c ADJUST 5.3 the following behavior is not desired and has
c been neutralized 4/2/08 Mike Potterfield
**** v 4.30VF  12-22-03 **********
c      ELSEIF((BBID.EQ.'91').OR.(BBID.EQ.'92').OR.(BBID.EQ.'93'))THEN
c       WRITE (LUNIT,867) BCARD
c  867   FORMAT (/,1X,A80,/,' INPUT BLUE BOOK FILE *91*, *92* OR ',
c     &  '*93* RECORDS ENCOUNTERED -- FATAL') 
c        CALL ABORT2
************************
c End of changes 4/2/08
      ENDIF

      GO TO 10
  100 CLOSE (IUNIT)

*** ERROR CHECK

      IF (NSTA .LE. 0) THEN
        WRITE (LUNIT,667)
  667   FORMAT ('0NO *80* CONTROL POINT RECORDS ENCOUNTERED--FATAL')
        CALL ABORT2
      ENDIF
      IF (FATAL) THEN
        WRITE (LUNIT,668)
  668   FORMAT ('0FATAL ERROR FLAG DUE TO PREVIOUS ERRORS')
        CALL ABORT2
      ENDIF

*** PROCESS GFILE TO DETERMINE MAX VECTOR SIZE AND 2.5 DIM STATIONS


***v 4.28k
         IF ( .NOT. LGPS) THEN
          if(gfile .ne. 'nogfile') then
          if(gfile .ne. 'NOGFILE') then
*************

        IGPS = 13
       
        OPEN (IGPS,ERR=920,STATUS='OLD',FILE=GFILE)

*** READ THE MAX NUMBER OF GPS VECTORS IN A GROUP

        NOB = 0
        DO 5 I = 1, MAXB
        NCCNT(I) = 0
    5   CONTINUE
        MAXVEC = 0
   20   READ (IGPS,23,END=200) GCARD
   23   FORMAT (A80)
**v6.3b
        if (GCARD(1:1) == ' ') goto 20
**so far for v6.3b

        CC1 = GCARD(1:1)

**v6.2.2
        IF (CC1 .EQ. 'A' .and. GCARD(79:80) .EQ. 'ZT') THEN
           WRITE (LUNIT,'(/)')
           WRITE (LUNIT,*)
     &'***************************************************************'

           WRITE (LUNIT, 227)
 227       FORMAT
     &(' NOTE:THE INPUT GPS VECTORS FILE (GFILE) WAS TRANSFORMED BY HTDP
     &')
           WRITE (LUNIT,*)
     &'***************************************************************'
           WRITE (LUNIT,'(/)')

        ENDIF
**so far for v6.2.2

        IF (CC1 .EQ. 'B') THEN
          NOD = 0
          NOC = 0
          NOB = NOB + 1

*****12-02-03*******
         IF(NOB .GE. MAXB) THEN
           WRITE (LUNIT, 228) MAXB
 228       FORMAT(/,'NUMBER OF SESSIONS EXCEED ',I8,' FATAL ERROR')
           CALL ABORT2
          ENDIF
************************
          READ (GCARD,21) MAX
   21     FORMAT (25X, I2)
          IF (MAX .GT. MAXVEC) MAXVEC = MAX
        ENDIF

*** COUNT NUMBER OF C/F REC IN EACH GRP AND STORE IN NCCNT ARRAY***
        
*****3-28*************************************************************
***del           IF(CC1 .EQ. 'C' .OR. CC1 .EQ. 'F') NOC = NOC + 1
***add
             IF(CC1 .EQ. 'C' .OR. CC1 .EQ. 'F') THEN 
                  NOC = NOC + 1
                  NCCNT(NOB) = NOC
                  IF(NOC .GT. MAXVEC) MAXVEC = NOC
             ENDIF
**********************************************************************

*** IF 2.5DIM ADJUSTMENT, SET ISN'S IN I2HLF TO 1

* v 4.28j***********
*       IF ( L2HLF  .AND.  CC1 .EQ. 'C' ) THEN
        if (l2hlf) then
         if( cc1 .eq. 'F' .or. cc1 .eq. 'C') then
*********************
          READ (GCARD,27) ISSN, JSSN
   27     FORMAT (1X, 2I4)
          IF ( GETSSN(ISSN, ISN)  .AND.  GETSSN(JSSN, JSN) ) THEN
            I2HLF(ISN) = 1
            I2HLF(JSN) = 1
          ELSE
            WRITE (LUNIT, 28) GCARD
   28       FORMAT (/, 'ILLEGAL SSN ON FOLLOWING CARD, FATAL ERROR',
     &              /, A80)
            CALL ABORT2
          ENDIF
         endif
        ENDIF
*******3-28 del**************************************
*****       IF(CC1 .EQ. 'D' .OR. CC1 .EQ. 'E') THEN
******         IF(NOD .EQ. 0) THEN
******             NCCNT(NOB) = NOC
*****             NOD = NOD + 1
*****             IF (NOC .GT. MAXVEC) MAXVEC = NOC
****          ELSE
****             NOD = NOD + 1
****          ENDIF
*****        ENDIF
**********************************************

        GO TO 20

*** ERROR TEST ON MAXVEC

  200   IF (MAXVEC .EQ. 0  .OR.  MAXVEC .GT. NVECS) THEN
          WRITE (LUNIT,22) MAXVEC, NVECS
   22     FORMAT (' MAXVEC =', I3, ' MAX GPS VECTORS MUST BE ',
     &            'GREATER THAN ZERO AND LESS THAN', I3, '.')
          CALL ABORT2
        ENDIF
  210   CLOSE (IGPS)
      ENDIF
***v 4.28k
      endif
      endif
***************

*** PROCESS DFILE TO DETERMINE 2.5 DIM STATIONS

  290 IF (.NOT. LDOP) THEN
      if(dfile .ne. 'nodfile') then
      if(dfile .ne. 'NODFILE') then
        IDOP = 14
        OPEN (IDOP,ERR=930,STATUS='OLD',FILE=DFILE)
        IF (L2HLF) THEN
  120     READ (IDOP,123,END=300) DCARD
  123     FORMAT (A89)
          CC2 = DCARD(1:2)
          IF ( CC2 .EQ. 'DP' ) THEN
            READ (DCARD,127) ISSN
  127       FORMAT (3X, I5)
            IF ( GETSSN(ISSN, ISN) ) THEN
              I2HLF(ISN) = 1
            ELSE
              WRITE (LUNIT, 128) DCARD
  128         FORMAT (/, 'ILLEGAL SSN ON FOLLOWING CARD, FATAL ERROR',
     &                /, A89)
              CALL ABORT2
            ENDIF
          ENDIF
          GOTO 120
        ENDIF
  300   CLOSE (IDOP)
      ENDIF
***v 4.28k
      endif
      endif
***************

*** DETERMINE NUMBER OF DUAL HEIGHT STATIONS

  330 N2HLF = 0
      IF (L2HLF) THEN
        DO 320 I = 1, NSTA
          IF (I2HLF(I) .EQ. 1) N2HLF = N2HLF + 1
  320   CONTINUE
        IF (N2HLF .EQ. 0) THEN
          CALL LINE (3)
          WRITE (LUNIT,310)
  310     FORMAT (/, ' WARNING - NO DUAL HEIGHT STATIONS', /)
        ENDIF
      ENDIF

*** CONVERT GLAT, GLON, ELLIP HT (SUM OF ORTHOMETRIC HT AND GEOID HT)
*** TO EARTH CENTERED (E.C.F.) X,Y,Z.  IF 2.5DIM ADJUSTMENT ALSO CONVERT
*** GLAT, GLON, AND ELLIP HT (DUAL HT) TO X,Y,Z IN SUBROUTINE ALLXYZ.
*** OBTAIN AVERAGE (E.G. WEIGHTED CENTER) LOCATION FOR GPS AND
*** DOPPLER ROTATIONS

  664 CALL ALLXYZ (B, GLAAVE, GLOAVE)
      IF (NGRT .GT. 0) THEN
        GLAT0 = GLAAVE
        GLON0 = GLOAVE
        SLAT = DSIN(GLAT0)
        CLAT = DCOS(GLAT0)
        SLON = DSIN(GLON0)
        CLON = DCOS(GLON0)

        R0(1,1) = -SLAT * CLON
        R0(1,2) = -SLAT * SLON
        R0(1,3) = CLAT

        R0(2,1) = -SLON
        R0(2,2) = CLON
        R0(2,3) = 0.D0

        R0(3,1) = CLAT * CLON
        R0(3,2) = CLAT * SLON
        R0(3,3) = SLAT
      ENDIF

      RETURN

*** NO BLUE BOOK -- FATAL ERROR

  666 WRITE (LUNIT,669)
  669 FORMAT ('0NO BLUE BOOK FOUND--FATAL ERROR')
      CALL ABORT2
      RETURN

*** ERROR MESSAGES - NO GFILE OR NO DFILE

  920 WRITE (LUNIT,925)
***v 4.28k
****  925 FORMAT ('0ERROR - NO GFILE FOUND IN SUBROUTINE FIRST')
  925 FORMAT ('0NO GFILE FOUND IN SUBROUTINE FIRST--FATAL ERROR')
      call abort2
      return
***     CLOSE (IGPS)
***     LGPS = .TRUE.
***     GOTO 290

  930 WRITE (LUNIT,935)
***  935 FORMAT ('0ERROR - NO DFILE FOUND IN SUBROUTINE FIRST')
  935 FORMAT ('0NO DFILE FOUND IN SUBROUTINE FIRST--FATAL ERROR')
      call abort2
      return
***      CLOSE (IDOP)
***      LDOP = .TRUE.
***      GOTO 330
************************
      END
C--------------------------------------------------------------------------------------------
      LOGICAL FUNCTION FIXVF (I)

*** RETURN TRUE IF FIXED VARIANCE FACTOR

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXVF = 40 )
      LOGICAL LFIXS
      COMMON /IVFTB1/ NIVFS, ICODES(MXVF), IOLDS(MXVF), INEWS(MXVF),
     &                LFIXS(MXVF)
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NVFTOT) THEN
        WRITE (LUNIT,1) I,NVFTOT
    1   FORMAT ('0********* ILLEGAL IVF=',I10,' FOR N=',I10,' IN FIXVF')
        CALL ABORT2
      ENDIF
      FIXVF = LFIXS(I)

      RETURN
      END
      SUBROUTINE FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)

*** COMPUTE THE OBS EQUATION COEFFICIENTS FROM THE PARAMETERS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, LENC = 10 )
      LOGICAL L2HLF, LEHT
      DIMENSION RI(3,3), RJ(3,3), RT(3,3)
      DIMENSION B(*)
      DIMENSION C(LENC), CL(LENC-4), CR(LENC-4)
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT
      SAVE RI, RJ

*** KIND 0 IS AN AUXILIARY PARAMETER CONSTRAINT

      IF (KIND .EQ. 0) THEN
        C(1) = 1.D0

*** KIND 1 TO 3 ARE SIMPLE COORDINATE CONSTRAINTS

      ELSEIF (KIND .GE. 1  .AND.  KIND .LE. 3) THEN
        C(1) = 1.D0

*** KIND 4 IS DIFFERENTIAL X (GPS)

      ELSEIF (KIND .EQ. 4) THEN
        CALL GETRI (RI, ISN, B)
        CALL GETRJ (RJ, JSN, B)
        IF ( IAUX .GT. 0  .OR.  IGRT .GT. 0 ) THEN
          IF (L2HLF) THEN
            CALL GETEHX (XI, ISN, B)
            CALL GETEHX (XJ, JSN, B)
          ELSE
            CALL GETECX (XI, ISN, B)
            CALL GETECX (XJ, JSN, B)
          ENDIF
          DELX = XJ - XI
          IF (IGRT .GT. 0) THEN
            IF (L2HLF) THEN
              CALL GETEHY (YI, ISN, B)
              CALL GETEHY (YJ, JSN, B)
              CALL GETEHZ (ZI, ISN, B)
              CALL GETEHZ (ZJ, JSN, B)
            ELSE
              CALL GETECY (YI, ISN, B)
              CALL GETECY (YJ, JSN, B)
              CALL GETECZ (ZI, ISN, B)
              CALL GETECZ (ZJ, JSN, B)
            ENDIF
            DELY = YJ - YI
            DELZ = ZJ - ZI
            CALL GETRT (DELX, DELY, DELZ, RT)
          ENDIF
        ENDIF
        IF (IDIM .EQ. 1) THEN
          C(1) = -RI(3,1)
          C(2) = +RJ(3,1)
          IF (IAUX .GT. 0) THEN
            C(3) = -DELX
            IF (IGRT .GT. 0) THEN
              C(4) = RT(1,1)
              C(5) = RT(1,2)
              C(6) = RT(1,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(3) = RT(1,1)
              C(4) = RT(1,2)
              C(5) = RT(1,3)
            ENDIF
          ENDIF
        ELSEIF (IDIM .EQ. 2  .AND.  ( .NOT. L2HLF) ) THEN
          C(1) = -RI(1,1)
          C(2) = -RI(2,1)
          C(3) = +RJ(1,1)
          C(4) = +RJ(2,1)
          IF (IAUX .GT. 0) THEN
            C(5) = -DELX
            IF (IGRT .GT. 0) THEN
              C(6) = RT(1,1)
              C(7) = RT(1,2)
              C(8) = RT(1,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(5) = RT(1,1)
              C(6) = RT(1,2)
              C(7) = RT(1,3)
            ENDIF
          ENDIF
        ELSE
          C(1) = -RI(1,1)
          C(2) = -RI(2,1)
          C(3) = -RI(3,1)
          C(4) = +RJ(1,1)
          C(5) = +RJ(2,1)
          C(6) = +RJ(3,1)
          IF (IAUX .GT. 0) THEN
            C(7) = -DELX
            IF (IGRT .GT. 0) THEN
              C(8) = RT(1,1)
              C(9) = RT(1,2)
              C(10) = RT(1,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(7) = RT(1,1)
              C(8) = RT(1,2)
              C(9) = RT(1,3)
            ENDIF
          ENDIF
        ENDIF

*** KIND 5 IS DIFFERENTIAL Y (GPS)

      ELSEIF (KIND .EQ. 5) THEN
        CALL GETRI (RI, ISN, B)
        CALL GETRJ (RJ, JSN, B)
        IF ( IAUX .GT. 0  .OR.  IGRT .GT. 0 ) THEN
          IF (L2HLF) THEN
            CALL GETEHY (YI, ISN, B)
            CALL GETEHY (YJ, JSN, B)
          ELSE
            CALL GETECY (YI, ISN, B)
            CALL GETECY (YJ, JSN, B)
          ENDIF
          DELY = YJ - YI
          IF (IGRT .GT. 0) THEN
            IF (L2HLF) THEN
              CALL GETEHX (XI, ISN, B)
              CALL GETEHX (XJ, JSN, B)
              CALL GETEHZ (ZI, ISN, B)
              CALL GETEHZ (ZJ, JSN, B)
            ELSE
              CALL GETECX (XI, ISN, B)
              CALL GETECX (XJ, JSN, B)
              CALL GETECZ (ZI, ISN, B)
              CALL GETECZ (ZJ, JSN, B)
            ENDIF
            DELX = XJ - XI
            DELZ = ZJ - ZI
            CALL GETRT (DELX, DELY, DELZ, RT)
          ENDIF
        ENDIF
        IF (IDIM .EQ. 1) THEN
          C(1) = -RI(3,2)
          C(2) = +RJ(3,2)
          IF (IAUX .GT. 0) THEN
            C(3) = -DELY
            IF (IGRT .GT. 0) THEN
              C(4) = RT(2,1)
              C(5) = RT(2,2)
              C(6) = RT(2,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(3) = RT(2,1)
              C(4) = RT(2,2)
              C(5) = RT(2,3)
            ENDIF
          ENDIF
        ELSEIF (IDIM .EQ. 2  .AND.  ( .NOT. L2HLF) ) THEN
          C(1) = -RI(1,2)
          C(2) = -RI(2,2)
          C(3) = +RJ(1,2)
          C(4) = +RJ(2,2)
          IF (IAUX .GT. 0) THEN
            C(5) = -DELY
            IF (IGRT .GT. 0) THEN
              C(6) = RT(2,1)
              C(7) = RT(2,2)
              C(8) = RT(2,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(5) = RT(2,1)
              C(6) = RT(2,2)
              C(7) = RT(2,3)
            ENDIF
          ENDIF
        ELSE
          C(1) = -RI(1,2)
          C(2) = -RI(2,2)
          C(3) = -RI(3,2)
          C(4) = +RJ(1,2)
          C(5) = +RJ(2,2)
          C(6) = +RJ(3,2)
          IF (IAUX .GT. 0) THEN
            C(7) = -DELY
            IF (IGRT .GT. 0) THEN
              C(8) = RT(2,1)
              C(9) = RT(2,2)
              C(10) = RT(2,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(7) = RT(2,1)
              C(8) = RT(2,2)
              C(9) = RT(2,3)
            ENDIF
          ENDIF
        ENDIF

*** KIND 6 IS DIFFERENTIAL Z (GPS)

      ELSEIF (KIND .EQ. 6) THEN
        CALL GETRI (RI, ISN, B)
        CALL GETRJ (RJ, JSN, B)
        IF ( IAUX .GT. 0  .OR.  IGRT .GT. 0 ) THEN
          IF (L2HLF) THEN
            CALL GETEHZ (ZI, ISN, B)
            CALL GETEHZ (ZJ, JSN, B)
          ELSE
            CALL GETECZ (ZI, ISN, B)
            CALL GETECZ (ZJ, JSN, B)
          ENDIF
          DELZ = ZJ - ZI
          IF (IGRT .GT. 0) THEN
            IF (L2HLF) THEN
              CALL GETEHX (XI, ISN, B)
              CALL GETEHX (XJ, JSN, B)
              CALL GETEHY (YI, ISN, B)
              CALL GETEHY (YJ, JSN, B)
            ELSE
              CALL GETECX (XI, ISN, B)
              CALL GETECX (XJ, JSN, B)
              CALL GETECY (YI, ISN, B)
              CALL GETECY (YJ, JSN, B)
            ENDIF
            DELX = XJ - XI
            DELY = YJ - YI
            CALL GETRT (DELX, DELY, DELZ, RT)
          ENDIF
        ENDIF
        IF (IDIM .EQ. 1) THEN
          C(1) = -RI(3,3)
          C(2) = +RJ(3,3)
          IF (IAUX .GT. 0) THEN
            C(3) = -DELZ
            IF (IGRT .GT. 0) THEN
              C(4) = RT(3,1)
              C(5) = RT(3,2)
              C(6) = RT(3,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(3) = RT(3,1)
              C(4) = RT(3,2)
              C(5) = RT(3,3)
            ENDIF
          ENDIF
        ELSEIF (IDIM .EQ. 2  .AND.  ( .NOT. L2HLF) ) THEN
          C(1) = -RI(1,3)
          C(2) = -RI(2,3)
          C(3) = +RJ(1,3)
          C(4) = +RJ(2,3)
          IF (IAUX .GT. 0) THEN
            C(5) = -DELZ
            IF (IGRT .GT. 0) THEN
              C(6) = RT(3,1)
              C(7) = RT(3,2)
              C(8) = RT(3,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(5) = RT(3,1)
              C(6) = RT(3,2)
              C(7) = RT(3,3)
            ENDIF
          ENDIF
        ELSE
          C(1) = -RI(1,3)
          C(2) = -RI(2,3)
          C(3) = -RI(3,3)
          C(4) = +RJ(1,3)
          C(5) = +RJ(2,3)
          C(6) = +RJ(3,3)
          IF (IAUX .GT. 0) THEN
            C(7) = -DELZ
            IF (IGRT .GT. 0) THEN
              C(8) = RT(3,1)
              C(9) = RT(3,2)
              C(10) = RT(3,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(7) = RT(3,1)
              C(8) = RT(3,2)
              C(9) = RT(3,3)
            ENDIF
          ENDIF
        ENDIF

*** KIND 7 IS MARK TO MARK DISTANCE

      ELSEIF (KIND .EQ. 7) THEN
        CALL GETLAH (ISN, JSN, B, P1, Q1, R1, T1, P2, Q2, R2, T2, S)
        IF (IDIM .EQ. 1) THEN
          C(1) = -T1/S
          C(2) = -T2/S
          IF (IAUX .GT. 0) C(3) = -S
        ELSEIF (IDIM .EQ. 2) THEN
          C(1) = -P1/S
          C(2) = -Q1/S
          C(3) = -P2/S
          C(4) = -Q2/S
          IF (IAUX .GT. 0) C(5) = -S
        ELSE
          C(1) = -P1/S
          C(2) = -Q1/S
          C(3) = -T1/S
          C(4) = -P2/S
          C(5) = -Q2/S
          C(6) = -T2/S
          IF (IAUX .GT. 0) C(7) = -S
        ENDIF

*** KIND 8 IS AN ASTRONOMIC AZIMUTH

      ELSEIF (KIND .EQ. 8) THEN
        CALL GETLAH (ISN, JSN, B, P1, Q1, R1, T1, P2, Q2, R2, T2, S)
        CALL GETGLA (GLAI, ISN, B)
        CALL GETGLO (GLOI, ISN, B)
        CALL GETGLA (GLAJ, JSN, B)
        CALL GETGLO (GLOJ, JSN, B)
        DLAM = GLOJ - GLOI

        SLAI = DSIN(GLAI)
        SLAJ = DSIN(GLAJ)
        CLAI = DCOS(GLAI)
        CLAJ = DCOS(GLAJ)
        SDL = DSIN(DLAM)
        CDL = DCOS(DLAM)
        RR1 = R1*R1

        IF (IDIM .EQ. 1) THEN
          C(1) = 0.D0
          C(2) = 0.D0
        ELSEIF (IDIM .EQ. 2) THEN
          C(1) = +Q1/RR1
          C(2) = -P1/RR1
          C(3) = -(Q1*(CLAI*CLAJ + SLAI*SLAJ*CDL) + P1*SLAJ*SDL)/RR1
          C(4) = +(P1*CDL - Q1*SLAI*SDL)/RR1
        ELSE
          C(1) = +Q1/RR1
          C(2) = -P1/RR1
          C(3) = 0.D0
          C(4) = -(Q1*(CLAI*CLAJ + SLAI*SLAJ*CDL) + P1*SLAJ*SDL)/RR1
          C(5) = +(P1*CDL - Q1*SLAI*SDL)/RR1
          C(6) = 0.D0
        ENDIF
        CALL GETLGH (ISN, JSN, B, C)

*** KIND 9 IS A ZENITH DISTANCE

      ELSEIF (KIND .EQ. 9) THEN
        CALL GETLAH (ISN, JSN, B, P1, Q1, R1, T1, P2, Q2, R2, T2, S)
        CALL GETGLA (GLAI, ISN, B)
        CALL GETGLO (GLOI, ISN, B)
        CALL GETGLA (GLAJ, JSN, B)
        CALL GETGLO (GLOJ, JSN, B)
        DLAM = GLOJ - GLOI

        SLAI = DSIN(GLAI)
        SLAJ = DSIN(GLAJ)
        CLAI = DCOS(GLAI)
        CLAJ = DCOS(GLAJ)
        SDL = DSIN(DLAM)
        CDL = DCOS(DLAM)
        SS2 = S*S

        IF (IDIM .EQ. 1) THEN
          C(1) = +R1/SS2
          C(2) = -(CLAI*CLAJ*CDL + SLAI*SLAJ + (T1*T2)/SS2)/R1
          IF (IAUX .GT. 0) C(3) = -R1
        ELSEIF (IDIM .EQ. 2) THEN
          C(1) = -(P1*T1)/(R1*SS2)
          C(2) = -(Q1*T1)/(R1*SS2)
          C(3) = -(-CLAI*SLAJ*CDL + SLAI*CLAJ + (T1*P2)/SS2)/R1
          C(4) = -(-CLAI*SDL + (T1*Q2)/SS2)/R1
          IF (IAUX .GT. 0) THEN
            CALL GETAUX (VAL, IAUX, B)
            RP = VAL*P1/R1
            RQ = VAL*Q1/R1
            C(1) = C(1) + RP
            C(2) = C(2) + RQ
            C(3) = C(3) - RP
            C(4) = C(4) - RQ
            C(5) = -R1
          ENDIF
        ELSE
          C(1) = -(P1*T1)/(R1*SS2)
          C(2) = -(Q1*T1)/(R1*SS2)
          C(3) = +R1/SS2
          C(4) = -(-CLAI*SLAJ*CDL + SLAI*CLAJ + (T1*P2)/SS2)/R1
          C(5) = -(-CLAI*SDL + (T1*Q2)/SS2)/R1
          C(6) = -(CLAI*CLAJ*CDL + SLAI*SLAJ + (T1*T2)/SS2)/R1
          IF (IAUX .GT. 0) THEN
            CALL GETAUX (VAL, IAUX, B)
            RP = VAL*P1/R1
            RQ = VAL*Q1/R1
            C(1) = C(1) + RP
            C(2) = C(2) + RQ
            C(4) = C(4) - RP
            C(5) = C(5) - RQ
            C(7) = -R1
          ENDIF
        ENDIF
        CALL GETLGH (ISN, JSN, B, C)

*** KIND 10 IS HORIZONTAL ANGLE (MARK-TO-MARK)

      ELSEIF (KIND .EQ. 10) THEN
        KSN = IAUX
        CALL GETLAH (ISN, JSN, B, P1L, Q1L, R1L, T1L, P2, Q2, R2, T2,S2)
        CALL GETLAH (ISN, KSN, B, P1R, Q1R, R1R, T1R, P3, Q3, R3, T3,S3)
        CALL GETGLA (GLAT1, ISN, B)
        CALL GETGLO (GLON1, ISN, B)
        CALL GETGLA (GLAT2, JSN, B)
        CALL GETGLO (GLON2, JSN, B)
        CALL GETGLA (GLAT3, KSN, B)
        CALL GETGLO (GLON3, KSN, B)
        SLAT1 = DSIN(GLAT1)
        SLAT2 = DSIN(GLAT2)
        SLAT3 = DSIN(GLAT3)
        CLAT1 = DCOS(GLAT1)
        CLAT2 = DCOS(GLAT2)
        CLAT3 = DCOS(GLAT3)
        SLONXL = DSIN(GLON2 - GLON1)
        CLONXL = DCOS(GLON2 - GLON1)
        SLONXR = DSIN(GLON3 - GLON1)
        CLONXR = DCOS(GLON3 - GLON1)
        DIVR1L = 1.D0/(R1L*R1L)
        DIVR1R = 1.D0/(R1R*R1R)
        IF (IDIM .EQ. 1) THEN
          C(1) = 0.D0
          RETURN
        ELSEIF (IDIM .EQ. 2) THEN
          CL(1) = Q1L*DIVR1L
          CL(2) = -P1L*DIVR1L
          CL(3) = -(Q1L*(CLAT1*CLAT2 + SLAT1*SLAT2*CLONXL) +
     &             P1L*SLAT2*SLONXL)*DIVR1L
          CL(4) = (P1L*CLONXL - Q1L*SLAT1*SLONXL)*DIVR1L
          CALL GETLGH (ISN, JSN, B, CL)
          CR(1) = Q1R*DIVR1R
          CR(2) = -P1R*DIVR1R
          CR(3) = -(Q1R*(CLAT1*CLAT3 + SLAT1*SLAT3*CLONXR) +
     &             P1R*SLAT3*SLONXR)*DIVR1R
          CR(4) = (P1R*CLONXR - Q1R*SLAT1*SLONXR)*DIVR1R
          CALL GETLGH (ISN, KSN, B, CR)
        ELSE
          CL(1) = Q1L*DIVR1L
          CL(2) = -P1L*DIVR1L
          CL(3) = 0.D0
          CL(4) = -(Q1L*(CLAT1*CLAT2 + SLAT1*SLAT2*CLONXL) +
     &             P1L*SLAT2*SLONXL)*DIVR1L
          CL(5) = (P1L*CLONXL - Q1L*SLAT1*SLONXL)*DIVR1L
          CL(6) = 0.D0
          CALL GETLGH (ISN, JSN, B, CL)
          CR(1) = Q1R*DIVR1R
          CR(2) = -P1R*DIVR1R
          CR(3) = 0.D0
          CR(4) = -(Q1R*(CLAT1*CLAT3 + SLAT1*SLAT3*CLONXR) +
     &             P1R*SLAT3*SLONXR)*DIVR1R
          CR(5) = (P1R*CLONXR - Q1R*SLAT1*SLONXR)*DIVR1R
          CR(6) = 0.D0
          CALL GETLGH (ISN, KSN, B, CR)
        ENDIF
        IF (IDIM .EQ. 2) THEN
          C(1) = CR(1) - CL(1)
          C(2) = CR(2) - CL(2)
          C(3) = -CL(3)
          C(4) = -CL(4)
          C(5) = CR(3)
          C(6) = CR(4)
        ELSE
          C(1) = CR(1) - CL(1)
          C(2) = CR(2) - CL(2)
          C(3) = -CL(4)
          C(4) = -CL(5)
          C(5) = CR(4)
          C(6) = CR(5)
      ENDIF

*** KIND 11 IS A HORIZONTAL DIRECTION

      ELSEIF (KIND .EQ. 11) THEN
        CALL GETLAH (ISN, JSN, B, P1, Q1, R1, T1, P2, Q2, R2, T2, S)
        CALL GETGLA (GLAI, ISN, B)
        CALL GETGLO (GLOI, ISN, B)
        CALL GETGLA (GLAJ, JSN, B)
        CALL GETGLO (GLOJ, JSN, B)
        DLAM = GLOJ - GLOI

        SLAI = DSIN(GLAI)
        SLAJ = DSIN(GLAJ)
        CLAI = DCOS(GLAI)
        CLAJ = DCOS(GLAJ)
        SDL = DSIN(DLAM)
        CDL = DCOS(DLAM)
        RR1 = R1*R1

        IF (IDIM .EQ. 1) THEN
          C(1) = 0.D0
          C(2) = 0.D0
          C(3) = -1.D0
        ELSEIF (IDIM .EQ. 2) THEN
          C(1) = +Q1/RR1
          C(2) = -P1/RR1
          C(3) = -(Q1*(CLAI*CLAJ + SLAI*SLAJ*CDL) + P1*SLAJ*SDL)/RR1
          C(4) = +(P1*CDL - Q1*SLAI*SDL)/RR1
          C(5) = -1.D0
        ELSE
          C(1) = +Q1/RR1
          C(2) = -P1/RR1
          C(3) = 0.D0
          C(4) = -(Q1*(CLAI*CLAJ + SLAI*SLAJ*CDL) + P1*SLAJ*SDL)/RR1
          C(5) = +(P1*CDL - Q1*SLAI*SDL)/RR1
          C(6) = 0.D0
          C(7) = -1.D0
        ENDIF
        CALL GETLGH (ISN, JSN, B, C)

*** KIND 12 IS A CONSTRAINED AZIMUTH

      ELSEIF (KIND .EQ. 12) THEN
        CALL GETLAH (ISN, JSN, B, P1, Q1, R1, T1, P2, Q2, R2, T2, S)
        CALL GETGLA (GLAI, ISN, B)
        CALL GETGLO (GLOI, ISN, B)
        CALL GETGLA (GLAJ, JSN, B)
        CALL GETGLO (GLOJ, JSN, B)
        DLAM = GLOJ - GLOI

        SLAI = DSIN(GLAI)
        SLAJ = DSIN(GLAJ)
        CLAI = DCOS(GLAI)
        CLAJ = DCOS(GLAJ)
        SDL = DSIN(DLAM)
        CDL = DCOS(DLAM)
        RR1 = R1*R1

        IF (IDIM .EQ. 1) THEN
          C(1) = 0.D0
          C(2) = 0.D0
        ELSEIF (IDIM .EQ. 2) THEN
          C(1) = +Q1/RR1
          C(2) = -P1/RR1
          C(3) = -(Q1*(CLAI*CLAJ + SLAI*SLAJ*CDL) + P1*SLAJ*SDL)/RR1
          C(4) = +(P1*CDL - Q1*SLAI*SDL)/RR1
        ELSE
          C(1) = +Q1/RR1
          C(2) = -P1/RR1
          C(3) = 0.D0
          C(4) = -(Q1*(CLAI*CLAJ + SLAI*SLAJ*CDL) + P1*SLAJ*SDL)/RR1
          C(5) = +(P1*CDL - Q1*SLAI*SDL)/RR1
          C(6) = 0.D0
        ENDIF
        CALL GETLGH (ISN, JSN, B, C)

*** KIND 13 IS A CONSTRAINED MARK TO MARK DISTANCE

      ELSEIF (KIND .EQ. 13) THEN
        CALL GETLAH (ISN, JSN, B, P1, Q1, R1, T1, P2, Q2, R2, T2, S)
        IF (IDIM .EQ. 1) THEN
          C(1) = -T1/S
          C(2) = -T2/S
        ELSEIF (IDIM .EQ. 2) THEN
          C(1) = -P1/S
          C(2) = -Q1/S
          C(3) = -P2/S
          C(4) = -Q2/S
        ELSE
          C(1) = -P1/S
          C(2) = -Q1/S
          C(3) = -T1/S
          C(4) = -P2/S
          C(5) = -Q2/S
          C(6) = -T2/S
        ENDIF

*** KIND 14 IS A CONSTRAINED ZENITH DISTANCE

      ELSEIF (KIND .EQ. 14) THEN
        CALL GETLAH (ISN, JSN, B, P1, Q1, R1, T1, P2, Q2, R2, T2, S)
        CALL GETGLA (GLAI, ISN, B)
        CALL GETGLO (GLOI, ISN, B)
        CALL GETGLA (GLAJ, JSN, B)
        CALL GETGLO (GLOJ, JSN, B)
        DLAM = GLOJ - GLOI

        SLAI = DSIN(GLAI)
        SLAJ = DSIN(GLAJ)
        CLAI = DCOS(GLAI)
        CLAJ = DCOS(GLAJ)
        SDL = DSIN(DLAM)
        CDL = DCOS(DLAM)
        SS2 = S*S

        IF (IDIM .EQ. 1) THEN
          C(1) = +R1/SS2
          C(2) = -(CLAI*CLAJ*CDL + SLAI*SLAJ + (T1*T2)/SS2)/R1
        ELSEIF (  IDIM .EQ. 2  .AND.
     &            ( I2HLF(ISN) .NE. 1  .OR.  I2HLF(JSN) .NE. 1 )  ) THEN
          C(1) = -(P1*T1)/(R1*SS2)
          C(2) = -(Q1*T1)/(R1*SS2)
          C(3) = -(-CLAI*SLAJ*CDL + SLAI*CLAJ + (T1*P2)/SS2)/R1
          C(4) = -(-CLAI*SDL + (T1*Q2)/SS2)/R1
        ELSE
          C(1) = -(P1*T1)/(R1*SS2)
          C(2) = -(Q1*T1)/(R1*SS2)
          C(3) = +R1/SS2
          C(4) = -(-CLAI*SLAJ*CDL + SLAI*CLAJ + (T1*P2)/SS2)/R1
          C(5) = -(-CLAI*SDL + (T1*Q2)/SS2)/R1
          C(6) = -(CLAI*CLAJ*CDL + SLAI*SLAJ + (T1*T2)/SS2)/R1
        ENDIF
        CALL GETLGH (ISN, JSN, B, C)

*** KIND 15 THRU 17 ARE CONSTRAINED HEIGHT DIFFERENCES

      ELSEIF (KIND .GE. 15  .AND.  KIND .LE. 17) THEN
        IF (  IDIM .EQ. 1  .OR.  IDIM .EQ. 3  .OR.
     &       (IDIM .EQ. 2  .AND.
     &        I2HLF(ISN) .EQ. 1  .AND.  I2HLF(JSN) .EQ. 1 )
     &     ) THEN
          C(1) = -1.D0
          C(2) = +1.D0
        ELSE
          C(1) = 0.D0
        ENDIF

*** KIND 18 IS E.C.F. X (DOPPLER)

      ELSEIF (KIND .EQ. 18) THEN
        CALL GETRI (RI, ISN, B)
        IF ( IAUX .GT. 0  .OR.  IGRT .GT. 0 ) THEN
          IF (L2HLF) THEN
            CALL GETEHX (XI, ISN, B)
          ELSE
            CALL GETECX (XI, ISN, B)
          ENDIF
          DELX = +XI
          IF (IGRT .GT. 0) THEN
            IF (L2HLF) THEN
              CALL GETEHY (YI, ISN, B)
              CALL GETEHZ (ZI, ISN, B)
            ELSE
              CALL GETECY (YI, ISN, B)
              CALL GETECZ (ZI, ISN, B)
            ENDIF
            DELY = +YI
            DELZ = +ZI
            CALL GETRT (DELX, DELY, DELZ, RT)
          ENDIF
        ENDIF
        IF (IDIM .EQ. 1) THEN
          C(1) = +RI(3,1)
          IF (IAUX .GT. 0) THEN
            C(2) = -DELX
            IF (IGRT .GT. 0) THEN
              C(3) = RT(1,1)
              C(4) = RT(1,2)
              C(5) = RT(1,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(2) = RT(1,1)
              C(3) = RT(1,2)
              C(4) = RT(1,3)
            ENDIF
          ENDIF
        ELSEIF (IDIM .EQ. 2  .AND.  ( .NOT. L2HLF) ) THEN
          C(1) = +RI(1,1)
          C(2) = +RI(2,1)
          IF (IAUX .GT. 0) THEN
            C(3) = -DELX
            IF (IGRT .GT. 0) THEN
              C(4) = RT(1,1)
              C(5) = RT(1,2)
              C(6) = RT(1,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(3) = RT(1,1)
              C(4) = RT(1,2)
              C(5) = RT(1,3)
            ENDIF
          ENDIF
        ELSE
          C(1) = +RI(1,1)
          C(2) = +RI(2,1)
          C(3) = +RI(3,1)
          IF (IAUX .GT. 0) THEN
            C(4) = -DELX
            IF (IGRT .GT. 0) THEN
              C(5) = RT(1,1)
              C(6) = RT(1,2)
              C(7) = RT(1,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(4) = RT(1,1)
              C(5) = RT(1,2)
              C(6) = RT(1,3)
            ENDIF
          ENDIF
        ENDIF

*** KIND 19 IS E.C.F. Y (DOPPLER)

      ELSEIF (KIND .EQ. 19) THEN
        CALL GETRI (RI, ISN, B)
        IF ( IAUX .GT. 0  .OR.  IGRT .GT. 0 ) THEN
          IF (L2HLF) THEN
            CALL GETEHY (YI, ISN, B)
          ELSE
            CALL GETECY (YI, ISN, B)
          ENDIF
          DELY = +YI
          IF (IGRT .GT. 0) THEN
            IF (L2HLF) THEN
              CALL GETEHX (XI, ISN, B)
              CALL GETEHZ (ZI, ISN, B)
            ELSE
              CALL GETECX (XI, ISN, B)
              CALL GETECZ (ZI, ISN, B)
            ENDIF
            DELX = +XI
            DELZ = +ZI
            CALL GETRT (DELX, DELY, DELZ, RT)
          ENDIF
        ENDIF
        IF (IDIM .EQ. 1) THEN
          C(1) = +RI(3,2)
          IF (IAUX .GT. 0) THEN
            C(2) = -DELY
            IF (IGRT .GT. 0) THEN
              C(3) = RT(2,1)
              C(4) = RT(2,2)
              C(5) = RT(2,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(2) = RT(2,1)
              C(3) = RT(2,2)
              C(4) = RT(2,3)
            ENDIF
          ENDIF
        ELSEIF (IDIM .EQ. 2  .AND.  ( .NOT. L2HLF) ) THEN
          C(1) = +RI(1,2)
          C(2) = +RI(2,2)
          IF (IAUX .GT. 0) THEN
            C(3) = -DELY
            IF (IGRT .GT. 0) THEN
              C(4) = RT(2,1)
              C(5) = RT(2,2)
              C(6) = RT(2,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(2) = RT(2,1)
              C(3) = RT(2,2)
              C(4) = RT(2,3)
            ENDIF
          ENDIF
        ELSE
          C(1) = +RI(1,2)
          C(2) = +RI(2,2)
          C(3) = +RI(3,2)
          IF (IAUX .GT. 0) THEN
            C(4) = -DELY
            IF (IGRT .GT. 0) THEN
              C(5) = RT(2,1)
              C(6) = RT(2,2)
              C(7) = RT(2,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(4) = RT(2,1)
              C(5) = RT(2,2)
              C(6) = RT(2,3)
            ENDIF
          ENDIF
        ENDIF

*** KIND 20 IS E.C.F. Z (DOPPLER)

      ELSEIF (KIND .EQ. 20) THEN
        CALL GETRI (RI, ISN, B)
        IF ( IAUX .GT. 0  .OR.  IGRT .GT. 0 ) THEN
          IF (L2HLF) THEN
            CALL GETEHZ (ZI, ISN, B)
          ELSE
            CALL GETECZ (ZI, ISN, B)
          ENDIF
          DELZ = +ZI
          IF (IGRT .GT. 0) THEN
            IF (L2HLF) THEN
              CALL GETEHX (XI, ISN, B)
              CALL GETEHY (YI, ISN, B)
            ELSE
              CALL GETECX (XI, ISN, B)
              CALL GETECY (YI, ISN, B)
            ENDIF
            DELX = +XI
            DELY = +YI
            CALL GETRT (DELX, DELY, DELZ, RT)
          ENDIF
        ENDIF
        IF (IDIM .EQ. 1) THEN
          C(1) = +RI(3,3)
          IF (IAUX .GT. 0) THEN
            C(2) = -DELZ
            IF (IGRT .GT. 0) THEN
              C(3) = RT(3,1)
              C(4) = RT(3,2)
              C(5) = RT(3,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(2) = RT(3,1)
              C(3) = RT(3,2)
              C(4) = RT(3,3)
            ENDIF
          ENDIF
        ELSEIF (IDIM .EQ. 2  .AND.  ( .NOT. L2HLF) ) THEN
          C(1) = +RI(1,3)
          C(2) = +RI(2,3)
          IF (IAUX .GT. 0) THEN
            C(3) = -DELZ
            IF (IGRT .GT. 0) THEN
              C(4) = RT(3,1)
              C(5) = RT(3,2)
              C(6) = RT(3,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(3) = RT(3,1)
              C(4) = RT(3,2)
              C(5) = RT(3,3)
            ENDIF
          ENDIF
        ELSE
          C(1) = +RI(1,3)
          C(2) = +RI(2,3)
          C(3) = +RI(3,3)
          IF (IAUX .GT. 0) THEN
            C(4) = -DELZ
            IF (IGRT .GT. 0) THEN
              C(5) = RT(3,1)
              C(6) = RT(3,2)
              C(7) = RT(3,3)
            ENDIF
          ELSE
            IF (IGRT .GT. 0) THEN
              C(4) = RT(3,1)
              C(5) = RT(3,2)
              C(6) = RT(3,3)
            ENDIF
          ENDIF
        ENDIF

*** KIND 21 IS A CONSTRAINED NORTH COORD DIFFERENCE
*** KIND 22 IS A CONSTRAINED EAST COORD DIFFERENCE

      ELSEIF(KIND.EQ.21.OR.KIND.EQ.22) THEN
        IF(IDIM.EQ.2.OR.IDIM.EQ.3) THEN
          C(1) = -1.D0
          C(2) = +1.D0
        ELSE
          C(1) = 0.D0
        ENDIF

      ELSE
        WRITE (LUNIT,1) KIND
    1   FORMAT ('0ILLEGAL KIND IN FORMC =', I5)
        CALL ABORT2
      ENDIF

      RETURN
      END
      SUBROUTINE FORMG (IUO, IUO2, NVEC, NR, G, B, FATAL, IVF,IAUX,IGRT)

*** ROUTINE TO RE-COMPUTE OBS.EQ. FOR A GPS CLUSTER

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( NVECS = 700, MXSSN = 9999 )
      DIMENSION ICM((NVECS+1)*3+1+3)
      DIMENSION KINDS(NVECS*3), ISNS(NVECS*3), JSNS(NVECS*3)
      DIMENSION LOBS(NVECS*3)
      LOGICAL FATAL, NOBIGV
      LOGICAL L2HLF, LEHT
      DIMENSION B(*), G(NR,*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

*** ABORT/PRINT BIG RESIDUALS

      NOBIGV = .FALSE.

      IF ( .NOT. L2HLF) THEN
        LENG = (NVEC+1)*IDIM + 1
      ELSE
        LENG = (NVEC+1)*3 + 1
      ENDIF
      IF (NGRT .GT. 0) LENG = LENG + 3
*     NC = NR + 3 + LENG
      NC = NR + 5 + LENG


*** READ SUPPORTING INDICIES

      READ (IUO,END=666) ICM, NICM, KINDS, ISNS, JSNS, LOBS
      WRITE (IUO2) ICM, NICM, KINDS, ISNS, JSNS, LOBS

*** LOAD WORK SPACE (G)

      DO 1 I = 1, NR
        READ (IUO,END=666) (G(I,J), J = 1, NC)
    1 CONTINUE

*** COMPUTE AND DECORRELATE OBS. EQ. AND MISCLOSURE

      CALL GOBSEQ (G, NR, NC, NVEC, LENG, B, ICM, NICM, KINDS, ISNS,
     &             JSNS, LOBS, IAUX, IVF, FATAL, NOBIGV, N4, IGRT)

*** UNLOAD WORK SPACE

      DO 500 I = 1, NR
        WRITE (IUO2) (G(I,J), J = 1, NC)
  500 CONTINUE

      RETURN

*** PREMATURE END OF FILE

  666 WRITE (LUNIT,667) NVEC
  667 FORMAT ('0PREMATURE FILE END IN FORMG -- NVEC = ',I5)
      CALL ABORT2
      RETURN
      END
      SUBROUTINE FORMG2 (IUO, IUO2, NVEC, NR, G, B, A, NX, FATAL, IVF,
     &                   IAUX, IGRT)

*** ROUTINE TO RE-COMPUTE OBS.EQ. FOR A GPS CLUSTER

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( NVECS = 700, MXSSN = 9999 )
      DIMENSION ICM((NVECS+1)*3+1+3)
      DIMENSION KINDS(NVECS*3), ISNS(NVECS*3), JSNS(NVECS*3)
      DIMENSION LOBS(NVECS*3)
      LOGICAL FATAL
      LOGICAL L2HLF, LEHT
      DIMENSION B(*), G(NR,*), A (*), NX(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT

      IF ( .NOT. L2HLF) THEN
        LENG = (NVEC+1)*IDIM + 1
      ELSE
        LENG = (NVEC+1)*3 + 1
      ENDIF
      IF (NGRT .GT. 0) LENG = LENG + 3
*     NC = NR + 3 + LENG
      NC = NR + 5 + LENG

*** READ SUPPORTING INDICIES

      READ (IUO) ICM, NICM, KINDS, ISNS, JSNS, LOBS
      WRITE (IUO2) ICM, NICM, KINDS, ISNS, JSNS, LOBS

*** LOAD WORK SPACE (G)

      DO 1 I = 1, NR
        READ (IUO) (G(I,J), J = 1, NC)
    1 CONTINUE

*** COMPUTE AND DECORRELATE OBS. EQ. AND MISCLOSURE

      CALL VFGPS (G, NR, NC, NVEC, LENG, B, ICM, NICM, KINDS, ISNS,
     &            JSNS, LOBS, IAUX, IVF, FATAL, A, NX, N2, N3, N4, IGRT)

*** UNLOAD WORK SPACE

      DO 500 I = 1, NR
        WRITE (IUO2) (G(I,J), J = 1, NC)
  500 CONTINUE

      RETURN
      END
      SUBROUTINE FORMIC (KIND, I, J, IAUX, IC, LENG, IGRT)

*** LOAD UNKNOWN NUMBER VECTOR WITH UNKNOWN NUMBERS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999, LENC = 10 )
      LOGICAL L2HLF, LEHT
      DIMENSION IC(LENC)
      DIMENSION LENGS3(22), LENGS2(22), LENGS1(22)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

      DATA LENGS3 /1,1,1,6,6,6,6,6,6,6,7,6,6,6,2,2,2,3,3,3,2,2/
      DATA LENGS2 /1,1,1,4,4,4,4,4,4,6,5,4,4,4,1,1,1,2,2,2,2,2/
      DATA LENGS1 /1,1,1,2,2,2,2,2,2,1,3,2,2,2,2,2,2,1,1,1,1,1/

*** KIND 0 IS AN AUXILIARY PARAMETER CONSTRAINT (I = IAUX)

      IF (KIND .EQ. 0) THEN
        IC(1) = IUNAUX(I)
        LENG = 1

*** KIND 1 TO 3 ARE SIMPLE COORDINATE CONSTRAINTS (I = ISTA)

      ELSEIF ( KIND .GE. 1  .AND.  KIND .LE. 3  .AND.
     &            I .GE. 1  .AND.     I .LE. NSTA )THEN
        IC(1) = IUNSTA(I, KIND)
        LENG = LENGS3(KIND)

*** KIND 4 THRU 6 IS DIFFERENTIAL ECF X,Y,Z (GPS)

      ELSEIF ( KIND .GE. 4  .AND.  KIND .LE. 6  .AND.
     &            I .GE. 1  .AND.  I .LE. NSTA  .AND.
     &            J .GE. 1  .AND.  J .LE. NSTA) THEN
        IF (IDIM .EQ. 1) THEN
          LENG = LENGS1(KIND)
          IC(1) = IUNSTA(I, 3)
          IC(2) = IUNSTA(J, 3)
        ELSEIF (IDIM .EQ. 2  .AND.  (.NOT.L2HLF) ) THEN
          LENG = LENGS2(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(J, 1)
          IC(4) = IUNSTA(J, 2)
        ELSE
          LENG = LENGS3(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(I, 3)
          IC(4) = IUNSTA(J, 1)
          IC(5) = IUNSTA(J, 2)
          IC(6) = IUNSTA(J, 3)
        ENDIF
        IF (IAUX .GT. 0) THEN
          IF (IGRT .GT. 0) THEN
            LENG = LENG + 4
            IC(LENG-3) = IUNAUX(IAUX)
            IC(LENG-2) = IUNGRT(IGRT, 1)
            IC(LENG-1) = IUNGRT(IGRT, 2)
            IC(LENG) = IUNGRT(IGRT, 3)
          ELSE
            LENG = LENG + 1
            IC(LENG) = IUNAUX(IAUX)
          ENDIF
        ELSE
          IF (IGRT .GT. 0) THEN
            LENG = LENG + 3
            IC(LENG-2) = IUNGRT(IGRT, 1)
            IC(LENG-1) = IUNGRT(IGRT, 2)
            IC(LENG) = IUNGRT(IGRT, 3)
          ENDIF
        ENDIF

*** KIND 7 IS MARK TO MARK DISTANCE

      ELSEIF ( KIND .EQ. 7  .AND.  I .GE. 1  .AND.  I .LE. NSTA  .AND.
     &                             J .GE. 1  .AND.  J .LE. NSTA) THEN
        IF (IDIM .EQ. 1) THEN
          LENG = LENGS1(KIND)
          IC(1) = IUNSTA(I, 3)
          IC(2) = IUNSTA(J, 3)
        ELSEIF (IDIM .EQ. 2) THEN
          LENG = LENGS2(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(J, 1)
          IC(4) = IUNSTA(J, 2)
        ELSE
          LENG = LENGS3(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(I, 3)
          IC(4) = IUNSTA(J, 1)
          IC(5) = IUNSTA(J, 2)
          IC(6) = IUNSTA(J, 3)
        ENDIF
        IF (IAUX .GT. 0) THEN
          LENG = LENG + 1
          IC(LENG) = IUNAUX(IAUX)
        ENDIF

***  KIND 8 IS AN ASTRONOMIC AZIMUTH

      ELSEIF (KIND .EQ. 8  .AND.  I .GE. 1  .AND.  I .LE. NSTA  .AND.
     &                            J .GE. 1  .AND.  J .LE. NSTA) THEN
        IF (IDIM .EQ. 1) THEN
          LENG = LENGS1(KIND)
          IC(1) = IUNSTA(I, 3)
          IC(2) = IUNSTA(J, 3)
        ELSEIF (IDIM .EQ. 2) THEN
          LENG = LENGS2(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(J, 1)
          IC(4) = IUNSTA(J, 2)
        ELSE
          LENG = LENGS3(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(I, 3)
          IC(4) = IUNSTA(J, 1)
          IC(5) = IUNSTA(J, 2)
          IC(6) = IUNSTA(J, 3)
        ENDIF

***  KIND 9 IS A ZENITH DISTANCE

      ELSEIF (KIND .EQ. 9  .AND.  I .GE. 1  .AND.  I .LE. NSTA  .AND.
     &                            J .GE. 1  .AND.  J .LE. NSTA) THEN
        IF (IDIM .EQ. 1) THEN
          LENG = LENGS1(KIND)
          IC(1) = IUNSTA(I, 3)
          IC(2) = IUNSTA(J, 3)
        ELSEIF (IDIM .EQ. 2) THEN
          LENG = LENGS2(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(J, 1)
          IC(4) = IUNSTA(J, 2)
        ELSE
          LENG = LENGS3(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(I, 3)
          IC(4) = IUNSTA(J, 1)
          IC(5) = IUNSTA(J, 2)
          IC(6) = IUNSTA(J, 3)
        ENDIF
        IF (IAUX .GT. 0) THEN
          LENG = LENG + 1
          IC(LENG) = IUNAUX(IAUX)
        ENDIF

*** KIND 10 IS HORIZONTAL ANGLE

      ELSEIF ( KIND.EQ.10 .AND.    I.GE.1  .AND.     I.LE.NSTA  .AND.
     &                             J.GE.1  .AND.     J.LE.NSTA  .AND.
     &                          IAUX.GE.1  .AND.  IAUX.LE.NSTA) THEN
        K = IAUX
        IF (IDIM .EQ. 1) THEN
          LENG = LENGS1(KIND)
          IC(1) = IUNSTA(I, 3)
        ELSE
          LENG = LENGS2(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(J, 1)
          IC(4) = IUNSTA(J, 2)
          IC(5) = IUNSTA(K, 1)
          IC(6) = IUNSTA(K, 2)
        ENDIF

***  KIND 11 IS A HORIZONTAL DIRECTION

      ELSEIF ( KIND .EQ. 11  .AND.  I .GE. 1  .AND.  I .LE. NSTA  .AND.
     &                              J .GE. 1  .AND.  J .LE. NSTA) THEN
        IF (IDIM .EQ. 1) THEN
          LENG = LENGS1(KIND)
          IC(1) = IUNSTA(I, 3)
          IC(2) = IUNSTA(J, 3)
          IC(3) = IUNROT(IAUX)
        ELSEIF (IDIM .EQ. 2) THEN
          LENG = LENGS2(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(J, 1)
          IC(4) = IUNSTA(J, 2)
          IC(5) = IUNROT(IAUX)
        ELSE
          LENG = LENGS3(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(I, 3)
          IC(4) = IUNSTA(J, 1)
          IC(5) = IUNSTA(J, 2)
          IC(6) = IUNSTA(J, 3)
          IC(7) = IUNROT(IAUX)
        ENDIF

***  KIND 12 IS A CONSTRAINED AZIMUTH

      ELSEIF (KIND .EQ. 12  .AND.  I .GE. 1  .AND.  I .LE. NSTA  .AND.
     &                             J .GE. 1  .AND.  J .LE. NSTA) THEN
        IF (IDIM .EQ. 1) THEN
          LENG = LENGS1(KIND)
          IC(1) = IUNSTA(I, 3)
          IC(2) = IUNSTA(J, 3)
        ELSEIF (IDIM .EQ. 2) THEN
          LENG = LENGS2(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(J, 1)
          IC(4) = IUNSTA(J, 2)
        ELSE
          LENG = LENGS3(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(I, 3)
          IC(4) = IUNSTA(J, 1)
          IC(5) = IUNSTA(J, 2)
          IC(6) = IUNSTA(J, 3)
        ENDIF

*** KIND 13 IS A CONSTRAINED MARK TO MARK DISTANCE

      ELSEIF (KIND .EQ. 13  .AND.  I .GE. 1  .AND.  I .LE. NSTA  .AND.
     &                             J .GE. 1  .AND.  J .LE. NSTA) THEN
        IF (IDIM .EQ. 1) THEN
          LENG = LENGS1(KIND)
          IC(1) = IUNSTA(I, 3)
          IC(2) = IUNSTA(J, 3)
        ELSEIF (IDIM .EQ. 2) THEN
          LENG = LENGS2(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(J, 1)
          IC(4) = IUNSTA(J, 2)
        ELSE
          LENG = LENGS3(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(I, 3)
          IC(4) = IUNSTA(J, 1)
          IC(5) = IUNSTA(J, 2)
          IC(6) = IUNSTA(J, 3)
        ENDIF

***  KIND 14 IS A CONSTRAINED ZENITH DISTANCE

      ELSEIF (KIND .EQ. 14  .AND.  I .GE. 1  .AND.  I .LE. NSTA  .AND.
     &                             J .GE. 1  .AND.  J .LE. NSTA) THEN
        IF (IDIM .EQ. 1) THEN
          LENG = LENGS1(KIND)
          IC(1) = IUNSTA(I, 3)
          IC(2) = IUNSTA(J, 3)
        ELSEIF (  IDIM .EQ. 2  .AND.
     &            ( I2HLF(I) .NE. 1  .OR.  I2HLF(J) .NE. 1 )  ) THEN
          LENG = LENGS2(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(J, 1)
          IC(4) = IUNSTA(J, 2)
        ELSE
          LENG = LENGS3(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(I, 3)
          IC(4) = IUNSTA(J, 1)
          IC(5) = IUNSTA(J, 2)
          IC(6) = IUNSTA(J, 3)
        ENDIF

*** KIND 15 THRU 17 ARE CONSTRAINED HEIGHT DIFFERENCES

      ELSEIF (KIND .GE. 15  .AND.  KIND .LE. 17    .AND.
     &           I .GE. 1   .AND.     I .LE. NSTA  .AND.
     &           J .GE. 1   .AND.     J .LE. NSTA) THEN
        IF (  IDIM .EQ. 1  .OR.  IDIM .EQ. 3  .OR.
     &       (IDIM .EQ. 2  .AND.
     &        I2HLF(I) .EQ. 1  .AND.  I2HLF(J) .EQ. 1 )
     &     ) THEN
          LENG = LENGS1(KIND)
          IC(1) = IUNSTA(I, 3)
          IC(2) = IUNSTA(J, 3)
        ELSE
          LENG = LENGS2(KIND)
          IC(1) = IUNSTA(I, 1)
        ENDIF

*** KIND 18 THRU 20 ARE E.C.F. XYZ OBSERVATIONS (DOPPLER)

      ELSEIF ( KIND .GE. 18  .AND.  KIND .LE. 20  .AND.
     &            I .GE. 1   .AND.  I .LE. NSTA ) THEN
        IF (IDIM .EQ. 1) THEN
          LENG = LENGS1(KIND)
          IC(1) = IUNSTA(I, 3)
        ELSEIF (IDIM .EQ. 2  .AND.  ( .NOT. L2HLF) ) THEN
          LENG = LENGS2(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
        ELSE
          LENG = LENGS3(KIND)
          IC(1) = IUNSTA(I, 1)
          IC(2) = IUNSTA(I, 2)
          IC(3) = IUNSTA(I, 3)
        ENDIF
        IF (IAUX .GT. 0) THEN
          IF (IGRT .GT. 0) THEN
            LENG = LENG + 4
            IC(LENG-3) = IUNAUX(IAUX)
            IC(LENG-2) = IUNGRT(IGRT, 1)
            IC(LENG-1) = IUNGRT(IGRT, 2)
            IC(LENG) = IUNGRT(IGRT, 3)
          ELSE
            LENG = LENG + 1
            IC(LENG) = IUNAUX(IAUX)
          ENDIF
        ELSE
          IF (IGRT .GT. 0) THEN
            LENG = LENG + 3
            IC(LENG-2) = IUNGRT(IGRT, 1)
            IC(LENG-1) = IUNGRT(IGRT, 2)
            IC(LENG) = IUNGRT(IGRT, 3)
          ENDIF
        ENDIF

*** KIND 21 IS A CONSTRAINED NORTH COORD DIFFERENCE

      ELSEIF(KIND.EQ.21) THEN
        IF(IDIM.EQ.2.OR.IDIM.EQ.3) THEN
          LENG  = LENGS3(KIND)
          IC(1) = IUNSTA(I,1)
          IC(2) = IUNSTA(J,1)
        ELSE
          LENG  = LENGS1(KIND)
          IC(1) = IUNSTA(I,3)
        ENDIF

*** KIND 22 IS A CONSTRAINED EAST COORD DIFFERENCE

      ELSEIF(KIND.EQ.22) THEN
        IF(IDIM.EQ.2.OR.IDIM.EQ.3) THEN
          LENG  = LENGS3(KIND)
          IC(1) = IUNSTA(I,2)
          IC(2) = IUNSTA(J,2)
        ELSE
          LENG  = LENGS1(KIND)
          IC(1) = IUNSTA(I,3)
        ENDIF

      ELSE
        WRITE (LUNIT,1) KIND, I, J
    1   FORMAT ('0ILLEGAL VALUES IN FORMIC =', 3I10)
        CALL ABORT2
      ENDIF

      RETURN
      END
      SUBROUTINE FORMOB (IUO,IUO2,B,G)

*** REFORM OBS EQUATIONS USING MOST RECENT PARAMETERS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL FATAL
      DIMENSION B(*),G(*)
      DIMENSION IC(LENC),C(LENC)
      DIMENSION COVECF(3,3)
      COMMON /OPT/   AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &               LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /CONST/  PI, PI2, RAD
      COMMON/UNITS/ LUNIT

*** LOOP OVER THE OBSERVATIONS

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
          CALL BIGL (KIND,ISN,JSN,IOBS,IVF,CMO,SD,FATAL)
          WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                 IVF,IAUX,IGRT

        ELSEIF (KIND .EQ. 18) THEN

*** CORRELATED TYPE (DOPPLER)

          CALL FORMC (KIND,C,B,ISN,JSN,IAUX,IGRT)
          CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
          CMO = OBS0 - OBSB
          CALL BIGL (KIND,ISN,JSN,IOBS,IVF,CMO,SD,FATAL)
          WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                 IVF,IAUX,IGRT
          READ (IUO,END=777) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                       IVF,IAUX,IGRT
          CALL FORMC (KIND,C,B,ISN,JSN,IAUX,IGRT)
          CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
          CMO = OBS0 - OBSB
          CALL BIGL (KIND,ISN,JSN,IOBS,IVF,CMO,SD,FATAL)
          WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                 IVF,IAUX,IGRT
          READ (IUO,END=777) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                       IVF,IAUX,IGRT
          CALL FORMC (KIND,C,B,ISN,JSN,IAUX,IGRT)
          CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,IAUX,IGRT)
          CMO = OBS0 - OBSB
          CALL BIGL (KIND,ISN,JSN,IOBS,IVF,CMO,SD,FATAL)
          WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                 IVF,IAUX,IGRT
          READ (IUO,END=777) ISIGNL, ( (COVECF(I,J), J = 1,3), I = 1,3 )
          WRITE (IUO2) ISIGNL, ( (COVECF(I,J), J = 1,3), I = 1,3 )
        ENDIF
      ELSE

*** DECORRELATED TYPE (GPS)

        NVEC = ISN
        IAUX = JSN
        NR = 3*NVEC
        WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &               IVF,IAUX,IGRT
        CALL FORMG (IUO,IUO2,NVEC,NR,G,B,FATAL,IVF,IAUX,IGRT)
      ENDIF
      GO TO 100

*** ABORT DUE TO LARGE MISCLOSURES

  777 IF (FATAL) THEN
        CALL LINE (3)
        WRITE (LUNIT,3) VM
    3   FORMAT ('0TERMINATED DUE TO MISCLOSURES (C-O)/SD EXCEEDING ',
     &          F7.1/)
        CALL ABORT2
      ENDIF

      REWIND IUO
      REWIND IUO2

*** EXCHANGE PRIMARY/SECONDARY OBS EQ FILE INDICATOR

      ITEMP = IUO
      IUO = IUO2
      IUO2 = ITEMP

      RETURN
      END
      
      SUBROUTINE FROM83 (I83, X83, Y83, Z83, X, Y, Z)

*** CONVERT X,Y,Z FROM ADJUSTMENT DATUM

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /UNITS/ LUNIT
      COMMON /XLATE/  RX(37), RY(37), RZ(37), SCALE(37),
     +             DELX(37), DELY(37), DELZ(37), COSRX(37), SINRX(37),
     +             COSRY(37), SINRY(37), COSRZ(37), SINRZ(37)
     
*** Abort program if reference system code greater than parameter array size.
      IF (I83 .GT. 37) THEN
        WRITE (LUNIT,1)
    1   FORMAT (1X, 130('*'))
        WRITE (LUNIT,2) I83
    2   FORMAT (' ERROR: UNDEFINED COORDINATE REFERENCE SYSTEM CODE = ',
     &          I2,' IN G-FILE.')
        WRITE (LUNIT,1)
        CALL ABORT2
      ENDIF

*** DO THE TRANSLATION

      X = X83 - DELX(I83)
      Y = Y83 - DELY(I83)
      Z = Z83 - DELZ(I83)

*** ROTATE ABOUT THE X AXIS

      XRX = X
      YRX = Y*COSRX(I83) - Z*SINRX(I83)
      ZRX = Z*COSRX(I83) + Y*SINRX(I83)

*** ROTATE ABOUT THE Y AXIS

      XRXY = XRX*COSRY(I83) + ZRX*SINRY(I83)
      YRXY = YRX
      ZRXY = ZRX*COSRY(I83) - XRX*SINRY(I83)

*** ROTATE ABOUT THE Z AXIS

      X = XRXY*COSRZ(I83) - YRXY*SINRZ(I83)
      Y = YRXY*COSRZ(I83) + XRXY*SINRZ(I83)
      Z = ZRXY

*** APPLY THE SCALE CHANGE

      XNEW = X - SCALE(I83)*X
      YNEW = Y - SCALE(I83)*Y
      ZNEW = Z - SCALE(I83)*Z

** 7/06/2005
      X = XNEW
      Y = YNEW
      Z = ZNEW
**************

      RETURN
      END
      
      SUBROUTINE GET83 (I83, XOLD, YOLD, ZOLD, XNEW, YNEW, ZNEW)

*** CONVERT X,Y,Z TO ADJUSTMENT DATUM

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /UNITS/ LUNIT
      COMMON /XLATE/  RX(37), RY(37), RZ(37), SCALE(37),
     +             DELX(37), DELY(37), DELZ(37), COSRX(37), SINRX(37),
     +             COSRY(37), SINRY(37), COSRZ(37), SINRZ(37)

*** Abort program if reference system code greater than parameter array size.
      IF (I83 .GT. 37) THEN
        WRITE (LUNIT,1)
    1   FORMAT (1X, 130('*'))
        WRITE (LUNIT,2) I83
    2   FORMAT (' ERROR: UNDEFINED COORDINATE REFERENCE SYSTEM CODE = ',
     &          I2,' IN G-FILE.')
        WRITE (LUNIT,1)
        CALL ABORT2
      ENDIF
      
*** DO THE TRANSLATION

      X = XOLD + DELX(I83)
      Y = YOLD + DELY(I83)
      Z = ZOLD + DELZ(I83)

*** ROTATE ABOUT THE X AXIS

      XRX = X
      YRX = Y*COSRX(I83) + Z*SINRX(I83)
      ZRX = Z*COSRX(I83) - Y*SINRX(I83)

*** ROTATE ABOUT THE Y AXIS

      XRXY = XRX*COSRY(I83) - ZRX*SINRY(I83)
      YRXY = YRX
      ZRXY = ZRX*COSRY(I83) + XRX*SINRY(I83)

*** ROTATE ABOUT THE Z AXIS

      X = XRXY*COSRZ(I83) + YRXY*SINRZ(I83)
      Y = YRXY*COSRZ(I83) - XRXY*SINRZ(I83)
      Z = ZRXY

*** APPLY THE SCALE CHANGE

      XNEW = X + SCALE(I83)*X
      YNEW = Y + SCALE(I83)*Y
      ZNEW = Z + SCALE(I83)*Z

      RETURN
      END
      SUBROUTINE GETALA (ALA,I,B)

*** ROUTINE TO RETRIEVE ASTRO LAT FROM CONTROL POINT DATA BLOCK

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
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETALA')
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      ALA = B(1+N)

      RETURN
      END
      SUBROUTINE GETALO (ALO,I,B)

*** ROUTINE TO RETRIEVE ASTRO LON FROM CONTROL POINT DATA BLOCK

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
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETALO')
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      ALO = B(2+N)

      RETURN
      END
      SUBROUTINE GETAUX (AUX, I, B)

*** ROUTINE TO RETRIEVE AUX. PARM. FROM CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I .LE. 0  .OR.  I .GT. NAUX) THEN
        WRITE (LUNIT,1) I, NAUX
    1   FORMAT ('0****ILLEGAL IAUX', I20, ' FOR NAUX=', I5,
     &          ' IN GETAUX ')
        CALL ABORT2
      ENDIF

      AUX = B(I)

      RETURN
      END
      SUBROUTINE GETDMS (VAL,ID,IM,S,ISIGN)

*** CONVERT RADIANS TO DEG, MIN, SEC

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
      SUBROUTINE GETECX (ECX,I,B)

*** ROUTINE TO RETRIEVE E.C.F. X FROM CONTROL POINT DATA BLOCK

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
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETECX')
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      ECX = B(5+N)

      RETURN
      END
      SUBROUTINE GETECY (ECY,I,B)

*** ROUTINE TO RETRIEVE E.C.F. Y FROM CONTROL POINT DATA BLOCK

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
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETECY')
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      ECY = B(6+N)

      RETURN
      END
      SUBROUTINE GETECZ (ECZ,I,B)

*** ROUTINE TO RETRIEVE E.C.F. Z FROM CONTROL POINT DATA BLOCK

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
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETECZ')
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      ECZ = B(7+N)

      RETURN
      END
      SUBROUTINE GETEHT (EHT,I,B)

*** ROUTINE TO RETRIEVE ELLIPSOIDAL HT. FROM CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETEHT')
        CALL ABORT2
      ENDIF
      N = NAUX+NZ+(I-1)*13
      EHT = B(10+N)

      RETURN
      END
      SUBROUTINE GETEHX (EHX,I,B)

*** ROUTINE TO RETRIEVE ELLIP HT X FROM CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETEHX')
        CALL ABORT2
      ENDIF
      N = NAUX+NZ+(I-1)*13
      EHX = B(11+N)

      RETURN
      END
      SUBROUTINE GETEHY (EHY,I,B)

*** ROUTINE TO RETRIEVE ELLIP. HT. Y FROM CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETEHY')
        CALL ABORT2
      ENDIF
      N = NAUX+NZ+(I-1)*13
      EHY = B(12+N)

      RETURN
      END
      SUBROUTINE GETEHZ (EHZ,I,B)

*** ROUTINE TO RETRIEVE ELLIP. HT. Z FROM CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NSTA) THEN
        WRITE (LUNIT,1) I,NSTA
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETEHZ')
        CALL ABORT2
      ENDIF
      N = NAUX+NZ+(I-1)*13
      EHZ = B(13+N)

      RETURN
      END
      SUBROUTINE GETGH (GH,I,B)

*** ROUTINE TO RETRIEVE GEOID HT. FROM CONTROL POINT DATA BLOCK

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
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETGH ')
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      GH = B(9+N)

      RETURN
      END
      SUBROUTINE GETGIC (ISN, JSN, IAUX, ICM, IC, N2, NICM, IGRT)

*** GET COLUMN INDEX IN WORK ARRAY FROM STATION NUMBERS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10, NVECS = 700 )
      DIMENSION ICM((NVECS+1)*3+1+3)
      LOGICAL GETICM
      DIMENSION IC(LENC)
      COMMON/UNITS/ LUNIT

      I4 = 4
      CALL FORMIC (I4, ISN, JSN, IAUX, IC, LENGIC, IGRT)

      DO 1 I = 1, LENGIC
        IF ( GETICM(IC(I), J, ICM, NICM) ) THEN
          IC(I) = N2 + J
        ELSE
          WRITE (LUNIT,666) ISN, JSN, ICM, IC
  666     FORMAT ('0ERROR IN GETGIC--',40I6)
          CALL ABORT2
        ENDIF
    1 CONTINUE

      RETURN
      END
      SUBROUTINE GETGLA (GLA,I,B)

*** ROUTINE TO RETRIEVE GEOD. LAT FROM CONTROL POINT DATA BLOCK

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
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETGLA')
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      GLA = B(3+N)

      RETURN
      END
      SUBROUTINE GETGLO (GLO,I,B)

*** ROUTINE TO RETRIEVE GEOD. LON FROM CONTROL POINT DATA BLOCK

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
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETGLO')
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      GLO = B(4+N)

      RETURN
      END
      LOGICAL FUNCTION GETGRT (ICODE,ITIME,IPARM)

*** TABLE LOOKUP

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PRMTB2/ NPARMS, ICODES(5), IOLDS(5), INEWS(5)

      IPARM = 1
  100 IF (IPARM.LE.NPARMS) THEN
        IF (ICODE.EQ.ICODES(IPARM) .AND.
     &      ITIME.GE.IOLDS(IPARM)  .AND.
     &      ITIME.LE.INEWS(IPARM)) THEN
          GETGRT = .TRUE.
          RETURN
        ELSE
          IPARM = IPARM+1
        ENDIF
        GO TO 100
      ENDIF

*** FELL THRU TABLE

      GETGRT = .FALSE.

      RETURN
      END
      LOGICAL FUNCTION GETICM (ICM, IICM, ICMS, NICMS)

*** TABLE LOOKUP

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION ICMS(*)

      IICM = 1
  100 IF (IICM .LE. NICMS) THEN
        IF ( ICM .EQ. ICMS(IICM) ) THEN
          GETICM = .TRUE.
          RETURN
        ELSE
          IICM = IICM + 1
        ENDIF
        GOTO 100
      ENDIF

*** FELL THRU TABLE

      GETICM = .FALSE.

      RETURN
      END
      LOGICAL FUNCTION GETIVF (ICODE,ITIME,IVF)

*** TABLE LOOKUP

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXVF = 40 )
      LOGICAL LFIXS
      COMMON /IVFTB1/ NIVFS, ICODES(MXVF), IOLDS(MXVF), INEWS(MXVF),
     &                LFIXS(MXVF)

      IVF = 1
  100 IF (IVF.LE.NIVFS) THEN
        IF (ICODE.EQ.ICODES(IVF) .AND.
     *      ITIME.GE.IOLDS(IVF)  .AND.
     *      ITIME.LE.INEWS(IVF)) THEN
          GETIVF = .TRUE.
          RETURN
        ELSE
          IVF = IVF+1
        ENDIF
        GO TO 100
      ENDIF

*** FELL THRU TABLE

      GETIVF = .FALSE.

      RETURN
      END
      SUBROUTINE GETLAH (I,J,B,P1,Q1,R1,T1,P2,Q2,R2,T2,S)

*** RETRIEVE THE LOCAL ASTRO HORIZON SYSTEM(P,Q,R,S,T)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)

*** RETRIEVE THE EARTH-CENTERED COORDS FOR STATIONS I & J
*** NOTE: THIS SUBROUTINE NOT CALLED BY GPS OR DOPPLER OBSERVATIONS
*** THUS, DO NOT NEED TO ACCOUNT FOR DUAL HEIGHT IN 2.5DIM ADJUSTMENT

      CALL GETECX (XI,I,B)
      CALL GETECX (XJ,J,B)
      CALL GETECY (YI,I,B)
      CALL GETECY (YJ,J,B)
      CALL GETECZ (ZI,I,B)
      CALL GETECZ (ZJ,J,B)

*** FORMULAS FROM HAVAGO(TECH. MEMO #17)

      DX = XJ-XI
      DY = YJ-YI
      DZ = ZJ-ZI

      CALL GETALA (ALA1,I,B)
      CALL GETALA (ALA2,J,B)
      CALL GETALO (ALO1,I,B)
      CALL GETALO (ALO2,J,B)

      SLA1 = DSIN(ALA1)
      SLO1 = DSIN(ALO1)
      SLA2 = DSIN(ALA2)
      SLO2 = DSIN(ALO2)
      CLA1 = DCOS(ALA1)
      CLO1 = DCOS(ALO1)
      CLA2 = DCOS(ALA2)
      CLO2 = DCOS(ALO2)

      P1 = -SLA1*(DX*CLO1+DY*SLO1)+DZ*CLA1
      Q1 = -DX*SLO1+DY*CLO1
      R1 = DSQRT(P1*P1+Q1*Q1)
      T1 = +CLA1*(DX*CLO1+DY*SLO1)+DZ*SLA1
      P2 = +SLA2*(DX*CLO2+DY*SLO2)-DZ*CLA2
      Q2 = +DX*SLO2-DY*CLO2
      R2 = DSQRT(P2*P2+Q2*Q2)
      T2 = -CLA2*(DX*CLO2+DY*SLO2)-DZ*SLA2
      S = DSQRT(DX*DX+DY*DY+DZ*DZ)

      RETURN
      END
      SUBROUTINE GETLGH (I,J,B,C)

*** TRANSFORM COEFFICIENTS TO LOCAL GEODETIC HORIZON

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      DIMENSION B(*),C(LENC)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT

      CALL GETALA (ALA1,I,B)
      CALL GETALO (ALO1,I,B)
      CALL GETGLA (GLA1,I,B)
      CALL GETGLO (GLO1,I,B)

      CALL GETALA (ALA2,J,B)
      CALL GETALO (ALO2,J,B)
      CALL GETGLA (GLA2,J,B)
      CALL GETGLO (GLO2,J,B)

      XM1 = DSIN(ALA1)*(ALO1-GLO1)
      XL1 = ALA1-GLA1
      XN1 = DCOS(ALA1)*(ALO1-GLO1)

      XM2 = DSIN(ALA2)*(ALO2-GLO2)
      XL2 = ALA2-GLA2
      XN2 = DCOS(ALA2)*(ALO2-GLO2)

      IF (IDIM.EQ.1) THEN
        CONTINUE
      ELSEIF (IDIM.EQ.2) THEN
        C(1) = +C(1)+XM1*C(2)
        C(2) = -XM1*C(1)+C(2)
        C(3) = +C(3)+XM2*C(4)
        C(4) = -XM2*C(3)+C(4)
      ELSE
        C(1) = +C(1)+XM1*C(2)+XL1*C(3)
        C(2) = -XM1*C(1)+C(2)+XN1*C(3)
        C(3) = -XL1*C(1)-XN1*C(2)+C(3)
        C(4) = +C(4)+XM2*C(5)+XL2*C(6)
        C(5) = -XM2*C(4)+C(5)+XN2*C(6)
        C(6) = -XL2*C(4)-XN2*C(5)+C(6)
      ENDIF

      RETURN
      END
      SUBROUTINE GETMSL (GMSL,I,B)

*** ROUTINE TO RETRIEVE MSL FROM CONTROL POINT DATA BLOCK

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
    1   FORMAT ('0****ILLEGAL ISN',I5,' FOR NSTA=',I5,' IN GETMSL')
        CALL ABORT2
      ENDIF
      IF (.NOT.L2HLF) THEN
        N = NAUX+NZ+(I-1)*9
      ELSE
        N = NAUX+NZ+(I-1)*13
      ENDIF
      GMSL = B(8+N)

      RETURN
      END
      
     
      LOGICAL FUNCTION GETPRM (ICODE, ITIME, IPARM)

*** TABLE LOOKUP

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXPRM = 40 )
      COMMON /PRMTB1/ NPARMS, ICODES(MXPRM), IOLDS(MXPRM), INEWS(MXPRM)

      IPARM = 1
  100 IF (IPARM .LE. NPARMS) THEN
        IF ( ICODE .EQ. ICODES(IPARM) .AND.
     &       ITIME .GE. IOLDS(IPARM)  .AND.
     &       ITIME .LE. INEWS(IPARM) ) THEN
          GETPRM = .TRUE.
          RETURN
        ELSE
          IPARM = IPARM + 1
        ENDIF
        GO TO 100
      ENDIF

*** FELL THRU TABLE

      GETPRM = .FALSE.

      RETURN
      END
      SUBROUTINE GETQ (SIGOBS, COVLA, Q)

***  GET Q MATRIX (CHOLESKY FACTORIZATION FROM HOG - FUNCTION ADDCOR)
***  SIGOBS IS CHANGED TO CHOLESKY MATRIX (UPPER)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION SIGOBS(3,3), COVLA(3,3), Q(3,3)
      DIMENSION DUMVEC(3)
      DATA TOL /1.D-10/
      COMMON/UNITS/ LUNIT

*** CHOLESKY FACTOR WEIGHT MATRIX

      DO 1 I = 1, 3
        IM1 = I - 1
        IF (IM1 .GT. 0) THEN
          DO 2 K = 1, IM1
            SIGOBS(I,I) = SIGOBS(I,I) - SIGOBS(K,I)*SIGOBS(K,I)
    2     CONTINUE
        ENDIF

*** SINGULARITY TEST

        IF (SIGOBS(I,I) .LT. TOL) THEN
          WRITE (LUNIT,661)
  661     FORMAT (' PROGRAMMER ERROR IN GETQ, DOPPLER COV MATRIX',
     &            ' SIGOBS IS NOT POSITIVE DEFINITE')
          CALL ABORT2
          RETURN
        ENDIF

*** RESUME CHOLESKY FACTOR

        SIGOBS(I,I) = DSQRT( SIGOBS(I,I) )
        IF (I+1 .LE. 3) THEN
          DO 6 J = I+1, 3
            IF (IM1 .GT. 0) THEN
              DO 4 K = 1, IM1
                SIGOBS(I,J) = SIGOBS(I,J) - SIGOBS(K,I)*SIGOBS(K,J)
    4         CONTINUE
            ENDIF
            SIGOBS(I,J) = SIGOBS(I,J)/SIGOBS(I,I)
            SIGOBS(J,I) = 0.D0
    6     CONTINUE
        ENDIF
    1 CONTINUE

      CALL ABAT (SIGOBS, COVLA, Q, DUMVEC, 3, 3)
      DO 60 I = 1, 3
        Q(I,I) = 1.D0 - Q(I,I)
   60 CONTINUE

      RETURN
      END
      SUBROUTINE GETRAD (ID,IM,IS,SIGN,VAL)

*** CONVERT DEG, MIN, SEC TO RADIANS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /CONST/ PI, PI2, RAD

      S = DBLE(IS)/1.D5

      VAL = (ID+IM/60.D0+S/3600.D0)/RAD
      VAL = DSIGN(VAL,SIGN)

      RETURN
      END
      SUBROUTINE GETRHS (G, NR, NC, N, M)

*** COPY M-TH COLUMN TO N-TH COLUMN

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION G(NR,NC)

      DO 1 I = 1, NR
        G(I,N) = G(I,M)
    1 CONTINUE

      RETURN
      END
      SUBROUTINE GETRI (RI,ISN,B)

*** BUILD A NEW ROTATION MATRIX IF NEW ISN

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*),RI(3,3)
      COMMON /LASTSN/ ISNX, JSNX

      IF (ISN.NE.ISNX) THEN
        CALL BUILDR (RI,ISN,B)
        ISNX = ISN
      ENDIF

      RETURN
      END
      SUBROUTINE GETRJ (RJ,JSN,B)

*** BUILD A NEW ROTATION MATRIX IF NEW JSN

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*),RJ(3,3)
      COMMON /LASTSN/ ISNX, JSNX

      IF (JSN.NE.JSNX) THEN
        CALL BUILDR (RJ,JSN,B)
        JSNX = JSN
      ENDIF

      RETURN
      END
      SUBROUTINE GETROT (ROT,I,B)

*** ROUTINE TO RETRIEVE ROT. PARM. FROM CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NZ) THEN
        WRITE (LUNIT,1) I,NZ
    1   FORMAT ('0****ILLEGAL IZ',I20,' FOR NZ=',I5,' IN GETROT')
        CALL ABORT2
      ENDIF

      ROT = B(I+NAUX)

      RETURN
      END
      SUBROUTINE GETRT (DELX,DELY,DELZ,RT)

*** RETURN THE PRODUCT OF THE GPS VECTOR DIFFERENCES TIMES THE R0 MATRIX
*** R0 IS THE TRANSPOSE OF THE GEODETIC HORIZON MATRIX FOR GLAT0, GLON0

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION RT(3,3)
      COMMON /GRTTB1/ GLAT0, GLON0, R0(3,3)


      RT(1,1) = -DELZ * R0(1,2) + DELY * R0(1,3)
      RT(1,2) = -DELZ * R0(2,2) + DELY * R0(2,3)
      RT(1,3) = -DELZ * R0(3,2) + DELY * R0(3,3)

      RT(2,1) =  DELZ * R0(1,1) - DELX * R0(1,3)
      RT(2,2) =  DELZ * R0(2,1) - DELX * R0(2,3)
      RT(2,3) =  DELZ * R0(3,1) - DELX * R0(3,3)

      RT(3,1) = -DELY * R0(1,1) + DELX * R0(1,2)
      RT(3,2) = -DELY * R0(2,1) + DELX * R0(2,2)
      RT(3,3) = -DELY * R0(3,1) + DELX * R0(3,2)

      RETURN
      END
      SUBROUTINE GETRTG (GRTX,GRTY,GRTZ,I)

*** ROUTINE TO RETRIEVE AUX GPS AND DOPPLER ROTATION PARAMETERS
*** FROM CONTROL POINT DATA BLOCK

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /GRTTB2/ BGRT(5,3)
      COMMON/UNITS/ LUNIT

      IF (I.LE.0 .OR. I.GT.NGRT ) THEN
        WRITE (LUNIT,1) I,NGRT
    1   FORMAT ('0****ILLEGAL IGRT',I20,' FOR NGRT=',I5,' IN GETRTG ')
        CALL ABORT2
      ENDIF

      GRTX = BGRT(I,1)
      GRTY = BGRT(I,2)
      GRTZ = BGRT(I,3)

      RETURN
      END
C---------------------------------------------------------------------------
      LOGICAL FUNCTION GETSSN (ISSN,I)

*** TABLE LOOKUP

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      COMMON /SSNTBL/ NSSN, ISSISN(MXSSN)
      COMMON/UNITS/ LUNIT

      IF (ISSN.LE.0 .OR. ISSN.GT.MXSSN) THEN
        WRITE (LUNIT,1) ISSN
    1   FORMAT ('0ILLEGAL ISSN ',I9)
        CALL ABORT2
      ENDIF

      I = ISSISN(ISSN)
      IF (I.EQ.0) THEN
        GETSSN = .FALSE.
      ELSE
        GETSSN = .TRUE.
      ENDIF

      RETURN
      END
C-----------------------------------------------------------------------------
      SUBROUTINE GETYR(RECORD,IYR)
      
*** DETERMINE 4-DIGIT B-FILE CLASSICAL OBS FROM *12* REC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 RECORD
      COMMON /YEARS/ ICBEG,IYBEG,ICEND,IYEND
      COMMON/UNITS/ LUNIT

      IF (ICBEG.EQ.ICEND) THEN
         IYR = ICBEG*100 + IYR
      ELSE
         IB = IYR-IYBEG
         IF(IB.LT.0) IB = -IB
         IE = IYR-IYEND
         IF(IE.LT.0) IE = -IE
         IF(IB.LE.IE) THEN
           IYR = ICBEG*100 + IYR
         ELSE
           IYR = ICEND*100 + IYR
         ENDIF
      ENDIF
      IF((IYR.LT.(ICBEG*100+IYBEG)).OR.
     *(IYR.GT.(ICEND*100+IYEND))) THEN
      WRITE(LUNIT,1) RECORD
 1    FORMAT(/' FOR B-FILE OBSERVATION:'/' ',A80/)
      WRITE(LUNIT,10) IYR,ICBEG*100+IYBEG,ICEND*100+IYEND
 10   FORMAT(' THE DERIVED 4-DIGIT YEAR ',I4/' IS NOT BETWEEN *12*',
     *' YEAR OBS BEGAN ',I4,' AND YEAR OBS END ',I4/)
      CALL ABORT2
      ENDIF
 
      RETURN
      END
      
      LOGICAL FUNCTION GETZ (ISSN,LIST,IZ)

*** TABLE LOOKUP

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MAXZ =  8000 )
      COMMON /ROTTB1/ IZNS(MAXZ)
      COMMON /ROTTB2/ LOOK(977), IPTRS(MAXZ), NKEY, MAX, IFREE
      LOGICAL LOKZ

      IZN = ISSN*1000+LIST
      IF (LOKZ(IZN,ILOOK,IPTR)) THEN
        IZ = IPTR
        GETZ = .TRUE.
      ELSE
        IZ = IFREE
        GETZ = .FALSE.
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION GET82(ISSN,I)

*** TABLE LOOKUP OF *82* RECORDS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      COMMON /CHILD/ N82,ISS82(MXSSN)
      COMMON/UNITS/ LUNIT

      IF (ISSN.LE.0 .OR. ISSN.GT.MXSSN) THEN
        WRITE(LUNIT,1) ISSN
    1   FORMAT ('0ILLEGAL ISSN ',I9)
        CALL ABORT2
      ENDIF

      I = ISS82(ISSN)
      IF (I.EQ.0) THEN
        GET82 = .FALSE.
      ELSE
        GET82 = .TRUE.
      ENDIF

      RETURN
      END
      SUBROUTINE GLOCAT (N1, N2, N3, N4, N5, NVEC, LENG)

*** DETERMINE WORK ARRAY COLUMN POINTERS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      N1 = 3*NVEC
      N2 = N1 + 1
      N3 = N2 + LENG
      N4 = N3 + 1
      N5 = N4 + 1

      RETURN
      END
      SUBROUTINE GOBSEQ (G, NR, NC, NVEC, LENG, B, ICM, NICM,KINDS,ISNS,
     &                   JSNS, LOBS, IAUX, IVF, FATAL, NOBIGV, N4, IGRT)

*** COMPUTE AND DECORRELATE OBS. EQ. AND MISCLOSURE

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10, MXSSN = 9999 )
      LOGICAL FATAL, NOBIGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      DIMENSION KINDS(*), ISNS(*), JSNS(*), LOBS(*)
      DIMENSION ICM(*)
      DIMENSION B(*), G(NR,NC)
      DIMENSION IC(LENC), C(LENC)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /LASTSN/ ISNX, JSNX
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT

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

*** COMPUTE COEFFICIENTS AND OBSERVATIONS

      DO 10 I = 1, NR
        KIND = KINDS(I)
        ISN = ISNS(I)
        JSN = JSNS(I)
        IOBS = LOBS(I)
        CALL COMPOB (KIND, OBS0, B, ADUMMY, ISN, JSN, IAUX, IGRT)
        CALL FORMC (KIND, C, B, ISN, JSN, IAUX, IGRT)
        CALL GETGIC (ISN, JSN, IAUX, ICM, IC, N2, NICM, IGRT)

*** LOAD WORK ARRAY

        CMO = OBS0 - G(I,NC)
        IF (IMODE .EQ. 0) CMO = 0.D0
        G(I,N2) = CMO

        DO 2 J = 1, LENG
          K = J + N2
          G(I,K) = 0.D0
    2   CONTINUE

        DO 20 J = 1, LENGC
          G(I,IC(J) ) = C(J)
   20   CONTINUE
   10 CONTINUE

*** DECORRELATE MISCLOSURE

      CALL COMRHS (G, NR)
      CALL PUTRHS (G, NR, NC, N2, N4)

*** TEST FOR LARGE MISCLOSURES

      IF ( .NOT. NOBIGV) THEN
        SD = 1.D0
        DO 15 I = 1, NR
          CMO = G(I,N2)
          CALL BIGV (KINDS(I), ISNS(I), JSNS(I), LOBS(I), IVF, CMO, SD,
     &               FATAL)
   15   CONTINUE
      ENDIF

*** DECORRELATE COEFFICIENTS

      DO 30 I = 1, LENG
        CALL GETRHS (G, NR, NC, N2, N2+I)
        CALL COMRHS (G, NR)
        CALL PUTRHS (G, NR, NC, N2, N2+I)
   30 CONTINUE

      RETURN
      END
      SUBROUTINE GPSSDV (G, NR, NC, NVEC, LENG, B, ICM, NICM, KINDS,
     &                   ISNS, JSNS, LOBS, IAUX, IVF, FATAL, NOBIGV,
     &                   A, NX, N2, N3, N4, GMDES, IGRT)

*** COMPUTE AND DECORRELATE OBS. EQ. AND MISCLOSURE

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( NVECS = 700, LENC = 10, MXVF = 40, MXSSN = 9999 )
      DIMENSION ICM((NVECS+1)*3+1+3)
      DIMENSION KINDS(NVECS*3), ISNS(NVECS*3), JSNS(NVECS*3)
      DIMENSION LOBS(NVECS*3)
      LOGICAL FATAL, NOBIGV, PROP, PROPCV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      DIMENSION CM((NVECS+1)*3+1+3), CM2((NVECS+1)*3+1+3), GMDES(*)
      DIMENSION B(*), G(NR,NC), A(*), NX(*)
      DIMENSION IC(LENC), C(LENC)
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON /VFTB2/  VFS(MXVF), VTV(MXVF), VFRN(MXVF), VFSS(MXVF)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /LASTSN/ ISNX, JSNX
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
** v 4.32vf
       LOGICAL IVCGPS,VSlogic
      COMMON/VCREC/SIGH,SIGU, VFHOR,VFUP, IVCGPS
      CHARACTER*13 VSPROJ
      PARAMETER (MAXPRJ=2500)
      COMMON/VSREC/SIGHS(MAXPRJ), SIGUS(MAXPRJ), ISETHU, VSPROJ(MAXPRJ),
     &             VSlogic
********
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

      NB = NUNK + 1
      NL = NR*(NR+1)/2
      NB1 = NX(NB)
      NB2 = NB1 + NL
      NB3 = NB2 + NL

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
***            SCALE BY VARIANCE FACTOR AS NEEDED
***            (NOTE: COV. OF PARMS IS SCALED BY FACTORS IN SUB. NORMAL)

      DO 60 IROW = 1, NR
        DO 50 I = 1, NICM
          CM(I) = G(IROW,N2+I)
   50   CONTINUE

        IF ( .NOT. PROP(CM, ICM, NICM, VARL0, A, NX, IFLAG) ) THEN
          IF (IFLAG .EQ. 1) THEN
            WRITE (LUNIT,666)
  666       FORMAT ('0STATE ERROR IN FIRST PROP OF GPSSDV')
            CALL ABORT2
          ELSE
            WRITE (LUNIT,667)
  667       FORMAT ('0PROFILE ERROR IN FIRST PROP OF GPSSDV')
            CALL ABORT2
          ENDIF
        ENDIF

        IDEX = INX(IROW,NR)
        A(IDEX+NB1) = VARL0
** v 4.30vf **************************************
C       IF (IROW .EQ. 1) THEN
C         VARLB = G(1,1)*G(1,1)
C       ELSE
C         VARLB = G(IROW,1)
C       ENDIF
        VARLB=G(IROW,N5)**2
*****     IF (ISETHU.OR.IVCGPS) THEN

          IF(ISETHU .GT. 0 .OR. IVCGPS) THEN

C  get variance of observation from the Cholesky factor
          VARLB=0.0
            DO 55 I=1,IROW
   55     VARLB=VARLB+G(I,IROW)**2
         ENDIF
***************************************************
        IF (NVFTOT .GT. 0) THEN
          IF (IVF .LT. 0  .OR.  IVF .GT. NVFTOT) THEN
            WRITE (LUNIT,663) IVF, NVFTOT
  663       FORMAT ('0ILLEGAL IVF=',I10,' FOR N=',I10,' IN GPSSDV')
            CALL ABORT2
          ELSEIF (IVF .NE. 0) THEN
            VFACTR = VFS(IVF)
          ELSE
            VFACTR = 1.D0
          ENDIF
        ELSE
          VFACTR = 1.D0
        ENDIF
        VARLB = VARLB*VFACTR

        DIFF = VARLB - VARL0
        IF (DIFF .LT. 1.D-10) DIFF = 0.D0
        SDV = DSQRT(DIFF)
        G(IROW,N4) = SDV
   60 CONTINUE

*** (STEP 3.1)  LOAD ADJ. OBS. COV. MATRIX INTO 1ST SCRATCH SPACE

      DO 30 I = 1, NR
        DO 31 K = 1, NICM
          CM(K) = G(I,N2+K)
   31   CONTINUE
        IDEX = INX(I,NR) - I
        DO 30 J = I, NR
          IF (I .NE. J) THEN
            DO 32 K = 1, NICM
              CM2(K) = G(J,N2+K)
   32       CONTINUE
            IF ( .NOT. PROPCV(CM, ICM, NICM, CM2, ICM, NICM,
     &                       COV, A, NX, IFLAG) ) THEN
              IF (IFLAG .EQ. 1) THEN
                WRITE (LUNIT,668)
  668           FORMAT ('0STATE ERROR IN PROPCV OF GPSSDV')
                CALL ABORT2
              ELSE
                WRITE (LUNIT,669)
  669           FORMAT ('0PROFILE ERROR IN PROPCV OF GPSSDV')
                CALL ABORT2
              ENDIF
            ENDIF
            A(NB1+IDEX+J) = COV
          ENDIF
   30 CONTINUE


*** (STEP 3.2)  COPY OBS. COV. CHOLESKY FACTOR INTO 2ND SCRATCH SPACE

      IPOINT = NB2
      DO 40 I = 1, NR
        DO 40 J = I, NR
          IPOINT = IPOINT + 1
          A(IPOINT) = G(I,J)
   40 CONTINUE

*** (STEP 3.3)  INVERT IMAGE OF OBS. COV. CHOLESKY (A(NB3) IS SCRATCH)

      CALL AVERT (A(NB2+1), A(NB3+1), NR)

*** (STEP 4)  COMPUTE MARGINALLY DETECTABLE ERRORS (1 SIGMA)
***           TABLE GMDES(*) IS PART OF PATCH WORK

      IF (IMODE .EQ. 0  .OR.  IMODE .EQ. 3) THEN
        DO 70 IROW = 1, NR
          DO 75 JCOL = 1, NICM
            CM(JCOL) = 0.D0
            DO 75 K = 1, NR
              IF (IROW .LE. K) THEN
                INDEX = INX(IROW,NR) - IROW + K
              ELSE
                INDEX = INX(K,NR) - K + IROW
              ENDIF
              CM(JCOL) = CM(JCOL) + G(K,N2+JCOL)*A(NB2+INDEX)
   75     CONTINUE
          IF ( .NOT. PROP(CM, ICM, NICM, VAR, A, NX, IFLAG) ) THEN
            IF (IFLAG .EQ. 1) THEN
              WRITE (LUNIT,664)
  664         FORMAT ('0STATE ERROR IN SECOND PROP OF GPSSDV')
              CALL ABORT2
            ELSE
              WRITE (LUNIT,665)
  665         FORMAT ('0PROFILE ERROR IN SECOND PROP OF GPSSDV')
              CALL ABORT2
            ENDIF
          ENDIF
          IF (NVFTOT .GT. 0) THEN
            IF (IVF .LT. 0  .OR.  IVF .GT. NVFTOT) THEN
              WRITE (LUNIT,663) IVF, NVFTOT
              CALL ABORT2
            ELSEIF (IVF .NE. 0) THEN
              VFACTR = VFS(IVF)
            ELSE
              VFACTR = 1.D0
            ENDIF
          ELSE
            VFACTR = 1.D0
          ENDIF
          DEINV = A( NB2+INX(IROW,NR) )*VFACTR - VAR
          IF (DEINV .LE. 1.D-6) THEN
            GMDES(IROW) = 0.D0
          ELSE
            GMDES(IROW) = 1.D0/DEINV
          ENDIF
   70   CONTINUE
      ENDIF

*** (STEP 5.)  COMPUTE OBSERVATIONS AND MISCLOSURES

      IF (IMODE .NE. 0) THEN
        DO 90 I = 1, NR
          KIND = KINDS(I)
          ISN = ISNS(I)
          JSN = JSNS(I)
          CALL COMPOB (KIND, OBS0, B, ADUMMY, ISN, JSN, IAUX, IGRT)
          CMO = OBS0 - G(I,NC)
          G(I,N2) = CMO
   90   CONTINUE
      ENDIF

*** (STEP 6.)  DECORRELATE MISCLOSURE

      IF (IMODE .NE. 0) THEN
        CALL COMRHS (G, NR)

        IF (NVFTOT .GT. 0  .AND.  IVF .GT. 0) THEN
          SFACTR = DSQRT(VFACTR)
          DO 15 I = 1, NR
            G(I,N2) = G(I,N2)/SFACTR
   15     CONTINUE
        ENDIF

*** TEST FOR LARGE MISCLOSURES

        IF ( .NOT. NOBIGV) THEN
          SD = 1.D0
          DO 16 I = 1, NR
            CMO = G(I,N2)
            CALL BIGL (KINDS(I), ISNS(I), JSNS(I), LOBS(I),
     &                 IVF, CMO, SD, FATAL)
   16     CONTINUE
        ENDIF
      ENDIF

*** (STEP 7)  COMPUTE REDUNDANCY NUMBERS USING IMAGES
***           CAUTION--THIS STEP DAMAGES THE COEFFICIENTS IN G

      DO 80 I = 1, NR
        W = 0.D0
        DO 85 J = 1, NR
          IF (I .LE. J) THEN
            K = INX(I,NR) - I + J
          ELSE
            K = INX(J,NR) + I - J
          ENDIF
          SIGA = A(NB1+K)
          SIGL = A(NB2+K)
            W = W + SIGA*SIGL
   85   CONTINUE
        W = W/VFACTR
        RN = 1.D0 - W
        G(I,N3) = RN
   80 CONTINUE

      RETURN
      END
      SUBROUTINE HEAD2

*** PRINT A HEADING AND SUBHEADING #2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON/UNITS/ LUNIT

      CALL HEAD

      IF (IMODE.EQ.0) THEN
        WRITE (LUNIT,1)
** v 4.29f
*1      FORMAT ('0SIMULATION'/T18,'COMPUTED',T30,'MDE(3-SIGMA)',
*    &           T50,'RN',T57,'FROM STATION',T89,'TO STATION', 
*    &           /T57,'SESSION ID')
 1      FORMAT ('0SIMULATION'/T18,'COMPUTED',T30,'MDE(3-SIGMA)',
     &           T51,'RN',T58,'FROM STATION',T96,'TO STATION', 
     &           /T58,'SESSION ID')
*********************************************
     
      ELSEIF (IMODE.EQ.3) THEN
        WRITE (LUNIT,2)
 2      FORMAT ('0NORMALIZED RESIDUALS'/T18,'COMPUTED',T34,'OBSERVED',
     &       T51,'V=C-O',T69,'SDV     V/SDV',T87,'MDE',T96,'RN',T100,
     &       'FROM STATION'/T49,'SEC',T57,'METER',T85,'3-SIGMA',
** v 4.28
*    &       T102,'TO STATION(S)')
     &       T100,'TO STATION(S)',/T100,'SESSION ID')
      ELSE
        WRITE (LUNIT,3)
**v 4.28
*3      FORMAT ('0QUASI-NORMALIZED RESIDUALS'/T17,'COMPUTED',T35,
*    &          'OBSERVED',T53,'V=C-O',T68,'V/SD'/T51,'SEC',T59,'METER')
** v 4.29f
*3      FORMAT ('0QUASI-NORMALIZED RESIDUALS'/T17,'COMPUTED',T35,
*    &          'OBSERVED',T53,'V=C-O',T68,'V/SDV',T74,'FROM STATION',
*    &         T105,'TO STATION',/T51,'SEC',T59,'METER',T74,'SESSID')  
 3      FORMAT ('0QUASI-NORMALIZED RESIDUALS'/T17,'COMPUTED',T35,
     &          'OBSERVED',T53,'V=C-O',T68,'V/SDV',T74,'FROM STATION',
     &         T111,'TO STATION',/T51,'SEC',T59,'METER',T74,'SESSID')  
************************************     
     
      ENDIF

      IF (IMODE.NE.3) THEN
        ILINE = 10
      ELSE
        ILINE = 12
      ENDIF

      RETURN
      END
      SUBROUTINE HEAD3

*** PRINT A HEADING AND SUBHEADING #3

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE
      COMMON/UNITS/ LUNIT

      CALL HEAD

      WRITE (LUNIT,1)
    1 FORMAT ('0ADJUSTED POSITIONS', //,
     +        8X, 'SSN', 1X, 'NAME', 29X, 'LATITUDE', 9X, 'LONGITUDE',
**     +        10X, 'M.S.L.',  3X, 'G. HT.', 4X, 'E. HT.')
**   +        10X, 'ORTHO HT',  1X, 'G. HT.', 4X, 'E. HT.')
     +        6X, 'ORTHO. HT.',  3X, 'G. HT.', 4X, 'E. HT.')
      ILINE = 10

      RETURN
      END
      SUBROUTINE HEAD4

*** PRINT A HEADING AND SUBHEADING #4

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE
      COMMON/UNITS/ LUNIT

      CALL HEAD

      IF (LSS) THEN
        WRITE (LUNIT,1)
    1   FORMAT ('0LENGTH RELATIVE ACCURACIES',
     *          3X,'(USING A-POSTERIORI WEIGHTS)')
      ELSE
        WRITE (LUNIT,2)
    2   FORMAT ('0LENGTH RELATIVE ACCURACIES',
     *          3X,'(USING A-PRIORI WEIGHTS)')
      ENDIF

      ILINE = 7

      RETURN
      END
      SUBROUTINE HEAD5

*** PRINT A HEADING AND SUBHEADING #5

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE
      COMMON/UNITS/ LUNIT

      CALL HEAD

      WRITE (LUNIT,1)
    1 FORMAT ('0ADJUSTED AUXILIARY PARAMETERS'//T2,'NUM',T34,'VALUE')

      ILINE = 9

      RETURN
      END
      SUBROUTINE HEAD6

*** PRINT A HEADING AND SUBHEADING #6

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      LOGICAL L2HLF, LEHT
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON/UNITS/ LUNIT

      CALL HEAD

*** ISSUE WARNING FOR DUAL HT STATIONS

      IF (IDIM .EQ. 2  .AND.  L2HLF) THEN
        ILINE = ILINE + 2
        WRITE (LUNIT,2)
   2    FORMAT ('0*** OBSERVATIONAL SUMMARY ***'
     &           /, ' WARNING - FOR THE OBSERVATION SUMMARY',
     &              ' COMPUTATIONS, DUAL HEIGHT STATIONS',
     &              ' ARE TREATED AS 2-DIM STATIONS',
     &         //,  3X, 'SSN  CMP  STATION NAME',
     &             21X, 'DIR', 9X, 'ANG', 9X, 'AZI', 9X, 'DIS',
     &              9X, 'ZD', 10X, 'GPS', 7X, 'DOP',
     &          /, 44X, 'FRM   TO', 4X, 'FRM   TO', 4X, 'FRM   TO',
     &              4X, 'FRM   TO', 4X, 'FRM   TO', 4X, 'FRM   TO', /)
        ILINE = 13
      ELSE

        WRITE (LUNIT,1)
   1    FORMAT ('0*** OBSERVATIONAL SUMMARY ***'
     &         //,  3X, 'SSN  CMP  STATION NAME',
     &             21X, 'DIR', 9X, 'ANG', 9X, 'AZI', 9X, 'DIS',
     &              9X, 'ZD', 10X, 'GPS', 7X, 'DOP',
     &          /, 44X, 'FRM   TO', 4X, 'FRM   TO', 4X, 'FRM   TO',
     &              4X, 'FRM   TO', 4X, 'FRM   TO', 4X, 'FRM   TO', /)

        ILINE = 10
      ENDIF

      RETURN
      END
      SUBROUTINE HEAD7

*** PRINT A HEADING AND SUBHEADING #7

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON/UNITS/ LUNIT

      CALL HEAD

      IF (IMODE.EQ.3) THEN
        WRITE (LUNIT,1)
 1      FORMAT ('0NORMALIZED RESIDUALS OF HORIZONTAL OBSERVATIONS',
     &          ' GROUPED AROUND INTERSECTION STATIONS', /,
     &          13X, 'COMPUTED', T31, 'OBSERVED', T50, 'V=C-O',
     &          T64, 'V/SDV RN', /,
     &          48X, 'SEC', T57, 'METER')
      ELSE
        WRITE (LUNIT,2)
 2      FORMAT ('0QUASI-NORMALIZED RESIDUALS OF HORIZONTAL OBSERVATION',
     &          ' GROUPED AROUND INTERSECTION STATIONS', /,
     &          12X, 'COMPUTED', T32, 'OBSERVED', T53,'V=C-O',
     &          T67, 'V/SDO', /,
     &          50X, 'SEC', T59, 'METER')
      ENDIF

      ILINE = 8

      RETURN
      END
      SUBROUTINE HEAD8

*** PRINT A HEADING AND SUBHEADING #8

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*1 ADIR1,ADIR2
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE
      COMMON /GRTTB1/ GLAT0, GLON0, R0(3,3)
      COMMON/UNITS/ LUNIT

      CALL HEAD

      CALL GETDMS (GLAT0,ID1,IM1,S1,ISIGN)
      IF (ISIGN.GT.0) THEN
        ADIR1 = 'N'
      ELSE
        ADIR1 = 'S'
      ENDIF
      CALL GETDMS (GLON0,ID2,IM2,S2,ISIGN)
      IF (ISIGN.GT.0) THEN
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

      WRITE (LUNIT,1) ID1,IM2,S1,ADIR1,ID2,IM2,S2,ADIR2
    1 FORMAT ('0ADJUSTED AUXILIARY GPS AND DOPPLER ROTATION PARAMETERS',
     &        //,' ROTATION ORIGIN IS ',2I3,F9.5,A1,3X,I4,I3,F9.5,A1,
     &        //,' NUM',17X,'X ROTATION',15X,'Y ROTATION',
     &                  15X,'Z ROTATION')

      ILINE = 12

      RETURN
      END
      SUBROUTINE HEAD9

*** PRINT A HEADING AND SUBHEADING #9

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE
      COMMON/UNITS/ LUNIT

      CALL HEAD

      WRITE (LUNIT,1)
    1 FORMAT ( '0DUAL HEIGHT DIFFERENCES', //,
     &          48X, 'ORIGINAL', 24X, 'ADJUSTED', /,
     &          34X, 2( 3X, '------------------------------'), /,
     &          3X, 'SSN NAME', 24X,
**    &          2( 5X, 'M.S.L.    G. HT.     E. HT.'), 6X, 'DIFF.')
     &          2( 5X, 'ORTHO HT G. HT.     E. HT.'), 6X, 'DIFF.')

      ILINE = 11

      RETURN
      END
      SUBROUTINE HORANG (BCARD,IUO,IOBS,B,NX,FATAL,LSN)

*** OBSERVATION EQUATIONS FOR HORIZONTAL ANGLES

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( LENC = 10 )
      CHARACTER*80 BCARD
      CHARACTER*3 ASS
      CHARACTER*1 TC
      LOGICAL GETSSN,ADDCON,GETIVF
      LOGICAL FATAL,LSN
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL REJECT,GET82
      DIMENSION IC(LENC),C(LENC)
      DIMENSION B(*),NX(*)
      COMMON /CONST/  PI, PI2, RAD
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON/UNITS/ LUNIT
      SAVE IYR, IMO, IDY

      READ (BCARD,1) IRT,ISSN,JSSN,ID,IM,ASS,KSSN
 1    FORMAT (7X,I2, 1X,I4,36X,I4, 9X,I3,I2,A3,I4)
      CALL NBLANK (ASS,3,IBLK)
      READ (ASS,3) SS
 3    FORMAT (F3.1)
      IF (IRT .EQ. 30) THEN
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
      ELSEIF (IRT .EQ. 32) THEN
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

      IF (GET82(ISSN,I)) return
      IF (GET82(JSSN,I)) return
      IF (GET82(KSSN,I)) RETURN

        IF ( .NOT. GETSSN(ISSN,ISN)) THEN
          CALL LINE (3)
          WRITE (LUNIT,2) BCARD
 2        FORMAT ('0NO *80* RECORD FOR--',A80/)
        ELSEIF ( .NOT. GETSSN(JSSN,JSN)) THEN
          CALL LINE (3)
          WRITE (LUNIT,2) BCARD
        ELSEIF ( .NOT. GETSSN(KSSN,KSN)) THEN
          CALL LINE (3)
          WRITE (LUNIT,2) BCARD
        ELSE

*** RETRIEVE THE STD DEV

          CALL STDDEV (BCARD,IRT,SD,REJECT)

*** KIND = 10 HORIZONTAL ANGLES

          KIND = 10
          IOBS = IOBS+1
          NOBS = NOBS+1
          NANG = NANG+1
          LSN = .TRUE.
          OBSB = (ID+IM/60.D0+SS/3600.D0)/RAD
          IF ( .NOT. GETIVF(30,ITIME,IVF)) IVF = 0
          IGRT = 0
          CALL FORMIC (KIND,ISN,JSN,KSN,IC,LENG,IGRT)
          CALL FORMC (KIND,C,B,ISN,JSN,KSN,IGRT)
          CALL COMPOB (KIND,OBS0,B,OBSB,ISN,JSN,KSN,IGRT)
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
 667        FORMAT ('0INSUFFICIENT STORAGE FOR HOR. ANGLES'/)
            CALL ABORT2
          ENDIF

          WRITE (IUO) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                IVF,KSN,IGRT

        ENDIF

      RETURN
      END
