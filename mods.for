C %P%
* file $RCSfile: mods.for,v $  v_$Revision: 87172 $

      SUBROUTINE ADJST (B, GOOGE, SHIFTS, G, NX, A, LAWORK, LNWORK,
     &                  IUO, IUO2)
********************************************************************************
* ADJUSTMENT ROUTINE
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*26 TBUF
      CHARACTER*99 SCCSID

      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
** v 4.30vf
      LOGICAL CONVRG, VFCVRG,VFCGPS

      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      DIMENSION B(*), GOOGE(*), SHIFTS(*), G(*), NX(*), A(*)
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON/UNITS/ LUNIT
      LOGICAL NLMODE,NLSCL
      COMMON /NLOPT/ NLMODE,NLSCL,NLUN1,NLUN2,NLUN3,NLUN4,NLUN5
***10-24-03
       common /opt2/ itmax2
******

**     2/7/2003  ** override input covariance matrix  by scaling horizontal and vertical components
** v 4.32vf
       LOGICAL IVCGPS,VSlogic
      COMMON/VCREC/SIGH,SIGU, VFHOR,VFUP, IVCGPS
      CHARACTER*13 VSPROJ
      PARAMETER (MAXPRJ=2500)
      COMMON/VSREC/SIGHS(MAXPRJ),SIGUS(MAXPRJ),ISETHU,VSPROJ(MAXPRJ),
     &       VSlogic

**v6.3

      parameter  (MXSSN = 9999)
      LOGICAL    overwrite_const_coord_in_newbb
      character  CC_records*80
      COMMON/MM6/overwrite_const_coord_in_newbb
      COMMON/MM6_array/CC_records(MXSSN)

**So far for v6.3

C
      EXTERNAL SECOND

      DATA ITYPE / 1 /

C      SCCSID='$Id: mods.for 87172 2016-01-21 15:33:05Z jarir.saleh $	20$Date: 2007/11/30 13:11:56 $ NGS'

*
*** GET STATION CONNECTION MATRIX FOR USE IN COMPUTING LOCAL ACCURACIES
*
      IF(NLMODE) CALL NLLIST(B, G, NX, LNWORK)

*** GET MINUTES ELAPSED SINCE START OF PROGRAM

      CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
      CALL LINE (2)
      WRITE (LUNIT,55) DIFM0
   55 FORMAT (/, ' MINUTES SINCE START OF PROGRAM IS', F7.1)

*** READ DATA AND FORM OBS EQS AND CONNECTION MATRIX
*** REORDER FOR MINIMUM PROFILE

c Mike Potterfield 2/10/06
c Open a new scratch file IUO3 for reject codes for each vector
c This scratch file is written in SECGPS and NEWGRP and
c read in OBSSUM and RGPS.  These are the only uses for it.

      IUO3 = 10
      OPEN (IUO3, FORM='UNFORMATTED', STATUS='SCRATCH')
      REWIND IUO3
      ENDFILE IUO3
      REWIND IUO3


      CALL SECOND (IUO, B, G, NX, LNWORK, IUO3)

*** INITIALIZE SHIFTS

      IF (IMODE .NE. 0) THEN
        IF (LPS) THEN
          DO 30 I = 1, NSTA
            I1 = IUNSHF(I, 1)
            I2 = IUNSHF(I, 2)
            I3 = IUNSHF(I, 3)
            SHIFTS(I1) = 0.D0
            SHIFTS(I2) = 0.D0
            SHIFTS(I3) = 0.D0
   30     CONTINUE
        ENDIF

*** LIST THE OBSERVATIONAL SUMMARY

c Mike Potterfield 2/10/06
c Add the new scratch file IUO3 to OBSSUM
        REWIND IUO3
        CALL OBSSUM (IUO, G,IUO3)
        REWIND IUO3

        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL LINE (2)
        WRITE (LUNIT,65) DIFM0
   65   FORMAT (/, ' TOTAL MINUTES TO COMMENCMENT OF ADJUSTMENT', F7.1)

*** FORM NORMALS AND SOLVE

        CALL HEAD
        CALL LINE (2)
        WRITE (LUNIT,2)
    2   FORMAT ('0', 11('*'), ' COMMENCING ADJUSTMENT ', 11('*') )
      ENDIF

      ITER = 0
  100 CALL NORMAL (IUO, ITER, B, A, G, NX, LAWORK, GOOGE)

*** UPDATE UNKNOWNS AND TEST CONVERGENCE

      IF (IMODE .EQ. 0) THEN
** v 4.30vf
C        CALL FINAL (IUO, A, NX, B, SHIFTS, GOOGE, G, SIGUWT, IUO3)

      ELSEIF ( CONVRG(A, NX, B, SHIFTS, SIGUWT, ITER) ) THEN
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL LINE (2)
        WRITE (LUNIT,75) ITER, DIFM
   75   FORMAT (/, ' MINUTES NEEDED FOR ITERATION', I3, ' =', F7.1)
        CALL LINE (2)
        WRITE (LUNIT,3)
    3   FORMAT ('0*********** ADJUSTMENT CONVERGED *************')
** v 4.30vf
********modified 4/8/2003
*        IF (NVFREE .LE. 0) THEN
*          CALL FINAL (IUO, A, NX, B, SHIFTS, GOOGE, G, SIGUWT)

*        ELSE
          IF (LUPI) THEN
            CALL UPAFIL (B, ITER)
            CALL UPDAT (B)
          ENDIF
*       ENDIF
***********************************
      ELSE
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL LINE (2)
        WRITE (LUNIT,75) ITER, DIFM
        IF (ITER .GE. ITMAX) THEN
          CALL LINE (2)
          WRITE (LUNIT,1)
    1     FORMAT ('0SLOWLY CONVERGING SOLUTION')
          IF (LDR) THEN
            CALL FINAL (IUO, A, NX, B, SHIFTS, GOOGE, G, SIGUWT,
     &          IUO3)
            CALL ABORT2
          ELSE
            CALL ABORT2
          ENDIF
        ELSE
          IF (LUPI) THEN
            CALL UPAFIL (B, ITER)
            CALL UPDAT (B)
          ENDIF
          CALL FORMOB (IUO, IUO2, B, G)
          ITER = ITER + 1
          GO TO 100
        ENDIF
      ENDIF

*** VARIANCE FACTOR ITERATION
*** VFCVRG REFORMS OBS. EQS., COMPUTES INVERSE, AND TESTS V.F. CONVERGE

      IF (NVFREE .GT. 0) THEN
  200   IF (VFCVRG(IUO, IUO2, B, G, A, NX, GOOGE) ) THEN
          CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
          CALL LINE (2)
          WRITE (LUNIT,75) ITER, DIFM
          CALL LINE (2)
          WRITE (LUNIT,5)
    5     FORMAT ('0***** VARIANCE FACTOR ADJUSTMENT CONVERGED *******')
** v 4.30vf
C          CALL FINAL (IUO, A, NX, B, SHIFTS, GOOGE, G, SIGUWT)
******************
        ELSEIF (ITER .GE. ITMAX) THEN
          CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
          CALL LINE (2)
          WRITE (LUNIT,75) ITER, DIFM
          CALL LINE (2)
          WRITE (LUNIT,4)
    4    FORMAT ('0SLOWLY CONVERGING VARIANCE FACTOR SOLUTION ',
     &           50('*') )
          IF (LDR) THEN
            CALL FINAL (IUO, A, NX, B, SHIFTS, GOOGE, G, SIGUWT,
     &         IUO3)
            CALL ABORT2
          ELSE
            CALL ABORT2
          ENDIF
        ELSE
          CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
          CALL LINE (2)
          WRITE (LUNIT,75) ITER, DIFM
          IF (LUPI) THEN
            CALL UPAFIL (B, ITER)
            CALL UPDAT (B)
          ENDIF
          CALL NORMAL (IUO, ITER, B, A, G, NX, LAWORK, GOOGE)
          ITER = ITER + 1
          IF ( CONVRG(A, NX, B, SHIFTS, SIGUWT, ITER) ) THEN
            CONTINUE
          ENDIF
          GO TO 200
        ENDIF
      ENDIF
** v 4.30vf
****ADDED 4/07/2003
*
      IF (IVCGPS ) THEN
  300   IF (VFCGPS(IUO, IUO2, B, G, A, NX, GOOGE) ) THEN
          CALL LINE (2)
          WRITE (LUNIT,305)
  305     FORMAT ('0***** GPS HORIZONTAL/VERTICAL VARIANCE FACTOR',
     &      ' ADJUSTMENT CONVERGED *******')
C          CALL FINAL (IUO, A, NX, B, SHIFTS, GOOGE, G, SIGUWT)

*** 10-24-03
*       ELSEIF (ITER .GE. ITMAX) THEN
        ELSEIF (ITER .GE. itmax2) THEN
**********
          CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
          CALL LINE (2)
          WRITE (LUNIT,75) ITER, DIFM
          CALL LINE (2)
          WRITE (LUNIT,304)
  304 FORMAT ('0SOLUTION FOR GPS HORIZONTAL/VERTICAL VARIANCE FACTORS ',
     &   'HAS NOT CONVERGED AFTER MAXIMUM ITERATIONS. ABORTING*******')
          IF (LDR) THEN
            CALL FINAL (IUO, A, NX, B, SHIFTS, GOOGE, G, SIGUWT, IUO3)
            CALL ABORT2
          ELSE
            CALL ABORT2
          ENDIF
        ELSE
          IF (LUPI) THEN
C
            CALL UPAFIL (B, ITER)
            CALL UPDAT (B)
          ENDIF
C
C   UPDATE GPS COVARIANCE MATRICES
          ITER = ITER + 1
          CALL UPCOV(IUO, IUO2, B, G, A)
          CALL NORMAL (IUO, ITER, B, A, G, NX, LAWORK, GOOGE)
          IF ( CONVRG(A, NX, B, SHIFTS, SIGUWT, ITER) ) THEN
            CONTINUE
          ENDIF
          CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
          CALL LINE (2)
          WRITE (LUNIT,75) ITER, DIFM
          GO TO 300
        ENDIF
      ENDIF
C
      CALL FINAL (IUO, A, NX, B, SHIFTS, GOOGE, G, SIGUWT, IUO3)
**v 4.30vf - end insert  4/07/03********

c Mike Potterfield 2/13/06
c Now nuke IUO3, it's no longer being used
        CLOSE (IUO3)

      RETURN
      END
C---------------------------------------------------------------------------------------------------

      SUBROUTINE AFIL (A)
********************************************************************************
* ROUTINE TO PROCESS THE A FILE
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ACARD
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7 ADJFIL, NAFILE
      CHARACTER*2 CC12
      LOGICAL FATAL
      DIMENSION A(*)
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      COMMON/UNITS/ LUNIT

*** INITIALIZE DEFAULTS

      CALL AFINIT
      FATAL = .FALSE.

*** HEADING

      CALL LINE (3)
      WRITE (LUNIT,1)
    1 FORMAT ('0', 32('*'), ' A-FILE CONTENTS ', 32('*')/)

*** OPEN ADJUSTMENT FILE (IUNIT=11)

***v 4.28k
****      IF (AFILE .EQ. 'NOAFILE') GOTO 200
      
      IUNIT = 11
      OPEN (IUNIT,ERR=200,STATUS='OLD',FILE=AFILE,IOSTAT=IOS)

*** READ A FILE -- PROCESS RECORDS WITHOUT SSN FIELDS

   10 READ (IUNIT,11,END=100) ACARD
   11 FORMAT (A80)
   
      READ (ACARD,14) CC12
   14 FORMAT (A2)
      CALL LINE (1)
      WRITE (LUNIT,12) ACARD
   12 FORMAT (5X, A80)

      IF (CC12 .EQ. 'AA') THEN
        CALL AFAA (ACARD)
      ELSEIF (CC12 .EQ. 'BB') THEN
        CALL AFBB (ACARD)
      ELSEIF (CC12 .EQ. 'DD') THEN
        CALL AFDD (ACARD)
      ELSEIF (CC12 .EQ. 'EE') THEN
        CALL AFEE (ACARD)
      ELSEIF (CC12 .EQ. 'GG') THEN
        CALL AFGG (ACARD)
      ELSEIF (CC12 .EQ. 'HD') THEN
        CALL AFHD (ACARD)
      ELSEIF (CC12 .EQ. 'II') THEN
        CALL AFII (ACARD)
      ELSEIF (CC12 .EQ. 'MM') THEN
        CALL AFMM (ACARD)
	
*************************************    
** v 4.16	
*     ELSEIF (CC12 .EQ. 'NL') THEN
*       CONTINUE
** v 4.26 
      ELSEIF (CC12 .EQ. 'NL') THEN
        CALL AFNL (ACARD)
*************************************

      ELSEIF (CC12 .EQ. 'PP') THEN
        CALL AFPP (ACARD)
      ELSEIF (CC12 .EQ. 'RR') THEN
        CALL AFRR (ACARD, FATAL)
      ELSEIF (CC12 .EQ. 'SS') THEN
        CALL AFSS (ACARD, A, FATAL)
      ELSEIF (CC12 .EQ. 'VV') THEN
        CALL AFVV (ACARD)
** v 4.30vf	
***    added feb 10, 2003
	ELSEIF(CC12 .EQ.'VS') THEN
          CALL AFVS (ACARD)
***     end insert 2/10/2003
      ELSEIF (CC12 .NE. 'CC'  .AND.  CC12 .NE. 'HC'  .AND.
     &        CC12 .NE. 'QQ'  .AND.  CC12 .NE. 'CA'  .AND.
     &        CC12 .NE. 'CD'  .AND.  CC12 .NE. 'CZ'  .AND.
     &        CC12 .NE. 'CH'  .AND.  CC12 .NE. '  '  .AND.
     &        CC12 .NE. '**'  .AND.  CC12 .NE. 'CP') THEN
        CALL LINE (3)
        WRITE (LUNIT,13) CC12
   13   FORMAT ('0*** ILLEGAL ADJUSTMENT FILE FORMAT TYPE -- ', A2/)
      ENDIF
      GO TO 10

*** END OF A FILE PROCESSING

  100 CLOSE (IUNIT)
      CALL EFINIT
      CALL LINE (2)
      WRITE (LUNIT,101)
  101 FORMAT ('0', 33('*'), ' END OF A-FILE ', 33('*') )
      IF (FATAL) then
        CALL ABORT2
      endif
      RETURN

*** A FILE NOT PRESENT

  200 CALL EFINIT
      CALL LINE (3)
***v 4.28k
***   WRITE (LUNIT,201) IOS
***  201 FORMAT ('0', 30('*'), ' NO A-FILE -- ALL DEFAULTS ACTIVE ',
***     &        30('*'), I10/)
      write(lunit,201)
  201 format ('0', 33('*'), ' A-FILE NOT FOUND--FATAL ERROR')     
      call abort2

      RETURN
      END


      SUBROUTINE AFINIT
********************************************************************************
* INITIALIZE THE DEFAULTS
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )

      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP

      LOGICAL L2HLF, LEHT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT

      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
     
*v 4.29      
      COMMON /OPRIN2/ CRIT2

      LOGICAL LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      LOGICAL NLMODE,NLSCL
      COMMON /NLOPT/ NLMODE,NLSCL,NLUN1,NLUN2,NLUN3,NLUN4,NLUN5
      LOGICAL LEB, LLB, LEG, LLG, LED, LLD
      COMMON /PECHO/  LEB, LLB, LEG, LLG, LED, LLD

      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      
      COMMON /YEARS/ ICBEG,IYBEG,ICEND,IYEND
** v 4.32vf
       LOGICAL IVCGPS,VSlogic
      COMMON/VCREC/SIGH,SIGU, VFHOR,VFUP, IVCGPS
      CHARACTER*13 VSPROJ
      PARAMETER (MAXPRJ=2500)
      COMMON/VSREC/SIGHS(MAXPRJ), SIGUS(MAXPRJ), ISETHU, VSPROJ(MAXPRJ),
     &       VSlogic
      
***10-24-03
      common /opt2/ itmax2
*******



      LDIR = .FALSE.
      LANG = .FALSE.
      LZEN = .FALSE.
      LDIS = .FALSE.
      LAZI = .FALSE.
      LGPS = .FALSE.
      LDOP = .FALSE.

      L2HLF = .FALSE.
      LEHT  = .FALSE.

      CRIT = 0.D0
      CRIT2 = 0.D0
      
      LBB  = .TRUE.
      LGF  = .TRUE.
      LCS  = .TRUE.
      LVD  = .TRUE.
      LVA  = .TRUE.
      LVZ  = .TRUE.
      LVS  = .TRUE.
      LVR  = .TRUE.
      LVG  = .TRUE.
      LVC  = .TRUE.
      LIS  = .TRUE.
      LPS  = .TRUE.
      LPG  = .TRUE.
      LDR  = .FALSE.
      LOS  = .TRUE.
      LAP  = .TRUE.
      LDF  = .TRUE.
      LVX  = .TRUE.
      LGV  = .TRUE.

      AX    = 6378137.D0
      E2    = 6.6943800229034156D-3
      DMSL  = 0.D0
      DGH   = 0.D0
      VM    = 300.D0
      VP    = 30.0
      CTOL  = 0.003D0
      ITMAX = 5
      IMODE = 1
      LMSL  = .TRUE.
      LSS   = .FALSE.
      LUP   = .FALSE.
      LABS  = .FALSE.
      LADJ  = .FALSE.
      LUPI  = .FALSE.

      LEB = .FALSE.
      LLB = .FALSE.
      LEG = .FALSE.
      LLG = .FALSE.
      LED = .FALSE.
      LLD = .FALSE.
      NLMODE = .FALSE.
      NLSCL  = .FALSE.

      NSTA  = 0
      NAUX  = 0
      NUNK  = 0
      IDIM  = 3
      NSTAS = 0
      NOBS  = 0
      NCON  = 0
      NZ    = 0
      NGRT  = 0
      
      ICBEG = 0
      IYBEG = 0
      ICEND = 0
      IYEND = 0
** v 4.32vf
      ISETHU=0
      IVCGPS=.FALSE.
      SIGH=1.0
      SIGU=1.0
*************
***10-24-03
      itmax2 = 5
*****
      CALL NEWPRM
      CALL NEWGRT
      CALL NEWIVF

      RETURN
      END


      SUBROUTINE AFPRNT
********************************************************************************
* ECHO THE ADJUSTMENT OPTIONS
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXSSN = 9999 )
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7 ADJFIL, NAFILE
      LOGICAL LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      LOGICAL LEB, LLB, LEG, LLG, LED, LLD
      COMMON /PECHO/  LEB, LLB, LEG, LLG, LED, LLD
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
* v 4.29     
      COMMON /OPRIN2/ CRIT2
     
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      LOGICAL NLMODE,NLSCL
      COMMON /NLOPT/ NLMODE,NLSCL,NLUN1,NLUN2,NLUN3,NLUN4,NLUN5

      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /BYPASS/ LDIR, LANG, LZEN, LDIS, LAZI, LGPS, LDOP
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT
** v 4.32vf , 6.1 
      LOGICAL IVCGPS,VSlogic
      CHARACTER*13 VSPROJ
      COMMON/VCREC/SIGH,SIGU, VFHOR,VFUP, IVCGPS
** v 6.3
      LOGICAL overwrite_const_coord_in_newbb
      COMMON/MM6/overwrite_const_coord_in_newbb 
      
      PARAMETER (MAXPRJ=2500)
      COMMON/VSREC/SIGHS(MAXPRJ), SIGUS(MAXPRJ), ISETHU, VSPROJ(MAXPRJ),
     &       VSlogic

      CALL HEAD
      CALL LINE (46)
      WRITE (LUNIT,7777)
 7777 FORMAT (' ***************** ADJUSTMENT FILE OPTIONS **********'/)

      WRITE (LUNIT,1) AX
    1 FORMAT (' ELLIPSOID SEMI-MAJOR AXIS = ', F12.3, ' METERS')

      WRITE (LUNIT,2) E2
    2 FORMAT (' ECCENTRICITY SQUARED = ', F20.18)

      WRITE (LUNIT,3) DMSL
    
    
** v 4.29f    
*   3 FORMAT (' DEFAULT MEAN SEA LEVEL = ', F8.3, ' METERS')
    3 FORMAT (' DEFAULT ORTHOMETRIC HEIGHT = ', F8.3, ' METERS')
**************************    

      WRITE (LUNIT,4) DGH
    4 FORMAT (' DEFAULT GEOID HEIGHT = ', F8.3, ' METERS')

      IF (LMSL) THEN
        WRITE (LUNIT,55)
** v 4.29f	
*  55   FORMAT (' ADJUST ORTHOMETRIC ELEVATIONS')
   55   FORMAT (' ADJUST ORTHOMETRIC HEIGHTS')
**********************
   
      ELSE
        WRITE (LUNIT,5)
    5   FORMAT (' ADJUST GEOID HEIGHTS')
      ENDIF

      IF (LSS) THEN
        WRITE (LUNIT,155)
  155   FORMAT (' SCALE SIGMAS BY A-POSTERIORI SIGMA OF UNIT WEIGHT')
      ELSE
        WRITE (LUNIT,156)
  156   FORMAT (' DO NOT SCALE SIGMAS BY A-POSTERIORI SIGMA')
      ENDIF

      IF (LABS) THEN
        WRITE (LUNIT,153)
  153   FORMAT (' ABORT IF SINGULARITIES')
      ELSE
        WRITE (LUNIT,154)
  154   FORMAT (' DO NOT ABORT IF SINGULARITIES')
      ENDIF

      IF (LUP) THEN
        WRITE (LUNIT,157) BBNAM
* 157   FORMAT (' UPDATE *80* RECORDS INTO NEW BBOOK FILENAME--', A40)
  157   FORMAT (' UPDATE *80* RECORDS INTO NEW BBOOK FILENAME--', A80)
      ELSE
        WRITE (LUNIT,158)
  158   FORMAT (' DO NOT UPDATE *80* RECORDS')
      ENDIF

**v6.3

      if (overwrite_const_coord_in_newbb) then
        WRITE (LUNIT,159) 
  159   FORMAT (
     &' UPDATE OUTPUT BLUE BOOK WITH CONSTRAINED VALUES FROM A-FILE')
      else
        WRITE (LUNIT,160) 
  160   FORMAT (
     &" DO NOT UPDATE OUTPUT B-BOOK WITH CONSTRAINED VALUES FROM A-FILE"
     &)
      endif
**So far for v6.3
      
      IF (LADJ) THEN
        WRITE (LUNIT,200) ADJFIL
  200   FORMAT (' CREATE ADJUSTED POSITIONS FILE WITH FILENAME--', A7)
      ELSE
        WRITE (LUNIT,201)
  201   FORMAT (' DO NOT CREATE ADJUSTED POSITIONS FILE')
      ENDIF

      IF (LUPI) THEN
** v 4.28
*       WRITE (LUNIT,167)  BBNAM, NAFILE
        WRITE (LUNIT,167)  NAFILE, BBNAM
* 167   FORMAT (' UPDATE BLUE BOOK AND ADJUSTMENT FILE AT THE END OF',
*    &          ' EACH ITERATION INTO FILES ', A7, ' AND ', A7)
* 167   FORMAT (' UPDATE ADJUSTMENT AND BLUE BOOK FILES AT THE END OF',
*    &          ' EACH ITERATION INTO FILES ', A7, ' AND ', A40)
  167   FORMAT (' UPDATE ADJUSTMENT AND BLUE BOOK FILES AT THE END OF',
     &          ' EACH ITERATION INTO FILES ', A7, ' AND ', A80)
************
      ELSE
        WRITE (LUNIT,168)
  168   FORMAT (' DO NOT UPDATE BLUE BOOK AND ADJUSTMENT FILE AT THE',
     &          ' END OF EACH ITERATION')
      ENDIF

      IF ( .NOT. L2HLF) THEN
        WRITE (LUNIT,6) IDIM
    6   FORMAT (' COMPUTE A', I2, '-DIMENSIONAL ADJUSTMENT')
      ELSE
        WRITE (LUNIT,61)
   61   FORMAT (' COMPUTE A 2.5-DIMENSIONAL ADJUSTMENT')
   
*** BIGADJUST ONLY
*       IF (LEHT) THEN
*         WRITE (LUNIT,62)
*  62   FORMAT (' OBTAIN ELLIPSOID HEIGHTS FROM EXTENDED *80* RECORDS')
*       ELSE
*         WRITE (LUNIT,63)
*  63   FORMAT (' OBTAIN ELLIPSOID HEIGHTS FROM ORTHOMETRIC AND GEOID',
*    &          ' HEIGHTS')
*       ENDIF
      ENDIF
      
      WRITE (LUNIT,7) ITMAX
    7 FORMAT (' COMPUTE NO MORE THAN', I5, ' ITERATIONS')

      IF (LDR) THEN
        WRITE (LUNIT,127)
  127   FORMAT (' DISPLAY STATISTICS IF SOLUTION SLOWLY CONVERGES')
      ELSE
        WRITE (LUNIT,128)
  128   FORMAT (' DO NOT DISPLAY STATISTICS IF SOLUTION SLOWLY',
     &          ' CONVERGES')
      ENDIF

      WRITE (LUNIT,8) VM
    8 FORMAT (' ABORT IF MISCLOSURE EXCEEDS', G8.1, ' SIGMA')

      WRITE (LUNIT,81) VP
   81 FORMAT (' PRINT WHEN MISCLOSURES EXCEED', G8.1, ' SIGMA')

      WRITE (LUNIT,9) CTOL
    9 FORMAT (' CONVERGE IF RMS SUM OF SHIFTS BELOW', F8.3, ' METERS')

      IF (IMODE .EQ. 0) THEN
        WRITE (LUNIT,90)
   90   FORMAT (' PERFORM A SIMULATION ***************************')
      ELSEIF (IMODE .EQ. 1) THEN
        WRITE (LUNIT,91)
   91   FORMAT (' COMPUTE QUASI-NORMALIZED RESIDUALS--NO INVERSE')
      ELSEIF (IMODE .EQ. 2) THEN
        WRITE (LUNIT,92)
   92   FORMAT (' COMPUTE QUASI-NORMALIZED RESIDUALS--COMPUTE INVERSE')
      ELSEIF (IMODE .EQ. 3) THEN
        WRITE (LUNIT,93)
   93   FORMAT (' COMPUTE NORMALIZED RESIDUALS AND INVERSE')
      ENDIF
    
      IF (NLMODE ) THEN
        WRITE(LUNIT,97)
   97   FORMAT (' COMPUTE NETWORK AND LOCAL ACCURACIES')
        IF ( IMODE .LE. 1) WRITE (LUNIT, 98)
   98   FORMAT (' THE ABOVE OPTION WILL BE IGNORED - IT REQUIRES ',
     1 'THAT THE MODE INDICATOR ON THE MM RECORD BE GREATER THAN 1')
        IF(NLSCL) THEN
          WRITE(LUNIT,99)
        ELSE 
          WRITE(LUNIT,100)
   99   FORMAT(' NETWORK AND LOCAL ACCURACIES WILL BE SCALED')
  100   FORMAT(' NETWORK AND LOCAL ACCURACIES WILL NOT BE SCALED')
        ENDIF
      ENDIF

      IF (LBB .AND. (.NOT. LEB) .AND. (.NOT. LLB)) WRITE (LUNIT,103)    
  103   FORMAT (' ECHO BLUE BOOK FILE')
      IF (LEB) WRITE (LUNIT,171)
  171   FORMAT (' ECHO BLUE BOOK OBSERVATIONS ONLY')
      IF (LLB) WRITE (LUNIT,172)
  172   FORMAT (' ECHO LARGE BLUE BOOK MISCLOSURES ONLY')
      IF ( .NOT. LBB) WRITE (LUNIT,104)
  104   FORMAT (' DO NOT ECHO BLUE BOOK FILE')

      IF (LDIR) WRITE (LUNIT,181)
  181 FORMAT (' BYPASS ALL HORIZONTAL DIRECTIONS')
      IF (LANG) WRITE (LUNIT,182)
  182 FORMAT (' BYPASS ALL HORIZONTAL ANGLES')
      IF (LZEN) WRITE (LUNIT,183)
  183 FORMAT (' BYPASS ALL ZENITH DISTANCES')
      IF (LDIS) WRITE (LUNIT,184)
  184 FORMAT (' BYPASS ALL DISTANCES')
      IF (LAZI) WRITE (LUNIT,185)
  185 FORMAT (' BYPASS ALL ASTRONOMIC AZIMUTHS')

      IF (LGF  .AND.  ( .NOT. LGPS) ) THEN
        WRITE (LUNIT,101)
  101   FORMAT (' ECHO GPS DATA TRANSFER FILE')
        IF (LEG) WRITE (LUNIT,173)
  173   FORMAT (' ECHO G-FORMAT OBSERVATIONS ONLY')
        IF (LLG) WRITE (LUNIT,174)
  174   FORMAT (' ECHO LARGE G-FORMAT MISCLOSURES ONLY')
      ELSE
        WRITE (LUNIT,102)
  102   FORMAT (' DO NOT ECHO GPS DATA TRANSFER FILE')
      ENDIF

      IF (LGPS) WRITE (LUNIT,186)
  186 FORMAT (' BYPASS ALL GPS DATA')

      IF (LDF  .AND.  ( .NOT. LDOP) ) THEN
c       WRITE (LUNIT,191)
  191   FORMAT (' ECHO DOPPLER DATA TRANSFER FILE')
        IF (LED) WRITE (LUNIT,193)
  193   FORMAT (' ECHO D-FORMAT OBSERVATIONS ONLY')
        IF (LLD) WRITE (LUNIT,194)
  194   FORMAT (' ECHO LARGE D-FORMAT MISCLOSURES ONLY')
      ELSE
c       WRITE (LUNIT,192)
  192   FORMAT (' DO NOT ECHO DOPPLER DATA TRANSFER FILE')
      ENDIF

c     IF (LDOP) WRITE (LUNIT,196)
c 196 FORMAT (' BYPASS ALL DOPPLER DATA')

      IF (LCS) THEN
        WRITE (LUNIT,111)
  111   FORMAT (' DISPLAY CONSTRAINTS')
      ELSE
        WRITE (LUNIT,112)
  112   FORMAT (' DO NOT DISPLAY CONSTRAINTS')
      ENDIF

* v 4.29
      WRITE (LUNIT,131) CRIT2
  131 FORMAT (' DISPLAY ALL GPS RESIDUALS GREATER OR EQUAL TO ',
     *F4.1,' MM')
  
* v 4.29  
      WRITE (LUNIT,1310) CRIT
 1310 FORMAT (' DISPLAY ALL NON-GPS RESIDUALS/SD GREATER OR EQUAL TO ',
     *F4.1,' SIGMA')

      IF (LVD) THEN
        IF ( .NOT. LDIR) THEN
          WRITE (LUNIT,113)
  113     FORMAT (' DISPLAY DIRECTION RESIDUALS')
        ELSE
          WRITE (LUNIT,114)
  114     FORMAT (' DO NOT DISPLAY DIRECTION RESIDUALS')
        ENDIF
      ENDIF

      IF (LVA) THEN
        IF ( .NOT. LANG) THEN
          WRITE (LUNIT,115)
  115     FORMAT (' DISPLAY ANGLE RESIDUALS')
        ELSE
          WRITE (LUNIT,116)
  116     FORMAT (' DO NOT DISPLAY ANGLE RESIDUALS')
        ENDIF
      ENDIF

      IF (LVZ) THEN
        IF ( .NOT. LZEN) THEN
          WRITE (LUNIT,117)
  117     FORMAT (' DISPLAY ZENITH DISTANCE RESIDUALS')
        ELSE
          WRITE (LUNIT,118)
  118     FORMAT (' DO NOT DISPLAY ZENITH DISTANCE RESIDUALS')
        ENDIF
      ENDIF

      IF (LVS) THEN
        IF ( .NOT. LDIS) THEN
          WRITE (LUNIT,119)
  119     FORMAT (' DISPLAY DISTANCE RESIDUALS')
        ELSE
          WRITE (LUNIT,120)
  120     FORMAT (' DO NOT DISPLAY DISTANCE RESIDUALS')
        ENDIF
      ENDIF

      IF (LVR) THEN
        IF ( .NOT. LAZI) THEN
          WRITE (LUNIT,121)
  121     FORMAT (' DISPLAY ASTRO-AZIMUTH RESIDUALS')
        ELSE
          WRITE (LUNIT,122)
  122     FORMAT (' DO NOT DISPLAY ASTRO-AZIMUTH RESIDUALS')
        ENDIF
      ENDIF

      IF (LVG) THEN
        IF ( .NOT. LGPS) THEN
          WRITE (LUNIT,123)
  123     FORMAT (' DISPLAY GPS RESIDUALS')
        ELSE
          WRITE (LUNIT,124)
  124     FORMAT (' DO NOT DISPLAY GPS RESIDUALS')
        ENDIF
      ENDIF

      IF (LVX) THEN
        IF ( .NOT. LDOP) THEN
          WRITE (LUNIT,197)
  197     FORMAT (' DISPLAY DOPPLER RESIDUALS')
        ELSE
c         WRITE (LUNIT,198)
  198     FORMAT (' DO NOT DISPLAY DOPPLER RESIDUALS')
        ENDIF
      ENDIF

      IF (LVC) THEN
        WRITE (LUNIT,125)
  125   FORMAT (' DISPLAY CONSTRAINED RESIDUALS')
      ELSE
        WRITE (LUNIT,126)
  126   FORMAT (' DO NOT DISPLAY CONSTRAINED RESIDUALS')
      ENDIF

      IF (LIS) THEN
        WRITE (LUNIT,129)
  129   FORMAT (' DISPLAY RESIDUALS GROUPED AROUND INTERSECTION STAS')
      ENDIF

      IF (LPS) THEN
        WRITE (LUNIT,107)
  107   FORMAT (' DISPLAY POSITION SHIFTS')
      ELSE
        WRITE (LUNIT,108)
  108   FORMAT (' DO NOT DISPLAY POSITION SHIFTS')
      ENDIF

      IF (LPG) THEN
        WRITE (LUNIT,109)
  109   FORMAT (' DISPLAY POSITION GOOGE NUMBERS')
      ELSE
        WRITE (LUNIT,110)
  110   FORMAT (' DO NOT DISPLAY POSITION GOOGE NUMBERS')
      ENDIF

      IF ( .NOT. LOS) WRITE (LUNIT,301)
  301 FORMAT (' DISPLAY ONLY AN ABBREVIATED OBSERVATIONAL SUMMARY')
      IF ( .NOT. LAP) WRITE (LUNIT,302)
  302 FORMAT (' DO NOT DISPLAY THE ADJUSTED POSITIONS')
** v 4.32vf
*** inserted 2/10/3003
	IF(ISETHU.GT.0) THEN
        WRITE(LUNIT,310) 
  310	FORMAT(' GPS SESSION COVARIANCE MATRICES INPUT IN THE G-FILE',
     &' WILL BE SCALED AS FOLLOWS:'/
     &'  PROJECT    HORIZONTAL  UP')
      DO 312 I=1,ISETHU
      WRITE(LUNIT,311) VSPROJ(I), SIGHS(I), SIGUS(I)
  311 FORMAT(1X,A13,2F7.3)
  312 CONTINUE
      WRITE(LUNIT,313)
  313 FORMAT()
      ENDIF
   
      IF(IVCGPS) WRITE(LUNIT,315)
  315 FORMAT(' VARIANCE FACTORS FOR GPS HORIZONTAL AND VERTICAL ',
     &  'COMPONENTS WILL BE COMPUTED')
***     end insert		

      RETURN
      END


      SUBROUTINE FINAL (IUO,A,NX,B,SHIFTS,GOOGE,G,SIGUWT,IUO3)
********************************************************************************
* LIST THE ADJUSTMENT RESULTS
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( MXFOT = 3500, MXSSN = 9999 )
      CHARACTER*26 TBUF
      LOGICAL INVERT,GETA
      LOGICAL LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &        LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL L2HLF, LEHT
      DIMENSION A(*),NX(*),B(*),SHIFTS(*),GOOGE(*),G(*)
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      COMMON /OPRINT/ CRIT, LBB, LGF, LCS, LVD, LVA, LVZ, LVS, LVR, LVG,
     &                LVC, LIS, LPS, LPG, LDR, LOS, LAP, LDF, LVX, LGV
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /NUMVFS/ NVFTOT, NVFREE
      COMMON /STATCT/ N84, N85, N86, NDIR, NANG, NGPS, NZD, NDS, NAZ,
     &                NQQ, NREJ, NGPSR, NDOP
      COMMON /FOURTH/ NFOT, NFOTS(MXFOT)
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      LOGICAL NLMODE,NLSCL
      COMMON /NLOPT/ NLMODE,NLSCL,NLUN1,NLUN2,NLUN3,NLUN4,NLUN5
      COMMON/UNITS/ LUNIT
** v 4.32vf
       LOGICAL IVCGPS,VSlogic
      COMMON/VCREC/SIGH,SIGU, VFHOR,VFUP, IVCGPS
      CHARACTER*13 VSPROJ
      PARAMETER (MAXPRJ=2500)
      COMMON/VSREC/SIGHS(MAXPRJ), SIGUS(MAXPRJ), ISETHU, VSPROJ(MAXPRJ),
     &             VSlogic
**v6.3

      integer*4  CC_counter
      LOGICAL    overwrite_const_coord_in_newbb
      character  CC_records*80
      COMMON/MM6/overwrite_const_coord_in_newbb
      COMMON/MM6_array/CC_records(MXSSN)

**So far for v6.3

****
      DATA  ITYPE / 1 /

***  COMPLETE THE COMPUTATION OF THE GOOGE NUMBERS, IF VARIANCE
***  FACTORS THEN THE GOOGE NUMBS ALREADY COMPUTED IN SUBR. VFCVRG

** v 4.30vf 
      IF (NVFREE .EQ. 0 .AND. .NOT.IVCGPS) THEN
******
*** WRITE OUT THE VERY SMALL GOOGE NUMBERS

        GTOL = 1.D-3
        CALL LINE (3)
        WRITE (LUNIT,101) GTOL
  101   FORMAT (//, ' THE FOLLOWING UNKNOWNS HAVE GOOGE NUMBERS',
     &              ' LESS THAN', 1PD9.2)
        IF ( .NOT. LOS) THEN
          CALL LINE (5)
          WRITE (LUNIT, 201)
  201     FORMAT (/, ' *** STATION CONNECTIONS FOR SMALL GOOGE',
     &               ' NUMBERS ***',
     &          /, 9X, 'SSN',
     &             7X, 'DIR', 9X, 'ANG', 9X, 'AZI',
     &             9X, 'DIS', 9X, 'ZD', 10X, 'GPS', 7X, 'DOP',
     &          /, 9X, 'CMP',
     &             5X, 'FRM   TO', 4X, 'FRM   TO', 4X, 'FRM   TO',
     &             4X, 'FRM   TO', 4X, 'FRM   TO', 4X, 'FRM   TO', /)
        ENDIF

        DO 4 I = 1,NUNK
          IF ( .NOT. GETA(I,I,VAL,A,NX)) THEN
            CALL INVIUN (I,ISN,IUTYPE,IUCODE)
            WRITE (LUNIT,5) I,ISN,IUTYPE,IUCODE
    5       FORMAT ('0 FATAL ERROR IN GOOGE COMPUTATION IN FINAL',4I8)
            CALL ABORT2
          ENDIF

***  ALLOW FOR GOOGE OF ZERO IF SOLVING A SINGULAR SYSTEM

          VAL2 = VAL*VAL
          GI = DIVID( VAL2, GOOGE(I) )
          GOOGE(I) = GI
          IF ( GI  .LT.  GTOL  .AND.  GI .GT. 0.D0 ) THEN
            CALL INVIUN (I,ISN,IUTYPE,IUCODE)
            CALL TINYG (ISN, IUTYPE, IUCODE, GI )
          ENDIF
    4   CONTINUE
      ENDIF

*** INVERT WITHIN PROFILE

** v 4.30vf
      IF (IMODE .NE. 1  .AND.  NVFREE .EQ. 0 .AND. .NOT.IVCGPS) THEN
********
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        IF ( .NOT. INVERT(A,NX)) THEN
          WRITE (LUNIT,666)
  666     FORMAT ('0STATE ERROR IN FINAL')
          CALL ABORT2
        ENDIF
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL LINE (2)
        WRITE (LUNIT,15) DIFM
   15   FORMAT (/,' MINUTES TO INVERT', F7.1)
      ENDIF

*** LIST THE JOB STATISTICS

      CALL JOBSTT

*** LIST ADJUSTED OBSERVATIONS AND RESIDUALS

      CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)

c Mike Potterfield 2/13/06
c Pass IUO3 to RGPS to identify project ids and rejects
      CALL RESID (IUO,A,B,NX,G,SIGUWT,IUO3)

      CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
      CALL LINE (2)
      WRITE (LUNIT,25) DIFM
   25 FORMAT (/,' MINUTES TO LIST RESIDUALS', F7.1)

*** LIST ADJUSTED AUXILIARY PARAMETERS

      IF (NAUX .GT. 0) THEN
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL ADJAUX (A,NX,GOOGE,B,SIGUWT)
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL LINE (2)
        WRITE (LUNIT,35) DIFM
   35   FORMAT (/,' MINUTES TO LIST AUXILIARY PARAMETERS', F7.1)
      ENDIF

*** LIST ADJUSTED AUXILIARY GPS AND DOPPLER ROTATION PARAMETERS

      IF (NGRT .GT. 0) THEN
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL ADJGRT (A,NX,GOOGE,SIGUWT)
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL LINE (2)
        WRITE (LUNIT,45) DIFM
   45   FORMAT (/,' MINUTES TO LIST ROTATION PARAMETERS', F7.1)
      ENDIF

*** LIST THE RESIDUALS GROUPED AROUND INTERSECTION STATIONS

      IF (LIS  .AND.  IMODE  .NE.  0  .AND.  NFOT .GT. 0) THEN
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL RESID2 (IUO,A,B,NX,G,SIGUWT)
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL LINE (2)
        WRITE (LUNIT,55) DIFM
   55   FORMAT (/,' MINUTES TO LIST RESIDUALS GROUPED AROUND',
     &            ' INTERSECTION STATIONS', F7.1)
      ENDIF

*** LIST ADJUSTED POSITIONS AND SHIFTS

      CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
      CALL ADJPOS (A,B,NX,SHIFTS,GOOGE,SIGUWT)
      CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
      CALL LINE (2)
      WRITE (LUNIT,65) DIFM
   65 FORMAT (/,' MINUTES TO LIST ADJUSTED POSITIONS', F7.1)

*** LIST SHIFTS FOR CONSTRAINED STATIONS ONLY

**v6.3

      call constrained_stations_shifts (A, B, NX)

**v6.3 so far for v6.3

*** LIST DUAL HEIGHT DIFFERENCES

      IF (L2HLF) THEN
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL DUALHT (A,B,NX,SHIFTS,GOOGE,SIGUWT)
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL LINE (2)
        WRITE (LUNIT,70) DIFM
   70   FORMAT (/,' MINUTES TO LIST DUAL HEIGHTS', F7.1)
      ENDIF

*** COMPUTE LOCAL AND NETWORK ACCURACIES
      IF (NLMODE) CALL NLACUR(A,B,NX,SHIFTS,GOOGE,SIGUWT)
*** COMPUTE ACCURACIES

      IF (NQQ .GT. 0) THEN
        IF (LPS) THEN
          CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
          CALL ACCUR (A,NX,B,SHIFTS,SIGUWT)
          CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
          CALL LINE (2)
          WRITE (LUNIT,75) DIFM
   75     FORMAT (/,' MINUTES TO COMPUTE ACCURACIES', F7.1)
        ELSE
          WRITE (LUNIT, 175)
  175     FORMAT (/, ' ****** ERROR - ACCURACIES CANNOT BE COMPUTED',
     +               ' UNLESS SHIFTS ARE ENABLED')
          CALL ABORT2
        ENDIF
      ENDIF

*** UPDATE CONTROL POINT RECORDS

      IF (LUP) THEN
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL UPDAT (B)
        CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
        CALL LINE (2)
        WRITE (LUNIT,85) DIFM
   85   FORMAT (/,' MINUTES TO UPDATE CONTROL POINTS', F7.1)
      ENDIF

      RETURN
      END


      SUBROUTINE HEAD
********************************************************************************
* GO TO A NEW PAGE AND PRINT A HEADING
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE
      COMMON /VERSN/  PGMVER
      COMMON /UNITS/  LUNIT

      WRITE (LUNIT,1) IPAGE, PGMVER
    1 FORMAT (/,'1', T54, 'NATIONAL GEODETIC SURVEY',
     &          /, ' PROGRAM ADJUST', T56, ' ADJUSTMENT PROGRAM',
     &             T122, 'PAGE', I4,
     &          /, T59, 'VERSION   ',A5, //)

*   1 FORMAT ('1', T54, 'NATIONAL GEODETIC SURVEY',
*    &          /, ' PROGRAM ADJUST', T56, ' ADJUSTMENT PROGRAM',
*    &             T122, 'PAGE', I4,

      IPAGE = IPAGE+1
      ILINE = 5

      RETURN
      END
C--------------------------------------------------------------------------------------
      SUBROUTINE UPDAT (B)

* UPDATE CONTROL PT RECORDS

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      DIMENSION B(*)
      CHARACTER*80 BCARD
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, BBNAM
      CHARACTER*7  ADJFIL, NAFILE
      CHARACTER*2  BBID,JCODE
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /OPT/    AX, E2, DMSL, DGH, VM, VP, CTOL, ITMAX, IMODE,
     &                LMSL, LSS, LUP, LABS, LADJ, LUPI
      LOGICAL NLMODE,NLSCL
      COMMON /NLOPT/ NLMODE,NLSCL,NLUN1,NLUN2,NLUN3,NLUN4,NLUN5
      COMMON/UNITS/ LUNIT

*** 10-31-03
*       common/vcrec2/vtvhor, vtvup, vtvcor, rnhor, rnup
** v 4.32vf
      LOGICAL IVCGPS,VSlogic
      COMMON/VCREC/SIGH,SIGU, VFHOR,VFUP, IVCGPS
      CHARACTER*13 VSPROJ
      PARAMETER (MAXPRJ=2500)
      COMMON/VSREC/SIGHS(MAXPRJ), SIGUS(MAXPRJ), ISETHU, VSPROJ(MAXPRJ),
     &             VSlogic

**v6.3

      parameter  (MXSSN = 9999)
      LOGICAL    overwrite_const_coord_in_newbb
      LOGICAL    un_constrained_Hz,un_constrained_HT
      character  CC_records*80,lat_lon_char*25,HT_char*7,HT_type*1
      COMMON/MM6/overwrite_const_coord_in_newbb
      COMMON/MM6_array/CC_records(MXSSN)

**So far for v6.3


**************

*** OPEN OLD BLUE BOOK

      IOLD = 14
      OPEN (IOLD,ERR=666,STATUS='OLD',FILE=BBOOK)
*** GET THE JOB CODE
      READ(IOLD,1,END=668) BCARD
      READ(BCARD,4) JCODE

*** OPEN NEW BLUE BOOK

      INEW = 15
      OPEN (INEW,ERR=667,STATUS='UNKNOWN',FILE=BBNAM)
      WRITE(INEW,2) BCARD

*** READ THE REST OF THE BLUE BOOK RECORDS

  100 READ (IOLD,1,END=668) BCARD
    1 FORMAT (A80)
      READ (BCARD,4) BBID
    4 FORMAT (7X,A2)

      IF (BBID .EQ. '80') THEN
C        BCARD(79:80) = 'BA'                  !correct ellipsoidal height order and class
        CALL UP80 (BCARD, B)
      ELSEIF (BBID .EQ. '84') THEN
*** v 4.29k
*       IF (.NOT. LMSL) CALL UP84 (BCARD, B)
        CONTINUE   
***	
      ELSEIF (BBID .EQ. '86') THEN
C        BCARD(54:55) = '41'                  !Correct ellipsoidal height order and class
        CALL UP86 (BCARD, B)
      ELSEIF (BBID .EQ. JCODE) THEN
        IF (NLMODE) CALL UPNL(INEW)
***10-31-03
*            if(ivcgps) then
*               write(inew,200) vtvhor,rnhor,vfhor,sigh,vtvup,rnup,
*     *           vfup,sigu,vtvcor
* 200           format(6x,'*93*',2f8.3,2f7.3,2f8.3,2f7.3,f8.3)
*            endif
************
***12-04-03*********
        if(ivcgps) then
          write(inew,200) sigh,sigu
 200      format(6x,'*93*',2f8.3)
        endif
***********************

        if(VSlogic) then
          do i=1,ISETHU
            write(inew,201) sighs(i),sigus(i),vsproj(i)
 201        format(6x,'*93*',2f8.3,1x,a13)
          enddo
        endif
      ENDIF

C  Changes to the newbb file

      if (BBID == '80')  BCARD(70:76) = '       '           !I added this
      if(BBID/='91'.and.BBID/='92'.and.BBID/='93') WRITE(INEW,2) BCARD
    2 FORMAT (a)
      GO TO 100

  668 CLOSE (IOLD)
      CLOSE (INEW)
      CALL LINE (2)
      WRITE (LUNIT,'(//)') 
      WRITE (LUNIT,5) BBNAM
*    5 FORMAT ('0UPDATED CONTROL POINT RECORDS IN FILE -- ', A40)
    5 FORMAT ('0UPDATED CONTROL POINT RECORDS IN FILE -- ', A80)
    
**v6.3

C  Now that the output bfile is completed, rewind it and overwrite it with the constrained values
C  for constrained points. Keep it the same for unconstrained points.

      if (overwrite_const_coord_in_newbb) then

*** OPEN NEW BLUE BOOK again

        INEW = 15
        OPEN (INEW,ERR=667,STATUS='UNKNOWN',FILE=BBNAM)
        ITMP = 14
        open (ITMP,ERR=667,STATUS='UNKNOWN',FILE='TMP')

C  Read NEWBB records one by one

  101   READ (INEW,'(a80)',END=6681) BCARD
        READ (BCARD,4) BBID

        IF (BBID .EQ. '80') THEN
          READ (BCARD,'(10x,i4)') ISSN
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
          write (ITMP,'(a80)') BCARD
        ELSEIF (BBID .EQ. '86') THEN
          READ (BCARD,'(10x,i4)') ISSN
          do ii=1,MXSSN
            read (CC_records(ii)(11:14),'(i4)') iissn
            if (ISSN == iissn) then
              read(CC_records(ii)(70:76),'(a7)') HT_char
              read(CC_records(ii)(77:77),'(a1)') HT_type
              un_constrained_HT = (HT_char == ' ')
              if((.not. un_constrained_HT) .and. HT_type== "E") then          !The station is constrained includin the height
                BCARD(46:52) = CC_records(ii)(70:76)
                exit
              elseif((.not. un_constrained_HT).and.HT_type==" ") then          !The station is constrained includin the height
                BCARD(17:23) = CC_records(ii)(70:76)
                exit
              endif
            endif
          enddo
          write (ITMP,'(a80)') BCARD
        ELSE
          write (ITMP,'(a80)') BCARD
        ENDIF
 
        GO TO 101

 6681   CLOSE (INEW)
        CLOSE (ITMP)

C  Now just overwrite INEW by ITMP

        INEW = 15
        OPEN (INEW,ERR=667,STATUS='UNKNOWN',FILE=BBNAM)
        ITMP = 14
        open (ITMP,ERR=667,STATUS='UNKNOWN',FILE='TMP')

c  Read TMP records one by one and overwrite NEWBB

  102   READ (ITMP,'(a80)',END=6682) BCARD
        write(INEW,'(a80)'         ) BCARD
        goto 102
 6682   CLOSE (INEW)
c       CLOSE (ITMP)
        CLOSE (ITMP,status="DELETE")
      endif

***so far for v6.3

      RETURN

*** NO OLD BLUE BOOK FOUND

  666 WRITE (LUNIT,699)
  699 FORMAT ('0NO OLD BLUE BOOK FOUND', /)
      CALL ABORT2
      RETURN

*** NOT ABLE TO OPEN NEW BB FILE

  667 WRITE (LUNIT,698)
  698 FORMAT ('0NOT ABLE TO OPEN NEW BLUE BOOK FILE', /)
      CALL ABORT2

      RETURN
      END

