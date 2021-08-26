C %P%

      SUBROUTINE UPCOV (IUO,IUO2,B,G,A)
********************************************************************************
* RECOMPUTE COVARIANCE MATRICES FOR GPS OBSERVATIONS
********************************************************************************
      PARAMETER ( NVECS = 700, LENC = 10, MXVF = 40, MXSSN = 9999 )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL LSS, LUP, LMSL, LABS, LADJ, LUPI
      LOGICAL FATAL,FATAL2,FIXVF,PROP,INVERT
      LOGICAL L2HLF, LEHT
      DIMENSION B(*),G(*),A(*)
	DIMENSION COVECF(3,3)
      DIMENSION IC(LENC), C(LENC)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      COMMON/UNITS/ LUNIT

** v 4.32vf
      LOGICAL IVCGPS,VSlogic
      COMMON/VCREC/SIGH,SIGU, VFHOR,VFUP, IVCGPS
      CHARACTER*13 VSPROJ
      PARAMETER (MAXPRJ=2500)
      COMMON/VSREC/SIGHS(MAXPRJ), SIGUS(MAXPRJ), ISETHU, VSPROJ(MAXPRJ),
     &             VSlogic
      CHARACTER*99 SCCSID
C
C
C      SCCSID='$Id: upcov.for 65375 2012-09-20 15:10:51Z jarir.saleh $	20$Date: 2008/08/25 16:04:06 $ NGS'

      WRITE (LUNIT, 6010) SIGH, SIGU
 6010 FORMAT(/' RESCALING UNCERTAINTIES OF GPS HORIZONTAL OBSERVATIONS',        
     & ' BY', F10.3/' AND                              VERTICAL',
     & ' OBSERVATIONS BY', F10.3)

  100 READ (IUO,END=777) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                   IVF,IAUX,IGRT
      IF (KIND .LE. 999) THEN
        IF (KIND .LT. 18  .OR.  KIND .GT. 20) THEN
          WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                 IVF,IAUX,IGRT

        ELSEIF (KIND .GE. 18 .AND. KIND.LE. 20) THEN
*** CORRELATED TYPE (DOPPLER)
          WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &                 IVF,IAUX,IGRT
          READ (IUO,END=777) ISIGNL, ( (COVECF(I,J), J = 1,3), I = 1,3 )
          WRITE (IUO2) ISIGNL, ( (COVECF(I,J), J = 1,3), I = 1,3 )
        ENDIF
C
      ELSE
*** GPS OBSERVATIONS
        NVEC = ISN
        IAUX = JSN
        NR = 3*NVEC
        WRITE (IUO2) KIND,ISN,JSN,IC,C,LENG,CMO,OBSB,SD,IOBS,
     &               IVF,IAUX,IGRT
C
      IF ( .NOT. L2HLF) THEN
        LENG = (NVEC+1)*IDIM + 1
      ELSE
        LENG = (NVEC+1)*3 + 1
      ENDIF
      IF (NGRT .GT. 0) LENG = LENG + 3

*     NC = NR + 3 + LENG
      NC = NR + 5 + LENG

	  CALL UPDATG (IUO, IUO2, NVEC, NR, NC, LENG, G, B, A)	 
      ENDIF
      GO TO 100

*** ABORT DUE TO LARGE MISCLOSURES

  777 CONTINUE
      REWIND IUO
      REWIND IUO2

*** EXCHANGE PRIMARY/SECONDARY OBS EQ FILE INDICATOR

      ITEMP = IUO
      IUO = IUO2
      IUO2 = ITEMP

*** HEADING

      CALL LINE (3)

      RETURN
      END
      SUBROUTINE UPDATG (IUO, IUO2, NVEC, NR, NC, LENG, G, B, A)

*** ROUTINE TO RE-COMPUTE OBS.EQ. FOR A GPS SESSION

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( NVECS = 700, MXSSN = 9999 )
      DIMENSION B(*), G(NR,NC), A (*)
      COMMON /STRUCT/ NSTA, NAUX, NUNK, IDIM, NSTAS, NOBS, NCON, NZ,NGRT
      COMMON /DUAL/   N2HLF, I2HLF(MXSSN), L2HLF, LEHT
      DIMENSION ICM((NVECS+1)*3+1+3)
      DIMENSION KINDS(NVECS*3), ISNS(NVECS*3), JSNS(NVECS*3)
      DIMENSION LOBS(NVECS*3)
** v 4.32vf
       LOGICAL IVCGPS,VSlogic
      COMMON/VCREC/SIGH,SIGU, VFHOR,VFUP, IVCGPS
      CHARACTER*13 VSPROJ
      PARAMETER (MAXPRJ=2500)
      COMMON/VSREC/SIGHS(MAXPRJ), SIGUS(MAXPRJ), ISETHU, VSPROJ(MAXPRJ),
     &             VSlogic
      COMMON/UNITS/ LUNIT

      CALL GLOCAT(N1,N2,N3,N4,N5,NVEC,LENG)

*** READ SUPPORTING INDICIES

      READ (IUO) ICM, NICM, KINDS, ISNS, JSNS, LOBS
      WRITE (IUO2) ICM, NICM, KINDS, ISNS, JSNS, LOBS

*** LOAD WORK SPACE (G)

      DO 1 I = 1, NR
        READ (IUO) (G(I,J), J = 1, NC)
    1 CONTINUE
C
C   RECONSTRUCT ORIGINAL COVARIANCE MATRIX OF THE OBSERVATIONS
C    BY USING STANDARD DEVIATIONS IN COLUMN N5 OF G
C    AND COVARIANCES BELOW THE DIAGONAL
C
      DO 20 I=1,NR
        DO 10 J=I,NR
          IF (I.EQ.J) GO TO 10
            G(I,J)=G(J,I)
   10   CONTINUE
        G(I,I)=G(I,N5)**2
   20 CONTINUE
C
C   THEN SCALE IT BY NEW HORIZONTAL AND VERTICAL FACTORS
C
      CALL SCALEG (G, NR, NC, LENG,N1,N2,N3,N4,N5, NVEC, B,
     &             KINDS, ISNS, JSNS, VFHOR, VFUP)

      DO 95 I=1,NR
C   MARK REJECTED OBSERVATIONS
      IF(G(I,N5).EQ.100.) THEN
        DO 93 J=1,NR
          G(J,I)=0.0
   93     G(I,J)=0.0
        G(I,I)=10000.
      ELSE
C  SAVE THE DIAGONAL TERM IN COLUMN N5
C
        G(I,N5)=SQRT(G(I,I))
      ENDIF
   95 CONTINUE

******************************************
c 7/17/05 - 7/25/05 mike potterfield 
c try to update the correlation coefficient here
      DO 301 J = 1,NVEC
         I = 3 * (J-1) + 1
c         G(I,N5) = G(I,I)
c         G(I+1,N5) = G(I+1,I+1)
c         G(I+2,N5) = G(I+2,I+2)
c next is correlation XY         
         IF (G(I,I) .EQ. 10000.) THEN
         G(I,N5+1 ) = 0.
         ELSE
         G(I,N5+1) = G(I,I+1)/(SQRT(G(I,I))*SQRT(G(I+1,I+1)))
         ENDIF
c next is correlation YZ         
         IF (G(I+1,I+1) .EQ. 10000.) THEN
         G(I+1,N5+1) = 0.
         ELSE
         G(I+1,N5+1) = G(I+1,I+2)/(SQRT(G(I+1,I+1))*SQRT(G(I+2,I+2)))
         ENDIF
c next is correlation XZ          
         IF (G(I+2,I+2) .EQ. 10000.) THEN
         G(I+2,N5+1) = 0.
         ELSE
         G(I+2,N5+1) = G(I,I+2)/(SQRT(G(I,I))*SQRT(G(I+2,I+2)))
         ENDIF
  301 CONTINUE   

************************************************
************   end insert for overridding input covariance matrix

*** CHOLESKY FACTOR COVARIANCE MATRIX IN WORK SPACE

      DO 101 I = 1, NR
        I1 = I - 1
        IF (I1 .GT. 0) THEN
          DO 102 K = 1, I1
            G(I,I) = G(I,I) - G(K,I)*G(K,I)
  102     CONTINUE
        ENDIF

*** SINGULARITY TEST

****************************3-27
*       WRITE  (LUNIT,999) G(I,I),I
*999    FORMAT ( ' G= ',1PD9.1,' ROW/COL= ',I5)
**********************************

        IF (G(I,I) .LT. 1.D-10) THEN
          WRITE (LUNIT,7) I
    7     FORMAT ('0CORRELATION MATRIX SINGULARITY IN UPDATG--',I5)
          CALL ABORT2
        ENDIF

*** RESUME CHOLESKY FACTOR

        G(I,I) = DSQRT( G(I,I) )
        IF (I+1 .LE. NR) THEN
          DO 106 J = I+1, NR
            IF (I1 .GT. 0) THEN
              DO 104 K = 1, I1
                G(I,J) = G(I,J) - G(K,I)*G(K,J)
  104         CONTINUE
            ENDIF
            G(I,J) = G(I,J)/G(I,I)
  106     CONTINUE
        ENDIF
  101 CONTINUE

***  DECORELLATE COEFFICIENTS OF OBS EQ IN PREPARATION FOR NEXT ITERATION
C
      DO 30 I = 1, LENG
        CALL CMRHS2(G,NR,N2+I)
   30 CONTINUE
C
C
******* ALSO  DECORRELATE MISCLOSURES
      CALL CMRHS2(G,NR,N4)
C
*** REWRITE G
      DO 500 I = 1, NR
        WRITE (IUO2) (G(I,J), J = 1, NC)
  500 CONTINUE
C
      RETURN
      END
