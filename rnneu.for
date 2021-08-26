C %P%

      SUBROUTINE RNNEU (G, NR, NC, NVEC, LENG, B, ICM, NICM, KINDS,
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
      COMMON /GPSHUS/ GPSVS(NVECS*3,2), GPSRNS(NVECS*3)
      COMMON/UNITS/ LUNIT
      CHARACTER*99 SCCSID

****  2/20/2003
      DIMENSION WORK1(3,3),WORK2(3,3),RM(3,3)
*********

C      SCCSID='$Id: rnneu.for 44197 2010-07-13 12:26:03Z bruce.tran $	20$Date: 2008/08/25 16:03:29 $ NGS'

*** LENGTH OF COEFF. ARRAYS

      IF ( .NOT. L2HLF) THEN
        LENGC = 2*IDIM
      ELSE
        LENGC = 2*3
      ENDIF
      IF (IAUX .GT. 0) LENGC = LENGC + 1
      IF (IGRT .GT. 0) LENGC = LENGC + 3


*** DETERMINE WORK ARRAY COLUMN POINTERS

      CALL GLOCAT (N1, N2, N3, N4, N5, NVEC, LENG)

      NB = NUNK + 1
      NL = NR*(NR+1)/2
      NB1 = NX(NB)
      NB2 = NB1 + NL
      NB3 = NB2 + NL


*** (STEP 2.)  COMPUTE VARIANCES OF ADJUSTED OBSERVATIONS
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
   60 CONTINUE

*** (STEP 3.1)  GET COVARIANCE TERMS OF ADJUSTED OBSERVATIONS

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
****       TO GET WEIGHT MATRIX 

      CALL AVERT (A(NB2+1), A(NB3+1), NR)


*** 2/20/2003  STEP 8 - COMPUTE REDUNDANCY NUMBERS FOR E,N,U OBSERVATIONS
      DO 250 IVEC=1,NVEC
        ISN=ISNS(3*IVEC)
        JSN=JSNS(3*IVEC)
        CALL GETGLA (GLA1, ISN, B)
        CALL GETGLO (GLO1, ISN, B)
        CALL GETGLA (GLA2, JSN, B)
        CALL GETGLO (GLO2, JSN, B)

*** COMPUTE MEAN ORIENTATION
        GLA = (GLA1 + GLA2)/2.D0
        GLO = (GLO1 + GLO2)/2.D0

        SB = DSIN(GLA)
        CB = DCOS(GLA)
        SL = DSIN(GLO)
        CL = DCOS(GLO)

*** PUT TOGETHER ROTATION MATRIX --- CMP
        RM(1,1)=-SB*CL
        RM(1,2)=-SL
        RM(1,3)=CB*CL
        RM(2,1)=-SB*SL
        RM(2,2)=CL 
        RM(2,3)=CB*SL
        RM(3,1)=CB
        RM(3,2)=0 
        RM(3,3)=SB
C
C   ROTATE COVARIANCE MATRIX OF THE ADJUSTED OBSERVATIONS 
C         (WHOSE UPPER TRIANGULAR PART IS STORED IN A(NB1)
C   AND THE COVARIANCE MATRIX OF THE OBSERVATIONS
C         (WHOSE UPPER TRIANGULAR PART IS STORED IN A(NB2)
C   INTO NEU COORDINATE SYSTEM (AVERAGE LOCAL NEU FOR EACH VECTOR)
C
C   COMPUTE R*Q-SUB-V*R(TRANSPOSE), STORE IN A(NB1)
C   COMPUTE THE DIAGONAL ELEMENT FIRST
C
        DO 120 I=1,3
          II=3*(IVEC-1)
          IB=II+I
          DO 120 J=1,3
            JB=3*(IVEC-1)+J
            WORK1(I,J)=0.0
            WORK2(I,J)=0.0
            DO 120 L=1,3
              IF (II+L .LE. JB) THEN
                K = INX(II+L,NR) - (II+L)+ JB
              ELSE
                K = INX(JB,NR) + II+L - JB
              ENDIF
              TERM1= RM(I,L)*A(NB1+K)
              WORK1(I,J)=WORK1(I,J)+ TERM1
  120   WORK2(I,J)=WORK2(I,J)+RM(I,L)*A(NB2+K)
        DO 121 I=1,3
          DO 121 J=I,3
            K=INX(II+I,NR)-(II+I) + 3*(IVEC-1) +J
            A(NB1+K)=0.0
            A(NB2+K)=0.0
            DO 121 L=1,3
              A(NB1+K)=A(NB1+K)+WORK1(I,L)*RM(J,L)
  121   A(NB2+K)=A(NB2+K)+WORK2(I,L)*RM(J,L)
C
C   (JVEC>IVEC)  ALL ELEMENTS ARE STORED IN A
        IF (IVEC.LT.NVEC) THEN
C   COMPUTE R*Q FOR THE OTHER SUBMATRICES IN ROW IVEC
          DO 140 JVEC=IVEC+1,NVEC
            DO 130 I=1,3
              DO 130 J=1,3
                K=INX(II+I,NR)-(II+I) + 3*(JVEC-1) +J
                WORK1(I,J)=A(NB1+K)
  130       WORK2(I,J)=A(NB2+K)
            DO 135 I=1,3
              DO 135 J=1,3
                K=INX(II+I,NR)-(II+I) + 3*(JVEC-1) +J
                A(NB1+K)=0.0
                A(NB2+K)=0.0
                DO 135 L=1,3
                  A(NB1+K)=A(NB1+K)+RM(I,L)*WORK1(L,J)
  135       A(NB2+K)=A(NB2+K)+RM(I,L)*WORK2(L,J)
  140     CONTINUE
        ENDIF
C
C   COMPUTE Q*R(TRANSPOSE) FOR THE OTHER SUBMATRICES IN COLUMN IVEC
C
        IF (IVEC.GT.1) THEN
          DO 160 JVEC=1,IVEC-1
            DO 150 I=1,3
              DO 150 J=1,3
                IROW=3*(JVEC-1)+I
                JCOL=3*(IVEC-1)+J
                K=INX(IROW,NR)-IROW+JCOL
                WORK1(I,J)=A(NB1+K)
  150       WORK2(I,J)=A(NB2+K)
C
            DO 155 I=1,3
              DO 155 J=1,3
                IROW=3*(JVEC-1)+I
                JCOL=3*(IVEC-1)+J
                K=INX(IROW,NR)-IROW+JCOL
                A(NB1+K)=0.0
                A(NB2+K)=0.0
                DO 155 L=1,3
                  A(NB1+K)=A(NB1+K) + WORK1(I,L)*RM(J,L)
  155       A(NB2+K)=A(NB2+K) + WORK2(I,L)*RM(J,L)
  160     CONTINUE
        ENDIF

  250 CONTINUE

C
*** (STEP 9)  COMPUTE REDUNDANCY NUMBERS FOR NEU COMPONENTS
***          STORE IN GPSRN
      DO 280 I = 1, NR
        W = 0.D0
        DO 285 J = 1, NR
          IF (I .LE. J) THEN
            K = INX(I,NR) - I + J
          ELSE
            K = INX(J,NR) + I - J
          ENDIF
          SIGA = A(NB1+K)
          SIGL = A(NB2+K)
          W = W + SIGA*SIGL
  285   CONTINUE
        RN = 1.D0 - W
        GPSRNS(I) = RN
  280 CONTINUE
      RETURN
      END
