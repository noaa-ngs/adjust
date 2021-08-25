C %P%

      SUBROUTINE OBSEQW (IUO, FULL, G, B, NX, NVEC, NR, NC, LENG,
     &                   NICM, ICM, KINDS, ISNS, JSNS, LOBS, IAUX,
     &                   IVF, FATAL, IOBS, IGRT,CVFLAG,PROJID)

*** COMPUTE & DECORRELATE OBS.EQ. FROM WORK AREA

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( NVECS = 700, LENC = 10 )
      DIMENSION ICM((NVECS+1)*3+1+3)
      DIMENSION KINDS(NVECS*3), ISNS(NVECS*3), JSNS(NVECS*3)
      DIMENSION LOBS(NVECS*3)
      LOGICAL ADDCON, FULL, FATAL, NOBIGV, CVFLAG
      DIMENSION G(NR,NC), B(*), NX(*)
      DIMENSION SXYZ(3,3),SNEU(3,3)
      DIMENSION RM(3,3),W(3)
      DIMENSION IC(LENC), C(LENC)
      COMMON/UNITS/ LUNIT
**     2/7/2003  ** override input covariance matrix
** v 4.32vf
      LOGICAL IVCGPS,VSlogic
      COMMON/VCREC/SIGH,SIGU, VFHOR,VFUP, IVCGPS
      CHARACTER*13 VSPROJ
      PARAMETER (MAXPRJ=2500)
      COMMON/VSREC/SIGHS(MAXPRJ), SIGUS(MAXPRJ), ISETHU, VSPROJ(MAXPRJ),
     &             VSlogic
      CHARACTER*13 PROJID
      CHARACTER*99 SCCSID
C
****************** end 2/7/2003 insert
C
C      SCCSID='$Id: obseqw.for 65375 2012-09-20 15:10:51Z jarir.saleh $	20$Date: 2008/08/25 16:02:55 $ NGS'

*** ABORT/PRINT LARGE RESIDUALS

      NOBIGV = .FALSE.

*** DEFINE SIGNAL FOR A GPS OBS. EQ. SECTION

      ISIGNL = 1000


***v 4.27
*** STORE STANDARD ERROR OF OBSERVATIONS AND CORRELATION IN G-STRUCTURE

      CALL GLOCAT(N1,N2,N3,N4,N5,NVEC,LENG)
      DO 100 J = 1,NVEC
        I = 3 * (J-1) + 1
        G(I,N5) = G(I,I)
        G(I+1,N5) = G(I+1,I+1)
        G(I+2,N5) = G(I+2,I+2)
        G(I,N5+1) = G(I,I+1)
        G(I+1,N5+1) = G(I+1,I+2)
        G(I+2,N5+1) = G(I,I+2)
  100 CONTINUE	 

******
*******************3-27*****************
*       WRITE  (LUNIT,990) NVEC     
*990    FORMAT ( ' NVEC= ' ,I5,/)

*       DO 997 I = 1,NR

*       WRITE  (LUNIT,998) G(I,I),I
*998    FORMAT ( ' G1= ',1PD9.1,' ROW/COL= ',I5)
*997    CONTINUE
************************************************

*** TRANSFORM CORR-TYPE MATRIX TO COVARIANCES
      IF(.NOT.CVFLAG) THEN
        DO 3 I = 1, NR
          GII = G(I,I)
          DO 3 J = 1, NR
            IF (I .NE. J) G(I,J) = G(I,J)*GII*G(J,J)
    3   CONTINUE
      ENDIF
      DO 5 I = 1, NR
        G(I,I) = G(I,I)*G(I,I)
C        G(I,1) = G(I,I)
    5 CONTINUE
**
**     2/7/2003  ** override input covariance matrix
C
      IF (ISETHU.GT.0) THEN
        DO 85 I=1,ISETHU
          IF(PROJID.EQ.VSPROJ(I)) THEN
 	    CALL SCALEG (G, NR, NC, LENG,N1,N2,N3,N4,N5, NVEC, B,
     &                   KINDS, ISNS, JSNS, SIGHS(I),SIGUS(I))
            GO TO 86
          ENDIF
  85    CONTINUE
  86    CONTINUE
      ENDIF
C
C
      DO 95 I=1,NR
C MARK REJECTED OBSERVATIONS
      IF(G(I,N5).EQ.100.) THEN
        DO 93 J=1,NR
          G(J,I)=0.0
   93     G(I,J)=0.0
          G(I,I)=10000.

** Oct. 15, 2003
        ELSE
          G(I,N5)=SQRT(G(I,I))
** end insert

        ENDIF
   95 CONTINUE

******************************************
c 7/17/05 - 7/25/05 mike potterfield 
c try to update the correlation coefficient here
      DO 101 J = 1,NVEC
         I = 3 * (J-1) + 1
c         G(I,N5) = G(I,I)
c         G(I+1,N5) = G(I+1,I+1)
c         G(I+2,N5) = G(I+2,I+2)
c next is correlation XY         
         IF (G(I,I) .EQ. 10000.) THEN
           G(I,N5+1) = 0.
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
  101 CONTINUE   

*************************************
*** CHOLESKY FACTOR COVARIANCE MATRIX IN WORK SPACE

      DO 1 I = 1, NR
        I1 = I - 1
        IF (I1 .GT. 0) THEN
          DO 2 K = 1, I1
            G(I,I) = G(I,I) - G(K,I)*G(K,I)
    2     CONTINUE
        ENDIF

*** SINGULARITY TEST

****************************3-27
*       WRITE  (LUNIT,999) G(I,I),I
*999    FORMAT ( ' G= ',1PD9.1,' ROW/COL= ',I5)
**********************************

        IF (G(I,I) .LT. 1.D-10) THEN
          WRITE (LUNIT,7) I
    7     FORMAT ('0CORRELATION MATRIX SINGULARITY IN OBSEQW--',I5)
          CALL ABORT2
        ENDIF

*** RESUME CHOLESKY FACTOR

        G(I,I) = DSQRT( G(I,I) )
        IF (I+1 .LE. NR) THEN
          DO 6 J = I+1, NR
            IF (I1 .GT. 0) THEN
              DO 4 K = 1, I1
                G(I,J) = G(I,J) - G(K,I)*G(K,J)
    4         CONTINUE
            ENDIF
            G(I,J) = G(I,J)/G(I,I)
  6       CONTINUE
        ENDIF
    1 CONTINUE

*** UPDATE CONNECTIVITY AND WRITE GPS OBS.EQ. SECTION

      IF ( .NOT. ADDCON(ICM, NICM, NX) ) THEN
        WRITE (LUNIT,666) NICM, (ICM(I), I = 1, NICM)
  666   FORMAT ('0INSUFFICIENT MEMORY GPS VECTORS',20I3)
        CALL ABORT2
      ENDIF

      LENGX = 0
      CMOX = 0.D0
      OBX = 0.D0
      SDXX = 0.D0

      WRITE (IUO) ISIGNL, NVEC, IAUX, IC, C, LENGX, CMOX, OBX, SDXX,
     &            IOBS, IVF, IAUX, IGRT
      WRITE (IUO) ICM, NICM, KINDS, ISNS, JSNS, LOBS
      CALL GOBSEQ (G, NR, NC, NVEC, LENG, B, ICM, NICM, KINDS, ISNS,
     &             JSNS, LOBS, IAUX, IVF, FATAL, NOBIGV, N4, IGRT)
      DO 10 I = 1, NR
        WRITE (IUO) (G(I,J), J = 1, NC)
   10 CONTINUE

      FULL = .FALSE.
      CALL NEWICM (NICM)

      RETURN
      END

      SUBROUTINE SCALEG (G, NR, NC, LENG,N1,N2,N3,N4,N5, NVEC, B,
     &   KINDS,ISNS, JSNS, SFHOR,SFUP)
********************************************************************************
*** SCALE AN INPUT COVARIANCE MATRIX (FOR A GPS SESSION)
***    BY APPLYING SCALE FACTORS FOR HORIZONTAL AND VERTICAL COMPONENETS
***    THE MATRIX IS STORED IN G
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( NVECS = 700, LENC = 10 )
      DIMENSION KINDS(NVECS*3), ISNS(NVECS*3), JSNS(NVECS*3)
      LOGICAL ADDCON, FULL, FATAL, NOBIGV, CVFLAG
      DIMENSION G(NR,NC), B(*)
      DIMENSION RM(3,3),W(3)
*      COMMON/UNITS/ LUNIT
**     2/7/2003  ** override input covariance matrix
C
C
      DO 50 IVEC=1,NVEC
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
C  PROPAGATE COVARIANCE MATRIX FROM XYZ TO LOCAL HORIZON
C
C      COMPUTE R*C
      II=3*(IVEC-1)
      DO 24 JVEC=1,NVEC
        DO 23 JCOL=1,3
          JJ=3*(JVEC-1)+JCOL
          DO 21 I=1,3
            W(I)=0.0
            DO 20 J=1,3
   20         W(I)=W(I)+RM(J,I)*G(II+J,JJ)
   21     CONTINUE
          DO 22 I=1,3
   22       G(II+I,JJ)=W(I)
   23   CONTINUE
   24 CONTINUE
C
C     COMPUTE R*C*R(TRANSPOSE)
      DO 34 JVEC=1,NVEC
        DO 33 IROW=1,3
          JJ=3*(JVEC-1)+IROW
          DO 31 I=1,3
            W(I)=0.0
            DO 30 J=1,3
   30         W(I)=W(I)+G(JJ,II+J)*RM(J,I)
   31     CONTINUE
          DO 32 I=1,3
   32       G(JJ,II+I)=W(I)
   33   CONTINUE
   34 CONTINUE
C
   50 CONTINUE
C
C     NOW SCALE ALL NUMBERS
C
      DO 55 I=1,NR
        IR=MOD(I-1,3)+1
        DO 55 J=1,NR
          JC=MOD(J-1,3)+1
          IF (IR .EQ. 1 .OR. IR .EQ.2) THEN 
	    G(I,J)=G(I,J)*SFHOR
          ELSE
            G(I,J)=G(I,J)*SFUP
          ENDIF
          IF (JC .EQ. 1 .OR. JC .EQ.2) THEN 
            G(I,J)=G(I,J)*SFHOR
          ELSE
            G(I,J)=G(I,J)*SFUP
          ENDIF
   55 CONTINUE
C
C
C     PROPAGATE COVARIANCE MATRIX FROM LOCAL HORIZON SYSTEM BACK TO XYZ
C
      DO 90 IVEC=1,NVEC
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
C  PROPAGATE COVARIANCE MATRIX FROM LOCAL HORIZON TO XYZ
C
C  COMPUTE R*C
        II=3*(IVEC-1)
        DO 64 JVEC=1,NVEC
          DO 63 JCOL=1,3
          JJ=3*(JVEC-1)+JCOL
            DO 61 I=1,3
              W(I)=0.0
              DO 60 J=1,3
   60           W(I)=W(I)+RM(I,J)*G(II+J,JJ)
   61       CONTINUE
            DO 62 I=1,3
   62         G(II+I,JJ)=W(I)
   63     CONTINUE
   64   CONTINUE
C
C  COMPUTE R*C*R(TRANSPOSE)
        DO 74 JVEC=1,NVEC
          DO 73 IROW=1,3
            JJ=3*(JVEC-1)+IROW
            DO 71 I=1,3
              W(I)=0.0
              DO 70 J=1,3
   70           W(I)=W(I)+G(JJ,II+J)*RM(I,J)
   71       CONTINUE
            DO 72 I=1,3
   72         G(JJ,II+I)=W(I)
   73     CONTINUE
   74   CONTINUE
C
   90 CONTINUE
C
C
************   end insert for overridding input covariance matrix

      RETURN
      END
