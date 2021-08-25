C %P%

      SUBROUTINE CMRHS2 (G, N, N1)
********************************************************************************
* SINGLE SUBSTITUTION FOR DECORRELATION FROM A CHOL. FACTOR
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION G(N,*)

      CHARACTER*99 SCCSID

C      SCCSID='$Id: cmrhs2.for 44197 2010-07-13 12:26:03Z bruce.tran $	20$Date: 2007/11/20 15:22:13 $ NGS'

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


      SUBROUTINE CMRHS3 (G, N, H, NH, ICOL)
********************************************************************************
* COMPUTE RT INVERSE TIMES V AND STORE IT BACK IN V
* WHERE RT IS IN THE UPPER TRIANGLE OF G AND HAS N ROWS
* V IS IN COLUME ICOL OF H, WHICH HAS NH ROWS
* 
*** SINGLE SUBSTITUTION FOR DECORRELATION FROM A CHOL. FACTOR
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION G(N,*), H(NH,*)

      DO 1 I = 1, N
        I1 = I - 1
        TEMP = 0.D0
        IF (I1 .GT. 0) THEN
          DO 2 K = 1, I1
            TEMP = TEMP + G(K,I)*H(K,ICOL)
    2     CONTINUE
        ENDIF
        H(I,ICOL) = (H(I,ICOL)-TEMP)/G(I,I)
    1 CONTINUE

      RETURN
      END

