C %P%

      LOGICAL FUNCTION OPENN(A,NX,N,LENGTH)

*** INITIALIZE NORMAL EQUATION MATRIX

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*99 SCCSID

      DIMENSION A(*),NX(*)
      COMMON/DENNIS/NR,N1,N2,ISTATE

C      SCCSID='$Id: hogk.for 81724 2014-12-22 16:17:53Z jarir.saleh $	20$Date: 2007/11/20 15:23:51 $ NGS'

*** SET STATE TO UNREDUCED

      ISTATE=0
      NR=N
      N1=N+1
      N2=N1+N

*** ZERO THE MATRIX

      IF(NX(N1).LE.LENGTH) THEN
        OPENN=.TRUE.
        DO 1 I=1,NX(N1)
    1   A(I)=0.D0
      ELSE
        OPENN=.FALSE.
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION GETA(II,JJ,VAL,A,NX)

*** RETURN ELEMENT OF MATRIX

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),NX(*)
      COMMON/DENNIS/N,N1,N2,ISTATE

*** CONVERT EXTERNAL NUMBERS TO INTERNAL

      IF(II.LT.N1) THEN
        I=NX(N2+II)
      ELSE
        I=N1
      ENDIF
      IF(JJ.LT.N1) THEN
        J=NX(N2+JJ)
      ELSE
        J=N1
      ENDIF

      IF(I.LE.J) THEN
        IF(J.EQ.N1) THEN
          VAL=A(NX(N)+I)
          GETA=.TRUE.
        ELSEIF(I.EQ.N) THEN
          VAL=A(NX(N))
          GETA=.TRUE.
        ELSE
          INDEX=NX(I)-I+J
          IF(INDEX.GE.NX(I+1)) THEN
            VAL=0.D0
            GETA=.FALSE.
          ELSE
            VAL=A(INDEX)
            GETA=.TRUE.
          ENDIF
        ENDIF
      ELSE
        IF(J.EQ.N) THEN
          VAL=A(NX(N))
          GETA=.TRUE.
        ELSEIF(I.EQ.N1) THEN
          VAL=A(NX(N)+J)
          GETA=.TRUE.
        ELSE
          INDEX=NX(J)-J+I
          IF(INDEX.GE.NX(J+1)) THEN
            VAL=0.D0
            GETA=.FALSE.
          ELSE
            VAL=A(INDEX)
            GETA=.TRUE.
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION PUTA(II,JJ,VAL,A,NX)

*** STORE ELEMENT OF MATRIX

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),NX(*)
      COMMON/DENNIS/N,N1,N2,ISTATE

*** CONVERT EXTERNAL NUMBERS TO INTERNAL

      IF(II.LT.N1) THEN
        I=NX(N2+II)
      ELSE
        I=N1
      ENDIF
      IF(JJ.LT.N1) THEN
        J=NX(N2+JJ)
      ELSE
        J=N1
      ENDIF

      IF(I.LE.J) THEN
        IF(J.EQ.N1) THEN
          A(NX(N)+I)=VAL
          PUTA=.TRUE.
        ELSEIF(I.EQ.N) THEN
          A(NX(N))=VAL
          PUTA=.TRUE.
        ELSE
          INDEX=NX(I)-I+J
          IF(INDEX.GE.NX(I+1)) THEN
            PUTA=.FALSE.
          ELSE
            A(INDEX)=VAL
            PUTA=.TRUE.
          ENDIF
        ENDIF
      ELSE
        IF(J.EQ.N) THEN
          A(NX(N))=VAL
          PUTA=.TRUE.
        ELSEIF(I.EQ.N1) THEN
          A(NX(N)+J)=VAL
          PUTA=.TRUE.
        ELSE
          INDEX=NX(J)-J+I
          IF(INDEX.GE.NX(J+1)) THEN
            PUTA=.FALSE.
          ELSE
            A(INDEX)=VAL
            PUTA=.TRUE.
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION ADDOBS(C,IC,LEN,OBS,WT,A,NX)

*** ACCUMULATE OBSERVATIONS INTO NORMALS

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL IN
      DIMENSION C(LEN),IC(LEN)
      DIMENSION A(*),NX(*)
      COMMON/DENNIS/N,N1,N2,ISTATE

*** LOOP OVER COEFFICIENTS

      A(NX(N1))=A(NX(N1))+OBS*WT*OBS
      DO 2 I=1,LEN
      II=NX(N2+IC(I))
      CIWT=C(I)*WT
      A(NX(N)+II)=A(NX(N)+II)+CIWT*OBS
      DO 2 J=I,LEN
      JJ=NX(N2+IC(J))

*** CHECK IF COEFFICIENTS WITHIN PROFILE

      IF(.NOT.IN(II,JJ,NX,INDEX)) THEN
        ADDOBS=.FALSE.
        RETURN
      ENDIF

*** ACCUMULATE

  2   A(INDEX)=A(INDEX)+CIWT*C(J)

*** ACCUMULATION SUCCESSFUL

      ADDOBS=.TRUE.

      RETURN
      END
      LOGICAL FUNCTION ADDCOR(C,IC,OBS,SIGOBS,LUNK,LOBS,A,NX,IFLAG)

*** ACCUMULATE CORRELATED OBSERVATIONS INTO NORMALS

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL IN
      DIMENSION C(*),IC(LUNK),SIGOBS(LOBS,LOBS),OBS(LOBS)
      DIMENSION A(*),NX(*)
      COMMON/DENNIS/N,N1,N2,ISTATE
      DATA TOL/1.D-10/

*** CHOLESKY FACTOR WEIGHT MATRIX

      DO 1 I=1,LOBS
      I1=I-1
      IF(I1.GT.0) THEN
        DO 2 K=1,I1
    2   SIGOBS(I,I)=SIGOBS(I,I)-SIGOBS(K,I)*SIGOBS(K,I)
      ENDIF

*** SINGULARITY TEST

      IF(SIGOBS(I,I).LT.TOL) THEN
        ADDCOR=.FALSE.
        IFLAG=1
        RETURN
      ENDIF

*** RESUME CHOLESKY FACTOR

      SIGOBS(I,I)=DSQRT(SIGOBS(I,I))
      IF(I+1.LE.LOBS) THEN
        DO 6 J=I+1,LOBS
        IF(I1.GT.0) THEN
          DO 4 K=1,I1
    4     SIGOBS(I,J)=SIGOBS(I,J)-SIGOBS(K,I)*SIGOBS(K,J)
        ENDIF
        SIGOBS(I,J)=SIGOBS(I,J)/SIGOBS(I,I)
    6   SIGOBS(J,I)=SIGOBS(I,J)
      ENDIF
    1 CONTINUE

*** TRANSFORM OBSERVATIONS

      CALL CHSOLV(SIGOBS,OBS,LOBS)

*** TRANSFORM COEFFICIENTS

      DO 3 I=1,LUNK
      J=(I-1)*LOBS+1
    3 CALL CHSOLV(SIGOBS,C(J),LOBS)

*** ACCUMULATE TRANSFORMED OBSERVATIONS -- UNIT WEIGHT

      DO 100 K=1,LOBS

*** LOOP OVER COEFFICIENTS

        A(NX(N1))=A(NX(N1))+OBS(K)*OBS(K)
        DO 5 I=1,LUNK
        II=NX(N2+IC(I))
        CI=C((I-1)*LOBS+K)
        A(NX(N)+II)=A(NX(N)+II)+CI*OBS(K)
        DO 5 J=I,LUNK
        JJ=NX(N2+IC(J))

*** CHECK IF COEFFICIENTS WITHIN PROFILE

        IF(.NOT.IN(II,JJ,NX,INDEX)) THEN
          ADDCOR=.FALSE.
          IFLAG=2
          RETURN
        ENDIF

*** ACCUMULATE

    5   A(INDEX)=A(INDEX)+CI*C((J-1)*LOBS+K)
  100 CONTINUE

*** ACCUMULATION SUCCESSFUL

      ADDCOR=.TRUE.
      IFLAG=0

      RETURN
      END
      SUBROUTINE CHSOLV(A,X,LENG)

*** SOLVE AY=X FOR Y AND STORE RESULT IN X (GIVEN CHOLESKY FACTOR, A)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(LENG,LENG),X(LENG)

*** SUBSTITUTION

      DO 5 I=1,LENG
      I1=I-1
      TEMP=0.D0
      IF(I1.GT.0) THEN
        DO 6 K=1,I1
    6   TEMP=TEMP+A(K,I)*X(K)
      ENDIF
    5 X(I)=(X(I)-TEMP)/A(I,I)

      RETURN
      END
      LOGICAL FUNCTION SOLVE(A,NX,TOL,ISING,GSING,LSING)

*** FACTOR AND SOLVE A LINEAR SYSTEM

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL IN
      DIMENSION ISING(*),GSING(*)
      DIMENSION A(*),NX(*)
      COMMON/DENNIS/N,N1,N2,ISTATE
      DATA BIG,TINY/1.D30,1.D-30/

      NM1=N-1
      LSING=0
      IF(ISTATE.EQ.0) THEN
        ISTATE=1
        SOLVE=.TRUE.

*** CHOLESKY FACTOR

        DO 1 I=N,1,-1
          IDIAG=NX(I)
          DIAG=A(IDIAG)
          IF(I.NE.N) THEN
            ITOP=NX(I+1)-1
            ILEN=ITOP-IDIAG
            IF(ILEN.GT.0) THEN
              ICOL=ILEN+I
              DO 2 J=I-1,1,-1
                JDIAG=NX(J)
                JTOP=NX(J+1)-1
                JLEN=JTOP-JDIAG
                JCOL=JLEN+J
                IF(JCOL.GT.I) THEN
                  JI=JDIAG-J+I
                  IF(ICOL.LT.JCOL) THEN
                    NP=ILEN
                    IP=ITOP-NP+1
                    JP=JTOP-(JCOL-ICOL)-NP+1
                  ELSE
                    NP=JLEN+J-I
                    JP=JTOP-NP+1
                    IP=ITOP-(ICOL-JCOL)-NP+1
                  ENDIF
                  CALL INNERP(A(IP),A(JP),NP,PROD)
                  A(JI)=A(JI)-PROD
                ENDIF
    2         CONTINUE
            ENDIF
            IA=IDIAG+1
            IB=NX(I+1)-1
            DO 4 JI=IA,IB
    4       A(IDIAG)=A(IDIAG)-A(JI)*A(JI)
          ENDIF

*** SINGULARITY TEST

          IF(DIAG.LE.0.D0) THEN
            GOOGE=0.D0
          ELSE
            GOOGE=A(IDIAG)/DIAG
          ENDIF
          IF(GOOGE.LT.TOL) THEN
            LSING=LSING+1
            ISING(LSING)=NX(N1+I)
            IF(GOOGE.LE.0.D0) THEN
              GSING(LSING)=0.D0
            ELSE
              GSING(LSING)=GOOGE
            ENDIF
            SOLVE=.FALSE.
            A(IDIAG)=BIG
          ENDIF

*** RESUME CHOLESKY FACTOR

          DIAG=DSQRT(A(IDIAG))
          A(IDIAG)=DIAG
          RDIAG=1.D0/DIAG
          DO 5 J=1,I-1
            JI=NX(J)-J+I
            IF(JI.LT.NX(J+1)) THEN
              TEMP=A(JI)*RDIAG
              IF(DABS(TEMP).LT.TINY) THEN
                A(JI)=0.D0
              ELSE
                A(JI)=TEMP
              ENDIF
            ENDIF
    5     CONTINUE
    1   CONTINUE

*** COMPLETE CHOLESKY FACTOR

*** FORWARD SUBSTITUTION

        NM1=N-1
        A(NX(N1)+N)=A(NX(N)+N)/A(NX(N))
        DO 7 I=NM1,1,-1
          TEMP=0.D0
          IA=NX(I)+1
          IB=NX(I+1)-1
          K=I
          DO 8 IK=IA,IB
            K=K+1
            TEMP=TEMP+A(IK)*A(NX(N1)+K)
    8     CONTINUE
    7   A(NX(N1)+I)=(A(NX(N)+I)-TEMP)/A(NX(I))

*** BASEMENT WINDOW

        TEMP=A(NX(N1))
        DO 10 I=1,N
          Y=A(NX(N1)+I)
   10   TEMP=TEMP-Y*Y
        A(NX(N1))=TEMP

*** BACKWARD SOLUTION

        DO 9 I=1,N
        I1=I-1
        TEMP=0.D0
        DO 6 K=1,I1
          KI=NX(K)-K+I
          IF(KI.LT.NX(K+1)) TEMP=TEMP+A(KI)*A(NX(N)+K)
    6   CONTINUE
    9   A(NX(N)+I)=(A(NX(N1)+I)-TEMP)/A(NX(I))
      ELSE
        SOLVE=.FALSE.
      ENDIF

      RETURN
      END
      SUBROUTINE INNERP(X,Y,N,PROD)

*** COMPUTE INNER PRODUCT

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(*),Y(*)

      PROD=0.D0
      DO 1 I=1,N
    1 PROD=PROD+X(I)*Y(I)

      RETURN
      END
      LOGICAL FUNCTION RESOLV(A,NX)

*** SOLVE THE FACTORED NORMAL EQUATIONS FOR A NEW R.H. SIDE

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL IN
      DIMENSION A(*),NX(*)
      COMMON/DENNIS/N,N1,N2,ISTATE

      IF(ISTATE.NE.1) THEN
        RESOLV=.FALSE.
      ELSE

*** FORWARD SUBSTITUTION

        RESOLV=.TRUE.

        NM1=N-1
        A(NX(N1)+N)=A(NX(N)+N)/A(NX(N))
        DO 7 I=NM1,1,-1
          TEMP=0.D0
          IA=NX(I)+1
          IB=NX(I+1)-1
          K=I
          DO 8 IK=IA,IB
            K=K+1
            TEMP=TEMP+A(IK)*A(NX(N1)+K)
    8     CONTINUE
    7   A(NX(N1)+I)=(A(NX(N)+I)-TEMP)/A(NX(I))

*** BASEMENT WINDOW

        TEMP=A(NX(N1))
        DO 10 I=1,N
          Y=A(NX(N1)+I)
   10   TEMP=TEMP-Y*Y
        A(NX(N1))=TEMP

*** BACKWARD SOLUTION

        DO 9 I=1,N
        I1=I-1
        TEMP=0.D0
        DO 6 K=1,I1
          KI=NX(K)-K+I
          IF(KI.LT.NX(K+1)) TEMP=TEMP+A(KI)*A(NX(N)+K)
    6   CONTINUE
    9   A(NX(N)+I)=(A(NX(N1)+I)-TEMP)/A(NX(I))
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION INVERT(A,NX)

*** INVERT THE CHOLESKY FACTOR

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL IN
      DIMENSION A(*),NX(*)
      COMMON/DENNIS/N,N1,N2,ISTATE

      NM1=N-1
      IF(ISTATE.NE.1) THEN
        INVERT=.FALSE.
      ELSE
        ISTATE=2
        INVERT=.TRUE.

*** PUT CHOLESKY RECIPROCALS ON DIAGONAL

        DO 1 I=1,N
    1   A(NX(I))=1.D0/A(NX(I))

*** COMPUTE GAUSSIAN FACTOR FROM CHOLESKY FACTOR

        DO 2 I=N,2,-1
          DIAG=A(NX(I))
          DO 7 J=1,I
            JI=NX(J)-J+I
            IF(JI.LT.NX(J+1)) A(JI)=DIAG*A(JI)
    7     CONTINUE
    2   CONTINUE
        A(1)=A(1)*A(1)

*** PROCEED BY COLUMNS--TOP TO BOTTOM

        DO 6 I=2,N

*** SAVE I-TH COL IN SCRATCH SPACE

          DO 3 J=1,I-1
            JI=NX(J)-J+I
            IF(JI.LT.NX(J+1)) A(NX(N1)+J)=-A(JI)
    3     CONTINUE

*** GET INVERSE OF OFF DIAGONAL ELEMENTS OF I-TH COL

          DO 5 J=1,I-1
            JI=NX(J)-J+I
            IF(JI.LT.NX(J+1)) THEN
              A(JI)=0.D0
              DO 4 K=1,I-1
                KI=NX(K)-K+I
                IF(KI.LT.NX(K+1)) THEN
                  YK=A(NX(N1)+K)
                  IF(K.LT.J) THEN
                    KJ=NX(K)-K+J
                    IF(KJ.LT.NX(K+1)) A(JI)=A(JI)+YK*A(KJ)
                  ELSE
                    JK=NX(J)-J+K
                    IF(JK.LT.NX(J+1)) A(JI)=A(JI)+YK*A(JK)
                  ENDIF
                ENDIF
    4         CONTINUE
            ENDIF
    5     CONTINUE

*** GET INVERSE OF DIAGONAL AT I-TH COL

          DO 8 J=1,I-1
            JI=NX(J)-J+I
            IF(JI.LT.NX(J+1)) THEN
              A(NX(I))=A(NX(I))+A(JI)*A(NX(N1)+J)
            ENDIF
    8     CONTINUE
    6   CONTINUE
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION PROP(C,IC,LEN,VAR,A,NX,IFLAG)

*** COMPUTE VARIANCE BY LINEAR ERROR PROPAGATION

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL IN
      DIMENSION C(LEN),IC(LEN)
      DIMENSION A(*),NX(*)
      COMMON/DENNIS/N,N1,N2,ISTATE

*** VERIFY PARAMETER COVARIANCE AVAILABLE

      IF(ISTATE.NE.2) THEN
        IFLAG=1
        PROP=.FALSE.
      ELSE
        VAR=0.D0

        DO 20 I=1,LEN
        II=NX(N2+IC(I))
        VAR=VAR+C(I)*C(I)*A(NX(II))
        I1=I+1
        IF(I1.LE.LEN) THEN
          C2=C(I)+C(I)
          DO 10 J=I1,LEN
          JJ=NX(N2+IC(J))
          IF(.NOT.IN(II,JJ,NX,INDEX)) THEN
            IFLAG=2
            PROP=.FALSE.
            RETURN
          ENDIF
   10     VAR=VAR+C2*C(J)*A(INDEX)
        ENDIF
   20   CONTINUE

        IFLAG=0
        PROP=.TRUE.
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION PROPCV(C1,IC1,LEN1,C2,IC2,LEN2,COV,A,NX,IFLAG)

*** COMPUTE COVARIANCE BY LINEAR ERROR PROPAGATION

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL GETA
      DIMENSION C1(LEN1),IC1(LEN1)
      DIMENSION C2(LEN2),IC2(LEN2)
      DIMENSION A(*),NX(*)
      COMMON/DENNIS/N,N1,N2,ISTATE

*** VERIFY PARAMETER COVARIANCES AVAILABLE

      IF(ISTATE.NE.2) THEN
        IFLAG=1
        PROPCV=.FALSE.
      ELSE
        COV=0.D0
        DO 20 I=1,LEN1
        CI=C1(I)
        II=IC1(I)
        DO 10 J=1,LEN2
          JJ=IC2(J)
          IF(.NOT.GETA(II,JJ,SIG,A,NX)) THEN
            IFLAG=2
            PROPCV=.FALSE.
            RETURN
          ENDIF
   10     COV=COV+CI*C2(J)*SIG
   20   CONTINUE

        IFLAG=0
        PROPCV=.TRUE.
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION CORR(A,NX)

*** COMPUTE CORRELATION MATRIX FROM COVARIANCE

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL IN
      DIMENSION A(*),NX(*)
      COMMON/DENNIS/N,N1,N2,ISTATE

      NM1=N-1
      IF(ISTATE.NE.2) THEN
        CORR=.FALSE.
      ELSE
        ISTATE=3
        CORR=.TRUE.

*** CONVERT VARIANCES TO STANDARD DEVIATIONS

        DO 1 I=1,N
    1   A(NX(I))=DSQRT(A(NX(I)))

*** DIVIDE EACH OFF-DIAGONAL ROW ELEMENT BY DIAGONAL

        DO 5 I=1,NM1
        RDIAG=1.D0/A(NX(I))
        IA=NX(I)+1
        IB=NX(I+1)-1
        IF(IA.LE.IB) THEN
          DO 2 IJ=IA,IB
    2     A(IJ)=A(IJ)*RDIAG
        ENDIF
    5   CONTINUE

*** DIVIDE EACH OFF-DIAGONAL COLUMN ELEMENT BY DIAGONAL

        DO 3 J=N,2,-1
        RDIAG=1.D0/A(NX(J))
        DO 3 I=1,J-1
        IF(IN(I,J,NX,IJ)) A(IJ)=A(IJ)*RDIAG
    3   CONTINUE

*** PLACE ONE'S ALONG DIAGONAL

        DO 4 I=1,N
    4   A(NX(I))=1.D0
      ENDIF

      RETURN
      END
      LOGICAL FUNCTION IN(I,J,NX,INDEX)

*** RETURN ELEMENT INDEX OF MATRIX PROFILE
*** CAUTION: I AND J ARE INTERNAL NUMBERS !!

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION NX(*)
      COMMON/DENNIS/N,N1,N2,ISTATE

      IF(I.LE.J) THEN
        IF(J.EQ.N1) THEN
          INDEX=NX(N)+I
          IN=.TRUE.
        ELSEIF(I.EQ.N) THEN
          INDEX=NX(N)
          IN=.TRUE.
        ELSE
          INDEX=NX(I)-I+J
          IF(INDEX.GE.NX(I+1)) THEN
            IN=.FALSE.
          ELSE
            IN=.TRUE.
          ENDIF
        ENDIF
      ELSE
        IF(J.EQ.N) THEN
          INDEX=NX(N)
          IN=.TRUE.
        ELSEIF(I.EQ.N1) THEN
          INDEX=NX(N)+J
          IN=.TRUE.
        ELSE
          INDEX=NX(J)-J+I
          IF(INDEX.GE.NX(J+1)) THEN
            IN=.FALSE.
          ELSE
            IN=.TRUE.
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END
      subroutine openg(n,nx,nxleng)

*** initialize linked neighbor list

      dimension nx(*)
      common/kathy/nr,n1,n2,iemx,ifree

*** compute indicies

      nr=n
      n1=n+1
      n2=n1+n

*** compute maximum number of edges
*** insure iemx (2*edges) is even

      iemx=((nxleng-n-n-1)/6)*2

      do 1 i=1,n
    1 nx(i)=-i
      ifree=1

      return
      end
      logical function addcon(ic,len,nx)

*** add observation equation connections to neighbor list

      logical addg
      dimension ic(*),nx(*)
      common/kathy/n,n1,n2,iemx,ifree

      if(len.le.1) then
c       write(*,2) len
    2   format('len error in addcon',i5)
        stop
      endif
      l1=len-1

*** loop over all connections

      do 1 i=1,l1
      do 1 j=i+1,len
      if(.not.addg(ic(i),ic(j),nx(1),nx(n1),nx(n1+iemx))) then
        addcon=.false.
        return
      endif
    1 continue
      addcon=.true.

      return
      end
      logical function addg(i,j,ihead,nbr,link)

*** add a connection

      dimension ihead(*),nbr(*),link(*)
      logical conect
      common/kathy/n,n1,n2,iemx,ifree

      addg=.true.

*** add connection if it does not exist

      if(.not.conect(i,j,ipoint,ihead,nbr,link)) then
        if(ipoint.eq.0) then
          ihead(i)=ifree
        else
          if(link(ipoint).ne.-i) then
c           write(*,1) i,j,ipoint,link(ipoint)
    1       format(4i5,' -- graph subroutine link error')
            stop
          else
            link(ipoint)=ifree
          endif
        endif
        nbr(ifree)=j
        link(ifree)=-i
        ifree=ifree+1
        if(ifree.gt.iemx) then
          addg=.false.
          return
        endif
      endif

*** add connection in other direction

      if(.not.conect(j,i,ipoint,ihead,nbr,link)) then
        if(ipoint.eq.0) then
          ihead(j)=ifree
        else
          if(link(ipoint).ne.-j) then
c           write(*,1) i,j,ipoint,link(ipoint)
            stop
          else
            link(ipoint)=ifree
          endif
        endif
        nbr(ifree)=i
        link(ifree)=-j
        ifree=ifree+1
        if(ifree.gt.iemx) then
          addg=.false.
          return
        endif
      endif

      return
      end
      logical function conect(i,j,ipntx,ihead,nbr,link)

*** is connection in the graph?
*** if not, find pointer to end of linked list

      dimension ihead(*),nbr(*),link(*)

      ipntx=0
      ipoint=ihead(i)

    1 if(ipoint.lt.0) then
        conect=.false.
        return
      else
        if(nbr(ipoint).eq.j) then
          conect=.true.
          return
        endif
        ipntx=ipoint
        ipoint=link(ipoint)
        go to 1
      endif

      return
      end
      subroutine reordr(nx)

*** reorder unknowns into king order
*** for minimum profile, row storage, upper triangle

      dimension nx(*)
      common/kathy/n,n1,n2,iemx,ifree

*** produce adjacency structure from linked list

      call trnslt(nx(1),nx(n1),nx(n1+iemx),nx(n1+3*iemx),nx(n1+2*iemx))

*** produce king order

      call gking(n,nx(n1+3*iemx),nx(n1+2*iemx),nx(n2+1),nx(n1+1),
     *            nx(3*n+2))

*** reverse order (for row profile)

      do 776 i=1,n
  776 nx(n2+1-i)=nx(n2+i)

*** invert order

      do 777 i=1,n
      j=nx(n1+i)
  777 nx(n2+j)=i

*** produce profile

      call profil(nx(1),nx(n1+3*iemx),nx(n1+2*iemx))

      return
      end
      subroutine trnslt(ihead,nbr,link,nbrpnt,nbrlst)

*** translate linked neighbor list into adjacency structure

      dimension ihead(*),nbr(*),link(*),nbrpnt(*),nbrlst(*)
      common/kathy/n,n1,n2,iemx,ifree

      ifree=1
      do 100 i=1,n
      ipoint=ihead(i)
      nbrpnt(i)=ifree
    1 if(ipoint.gt.0) then
        nbrlst(ifree)=nbr(ipoint)
        ipoint=link(ipoint)
        ifree=ifree+1
        go to 1
      endif
  100 continue
      nbrpnt(n1)=ifree

      return
      end
      subroutine profil(nx,nbrpnt,nbrlst)

*** find variable column profile from adjacency structure

      dimension nx(*),nbrpnt(*),nbrlst(*)
      common/kathy/n,n1,n2,iemx,ifree

*** store row widths in nx(2 thru n1)

      do 20 ir=1,n
      i=nx(n2+ir)
      lbegin=nbrpnt(ir)
      lend=nbrpnt(ir+1)-1
      lrow=0
      if(lbegin.le.lend) then
        do 10 ipt=lbegin,lend
        jc=nbrlst(ipt)
        j=nx(n2+jc)
        if(j.gt.i) then
          len=j-i
          if(len.gt.lrow) lrow=len
        endif
   10   continue
      endif
   20 nx(i+1)=lrow+1

*** convert row widths to profile

      nx(1)=1
      do 30 i=2,n
   30 nx(i)=nx(i)+nx(i-1)
      nx(n1)=nx(n)+n1

      return
      end
      subroutine gking(neqns,nbrpnt,nbrlst,neword,mask,levptr)

*** find king ordering for a general graph
*** for each component, ordering is found by subroutine king

      dimension nbrlst(*),mask(*),neword(*),levptr(*),nbrpnt(*)

*** initialize component mask

      do 100 i=1,neqns
  100 mask(i)=1
      num=1

*** loop over all equations

      do 200 i=1,neqns
      if(mask(i).eq.0) go to 200
      iroot=i

*** find pseudo-peripheral node

      call fnroot(iroot,nbrpnt,nbrlst,mask,nlvl,levptr,neword(num))

*** number component in king order

      call king(iroot,nbrpnt,nbrlst,mask,neword(num),lcomp,levptr)
      num=num+lcomp
      if(num.gt.neqns) return
  200 continue

      return
      end
      subroutine fnroot(iroot,nbrpnt,nbrlst,mask,nlvl,levptr,ls)

*** find pseudo-peripheral node

      dimension nbrlst(*),ls(*),mask(*),levptr(*),nbrpnt(*)

*** determine level structure rooted at iroot

      call rootls(iroot,nbrpnt,nbrlst,mask,nlvl,levptr,ls)
      lcomp=levptr(nlvl+1)-1
      if(nlvl.eq.1.or.nlvl.eq.lcomp) return

*** pick a node with minimum degree from the last level

  100 jstrt=levptr(nlvl)
      mindeg=lcomp
      iroot=ls(jstrt)
      if(lcomp.eq.jstrt) go to 400
        do 300 j=jstrt,lcomp
        node=ls(j)
        ndeg=0
        kstrt=nbrpnt(node)
        kstop=nbrpnt(node+1)-1
        do 200 k=kstrt,kstop
        nabor=nbrlst(k)
        if(mask(nabor).gt.0) ndeg=ndeg+1
  200   continue
        if(ndeg.ge.mindeg) go to 300
        iroot=node
        mindeg=ndeg
  300   continue

*** generate its rooted level structure

  400 call rootls(iroot,nbrpnt,nbrlst,mask,nunlvl,levptr,ls)
      if(nunlvl.le.nlvl) return
      nlvl=nunlvl
      if(nlvl.lt.lcomp) go to 100

      return
      end
      subroutine rootls(iroot,nbrpnt,nbrlst,mask,nlvl,levptr,ls)

*** generate level structure rooted at iroot

      dimension nbrlst(*),ls(*),mask(*),levptr(*),nbrpnt(*)

*** initialization

      mask(iroot)=0
      ls(1)=iroot
      nlvl=0
      lvlend=0
      lcomp=1

*** lbegin is pointer to beginning of current level
*** lvlend is pointer to end of current level

  200 lbegin=lvlend+1
      lvlend=lcomp
      nlvl=nlvl+1
      levptr(nlvl)=lbegin

*** generate next level by finding all masked neighbors of nodes in
*** the current level

      do 400 i=lbegin,lvlend
      node=ls(i)
      jstrt=nbrpnt(node)
      jstop=nbrpnt(node+1)-1
      if(jstop.lt.jstrt) go to 400
      do 300 j=jstrt,jstop
      nbr=nbrlst(j)
      if(mask(nbr).eq.0) go to 300
      lcomp=lcomp+1
      ls(lcomp)=nbr
      mask(nbr)=0
  300 continue
  400 continue

*** compute current level width--if nonzero, generate next level

      lvsize=lcomp-lvlend
      if(lvsize.gt.0) go to 200

*** reset mask to 1 for nodes in the level structure

      levptr(nlvl+1)=lvlend+1
      do 500 i=1,lcomp
      node=ls(i)
  500 mask(node)=1

      return
      end
      subroutine king(iroot,xadj,adj,mask,neword,nord,nlist)

*** number connected components using king algorithm

      implicit integer(a-z)
      dimension xadj(*),adj(*),mask(*),neword(*),nlist(*)
      common/kathy/nr,n1,n2,iemx,ifree

*** start with pseudo-periperial node, iroot
*** mark its unordered neighbors with -1, and maintain a neighbor count

      nord=1
      neword(nord)=iroot
      mask(iroot)=0
      ncount=0
      istrt=xadj(iroot)
      istop=xadj(iroot+1)-1
      do 2 j=istrt,istop
        k=adj(j)
        if(mask(k).eq.+1) then
          mask(k)=-1
          ncount=ncount+1
          nlist(ncount)=k
        endif
    2 continue
      if(ncount.le.0) return

*** main loop, execute until no more neighbors
*** locate the node (mnode) with the minumum adjacency count (minadj)

  100 minadj=nr+1
      mnode=0
      do 1 i=1,ncount
          inode=nlist(i)

*** pretend the unordered neighbor is ordered and get a trial count (iadj)

          iadj=ncount-1
          istrt=xadj(inode)
          istop=xadj(inode+1)-1
          do 4 j=istrt,istop
            if(mask(adj(j)).eq.+1) iadj=iadj+1
    4     continue

*** accumulate the minimum count

          if(iadj.lt.minadj) then
            minadj=iadj
            mnode=inode
          endif
    1 continue
      if(mnode.eq.0) stop 33333

*** order the winning neighbor & update the unordered neighbors
*** condense the neighbor list and add new neighbors

      nord=nord+1
      neword(nord)=mnode
      mask(mnode)=0
      j=0
      do 5 i=1,ncount
        if(nlist(i).ne.mnode) then
          j=j+1
          nlist(j)=nlist(i)
        endif
    5 continue
      ncount=ncount-1
      if(ncount.ne.j) stop 33433
      istrt=xadj(mnode)
      istop=xadj(mnode+1)-1
      do 3 j=istrt,istop
        k=adj(j)
        if(mask(k).eq.+1) then
          mask(k)=-1
          ncount=ncount+1
          nlist(ncount)=k
        endif
    3 continue

*** loop until no more unmarked neighors

      if(ncount.gt.0) go to 100

      return
      end
