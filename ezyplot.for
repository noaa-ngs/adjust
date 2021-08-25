C %P%

      SUBROUTINE SETRNG (XSTRT, XSTOP, YSTRT, YSTOP)
********************************************************************************
* SET UNITS FOR PLOTTING
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ROW(48)
      CHARACTER*99 SCCSID
      COMMON /ROWS/ ROW
      COMMON /GRID/ X0, DX, Y0, DY, RDX, RDY, NROWS, NCOLS

C     SCCSID='$Id: ezyplot.for 44197 2010-07-13 12:26:03Z bruce.tran $	20$Date: 2007/11/20 15:22:44 $ NGS'
      NCOLS = 80
      NROWS = 48

      IF (XSTOP .LE. XSTRT) STOP 888
      IF (YSTOP .LE. YSTRT) STOP 889

      X0 = XSTRT
      DX = (NCOLS - 1)/(XSTOP - XSTRT)
      RDX = 1.D0/DX

      Y0 = YSTRT
      DY = (NROWS - 1)/(YSTOP - YSTRT)
      RDY = 1.D0/DY

      RETURN
      END


      SUBROUTINE CLRPLT
********************************************************************************
* CLEAR THE PLOT
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ROW(48), B80
      COMMON /ROWS/ ROW
      COMMON /GRID/ X0, DX, Y0, DY, RDX, RDY, NROWS, NCOLS

      DATA B80 /'                                                       
     &                         '/

      DO 1 IY = 1, NROWS
        ROW(IY) = B80
    1 CONTINUE

      RETURN
      END


      SUBROUTINE PUTXAX (Y, CHR)
********************************************************************************
* PLOT AN X AXIS AT Y ORDINATE
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*1 CHR
      CHARACTER*80 ROW(48)
      COMMON /ROWS/ ROW
      COMMON /GRID/ X0, DX, Y0, DY, RDX, RDY, NROWS, NCOLS

      IY = (Y - Y0)*DY + 1.4999D0

      IF (IY .LT. 1)     RETURN
      IF (IY .GT. NROWS) RETURN

      DO 1 IX = 1, NCOLS
        ROW(IY)(IX:IX) = CHR
    1 CONTINUE

      RETURN
      END


      SUBROUTINE PUTYAX (X, CHR)
********************************************************************************
* PLOT A Y AXIS AT X ORDINATE
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*1 CHR
      CHARACTER*80 ROW(48)
      COMMON /ROWS/ ROW
      COMMON /GRID/ X0, DX, Y0, DY, RDX, RDY, NROWS, NCOLS

      IX = (X - X0)*DX + 1.4999D0

      IF (IX .LT. 1)     RETURN
      IF (IX .GT. NCOLS) RETURN

      DO 1 IY = 1, NROWS
        ROW(IY)(IX:IX) = CHR
    1 CONTINUE

      RETURN
      END


      SUBROUTINE PUTXTK (YHOR, XST, XINC, CHR)
********************************************************************************
* PLOT TICKS ON X AXIS
* START AT XST, INCREMENT BY XINC, BOTH WAYS
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*1 CHR
      CHARACTER*80 ROW(48)
      COMMON /ROWS/ ROW
      COMMON /GRID/ X0, DX, Y0, DY, RDX, RDY, NROWS, NCOLS

      IF (XINC .LE. 0.D0) RETURN
      IF ( DINT(XINC*DX) .LT. 1 ) XINC = RDX
      IY = (YHOR - Y0)*DY + 1.4999D0
      IF (IY .LT. 1)     RETURN
      IF (IY .GT. NROWS) RETURN

*** LOOP UNTIL MAX RANGE EXCEEDED

      X = XST
  100 IX = (X - X0)*DX + 1.4999D0
      IF (IX .GE .1) THEN
        IF (IX .GT. NCOLS) GO TO 199
        ROW(IY)(IX:IX) = CHR
      ENDIF
      X = X + XINC
      GO TO 100

*** LOOP UNTIL MIN RANGE EXCEEDED

  199 X = XST
  200 IX = (X - X0)*DX + 1.4999D0
      IF (IX .LE. NCOLS) THEN
        IF (IX .LT. 1) GO TO 300
        ROW(IY)(IX:IX) = CHR
      ENDIF
      X = X - XINC
      GO TO 200

  300 RETURN
      END


      SUBROUTINE PUTYTK (XVERT, YST, YINC, CHR)
********************************************************************************
* PLOT TICKS ON Y AXIS
* START AT YST, INCREMENT BY YINC, BOTH WAYS
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*1 CHR
      CHARACTER*80 ROW(48)
      COMMON /ROWS/ ROW
      COMMON /GRID/ X0, DX, Y0, DY, RDX, RDY, NROWS, NCOLS

      IF (YINC .LE. 0.D0) RETURN
      IF ( DINT(YINC*DY) .LT. 1 ) YINC = RDY
      IX = (XVERT - X0)*DX + 1.4999D0
      IF (IX .LT. 1)     RETURN
      IF (IX .GT. NCOLS) RETURN

*** LOOP UNTIL MAX RANGE EXCEEDED

      Y = YST
  100 IY = (Y - Y0)*DY + 1.4999D0
      IF (IY .GE. 1) THEN
        IF (IY .GT. NROWS) GO TO 199
        ROW(IY)(IX:IX) = CHR
      ENDIF
      Y = Y + YINC
      GO TO 100

*** LOOP UNTIL MIN RANGE EXCEEDED

  199 Y = YST
  200 IY = (Y - Y0)*DY + 1.4999D0
      IF (IY .LE. NROWS) THEN
        IF (IY .LT. 1) GO TO 300
        ROW(IY)(IX:IX) = CHR
      ENDIF
      Y = Y - YINC
      GO TO 200

  300 RETURN
      END


      SUBROUTINE PUTXY (X, Y, CHR)
********************************************************************************
* PLOT CHARACTER AT X,Y
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*1 CHR
      CHARACTER*80 ROW(48)
      COMMON /ROWS/ ROW
      COMMON /GRID/ X0, DX, Y0, DY, RDX, RDY, NROWS, NCOLS

      IX = (X - X0)*DX + 1.4999D0
      IY = (Y - Y0)*DY + 1.4999D0

      IF (IX .LT. 1)     RETURN
      IF (IX .GT. NCOLS) RETURN
      IF (IY .LT. 1)     RETURN
      IF (IY .GT. NROWS) RETURN

      ROW(IY)(IX:IX) = CHR

      RETURN
      END


      SUBROUTINE GETPLT (IUNIT)
********************************************************************************
* DUMP THE PLOT
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*80 ROW(48)
      COMMON /ROWS/ ROW
      COMMON /GRID/ X0, DX, Y0, DY, RDX, RDY, NROWS, NCOLS
      COMMON/UNITS/ LUNIT

      IF (IUNIT .GT. 0) THEN
        DO 3 IY = NROWS, 1, -1
          WRITE (IUNIT, 4) ROW(IY)
    3   CONTINUE
    4   FORMAT (' ', A80)
      ELSE
        DO 5 IY = NROWS, 1, -1
          WRITE (LUNIT, 4) ROW(IY)
    5   CONTINUE
      ENDIF

      RETURN
      END


      SUBROUTINE PUTLIN (XSTART, YSTART, XSTOP, YSTOP, CHR)
********************************************************************************
* PLOT A LINE
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*1 CHR
      CHARACTER*80 ROW(48)
      COMMON /ROWS/ ROW
      COMMON /GRID/ X0, DX, Y0, DY, RDX, RDY, NROWS, NCOLS

      IF (YSTOP .EQ. YSTART) THEN
        IF (XSTOP .EQ. XSTART) THEN
          CALL PUTXY (XSTART, YSTART, CHR)
          RETURN
        ELSE
          DYDX = 0.D0
          B = YSTART
          IX1 = (XSTART - X0)*DX + 1.4999D0
          IX2 = (XSTOP  - X0)*DX + 1.4999D0
          CALL COLLOP (IX1, IX2, DYDX, B, CHR)
        ENDIF
      ELSE
        IF (XSTOP .EQ. XSTART) THEN
          DXDY = 0.D0
          A = XSTART
          IY1 = (YSTART - Y0)*DY + 1.4999D0
          IY2 = (YSTOP  - Y0)*DY + 1.4999D0
          CALL ROWLOP (IY1, IY2, DXDY, A, CHR)
        ELSE
          DXDY = (XSTOP - XSTART)/(YSTOP - YSTART)
          RDXRDY = RDX/RDY
          IF ( DABS(DXDY) .LT. RDXRDY ) THEN
            A = XSTART - DXDY*YSTART
            IY1 = (YSTART - Y0)*DY + 1.4999D0
            IY2 = (YSTOP  - Y0)*DY + 1.4999D0
            CALL ROWLOP (IY1, IY2, DXDY, A, CHR)
          ELSE
            DYDX = 1.D0/DXDY
            B = YSTART - DYDX*XSTART
            IX1 = (XSTART - X0)*DX + 1.4999D0
            IX2 = (XSTOP  - X0)*DX + 1.4999D0
            CALL COLLOP (IX1, IX2, DYDX, B, CHR)
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END


      SUBROUTINE ROWLOP (IY1, IY2, DXDY, A, CHR)
********************************************************************************
* DRAW A LINE--LOOP BY ROWS
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*1 CHR
      CHARACTER*80 ROW(48)
      COMMON /ROWS/ ROW
      COMMON /GRID/ X0, DX, Y0, DY, RDX, RDY, NROWS, NCOLS

*** ORDER AND TEST ROW LIMITS

      IF (IY2 .LT. IY1) THEN
        TEMP = IY1
        IY1 = IY2
        IY2 = TEMP
      ENDIF
      IF (IY1 .GT. NROWS) RETURN
      IF (IY2 .LT. 1)     RETURN
      IF (IY1 .LT. 1)     IY1 = 1
      IF (IY2 .GT. NROWS) IY2 = NROWS

*** PLOT LINE

      DO 1 IY = IY1, IY2
        Y = (IY - 1)*RDY + Y0
        X = DXDY*Y + A
        IX = (X - X0)*DX + 1.4999D0
        IF (IX .LT. 1)     GO TO 1
        IF (IX .GT. NCOLS) GO TO 1
        ROW(IY)(IX:IX) = CHR
    1 CONTINUE

      RETURN
      END


      SUBROUTINE COLLOP (IX1, IX2, DYDX, B, CHR)
********************************************************************************
* DRAW A LINE--LOOP BY COLUMNS
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*1 CHR
      CHARACTER*80 ROW(48)
      COMMON /ROWS/ ROW
      COMMON /GRID/ X0, DX, Y0, DY, RDX, RDY, NROWS, NCOLS

*** ORDER AND TEST COLUMN LIMITS

      IF (IX2 .LT. IX1) THEN
        TEMP = IX1
        IX1 = IX2
        IX2 = TEMP
      ENDIF
      IF (IX1 .GT. NCOLS) RETURN
      IF (IX2 .LT. 1)     RETURN
      IF (IX1 .LT. 1)     IX1 = 1
      IF (IX2 .GT. NCOLS) IX2 = NCOLS

*** PLOT LINE

      DO 1 IX = IX1, IX2
        X = (IX - 1)*RDX + X0
        Y = DYDX*X + B
        IY = (Y - Y0)*DY + 1.4999D0
        IF (IY .LT. 1)     GO TO 1
        IF (IY .GT. NROWS) GO TO 1
        ROW(IY)(IX:IX) = CHR
    1 CONTINUE

      RETURN
      END


      SUBROUTINE COMINC (ZMIN, ZMAX, ZINC)
********************************************************************************
* COMPUTE AN INCREMENT FOR TICK MARKS
********************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      DEL = ZMAX - ZMIN

*** AIM AT 10 TICKS PER AXIS, (MAY BE FEWER)

      TICK = DEL/10.D0

*** SCALE FROM 1.0 TO 9.999999

      VAL = DLOG10(TICK)
      IF (VAL .GE. 0.D0) THEN
        ILOG = VAL
      ELSE
        ILOG = VAL - 1.D0
      ENDIF
      STICK = TICK/10.D0**ILOG

*** SELECT A GRANULARITY

      IF (STICK .LE. 1.1D0) THEN
        ZINC = 10.D0**ILOG
      ELSEIF (STICK .LE. 2.D0) THEN
        ZINC = 2.D0*10.D0**ILOG
      ELSEIF (STICK .LE. 5.D0) THEN
        ZINC = 5.D0*10.D0**ILOG
      ELSE
        ZINC = 10.D0**(ILOG + 1)
      ENDIF

      RETURN
      END
