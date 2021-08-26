* /ngslib/source/Adjust/SCCS/s.adjust.for
* file adjust.for
      PROGRAM ADJUST
***********************************************************************
*                                                                     *
* PROGRAM :   ADJUST                                                  *
*                                                                     *
* PURPOSE:    WEIGHTED LEAST SQUARES ADJUSTMENT OF HORIZONTAL FIELD   *
*             OBSERVATIONS, EARTH CENTERED DIFFERENTIAL XYZ DATA      *
*             (I.E. GLOBAL POSITIONING SYSTEM (GPS) DATA)             *
*             AND XYZ COORDINATE OBSERVATIONAL DATA (DOPPLER)         *
*             IN THE NATIONAL GEODETIC SURVEY DATA BASE               *
*             INPUT FORMATS (BLUE-BOOK, REVISED, JAN 1989 VERSION).   *
*                                                                     *
* VERSION CODE:  6.4.2 
*                                                                     *
* VERSION DATE:  (ccyy/mm/dd)  2017/11/06
*                                                                     *
*       AUTHORS:   ALICE R. DREW                                      *
*                  N/CG12X2                                           *
*                  NATIONAL GEODETIC SURVEY, NOS, NOAA                *
*                  ROCKVILLE, MD   20852                              *
*                                                                     *
*                  DENNIS G. MILBERT                                  *
*                  N/CG113                                            *
*                  NATIONAL GEODETIC SURVEY, NOS, NOAA                *
*                  ROCKVILLE, MD    20852                             *
*                                                                     *
*                  WILLIAM G. KASS                                    *
*                  N/CG121                                            *
*                  NATIONAL GEODETIC SURVEY, NOS, NOAA                *
*                  ROCKVILLE, MD   20852                              *
*                                                                     *
*                                                                     *
*                  DISCLAIMER                                         *
*                                                                     *
*   THIS PROGRAM AND SUPPORTING INFORMATION IS FURNISHED BY THE       *
* GOVERNMENT OF THE UNITED STATES OF AMERICA, AND IS ACCEPTED AND     *
* USED BY THE RECIPIENT WITH THE UNDERSTANDING THAT THE UNITED STATES *
* GOVERNMENT MAKES NO WARRANTIES, EXPRESS OR IMPLIED, CONCERNING THE  *
* ACCURACY, COMPLETENESS, RELIABILITY, OR SUITABILITY OF THIS         *
* PROGRAM, OF ITS CONSTITUENT PARTS, OR OF ANY SUPPORTING DATA.       *
*                                                                     *
*   THE GOVERNMENT OF THE UNITED STATES OF AMERICA SHALL BE UNDER NO  *
* LIABILITY WHATSOEVER RESULTING FROM ANY USE OF THIS PROGRAM.  THIS  *
* PROGRAM SHOULD NOT BE RELIED UPON AS THE SOLE BASIS FOR SOLVING A   *
* PROBLEM WHOSE INCORRECT SOLUTION COULD RESULT IN INJURY TO PERSON   *
* OR PROPERTY.                                                        *
*                                                                     *
*   THIS PROGRAM IS PROPERTY OF THE GOVERNMENT OF THE UNITED STATES   *
* OF AMERICA.  THEREFORE, THE RECIPIENT FURTHER AGREES NOT TO ASSERT  *
* PROPRIETARY RIGHTS THEREIN AND NOT TO REPRESENT THIS PROGRAM TO     *
* ANYONE AS BEING OTHER THAN A GOVERNMENT PROGRAM.                    *
*                                                                     *
***********************************************************************
*
*** PARAMETER STATEMENTS ARE SET-UP THRU OUT THE CODE TO DIMENSION
*** THE SIZE OF VECTORS.
*
* PARAMETER  DEFINITION
* --------------------------------------------------------------------
*   LENA     WORK ARRAY SIZE IN SINGLE PRECISION WORDS
*   LENIW    WORK ARRAY SIZE IN DOUBLE PRECISION WORDS =2*LENA
*            IN GENERAL, IF THIS PROGRAM NEEDS ADDITIONAL STORAGE
*            ALLOCATION, ONLY THESE TWO PARAMETERS NEED TO BE INCREASED.
*            UNFORTUNATELY, ONE MUST THEN RE-COMPILE AND LINK
*            BEFORE RE-EXECUTING.
*   MXSSN    MAXINUM NUMBER OF SSN'S
*   NVECS    MAXIMUM NUMBER OF GPS VECTORS IN A GROUP OF GPS OBS (GROUP = session)
*   MAXZ     MAXINUM NUMBER OF DIRECTION ROTATION UNKNOWNS
*            (I.E., THE MAXIMUM NUMBER OF *20* BLUE BOOK RECORDS
*             WHICH IS, THE MAXUMUM NUMBER OF DIRECTION LISTS)
*   MXPRM    MAXIMUM NUMBER OF AUXILLARY PARAMETERS
*   MXVF     MAXIMUM NUMBER OF VARIANCE FACTORS
*   MXRD     MAXIMUM NUMBER OF 'GREATEST RESIDUALS' LISTED IN OUTPUT
*   LENC     LENGTH OF THE C (PARTIAL DERIVITIVE) AND THE IC VECTORS
*   MXFOT    MAXIMUM NUMBER OF FOURTH ORDER STATIONS (LANDMARKS)
*   NSING    SIZE IF THE SINGULARITY ARRAY
*   MXALL    TOTAL POSIBLE PARAMETERS = MXPRM + 15 + MAXZ + MXSSN*3
*            NOTE THAT THE MAXIMUM NUMBER OF GPS ROTATION PARAMETER
*            SETS IS 5.  AND THERE ARE 3 PARAMETERS PER SET
*                                                                     *
***********************************************************************

*** SIMULTANEOUS NETWORK ADJUSTMENT OF BLUE BOOK (GPS)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

***   PARAMETER ( LENA =  2000000 )
***   PARAMETER ( LENA = 10000000 )
      PARAMETER ( LENA = 60000000 )
******
******
      PARAMETER ( LENIW = LENA+LENA )
      DIMENSION IW(LENIW), A(LENA)
      EQUIVALENCE ( A(1), IW(1) )

      CHARACTER*99 SCCSID
      CHARACTER*5  PGMVER
      CHARACTER*15 PGMDAT
      CHARACTER*80 BBOOK, AFILE, GFILE, DFILE, FILE, ADJOUT, BBNAM
      CHARACTER*26 TBUF
      CHARACTER*7  ADJFIL, NAFILE
      CHARACTER*1 RESP  
      COMMON /NAME/   BBOOK, AFILE, GFILE, DFILE, BBNAM, ADJFIL, NAFILE
      COMMON /PAGEIT/ MAXLIN, IPAGE, ILINE
      COMMON /VERSN/  PGMVER
      COMMON /CONST/  PI, PI2, RAD
      LOGICAL NLMODE,NLSCL
      COMMON /NLOPT/ NLMODE,NLSCL,NLUN1,NLUN2,NLUN3,NLUN4,NLUN5
      COMMON/UNITS/ LUNIT

*** THE FOLLOWING COMMON BLOCK IS HERE TO FORCE PERMANENT
*** STORAGE -- NECESSARY IN SOME MACHINES

      COMMON /BLOCK/ A

      LNWORK = LENIW
      LAWORK = LENA

*** DEFINE CONSTANTS
      SCCSID='@(#)adjust.for	6.4.2 	2017/11/06 NGS'
      PGMVER='6.4.2'
      PGMDAT='2017/11/06'

      PI2 = 2.D0*DATAN(1.D0)
      PI = 2.D0*PI2
      RAD = 180.D0/PI
      LUNIT = 16

*** DEFINE USER FILES; BLUE BOOK, ADJUSTMENT OPTION FILE, GPS FILE,
*** AND DOPPLER FILE

      BBOOK = 'BBOOK'
      AFILE = 'AFILE'
      GFILE = 'GFILE'
      DFILE = 'DFILE'

*** use with all compilers
      WRITE (*,2) PGMVER,PGMDAT
   2  FORMAT(' WELCOME TO ADJUST VERSION ',A5,
     *       '   DATE(ccyy/mm/dd) ',A15,/)

      WRITE (*,10)
   10 FORMAT (' ENTER INPUT BLUE BOOK FILENAME (DEFAULT=''BBOOK''):',/)
      READ (*,20,END=99,ERR=99) FILE
   20 FORMAT (A80)
      IF (FILE .NE. ' ') BBOOK = FILE
      WRITE (*,20) BBOOK

      WRITE (*,30)
***v 4.28k
***   30 FORMAT (' ENTER ADJUSTMENT FILE FILENAME', /,
***     *        ' (DEFAULT=''AFILE'', IF THERE ISNT ONE,',
***     *        ' ENTER: ''NOAFILE''):', /)
   30 FORMAT (' ENTER ADJUSTMENT FILE FILENAME',
     *        ' (DEFAULT=''AFILE''):'/) 
******************************************     
     
      READ (*,20,END=99,ERR=99) FILE
      IF (FILE .NE. ' ') AFILE = FILE
      WRITE (*,20) AFILE

      WRITE (*,40)
   40 FORMAT(' ENTER GPS FILE FILENAME (DEFAULT=''GFILE''):',/)
      READ (*,20,END=99,ERR=99) FILE
      IF (FILE .NE. ' ') GFILE = FILE
      WRITE (*,20) GFILE

c     WRITE (*,50)
   50 FORMAT (' ENTER DOPPLER FILE FILENAME', /,
     *        ' (DEFAULT=''DFILE'', IF THERE ISNT ONE,',
     *        ' ENTER: ''NODFILE''):', /)
c     READ (*,20,END=99,ERR=99) FILE
      FILE = 'NODFILE'
      IF (FILE .NE. ' ') DFILE = FILE
c     WRITE (*,20) DFILE

      WRITE (*,55)
   55 FORMAT (' ENTER THE DESIRED NAME FOR THE PRINTED OUTPUT : '/,
     *        ' (Make sure it is not an existing filename.)')
     
      READ(*,20) ADJOUT
      WRITE (*,20) ADJOUT
      
**v 4.28

      WRITE (*,56)
   56 FORMAT (' ENTER THE DESIRED NAME FOR THE OUTPUT BLUE BOOK',
     *        ' FILENAME : '/,
     *        ' (Make sure it is not an existing filename.)')
     
      READ(*,20) BBNAM
      WRITE (*,20) BBNAM

   99 CONTINUE
      CLOSE (UNIT=6)

********************************************************************************
*** use with PC compiler
*** The following statement is non-standard fortran!!
*** Use it only with the WATCOM Compiler!!
* 	!DEC$ IF DEFINED (WATCOM)
*     OPEN ( UNIT=LUNIT, CARRIAGECONTROL='YES', FILE=ADJOUT)
*	!DEC$ ELSE IF DEFINED (DF)
*     OPEN ( UNIT=LUNIT, CARRIAGECONTROL='FORTRAN', FILE=ADJOUT)
*	!DEC$ ENDIF

*** use with sun compiler
*      OPEN (LUNIT,FILE=ADJOUT,STATUS='NEW',
*     *ACCESS='SEQUENTIAL',FORM='FORMATTED',BLANK='ZERO')
********************************************************************************

      CALL GET_OFIL( LUNIT, ADJOUT )

*** DEFINE SCRATCH FILES

      IUO = 8
      OPEN (IUO, FORM='UNFORMATTED', STATUS='SCRATCH')
      IUO2 = 9
      OPEN (IUO2, FORM='UNFORMATTED', STATUS='SCRATCH')
      REWIND IUO
      ENDFILE IUO
      REWIND IUO
      REWIND IUO2
      ENDFILE IUO2
      REWIND IUO2

***  DEFINE UNITS FOR NETWORK AND LOCAL ACCURACY COMPUTATIONS
	NLUN1=18
	NLUN2=19
	NLUN3=20 
	NLUN4=21
	NLUN5=22 

*** SET PAGE PARAMETERS AND PRINT FIRST HEADING

      MAXLIN = 58
      IPAGE = 1
      CALL HEAD

*** WRITE FILENAMES TO OUTPUT

c     WRITE (LUNIT,100) BBOOK,AFILE,GFILE,DFILE,ADJOUT,BBNAM
      WRITE (LUNIT,100) BBOOK,AFILE,GFILE,ADJOUT,BBNAM
* 100  FORMAT (' BLUE BOOK FILE = ',A40/,' ADJUSTMENT FILE = ',
*     *A40/,' GPS FILE = ',A40/,' DOPPLER FILE = ',A40/,
*     *' OUTPUT FILE = ',A40/,' OUTPUT BLUE BOOK FILE = ',A40/)
 100  FORMAT (' BLUE BOOK FILE = ',A80/,' ADJUSTMENT FILE = ',
c    *A80/,' GPS FILE = ',A80/,' DOPPLER FILE = ',A80/,
     *A80/,' GPS FILE = ',A80/,
     *' OUTPUT FILE = ',A80/,' OUTPUT BLUE BOOK FILE = ',A80/)
     
      ILINE = ILINE + 6
            
*** GET SYSTEM TIME

      ITYPE = 0
      CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
      CALL LINE (2)
      WRITE (LUNIT, 1) TBUF(1:24)
    1 FORMAT ('0SYSTEM TIME IS ', A24)

*** READ ADJUSTMENT FILE AND ECHO OPTIONS

      CALL AFIL (A)
      CALL AFPRNT

*** FIRST PASS OF DATA

      CALL FIRST (A)

*** ALLOCATE STORAGE

      CALL ALOCAT (ID1, ID2, ID3, II4, ID5, LAWORK, LNWORK)

*** READ ADJUSTMENT FILE AGAIN FOR REMAINING RECORDS

      CALL AFILE2

*** PERFORM ADJUSTMENT

      CALL ADJST (A(1), A(ID1+1), A(ID2+1), A(ID3+1), IW(II4+1),
     *            A(ID5+1), LAWORK, LNWORK, IUO, IUO2)

*** END OF ADJUSTMENT

      ITYPE = 0
      CALL SYSTIM (ITYPE, TBUF, DIFM, DIFM0)
      CALL LINE (4)
      WRITE (LUNIT, 3) TBUF(1:24)
    3 FORMAT ('0END OF ADJUST PROCESSING',
     &        /, ' SYSTEM TIME IS ', A24,
     &        /, ' HAVE A NICE DAY!')
      CLOSE (IUO)
      CLOSE (IUO2)
      STOP
      END
