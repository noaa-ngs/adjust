C %P%

      SUBROUTINE TOLGH2(VX,VY,VZ,VN,VE,VU,ISN, JSN, B, G, IVEC,NR,NC)
********************************************************************************
* COMPUTE RESIDUALS IN E,N,U. FOR A GPS VECTOR
*   PUT V-TILDE H IN COLUMN N2+1 OF G
*   PUT V-TILDE U IN COLUMN N2+2 OF G
********************************************************************************
      PARAMETER(NVECS=700)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(*)
      DIMENSION RM(3,3)
      DIMENSION G(NR,NC)
      COMMON/UNITS/ LUNIT
      COMMON /GPSHUS/ GPSVS(NVECS*3,2), GPSRNS(NVECS*3)
      CHARACTER*99 SCCSID

C      SCCSID='$Id: tolgh2.for 44197 2010-07-13 12:26:03Z bruce.tran $	20$Date: 2008/08/25 16:03:50 $ NGS'

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
      RM(2,1)=-SL
      RM(3,1)=CB*CL
      RM(1,2)=-SB*SL
      RM(2,2)=CL 
      RM(3,2)=CB*SL
      RM(1,3)=CB
      RM(2,3)=0
      RM(3,3)=SB

*** 
      VN=RM(1,1)*VX+RM(1,2)*VY+RM(1,3)*VZ
      VE=RM(2,1)*VX+RM(2,2)*VY+RM(2,3)*VZ
      VU=RM(3,1)*VX+RM(3,2)*VY+RM(3,3)*VZ


C 2/4/03

C  save V-tilde sub H
      I=3*IVEC

      GPSVS(I-2,1)=RM(1,1)*VN+RM(2,1)*VE
      GPSVS(I-1,1)=RM(1,2)*VN+RM(2,2)*VE
      GPSVS(I  ,1)=RM(1,3)*VN+RM(2,3)*VE
C  save V-tilde sub U
      GPSVS(I-2,2)=RM(3,1)*VU
      GPSVS(I-1,2)=RM(3,2)*VU
      GPSVS(I  ,2)=RM(3,3)*VU
C

      RETURN
      END

