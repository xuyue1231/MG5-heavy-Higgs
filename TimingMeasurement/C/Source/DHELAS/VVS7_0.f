C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1,2)*P(2,2) - P(-1,2)**2*Metric(1,2)
C     
      SUBROUTINE VVS7_0(V1, V2, S3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 S3(*)
      COMPLEX*16 COUP
      COMPLEX*16 TMP10
      REAL*8 P2(0:3)
      COMPLEX*16 TMP17
      COMPLEX*16 TMP15
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP25
      COMPLEX*16 V1(*)
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      TMP15 = (P2(0)*V1(3)-P2(1)*V1(4)-P2(2)*V1(5)-P2(3)*V1(6))
      TMP25 = (P2(0)*P2(0)-P2(1)*P2(1)-P2(2)*P2(2)-P2(3)*P2(3))
      TMP17 = (P2(0)*V2(3)-P2(1)*V2(4)-P2(2)*V2(5)-P2(3)*V2(6))
      TMP10 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      VERTEX = COUP*S3(3)*(-CI*(TMP15*TMP17)+CI*(TMP10*TMP25))
      END

