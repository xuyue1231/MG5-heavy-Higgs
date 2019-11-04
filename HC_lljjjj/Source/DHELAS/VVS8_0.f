C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1,1)*P(2,1) + P(1,2)*P(2,2) - P(-1,1)**2*Metric(1,2) -
C      P(-1,2)**2*Metric(1,2)
C     
      SUBROUTINE VVS8_0(V1, V2, S3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 S3(*)
      REAL*8 P1(0:3)
      COMPLEX*16 TMP10
      REAL*8 P2(0:3)
      COMPLEX*16 TMP19
      COMPLEX*16 TMP17
      COMPLEX*16 TMP16
      COMPLEX*16 TMP4
      COMPLEX*16 VERTEX
      COMPLEX*16 COUP
      COMPLEX*16 TMP9
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP18
      P1(0) = DBLE(V1(1))
      P1(1) = DBLE(V1(2))
      P1(2) = DIMAG(V1(2))
      P1(3) = DIMAG(V1(1))
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      TMP9 = (P2(0)*V1(3)-P2(1)*V1(4)-P2(2)*V1(5)-P2(3)*V1(6))
      TMP18 = (P1(0)*P1(0)-P1(1)*P1(1)-P1(2)*P1(2)-P1(3)*P1(3))
      TMP4 = (P1(0)*V2(3)-P1(1)*V2(4)-P1(2)*V2(5)-P1(3)*V2(6))
      TMP17 = (P2(0)*V2(3)-P2(1)*V2(4)-P2(2)*V2(5)-P2(3)*V2(6))
      TMP16 = (P1(0)*V1(3)-P1(1)*V1(4)-P1(2)*V1(5)-P1(3)*V1(6))
      TMP19 = (P2(0)*P2(0)-P2(1)*P2(1)-P2(2)*P2(2)-P2(3)*P2(3))
      TMP10 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      VERTEX = COUP*S3(3)*(TMP10*(+CI*(TMP18+TMP19))+(-CI*(TMP4*TMP16
     $ +TMP9*TMP17)))
      END

