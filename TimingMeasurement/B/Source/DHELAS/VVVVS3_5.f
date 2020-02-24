C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)
C     
      SUBROUTINE VVVVS3_5(V1, V2, V3, V4, COUP, M5, W5,S5)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP12
      COMPLEX*16 V3(*)
      COMPLEX*16 S5(3)
      COMPLEX*16 V4(*)
      COMPLEX*16 TMP15
      REAL*8 M5
      COMPLEX*16 DENOM
      REAL*8 W5
      REAL*8 P5(0:3)
      COMPLEX*16 COUP
      COMPLEX*16 TMP19
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP13
      S5(1) = +V1(1)+V2(1)+V3(1)+V4(1)
      S5(2) = +V1(2)+V2(2)+V3(2)+V4(2)
      P5(0) = -DBLE(S5(1))
      P5(1) = -DBLE(S5(2))
      P5(2) = -DIMAG(S5(2))
      P5(3) = -DIMAG(S5(1))
      TMP15 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      TMP19 = (V3(3)*V1(3)-V3(4)*V1(4)-V3(5)*V1(5)-V3(6)*V1(6))
      TMP13 = (V3(3)*V4(3)-V3(4)*V4(4)-V3(5)*V4(5)-V3(6)*V4(6))
      TMP12 = (V2(3)*V4(3)-V2(4)*V4(4)-V2(5)*V4(5)-V2(6)*V4(6))
      DENOM = COUP/(P5(0)**2-P5(1)**2-P5(2)**2-P5(3)**2 - M5 * (M5 -CI
     $ * W5))
      S5(3)= DENOM*(-CI*(TMP13*TMP15)+CI*(TMP12*TMP19))
      END

