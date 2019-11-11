C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1,1)*P(2,1) - P(-1,1)*P(-1,1)*Metric(1,2)
C     
      SUBROUTINE VVS11_3(V1, V2, COUP, M3, W3,S3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP2
      COMPLEX*16 S3(5)
      COMPLEX*16 TMP1
      REAL*8 W3
      COMPLEX*16 TMP0
      REAL*8 M3
      COMPLEX*16 P3(0:3)
      COMPLEX*16 DENOM
      COMPLEX*16 P1(0:3)
      COMPLEX*16 COUP
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP3
      P1(0) = V1(1)
      P1(1) = V1(2)
      P1(2) = V1(3)
      P1(3) = V1(4)
      S3(1) = +V1(1)+V2(1)
      S3(2) = +V1(2)+V2(2)
      S3(3) = +V1(3)+V2(3)
      S3(4) = +V1(4)+V2(4)
      P3(0) = -S3(1)
      P3(1) = -S3(2)
      P3(2) = -S3(3)
      P3(3) = -S3(4)
      TMP1 = (P1(0)*V1(5)-P1(1)*V1(6)-P1(2)*V1(7)-P1(3)*V1(8))
      TMP0 = (P1(0)*V2(5)-P1(1)*V2(6)-P1(2)*V2(7)-P1(3)*V2(8))
      TMP3 = (P1(0)*P1(0)-P1(1)*P1(1)-P1(2)*P1(2)-P1(3)*P1(3))
      TMP2 = (V2(5)*V1(5)-V2(6)*V1(6)-V2(7)*V1(7)-V2(8)*V1(8))
      DENOM = COUP/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 -CI
     $ * W3))
      S3(5)= DENOM*(-CI*(TMP2*TMP3)+CI*(TMP0*TMP1))
      END

