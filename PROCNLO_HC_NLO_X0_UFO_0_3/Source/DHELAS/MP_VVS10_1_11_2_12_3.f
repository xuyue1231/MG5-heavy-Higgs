C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
Coup(1) * (Epsilon(1,2,-1,-2)*P(-1,1)*P(-2,2)) + Coup(2) * (Metric(1,2)) + Coup(3) * (P(1,1)*P(2,1) - P(-1,1)*P(-1,1)*Metric(1,2)) + Coup(4) * (P(1,2)*P(2,1) - P(-1,1)*P(-1,2)*Metric(1,2)) + Coup(5) * (P(1,2)*P(2,2) - P(-1,2)*P(-1,2)*Metric(1,2))
C     
      SUBROUTINE MP_VVS10_1_11_2_12_3(V1, V2, COUP1, COUP2, COUP3,
     $  COUP4, COUP5, M3, W3,S3)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 S3(5)
      COMPLEX*32 TMP1
      COMPLEX*32 P3(0:3)
      COMPLEX*32 COUP2
      COMPLEX*32 TMP2
      COMPLEX*32 P2(0:3)
      COMPLEX*32 COUP1
      COMPLEX*32 COUP3
      COMPLEX*32 COUP4
      COMPLEX*32 V2(*)
      COMPLEX*32 TMP12
      REAL*16 W3
      COMPLEX*32 COUP5
      COMPLEX*32 TMP7
      COMPLEX*32 DENOM
      COMPLEX*32 P1(0:3)
      COMPLEX*32 V1(*)
      COMPLEX*32 TMP5
      COMPLEX*32 TMP0
      REAL*16 M3
      COMPLEX*32 TMP4
      COMPLEX*32 TMP3
      COMPLEX*32 TMP8
      P1(0) = V1(1)
      P1(1) = V1(2)
      P1(2) = V1(3)
      P1(3) = V1(4)
      P2(0) = V2(1)
      P2(1) = V2(2)
      P2(2) = V2(3)
      P2(3) = V2(4)
      S3(1) = +V1(1)+V2(1)
      S3(2) = +V1(2)+V2(2)
      S3(3) = +V1(3)+V2(3)
      S3(4) = +V1(4)+V2(4)
      P3(0) = -S3(1)
      P3(1) = -S3(2)
      P3(2) = -S3(3)
      P3(3) = -S3(4)
      TMP8 = (P2(0)*P2(0)-P2(1)*P2(1)-P2(2)*P2(2)-P2(3)*P2(3))
      TMP5 = (P2(0)*V1(5)-P2(1)*V1(6)-P2(2)*V1(7)-P2(3)*V1(8))
      TMP4 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      TMP7 = (P2(0)*V2(5)-P2(1)*V2(6)-P2(2)*V2(7)-P2(3)*V2(8))
      TMP2 = (V2(5)*V1(5)-V2(6)*V1(6)-V2(7)*V1(7)-V2(8)*V1(8))
      TMP1 = (P1(0)*V1(5)-P1(1)*V1(6)-P1(2)*V1(7)-P1(3)*V1(8))
      TMP0 = (P1(0)*V2(5)-P1(1)*V2(6)-P1(2)*V2(7)-P1(3)*V2(8))
      TMP3 = (P1(0)*P1(0)-P1(1)*P1(1)-P1(2)*P1(2)-P1(3)*P1(3))
      TMP12 = (-1Q0)*(P1(0)*(P2(1)*(V2(7)*V1(8)-V2(8)*V1(7))+(P2(2)
     $ *(V2(8)*V1(6)-V2(6)*V1(8))+P2(3)*(V2(6)*V1(7)-V2(7)*V1(6))))
     $ +(P1(1)*(P2(0)*(V2(8)*V1(7)-V2(7)*V1(8))+(P2(2)*(V2(5)*V1(8)
     $ -V2(8)*V1(5))+P2(3)*(V2(7)*V1(5)-V2(5)*V1(7))))+(P1(2)*(P2(0)
     $ *(V2(6)*V1(8)-V2(8)*V1(6))+(P2(1)*(V2(8)*V1(5)-V2(5)*V1(8))
     $ +P2(3)*(V2(5)*V1(6)-V2(6)*V1(5))))+P1(3)*(P2(0)*(V2(7)*V1(6)
     $ -V2(6)*V1(7))+(P2(1)*(V2(5)*V1(7)-V2(7)*V1(5))+P2(2)*(V2(6)
     $ *V1(5)-V2(5)*V1(6)))))))
      DENOM = 1Q0/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 -CI*
     $  W3))
      S3(5)= DENOM*(TMP2*(-1Q0)*(-CI*(COUP2)+CI*(TMP3*COUP3+TMP4*COUP4
     $ +TMP8*COUP5))+(TMP0*(+CI*(TMP1*COUP3+TMP5*COUP4))+(+CI*(COUP1
     $ *TMP12+TMP5*TMP7*COUP5))))
      END

