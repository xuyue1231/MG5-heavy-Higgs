C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
Coup(1) * (Epsilon(1,2,-1,-2)*P(-1,1)*P(-2,2)) + Coup(2) * (Metric(1,2)) + Coup(3) * (P(1,1)*P(2,1) - P(-1,1)*P(-1,1)*Metric(1,2)) + Coup(4) * (P(1,2)*P(2,1) - P(-1,1)*P(-1,2)*Metric(1,2)) + Coup(5) * (P(1,2)*P(2,2) - P(-1,2)*P(-1,2)*Metric(1,2))
C     
      SUBROUTINE MP_VVS10_1_11_2_12_2(V1, S3, COUP1, COUP2, COUP3,
     $  COUP4, COUP5, M2, W2,V2)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 V2(8)
      COMPLEX*32 DENOM
      COMPLEX*32 S3(*)
      COMPLEX*32 TMP1
      COMPLEX*32 COUP2
      COMPLEX*32 COUP5
      REAL*16 W2
      COMPLEX*32 P2(0:3)
      COMPLEX*32 COUP1
      COMPLEX*32 TMP4
      COMPLEX*32 TMP5
      COMPLEX*32 P1(0:3)
      REAL*16 OM2
      REAL*16 M2
      COMPLEX*32 COUP3
      COMPLEX*32 COUP4
      COMPLEX*32 V1(*)
      COMPLEX*32 TMP3
      COMPLEX*32 TMP8
      P1(0) = V1(1)
      P1(1) = V1(2)
      P1(2) = V1(3)
      P1(3) = V1(4)
      OM2 = 0Q0
      IF (M2.NE.0Q0) OM2=1Q0/M2**2
      V2(1) = +V1(1)+S3(1)
      V2(2) = +V1(2)+S3(2)
      V2(3) = +V1(3)+S3(3)
      V2(4) = +V1(4)+S3(4)
      P2(0) = -V2(1)
      P2(1) = -V2(2)
      P2(2) = -V2(3)
      P2(3) = -V2(4)
      TMP5 = (P2(0)*V1(5)-P2(1)*V1(6)-P2(2)*V1(7)-P2(3)*V1(8))
      TMP4 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      TMP1 = (P1(0)*V1(5)-P1(1)*V1(6)-P1(2)*V1(7)-P1(3)*V1(8))
      TMP8 = (P2(0)*P2(0)-P2(1)*P2(1)-P2(2)*P2(2)-P2(3)*P2(3))
      TMP3 = (P1(0)*P1(0)-P1(1)*P1(1)-P1(2)*P1(2)-P1(3)*P1(3))
      DENOM = 1Q0/(P2(0)**2-P2(1)**2-P2(2)**2-P2(3)**2 - M2 * (M2 -CI*
     $  W2))
      V2(5)= DENOM*CI * S3(5)*(COUP1*(P1(1)*(P2(2)*V1(8)-P2(3)*V1(7))
     $ +(P1(2)*(P2(3)*V1(6)-P2(1)*V1(8))+P1(3)*(P2(1)*V1(7)-P2(2)*V1(6)
     $ )))+(P2(0)*(OM2*(TMP5*(COUP2-TMP3*COUP3)+TMP1*TMP4*COUP3)-TMP5
     $ *COUP5)+(V1(5)*(TMP3*COUP3+TMP4*COUP4+TMP8*COUP5-COUP2)-P1(0)
     $ *(TMP1*COUP3+TMP5*COUP4))))
      V2(6)= DENOM*(-CI )* S3(5)*(COUP1*(P1(0)*(P2(3)*V1(7)-P2(2)*V1(8)
     $ )+(P1(2)*(P2(0)*V1(8)-P2(3)*V1(5))+P1(3)*(P2(2)*V1(5)-P2(0)
     $ *V1(7))))+(P2(1)*(OM2*(TMP5*(TMP3*COUP3-COUP2)-TMP1*TMP4*COUP3)
     $ +TMP5*COUP5)+(V1(6)*(-1Q0)*(TMP3*COUP3+TMP4*COUP4+TMP8*COUP5
     $ -COUP2)+P1(1)*(TMP1*COUP3+TMP5*COUP4))))
      V2(7)= DENOM*(-CI )* S3(5)*(COUP1*(P1(0)*(P2(1)*V1(8)-P2(3)*V1(6)
     $ )+(P1(1)*(P2(3)*V1(5)-P2(0)*V1(8))+P1(3)*(P2(0)*V1(6)-P2(1)
     $ *V1(5))))+(P2(2)*(OM2*(TMP5*(TMP3*COUP3-COUP2)-TMP1*TMP4*COUP3)
     $ +TMP5*COUP5)+(V1(7)*(-1Q0)*(TMP3*COUP3+TMP4*COUP4+TMP8*COUP5
     $ -COUP2)+P1(2)*(TMP1*COUP3+TMP5*COUP4))))
      V2(8)= DENOM*(-CI )* S3(5)*(COUP1*(P1(0)*(P2(2)*V1(6)-P2(1)*V1(7)
     $ )+(P1(1)*(P2(0)*V1(7)-P2(2)*V1(5))+P1(2)*(P2(1)*V1(5)-P2(0)
     $ *V1(6))))+(P2(3)*(OM2*(TMP5*(TMP3*COUP3-COUP2)-TMP1*TMP4*COUP3)
     $ +TMP5*COUP5)+(V1(8)*(-1Q0)*(TMP3*COUP3+TMP4*COUP4+TMP8*COUP5
     $ -COUP2)+P1(3)*(TMP1*COUP3+TMP5*COUP4))))
      END

