C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFV2_0(F1, F2, V3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V3(*)
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 VERTEX
      COMPLEX*16 COUP
      COMPLEX*16 TMP9
      TMP9 = (F1(5)*(F2(7)*(V3(5)+V3(8))+F2(8)*(V3(6)+CI*(V3(7))))
     $ +F1(6)*(F2(7)*(V3(6)-CI*(V3(7)))+F2(8)*(V3(5)-V3(8))))
      VERTEX = COUP*(-CI * TMP9)
      END

