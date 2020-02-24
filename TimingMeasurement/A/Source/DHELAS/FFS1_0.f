C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma5(2,1)
C     
      SUBROUTINE FFS1_0(F1, F2, S3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 S3(*)
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP28
      COMPLEX*16 COUP
      TMP28 = (F2(5)*F1(5)+F2(6)*F1(6)-F2(3)*F1(3)-F2(4)*F1(4))
      VERTEX = COUP*(-CI * TMP28*S3(3))
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma5(2,1)
C     
      SUBROUTINE FFS1_2_0(F1, F2, S3, COUP1, COUP2,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP2
      COMPLEX*16 S3(*)
      COMPLEX*16 F1(*)
      COMPLEX*16 COUP1
      COMPLEX*16 F2(*)
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP
      CALL FFS1_0(F1,F2,S3,COUP1,VERTEX)
      CALL FFS2_0(F1,F2,S3,COUP2,TMP)
      VERTEX = VERTEX + TMP
      END

