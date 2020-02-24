ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_1 = -(MDL_EE*MDL_COMPLEXI)/3.000000D+00
      GC_2 = (2.000000D+00*MDL_EE*MDL_COMPLEXI)/3.000000D+00
      GC_16 = -(MDL_CA*MDL_COMPLEXI*MDL_GHZA*MDL_KHZA)
      GC_66 = MDL_COMPLEXI*MDL_GAZA*MDL_KAZA*MDL_SA
      GC_72 = (MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_73 = (MDL_CKM1X1*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_74 = (MDL_CKM1X2*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_75 = (MDL_CKM2X1*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_76 = (MDL_CKM2X2*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_77 = -(MDL_CW*MDL_EE*MDL_COMPLEXI)/(2.000000D+00*MDL_SW)
      GC_78 = (MDL_CW*MDL_EE*MDL_COMPLEXI)/(2.000000D+00*MDL_SW)
      GC_89 = -(MDL_EE*MDL_COMPLEXI*MDL_SW)/(6.000000D+00*MDL_CW)
      GC_90 = (MDL_EE*MDL_COMPLEXI*MDL_SW)/(2.000000D+00*MDL_CW)
      GC_99 = (MDL_CA*MDL_EE__EXP__2*MDL_COMPLEXI*MDL_KSM*MDL_VEV)
     $ /(2.000000D+00*MDL_SW__EXP__2)
      GC_100 = MDL_CA*MDL_EE__EXP__2*MDL_COMPLEXI*MDL_KSM*MDL_VEV
     $ +(MDL_CA*MDL_CW__EXP__2*MDL_EE__EXP__2*MDL_COMPLEXI*MDL_KSM
     $ *MDL_VEV)/(2.000000D+00*MDL_SW__EXP__2)+(MDL_CA*MDL_EE__EXP__2
     $ *MDL_COMPLEXI*MDL_KSM*MDL_SW__EXP__2*MDL_VEV)/(2.000000D+00
     $ *MDL_CW__EXP__2)
      GC_101 = -((MDL_CA*MDL_COMPLEXI*MDL_KHBB*MDL_YB)/MDL_SQRT__2)
      GC_102 = (MDL_KABB*MDL_SA*MDL_YB)/MDL_SQRT__2
      GC_105 = -((MDL_CA*MDL_COMPLEXI*MDL_KHLL*MDL_YTAU)/MDL_SQRT__2)
      GC_106 = (MDL_KALL*MDL_SA*MDL_YTAU)/MDL_SQRT__2
      GC_107 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM1X1)/(MDL_SW
     $ *MDL_SQRT__2)
      GC_108 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM1X2)/(MDL_SW
     $ *MDL_SQRT__2)
      GC_109 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM2X1)/(MDL_SW
     $ *MDL_SQRT__2)
      GC_110 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM2X2)/(MDL_SW
     $ *MDL_SQRT__2)
      GC_38 = -((MDL_CA*MDL_COMPLEXI*MDL_KHDA)/MDL_LAMBDA)
      GC_39 = -((MDL_CA*MDL_COMPLEXI*MDL_KHDW)/MDL_LAMBDA)
      END