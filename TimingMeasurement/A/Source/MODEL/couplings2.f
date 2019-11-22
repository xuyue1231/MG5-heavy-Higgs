ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP2()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_40 = -((MDL_CA*MDL_COMPLEXI*MDL_KHDZ)/MDL_LAMBDA)
      GC_42 = -((MDL_CA*MDL_COMPLEXI*MDL_KHZZ)/MDL_LAMBDA)
      GC_43 = (MDL_COMPLEXI*MDL_KL)/(4.000000D+00*MDL_LAMBDA)
      GC_45 = (MDL_COMPLEXI*MDL_KQ)/(4.000000D+00*MDL_LAMBDA)
      GC_49 = (MDL_COMPLEXI*MDL_KQ3)/(4.000000D+00*MDL_LAMBDA)
      GC_56 = -((MDL_COMPLEXI*MDL_KQ3*MDL_MB)/MDL_LAMBDA)
      GC_58 = -((MDL_COMPLEXI*MDL_KL*MDL_MTA)/MDL_LAMBDA)
      GC_68 = (MDL_COMPLEXI*MDL_KAZZ*MDL_SA)/(2.000000D+00*MDL_LAMBDA)
      END
