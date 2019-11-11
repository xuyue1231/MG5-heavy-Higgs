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
      R2_DXUW = ((MDL_CKM11*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2))
     $ *MDL_R2MIXEDFACTOR_FIN_
      GC_5 = MDL_COMPLEXI*G
      END
