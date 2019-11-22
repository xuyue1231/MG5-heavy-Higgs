ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP3()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_10 = -(MDL_CA*MDL_COMPLEXI*MDL_GHGG*MDL_KHGG)
      GC_13 = -(MDL_CA*MDL_COMPLEXI*MDL_GHHGG*MDL_KHHGG)
      GC_62 = (MDL_COMPLEXI*MDL_GAAGG*MDL_KAAGG*MDL_SA)/8.000000D+00
      GC_64 = (MDL_COMPLEXI*MDL_GAGG*MDL_KAGG*MDL_SA)/8.000000D+00
      END
