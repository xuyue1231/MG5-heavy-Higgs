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
      GC_6 = -G
      GC_7 = MDL_COMPLEXI*G
      GC_8 = MDL_COMPLEXI*MDL_G__EXP__2
      GC_10 = -(MDL_CA*MDL_COMPLEXI*MDL_GHGG*MDL_KHGG)
      GC_11 = -(MDL_CA*G*MDL_GHGG*MDL_KHGG)
      GC_12 = MDL_CA*MDL_COMPLEXI*MDL_G__EXP__2*MDL_GHGG*MDL_KHGG
      GC_13 = -(MDL_CA*MDL_COMPLEXI*MDL_GHHGG*MDL_KHHGG)
      GC_62 = (MDL_COMPLEXI*MDL_GAAGG*MDL_KAAGG*MDL_SA)/8.000000D+00
      GC_64 = (MDL_COMPLEXI*MDL_GAGG*MDL_KAGG*MDL_SA)/8.000000D+00
      GC_65 = (G*MDL_GAGG*MDL_KAGG*MDL_SA)/4.000000D+00
      END
