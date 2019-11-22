ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Parameters that should not be recomputed event by event.
C     
      IF(READLHA) THEN

        G = 2 * DSQRT(AS*PI)  ! for the first init

        MDL_COS__CABI = COS(MDL_CABI)

        MDL_CKM1X1 = MDL_COS__CABI

        MDL_SIN__CABI = SIN(MDL_CABI)

        MDL_CKM1X2 = MDL_SIN__CABI

        MDL_CKM2X1 = -MDL_SIN__CABI

        MDL_CKM2X2 = MDL_COS__CABI

        MDL_CA__EXP__2 = MDL_CA**2

        MDL_SA = SQRT(DCMPLX(1.000000D+00-MDL_CA__EXP__2))

        MDL_COMPLEXI = DCMPLX(0.000000D+00,1.000000D+00)

        MDL_KHDW = MDL_COMPLEXI*MDL_KHDWI+MDL_KHDWR

        MDL_MZ__EXP__2 = MDL_MZ**2

        MDL_MZ__EXP__4 = MDL_MZ**4

        MDL_SQRT__2 = SQRT(DCMPLX(2.000000D+00))

        MDL_NB__2__EXP__0_75 = 2.000000D+00**7.500000D-01

        MDL_MX0__EXP__2 = MDL_MX0**2

        MDL_CA__EXP__4 = MDL_CA**4

        MDL_CA__EXP__3 = MDL_CA**3

        MDL_CONJG__CKM1X1 = CONJG(DCMPLX(MDL_CKM1X1))

        MDL_CONJG__CKM1X2 = CONJG(DCMPLX(MDL_CKM1X2))

        MDL_CONJG__CKM2X1 = CONJG(DCMPLX(MDL_CKM2X1))

        MDL_CONJG__CKM2X2 = CONJG(DCMPLX(MDL_CKM2X2))

        MDL_CONJG__KHDW = CONJG(DCMPLX(MDL_KHDW))

        MDL_AEW = 1.000000D+00/AEWM1

        MDL_MW = SQRT(DCMPLX(MDL_MZ__EXP__2/2.000000D+00
     $   +SQRT(DCMPLX(MDL_MZ__EXP__4/4.000000D+00-(MDL_AEW*PI
     $   *MDL_MZ__EXP__2)/(MDL_GF*MDL_SQRT__2)))))

        MDL_SQRT__AEW = SQRT(DCMPLX(MDL_AEW))

        MDL_EE = 2.000000D+00*MDL_SQRT__AEW*SQRT(DCMPLX(PI))

        MDL_MW__EXP__2 = MDL_MW**2

        MDL_SW2 = 1.000000D+00-MDL_MW__EXP__2/MDL_MZ__EXP__2

        MDL_CW = SQRT(DCMPLX(1.000000D+00-MDL_SW2))

        MDL_SQRT__SW2 = SQRT(DCMPLX(MDL_SW2))

        MDL_SW = MDL_SQRT__SW2

        MDL_AD = (MDL_EE*(-5.000000D-01+(2.000000D+00*MDL_SW2)
     $   /3.000000D+00))/(2.000000D+00*MDL_CW*MDL_SW)

        MDL_AL = (MDL_EE*(-5.000000D-01+2.000000D+00*MDL_SW2))
     $   /(2.000000D+00*MDL_CW*MDL_SW)

        MDL_AN = MDL_EE/(4.000000D+00*MDL_CW*MDL_SW)

        MDL_AU = (MDL_EE*(5.000000D-01-(4.000000D+00*MDL_SW2)
     $   /3.000000D+00))/(2.000000D+00*MDL_CW*MDL_SW)

        MDL_BD = -MDL_EE/(4.000000D+00*MDL_CW*MDL_SW)

        MDL_BL = -MDL_EE/(4.000000D+00*MDL_CW*MDL_SW)

        MDL_BN = MDL_EE/(4.000000D+00*MDL_CW*MDL_SW)

        MDL_BU = MDL_EE/(4.000000D+00*MDL_CW*MDL_SW)

        MDL_GWWZ = -((MDL_CW*MDL_EE)/MDL_SW)

        MDL_G1 = MDL_EE/MDL_CW

        MDL_GW = MDL_EE/MDL_SW

        MDL_VEV = (2.000000D+00*MDL_MW*MDL_SW)/MDL_EE

        MDL_EE__EXP__2 = MDL_EE**2

        MDL_GAAA = MDL_EE__EXP__2/(3.000000D+00*PI**2*MDL_VEV)

        MDL_VEV__EXP__2 = MDL_VEV**2

        MDL_CW__EXP__2 = MDL_CW**2

        MDL_GAZA = ((-5.000000D+00+8.000000D+00*MDL_CW__EXP__2)
     $   *SQRT(DCMPLX(MDL_EE__EXP__2*MDL_GF*MDL_MZ__EXP__2)))
     $   /(6.000000D+00*MDL_NB__2__EXP__0_75*PI**2*MDL_VEV)

        MDL_GHAA = (4.700000D+01*MDL_EE__EXP__2)/(7.200000D+01*PI**2
     $   *MDL_VEV)

        MDL_GHZA = ((-1.300000D+01+9.400000D+01*MDL_CW__EXP__2)
     $   *SQRT(DCMPLX(MDL_EE__EXP__2*MDL_GF*MDL_MZ__EXP__2)))
     $   /(3.600000D+01*MDL_NB__2__EXP__0_75*PI**2*MDL_VEV)

        MDL_LAM = MDL_MX0__EXP__2/(2.000000D+00*MDL_VEV__EXP__2)

        MDL_YB = (MDL_YMB*MDL_SQRT__2)/MDL_VEV

        MDL_YT = (MDL_YMT*MDL_SQRT__2)/MDL_VEV

        MDL_YTAU = (MDL_YMTAU*MDL_SQRT__2)/MDL_VEV

        MDL_MUH = SQRT(DCMPLX(MDL_LAM*MDL_VEV__EXP__2))

        MDL_SW__EXP__2 = MDL_SW**2

      ENDIF
C     
C     Parameters that should be recomputed at an event by even basis.
C     
      AS = G**2/4/PI

      MDL_SQRT__AS = SQRT(DCMPLX(AS))

      MDL_G__EXP__2 = G**2

      MDL_GAAGG = MDL_G__EXP__2/(8.000000D+00*PI**2*MDL_VEV__EXP__2)

      MDL_GAGG = MDL_G__EXP__2/(8.000000D+00*PI**2*MDL_VEV)

      MDL_GHGG = -MDL_G__EXP__2/(1.200000D+01*PI**2*MDL_VEV)

      MDL_GHHGG = MDL_G__EXP__2/(1.200000D+01*PI**2*MDL_VEV__EXP__2)

C     
C     Definition of the EW coupling used in the write out of aqed
C     
      GAL(1) = 3.5449077018110318D0 / DSQRT(AEWM1)
      GAL(2) = 1D0

