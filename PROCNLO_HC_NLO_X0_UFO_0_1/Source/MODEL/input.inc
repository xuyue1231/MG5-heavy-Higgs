ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DOUBLE PRECISION MDL_SQRT__AS,MDL_G__EXP__2,MDL_G__EXP__4
     $ ,MDL_R2MIXEDFACTOR_FIN_,MDL_GWCFT_UV_B_1EPS_
     $ ,MDL_GWCFT_UV_T_1EPS_,MDL_BWCFT_UV_1EPS_,MDL_TWCFT_UV_1EPS_
     $ ,MDL_G__EXP__3,MDL_MU_R__EXP__2,MDL_GWCFT_UV_B_FIN_
     $ ,MDL_GWCFT_UV_T_FIN_,MDL_BWCFT_UV_FIN_,MDL_TWCFT_UV_FIN_
     $ ,MDL_GAGG,MDL_GHGG,MDL_CKM11,MDL_CONJG__CKM3X3,MDL_CONJG__CKM11
     $ ,MDL_LHV,MDL_CKM3X3,MDL_CONJG__CKM33,MDL_NCOL,MDL_CA,MDL_TF
     $ ,MDL_CF,MDL_MZ__EXP__2,MDL_MZ__EXP__4,MDL_SQRT__2
     $ ,MDL_MX0__EXP__2,MDL_COSA__EXP__2,MDL_SINA,MDL_NB__2__EXP__0_75
     $ ,MDL_NCOL__EXP__2,MDL_MB__EXP__2,MDL_MT__EXP__2,MDL_AEW
     $ ,MDL_SQRT__AEW,MDL_EE,MDL_MW__EXP__2,MDL_SW2,MDL_CW
     $ ,MDL_SQRT__SW2,MDL_SW,MDL_G1,MDL_GW,MDL_VEV,MDL_VEV__EXP__2
     $ ,MDL_LAM,MDL_YB,MDL_YT,MDL_YTAU,MDL_MUH,MDL_AXIALZUP
     $ ,MDL_AXIALZDOWN,MDL_VECTORZUP,MDL_VECTORZDOWN,MDL_VECTORAUP
     $ ,MDL_VECTORADOWN,MDL_VECTORWMDXU,MDL_AXIALWMDXU,MDL_VECTORWPUXD
     $ ,MDL_AXIALWPUXD,MDL_EE__EXP__2,MDL_GAAA,MDL_CW__EXP__2,MDL_GAZA
     $ ,MDL_GHAA,MDL_GHZA,MDL_GW__EXP__2,MDL_SW__EXP__2,AEWM1,MDL_GF
     $ ,AS,MDL_YMB,MDL_YMT,MDL_YMTAU,MDL_LAMBDA,MDL_COSA,MDL_KSM
     $ ,MDL_KHTT,MDL_KATT,MDL_KHBB,MDL_KABB,MDL_KHLL,MDL_KALL,MDL_KHAA
     $ ,MDL_KAAA,MDL_KHZA,MDL_KAZA,MDL_KHZZ,MDL_KAZZ,MDL_KHWW,MDL_KAWW
     $ ,MDL_KHDA,MDL_KHDZ,MDL_KHDWR,MDL_KHDWI

      COMMON/PARAMS_R/ MDL_SQRT__AS,MDL_G__EXP__2,MDL_G__EXP__4
     $ ,MDL_R2MIXEDFACTOR_FIN_,MDL_GWCFT_UV_B_1EPS_
     $ ,MDL_GWCFT_UV_T_1EPS_,MDL_BWCFT_UV_1EPS_,MDL_TWCFT_UV_1EPS_
     $ ,MDL_G__EXP__3,MDL_MU_R__EXP__2,MDL_GWCFT_UV_B_FIN_
     $ ,MDL_GWCFT_UV_T_FIN_,MDL_BWCFT_UV_FIN_,MDL_TWCFT_UV_FIN_
     $ ,MDL_GAGG,MDL_GHGG,MDL_CKM11,MDL_CONJG__CKM3X3,MDL_CONJG__CKM11
     $ ,MDL_LHV,MDL_CKM3X3,MDL_CONJG__CKM33,MDL_NCOL,MDL_CA,MDL_TF
     $ ,MDL_CF,MDL_MZ__EXP__2,MDL_MZ__EXP__4,MDL_SQRT__2
     $ ,MDL_MX0__EXP__2,MDL_COSA__EXP__2,MDL_SINA,MDL_NB__2__EXP__0_75
     $ ,MDL_NCOL__EXP__2,MDL_MB__EXP__2,MDL_MT__EXP__2,MDL_AEW
     $ ,MDL_SQRT__AEW,MDL_EE,MDL_MW__EXP__2,MDL_SW2,MDL_CW
     $ ,MDL_SQRT__SW2,MDL_SW,MDL_G1,MDL_GW,MDL_VEV,MDL_VEV__EXP__2
     $ ,MDL_LAM,MDL_YB,MDL_YT,MDL_YTAU,MDL_MUH,MDL_AXIALZUP
     $ ,MDL_AXIALZDOWN,MDL_VECTORZUP,MDL_VECTORZDOWN,MDL_VECTORAUP
     $ ,MDL_VECTORADOWN,MDL_VECTORWMDXU,MDL_AXIALWMDXU,MDL_VECTORWPUXD
     $ ,MDL_AXIALWPUXD,MDL_EE__EXP__2,MDL_GAAA,MDL_CW__EXP__2,MDL_GAZA
     $ ,MDL_GHAA,MDL_GHZA,MDL_GW__EXP__2,MDL_SW__EXP__2,AEWM1,MDL_GF
     $ ,AS,MDL_YMB,MDL_YMT,MDL_YMTAU,MDL_LAMBDA,MDL_COSA,MDL_KSM
     $ ,MDL_KHTT,MDL_KATT,MDL_KHBB,MDL_KABB,MDL_KHLL,MDL_KALL,MDL_KHAA
     $ ,MDL_KAAA,MDL_KHZA,MDL_KAZA,MDL_KHZZ,MDL_KAZZ,MDL_KHWW,MDL_KAWW
     $ ,MDL_KHDA,MDL_KHDZ,MDL_KHDWR,MDL_KHDWI


      DOUBLE COMPLEX MDL_COMPLEXI,MDL_KHDW,MDL_CONJG__KHDW,MDL_I1X33
     $ ,MDL_I2X33,MDL_I3X33,MDL_I4X33,MDL_VECTOR_TBGP,MDL_AXIAL_TBGP
     $ ,MDL_VECTOR_TBGM,MDL_AXIAL_TBGM

      COMMON/PARAMS_C/ MDL_COMPLEXI,MDL_KHDW,MDL_CONJG__KHDW,MDL_I1X33
     $ ,MDL_I2X33,MDL_I3X33,MDL_I4X33,MDL_VECTOR_TBGP,MDL_AXIAL_TBGP
     $ ,MDL_VECTOR_TBGM,MDL_AXIAL_TBGM

