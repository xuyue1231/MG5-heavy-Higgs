include ../../Source/make_opts
SHELL = /bin/bash
LIBDIR = ../../lib/
PROG   = check
PROG_SPLITORDERS = check_sa_born_splitOrders
LINKLIBS =  -L../../lib/ -ldhelas -lmodel
LIBS = $(LIBDIR)libdhelas.$(libext) $(LIBDIR)libmodel.$(libext)
PROCESS=  matrix.o
CHECK_SA=  check_sa.o
CHECK_SA_SPLITORDERS=  check_sa_born_splitOrders.o
 
$(PROG): $(PROCESS) $(CHECK_SA) makefile $(LIBS)
	$(FC) $(FFLAGS) -o $(PROG) $(PROCESS) $(CHECK_SA) $(LINKLIBS)

$(PROG_SPLITORDERS): $(PROCESS) $(CHECK_SA_SPLITORDERS) makefile $(LIBS)
	$(FC) $(FFLAGS) -o $(PROG) $(PROCESS) $(CHECK_SA_SPLITORDERS) $(LINKLIBS)

driver.f: nexternal.inc pmass.inc ngraphs.inc coupl.inc

# For python linking (require f2py part of numpy)
ifeq ($(origin MENUM),undefined)
  MENUM=2
endif

matrix$(MENUM)py.so: matrix.f makefile $(LIBS)
	touch __init__.py
	$(F2PY) $(LINKLIBS) -c  matrix.f -m matrix$(MENUM)py
# For python linking (require f2py part of numpy)
ifeq ($(origin MENUM),undefined)
  MENUM=2
endif
P1_uux_zhh_z_emep_hh_uuuxux/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_uuuxux;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_uuuxux/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_uuuxux;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_ucuxcx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_ucuxcx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_ucuxcx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_ucuxcx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_uduxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_uduxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_uduxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_uduxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_uduxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_uduxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_uduxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_uduxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_udcxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_udcxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_udcxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_udcxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_udcxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_udcxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_udcxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_udcxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_usuxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_usuxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_usuxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_usuxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_usuxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_usuxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_usuxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_usuxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_uscxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_uscxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_uscxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_uscxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_uscxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_uscxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_uscxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_uscxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_cccxcx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_cccxcx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_cccxcx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_cccxcx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_cduxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_cduxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_cduxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_cduxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_cduxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_cduxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_cduxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_cduxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_cdcxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_cdcxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_cdcxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_cdcxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_cdcxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_cdcxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_cdcxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_cdcxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_csuxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_csuxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_csuxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_csuxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_csuxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_csuxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_csuxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_csuxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_cscxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_cscxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_cscxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_cscxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_cscxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_cscxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_cscxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_cscxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_dddxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_dddxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_dddxdx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_dddxdx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_dsdxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_dsdxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_dsdxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_dsdxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_emep_hh_sssxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_emep_hh_sssxsx;make matrix$(MENUM)py.so
P1_uux_zhh_z_mummup_hh_sssxsx/matrix$(MENUM)py.so:
	cd P1_uux_zhh_z_mummup_hh_sssxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_uuuxux/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_uuuxux;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_uuuxux/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_uuuxux;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_ucuxcx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_ucuxcx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_ucuxcx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_ucuxcx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_uduxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_uduxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_uduxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_uduxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_uduxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_uduxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_uduxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_uduxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_udcxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_udcxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_udcxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_udcxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_udcxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_udcxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_udcxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_udcxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_usuxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_usuxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_usuxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_usuxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_usuxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_usuxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_usuxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_usuxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_uscxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_uscxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_uscxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_uscxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_uscxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_uscxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_uscxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_uscxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_cccxcx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_cccxcx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_cccxcx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_cccxcx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_cduxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_cduxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_cduxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_cduxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_cduxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_cduxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_cduxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_cduxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_cdcxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_cdcxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_cdcxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_cdcxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_cdcxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_cdcxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_cdcxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_cdcxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_csuxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_csuxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_csuxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_csuxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_csuxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_csuxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_csuxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_csuxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_cscxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_cscxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_cscxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_cscxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_cscxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_cscxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_cscxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_cscxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_dddxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_dddxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_dddxdx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_dddxdx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_dsdxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_dsdxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_dsdxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_dsdxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_emep_hh_sssxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_emep_hh_sssxsx;make matrix$(MENUM)py.so
P1_ddx_zhh_z_mummup_hh_sssxsx/matrix$(MENUM)py.so:
	cd P1_ddx_zhh_z_mummup_hh_sssxsx;make matrix$(MENUM)py.so
all:  P1_uux_zhh_z_emep_hh_uuuxux/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_uuuxux/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_ucuxcx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_ucuxcx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_uduxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_uduxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_uduxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_uduxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_udcxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_udcxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_udcxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_udcxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_usuxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_usuxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_usuxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_usuxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_uscxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_uscxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_uscxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_uscxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_cccxcx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_cccxcx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_cduxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_cduxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_cduxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_cduxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_cdcxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_cdcxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_cdcxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_cdcxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_csuxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_csuxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_csuxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_csuxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_cscxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_cscxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_cscxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_cscxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_dddxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_dddxdx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_dsdxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_dsdxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_emep_hh_sssxsx/matrix$(MENUM)py.so  P1_uux_zhh_z_mummup_hh_sssxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_uuuxux/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_uuuxux/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_ucuxcx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_ucuxcx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_uduxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_uduxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_uduxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_uduxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_udcxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_udcxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_udcxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_udcxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_usuxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_usuxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_usuxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_usuxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_uscxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_uscxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_uscxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_uscxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_cccxcx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_cccxcx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_cduxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_cduxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_cduxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_cduxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_cdcxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_cdcxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_cdcxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_cdcxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_csuxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_csuxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_csuxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_csuxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_cscxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_cscxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_cscxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_cscxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_dddxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_dddxdx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_dsdxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_dsdxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_emep_hh_sssxsx/matrix$(MENUM)py.so  P1_ddx_zhh_z_mummup_hh_sssxsx/matrix$(MENUM)py.so 
	echo 'done'