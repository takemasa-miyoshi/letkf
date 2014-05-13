include ../configure.user

PROGS   = ss2grd grd2ss sscycle ss2grdp grdctl

OBJS    = ssio_tools.o physcons.o

COMMON_DIR = ../common

COMMON_OBJS = $(COMMON_DIR)/SFMT.o            \
              $(COMMON_DIR)/common.o          \
              $(COMMON_DIR)/common_gfs.o      \
              $(COMMON_DIR)/common_gfs_pres.o

INCS    = -I$(SIGIO_INC) \
          -I$(SFCIO_INC)

LIBS    = -L$(BACIO_LIB) -lbacio_4 \
          -L$(SIGIO_LIB) -lsigio_4 \
          -L$(SFCIO_LIB) -lsfcio_4 \
          -L$(SP_LIB) -lsp_4

all:	$(PROGS)

ss2grd:	ss2grd.o $(OBJS) common_modules
	$(SFC) $(FOPTS) -o ss2grd ss2grd.o $(COMMON_OBJS) $(OBJS) $(LIBS)
grd2ss:	grd2ss.o $(OBJS) common_modules
	$(SFC) $(FOPTS) -o grd2ss grd2ss.o $(COMMON_OBJS) $(OBJS) $(LIBS)
sscycle:	sscycle.o $(OBJS) common_modules
	$(SFC) $(FOPTS) -o sscycle sscycle.o $(COMMON_OBJS) $(OBJS) $(LIBS)
ss2grdp:	ss2grdp.o $(OBJS) common_modules
	$(SFC) $(FOPTS) -o ss2grdp ss2grdp.o $(COMMON_OBJS) $(OBJS) $(LIBS)
grdctl:	grdctl.o $(OBJS) common_modules
	$(SFC) $(FOPTS) -o grdctl grdctl.o $(COMMON_OBJS) $(OBJS) $(LIBS)

ss2grd.o:	ss2grd.f90 $(OBJS) common_modules
	$(SFC) $(FOPTS) -c ss2grd.f90 -I$(COMMON_DIR) $(INCS)
grd2ss.o:	grd2ss.f90 $(OBJS) common_modules
	$(SFC) $(FOPTS) -c grd2ss.f90 -I$(COMMON_DIR) $(INCS)
sscycle.o:	sscycle.f90 $(OBJS) common_modules
	$(SFC) $(FOPTS) -c sscycle.f90 -I$(COMMON_DIR) $(INCS)
ss2grdp.o:	ss2grdp.f90 $(OBJS) common_modules
	$(SFC) $(FOPTS) -c ss2grdp.f90 -I$(COMMON_DIR) $(INCS)
grdctl.o:	grdctl.f90 $(OBJS) common_modules
	$(SFC) $(FOPTS) -c grdctl.f90 -I$(COMMON_DIR) $(INCS)

physcons.o:	common_modules
	$(SFC) $(FOPTS) -c physcons.f90 -I$(COMMON_DIR)
ssio_tools.o:	common_modules physcons.o
	$(SFC) $(FOPTS) -c ssio_tools.f90 -I$(COMMON_DIR)

common_modules:
	( cd $(COMMON_DIR) ; make )

clean:
	rm -f *.o *.mod $(PROGS)
