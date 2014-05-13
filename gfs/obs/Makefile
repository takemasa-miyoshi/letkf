include ../configure.user

PROGS   = superob dec_prepbufr dec_prcp readdiag_conv

OBJS    = obs_tools.o superob_tools.o obs_gfs_ext.o

COMMON_DIR = ../common

COMMON_OBJS = $(COMMON_DIR)/SFMT.o           \
              $(COMMON_DIR)/common.o         \
              $(COMMON_DIR)/common_mtx.o     \
              $(COMMON_DIR)/netlib.o         \
              $(COMMON_DIR)/common_letkf.o   \
              $(COMMON_DIR)/common_gfs.o     \
              $(COMMON_DIR)/common_obs_gfs.o

LIBS    = -L$(SP_LIB) -lsp_4 $(LIBBLAS)

all:	$(PROGS)

superob:	superob.o $(OBJS) common_modules
	$(SFC) $(FOPTS) $(FMCMODEL) -o superob superob.o $(COMMON_OBJS) $(OBJS) \
	$(LIBS) -L$(SIGIO_LIB) -lsigio_4
dec_prepbufr:	dec_prepbufr.o common_modules
	$(SFC) $(FOPTS) -o dec_prepbufr dec_prepbufr.o $(COMMON_OBJS) \
	$(LIBS) -L$(BUFR_LIB) -lbufrlib
dec_prcp:	dec_prcp.o common_modules
	$(SFC) $(FOPTS) -o dec_prcp dec_prcp.o $(COMMON_OBJS) $(LIBS)
readdiag_conv:	readdiag_conv.o common_modules
	$(SFC) $(FOPTS2) -o readdiag_conv readdiag_conv.o $(COMMON_OBJS) $(OBJS) \
	$(LIBS)

superob.o:      superob.f90 $(OBJS) common_modules
	$(SFC) $(FOPTS) $(FMCMODEL) -c superob.f90 -I$(COMMON_DIR) -I$(SIGIO_INC)
dec_prepbufr.o:	dec_prepbufr.f90 common_modules
	$(SFC) $(FOPTS) -c dec_prepbufr.f90 -I$(COMMON_DIR)
dec_prcp.o:	dec_prcp.f90 common_modules
	$(SFC) $(FOPTS) -c dec_prcp.f90 -I$(COMMON_DIR)
readdiag_conv.o:	readdiag_conv.f90 $(OBJS) common_modules
	$(SFC) $(FOPTS2) -c readdiag_conv.f90 -I$(COMMON_DIR)

obs_tools.o:	obs_tools.f90 common_modules
	$(SFC) $(FOPTS) -c obs_tools.f90 -I$(COMMON_DIR)
superob_tools.o:	superob_tools.f90 obs_tools.o obs_gfs_ext.o common_modules
	$(SFC) $(FOPTS) -c superob_tools.f90 -I$(COMMON_DIR)
obs_gfs_ext.o:	obs_gfs_ext.f90 common_modules
	$(SFC) $(FOPTS) -c obs_gfs_ext.f90 -I$(COMMON_DIR)

common_modules:
	( cd $(COMMON_DIR) ; make )

clean:
	rm -f *.o *.mod $(PROGS)
