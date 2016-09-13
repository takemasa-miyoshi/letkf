include ../configure.user

PROGS   = obsope obsmake #dec_prepbufr #superob dec_prcp readdiag_conv

OBJS    = obsope_tools.o #obs_tools.o superob_tools.o obs_gfs_ext.o

COMMON_DIR = ../common

COMMON_OBJS = $(COMMON_DIR)/SFMT.o           \
              $(COMMON_DIR)/common.o         \
              $(COMMON_DIR)/common_mpi.o     \
              $(COMMON_DIR)/common_mtx.o     \
              $(COMMON_DIR)/netlib.o         \
              $(COMMON_DIR)/common_letkf.o   \
              $(COMMON_DIR)/common_ncio.o    \
              $(COMMON_DIR)/common_scale.o     \
              $(COMMON_DIR)/common_obs_scale.o \
              $(COMMON_DIR)/common_mpi_scale.o \
              $(COMMON_DIR)/common_nml.o
#              $(COMMON_DIR)/mod_admin_time.o
#              $(COMMON_DIR)/common_scalelib.o
#              $(COMMON_DIR)/common_precip.o

INCS    = $(SCALE_INC)

LIBS    = $(SCALE_LIB) $(NETCDF_LIB) $(LAPACK_LIB)

ifeq ($(H08),T)
COMMON_OBJS += $(COMMON_DIR)/scale_H08_fwd.o
LIBS    += $(RT_INCLUDE) $(RT_LDFLAGS)
endif

all:	$(PROGS)

obsope:	obsope.o $(OBJS) common_modules
	$(FC) $(FOPTS) $(FMCMODEL) -o obsope obsope.o \
	$(COMMON_OBJS) $(OBJS) $(LIBS)
obsmake:	obsmake.o $(OBS) common_modules
	$(FC) $(FOPTS) $(FMCMODEL) -o obsmake obsmake.o $(COMMON_OBJS) $(OBJS) $(LIBS)
superob:	superob.o $(OBJS) common_modules
	$(SFC) $(FOPTS) $(FMCMODEL) -o superob superob.o $(COMMON_OBJS) $(OBJS) \
	$(LIBS) -L$(SIGIO_LIB) -lsigio_4
dec_prepbufr:	dec_prepbufr.o common_modules
	$(FC) $(FOPTS) -o dec_prepbufr dec_prepbufr.o $(COMMON_OBJS) \
	$(LIBS) $(BUFR_LIB)
dec_prcp:	dec_prcp.o common_modules
	$(SFC) $(FOPTS) -o dec_prcp dec_prcp.o $(COMMON_OBJS) $(LIBS)
readdiag_conv:	readdiag_conv.o common_modules
	$(SFC) $(FOPTS) -o readdiag_conv readdiag_conv.o $(COMMON_OBJS) $(OBJS) \
	$(LIBS)

obsope.o:	obsope.f90 $(OBJS) common_modules
	$(FC) $(FOPTS) -c obsope.f90 -I$(COMMON_DIR) $(INCS)
obsmake.o:	obsmake.f90 $(OBJS) common_modules
	$(FC) $(FOPTS) -c obsmake.f90 -I$(COMMON_DIR) $(INCS)
superob.o:      superob.f90 $(OBJS) common_modules
	$(SFC) $(FOPTS) $(FMCMODEL) -c superob.f90 -I$(COMMON_DIR) -I$(SIGIO_INC)
dec_prepbufr.o:	dec_prepbufr.f90 common_modules
	$(FC) $(FOPTS) -c dec_prepbufr.f90 -I$(COMMON_DIR) $(INCS)
dec_prcp.o:	dec_prcp.f90 common_modules
	$(SFC) $(FOPTS) -c dec_prcp.f90 -I$(COMMON_DIR)
readdiag_conv.o:	readdiag_conv.f90 $(OBJS) common_modules
	$(SFC) $(FOPTS) -c readdiag_conv.f90 -I$(COMMON_DIR)

obsope_tools.o:	obsope_tools.f90 common_modules
	$(FC) $(FOPTS) -c obsope_tools.f90 -I$(COMMON_DIR) $(INCS)
obs_tools.o:	obs_tools.f90 common_modules
	$(SFC) $(FOPTS) -c obs_tools.f90 -I$(COMMON_DIR)
superob_tools.o:	superob_tools.f90 obs_tools.o obs_gfs_ext.o common_modules
	$(SFC) $(FOPTS) -c superob_tools.f90 -I$(COMMON_DIR)
obs_gfs_ext.o:	obs_gfs_ext.f90 common_modules
	$(SFC) $(FOPTS) -c obs_gfs_ext.f90 -I$(COMMON_DIR)

common_modules:
	( cd $(COMMON_DIR) ; make )

clean:
	rm -f *.o *.mod *.lst $(PROGS) obsope
