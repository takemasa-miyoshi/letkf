include ../configure.user

OBJS    = SFMT.o common.o common_mpi.o common_mtx.o netlib.o netlibblas.o common_letkf.o \
	  common_gfs.o common_gfs_pres.o common_obs_gfs.o common_mpi_gfs.o common_precip.o

all:	$(OBJS) enssize datetime

enssize:        enssize.o SFMT.o common.o common_mtx.o netlib.o common_letkf.o
	$(SFC) $(FOPTS) -o enssize enssize.o SFMT.o common.o common_mtx.o netlib.o common_letkf.o \
	$(LIBBLAS) -L$(SP_LIB) -lsp_4
datetime:	datetime.o SFMT.o common.o
	$(SFC) $(FOPTS) -o datetime datetime.o SFMT.o common.o

common.o:	../../common/common.f90
	$(SFC) $(FOPTS) -c ../../common/common.f90
common_mpi.o:	../../common/common_mpi.f90 common.o
	$(FC) $(FOPTS) -c ../../common/common_mpi.f90
common_mtx.o:	../../common/common_mtx.f90 common.o
	$(SFC) $(FOPTS) -c ../../common/common_mtx.f90
SFMT.o: ../../common/SFMT.f90
	$(SFC) $(FOPTS) -c ../../common/SFMT.f90
netlib.o:	../../common/netlib.f
	$(SFC) $(FOPTS) -c ../../common/netlib.f
netlibblas.o:	../../common/netlibblas.f
	$(SFC) $(FOPTS) -c ../../common/netlibblas.f
common_letkf.o:	../../common/common_letkf.f90 common.o common_mtx.o
	$(SFC) $(FOPTS) -c ../../common/common_letkf.f90
common_gfs.o:	common_gfs.f90 common.o
	$(SFC) $(FOPTS) -c common_gfs.f90
common_gfs_pres.o:   common_gfs_pres.f90 common.o common_gfs.o
	$(SFC) $(FOPTS) -c common_gfs_pres.f90
common_obs_gfs.o:	common_obs_gfs.f90 common.o common_gfs.o
	$(SFC) $(FOPTS) -c common_obs_gfs.f90
common_mpi_gfs.o:	common_mpi_gfs.f90 common.o common_mpi.o common_gfs.o common_obs_gfs.o
	$(FC) $(FOPTS) -c common_mpi_gfs.f90
common_precip.o:	common_precip.f90 common.o common_gfs.o common_letkf.o
	$(SFC) $(FOPTS) -c common_precip.f90
enssize.o:	enssize.f90 common_letkf.o
	$(SFC) $(FOPTS) -c enssize.f90
datetime.o:	datetime.f90 common.o
	$(SFC) $(FOPTS) -c datetime.f90

clean:
	rm -f *.o *.mod enssize datetime
