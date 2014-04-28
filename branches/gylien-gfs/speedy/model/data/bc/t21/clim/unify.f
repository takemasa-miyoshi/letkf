      program leggiscrivi

      parameter(nlon=64,nlat=32,ntime=12)
      parameter(undef=-9.999e+19)

      real sst(nlon,nlat,ntime),ssw1(nlon,nlat,ntime),
     *     ssw2(nlon,nlat,ntime),ssw3(nlon,nlat,ntime)
      real lst(nlon,nlat,ntime),lsw1(nlon,nlat,ntime),
     *     lsw2(nlon,nlat,ntime),lsw3(nlon,nlat,ntime)
      real st(nlon,nlat,ntime),sw1(nlon,nlat,ntime),
     *     sw2(nlon,nlat,ntime),sw3(nlon,nlat,ntime)

      character*80 filename

c read sst

      filename = 'fcsfh8190clim.t21.sea'
      open(10,file=filename,form='unformatted')

      do it = 1,ntime
        read(10) ((sst(i,j,it),i=1,nlon),j=1,nlat)
        read(10) ((ssw1(i,j,it),i=1,nlon),j=1,nlat)
        read(10) ((ssw2(i,j,it),i=1,nlon),j=1,nlat)
        read(10) ((ssw3(i,j,it),i=1,nlon),j=1,nlat)
      end do

      close(10)

c read lst

      filename = 'fcsfh8190clim.t21.land'
      open(20,file=filename,form='unformatted')

      do it = 1,ntime
        read(20) ((lst(i,j,it),i=1,nlon),j=1,nlat)
        read(20) ((lsw1(i,j,it),i=1,nlon),j=1,nlat)
        read(20) ((lsw2(i,j,it),i=1,nlon),j=1,nlat)
        read(20) ((lsw3(i,j,it),i=1,nlon),j=1,nlat)
      end do

      close(20)

c unify fields

      do it = 1,ntime
        do j = 1,nlat
          do i = 1,nlon
            st(i,j,it)  = sst(i,j,it)
            sw1(i,j,it) = ssw1(i,j,it)
            sw2(i,j,it) = ssw2(i,j,it)
            sw3(i,j,it) = ssw3(i,j,it)
            if(st(i,j,it).eq.undef) st(i,j,it)=lst(i,j,it)
            if(sw1(i,j,it).eq.undef) sw1(i,j,it)=lsw1(i,j,it)
            if(sw2(i,j,it).eq.undef) sw2(i,j,it)=lsw2(i,j,it)
            if(sw3(i,j,it).eq.undef) sw3(i,j,it)=lsw3(i,j,it)
          end do
        end do
      end do


c write st

      filename = 'fcsfh8190clim.t21'
      open(30,file=filename,form='unformatted')

      do it = 1,ntime
        write(30) ((st(i,j,it),i=1,nlon),j=1,nlat)
        write(30) ((sw1(i,j,it),i=1,nlon),j=1,nlat)
        write(30) ((sw2(i,j,it),i=1,nlon),j=1,nlat)
        write(30) ((sw3(i,j,it),i=1,nlon),j=1,nlat)
      end do

      close(10)

      stop
      end
