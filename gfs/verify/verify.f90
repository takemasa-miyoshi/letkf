!===============================================================================
! verify: main program
!-------------------------------------------------------------------------------
!
!
!
!===============================================================================
program verify
!-------------------------------------------------------------------------------

  use common
  use common_gfs
  use common_gfs_pres
  use common_obs_gfs
  implicit none

  character(8) :: vrffile = 'fcst.grd'
  character(9) :: vrfpfile = 'fcstp.grd'
  character(9) :: obsfile  = 'obsXX.dat'
  character(9) :: anafile = 'anaXX.grd'

  character(12) :: vrfobsfile  = 'vrfobsXX.dat'
  character(12) :: vrfanafile = 'vrfanaXX.dat'


  integer, parameter :: nvrf_obs = 1
  integer, parameter :: nvrf_ana = 1

  integer, parameter :: narea = 4
!--------------------------------------------- Globe,     NH,     TR,     SH
  real(r_size), parameter :: vlon1(narea) = (/  0.d0,   0.d0,   0.d0,   0.d0/)
  real(r_size), parameter :: vlon2(narea) = (/360.d0, 360.d0, 360.d0, 360.d0/)
  real(r_size), parameter :: vlat1(narea) = (/-90.d0,  20.d0, -20.d0, -90.d0/)
  real(r_size), parameter :: vlat2(narea) = (/ 90.d0,  90.d0,  20.d0, -20.d0/)

  real(r_sngl), parameter :: missing = -9.99e33 ! missing value in grads files

  integer, parameter :: nv3dv = 4 ! u,v,t,q
  integer, parameter :: nv2dv = 1 ! ps
  
  real(r_size), parameter :: threshold_dz = 500.0d0 ! threshold of surface height correction (m)
  real(r_size), parameter :: threshold_dp = 0.5d0   ! threshold of pressure level (Pa)

  real(r_size) :: v3d(nlon,nlat,nlev,nv3dx)
  real(r_size) :: v2d(nlon,nlat,nv2dx)

  real(r_size) :: v3dp(nlon,nlat,nlevp,nv3dp)
  real(r_size) :: v2dp(nlon,nlat,nv2dp)

  integer :: num3d_obs(narea,nlevp,nv3dv)
  real(r_size) :: bias3d_obs(narea,nlevp,nv3dv)
  real(r_size) :: abse3d_obs(narea,nlevp,nv3dv)
  real(r_size) :: rmse3d_obs(narea,nlevp,nv3dv)
  real(r_sngl) :: vrf3d_obs(narea,nlevp,nv3dv*4)

  integer :: num2d_obs(narea,nv2dv)
  real(r_size) :: bias2d_obs(narea,nv2dv)
  real(r_size) :: abse2d_obs(narea,nv2dv)
  real(r_size) :: rmse2d_obs(narea,nv2dv)
  real(r_sngl) :: vrf2d_obs(narea,nv2dv*4)


  real(r_size) :: v3dpa(nlon,nlat,nlevp,nv3dp)
  real(r_size) :: v2dpa(nlon,nlat,nv2dp)


  integer :: num3d_ana(narea,nlevp,nv3dp)
  real(r_size) :: wei3d_ana(narea,nlevp,nv3dp)
  real(r_size) :: bias3d_ana(narea,nlevp,nv3dp)
  real(r_size) :: abse3d_ana(narea,nlevp,nv3dp)
  real(r_size) :: rmse3d_ana(narea,nlevp,nv3dp)
  real(r_sngl) :: vrf3d_ana(narea,nlevp,nv3dp*4)

  integer :: num2d_ana(narea,nv2dp)
  real(r_size) :: wei2d_ana(narea,nv2dp)
  real(r_size) :: bias2d_ana(narea,nv2dp)
  real(r_size) :: abse2d_ana(narea,nv2dp)
  real(r_size) :: rmse2d_ana(narea,nv2dp)
  real(r_sngl) :: vrf2d_ana(narea,nv2dp*4)


  real(r_size), allocatable :: obselm(:)
  real(r_size), allocatable :: obslon(:)
  real(r_size), allocatable :: obslat(:)
  real(r_size), allocatable :: obslev(:)
  real(r_size), allocatable :: obsdat(:)
  real(r_size), allocatable :: obserr(:)
  real(r_size), allocatable :: obstyp(:)

  real(r_size) :: ri, rj, rk
  real(r_size) :: dz, hx, dep, wei, latm1, latm2
  integer :: n, nobs, ivrf, k, vk, vid, iarea
  integer :: i, i1, i2, j, j1, j2
  
!  integer :: nobsused

  integer :: iunit, iolen, irec


  logical :: obtype_vrf(nobtype,nvrf_obs)

!-------------------------------------------------------------------------------

  obtype_vrf(:,:) = .false.

  obtype_vrf(1,1) = .true.
  obtype_vrf(8,1) = .true.
  obtype_vrf(9,1) = .true.
  obtype_vrf(10,1) = .true.

!  obtype_vrf(1,2) = .true.

!  obtype_vrf(1,3) = .true.

!  obtype_vrf(:,4) = .true.

!===============================================================================

  call set_common_gfs

!-------------------------------------------------------------------------------

  call read_grdx(vrffile,v3d,v2d)
  call read_grdp(vrfpfile,v3dp,v2dp)

!-------------------------------------------------------------------------------

  write (*, *)
  write (*, '(A)') '================================================'
  write (*, '(A)') ' Verification against observations'
  write (*, '(A)') '================================================'

!-------------------------------------------------------------------------------

  do ivrf = 1, nvrf_obs

    write (obsfile(4:5), '(I2.2)') ivrf
    write (vrfobsfile(7:8), '(I2.2)') ivrf

    write (*, '(A,I4)') 'obs set =', ivrf
    write (*, '(2A)') 'read obs file: ', obsfile

    call get_nobs(obsfile, 7, nobs)

    if (nobs == 0) cycle

    allocate ( obselm (nobs) )
    allocate ( obslon (nobs) )
    allocate ( obslat (nobs) )
    allocate ( obslev (nobs) )
    allocate ( obsdat (nobs) )
    allocate ( obserr (nobs) )
    allocate ( obstyp (nobs) )


    call read_obs(obsfile, nobs, &
                  obselm, obslon, obslat, obslev, obsdat, obserr, obstyp)


!    nobsused = 0

    num3d_obs = 0
    num2d_obs = 0
    bias3d_obs = 0.0d0
    bias2d_obs = 0.0d0
    abse3d_obs = 0.0d0
    abse2d_obs = 0.0d0
    rmse3d_obs = 0.0d0
    rmse2d_obs = 0.0d0
    

    do n = 1, nobs
    
      if (.not. obtype_vrf(nint(obstyp(n)),ivrf)) cycle

      select case (nint(obselm(n)))
        case (id_u_obs)
          vid = 1
        case (id_v_obs)
          vid = 2
        case (id_t_obs,id_tv_obs)
          vid = 3
        case (id_q_obs)
          vid = 4
        case (id_ps_obs)
          vid = 1
        case default
          cycle
      end select

      vk = 0
      if (nint(obselm(n)) <= 9999) then
        do k = 1, nlevp
          if (abs(obslev(n) - levp(k)) < threshold_dp) then
            vk = k
            cycle
          end if
        end do
        if (vk == 0) cycle
      end if

      if (obslon(n) <    0.0d0) obslon(n) = obslon(n) + 360.0d0
      if (obslon(n) >= 360.0d0) obslon(n) = obslon(n) - 360.0d0

      call phys2ijk(v3d(:,:,:,iv3d_p), obselm(n), &
                    obslon(n), obslat(n), obslev(n), ri, rj, rk)

      if (ceiling(ri) < 2 .or. nlon+1 < ceiling(ri)) cycle
      if (ceiling(rj) < 2 .or. nlat < ceiling(rj)) cycle

      if (ceiling(rk) > nlev) cycle
      if (ceiling(rk) < 2 .and. nint(obselm(n)) /= id_ps_obs) cycle
      
      if(nint(obselm(n)) == id_ps_obs .and. obsdat(n) < -100.0d0) cycle
      if(nint(obselm(n)) == id_ps_obs) then
        call itpl_2d(v2d(:,:,iv2d_orog),ri,rj,dz)
        rk = obslev(n) - dz
        if(abs(rk) > threshold_dz) then ! pressure adjustment threshold
!          print '(A)','* PS obs vertical adjustment beyond threshold'
!          print '(A,F10.2,A,F6.2,A,F6.2,A)','  dz=',rk,&
!           & ', (lon,lat)=(',obslon(n),',',obslat(n),')'
          cycle
        end if
      end if

      ! observational operator
      call Trans_XtoY(obselm(n), ri, rj, rk, v3d, v2d, hx)


      dep = hx - obsdat(n)

!write (30, '(2I6,6F10.3)') vk, nint(obselm(n)), ri, rj, rk, hx, obsdat(n), dep


      if (nint(obselm(n)) == id_ps_obs) dep = dep / 100.

!      nobsused = nobsused + 1

      do iarea = 1, narea
        if (vlon1(iarea) <= obslon(n) .and. obslon(n) <= vlon2(iarea) .and. &
            vlat1(iarea) <= obslat(n) .and. obslat(n) <= vlat2(iarea)) then
          if (vk == 0) then ! 2D variables
            num2d_obs(iarea,vid) = num2d_obs(iarea,vid) + 1
            bias2d_obs(iarea,vid) = bias2d_obs(iarea,vid) + dep
            abse2d_obs(iarea,vid) = abse2d_obs(iarea,vid) + abs(dep)
            rmse2d_obs(iarea,vid) = rmse2d_obs(iarea,vid) + dep*dep
          else
            num3d_obs(iarea,vk,vid) = num3d_obs(iarea,vk,vid) + 1
            bias3d_obs(iarea,vk,vid) = bias3d_obs(iarea,vk,vid) + dep
            abse3d_obs(iarea,vk,vid) = abse3d_obs(iarea,vk,vid) + abs(dep)
            rmse3d_obs(iarea,vk,vid) = rmse3d_obs(iarea,vk,vid) + dep*dep
          end if
        end if
      end do



    end do ! [ n = 1, nobs ]

    do vid = 1, nv3dv
      do k = 1, nlevp
        do iarea = 1, narea
          vrf3d_obs(iarea,k,(vid-1)*4+1) = real(num3d_obs(iarea,k,vid), r_sngl)
          if (num3d_obs(iarea,k,vid) == 0) then
            vrf3d_obs(iarea,k,(vid-1)*4+2) = missing
            vrf3d_obs(iarea,k,(vid-1)*4+3) = missing
            vrf3d_obs(iarea,k,(vid-1)*4+4) = missing
          else
            vrf3d_obs(iarea,k,(vid-1)*4+2) = real(bias3d_obs(iarea,k,vid)/num3d_obs(iarea,k,vid), r_sngl)
            vrf3d_obs(iarea,k,(vid-1)*4+3) = real(abse3d_obs(iarea,k,vid)/num3d_obs(iarea,k,vid), r_sngl)
            vrf3d_obs(iarea,k,(vid-1)*4+4) = real(sqrt(rmse3d_obs(iarea,k,vid)/num3d_obs(iarea,k,vid)), r_sngl)
          end if
        end do
      end do
    end do

    do vid = 1, nv2dv
      do iarea = 1, narea
        vrf2d_obs(iarea,(vid-1)*4+1) = real(num2d_obs(iarea,vid), r_sngl)
        if (num2d_obs(iarea,vid) == 0) then
          vrf2d_obs(iarea,(vid-1)*4+2) = missing
          vrf2d_obs(iarea,(vid-1)*4+3) = missing
          vrf2d_obs(iarea,(vid-1)*4+4) = missing
        else
          vrf2d_obs(iarea,(vid-1)*4+2) = real(bias2d_obs(iarea,vid)/num2d_obs(iarea,vid), r_sngl)
          vrf2d_obs(iarea,(vid-1)*4+3) = real(abse2d_obs(iarea,vid)/num2d_obs(iarea,vid), r_sngl)
          vrf2d_obs(iarea,(vid-1)*4+4) = real(sqrt(rmse2d_obs(iarea,vid)/num2d_obs(iarea,vid)), r_sngl)
        end if
      end do
    end do

    iunit=55
    inquire(iolength=iolen) iolen
    open(iunit,file=vrfobsfile,form='unformatted',access='direct',recl=narea*iolen)

    irec=1
    do n=1,nv3dv*4
      do k=1,nlevp
        write(iunit,rec=irec) vrf3d_obs(:,k,n)
        irec = irec + 1
      end do
    end do

    do n=1,nv2dv*4
      write(iunit,rec=irec) vrf2d_obs(:,n)
      irec = irec + 1
    end do

    close(iunit)


    deallocate ( obselm )
    deallocate ( obslon )
    deallocate ( obslat )
    deallocate ( obslev )
    deallocate ( obsdat )
    deallocate ( obserr )
    deallocate ( obstyp )

  end do ! [ ivrf = 1, nvrf_obs ]

!-------------------------------------------------------------------------------

  write (*, *)
  write (*, '(A)') '================================================'
  write (*, '(A)') ' Verification against model analyses'
  write (*, '(A)') '================================================'

!-------------------------------------------------------------------------------

  do ivrf = 1, nvrf_ana

    write (anafile(4:5), '(I2.2)') ivrf
    write (vrfanafile(7:8), '(I2.2)') ivrf

    write (*, '(A,I4)') 'anal set =', ivrf
    write (*, '(2A)') 'read ana file: ', anafile

    call read_grdp(anafile,v3dpa,v2dpa)

    num3d_ana = 0
    num2d_ana = 0
    wei3d_ana = 0.0d0
    wei2d_ana = 0.0d0
    bias3d_ana = 0.0d0
    bias2d_ana = 0.0d0
    abse3d_ana = 0.0d0
    abse2d_ana = 0.0d0
    rmse3d_ana = 0.0d0
    rmse2d_ana = 0.0d0


    do iarea = 1, narea


      do i = 1, nlon
        if (lon(i) >= vlon1(iarea)) then
          i1 = i
          exit
        end if
      end do
      do i = nlon, 1, -1
        if (lon(i) <= vlon2(iarea)) then
          i2 = i
          exit
        end if
      end do
      do j = 1, nlat
        if (lat(j) >= vlat1(iarea)) then
          j1 = j
          exit
        end if
      end do
      do j = nlat, 1, -1
        if (lat(j) <= vlat2(iarea)) then
          j2 = j
          exit
        end if
      end do


      do vid = 1, nv3dp
        do k = 1, nlevp
          do j = j1, j2
            if (j == 1) then
              latm1 = -90.0d0
            else
              latm1 = 0.5d0*(lat(j-1) + lat(j))
            end if
            if (j == nlat) then
              latm2 = 90.0d0
            else
              latm2 = 0.5d0*(lat(j) + lat(j+1))
            end if
            wei = abs(sin(latm2*pi/180.0d0) - sin(latm1*pi/180.0d0))
!            wei = abs(latm2-latm1) * cos(lat(j)*pi/180.0d0)

            do i = i1, i2
              if (v3dp(i,j,k,vid) /= missing .and. v3dpa(i,j,k,vid) /= missing) then
                dep = v3dp(i,j,k,vid) - v3dpa(i,j,k,vid)
                num3d_ana(iarea,k,vid) = num3d_ana(iarea,k,vid) + 1
                wei3d_ana(iarea,k,vid) = wei3d_ana(iarea,k,vid) + wei
                bias3d_ana(iarea,k,vid) = bias3d_ana(iarea,k,vid) + wei * dep
                abse3d_ana(iarea,k,vid) = abse3d_ana(iarea,k,vid) + wei * abs(dep)
                rmse3d_ana(iarea,k,vid) = rmse3d_ana(iarea,k,vid) + wei * dep*dep
              end if
            end do
          end do
        end do

        do j = j1, j2
          if (j == 1) then
            latm1 = -90.0d0
          else
            latm1 = 0.5d0*(lat(j-1) + lat(j))
          end if
          if (j == nlat) then
            latm2 = 90.0d0
          else
            latm2 = 0.5d0*(lat(j) + lat(j+1))
          end if
          wei = abs(sin(latm2*pi/180.0d0) - sin(latm1*pi/180.0d0))
!          wei = abs(latm2-latm1) * cos(lat(j)*pi/180.0d0)

          do i = i1, i2
            if (v2dp(i,j,vid) /= missing .and. v2dpa(i,j,vid) /= missing) then
              dep = v2dp(i,j,vid) - v2dpa(i,j,vid)
              num2d_ana(iarea,vid) = num2d_ana(iarea,vid) + 1
              wei2d_ana(iarea,vid) = wei2d_ana(iarea,vid) + wei
              bias2d_ana(iarea,vid) = bias2d_ana(iarea,vid) + wei * dep
              abse2d_ana(iarea,vid) = abse2d_ana(iarea,vid) + wei * abs(dep)
              rmse2d_ana(iarea,vid) = rmse2d_ana(iarea,vid) + wei * dep*dep
            end if
          end do
        end do
      end do


    end do ! [ iarea = 1, narea ]




    do vid = 1, nv3dp
      do k = 1, nlevp
        do iarea = 1, narea
          vrf3d_ana(iarea,k,(vid-1)*4+1) = real(num3d_ana(iarea,k,vid), r_sngl)
          if (num3d_ana(iarea,k,vid) == 0) then
            vrf3d_ana(iarea,k,(vid-1)*4+2) = missing
            vrf3d_ana(iarea,k,(vid-1)*4+3) = missing
            vrf3d_ana(iarea,k,(vid-1)*4+4) = missing
          else
            vrf3d_ana(iarea,k,(vid-1)*4+2) = real(bias3d_ana(iarea,k,vid)/wei3d_ana(iarea,k,vid), r_sngl)
            vrf3d_ana(iarea,k,(vid-1)*4+3) = real(abse3d_ana(iarea,k,vid)/wei3d_ana(iarea,k,vid), r_sngl)
            vrf3d_ana(iarea,k,(vid-1)*4+4) = real(sqrt(rmse3d_ana(iarea,k,vid)/wei3d_ana(iarea,k,vid)), r_sngl)
          end if
        end do
      end do
    end do

    do vid = 1, nv2dp
      do iarea = 1, narea
        vrf2d_ana(iarea,(vid-1)*4+1) = real(num2d_ana(iarea,vid), r_sngl)
        if (num2d_ana(iarea,vid) == 0) then
          vrf2d_ana(iarea,(vid-1)*4+2) = missing
          vrf2d_ana(iarea,(vid-1)*4+3) = missing
          vrf2d_ana(iarea,(vid-1)*4+4) = missing
        else
          vrf2d_ana(iarea,(vid-1)*4+2) = real(bias2d_ana(iarea,vid)/wei2d_ana(iarea,vid), r_sngl)
          vrf2d_ana(iarea,(vid-1)*4+3) = real(abse2d_ana(iarea,vid)/wei2d_ana(iarea,vid), r_sngl)
          vrf2d_ana(iarea,(vid-1)*4+4) = real(sqrt(rmse2d_ana(iarea,vid)/wei2d_ana(iarea,vid)), r_sngl)
        end if
      end do
    end do

    iunit=55
    inquire(iolength=iolen) iolen
    open(iunit,file=vrfanafile,form='unformatted',access='direct',recl=narea*iolen)

    irec=1
    do n=1,nv3dp*4
      do k=1,nlevp
        write(iunit,rec=irec) vrf3d_ana(:,k,n)
        irec = irec + 1
      end do
    end do

    do n=1,nv2dp*4
      write(iunit,rec=irec) vrf2d_ana(:,n)
      irec = irec + 1
    end do

    close(iunit)




  end do ! [ ivrf = 1, nvrf_ana ]

!-------------------------------------------------------------------------------
end program
!===============================================================================
