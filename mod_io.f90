
      module mod_io

      use mod_param
      use mod_mpimess
      use mod_varandtypes
      use mod_time

      implicit none
      private

      integer , dimension(4) , save :: icount , istart

      public :: read_config
      public :: read_init
      public :: out_init
      public :: qmax_init
      public :: chym_out
      public :: chym_wqmax
      public :: read_runoff
      public :: chym_rst
      public :: createfile , closefile
      public :: createfile_rst
      public :: createfile_qmax
      public :: add_timestep
      public :: createnc1
      public :: createnc2
      public :: createnc3
      public :: write_dynvar

  interface write_dynvar
    module procedure write_dynvar_real
    module procedure write_dynvar_integer2
    module procedure write_dynvar_integer4
  end interface write_dynvar

  integer :: year0, month0, day0, hour0
  integer :: year1, month1, day1, hour1
  character(len=6) , dimension(12) :: xmese
  character(len=10) :: amese
  character(len=4) :: xyear
  character(len=2) :: xmonth
  logical :: bisest
  character(len=32) :: isodate0 , isodate1 , source
  character(len=64) :: timeunit

      contains
!
      subroutine read_config(ifile)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      character(len=*) :: ifile
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      real :: fdum
      character(len=100) :: pname
!
!-----------------------------------------------------------------------
!     Read parameters
!-----------------------------------------------------------------------
!
      call read_rec(ifile, 'NLON', fdum, pname)
      nlon = int(fdum+0.001)
      call read_rec(ifile, 'NLAT', fdum, pname)
      nlat = int(fdum+0.001)
      nbc = nlat
      nlc = nlon
      call read_rec(ifile, 'CONVFAC', fdum, pname)
      convfac = fdum
      call read_rec(ifile, 'ISREAD', fdum, pname)
      isread = int(fdum+0.001)
      call read_rec(ifile, 'ISWRIT', fdum, pname)
      iswrit = int(fdum+0.001)
      call read_rec(ifile, 'TSDATE', fdum, tdate)
      read(tdate(1:4),'(i)') jahr1
      read(tdate(5:6),'(i)') jahr2
      read(tdate(7:8),'(i)') jahr3
      read(tdate(9:10),'(i)') jahr4
      read(tdate(1:10),'(i)') sdate
      call read_rec(ifile, 'TEDATE', fdum, tdate)
      read(tdate(1:10),'(i)') edate
      call read_rec(ifile, 'NDSTEP', fdum, pname)
      dstep = int(fdum+0.001)
      call read_rec(ifile, 'NSTEP', fdum, pname)
      nstep = int(fdum+0.001)
      call read_rec(ifile, 'NDTST', fdum, pname)
      step = int(fdum+0.001)
      call read_rec(ifile, 'NINIRUN', fdum, pname)
      inirun = int(fdum+0.001)
      call read_rec(ifile, 'NINYDAY', fdum, pname)
      yday = int(fdum+0.001)
!
      call read_rec(ifile, 'TDNCAL', fdum, calendario)
      call read_rec(ifile, 'TDNSIM', fdum, sim_name)
      call read_rec(ifile, 'TDNINP', fdum, filemrro)
      call read_rec(ifile, 'TDNRES', fdum, chym_restart)
      call read_rec(ifile, 'TDNINI', fdum, chym_initial)
      call read_rec(ifile, 'TDNOUT', fdum, chym_output)
      call read_rec(ifile, 'TDNQMAX', fdum, chym_qmax)
      call read_rec(ifile, 'TDNSTK', fdum, filestatic)
!
      time = sdate
      if (  trim(calendario) == "gregorian" &
       .or. trim(calendario) == "proleptic" &
       .or. trim(calendario) == "Gregorian" &
       .or. trim(calendario) == "standard"  &
       .or. trim(calendario) == "proleptic_gregorian") then

        ical = 0
      else if (trim(calendario) == "no_leap" .or. trim(calendario) == "365_day") then
        ical = 1
      else if (trim(calendario) == "360_day") then
        ical = 2
      else
        print*,"ERROR please define calendar!"
        call exit(1)
      end if

      end subroutine read_config


      subroutine read_init()
      use netcdf
      implicit none
      integer :: i, ios
      integer :: ihead(8)
      integer :: ncid, varid, ok
      if (.not. allocated(alfa)) allocate(alfa(nlc,nbc))
      if (.not. allocated(fmap)) allocate(fmap(nlc,nbc))
      if (.not. allocated(accl)) allocate(accl(nlc,nbc))
      if (.not. allocated(luse)) allocate(luse(nlc,nbc))
      if (.not. allocated(port)) allocate(port(nlc,nbc))
      if (.not. allocated(port_out)) allocate(port_out(nlc,nbc))
      if (.not. allocated(port_outs)) allocate(port_outs(nlc,nbc))
      if (.not. allocated(port_qmaxs)) allocate(port_qmaxs(nlc,nbc))
      if (.not. allocated(wkm1)) allocate(wkm1(nlc,nbc))
      if (.not. allocated(bwet)) allocate(bwet(nlc,nbc))
      if (.not. allocated(h2o)) allocate(h2o(nlc,nbc))
      if (.not. allocated(lat1)) allocate(lat1(nbc))
      if (.not. allocated(lon1)) allocate(lon1(nlc))
      if (.not. allocated(accl)) allocate(accl(nlc,nbc))
      if (.not. allocated(luse)) allocate(luse(nlc,nbc))
      if (.not. allocated(port)) allocate(port(nlc,nbc))
      if (.not. allocated(wkm1)) allocate(wkm1(nlc,nbc))
      if (.not. allocated(bwet)) allocate(bwet(nlc,nbc))
      if (.not. allocated(h2o)) allocate(h2o(nlc,nbc))
      if (.not. allocated(chym_area)) allocate(chym_area(nlc,nbc))
      if (.not. allocated(chym_drai)) allocate(chym_drai(nlc,nbc))
      if (.not. allocated(chym_dx)) allocate(chym_dx(nlc,nbc))
      if (.not. allocated(chym_lat)) allocate(chym_lat(nlc,nbc))
      if (.not. allocated(chym_lon)) allocate(chym_lon(nlc,nbc))
      if (.not. allocated(manning)) allocate(manning(lntypes))
      port=0; wkm1=0; bwet=0; h2o=0
      port_qmaxs = 0
!
!-----------------------------------------------------------------------
!     Set coefficients
!-----------------------------------------------------------------------
!
      !manning = 0.043
      manning( 1) = 10.0
      manning( 2) = 0.040
      manning( 3) = 0.150
      manning( 4) = 0.100
      manning( 5) = 0.100
      manning( 6) = 0.100
      manning( 7) = 0.650
      manning( 8) = 0.030
      manning( 9) = 0.030
      manning(10) = 0.080
      manning(11) = 0.030
      manning(12) = 0.030
      manning(13) = 0.080
      manning(14) = 0.035
      manning(15) = 0.035
      manning(16) = 0.080
      manning(17) = 0.080
      manning(18) = 0.120
      manning(19) = 0.150
      manning(20) = 0.150
      manning(21) = 0.150
      manning(22) = 0.120
      manning(23) = 0.100
      manning(24) = 0.110
      manning(25) = 0.110
      manning(26) = 0.110
      manning(27) = 0.150
      manning(28) = 0.180
      manning(29) = 0.180
      manning(30) = 0.075
      manning(31) = 0.250
      manning(32) = 0.100
      manning(33) = 0.180
      manning(34) = 0.100
      manning(35) = 0.075
      manning(36) = 0.055
      manning(37) = 0.055
      manning(38) = 0.063
      manning(39) = 0.060
      manning(40) = 0.080
      manning(41) = 0.080
      manning(42) = 0.040
      manning(43) = 0.100
      manning(44) = 0.080
      manning(45) = 0.080
      manning(46) = 0.080
      manning(47) = 0.080
      manning(48) = 0.100
      manning(49) = 0.030
      manning(50) = 0.030
      manning(51) = 0.030
      manning(52) = 0.030
      manning(53) = 0.030
      manning(54) = 0.120
      manning(55) = 0.075
      manning(56) = 0.100
      manning(57) = 0.100
      manning(58) = 0.080
      manning(58) = 0.080
      manning(59) = 0.080
      manning(60) = 0.100
      manning(61) = 0.100
      manning(62) = 0.150
      manning(63) = 0.075
      manning(64) = 0.080
      manning(65) = 0.080
      manning(66) = 0.080
      manning(67) = 0.080
      manning(68) = 0.080
      manning(69) = 0.030
      manning(70) = 0.030
      manning(71) = 0.035
      manning(72) = 0.075
      manning(73) = 0.035
      manning(74) = 0.035
      manning(75) = 0.035
      manning(76) = 0.035
      manning(77) = 0.150
      manning(78) = 0.120
      manning(79) = 0.120
      manning(80) = 0.035
      manning(81) = 0.035
      manning(82) = 0.035
      manning(83) = 0.035
      manning(84) = 0.035
      manning(85) = 0.035
      manning(86) = 0.030
      manning(87) = 0.030
      manning(88) = 0.050
      manning(89) = 0.100
      manning(91) = 0.100
      manning(92) = 0.080
      manning(93) = 0.045
      manning(94) = 0.040
      manning(95) = 0.120
      manning(96) = 0.100
      manning(97:lntypes) = 0.043
!
!-----------------------------------------------------------------------
!     Open, read and close restart file
!-----------------------------------------------------------------------
!
      pstep = 0
      if (isread /= 0) then
      if (myid == 0 ) then
        print*, "read chym restart data"
        call chym_inifile()
      end if
      end if

!
!-----------------------------------------------------------------------
!     Open, read and close static file
!-----------------------------------------------------------------------
!
      if (myid == 0 ) then
      print*,"Inizio lettura statico file:", trim(filestatic)
      call nio_check(nf90_open(trim(filestatic),                        &
         nf90_nowrite, ncid),1)
      call nio_check(nf90_inq_varid(ncid, 'lon', varid),2)
      call nio_check(nf90_get_var(ncid, varid, lon1(:)),3)

      call nio_check(nf90_inq_varid(ncid, 'lat', varid),4)
      call nio_check(nf90_get_var(ncid, varid, lat1(:)),5)

      call nio_check(nf90_inq_varid(ncid, 'fdm', varid),6)
      call nio_check(nf90_get_var(ncid, varid, fmap(:,:)),7)

      call nio_check(nf90_inq_varid(ncid, 'acc', varid),8)
      call nio_check(nf90_get_var(ncid, varid, accl(:,:)),9)

      call nio_check(nf90_inq_varid(ncid, 'lus', varid),10)
      call nio_check(nf90_get_var(ncid, varid, luse(:,:)),11)

      call nio_check(nf90_inq_varid(ncid, 'aer', varid),12)
      call nio_check(nf90_get_var(ncid, varid, chym_area(:,:)),13)

      call nio_check(nf90_inq_varid(ncid, 'dra', varid),14)
      call nio_check(nf90_get_var(ncid, varid, chym_drai(:,:)),15)

      call nio_check(nf90_close(ncid),16)

      end if

      call mpi_bcast(lon1(1),nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(lat1(1),nbc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(fmap(1,1),nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(accl(1,1),nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(luse(1,1),nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(chym_area(1,1),nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(chym_drai(1,1),nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_barrier(mycomm,mpierr)

      do i=1,nlc
        chym_lat(i,:) = lat1
      enddo
      do i=1,nbc
        chym_lon(:,i) = lon1
      enddo
      end subroutine read_init


  subroutine createfile(outname,ncdata)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: outname
    integer, intent(in) :: ncdata
    logical :: fcentury
    integer :: i,j

    fcentury = .true.
    chunksizes(1) = nlon/5
    chunksizes(2) = nlat/5
    chunksizes(3) = 1

!-----------------------------------------------------------------------
!     Create netCDF file
!-----------------------------------------------------------------------
!
    if (.not. allocated(chymout%dimid)) allocate(chymout%dimid(3))
    if (.not. allocated(chymout%varid)) allocate(chymout%varid(4))
!
    call nio_check(nf90_create(trim(outname), nf90_netcdf4,    &
                                chymout%ncid))
!

!
!-----------------------------------------------------------------------
!     Define dimensions
!-----------------------------------------------------------------------
!
    call nio_check(nf90_def_dim(chymout%ncid, 'lon',                  &
                                nlc, chymout%dimid(1)),101)
    call nio_check(nf90_def_dim(chymout%ncid, 'lat',                  &
                                nbc, chymout%dimid(2)),102)
    call nio_check(nf90_def_dim(chymout%ncid, 'time',                 &
                                nf90_unlimited, chymout%dimid(3)),103)
!
!-----------------------------------------------------------------------
!     Define dimension variables
!-----------------------------------------------------------------------
!
    call nio_check(nf90_def_var(chymout%ncid, 'lon', nf90_real,       &
                   chymout%dimid(1), chymout%varid(1)),104)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(1),       &
                   'long_name', 'Longitude'),105)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(1),       &
                   'units', 'degrees_east'),106)
!
    call nio_check(nf90_def_var(chymout%ncid, 'lat', nf90_real,       &
                   chymout%dimid(2), chymout%varid(2)),107)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(2),       &
                   'long_name', 'Latitude'),108)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(2),       &
                   'units', 'degrees_north'),109)
!
      write(amese,'(i10)') edate
      call gmafromindex(ncdata,hour0,day0,month0,year0)
      write(xyear,'(i4)')year0 ; write(xmonth,'(i2)')month0
      bisest = .false.
      if (MOD(year0,400)==0.or.(MOD(year0,4)==0.and.(.not. &
                 (MOD(year0,100)==0)))) then
        xmese = (/'013124','022924','033124','043024','053124','063024', &
                '073124','083124','093024','103124','113024','123124'/)
        bisest = .true.
      else
        xmese = (/'013124','022824','033124','043024','053124','063024', &
                '073124','083124','093024','103124','113024','123124'/)
      end if
      if ((xyear//xmonth=='199912').and.(edate<1000000) ) then
        fcentury = .false.
      end if
      if (xyear(3:4)//xmese(month0) > amese .and. fcentury) then
         call gmafromindex(edate,hour1,day1,month1,year1)
      else
         call gmafromindex(ncdata,hour1,day1,month1,year1)
        read(xmese(month0)(1:2),'(i2)')month1
        read(xmese(month0)(3:4),'(i2)')day1
        read(xmese(month0)(5:6),'(i2)')hour1
      end if

      write(isodate0,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           year0,'-',month0,'-',day0,' ',hour0,':00:00'
      write(isodate1,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           year1,'-',month1,'-',day1,' ',hour1,':00:00'
      write(timeunit,'(a,i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           'seconds since ',year0,'-',month0,'-',day0,' ',hour0, &
          ':00:00 UTC'

    call nio_check(nf90_put_att(chymout%ncid,nf90_global,'date_start',isodate0),111)
    call nio_check(nf90_put_att(chymout%ncid,nf90_global,'date_end',isodate1),112)
    call nio_check(nf90_def_var(chymout%ncid, 'time', nf90_double,       &
                   chymout%dimid(3), chymout%varid(3)),113)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(3),       &
                   'long_name', 'Time'),114)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(3),       &
                   'units', timeunit),115)
    call nio_check(nf90_put_att(chymout%ncid,chymout%varid(3),        &
         'calendar',trim(calendario)),116)
!
!-----------------------------------------------------------------------
!     Define variables
!-----------------------------------------------------------------------
!
      !call nio_check(nf90_def_var(chymout%ncid, 'dis', nf90_int,       &
      !               chymout%dimid, chymout%varid(4)),117)
      call nio_check(nf90_def_var(chymout%ncid, 'dis', nf90_float,      &
                     chymout%dimid, chymout%varid(4)),117)
      call nio_check(nf90_def_var_deflate(chymout%ncid,chymout%varid(4),&
           1,1,1),118)
      call nio_check(nf90_def_var_chunking(chymout%ncid,                &
           chymout%varid(4),0, chunksizes),119)
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'long_name', 'River Discharge'),120)
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'missing_value', 1.0e+36),121)
      !call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
      !               'missing_value', 2147483647),121)
      !call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
      !               'scale_factor', 0.00011641532188114492),122)
      !call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
      !               'add_offset', 250000),123)
!
!-----------------------------------------------------------------------
!     Exit define mode
!-----------------------------------------------------------------------
!
      call nio_check(nf90_enddef(chymout%ncid),124)
!
!-----------------------------------------------------------------------
!     Fill coordinate variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_put_var(chymout%ncid, chymout%varid(1), lon1),125)
      call nio_check(nf90_put_var(chymout%ncid, chymout%varid(2), lat1),126)
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymout%ncid),127)
      iostep = 0
   end subroutine createfile


      subroutine out_init()
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: i
      character(len=100) :: str
      chunksizes(1) = nlon/5
      chunksizes(2) = nlat/5
      chunksizes(3) = 1
!
!-----------------------------------------------------------------------
!     Create netCDF file
!-----------------------------------------------------------------------
!
      if (.not. allocated(chymout%dimid)) allocate(chymout%dimid(3))
      if (.not. allocated(chymout%varid)) allocate(chymout%varid(4))
!
      call nio_check(nf90_create(trim(chym_output), nf90_netcdf4,    &
                                  chymout%ncid))
!
!-----------------------------------------------------------------------
!     Define dimensions
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_dim(chymout%ncid, 'lon',                  &
                                  nlc, chymout%dimid(1)))
      call nio_check(nf90_def_dim(chymout%ncid, 'lat',                  &
                                  nbc, chymout%dimid(2)))
      call nio_check(nf90_def_dim(chymout%ncid, 'time',                 &
                                  nf90_unlimited, chymout%dimid(3)))
!
!-----------------------------------------------------------------------
!     Define dimension variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_var(chymout%ncid, 'lon', nf90_real,       &
                     chymout%dimid(1), chymout%varid(1)))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(1),       &
                     'long_name', 'Longitude'))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(1),       &
                     'units', 'degrees_east'))
!
      call nio_check(nf90_def_var(chymout%ncid, 'lat', nf90_real,       &
                     chymout%dimid(2), chymout%varid(2)))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(2),       &
                     'long_name', 'Latitude'))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(2),       &
                     'units', 'degrees_north'))
!
      call nio_check(nf90_def_var(chymout%ncid, 'time', nf90_int,       &
                     chymout%dimid(3), chymout%varid(3)))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(3),       &
                     'long_name', 'Time'))
      write(str,fmt='("days since ",I4,"-",I2.2,"-",I2.2," 00:00:00")') &
           jahr1, jahr2, jahr3
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(3),       &
                     'units', str))
      call nio_check(nf90_put_att(chymout%ncid,chymout%varid(3),        &
           'calendar','360_day'))
!
!-----------------------------------------------------------------------
!     Define variables
!-----------------------------------------------------------------------
!
      !call nio_check(nf90_def_var(chymout%ncid, 'dis', nf90_int,       &
      !               chymout%dimid, chymout%varid(4)))
      call nio_check(nf90_def_var(chymout%ncid, 'dis', nf90_float,      &
                     chymout%dimid, chymout%varid(4)))
      call nio_check(nf90_def_var_deflate(chymout%ncid,chymout%varid(4),&
           1,1,1))
      call nio_check(nf90_def_var_chunking(chymout%ncid,                &
           chymout%varid(4),0, chunksizes))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'long_name', 'River Discharge'))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'missing_value', 1.0e+36))
      !call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
      !               'missing_value', 2147483647))
      !call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
      !               'scale_factor', 0.00011641532188114492))
      !call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
      !               'add_offset', 250000))
!
!-----------------------------------------------------------------------
!     Exit define mode
!-----------------------------------------------------------------------
!
      call nio_check(nf90_enddef(chymout%ncid))
!
!-----------------------------------------------------------------------
!     Fill coordinate variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_put_var(chymout%ncid, chymout%varid(1), lon1))
      call nio_check(nf90_put_var(chymout%ncid, chymout%varid(2), lat1))
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymout%ncid))
      end subroutine out_init

      subroutine qmax_init()
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: i
      character(len=100) :: str
      chunksizes(1) = nlon/5
      chunksizes(2) = nlat/5
      chunksizes(3) = 1
!
!-----------------------------------------------------------------------
!     Create netCDF file
!-----------------------------------------------------------------------
!
      if (.not. allocated(chymqmax%dimid)) allocate(chymqmax%dimid(3))
      if (.not. allocated(chymqmax%varid)) allocate(chymqmax%varid(4))
!
      call nio_check(nf90_create(trim(chym_qmax), nf90_netcdf4,    &
                                  chymqmax%ncid))
!
!-----------------------------------------------------------------------
!     Define dimensions
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_dim(chymqmax%ncid, 'lon',                  &
                                  nlc, chymqmax%dimid(1)))
      call nio_check(nf90_def_dim(chymqmax%ncid, 'lat',                  &
                                  nbc, chymqmax%dimid(2)))
      call nio_check(nf90_def_dim(chymqmax%ncid, 'time',                 &
                                  nf90_unlimited, chymqmax%dimid(3)))
!
!-----------------------------------------------------------------------
!     Define dimension variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_var(chymqmax%ncid, 'lon', nf90_real,       &
                     chymqmax%dimid(1), chymqmax%varid(1)))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(1),       &
                     'long_name', 'Longitude'))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(1),       &
                     'units', 'degrees_east'))
!
      call nio_check(nf90_def_var(chymqmax%ncid, 'lat', nf90_real,       &
                     chymqmax%dimid(2), chymqmax%varid(2)))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(2),       &
                     'long_name', 'Latitude'))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(2),       &
                     'units', 'degrees_north'))
!
      call nio_check(nf90_def_var(chymqmax%ncid, 'time', nf90_int,       &
                     chymqmax%dimid(3), chymqmax%varid(3)))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(3),       &
                     'long_name', 'Time'))
      write(str,fmt='("years since ",I4,"-",I2.2,"-",I2.2," 00:00:00")') &
           jahr1, jahr2, jahr3
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(3),       &
                     'units', str))
      call nio_check(nf90_put_att(chymqmax%ncid,chymqmax%varid(3),        &
           'calendar','360_day'))
!
!-----------------------------------------------------------------------
!     Define variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_var(chymqmax%ncid, 'qmax', nf90_real,       &
                     chymqmax%dimid, chymqmax%varid(4)))
      call nio_check(nf90_def_var_deflate(chymqmax%ncid,chymqmax%varid(4),&
           1,1,1))
      call nio_check(nf90_def_var_chunking(chymqmax%ncid,                &
           chymqmax%varid(4),0, chunksizes))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(4),       &
                     'long_name', 'River Qmax'))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(4),       &
                     'missing_value', 1.0e20))
!
!-----------------------------------------------------------------------
!     Exit define mode
!-----------------------------------------------------------------------
!
      call nio_check(nf90_enddef(chymqmax%ncid))
!
!-----------------------------------------------------------------------
!     Fill coordinate variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_put_var(chymqmax%ncid, chymqmax%varid(1), lon1))
      call nio_check(nf90_put_var(chymqmax%ncid, chymqmax%varid(2), lat1))
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymqmax%ncid))
!
      end subroutine qmax_init

      subroutine createfile_qmax(restname,ncdata)
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
    character(len=*) , intent(in) :: restname
    integer, intent(in) :: ncdata
    logical :: fcentury
    integer :: i, j
    fcentury = .true.
      chunksizes(1) = nlon/5
      chunksizes(2) = nlat/5
      chunksizes(3) = 1
!
!-----------------------------------------------------------------------
!     Create netCDF file
!-----------------------------------------------------------------------
!
      if (.not. allocated(chymqmax%dimid)) allocate(chymqmax%dimid(3))
      if (.not. allocated(chymqmax%varid)) allocate(chymqmax%varid(5))
!
      call nio_check(nf90_create(trim(restname), nf90_netcdf4,chymqmax%ncid))
!
!-----------------------------------------------------------------------
!     Define dimensions
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_dim(chymqmax%ncid, 'lon', &
                                  nlc, chymqmax%dimid(1)))
      call nio_check(nf90_def_dim(chymqmax%ncid, 'lat', &
                                  nbc, chymqmax%dimid(2)))
      call nio_check(nf90_def_dim(chymqmax%ncid, 'time', &
                                  nf90_unlimited, chymqmax%dimid(3)))
      write(amese,'(i10)') edate
      call gmafromindex(ncdata,hour0,day0,month0,year0)
      write(xyear,'(i4)')year0 ; write(xmonth,'(i2)')month0
      bisest = .false.
      if (MOD(year0,400)==0.or.(MOD(year0,4)==0.and.(.not. &
                 (MOD(year0,100)==0)))) then
        xmese = (/'013124','022924','033124','043024','053124','063024', &
                '073124','083124','093024','103124','113024','123124'/)
        bisest = .true.
      else
        xmese = (/'013124','022824','033124','043024','053124','063024', &
                '073124','083124','093024','103124','113024','123124'/)
      end if
      if ((xyear//xmonth=='199912').and.(edate<1000000) ) then
        fcentury = .false.
      end if
      if (xyear(3:4)//xmese(month0) > amese .and. fcentury) then
         call gmafromindex(edate,hour1,day1,month1,year1)
      else
         call gmafromindex(ncdata,hour1,day1,month1,year1)
        read(xmese(month0)(1:2),'(i2)')month1
        read(xmese(month0)(3:4),'(i2)')day1
        read(xmese(month0)(5:6),'(i2)')hour1
      end if

      write(isodate0,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           year0,'-',month0,'-',day0,' ',hour0,':00:00'
      write(isodate1,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           year1,'-',month1,'-',day1,' ',hour1,':00:00'
      write(timeunit,'(a,i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           'seconds since ',year0,'-',month0,'-',day0,' ',hour0, &
          ':00:00 UTC'
    call nio_check(nf90_put_att(chymqmax%ncid,nf90_global,'date_start',isodate0))
    call nio_check(nf90_put_att(chymqmax%ncid,nf90_global,'date_end',isodate1))

!
!-----------------------------------------------------------------------
!     Define dimension variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_var(chymqmax%ncid, 'lon', nf90_real, &
                     chymqmax%dimid(1), chymqmax%varid(1)))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(1), &
                     'long_name', 'Longitude'))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(1), &
                     'units', 'degrees_east'))
!
      call nio_check(nf90_def_var(chymqmax%ncid, 'lat', nf90_real, &
                     chymqmax%dimid(2), chymqmax%varid(2)))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(2), &
                     'long_name', 'Latitude'))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(2), &
                     'units', 'degrees_north'))
!
    call nio_check(nf90_def_var(chymqmax%ncid, 'time', nf90_double,       &
                   chymqmax%dimid(3), chymqmax%varid(3)))
    call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(3),       &
                   'long_name', 'Time'))
    call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(3),       &
                   'units', timeunit))
    call nio_check(nf90_put_att(chymqmax%ncid,chymqmax%varid(3),        &
         'calendar',trim(calendario)))

!
!-----------------------------------------------------------------------
!     Define variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_var(chymqmax%ncid, 'qmax', nf90_real, &
                     chymqmax%dimid, chymqmax%varid(4)))
      call nio_check(nf90_def_var_deflate(chymqmax%ncid,chymqmax%varid(4),&
           1,1,1))
      call nio_check(nf90_def_var_chunking(chymqmax%ncid,                &
           chymqmax%varid(4),0, chunksizes))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(4), &
                     'long_name', 'River Qmax'))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(4), &
                     'missing_value', 1.0e20))
!
!-----------------------------------------------------------------------
!     Exit define mode
!-----------------------------------------------------------------------
!
      call nio_check(nf90_enddef(chymqmax%ncid))
!
!-----------------------------------------------------------------------
!     Fill coordinate variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_put_var(chymqmax%ncid, chymqmax%varid(1), lon1))
      call nio_check(nf90_put_var(chymqmax%ncid, chymqmax%varid(2), lat1))
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymqmax%ncid))
!
      iqstep = 0
      end subroutine createfile_qmax

      subroutine createfile_rst(restname,ncdata)
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
    character(len=*) , intent(in) :: restname
    integer, intent(in) :: ncdata
    logical :: fcentury
    integer :: i, j
    fcentury = .true.
      chunksizes(1) = nlon/5
      chunksizes(2) = nlat/5
      chunksizes(3) = 1
!
!-----------------------------------------------------------------------
!     Create netCDF file
!-----------------------------------------------------------------------
!
      if (.not. allocated(chymrst%dimid)) allocate(chymrst%dimid(3))
      if (.not. allocated(chymrst%varid)) allocate(chymrst%varid(5))
!
      call nio_check(nf90_create(trim(restname), nf90_netcdf4,chymrst%ncid))
!
!-----------------------------------------------------------------------
!     Define dimensions
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_dim(chymrst%ncid, 'lon', &
                                  nlc, chymrst%dimid(1)))
      call nio_check(nf90_def_dim(chymrst%ncid, 'lat', &
                                  nbc, chymrst%dimid(2)))
      call nio_check(nf90_def_dim(chymrst%ncid, 'time', &
                                  nf90_unlimited, chymrst%dimid(3)))
      write(amese,'(i10)') edate
      call gmafromindex(ncdata,hour0,day0,month0,year0)
      write(xyear,'(i4)')year0 ; write(xmonth,'(i2)')month0
      bisest = .false.
      if (MOD(year0,400)==0.or.(MOD(year0,4)==0.and.(.not. &
                 (MOD(year0,100)==0)))) then
        xmese = (/'013124','022924','033124','043024','053124','063024', &
                '073124','083124','093024','103124','113024','123124'/)
        bisest = .true.
      else
        xmese = (/'013124','022824','033124','043024','053124','063024', &
                '073124','083124','093024','103124','113024','123124'/)
      end if
      if ((xyear//xmonth=='199912').and.(edate<1000000) ) then
        fcentury = .false.
      end if
      if (xyear(3:4)//xmese(month0) > amese .and. fcentury) then
         call gmafromindex(edate,hour1,day1,month1,year1)
      else
         call gmafromindex(ncdata,hour1,day1,month1,year1)
        read(xmese(month0)(1:2),'(i2)')month1
        read(xmese(month0)(3:4),'(i2)')day1
        read(xmese(month0)(5:6),'(i2)')hour1
      end if

      write(isodate0,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           year0,'-',month0,'-',day0,' ',hour0,':00:00'
      write(isodate1,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           year1,'-',month1,'-',day1,' ',hour1,':00:00'
      write(timeunit,'(a,i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           'seconds since ',year0,'-',month0,'-',day0,' ',hour0, &
          ':00:00 UTC'
    call nio_check(nf90_put_att(chymrst%ncid,nf90_global,'date_start',isodate0))
    call nio_check(nf90_put_att(chymrst%ncid,nf90_global,'date_end',isodate1))

!
!-----------------------------------------------------------------------
!     Define dimension variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_var(chymrst%ncid, 'lon', nf90_real, &
                     chymrst%dimid(1), chymrst%varid(1)))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(1), &
                     'long_name', 'Longitude'))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(1), &
                     'units', 'degrees_east'))
!
      call nio_check(nf90_def_var(chymrst%ncid, 'lat', nf90_real, &
                     chymrst%dimid(2), chymrst%varid(2)))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(2), &
                     'long_name', 'Latitude'))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(2), &
                     'units', 'degrees_north'))
!
    call nio_check(nf90_def_var(chymrst%ncid, 'time', nf90_double,       &
                   chymrst%dimid(3), chymrst%varid(3)))
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(3),       &
                   'long_name', 'Time'))
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(3),       &
                   'units', timeunit))
    call nio_check(nf90_put_att(chymrst%ncid,chymrst%varid(3),        &
         'calendar',trim(calendario)))

!
!-----------------------------------------------------------------------
!     Define variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_var(chymrst%ncid, 'dis', nf90_real, &
                     chymrst%dimid, chymrst%varid(4)))
      call nio_check(nf90_def_var_deflate(chymrst%ncid,chymrst%varid(4),&
           1,1,1))
      call nio_check(nf90_def_var_chunking(chymrst%ncid,                &
           chymrst%varid(4),0, chunksizes))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(4), &
                     'long_name', 'River Discharge'))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(4), &
                     'missing_value', 1.0e20))
      call nio_check(nf90_def_var(chymrst%ncid, 'h2o', nf90_real, &
                     chymrst%dimid, chymrst%varid(5)))
      call nio_check(nf90_def_var_deflate(chymrst%ncid,chymrst%varid(5),&
           1,1,1))
      call nio_check(nf90_def_var_chunking(chymrst%ncid,                &
           chymrst%varid(5),0, chunksizes))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(5), &
                     'long_name', 'Total water'))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(5), &
                     'missing_value', 1.0e20))
!
!-----------------------------------------------------------------------
!     Exit define mode
!-----------------------------------------------------------------------
!
      call nio_check(nf90_enddef(chymrst%ncid))
!
!-----------------------------------------------------------------------
!     Fill coordinate variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(1), lon1))
      call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(2), lat1))
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymrst%ncid))
!
      irstep = 0
      end subroutine createfile_rst
  subroutine createnc1(hourstep)
    implicit none
    integer ora,giorno,mese,anno,yy,mm
    integer :: oldmese             ! = 01
    integer , intent(in) :: hourstep
    character(len=10) :: cfl
    character(len=256) :: mess
    logical :: first
    data first /.true./
    save first
    save mm,yy,oldmese

    call gmafromindex(time,ora,giorno,mese,anno)
    if ( first ) then
       read(tsdate(5:6),'(i2)') mm
       read(tsdate(1:4),'(i4)') yy
       oldmese = mese
       first = .false.
    end if
    if (mese /= oldmese .and. hourstep > 1) then
      mm = mm + 1
      if (mm >= 13) then
         mm =  1
           yy = yy + 1
      end if
      if (mm < 10) then
         write(cfl,'(i4,a,i1,a)') yy,'0',mm,'0100'
      else
         write(cfl,'(i4,i2,a)') yy,mm,'0100'
      end if
      call closefile(chymout%ncid)
      write (mess,'(15x,a,a)') 'Creating new output file:  '// &
         trim(sim_name),'_'//trim(cfl)//'.nc'
      write (6,'(a)') mess(1:len_trim(mess))
      call createfile(trim(sim_name)//'_'//trim(cfl)//'.nc',time)
    end if
    oldmese = mese
    return
  end subroutine createnc1

  subroutine createnc2(hourstep)
    implicit none
    integer ora,giorno,mese,anno,yy,mm
    integer :: oldanno             ! = 01
    integer , intent(in) :: hourstep
    character(len=10) :: cfl
    character(len=256) :: mess
    logical :: first
    data first /.true./
    save first
    save mm,yy,oldanno

    call gmafromindex(time,ora,giorno,mese,anno)
    if ( first ) then
       read(tsdate(5:6),'(i2)') mm
       read(tsdate(1:4),'(i4)') yy
       oldanno = anno
       first = .false.
    end if

    if ( anno /= oldanno .and. hourstep > 1 ) then
      yy = yy + 1
      write(cfl,'(i4,a)') yy,'010100'
      write (mess,'(15x,a,a)') 'Creating new Restart file:  '// &
         trim(sim_name),'_'//trim(cfl)//'_rst.nc'
      write (6,'(a)') mess(1:len_trim(mess))
      call createfile_rst(trim(sim_name)//'_'//trim(cfl)//'_rst.nc',time)
      call add_timestep(chymrst%ncid,chymrst%varid(3),irstep)
      call write_dynvar(chymrst%ncid,chymrst%varid(4),port,irstep)
      call write_dynvar(chymrst%ncid,chymrst%varid(5),h2o,irstep)
      call closefile(chymrst%ncid)
    end if
    oldanno = anno
    return
  end subroutine createnc2

  subroutine createnc3(hourstep)
    implicit none
    integer ora,giorno,mese,anno,yy,mm
    integer :: oldanno         ! = 01
    integer , intent(in) :: hourstep
    character(len=10) :: cfl
    character(len=256) :: mess
    logical :: first
    data first /.true./
    save first
    save mm,yy,oldanno

    call gmafromindex(time,ora,giorno,mese,anno)
    if ( first ) then
       read(tsdate(5:6),'(i2)') mm
       read(tsdate(1:4),'(i4)') yy
       oldanno = anno
       first = .false.
    end if
    if (anno /= oldanno .and. hourstep > 1) then
      yy = yy + 1
      write(cfl,'(i4,a)') yy-1,'000100'
      write (mess,'(15x,a,a)') 'Creating new Qmax file:  '// &
         trim(sim_name),'_'//trim(cfl)//'_qmax.nc'
      write (6,'(a)') mess(1:len_trim(mess))
      call createfile_qmax(trim(sim_name)//'_'//trim(cfl)//'_qmax.nc',time)
      call add_timestep(chymqmax%ncid,chymqmax%varid(3),iqstep)
      call write_dynvar(chymqmax%ncid,chymqmax%varid(4),port_qmaxs,iqstep)
      call closefile(chymqmax%ncid)
      port_qmaxs = 0
    end if
    oldanno = anno
    return
  end subroutine createnc3



      subroutine chym_rst(istep)
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      integer, intent(in) :: istep
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: start(4), count(4), len
!
!-----------------------------------------------------------------------
!     Write data to file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_inquire_dimension(chymrst%ncid,               &
                     chymrst%dimid(3), len=len))
      start = (/ len+1, 1, 1, 1 /)
      count = (/ 1, 1, 1, 1 /)
      call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(3),       &
                    (/ istep-1 /), start, count))
!
      start = (/ 1, 1, len+1, 1 /)
      count = (/ nlc, nbc, 1, 1 /)
      call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(4),       &
                       port, start, count))
!
      start = (/ 1, 1, len+1, 1 /)
      count = (/ nlc, nbc, 1, 1 /)
      call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(5),       &
                       h2o, start, count))
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymrst%ncid))
!
!-----------------------------------------------------------------------
!     Close file
!-----------------------------------------------------------------------
!
      if (istep == nstep) call nio_check(nf90_close(chymrst%ncid))
!
      end subroutine chym_rst
      subroutine chym_inifile()
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: ncid, varid, start(4), count(4)
!
!-----------------------------------------------------------------------
!     Open netCDF file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_open(trim(chym_initial), nf90_nowrite, ncid))
!
!-----------------------------------------------------------------------
!     Read variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_inq_varid(ncid, 'dis', varid))
      start = (/ 1, 1, 1, 1 /)
      count = (/ nlc, nbc, 1, 1 /)
      call nio_check(nf90_get_var(ncid, varid, port,                    &
                       start=start, count=count))
!
      call nio_check(nf90_inq_varid(ncid, 'h2o', varid))
      start = (/ 1, 1, 1, 1 /)
      count = (/ nlc, nbc, 1, 1 /)
      call nio_check(nf90_get_var(ncid, varid, h2o,                     &
                       start=start, count=count))
!
      call nio_check(nf90_inq_varid(ncid, 'time', varid))
      call nio_check(nf90_get_var(ncid, varid, pstep))
      pstep = pstep+1
!
!-----------------------------------------------------------------------
!     Close file
!-----------------------------------------------------------------------
!

      end subroutine chym_inifile


      subroutine read_runoff(step)
      use netcdf
      integer, intent(in) :: step
      integer, save :: ncid, varid, oldyear=-2
      logical, save :: first = .true.
      integer, save :: stepminus = 0
      integer :: stepr
      character(len=256) filerunoff
      if (myid == 0) then
      stepr = step - stepminus
        if (first == .false. .and. year == oldyear + 1) then
           call nio_check(nf90_close(ncid))
           first = .true.
           stepminus = step-1
           stepr = step - stepminus
        end if
        print*,"We are going to read step : ",stepr," of chym_runoff"
        if (first == .true.) then
          oldyear=year
          write(filerunoff,'(a,i4,a)') trim(filemrro)//'_',oldyear, &
              '.nc'
          print*,"we are opening file: ", trim(filerunoff)
          call nio_check(nf90_open(trim(filerunoff),                     &
             nf90_nowrite, ncid))
          if ( nf90_inq_varid(ncid, 'mrro', varid) /= nf90_noerr ) then
            call nio_check(nf90_inq_varid(ncid, 'ro',varid))
          end if
          first = .false.
        end if
      call nio_check(nf90_get_var(ncid, varid, chym_runoff(:,:),        &
           start = (/ 1, 1, stepr /), count = (/nlc , nbc, 1 /)))

      end if

      call mpi_bcast(chym_runoff(1,1),nbc*nlc,MPI_REAL,0,mycomm,mpierr)

      end subroutine read_runoff

      subroutine chym_out(istep)
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      integer, intent(in) :: istep
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: len, start(4), count(4)
!
!-----------------------------------------------------------------------
!     Write data to file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_inquire_dimension(chymout%ncid,               &
                     chymout%dimid(3), len=len))
!
      start = (/ len+1, 1, 1, 1 /)
      count = (/ 1, 1, 1, 1 /)
      call nio_check(nf90_put_var(chymout%ncid, chymout%varid(3),       &
                    (/ istep-1 /), start, count))
!
      start = (/ 1, 1, len+1, 1 /)
      count = (/ nlc, nbc, 1, 1 /)
      call nio_check(nf90_put_var(chymout%ncid, chymout%varid(4),       &
                     port_out, start, count))
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymout%ncid))
!
!-----------------------------------------------------------------------
!     Close file
!-----------------------------------------------------------------------
!
      if (istep == nstep) call nio_check(nf90_close(chymout%ncid))
!
      end subroutine chym_out

      subroutine chym_wqmax(istep)
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      integer, intent(in) :: istep
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: len, start(4), count(4)
!
!-----------------------------------------------------------------------
!     Write data to file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_inquire_dimension(chymqmax%ncid,               &
                     chymqmax%dimid(3), len=len))
!
      start = (/ len+1, 1, 1, 1 /)
      count = (/ 1, 1, 1, 1 /)
      call nio_check(nf90_put_var(chymqmax%ncid, chymqmax%varid(3),       &
                    (/ istep-1 /), start, count))
!
      start = (/ 1, 1, len+1, 1 /)
      count = (/ nlc, nbc, 1, 1 /)
      call nio_check(nf90_put_var(chymqmax%ncid, chymqmax%varid(4),       &
                     port_qmaxs, start, count))
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymqmax%ncid))
!
!-----------------------------------------------------------------------
!     Close file
!-----------------------------------------------------------------------
!
      if (istep == nstep) call nio_check(nf90_close(chymqmax%ncid))
!
      end subroutine chym_wqmax

!

      subroutine read_rec(ifile, key, value, pname)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      character(len=*), intent(in) :: ifile, key
      character(len=*), intent(inout) :: pname
      real, intent(inout) :: value
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: ios
!
!-----------------------------------------------------------------------
!     Open configuration file
!-----------------------------------------------------------------------
!
      open(lu, file=trim(ifile), access='sequential', form='formatted', &
           status='old', iostat=ios)
      if (ios /= 0) then
         write(*,*) '[error] -- file '//trim(ifile)//' not found!'
         stop
      endif
!
!-----------------------------------------------------------------------
!     Read and find the parameter
!-----------------------------------------------------------------------
!
      ios = 0
      do while (ios == 0)
        read(lu, fmt='(A90)', iostat=ios) pname
        if (index(pname, key) /= 0) then
          if (key(1:1) == 't' .or. key(1:1) == 'T') then
            read(lu, '(A90)') pname
            if (myid == 0) then
            write(*,*) 'Parameter: '//trim(key)//' = ', trim(pname)
            endif
          else
            read(lu, *) value
            if (myid == 0) then
            write(*,*) 'Parameter: '//trim(key)//' = ', value
            endif
          end if
          exit
        end if
      end do
!
!-----------------------------------------------------------------------
!     Close configuration file
!-----------------------------------------------------------------------
!
      close(lu, status='keep', iostat=ios)
      if (ios /= 0) then
        write(*,*) '[error] -- '//trim(ifile)//'is not closed!'
      endif
      end subroutine read_rec

      subroutine nio_check(status,nerr)
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      integer, intent(in) :: status
      integer,optional, intent(in) :: nerr

      if (status /= nf90_noerr) then
        print*, trim(nf90_strerror(status))
        if (present(nerr) == .True.) then
        print*, "ERRORRRR in NETCDF, err Number:",nerr
        else
        print*, "ERRORRRR in NETCDF, err Number: undefined"
        end if
        stop 2
      end if
!
      end subroutine nio_check
  subroutine closefile(ncid)
    use netcdf
    implicit none
    integer, intent(in) :: ncid
    call nio_check(nf90_close(ncid))
  end subroutine closefile
  subroutine add_timestep(ncid,time_varid,istep)
    use netcdf
    implicit none
    double precision , dimension(1) :: xtime
    integer, intent(in) :: ncid,time_varid
    integer, intent(inout) :: istep
    istep = istep + 1
    xtime(1) = dble(istep-1)*3600.0D0*dstep
    istart(1) = istep
    icount(1) = 1
    call nio_check(nf90_put_var(ncid,time_varid,xtime,istart(1:1), &
           icount(1:1)))
  end subroutine add_timestep
  subroutine write_dynvar_real(ncid,vname_id,rval,iistep)
    use netcdf
    implicit none
    integer , intent(in) :: vname_id,ncid,iistep
    real , dimension(:,:) , intent(in) :: rval
    istart(1) = 1
    istart(2) = 1
    istart(3) = iistep
    icount(1) = nlc
    icount(2) = nbc
    icount(3) = 1
    call nio_check(nf90_put_var(ncid,vname_id,rval,istart,icount))
  end subroutine write_dynvar_real

  subroutine write_dynvar_integer2(ncid,vname_id,ival,iistep)
    use netcdf
    implicit none
    integer , intent(in) :: vname_id,ncid,iistep
    integer(2) , dimension(:,:) , intent(in) :: ival
    istart(1) = 1
    istart(2) = 1
    istart(3) = iistep
    icount(1) = nlc
    icount(2) = nbc
    icount(3) = 1
    call nio_check(nf90_put_var(ncid,vname_id,ival,istart,icount))
  end subroutine write_dynvar_integer2

  subroutine write_dynvar_integer4(ncid,vname_id,ival,iistep)
    use netcdf
    implicit none
    integer , intent(in) :: vname_id,ncid,iistep
    integer(4) , dimension(:,:) , intent(in) :: ival
    integer :: indx
    istart(1) = 1
    istart(2) = 1
    istart(3) = iistep
    icount(1) = nlc
    icount(2) = nbc
    icount(3) = 1
    call nio_check(nf90_put_var(ncid,vname_id,ival,istart,icount))
  end subroutine write_dynvar_integer4


      end module mod_io
