
module mod_io

  use, intrinsic :: iso_fortran_env
  use mod_param
  use mod_mpimess
  use mod_varandtypes
  use mod_time

  implicit none
  private

  integer , dimension(4) , save :: icount , istart

  public :: read_config
  public :: read_init
  public :: read_runoff
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
  character(len=32) :: isodate0 , isodate1
  character(len=64) :: timeunit

  contains

    subroutine read_config(ifile)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      character(len=*) , intent(in) :: ifile
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: iretval
      integer :: lun

      namelist /iniparam/ convfac, thrriv, irloss, irmonfac, efficiency
      namelist /inputparam/ isread, inirun, yday, dstep, step, &
                            tdninp, tdnini, tdnstk
      namelist /outparam/ iorstfreq, tdnsim
      namelist /timeparam/ sdate, edate, calendar
!
!-----------------------------------------------------------------------
!     Read parameters
!-----------------------------------------------------------------------
!
      open(newunit=lun, file=ifile, status='old', &
           action='read', iostat=iretval)
      if ( iretval /= 0 ) then
        write(error_unit,*) 'Error opening namelist file'//trim(ifile)
        stop
      end if
      read(lun, nml=iniparam, iostat=iretval)
      if ( iretval /= 0 ) then
        write(error_unit,*) 'Error reading iniparam namelist'
        stop
      end if
      rewind(lun)
      read(lun, nml=inputparam, iostat=iretval)
      if ( iretval /= 0 ) then
        write(error_unit,*) 'Error reading inputparam namelist'
        stop
      end if
      rewind(lun)
      iorstfreq = 2 ! monthly
      read(lun, nml=outparam, iostat=iretval)
      if ( iretval /= 0 ) then
        write(error_unit,*) 'Error reading outparam namelist'
        stop
      end if
      rewind(lun)
      read(lun, nml=timeparam, iostat=iretval)
      if ( iretval /= 0 ) then
        write(error_unit,*) 'Error reading timeparam namelist'
        stop
      end if
      jahr1 = sdate/1000000
      jahr2 = (sdate-jahr1*1000000)/10000
      jahr3 = (sdate-(jahr1*1000000+jahr2*10000))/100
      jahr4 = sdate-jahr1*1000000+jahr2*10000+jahr3*100

      time = sdate
      if (  trim(calendar) == "gregorian" &
       .or. trim(calendar) == "proleptic" &
       .or. trim(calendar) == "Gregorian" &
       .or. trim(calendar) == "standard"  &
       .or. trim(calendar) == "proleptic_gregorian") then

        ical = 0
      else if (trim(calendar) == "no_leap" .or. trim(calendar) == "365_day") then
        ical = 1
      else if (trim(calendar) == "360_day") then
        ical = 2
      else
        print*,"ERROR please define calendar!"
        call exit(1)
      end if

      if ( myid == 0 ) then
        print *, 'Config read in.'
      end if
    end subroutine read_config

    subroutine read_init()
      use netcdf
      implicit none
      integer :: ncid, dimid, varid
#ifdef RUNOFF
#ifdef NILE
      integer :: ii, jj
#endif
      integer :: i, j, ilnd, idir
#endif

!
!-----------------------------------------------------------------------
!     Open, read and close static file
!-----------------------------------------------------------------------
!
      if (myid == 0 ) then
        print*,"Reading Static Data file:", trim(tdnstk)
        call nio_check(nf90_open(trim(tdnstk),                        &
           nf90_nowrite, ncid),1)

        call nio_check(nf90_inq_dimid(ncid,'lon',dimid),1)
        call nio_check(nf90_inquire_dimension(ncid,dimid,len=nlon),1)
        call nio_check(nf90_inq_dimid(ncid,'lat',dimid),1)
        call nio_check(nf90_inquire_dimension(ncid,dimid,len=nlat),1)
        print *, "DIMENSIONS : ", nlon,  ' x ', nlat
      end if

      call mpi_bcast(nlon,1,MPI_INT,0,mycomm,mpierr)
      call mpi_bcast(nlat,1,MPI_INT,0,mycomm,mpierr)

      nbc = nlat
      nlc = nlon
      if (myid == 0 ) then
        print *, 'Allocating space.'
      end if
      if (.not. allocated(alfa)) allocate(alfa(nlc,nbc))
      if (.not. allocated(fmap)) allocate(fmap(nlc,nbc))
      if (.not. allocated(accl)) allocate(accl(nlc,nbc))
      if (.not. allocated(luse)) allocate(luse(nlc,nbc))
      if (.not. allocated(farm)) allocate(farm(nlc,nbc))
      if (.not. allocated(port)) allocate(port(nlc,nbc))
      if (.not. allocated(port_out)) allocate(port_out(nlc,nbc))
#ifdef RUNOFF
      if (.not. allocated(roff_out)) allocate(roff_out(nlc,nbc))
#endif
      if (.not. allocated(port_qmaxs)) allocate(port_qmaxs(nlc,nbc))
      if (.not. allocated(wkm1)) allocate(wkm1(nlc,nbc))
      if (.not. allocated(bwet)) allocate(bwet(nlc,nbc))
      if (.not. allocated(h2o)) allocate(h2o(nlc,nbc))
      if (.not. allocated(wkm1)) allocate(wkm1(nlc,nbc))
      if (.not. allocated(bwet)) allocate(bwet(nlc,nbc))
      if (.not. allocated(h2o)) allocate(h2o(nlc,nbc))
      if (.not. allocated(chym_area)) allocate(chym_area(nlc,nbc))
      if (.not. allocated(chym_drai)) allocate(chym_drai(nlc,nbc))
      if (.not. allocated(chym_lsm)) allocate(chym_lsm(nlc,nbc))
      if (.not. allocated(chym_dx)) allocate(chym_dx(nlc,nbc))
      if (.not. allocated(chym_lat)) allocate(chym_lat(nlc,nbc))
      if (.not. allocated(chym_lon)) allocate(chym_lon(nlc,nbc))
      if (.not. allocated(corner_lon)) allocate(corner_lon(4,nlc,nbc))
      if (.not. allocated(corner_lat)) allocate(corner_lat(4,nlc,nbc))
      if (.not. allocated(manning)) allocate(manning(lntypes))

      port = 0
      wkm1 = 0
      bwet = 0
      h2o = 0
      port_qmaxs = 0

      if (myid == 0 ) then

        print *, 'Reading data from static file...'

        call nio_check(nf90_inq_varid(ncid, 'manning', varid),1)
        call nio_check(nf90_get_var(ncid, varid, manning),1)

        call nio_check(nf90_inq_varid(ncid, 'lon', varid),2)
        call nio_check(nf90_get_var(ncid, varid, chym_lon),3)

        call nio_check(nf90_inq_varid(ncid, 'corner_lon', varid),2)
        call nio_check(nf90_get_var(ncid, varid, corner_lon),3)

        call nio_check(nf90_inq_varid(ncid, 'lat', varid),4)
        call nio_check(nf90_get_var(ncid, varid, chym_lat),5)

        call nio_check(nf90_inq_varid(ncid, 'corner_lat', varid),4)
        call nio_check(nf90_get_var(ncid, varid, corner_lat),5)

        call nio_check(nf90_inq_varid(ncid, 'fdm', varid),6)
        call nio_check(nf90_get_var(ncid, varid, fmap(:,:)),7)

        call nio_check(nf90_inq_varid(ncid, 'acc', varid),8)
        call nio_check(nf90_get_var(ncid, varid, accl(:,:)),9)

        call nio_check(nf90_inq_varid(ncid, 'lus', varid),10)
        call nio_check(nf90_get_var(ncid, varid, luse(:,:)),11)

        call nio_check(nf90_inq_varid(ncid, 'alf', varid),12)
        call nio_check(nf90_get_var(ncid, varid, alfa(:,:)),13)

        call nio_check(nf90_inq_varid(ncid, 'aer', varid),14)
        call nio_check(nf90_get_var(ncid, varid, chym_area(:,:)),15)

        call nio_check(nf90_inq_varid(ncid, 'dra', varid),16)
        call nio_check(nf90_get_var(ncid, varid, chym_drai(:,:)),17)

        call nio_check(nf90_close(ncid),18)

        print *, 'Done!'

      end if
!
!-----------------------------------------------------------------------
!     Open, read and close restart file
!-----------------------------------------------------------------------
!
      call mpi_bcast(manning,lntypes,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(chym_lon,nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(chym_lat,nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(corner_lat,4*nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(corner_lon,4*nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(fmap,nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(accl,nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(alfa,nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(luse,nbc*nlc,MPI_INT, 0,mycomm,mpierr)
      call mpi_bcast(chym_area,nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(chym_drai,nbc*nlc,MPI_REAL, 0,mycomm,mpierr)
      call mpi_barrier(mycomm,mpierr)

#ifdef RUNOFF
      chym_lsm = 0.0
      do j = 2 , nbc-1
        do i = 2, nlc-1
          idir = fmap(i,j)
          if ( idir >= 1 .and. idir <= 8 ) then
            ilnd = luse(i+ir(idir),j+jr(idir))
            if ( chym_drai(i,j) > thrriv .and. ilnd == ocean ) then
              chym_lsm(i,j) = 1.0
            end if
          end if
#ifdef NILE
          if ( is_inbox(lat_damietta,lon_damietta, &
                corner_lat(:,i,j),corner_lon(:,i,j)) ) then
            call find_nearest_land(i,j,ii,jj)
            chym_lsm(ii,jj) = 1.0
          end if
          if ( is_inbox(lat_rosetta,lon_rosetta, &
                corner_lat(:,i,j),corner_lon(:,i,j)) ) then
            call find_nearest_land(i,j,ii,jj)
            chym_lsm(ii,jj) = 1.0
          end if
#endif
        end do
      end do
#endif

      do j = 2 , nbc-1
        do i = 2, nlc-1
          idir = fmap(i,j)
          chym_dx(i,j) = geodistance(chym_lat(i,j),chym_lon(i,j), &
                 chym_lat(i+ir(idir),j+jr(idir)),                 &
                 chym_lon(i+ir(idir),j+jr(idir)))
          if ( luse(i,j) == 30 .or. luse(i,j) == 31 .or. &
               luse(i,j) == 35 .or. luse(i,j) == 36 .or. &
               luse(i,j) == 37 .or. luse(i,j) == 38 .or. &
               luse(i,j) == 39 .or. luse(i,j) == 76 .or. &
               luse(i,j) == 92 .or. luse(i,j) == 93 .or. &
               luse(i,j) == 94 .or. luse(i,j) == 95 .or. &
               luse(i,j) == 96 ) then
            farm(i,j) = .true.
          else
            farm(i,j) = .false.
          end if
        end do
      end do

      if (isread /= 0) then
        if (myid == 0 ) then
          print*, "read chym restart data"
          call chym_inifile()
        end if
      end if
    end subroutine read_init


  subroutine createfile(outname,ncdata)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: outname
    integer, intent(in) :: ncdata
    logical :: fcentury

    fcentury = .true.
    chunksizes(1) = nlon/5
    chunksizes(2) = nlat/5
    chunksizes(3) = 1

!-----------------------------------------------------------------------
!     Create netCDF file
!-----------------------------------------------------------------------
!
    if (.not. allocated(chymout%dimid)) allocate(chymout%dimid(3))
#ifdef RUNOFF
    if (.not. allocated(chymout%varid)) allocate(chymout%varid(5))
#else
    if (.not. allocated(chymout%varid)) allocate(chymout%varid(4))
#endif
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
                   chymout%dimid(1:2), chymout%varid(1)),104)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(1),       &
                   'long_name', 'Longitude'),105)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(1),       &
                   'units', 'degrees_east'),106)
!
    call nio_check(nf90_def_var(chymout%ncid, 'lat', nf90_real,       &
                   chymout%dimid(1:2), chymout%varid(2)),107)
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
         'calendar',trim(calendar)),116)
!
!-----------------------------------------------------------------------
!     Define variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_var(chymout%ncid, 'dis', nf90_float,      &
                     chymout%dimid, chymout%varid(4)),117)
      call nio_check(nf90_def_var_deflate(chymout%ncid,chymout%varid(4),&
           1,1,1),118)
      call nio_check(nf90_def_var_chunking(chymout%ncid,                &
           chymout%varid(4),0, chunksizes),119)
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'long_name', 'River Discharge'),120)
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'units', 'm3 s-1'),120)
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'missing_value', 1.0e+36),121)
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'coordinates', "lat lon"))
#ifdef RUNOFF
      call nio_check(nf90_def_var(chymout%ncid, 'runoff', nf90_float,      &
                     chymout%dimid, chymout%varid(5)),117)
      call nio_check(nf90_def_var_deflate(chymout%ncid,chymout%varid(5),&
           1,1,1),118)
      call nio_check(nf90_def_var_chunking(chymout%ncid,                &
           chymout%varid(5),0, chunksizes),119)
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(5),       &
                     'units', 'm s-1'),120)
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(5),       &
                     'long_name', 'River runoff'),120)
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(5),       &
                     'missing_value', 1.0e+36),121)
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(5),       &
                     'coordinates', "lat lon"))
#endif
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
      call nio_check(nf90_put_var(chymout%ncid, chymout%varid(1), chym_lon),125)
      call nio_check(nf90_put_var(chymout%ncid, chymout%varid(2), chym_lat),126)
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymout%ncid),127)
      iostep = 0
   end subroutine createfile

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
                     chymqmax%dimid(1:2), chymqmax%varid(1)))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(1), &
                     'long_name', 'Longitude'))
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(1), &
                     'units', 'degrees_east'))
!
      call nio_check(nf90_def_var(chymqmax%ncid, 'lat', nf90_real, &
                     chymqmax%dimid(1:2), chymqmax%varid(2)))
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
         'calendar',trim(calendar)))

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
      call nio_check(nf90_put_att(chymqmax%ncid, chymqmax%varid(4), &
                     'coordinates', "lat lon"))
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
      call nio_check(nf90_put_var(chymqmax%ncid, chymqmax%varid(1), chym_lon))
      call nio_check(nf90_put_var(chymqmax%ncid, chymqmax%varid(2), chym_lat))
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
                     chymrst%dimid(1:2), chymrst%varid(1)))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(1), &
                     'long_name', 'Longitude'))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(1), &
                     'units', 'degrees_east'))
!
      call nio_check(nf90_def_var(chymrst%ncid, 'lat', nf90_real, &
                     chymrst%dimid(1:2), chymrst%varid(2)))
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
         'calendar',trim(calendar)))

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
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(4),       &
                     'coordinates', "lat lon"))
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
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(5),       &
                     'coordinates', "lat lon"))
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
      call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(1), chym_lon))
      call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(2), chym_lat))
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
       mm = mese
       yy = anno
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
         trim(tdnsim),'_'//trim(cfl)//'.nc'
      write (6,'(a)') mess(1:len_trim(mess))
      call createfile(trim(tdnsim)//'_'//trim(cfl)//'.nc',time)
    end if
    oldmese = mese
    return
  end subroutine createnc1

  subroutine createnc2(hourstep,itype)
    implicit none
    integer , intent(in) :: hourstep , itype
    integer ora,giorno,mese,anno,yy,mm
    integer :: old
    character(len=10) :: cfl
    character(len=256) :: mess
    logical :: first , donow
    data first /.true./
    save first
    save mm,yy,old

    call gmafromindex(time,ora,giorno,mese,anno)
    if ( first ) then
       mm = mese
       yy = anno
       if ( itype == 1 ) then
         old = anno
       else
         old = mese
       end if
       first = .false.
    end if

    if ( itype == 1 ) then
      donow = anno /= old
    else
      donow = mese /= old
    end if
    if ( donow .and. hourstep > 1 ) then
      mm = mm + 1
      if ( mm == 13 ) then
        yy = yy + 1
        mm = 1
      end if
      write(cfl,'(i0.4,i0.2,a)') yy,mm,'0100'
      write (mess,'(15x,a,a)') 'Creating new Restart file:  '// &
         trim(tdnsim),'_'//trim(cfl)//'_rst.nc'
      write (6,'(a)') mess(1:len_trim(mess))
      call createfile_rst(trim(tdnsim)//'_'//trim(cfl)//'_rst.nc',time)
      call add_timestep(chymrst%ncid,chymrst%varid(3),irstep)
      call write_dynvar(chymrst%ncid,chymrst%varid(4),port,irstep)
      call write_dynvar(chymrst%ncid,chymrst%varid(5),h2o,irstep)
      call closefile(chymrst%ncid)
      if ( itype == 1 ) then
        old = anno
      else
        old = mese
      end if
    end if
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
       mm = mese
       yy = anno
       oldanno = anno
       first = .false.
    end if
    if (anno /= oldanno .and. hourstep > 1) then
      yy = yy + 1
      write(cfl,'(i4,a)') yy-1,'000100'
      write (mess,'(15x,a,a)') 'Creating new Qmax file:  '// &
         trim(tdnsim),'_'//trim(cfl)//'_qmax.nc'
      write (6,'(a)') mess(1:len_trim(mess))
      call createfile_qmax(trim(tdnsim)//'_'//trim(cfl)//'_qmax.nc',time)
      call add_timestep(chymqmax%ncid,chymqmax%varid(3),iqstep)
      call write_dynvar(chymqmax%ncid,chymqmax%varid(4),port_qmaxs,iqstep)
      call closefile(chymqmax%ncid)
      port_qmaxs = 0
    end if
    oldanno = anno
    return
  end subroutine createnc3

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
      integer :: ncid, pstep, varid, start(4), count(4)
!
!-----------------------------------------------------------------------
!     Open netCDF file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_open(trim(tdnini), nf90_nowrite, ncid))
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
        if (.not. first .and. year == oldyear + 1) then
           call nio_check(nf90_close(ncid))
           first = .true.
           stepminus = step-1
           stepr = step - stepminus
        end if
        print*,"We are going to read step : ",stepr," of chym_runoff"
        if (first) then
          oldyear=year
          write(filerunoff,'(a,i4,a)') trim(tdninp)//'_',oldyear, &
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
        chym_runoff = chym_runoff*convfac
      end if

      call mpi_bcast(chym_runoff(1,1),nbc*nlc,MPI_REAL,0,mycomm,mpierr)

      end subroutine read_runoff

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
        if (present(nerr)) then
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
    istart(1) = 1
    istart(2) = 1
    istart(3) = iistep
    icount(1) = nlc
    icount(2) = nbc
    icount(3) = 1
    call nio_check(nf90_put_var(ncid,vname_id,ival,istart,icount))
  end subroutine write_dynvar_integer4

  real function geodistance(latt1,lonn1,latt2,lonn2)
    implicit none
    real, parameter :: rad = 6371000.0
    real, parameter :: dpi = 6.2831855
    real latt1,lonn1,latt2,lonn2,lt1,lt2,ln1,ln2,x,y
    lt1=latt1*dpi/360.
    lt2=latt2*dpi/360.
    ln1=lonn1*dpi/360.
    ln2=lonn2*dpi/360.
    if (abs(latt1-latt2).lt.0.2.and.abs(lonn1-lonn2).lt.0.2) then
        x=(rad*cos(lt1)*(ln1-ln2))*(rad*cos(lt2)*(ln1-ln2))
        y=(rad*(lt1-lt2))**2
        geodistance=sqrt(x+y)
    else
       x=sin(lt1)*sin(lt2)+cos(lt1)*cos(lt2)*cos((ln1)-(ln2))
       if (x.gt.1) x=1.0
       geodistance=acos(x)*rad
     endif
    if (geodistance.lt.0.1) geodistance=0.1
    return
  end function geodistance

end module mod_io
