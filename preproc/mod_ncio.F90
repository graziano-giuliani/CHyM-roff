!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP CHyM.
!
!    ICTP CHyM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP CHyM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP CHyM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_ncio
  use, intrinsic :: iso_fortran_env
  use netcdf

  implicit none

  private

  integer :: ncstatus

  integer, parameter :: maxdim = 8
  integer, parameter :: maxvar = 16

  type ncoutfile
    character(len=512) :: fname
    character(len=512) :: fcomment
    integer :: nx, ny, lntypes
    integer :: handler
    integer, dimension(maxdim) :: idimid
    integer, dimension(maxvar) :: ivarid
  end type ncoutfile

  integer, public, parameter :: latid = 1
  integer, public, parameter :: lonid = 2
  integer, public, parameter :: areaid = 3
  integer, public, parameter :: demid = 4
  integer, public, parameter :: lucid = 5
  integer, public, parameter :: fdmid = 6
  integer, public, parameter :: accid = 7
  integer, public, parameter :: draid = 8
  integer, public, parameter :: runid = 9
  integer, public, parameter :: alfid = 10
  integer, public, parameter :: ctrid = 11
  integer, public, parameter :: basid = 12
  integer, public, parameter :: mskid = 13

  type ncvariable
    character(len=16) :: vname
    integer :: vtype
    character(len=32) :: vstd
    character(len=64) :: vlong
    character(len=32) :: vunit
  end type ncvariable

  type(ncvariable), dimension(maxvar), parameter :: ncvars = &
      [ ncvariable('lat',nf90_float,'latitude','Latitude','degrees_north'),  &
        ncvariable('lon',nf90_float,'longitude','Longitude','degrees_east'), &
        ncvariable('aer',nf90_float,'cell_area','Cell Area','km2'),          &
        ncvariable('dem',nf90_float,'elevation','Elevation','m'),            &
        ncvariable('lus',nf90_int,'soil_category','Land Use Category','1'),  &
        ncvariable('fdm',nf90_int,'flow_direction','Flow direction','1'),    &
        ncvariable('acc',nf90_float,'slope',                                 &
                         'Tangent angle for flow direction','km2'),          &
        ncvariable('dra',nf90_float,'drainage_area','Drainage area','km2'),  &
        ncvariable('run',nf90_float,'runoff_time','Runoff Time','hours'),    &
        ncvariable('alf',nf90_float,'flow_velocity','Flow velocity','km/h'), &
        ncvariable('ctr',nf90_float,'flow_control','Flow control','m'),      &
        ncvariable('bas',nf90_int,'basin','Basin code','1'),                 &
        ncvariable('msk',nf90_int,'mask','Integration mask','1'),            &
        ncvariable('unk',nf90_float,'unknown','Unknown','1'),                &
        ncvariable('unk',nf90_float,'unknown','Unknown','1'),                &
        ncvariable('unk',nf90_float,'unknown','Unknown','1') ]

  public :: ncoutfile
  public :: create_outfile
  public :: dispose_outfile
  public :: add_variable
  public :: set_writemod
  public :: write_variable
  public :: outfile_attribute
  public :: grid_dimensions
  public :: grid_coordinates
  public :: ncfromfile

  interface outfile_attribute
    module procedure outfile_attribute_integer
    module procedure outfile_attribute_real
    module procedure outfile_attribute_text
  end interface

  interface ncfromfile
    module procedure integer_ncfromfile
    module procedure real_ncfromfile
  end interface ncfromfile

  interface write_variable
    module procedure write_variable_real_2d
    module procedure write_variable_integer_2d
  end interface write_variable

  interface read_variable
    module procedure read_real_var1d
    module procedure read_real_var2d
    module procedure read_real_var2d_from3d
    module procedure read_real_var2d_from4d
    module procedure read_real_var3d
    module procedure read_real_var4d
    module procedure read_integer_var2d
  end interface read_variable

  contains

  subroutine grid_dimensions(gridfile,nlat,nlon)
    implicit none
    character(len=*) , intent(in) :: gridfile
    integer, intent(out) :: nlat, nlon
    integer :: ncid
    ncstatus = nf90_open(gridfile,nf90_nowrite,ncid)
    call checkerror(__LINE__,'Cannot open',gridfile)
    nlon = dimlen(ncid,gridfile,'grid_xsize')
    nlat = dimlen(ncid,gridfile,'grid_ysize')
    ncstatus = nf90_close(ncid)
    call checkerror(__LINE__,'Cannot close',gridfile)
  end subroutine grid_dimensions

  subroutine grid_coordinates(gridfile,lat,lon,corner_lat,corner_lon,area)
    implicit none
    character(len=*) , intent(in) :: gridfile
    real, dimension(:,:), intent(out) :: lat, lon, area
    real, dimension(:,:,:), intent(out) :: corner_lat, corner_lon
    integer :: ncid
    ncstatus = nf90_open(gridfile,nf90_nowrite,ncid)
    call checkerror(__LINE__,'Cannot open',gridfile)
    call read_variable(ncid,gridfile,'grid_center_lon',lon)
    call read_variable(ncid,gridfile,'grid_center_lat',lat)
    call read_variable(ncid,gridfile,'grid_corner_lat',corner_lat)
    call read_variable(ncid,gridfile,'grid_corner_lon',corner_lon)
    call read_variable(ncid,gridfile,'cell_area',area)
    ncstatus = nf90_close(ncid)
    call checkerror(__LINE__,'Cannot close',gridfile)
  end subroutine grid_coordinates

  subroutine real_ncfromfile(filename,varname,thevar)
    implicit none
    character(len=*) , intent(in) :: filename, varname
    real, dimension(:,:) :: thevar
    integer :: ncid
    ncstatus = nf90_open(filename,nf90_nowrite,ncid)
    call checkerror(__LINE__,'Cannot open',filename)
    call read_variable(ncid,filename,varname,thevar)
    ncstatus = nf90_close(ncid)
    call checkerror(__LINE__,'Cannot close',filename)
  end subroutine real_ncfromfile

  subroutine integer_ncfromfile(filename,varname,thevar)
    implicit none
    character(len=*) , intent(in) :: filename, varname
    integer, dimension(:,:) :: thevar
    integer :: ncid
    ncstatus = nf90_open(filename,nf90_nowrite,ncid)
    call checkerror(__LINE__,'Cannot open',filename)
    call read_variable(ncid,filename,varname,thevar)
    ncstatus = nf90_close(ncid)
    call checkerror(__LINE__,'Cannot close',filename)
  end subroutine integer_ncfromfile

  integer function dimlen(ncid,fname,dname) result(dlen)
    implicit none
    integer , intent(in) :: ncid
    character(len=*) , intent(in) :: fname, dname
    integer :: idimid
    ncstatus = nf90_inq_dimid(ncid,dname,idimid)
    call checkerror(__LINE__,'Cannot read dim '//dname,fname)
    ncstatus = nf90_inquire_dimension(ncid,idimid,len=dlen)
    call checkerror(__LINE__,'Cannot read lenght of dim '//dname,fname)
  end function dimlen

  subroutine read_real_var1d(ncid,fname,vname,val)
    implicit none
    integer , intent(in) :: ncid
    character(len=*) , intent(in) :: fname, vname
    real, dimension(:), intent(out) :: val
    integer :: ivarid
    ncstatus = nf90_inq_varid(ncid,vname,ivarid)
    call checkerror(__LINE__,'Cannot find var '//vname,fname)
    ncstatus = nf90_get_var(ncid,ivarid,val)
    call checkerror(__LINE__,'Cannot read var '//vname,fname)
  end subroutine read_real_var1d

  subroutine read_real_var2d(ncid,fname,vname,val)
    implicit none
    integer , intent(in) :: ncid
    character(len=*) , intent(in) :: fname, vname
    real, dimension(:,:), intent(out) :: val
    integer :: ivarid
    ncstatus = nf90_inq_varid(ncid,vname,ivarid)
    call checkerror(__LINE__,'Cannot find var '//vname,fname)
    ncstatus = nf90_get_var(ncid,ivarid,val)
    call checkerror(__LINE__,'Cannot read var '//vname,fname)
  end subroutine read_real_var2d

  subroutine read_real_var2d_from3d(ncid,fname,vname,val,irec)
    implicit none
    integer , intent(in) :: ncid , irec
    character(len=*) , intent(in) :: fname, vname
    real, dimension(:,:), intent(out) :: val
    integer :: ivarid
    integer , dimension(3) :: istart , icount
    istart(1) = 1
    istart(2) = 1
    istart(3) = irec
    icount(1) = size(val,1)
    icount(2) = size(val,2)
    icount(3) = 1
    ncstatus = nf90_inq_varid(ncid,vname,ivarid)
    call checkerror(__LINE__,'Cannot find var '//vname,fname)
    ncstatus = nf90_get_var(ncid,ivarid,val,istart,icount)
    call checkerror(__LINE__,'Cannot read var '//vname,fname)
  end subroutine read_real_var2d_from3d

  subroutine read_real_var2d_from4d(ncid,fname,vname,val,irec1,irec2)
    implicit none
    integer , intent(in) :: ncid , irec1 , irec2
    character(len=*) , intent(in) :: fname, vname
    real, dimension(:,:), intent(out) :: val
    integer :: ivarid
    integer , dimension(4) :: istart , icount
    istart(1) = 1
    istart(2) = 1
    istart(3) = irec1
    istart(4) = irec2
    icount(1) = size(val,1)
    icount(2) = size(val,2)
    icount(3) = 1
    icount(4) = 1
    ncstatus = nf90_inq_varid(ncid,vname,ivarid)
    call checkerror(__LINE__,'Cannot find var '//vname,fname)
    ncstatus = nf90_get_var(ncid,ivarid,val,istart,icount)
    call checkerror(__LINE__,'Cannot read var '//vname,fname)
  end subroutine read_real_var2d_from4d

  subroutine read_real_var3d(ncid,fname,vname,val)
    implicit none
    integer , intent(in) :: ncid
    character(len=*) , intent(in) :: fname, vname
    real, dimension(:,:,:), intent(out) :: val
    integer :: ivarid
    ncstatus = nf90_inq_varid(ncid,vname,ivarid)
    call checkerror(__LINE__,'Cannot find var '//vname,fname)
    ncstatus = nf90_get_var(ncid,ivarid,val)
    call checkerror(__LINE__,'Cannot read var '//vname,fname)
  end subroutine read_real_var3d

  subroutine read_real_var4d(ncid,fname,vname,val)
    implicit none
    integer , intent(in) :: ncid
    character(len=*) , intent(in) :: fname, vname
    real, dimension(:,:,:,:), intent(out) :: val
    integer :: ivarid
    ncstatus = nf90_inq_varid(ncid,vname,ivarid)
    call checkerror(__LINE__,'Cannot find var '//vname,fname)
    ncstatus = nf90_get_var(ncid,ivarid,val)
    call checkerror(__LINE__,'Cannot read var '//vname,fname)
  end subroutine read_real_var4d

  subroutine read_integer_var2d(ncid,fname,vname,val)
    implicit none
    integer , intent(in) :: ncid
    character(len=*) , intent(in) :: fname, vname
    integer, dimension(:,:), intent(out) :: val
    integer :: ivarid
    ncstatus = nf90_inq_varid(ncid,vname,ivarid)
    call checkerror(__LINE__,'Cannot find var '//vname,fname)
    ncstatus = nf90_get_var(ncid,ivarid,val)
    call checkerror(__LINE__,'Cannot read var '//vname,fname)
  end subroutine read_integer_var2d

  subroutine create_outfile(finfo)
    implicit none
    type(ncoutfile), intent(inout) :: finfo
    character(len=64) :: dstring
    integer, dimension(8) :: tvals
    call date_and_time(values=tvals)
    write(dstring,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)') &
        tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,       &
        tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,             &
        ' : Created by CHyM preproc program'
    ncstatus = nf90_create(finfo%fname,&
                ior(nf90_netcdf4,nf90_clobber),finfo%handler)
    ncstatus = nf90_def_dim(finfo%handler,'corners',4,finfo%idimid(1))
    call checkerror(__LINE__,'Cannot create dimension lntypes',finfo%fname)
    call checkerror(__LINE__,'Cannot create file',finfo%fname)
    ncstatus = nf90_def_dim(finfo%handler,'lon',finfo%nx,finfo%idimid(2))
    call checkerror(__LINE__,'Cannot create dimension lon',finfo%fname)
    ncstatus = nf90_def_dim(finfo%handler,'lat',finfo%ny,finfo%idimid(3))
    call checkerror(__LINE__,'Cannot create dimension lon',finfo%fname)
    ncstatus = nf90_def_dim(finfo%handler,'lntypes',finfo%lntypes, &
                            finfo%idimid(4))
    call checkerror(__LINE__,'Cannot create dimension lntypes',finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,nf90_global,'Conventions','CF-1.6')
    call checkerror(__LINE__,'Cannot add basic attribute',finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,nf90_global,'history',trim(dstring))
    call checkerror(__LINE__,'Cannot add basic attribute',finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,nf90_global,'institution','ICTP')
    call checkerror(__LINE__,'Cannot add basic attribute',finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,nf90_global,'source','ChyM model')
    call checkerror(__LINE__,'Cannot add basic attribute',finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,nf90_global,'title', &
                       'ChyM static input basins dataset')
    call checkerror(__LINE__,'Cannot add basic attribute',finfo%fname)
    call add_variable(finfo,latid)
    call add_variable(finfo,lonid)
    call add_variable(finfo,areaid)
    call add_variable(finfo,mskid)
  end subroutine create_outfile

  subroutine set_writemod(finfo,manning,lat,lon,clat,clon,area,mask)
    implicit none
    type(ncoutfile), intent(inout) :: finfo
    real, dimension(:,:), intent(in) :: lat, lon, area
    integer, dimension(:,:), intent(in) :: mask
    real, dimension(:,:,:), intent(in) :: clon, clat
    real, dimension(:), intent(in) :: manning
    integer :: varid(3)
    ncstatus = nf90_def_var(finfo%handler,"corner_lat",      &
                            nf90_float, finfo%idimid(1:3), varid(1))
    call checkerror(__LINE__, 'Cannot add variable corner lat', finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,varid(1), &
                            'standard_name','latitude_bounds')
    call checkerror(__LINE__, 'Cannot add attribute to corner lat', finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,varid(1), &
                            'long_name','Latitude Bounds')
    call checkerror(__LINE__, 'Cannot add attribute to corner lat', finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,varid(1), 'units','degrees_north')
    call checkerror(__LINE__, 'Cannot add attribute to corner lat', finfo%fname)

    ncstatus = nf90_def_var(finfo%handler,"corner_lon",      &
                            nf90_float, finfo%idimid(1:3), varid(2))
    call checkerror(__LINE__, 'Cannot add variable corner lon', finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,varid(2), &
                            'standard_name','longitude_bounds')
    call checkerror(__LINE__, 'Cannot add attribute to corner lon', finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,varid(2), &
                            'long_name','Longitude Bounds')
    call checkerror(__LINE__, 'Cannot add attribute to corner lon', finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,varid(2), 'units','degrees_east')
    call checkerror(__LINE__, 'Cannot add attribute to corner lon', finfo%fname)

    ncstatus = nf90_def_var(finfo%handler,"manning",      &
                            nf90_float, finfo%idimid(4), &
                            varid(3))
    call checkerror(__LINE__, &
                  'Cannot add variable manning', finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,varid(3), &
                            'standard_name','manning_coefficient')
    call checkerror(__LINE__, &
                  'Cannot add standard name to manning', finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,varid(3), &
                            'long_name','Manning coefficient')
    call checkerror(__LINE__, &
                  'Cannot add long name to manning', finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,varid(3), &
                            'units','sm^1/3')
    call checkerror(__LINE__, &
                  'Cannot add units name to manning', finfo%fname)

    ncstatus = nf90_enddef(finfo%handler)
    call checkerror(__LINE__, 'Cannot put file in write mode', finfo%fname)

    ncstatus = nf90_put_var(finfo%handler,finfo%ivarid(latid),lat)
    call checkerror(__LINE__, 'Cannot write lat in file', finfo%fname)
    ncstatus = nf90_put_var(finfo%handler,finfo%ivarid(lonid),lon)
    call checkerror(__LINE__, 'Cannot write lon in file', finfo%fname)
    ncstatus = nf90_put_var(finfo%handler,finfo%ivarid(areaid),area)
    call checkerror(__LINE__, 'Cannot write area in file', finfo%fname)
    ncstatus = nf90_put_var(finfo%handler,finfo%ivarid(mskid),mask)
    call checkerror(__LINE__, 'Cannot write mask in file', finfo%fname)
    ncstatus = nf90_put_var(finfo%handler,varid(1),clat)
    call checkerror(__LINE__, 'Cannot write corner lat in file', finfo%fname)
    ncstatus = nf90_put_var(finfo%handler,varid(2),clon)
    call checkerror(__LINE__, 'Cannot write corner lon in file', finfo%fname)
    ncstatus = nf90_put_var(finfo%handler,varid(3),manning)
    call checkerror(__LINE__, 'Cannot write manning in file', finfo%fname)
  end subroutine set_writemod

  subroutine write_variable_real_2d(finfo,vid,val)
    implicit none
    type(ncoutfile), intent(inout) :: finfo
    integer, intent(in) :: vid
    real, dimension(:,:) :: val
    ncstatus = nf90_put_var(finfo%handler,finfo%ivarid(vid),val)
    call checkerror(__LINE__, &
                  'Cannot write variable '//ncvars(vid)%vname, &
                  finfo%fname)
  end subroutine write_variable_real_2d

  subroutine write_variable_integer_2d(finfo,vid,val)
    implicit none
    type(ncoutfile), intent(inout) :: finfo
    integer, intent(in) :: vid
    integer, dimension(:,:) :: val
    ncstatus = nf90_put_var(finfo%handler,finfo%ivarid(vid),val)
    call checkerror(__LINE__, &
                  'Cannot write variable '//ncvars(vid)%vname, &
                  finfo%fname)
  end subroutine write_variable_integer_2d

  subroutine add_variable(finfo,vid)
    implicit none
    type(ncoutfile), intent(inout) :: finfo
    integer, intent(in) :: vid
    ncstatus = nf90_def_var(finfo%handler,ncvars(vid)%vname,      &
                            ncvars(vid)%vtype, finfo%idimid(2:3), &
                            finfo%ivarid(vid))
    call checkerror(__LINE__, &
                  'Cannot add basic variable '//ncvars(vid)%vname, &
                  finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,finfo%ivarid(vid), &
                            'standard_name',ncvars(vid)%vstd)
    call checkerror(__LINE__, &
                  'Cannot add standard name to '//ncvars(vid)%vname, &
                   finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,finfo%ivarid(vid), &
                            'long_name',ncvars(vid)%vlong)
    call checkerror(__LINE__, &
                  'Cannot add long name to '//ncvars(vid)%vname, &
                  finfo%fname)
    ncstatus = nf90_put_att(finfo%handler,finfo%ivarid(vid), &
                            'units',ncvars(vid)%vunit)
    call checkerror(__LINE__, &
                  'Cannot add units to '//ncvars(vid)%vname, &
                  finfo%fname)
    if ( ncvars(vid)%vname /= "lon" .and. ncvars(vid)%vname /= "lat" ) then
      ncstatus = nf90_put_att(finfo%handler,finfo%ivarid(vid), &
                              'coordinates',"lat lon")
      call checkerror(__LINE__, &
                    'Cannot add coordinates to '//ncvars(vid)%vname, &
                    finfo%fname)
    end if
  end subroutine add_variable

  subroutine dispose_outfile(finfo)
    implicit none
    type(ncoutfile), intent(inout) :: finfo
    ncstatus = nf90_close(finfo%handler)
    call checkerror(__LINE__,'Cannot close file',finfo%fname)
    finfo%handler = -1
  end subroutine dispose_outfile

  subroutine outfile_attribute_integer(finfo,aname,aval)
    type(ncoutfile), intent(in) :: finfo
    character(len=*), intent(in) :: aname
    integer, intent(in) :: aval
    ncstatus = nf90_put_att(finfo%handler,nf90_global,aname,aval)
    call checkerror(__LINE__,'Cannot add attribute '//aname,finfo%fname)
  end subroutine outfile_attribute_integer

  subroutine outfile_attribute_real(finfo,aname,aval)
    type(ncoutfile), intent(in) :: finfo
    character(len=*), intent(in) :: aname
    real, intent(in) :: aval
    ncstatus = nf90_put_att(finfo%handler,nf90_global,aname,aval)
    call checkerror(__LINE__,'Cannot add attribute '//aname,finfo%fname)
  end subroutine outfile_attribute_real

  subroutine outfile_attribute_text(finfo,aname,aval)
    type(ncoutfile), intent(in) :: finfo
    character(len=*), intent(in) :: aname
    character(len=*), intent(in) :: aval
    ncstatus = nf90_put_att(finfo%handler,nf90_global,aname,aval)
    call checkerror(__LINE__,'Cannot add attribute '//aname,finfo%fname)
  end subroutine outfile_attribute_text

  subroutine checkerror(line,msg,iofile)
    implicit none
    integer , intent(in) :: line
    character(len=*) , intent(in) :: iofile, msg

    if ( ncstatus /= nf90_noerr ) then
      write(error_unit, *) 'Operating on ',trim(iofile)
      write(error_unit, *) 'At line ',line," : ",msg
      write(error_unit, *) 'NetCDF library error  : ',nf90_strerror(ncstatus)
      stop
    end if
  end subroutine checkerror

end module mod_ncio
