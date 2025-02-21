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
program preproc

  use, intrinsic :: iso_fortran_env
  use mod_param
  use mod_crtstat
  use mod_ncio

  implicit none

  character(len=8) :: VERSION = '6.09_cpl'
  real :: maxlon , minlon
  real :: maxlat , minlat

  type(ncoutfile) :: outf

  write (output_unit,'(7x,a,1x,a,1x,a)') &
      'Start of CHYM preprocessing version', VERSION, 'fortran code'
  write (output_unit,'(2x,a,a,a,a,a,a)') &
        ' GIT Revision: ', GIT_REV, ' compiled on: ', __DATE__ , &
        ', at ', __TIME__

  call setparam

  write(output_unit, '(5x,a)') 'Creating static fields.'

  write(output_unit, '(7x,a)') 'Reading data from input files.'
  call grid_dimensions(gridfile,nlat,nlon)
  call getspace( )
  call grid_coordinates(gridfile,lat,lon,area)
  area = area * 1.0e-6 ! Make it in km2
  minlat = minval(lat)
  maxlat = maxval(lat)
  minlon = minval(lon)
  maxlon = maxval(lon)
  if ( maxlon-minlon > 5.0 .or. maxlat-minlat > 5.0 ) ifactor_res = 600
  write(output_unit, '(12x,a,2f8.3)') 'Latitude in range  ', minlat, maxlat
  write(output_unit, '(12x,a,2f8.3)') 'Longitude in range ', minlon, maxlon
  call ncfromfile(demfile,'dem',dem)
  call ncfromfile(landfile,'luc',luc)
  call ncfromfile(maskfile,'mask',mask)
  write(output_unit, '(7x,a)') 'Done.'

  call buildflowdirmap( )
  call areamatrix(area,drai)
  call buildacclivitymap
  call reconnectdem
  call runoffspeed

  outf%fname = outfile
  outf%fcomment = 'None'
  outf%nx = nlon
  outf%ny = nlat
  call create_outfile(outf)
  call outfile_attribute(outf,'minimum_latitude',minlat)
  call outfile_attribute(outf,'maximum_latitude',maxlat)
  call outfile_attribute(outf,'minimum_longitude',minlon)
  call outfile_attribute(outf,'maximum_longitude',maxlon)

  call write_variable(outf,latid,lat)
  call write_variable(outf,lonid,lon)
  call write_variable(outf,areaid,area)

  call add_variable(outf,demid)
  call add_variable(outf,lucid)
  call add_variable(outf,fdmid)
  call add_variable(outf,accid)
  call add_variable(outf,draid)
  call add_variable(outf,runid)
  call add_variable(outf,alfid)
  call write_variable(outf,demid,dem)
  call write_variable(outf,lucid,luc)
  call write_variable(outf,fdmid,fmap)
  call write_variable(outf,accid,accl)
  call write_variable(outf,draid,drai)
  call write_variable(outf,runid,runt)
  call write_variable(outf,alfid,alfa)

  call dispose_outfile(outf)

  call freespace( )
end program preproc
