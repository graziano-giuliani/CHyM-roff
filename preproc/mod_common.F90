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

module mod_common

  use, intrinsic :: iso_fortran_env

  implicit none

  private

  integer, public, dimension(9), parameter :: ir = [-1, 0, 1, 1, 1, 0,-1,-1, 0]
  integer, public, dimension(9), parameter :: jr = [ 1, 1, 1, 0,-1,-1,-1, 0, 0]

  real, parameter :: erad = 6371000.0
  real, parameter :: pi = 4*atan(1.0)
  real, parameter :: rad2deg = 180.0/pi
  real, parameter :: deg2rad = pi/180.0
  real, parameter :: dpi = 2.0*pi

  integer, public :: nlat
  integer, public :: nlon
  integer, public :: n2df
  integer, public :: n4df
  integer, public :: ifactor_res = 300
  integer, public , parameter :: ocean = 15
  integer, public , parameter :: iwater = 14
  integer, public , parameter :: iriver = 100
  integer, public , parameter :: lntypes = 110

  real, dimension(:,:) , allocatable, public :: lat
  real, dimension(:,:) , allocatable, public :: lon
  real, dimension(:,:) , allocatable, public :: area
  real, dimension(:,:) , allocatable, public :: dem
  integer, dimension(:,:) , allocatable, public :: luc
  integer, dimension(:,:) , allocatable, public :: mask

  real, dimension(:,:) , allocatable, public :: dx
  integer, dimension(:,:) , allocatable, public :: fmap
  integer, dimension(:,:) , allocatable, public :: noflow
  real, dimension(:,:) , allocatable, public :: accl
  real, dimension(:,:) , allocatable, public :: drai
  real, dimension(:,:) , allocatable, public :: alfa
  real, dimension(:,:) , allocatable, public :: runt

  real, public, dimension(lntypes) :: manning

  character(len=512) , public :: gridfile
  character(len=512) , public :: demfile
  character(len=512) , public :: landfile
  character(len=512) , public :: maskfile
  character(len=512) , public :: outfile

  ! Return flow factor
  real , public :: cpar1  = 4.8e-07
  ! Alpha coefficients for hydraulic radius
  real , public :: cpar2  = 0.0015
  ! Beta coefficients for hydraulic radius
  real , public :: cpar3  = 0.050
  ! Melting temperature factor
  real , public :: cpar4  = 0.050
  ! Melting shortwave rad. factor
  real , public :: cpar5  = 0.0094
  ! River/land threshold (Km2)
  real , public :: cpar6  = 500.0
  ! Number of days to consider for return flow
  real , public :: cpar7  = 90.0
  ! Reduction of land/channel manning coefficient
  real , public :: cpar8  = 4.5
  ! River/land threshold (Km2) for returnflowcient
  real , public :: cpar9  = 200.0

  integer, public :: angiocycle = 1
  integer, public :: angionp = 40
  integer, public :: ncyc1 = 100
  integer, public :: ncyc2 = 1000
  integer, public :: ncyc3 = 21

  public :: read_namelist, getspace, freespace

  public :: distance

  contains

  subroutine read_namelist
    implicit none
    integer :: lun
    integer :: iretval
    character(len=512) :: namelistfile

    namelist /chymconfig/ gridfile, demfile, landfile, maskfile, outfile
    namelist /chymparam/ cpar1, cpar2, cpar3, cpar4, cpar5, &
                         cpar6, cpar7, cpar8, cpar9, &
                         ncyc1, ncyc2, ncyc3, &
                         angiocycle, angionp
    namelist /manningparam/ manning

    call getarg(1, namelistfile)
    open(newunit=lun, file=namelistfile, status='old', &
         action='read', iostat=iretval)
    if ( iretval /= 0 ) then
      write(error_unit,*) 'Error opening namelist file'//trim(namelistfile)
      stop
    end if
    write(output_unit,'(7x,a)') 'Reading model namelist'
    read(lun, nml=chymconfig, iostat=iretval)
    if ( iretval /= 0 ) then
      write(error_unit,*) 'Error reading chymconfig namelist'
      stop
    end if
    rewind(lun)
    read(lun, nml=chymparam, iostat=iretval)
    if ( iretval /= 0 ) then
      write(error_unit,*) 'Error reading chymparam namelist'
      stop
    end if
    rewind(lun)
    read(lun, nml=manningparam, iostat=iretval)
    if ( iretval /= 0 ) then
      write(error_unit,*) 'Error reading chymparam namelist'
      stop
    end if
  end subroutine read_namelist

  real function distance(lat1,lon1,lat2,lon2)
    implicit none
    real , intent(in) :: lat1 , lon1 , lat2 , lon2
    real :: lt1 , lt2 , ln1 , ln2 , x , y
    lt1 = lat1*deg2rad
    lt2 = lat2*deg2rad
    ln1 = lon1*deg2rad
    ln2 = lon2*deg2rad
    if ( abs(lat1-lat2) < 0.2 .and. abs(lon1-lon2) < 0.2 ) then
      x = (erad*cos(lt1)*(ln1-ln2))*(erad*cos(lt2)*(ln1-ln2))
      y = (erad*(lt1-lt2))**2
      distance = sqrt(x+y)
    else
      x = sin(lt1)*sin(lt2) + cos(lt1)*cos(lt2)*cos((ln1)-(ln2))
      if ( x > 1 ) x = 1.0
      distance = acos(x)*erad
    endif
    if ( distance < 0.1 ) distance = 0.1
  end function distance

  subroutine getspace
    implicit none
    allocate(lat(nlon,nlat))
    allocate(lon(nlon,nlat))
    allocate(area(nlon,nlat))
    allocate(dem(nlon,nlat))
    allocate(luc(nlon,nlat))
    allocate(mask(nlon,nlat))

    allocate(dx(nlon,nlat))
    allocate(noflow(nlon,nlat))
    allocate(fmap(nlon,nlat))
    allocate(accl(nlon,nlat))
    allocate(drai(nlon,nlat))
    allocate(alfa(nlon,nlat))
    allocate(runt(nlon,nlat))
  end subroutine getspace

  subroutine freespace
    implicit none
    deallocate(lat)
    deallocate(lon)
    deallocate(area)
    deallocate(dem)
    deallocate(luc)
    deallocate(mask)

    deallocate(dx)
    deallocate(noflow)
    deallocate(fmap)
    deallocate(accl)
    deallocate(drai)
    deallocate(alfa)
    deallocate(runt)
  end subroutine freespace

end module mod_common
