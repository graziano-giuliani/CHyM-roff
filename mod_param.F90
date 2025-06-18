module mod_param

  implicit none

  real, parameter :: chym_acella = 6307.744**2
  real :: convfac = 1.0  ! for mm/s, kg/m^2/s. Set to 1/3.6 for m/h
  integer :: iostep , irstep , iqstep, ivar , rvar , qvar
  integer :: nlc,nbc,nlon,nlat
  integer :: inirun,step,yday,time,oldyear, ical
  integer :: hourstep, sdate, edate ,dstep
  integer :: hour,day,month,year,now
  integer , dimension(12) :: mesi
  data mesi/31 , 28 , 31 , 30 , 31 , 30 , 31 , 31 , 30 , 31 , 30 , 31/
  character(len=256) :: tdnstk,tdninp
  character(len=256) :: tdnini,tdnsim
  character(len=32) :: calendar
  character(len=256) :: filename,filenamerst,filenameqmax
  integer :: jahr1, jahr2, jahr3, jahr4
  integer :: isread, iorstfreq
  integer, parameter :: lu = 11

  integer ir(9),jr(9)
  data ir /-1, 0, 1, 1, 1, 0,-1,-1,0/
  data jr / 1, 1, 1, 0,-1,-1,-1, 0,0/

  real, parameter :: irloss = 0.05
  real, dimension(12), parameter :: irmonfac = &
    [ 0.00, 0.00, 0.00, 0.02, 0.10, 0.15, 0.20, 0.15, 0.02, 0.00, 0.00, 0.00 ]

  real :: thrriv = 5400.0

  real :: deltat

  ! area of the chym grid cells
  real, allocatable :: chym_area(:,:)
  real, allocatable :: alfa(:,:)
  integer, allocatable :: fmap(:,:)
  real, allocatable :: accl(:,:)
  integer, allocatable :: luse(:,:)
  logical, allocatable :: farm(:,:)
  real, allocatable :: port_sub(:,:)
  real, allocatable :: h2o_sub(:,:)
  real, allocatable :: bwet_sub(:,:)
  real, allocatable :: wkm1_sub(:,:)
  real, allocatable :: port(:,:)
  real, allocatable :: port_out(:,:)
  real, allocatable :: roff_out(:,:)
  real, allocatable :: port_qmaxs(:,:)
  real, allocatable :: wkm1(:,:)
  real, allocatable :: bwet(:,:)
  real, allocatable :: h2o(:,:)
  real, allocatable :: port1d(:)
  real, allocatable :: wkm11d(:)
  real, allocatable :: bwet1d(:)
  real, allocatable :: h2o1d(:)
  real, allocatable :: portsub1d(:)
  real, allocatable :: wkm1sub1d(:)
  real, allocatable :: bwetsub1d(:)
  real, allocatable :: h2osub1d(:)
  real, allocatable :: port1d_n(:)
  real, allocatable :: wkm11d_n(:)
  real, allocatable :: bwet1d_n(:)
  real, allocatable :: h2o1d_n(:)
  real, allocatable :: portsub1d_n(:)
  real, allocatable :: wkm1sub1d_n(:)
  real, allocatable :: bwetsub1d_n(:)
  real, allocatable :: h2osub1d_n(:)
  real, allocatable :: manning(:)
  ! area of the chym grid cells
  real, allocatable :: chym_drai(:,:)
  ! CHyM land-sea mask
  real, allocatable :: chym_lsm(:,:)
  ! CHyM dx
  real, allocatable :: chym_dx(:,:)
  ! CHyM lat
  real, allocatable :: chym_lat(:,:)
  real, allocatable :: corner_lat(:,:,:)
  ! CHyM lon
  real, allocatable :: chym_lon(:,:)
  real, allocatable :: corner_lon(:,:,:)

  real, allocatable :: chym_runoff(:,:)

  integer, parameter :: lntypes = 110
  integer, parameter :: ocean = 15
  integer, parameter :: lake = 14
  integer , dimension(3) :: chunksizes

#ifdef NILE
  real, parameter, dimension(12) :: nile_fresh_flux = &
      [ 336,  396,  407,  399,  472,  615,  &
        634,  557,  412,  373,  372,  352]
  integer :: idamietta = -1
  integer :: jdamietta = -1
  integer :: irosetta = -1
  integer :: jrosetta = -1
  real , parameter :: lat_damietta = 31.52
  real , parameter :: lon_damietta = 31.83
  real , parameter :: lat_rosetta = 31.46
  real , parameter :: lon_rosetta = 30.36
#endif

!
!-----------------------------------------------------------------------
!     CHyM model netCDF data type
!-----------------------------------------------------------------------
!
  type CHyM_IO
    integer :: ncid
    integer, allocatable :: dimid(:)
    integer, allocatable :: varid(:)
  end type CHyM_IO
!
  type(CHyM_IO) :: chymout, chymrst, chymqmax

  contains

  logical function is_inbox(lat,lon,clat,clon)
    implicit none
    real, intent(in) :: lat, lon
    real, dimension(4), intent(in) :: clat, clon
    real :: m1, b1, m2, b2, m3, b3, m4, b4
    logical :: l1, l2, l3, l4
    is_inbox = .false.

    ! ASSUME NON PATOLOGICAL CASES, i.e. clon, clat proper quadrilater.
    if ( clat(1) /= clat(2) ) then
      m1 = (clat(1)-clat(2))/(clon(1)-clon(2))
      b1 = (clon(1)*clat(2)-clon(2)*clat(1))/(clon(1)-clon(2))
      l1 = (lat > m1*lon+b1)
    else
      l1 = lat > clat(1)
    end if
    if ( clat(3) /= clat(4) ) then
      m2 = (clat(3)-clat(4))/(clon(3)-clon(4))
      b2 = (clon(3)*clat(4)-clon(4)*clat(3))/(clon(3)-clon(4))
      l2 = lat < m2*lon+b2
    else
      l2 = lat < clat(3)
    end if
    if ( clon(1) /= clon(4) ) then
      m3 = (clon(1)-clon(4))/(clat(1)-clat(4))
      b3 = (clat(1)*clon(4)-clat(4)*clon(1))/(clat(1)-clat(4))
      l3 = lon > m3*lat+b3
    else
      l3 = lon > clon(1)
    end if
    if ( clon(2) /= clon(3) ) then
      m4 = (clon(2)-clon(3))/(clat(2)-clat(3))
      b4 = (clat(2)*clon(3)-clat(3)*clon(2))/(clat(2)-clat(3))
      l4 = lon < m4*lat+b4
    else
      l4 = lon < clon(3)
    end if
    is_inbox = l1 .and. l2 .and. l3 .and. l4
  end function is_inbox

  subroutine find_nearest_land(i,j,ii,jj)
    implicit none
    integer, intent(in) :: i, j
    integer, intent(out) :: ii, jj
    if ( luse(i,j) /= ocean ) then
      ii = i
      jj = j
    else
      if ( all(luse(i-1:i+1,j-1:j+1) == ocean) ) then
        print *, 'No land point found! Will modify landmask!'
        ii = i
        jj = j
        return
      end if
      do jj = j-1, j+1
        do ii = i-1, i+1
          if ( luse(ii,jj) /= ocean ) then
            return
          end if
        end do
      end do
    end if
  end subroutine find_nearest_land

  subroutine find_nearest_ocean(i,j,ii,jj)
    implicit none
    integer, intent(in) :: i, j
    integer, intent(out) :: ii, jj

    if ( luse(i,j) == ocean ) then
      ii = i
      jj = j
    else
      if ( all(luse(i-1:i+1,j-1:j+1) /= ocean) ) then
        print *, 'No ocean point found! Will modify landmask!'
        ii = i
        jj = j
        return
      end if
      do jj = j-1, j+1
        do ii = i-1, i+1
          if ( luse(ii,jj) == ocean ) then
            return
          end if
        end do
      end do
    end if
  end subroutine find_nearest_ocean
end module mod_param
