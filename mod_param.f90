      module mod_param

      implicit none

      real, parameter :: chym_acella = 6307.744**2
      integer :: iostep , irstep , iqstep, ivar , rvar , qvar
      integer :: nlc,nbc,step,nstep,nlon,nlat
      integer :: inirun,yday,time,oldyear, ical
      integer :: hourstep, sdate, edate ,dstep
      integer :: hour,day,month,year,now
      integer , dimension(12) :: mesi
      data mesi/31 , 28 , 31 , 30 , 31 , 30 , 31 , 31 , 30 , 31 , 30 , 31/
      character(len=256) :: filestatic,filemrro,chym_output
      character(len=256) :: chym_restart,chym_initial,tdate
      character(len=256) :: chym_qmax,tsdate,sim_name,calendario
      character(len=256) :: filename,filenamerst,filenameqmax
      integer :: jahr1, jahr2, jahr3, jahr4
      integer :: isread, iswrit
      integer, parameter :: lu = 11

      integer ir(9),jr(9)
      data ir /-1, 0, 1, 1, 1, 0,-1,-1,0/
      data jr / 1, 1, 1, 0,-1,-1,-1, 0,0/

      integer :: pstep
      integer :: deltat

      ! area of the chym grid cells
      real, allocatable :: chym_area(:,:)
      real, allocatable :: alfa(:,:)
      real, allocatable :: fmap(:,:)
      real, allocatable :: accl(:,:)
      real, allocatable :: luse(:,:)
      real, allocatable :: port_sub(:,:)
      real, allocatable :: h2o_sub(:,:)
      real, allocatable :: bwet_sub(:,:)
      real, allocatable :: wkm1_sub(:,:)
      real, allocatable :: port(:,:)
      integer, allocatable :: port_out(:,:)
      real, allocatable :: port_outs(:,:)
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
      ! CHyM lon
      real, allocatable :: chym_lon(:,:)
      real, allocatable :: lat1(:)
      ! CHyM lon1
      real, allocatable :: lon1(:)
      real, allocatable :: chym_runoff(:,:)

      integer, parameter :: lntypes = 110
      integer, parameter :: mare = 15
      integer, parameter :: lago = 14
      integer , dimension(3) :: chunksizes
      real,parameter :: scale_factor=0.00011641532188114492
      real,parameter :: add_offset=250000

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
      end module mod_param
