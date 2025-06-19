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

module mod_crtstat

  use, intrinsic :: iso_fortran_env
  use mod_common
  use mod_cellaut

  private

  public :: applymask
  public :: buildflowdirmap, buildacclivitymap
  public :: areamatrix, runoffspeed
  public :: reconnectdem

  contains

  subroutine applymask
    implicit none
    integer :: i , j
    integer , parameter :: fillin = 76
    do j = 1 , nlat
      do i = 1 , nlon
        if ( mask(i,j) < 2 ) luc(i,j) = 15
        if ( dem(i,j) <= 0.0 ) dem(i,j) = 1.0
      end do
    end do
    do j = 1 , nlat
      do i = 1 , nlon
        if ( luc(i,j) == 15 .and. mask(i,j) == 2 ) luc(i,j) = fillin
        if ( luc(i,j) == 14 ) luc(i,j) = fillin
      end do
    end do
    do j = 1 , nlat
      do i = 1 , nlon
        if ( luc(i,j) == 15 ) dem(i,j) = 0.0
      end do
    end do
  end subroutine applymask

  logical function inside(n,n1,n2)
    implicit none
    integer, intent(in) :: n , n1 , n2
    inside = .true.
    if ( n < n1 .or. n > n2 ) inside = .false.
  end function inside

  subroutine demholefilling
    implicit none
    real , dimension(8) :: hc
    integer :: i , j , ihole , k
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
        if ( luc(i,j) /= ocean ) then
          ihole = 0
          do k = 1 , 8
            hc(k) = dem(i+ir(k),j+jr(k))
            if ( dem(i,j) < dem(i+ir(k),j+jr(k)) ) ihole = ihole + 1
          end do
          if ( ihole == 8 ) dem(i,j) = minval(hc(1:8)) + 1
        end if
      end do
    end do
  end subroutine demholefilling

  subroutine buildflowdirmap
    implicit none
    character(len=2) , dimension(8) :: col
    integer :: i , ii , j , ncyc , nzero , nzero2
    integer , dimension(8) :: icol
    real , dimension(0:8) :: rk
    real , dimension(nlon,nlat) :: plot , work
    integer , dimension(nlon,nlat) :: icl
    character(len=80) :: title
    logical :: plotgeop , plotlake
    real , dimension(nlon,nlat) :: savedem
    data icol/5 , 7 , 9 , 11 , 3 , 17 , 4 , 2/
    data col/'NW' , 'N' , 'NE' , 'E' , 'SE' , 'S' , 'SW' , 'W'/
    data plotlake/.false./
    data plotgeop/.true./

    write(output_unit,'(/5x,a)') 'Building Flow Direction Map.'
    rk(0) = 0.0
    nzero = 0
    do i = 1 , nlon
      fmap(i,1) = -1
      fmap(i,nlat) = -1
      noflow(i,1) = 0
      noflow(i,nlat) = 0
    end do
    do j = 1 , nlat
      fmap(1,j) = -1
      fmap(nlon,j) = -1
      noflow(1,j) = 0
      noflow(nlon,j) = 0
    end do
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
        if ( luc(i,j) /= ocean ) then
          noflow(i,j) = 8
        else
          noflow(i,j) = 0
        end if
      end do
    end do

    savedem = dem

    ncyc = ncyc1
    write(output_unit,'(7x,a,i8,a)') &
        'DEM Smoothing by CA algorithm:', ncyc, ' cycles'
    do ii = 1 , ncyc
      call d2cellcycle(dem,noflow,plot,nlon,nlat,0.005)
    end do
    ncyc = ncyc2
    write(output_unit,'(7x,a,i8,a)') &
        'DEM Smoothing by filling algorithm:',ncyc,' cycles'
    do ii = 1 , ncyc
      call demholefilling
    end do

    call establishfowdir
    call noflowcounter(noflow,nzero)
    write(output_unit,'(7x,a,i8)') 'No-flow points are',nzero

    ii = 0
    do while ( ii < ncyc3 .and. nzero > 0 )
      ii = ii + 1
      call d2cellcycle(dem,noflow,plot,nlon,nlat,0.035)
      call establishfowdir
      call noflowcounter(noflow,nzero)
      if ( mod(ii,5) == 0 ) then
        write(title,'(i10,a)') nzero , ' No-flow points after'
        write(output_unit,'(a,i8,a)') title, ii,' CA Algorithm Correction'
      end if
    end do
    if ( nzero > 0 ) then
      write(output_unit,'(15x,a)') &
         'Correcting flow-dir with direction-algorithm.'
      call dircorrflow(icl)
      call dircorrflow(icl)
      call noflowcounter(icl,nzero2)
      write(output_unit,'(7x,a,i8)') 'No-flow points are now ',nzero2
    end if
    plot = 1.0
    call flowcheck(plot,work)
    write(output_unit,'(12x,a)') 'Done.'
  end subroutine buildflowdirmap

  subroutine areamatrix(area,drai)
    implicit none
    real, dimension(nlon,nlat), intent(in) :: area
    real, dimension(nlon,nlat), intent(inout) :: drai
    real, dimension(nlon,nlat) :: wrk1
    wrk1 = area
    call flowcheck(wrk1,drai)
  end subroutine areamatrix

  subroutine buildacclivitymap
    implicit none
    real :: thresh
    logical, save :: first = .true.
    integer :: i, j, k, ncor, nnofl
    real, dimension(nlon,nlat) :: wk, wk1
    if ( first ) then
      write(output_unit,'(/5x,a)') 'Building Incline Map.'
    else
      write(output_unit,'(/12x,a)') 'Correcting Incline Map.'
    end if
    accl(:,:) = 0.0
    wk1(:,:) = 0.0
    ncor = 0
    nnofl = 0
    do j = 2 , nlat - 1
      do i = 2 , nlon - 1
        k = fmap(i,j)
        if ( luc(i,j)/=ocean .and. k>0 .and. k<=8 ) then
          accl(i,j) = (dem(i,j)-dem(i+ir(k),j+jr(k)))  &
                      /distance(lat(i,j),lon(i,j),lat(i+ir(k),j+jr(k)), &
                      lon(i+ir(k),j+jr(k)))
          if ( accl(i,j) <= 1.0E-05 ) then
            noflow(i,j) = 8
            ncor = ncor + 1
            accl(i,j) = 0.0
          else
            noflow(i,j) = 0
          end if
        else if ( luc(i,j)/=ocean .and. luc(i,j)/=iwater .and. k==0 ) then
          nnofl = nnofl + 1
          noflow(i,j) = 8
          accl(i,j) = 0.0
        else if ( luc(i,j)==iriver .or. luc(i,j)==iwater ) then
          accl(i,j) = 0.0
          noflow(i,j) = 8
        else if ( luc(i,j)==ocean ) then
          accl(i,j) = 0.0
          noflow(i,j) = 0
        else
          write(output_unit,'(10x,a,4i8)') &
            ' Flux error inside buildacclivitymap.' , i , j , k , luc(i,j)
          nnofl = nnofl + 1
          noflow(i,j) = 8
          accl(i,j) = 0.0
        end if
        wk1(i,j) = fmap(i,j)
      end do
    end do
    do j = 1 , nlat
      noflow(1,j) = -1
      noflow(nlon,j) = -1
      accl(1,j) = 0.0
      accl(nlon,j) = 0.0
    end do
    do i = 1 , nlon
      noflow(i,1) = -1
      noflow(i,nlat) = -1
      accl(i,1) = 0.0
      accl(i,nlat) = 0.0
    end do
    write(output_unit,'(13x,i5,a)') &
        ncor , ' cells have acclivity lower than 1E-05.'
    write(output_unit,'(13x,i5,a)') nnofl , ' points have no flow direction.'
    write(output_unit,'(15x,a)') 'Correcting Incline field with CA algorithm.'

! The following CA cycles affect only the grid points for which
! acclivity has not been defined (noflow=8)
!    if ( mchym(13)/=10 ) then
      do i = 1 , 100
        call d2cellcycle(accl,noflow,wk,nlon,nlat,0.1)
      end do
!    end if

! This control has been eliminated, in few situation at least, it seems
! that increase the number of steps of integration is enough to
! solve numerica singularities along  the river path.
    ncor = 0
    thresh = 0.00030
    do i = 1 , nlon
      do j = 1 , nlat
        if ( accl(i,j)>thresh .and. drai(i,j)>200.0 ) then
          ncor = ncor + 1
!          accl(i,j)=thresh
        end if
      end do
    end do
    write(output_unit,'(13x,i8,a,f8.4)') &
      ncor,' cells had acclivity greater than ',thresh
    first = .false.
    write(output_unit,'(12x,a)') 'Done.'
  end subroutine buildacclivitymap

  subroutine reconnectdem
    implicit none
    integer :: i, iter, j
    logical :: modify
    integer :: nsec
    real , dimension(nlon,nlat) :: work , plot
    integer, dimension(nlon*nlat) :: isec , jsec
    character(len=60) :: str1 , str2
    write(output_unit,'(/12x,a)') 'Reconnecting severe DEM singularities.'
    modify = .false.
    do iter = 1 , angiocycle
      write(output_unit,'(15x,a,i1)') 'Iteration number ' , iter
      call findthemouths(drai,100.,0,nsec,isec,jsec)
      call basinpaint(dem,fmap,luc,work,nlon,nlat,1,isec,jsec,plot,0)
      do i = 2 , nlon - 1
        do j = 2 , nlat - 1
          if ( work(i,j) > 10.0 .and. nint(plot(i,j)) == -5 ) then
            write(str1,'(i5,a,i5)') i , '-' , j
            write(str2,'(a,f10.1,a)') 'CHyM Angioplasty to cell: '//   &
                     trim(str1)//'(' , work(i,j) ,')'
            write(output_unit,'(15x,a)') str2
            call angioplasty(dem,fmap,plot,nlon,nlat,i,j)
            modify = .true.
          end if
        end do
      end do
    end do
    write(output_unit,'(12x,a)') 'Done.'
    if ( modify ) then
      call areamatrix(area,drai)
      call buildacclivitymap
    end if
  end subroutine reconnectdem

  subroutine angioplasty(dem,fmap,wk1,nlon,nlat,i,j)
    implicit none
    integer, intent(in) :: i , j , nlat , nlon
    real, dimension(nlon,nlat), intent(inout) :: dem , wk1
    integer, dimension(nlon,nlat), intent(inout) :: fmap
    real :: h1 , h2 , slope
    integer :: i1 , i2 , ibet , idir , idist , ii , iskip , istep , j1 , j2 , &
               jbet , jdir , jj , jskip , jstep , mindist , ncal , np , nstep
    logical :: plot
    plot = .false.
    ncal = 1
    if ( ncal == 10 .and. plot ) then
      plot = .false.
      write(output_unit,'(12x,a)') &
          'Too many calls to angioplasty. No plot produced.'
    end if
    np = angionp
    i1 = i - np
    if ( i1<1 ) i1 = 1
    j1 = j - np
    if ( j1<1 ) j1 = 1
    i2 = i + np
    if ( i2>nlon ) i2 = nlon
    j2 = j + np
    if ( j2>nlat ) j2 = nlat
    ibet = 0
    jbet = 0
    mindist = 1000
    do ii = i1 , i2
      do jj = j1 , j2
        idist = iabs(i-ii) + iabs(j-jj)
        if ( dem(ii,jj)<dem(i,j) .and. &
             nint(wk1(ii,jj))/=-5 .and. idist<mindist ) then
          ibet = ii
          jbet = jj
          mindist = iabs(i-ii) + iabs(j-jj)
        end if
      end do
    end do
    if ( ibet==0 ) then
      write(output_unit,'(15x,a)') 'Now Applying less restrictive algorithm'
      do ii = i1 , i2
        do jj = j1 , j2
          idist = iabs(i-ii) + iabs(j-jj)
          if ( nint(wk1(ii,jj))/=-5 .and. idist<mindist ) then
            ibet = ii
            jbet = jj
            mindist = iabs(i-ii) + iabs(j-jj)
          end if
        end do
      end do
      if ( ibet/=0 ) dem(i,j) = dem(ibet,jbet) + mindist
    end if
    if ( ibet==0 ) then
      write(output_unit,'(12x,a)') 'Cannot solve singularity.'
      return
    end if
    if ( ibet<i ) then
      istep = -1
      idir = 8
      iskip = 1
    else
      istep = 1
      idir = 4
      iskip = -1
    end if
    if ( jbet<j ) then
      jstep = -1
      jdir = 6
      jskip = 1
    else
      jstep = 1
      jdir = 2
      jskip = -1
    end if
    nstep = 0
    h1 = dem(i,j)
    h2 = dem(ibet,jbet)
    slope = (dem(i,j)-dem(ibet,jbet))/mindist
    if ( ibet<i ) then
      fmap(i,j) = 8
    else
      fmap(i,j) = 4
    end if
    wk1(i,j) = dem(i,j)
    do ii = i + istep , ibet + iskip , istep
      nstep = nstep + 1
      fmap(ii,j) = idir
      dem(ii,j) = h1 - nstep*slope
      wk1(ii,j) = dem(ii,j)
    end do
    do jj = j , jbet + jskip , jstep
      nstep = nstep + 1
      fmap(ibet,jj) = jdir
      dem(ibet,jj) = h1 - nstep*slope
      wk1(ibet,jj) = dem(ibet,jj)
    end do
  end subroutine angioplasty

  subroutine findthemouths(w,cut,ilog,nsec,isec,jsec)
    implicit none
    real, dimension(nlon,nlat), intent(in) :: w
    real, intent(in) :: cut
    integer, intent(in) :: ilog
    integer, intent(out) :: nsec
    integer, dimension(nlon*nlat), intent(inout) :: isec , jsec
    integer :: i , idir , j
    if ( ilog == 1 ) then
      write(output_unit,'(21x,a)') 'Rivers mouths finding...'
      write(output_unit,'(8x,68(''-''))')
      write(output_unit,'(10x,a)') &
          '#    I   J     Lat      Lon   Area(Km2)   Name'
      write(output_unit,'(8x,68(''-''))')
    end if
    nsec = 0
    do j = nlat - 1 , 2 , -1
      do i = 2 , nlon - 1
        idir = fmap(i,j)
        if ( idir>=1 .and. idir<=8 ) then
          if ( luc(i+ir(idir),j+jr(idir))==ocean .and. w(i,j)>cut ) then
            nsec = nsec + 1
            isec(nsec) = i
            jsec(nsec) = j
            if ( ilog==1 ) &
              write(output_unit,'(9x,i2,2x,2i4,2f9.4,f11.1,3x,a25)') nsec ,  &
                i , j , lat(i,j) , lon(i,j) , w(i,j)
          end if
        end if
      end do
    end do
  end subroutine findthemouths

  subroutine establishfowdir
    implicit none
    integer :: i, j, k
    integer , dimension(1) :: iloc
    real , dimension(0:8) :: rk

    rk(0) = 0.0
    do j = 2 , nlat - 1
      do i = 2 , nlon - 1
        if ( luc(i,j) /= ocean ) then  ! land use chk
          do k = 1 , 8
            if ( dem(i+ir(k),j+jr(k)) < dem(i,j) ) then
              rk(k) = (dem(i,j)-dem(i+ir(k),j+jr(k))) / &
                distance(lat(i,j),lon(i,j), &
                         lat(i+ir(k),j+jr(k)),lon(i+ir(k),j+jr(k)))
            else
              rk(k) = -1.0
            end if
          end do
          iloc = maxloc(rk(0:8))
          fmap(i,j) = iloc(1) - 1
        end if
      end do
    end do
  end subroutine establishfowdir

  subroutine noflowcounter(cam,n_noflow)
    implicit none
    integer, dimension(nlon,nlat), intent(out) :: cam
    integer, intent(out) :: n_noflow
    integer :: i , j
    n_noflow = 0
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
        if ( luc(i,j)==ocean .or. luc(i,j)==iwater ) then
          cam(i,j) = 0
        else if ( fmap(i,j)==0 ) then
          n_noflow = n_noflow + 1
          cam(i,j) = 8
        else
          cam(i,j) = 0
        end if
      end do
    end do
  end subroutine noflowcounter

  subroutine dircorrflow(icl)
    implicit none
    integer, dimension(nlon,nlat), intent(inout) :: icl
    integer :: i , j , k
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
        if ( icl(i,j) == 8 ) then
          do k = 1 , 4
            if ( fmap(i+ir(k),j+jr(k)) == fmap(i+ir(k+4),j+jr(k+4)) ) then
              fmap(i,j) = fmap(i+ir(k),j+jr(k))
              icl(i,j) = 0
            end if
          end do
        end if
        if ( icl(i,j) == 8 ) then
          do k = 1 , 8
            if ( fmap(i+ir(k),j+jr(k)) == k ) then
              fmap(i,j) = k
              icl(i,j) = 0
            end if
          end do
        end if
      end do
    end do
  end subroutine dircorrflow

  subroutine flowcheck(water,wk1)
    implicit none
    real, dimension(nlon,nlat), intent(inout) :: water, wk1
    real, dimension(nlon,nlat) :: wk
    integer :: i, idir, ii, isa, j, jsa, niter
    real :: xmax

    wk1(:,:) = 0.0
    write(output_unit,'(18x,a)') 'Flow Check Module Called.'
    niter = (nlon+nlat)*2
    do ii = 1 , niter
      wk(:,:) = 0.0
      do j = 2 , nlat - 1
        do i = 2 , nlon - 1
          idir = fmap(i,j)
          if ( idir > 0 ) then
            wk(i+ir(idir),j+jr(idir)) = wk(i+ir(idir),j+jr(idir)) + water(i,j)
            water(i,j) = 0
          end if
        end do
      end do
      do j = 1 , nlat
        do i = 1 , nlon
          water(i,j) = water(i,j) + wk(i,j)
          wk1(i,j) = wk1(i,j) + wk(i,j)
        end do
      end do
    end do
    wk(1:nlon,1:nlat) = wk1(1:nlon,1:nlat)
    xmax = 0.0
    write(output_unit,'(18x,a)') 'Flow Check Module First part ended.'
    do j = 1 , nlat
      do i = 1 , nlon
        if ( wk1(i,j) > xmax ) then
          isa = i
          jsa = j
          xmax = wk1(i,j)
        end if
        if ( luc(i,j) /= ocean .and. luc(i,j) /= iwater ) then
          if ( wk1(i,j) > float((nlon+nlat)*ifactor_res) ) then
            write(output_unit,'(20x,a,2i5,f20.5)') &
                'Discarded: ' , i , j , wk1(i,j)
            noflow(i,j) = 8
          end if
        end if
      end do
    end do
    write(output_unit,'(21x,a,f20.5)') &
        'Maximum value of "Rolling Stones" is: ' , xmax
  end subroutine flowcheck

  subroutine runoffspeed
    implicit none
    real :: alfamin , alfamax , delta , xgamma , mann , tresh , vmax
    real :: enne , hrad
    real, dimension(nlon,nlat) :: wrk2 , wk
    integer :: i , idir , intstep , j , land
    !!!
    wrk2 = 0
    alfa = 0.0
    xgamma = 0.33
    delta = cpar8          ! Param. for land/channel flow
    tresh = cpar6
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
        idir = fmap(i,j)
        land = luc(i,j)
        mann = manning(land)
        if ( idir>=1 .and. idir<=8 .and. land/=ocean .and. land>0 ) then
          if ( land>100 .or. land<=0 ) then
            write(error_unit,'(10x,a,i5)') &
                'Wrong value for landuse code: ' , land
            stop 'flux error inside runoffspeed'
          end if
          dx(i,j) = distance(lat(i,j),lon(i,j), &
                             lat(i+ir(idir),j+jr(idir)), &
                             lon(i+ir(idir),j+jr(idir)))
          !enne = mann/(1.+(delta-1.)*(1.+max((drai(i,j)-tresh),0.0)/tresh))
          enne = mann/delta
          hrad = cpar2 + cpar3*(max(drai(i,j),tresh)**xgamma)
          alfa(i,j) = ((hrad**0.6666)*(accl(i,j)**0.5))/enne
        end if
      end do
    end do

    vmax = 0.0
    do i = 1 , nlon
      do j = 1 , nlat
        if ( alfa(i,j)>vmax .and. drai(i,j)>100.0 ) vmax = alfa(i,j)
      end do
    end do
    intstep = nint(vmax*6.0)
    if ( intstep<10 ) intstep = 10
    alfamin = 0.10
    do i = 1 , nlon
      do j = 1 , nlat
        if ( luc(i,j)==ocean ) then
          wrk2(i,j) = -10.0
        else if ( alfa(i,j)<alfamin ) then
          wrk2(i,j) = dx(i,j)/alfamin
        else
          wrk2(i,j) = dx(i,j)/alfa(i,j)
        end if
      end do
    end do
    call runofftime(wrk2,fmap,runt,nlon,nlat)
    wk = 1
    call rollingstones2(wk,fmap,wrk2,nlon,nlat)
    call rollingstones2(runt,fmap,wk,nlon,nlat)
    do i = 1 , nlon
      do j = 1 , nlat
        runt(i,j) = wk(i,j)/wrk2(i,j) - runt(i,j)
      end do
    end do
  end subroutine runoffspeed

  subroutine runofftime(rntime,fmap,rout,nlon,nlat)
    implicit none
    integer, intent(in) :: nlat , nlon
    integer, dimension(nlon,nlat), intent(in) :: fmap
    real, dimension(nlon,nlat), intent(in) :: rntime
    real, dimension(nlon,nlat), intent(inout) :: rout
    logical :: exceed
    integer :: i , idir , ist , j , jst , n
    integer, save :: ncycle = 10000

    exceed = .false.
    do i = 1 , nlon
      do j = 1 , nlat
        rout(i,j) = rntime(i,j)/2
        ist = i
        jst = j
        idir = fmap(ist,jst)
        do n = 1 , ncycle
          if ( i>1 .and. i<nlon .and. j>1 .and. j<nlat .and. &
               idir>=1 .and. idir<=8 ) then
            call chymmoveahead(idir,ist,jst)
            rout(i,j) = rout(i,j) + rntime(ist,jst)
            idir = fmap(ist,jst)
          else
            rout(i,j) = rout(i,j)/3600
            exit
          end if
        end do
        if ( n>=ncycle ) exceed = .true.
      end do
    end do
    if ( exceed ) &
      write(output_unit,'(/,10x,a/)') 'CHyM severe warning: runofftime'// &
               ' maximum number of cycles exceeded'
  end subroutine runofftime

  subroutine chymmoveahead(idir,i,j)
    implicit none
    integer, intent(inout) :: i , j
    integer, intent(in) :: idir
    integer :: di , dj
    call neighborhood(idir,di,dj)
    i = i + di
    j = j + dj
  end subroutine chymmoveahead

  subroutine neighborhood(id,i,j)
    implicit none
    integer, intent(in) :: id
    integer, intent(out) :: i , j

    integer, dimension(0:48), parameter :: ir = &
        [ 0,-1, 0, 1, 1, 1, 0,-1,-1,-2,-1, 0, 1, 2, 2, 2, 2, &
          2, 1, 0,-1,-2,-2,-2,-2,-3,-2,-1, 0, 1, 2, 3, 3, 3, &
          3, 3, 3, 3, 2, 1, 0,-1,-2,-3,-3,-3,-3,-3,-3 ]
    integer, dimension(0:48), parameter :: jr = &
        [ 0, 1, 1, 1, 0,-1,-1,-1, 0, 2, 2, 2, 2, 2, 1, 0,-1, &
         -2,-2,-2,-2,-2,-1, 0, 1, 3, 3, 3, 3, 3, 3, 3, 2, 1, &
          0,-1,-2,-3,-3,-3,-3,-3,-3,-3,-2,-1, 0, 1, 2 ]
    if ( id < 0 .or. id > 48 ) then
      write(error_unit,'(10x,a,i4,a)') &
         'Invalid id inside subroutine neighborhood: ', id , '. Exiting...'
      call exit(0)
    end if
    i = ir(id)
    j = jr(id)
  end subroutine neighborhood

  subroutine rollingstones2(rntime,fmap,cout,nlon,nlat)
    implicit none
    integer, intent(in) :: nlat , nlon
    integer, dimension(nlon,nlat), intent(in) :: fmap
    real, dimension(nlon,nlat), intent(in) :: rntime
    real, dimension(nlon,nlat), intent(inout) :: cout
    integer :: di , dj , i , idir , ist , j , jst , n
    logical :: exceed
    integer, save :: ncycle = 10000

    exceed = .false.
    cout = rntime
    do i = 1 , nlon
      do j = 1 , nlat
        ist = i
        jst = j
        idir = fmap(ist,jst)
        do n = 1 , ncycle
          if ( i > 1 .and. i < nlon .and. &
               j > 1 .and. j < nlat .and. &
               idir >= 1 .and. idir < 8 ) then
            call neighborhood(idir,di,dj)
            ist = ist + di
            jst = jst + dj
            cout(ist,jst) = cout(ist,jst) + rntime(i,j)
            idir = fmap(ist,jst)
          else
            exit
          end if
        end do
        if ( n>=ncycle ) exceed = .true.
      end do
    end do
    if ( exceed ) write (error_unit,'(/,10x,a/)')    &
        'CHyM severe warning: rollingstones2'//      &
        ' maximum number of cycles exceeded'
  end subroutine rollingstones2

  subroutine basinpaint(dem,fmap,luse,work,nlon,nlat,nisa,isa,jsa,wk1,flg)
    implicit none
    integer, intent(in) :: nlat, nlon, nisa, flg
    integer, dimension(nlon*nlat), intent(in) :: isa , jsa
    real, dimension(nlon,nlat), intent(in) :: dem
    integer, dimension(nlon,nlat), intent(in) :: fmap , luse
    real, dimension(nlon,nlat), intent(out) :: wk1
    real, dimension(nlon,nlat), intent(out) :: work
    integer :: i , idir , ii , iind , j , jj, iriv
    logical :: notyet

    if ( flg == 1 .or. flg == 0 ) then
      wk1 = -10
      work = 0.0
    end if
    do iriv = 1, nisa
      do j = 2 , nlat - 1
        do i = 2 , nlon - 1
          ii = i
          jj = j
          iind = 0
          notyet = .true.
          do while ( iind <= 2*(nlon+nlat) .and. notyet )
            iind = iind + 1
            idir = fmap(ii,jj)
            if ( idir > 0 ) then
              ii = ii + ir(idir)
              jj = jj + jr(idir)
              if ( ii == isa(iriv) .and. jj == jsa(iriv) ) then
                notyet = .false.
                if ( flg == 0 ) then
                  wk1(i,j) = dem(i,j)
                else if ( flg == -1 ) then
                  wk1(i,j) = -999
                else
                  wk1(i,j) = flg
                end if
              else if ( ii == 1 .or. ii == nlon .or. &
                        jj == 1 .or. jj == nlat ) then
                notyet = .false.
              end if
            else if ( luse(ii,jj) == ocean ) then
              notyet = .false.
            else                               ! noflow
              wk1(i,j) = -5.0
              if ( flg /= 0 ) wk1(i,j) = -10
              work(ii,jj) = work(ii,jj) + 1
              notyet = .false.
            end if
          end do
        end do
      end do
    end do
  end subroutine basinpaint

end module mod_crtstat
