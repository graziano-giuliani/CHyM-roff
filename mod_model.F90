module mod_model

  use mod_param
  use mod_mpimess
  use mod_varandtypes

  implicit none
  private

  public :: chymmodel

  contains

    subroutine chymmodel(chym_runoff,imon,iday)
      implicit none
      real, intent(in) :: chym_runoff(:,:)
      integer, intent(in) :: imon, iday
      real :: tmp
      logical :: first
      data first /.true./
      save first

      integer :: i,j,ii,jj,idir
      integer,save :: i1,i2,j1,j2
      real :: rainload,dm
      real, parameter :: maxfrac = 0.95

      if (first) then
         first=.false.
         call mpi_barrier(mycomm,mpierr)
         i1 = ide1gb+1
         i2 = ide2gb-1
         j1 = jde1gb+1
         j2 = jde2gb-1
         if (ma%has_bdytop) i2 = ide2gb
         if (ma%has_bdybottom) i1 = ide1gb
         if (ma%has_bdyleft) j1 = jde1gb
         if (ma%has_bdyright) j2 = jde2gb
#ifdef NILE
         idamietta = -1
         jdamietta = -1
         irosetta = -1
         jrosetta = -1
         do i = i1, i2
           do j = j1, j2
             if ( is_inbox(lat_damietta,lon_damietta, &
                    corner_lat(:,j,i),corner_lon(:,j,i)) ) then
               call find_nearest_land(j,i,jj,ii)
               idamietta = ii
               jdamietta = jj
               print *, 'Damietta is at ',jj,ii
             end if
             if ( is_inbox(lat_rosetta,lon_rosetta, &
                    corner_lat(:,j,i),corner_lon(:,j,i)) ) then
               call find_nearest_land(j,i,jj,ii)
               irosetta = ii
               jrosetta = jj
               print *, 'Rosetta is at ',jj,ii
             end if
           end do
         end do
#endif
      endif
      call exchange_bdy(h2o_sub(jde1gb:jde2gb,ide1gb:ide2gb),2,jde1,     &
         jde2,ide1,ide2,jde1gb,jde2gb,ide1gb,ide2gb)
      call exchange_bdy(bwet_sub(jde1gb:jde2gb,ide1gb:ide2gb),2,jde1,    &
         jde2,ide1,ide2,jde1gb,jde2gb,ide1gb,ide2gb)
      call exchange_bdy(port_sub(jde1gb:jde2gb,ide1gb:ide2gb),2,jde1,    &
         jde2,ide1,ide2,jde1gb,jde2gb,ide1gb,ide2gb)
      do i = i1, i2
        do j = j1, j2
          idir = fmap(j,i)
          if (luse(j,i) /= ocean .and. idir >= 1 .and. idir <= 8) then
            dm = port_sub(j,i)*deltat
            if ( dm > maxfrac * h2o_sub(j,i) ) then
              dm = maxfrac * h2o_sub(j,i)
            end if
            wkm1_sub(j,i) = wkm1_sub(j,i) - dm
            wkm1_sub(j+ir(idir),i+jr(idir))=                            &
                 wkm1_sub(j+ir(idir),i+jr(idir)) + dm
          end if
        end do
      end do
      call exchange_bdy(wkm1_sub(jde1gb:jde2gb,ide1gb:ide2gb),2,jde1,   &
           jde2,ide1,ide2,jde1gb,jde2gb,ide1gb,ide2gb)
      do i = i1, i2
        do j = j1, j2
          idir = fmap(j,i)
          if (luse(j,i) /= ocean .and. idir >= 1 .and. idir <= 8) then
             ! m3 of water recharge in the  grid cell
             rainload = chym_area(j,i)*1.0e+06*(chym_runoff(j,i))*deltat
             h2o_sub(j,i) = h2o_sub(j,i) + wkm1_sub(j,i) + &
                            convfac*rainload
             h2o_sub(j,i) = max(h2o_sub(j,i),0.0)
             if ( farm(j,i) .and. chym_lat(j,i) > 0.0 ) then
               ! Some water is not transmitted for irrigation
               bwet_sub(j,i) = (1.-irmonfac(imon)/deltat) * &
                   h2o_sub(j,i) / chym_dx(j,i)
               ! Evaporative loss of the water for irrigation
               h2o_sub(j,i) = (1.-irloss*irmonfac(imon)/deltat) * h2o_sub(j,i)
             else
               bwet_sub(j,i) = h2o_sub(j,i) / chym_dx(j,i)
             end if
             port_sub(j,i) = alfa(j,i) * bwet_sub(j,i)
          end if
        end do
      end do
#ifdef NILE
      tmp = 0.5 * mval(nile_fresh_flux,imon,iday)
      if ( idamietta > 0 .and. jdamietta > 0 ) then
        port_sub(jdamietta,idamietta) = tmp
      end if
      if ( irosetta > 0 .and. jrosetta > 0 ) then
        port_sub(jrosetta,irosetta) = tmp
      end if
#endif
    end subroutine chymmodel

#ifdef NILE

  real function mval(series,imon,iday)
    implicit none
    integer , intent(in) :: imon, iday
    real, dimension(12) :: series
    integer :: imonp, imonn
    real :: w1, w2, w3
    real :: fm, hm
    real, dimension(12), parameter :: mpd = &
     [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
    imonp = imon - 1
    imonn = imon + 1
    if ( imonp < 1 ) imonp = 12
    if ( imonn > 12 ) imonn = 1
    fm = mpd(imon)
    hm = fm/2.0
    if ( iday < hm ) then
      w1 = 0.5 - (iday-1+0.5)/fm
      w2 = 1.0 - w1
      w3 = 0.0
    else
      w1 = 0.0
      w2 = 1.0 - (iday-hm+0.5)/fm
      w3 = 1.0 - w2
    end if
    mval = ( w1 * series(imonp) + &
             w2 * series(imon)  + &
             w3 * series(imonn) )
  end function mval
#endif

end module mod_model
