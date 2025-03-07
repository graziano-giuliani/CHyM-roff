module mod_model

  use mod_param
  use mod_mpimess
  use mod_varandtypes

  implicit none
  private

  public :: chymmodel

  contains

    subroutine chymmodel(istep, chym_runoff)

      implicit none

      integer, intent(in) :: istep
      real, intent(in) :: chym_runoff(:,:)
      logical :: first
      data first /.true./
      save first

      integer :: i,j,idir
      integer,save :: i1,i2,j1,j2
      real :: rainload,dm

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
          if (luse(j,i) /= mare .and. idir >= 1 .and. idir <= 8) then
            dm = port_sub(j,i)*deltat
            if ( dm > h2o_sub(j,i) ) then
               dm = h2o_sub(j,i)
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
          if (luse(j,i) /= mare .and. idir >= 1 .and. idir <= 8) then
             ! m3 of water recharge in the  grid cell
             rainload = chym_area(j,i)*1.0e+06*(chym_runoff(j,i))*deltat
             h2o_sub(j,i) = h2o_sub(j,i) + wkm1_sub(j,i) + &
                            convfac*rainload
             bwet_sub(j,i) = h2o_sub(j,i) / chym_dx(j,i)
             port_sub(j,i) = alfa(j,i) * bwet_sub(j,i)
          end if
        end do
      end do
    end subroutine chymmodel

end module mod_model
