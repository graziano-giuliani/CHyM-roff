!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP CHyM.
!
!    ICTP CHyM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by
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

module mod_iface

  use mod_param
  use mod_io
  use mod_model
  use mod_mpimess
  use mod_varandtypes
  use mod_time

  implicit none
  private

  public chym_init
  public chym_run
  public chym_close

  contains

    subroutine chym_init()
      implicit none
      integer :: i, i1, j1
      integer :: displacem
      integer :: is,js,ks
      character(len=256) :: namelistfile
      character(len=10) :: tsdate

      call getarg(1, namelistfile)

      call read_config(trim(namelistfile))

      call read_init()

      if (.not. allocated(chym_runoff)) allocate(chym_runoff(nlc,       &
       nbc))

      i1 = nlc
      j1 = (nbc/nproc)+4
      iy = nbc
      jx = nlc

      call set_nproc
      call setup_model_indexes
      if (allocated(port_sub)) deallocate(port_sub)
      allocate (port_sub(jde1gb:jde2gb,ide1gb:ide2gb))
      if (allocated(wkm1_sub)) deallocate(wkm1_sub)
      allocate (wkm1_sub(jde1gb:jde2gb,ide1gb:ide2gb))
      if (allocated(h2o_sub)) deallocate(h2o_sub)
      allocate (h2o_sub(jde1gb:jde2gb,ide1gb:ide2gb))
      if (allocated(bwet_sub)) deallocate(bwet_sub)
      allocate (bwet_sub(jde1gb:jde2gb,ide1gb:ide2gb))
      port_sub = 0.
      wkm1_sub = 0.
      h2o_sub = 0.
      bwet_sub = 0.
      allocate(cartesian_np(nproc), cartesian_dis(nproc))
      allocate(cartesian_npoint(nproc), cartesian_displ(nproc))
      if (myid == 0) then
        allocate(ide1p(nproc), ide2p(nproc), jde1p(nproc), jde2p(nproc))
        allocate(iypp(nproc) , jxpp(nproc))
        allocate(ide1gbp(nproc), ide2gbp(nproc), jde1gbp(nproc),        &
           jde2gbp(nproc))
        allocate(iypgbp(nproc) , jxpgbp(nproc))
      else
        allocate(ide1p(1), ide2p(1), jde1p(1), jde2p(1))
        allocate(iypp(1) , jxpp(1))
        allocate(ide1gbp(1), ide2gbp(1), jde1gbp(1),        &
           jde2gbp(1))
        allocate(iypgbp(1) , jxpgbp(1))
      end if
      call mpi_reduce(nngb, nngbp, 1, mpi_integer, mpi_sum, 0, mycomm,  &
        mpierr)
      if (myid  == 0) then
        allocate(port1d(nngbp),bwet1d(nngbp),h2o1d(nngbp),wkm11d(nngbp))
        port1d=0.
        wkm11d=0.
        bwet1d=0.
        h2o1d=0.
      end if
      call mpi_gather(iyp,1,mpi_integer,                                &
                       iypp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jxp,1,mpi_integer,                                &
                       jxpp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(ide1,1,mpi_integer,                               &
                       ide1p,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(ide2,1,mpi_integer,                               &
                       ide2p,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jde1,1,mpi_integer,                               &
                       jde1p,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jde2,1,mpi_integer,                               &
                       jde2p,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(ide1gb,1,mpi_integer,                             &
                       ide1gbp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(ide2gb,1,mpi_integer,                             &
                       ide2gbp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jde1gb,1,mpi_integer,                             &
                       jde1gbp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jde2gb,1,mpi_integer,                             &
                       jde2gbp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jde2gb,1,mpi_integer,                             &
                       jde2gbp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(iypgb,1,mpi_integer,                              &
                       iypgbp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jxpgb,1,mpi_integer,                              &
                       jxpgbp,1,mpi_integer,0,mycomm,mpierr)
      allocate(h2osub1d(iypgb*jxpgb))
      allocate(portsub1d(iypgb*jxpgb))
      allocate(bwetsub1d(iypgb*jxpgb))
      allocate(wkm1sub1d(iypgb*jxpgb))
      allocate(portsub1d_n(iyp*jxp))
      allocate(h2osub1d_n(iyp*jxp))
      allocate(bwetsub1d_n(iyp*jxp))
      allocate(wkm1sub1d_n(iyp*jxp))
      allocate(port1d_n(jx*iy))
      allocate(h2o1d_n(jx*iy))
      allocate(bwet1d_n(jx*iy))
      allocate(wkm11d_n(jx*iy))
      h2osub1d=0.
      portsub1d=0.
      bwetsub1d=0.
      wkm1sub1d=0.
      portsub1d_n=0.
      h2osub1d_n=0.
      bwetsub1d_n=0.
      wkm1sub1d_n=0.
      port1d_n=0.
      h2o1d_n=0.
      bwet1d_n=0.
      wkm11d_n=0.

      if (myid == 0) then
         cartesian_dis(1) = 0
         cartesian_displ(1) = 0
         do i = 1,nproc
           cartesian_np(i) = iypp(i)*jxpp(i)
           cartesian_npoint(i) = iypgbp(i)*jxpgbp(i)
           if (i == nproc) exit
           cartesian_dis(i+1) = cartesian_dis(i) + cartesian_np(i)
           cartesian_displ(i+1) = cartesian_displ(i) +                  &
             cartesian_npoint(i)
         end do
      end if
      call mpi_bcast(cartesian_np,nproc,mpi_integer,0,mycomm,mpierr)
      call mpi_bcast(cartesian_dis,nproc,mpi_integer,0,mycomm,mpierr)
      if (isread.eq.1) then
        if (myid == 0) then
          displacem = 0
          cartesian_dis(1) = 0
          cartesian_displ(1) = 0
          do i = 1,nproc
            call mypack_real_grid(h2o,h2o1d,ide1gbp(i),ide2gbp(i),      &
                             jde1gbp(i),jde2gbp(i),displacem)
            call mypack_real_grid(port,port1d,ide1gbp(i),ide2gbp(i),    &
                             jde1gbp(i),jde2gbp(i),displacem)
            displacem = displacem + iypgbp(i)*jxpgbp(i)
            cartesian_np(i) = iypp(i)*jxpp(i)
            cartesian_npoint(i) = iypgbp(i)*jxpgbp(i)
            if (i == nproc) exit
            cartesian_dis(i+1) = cartesian_dis(i) + cartesian_np(i)
            cartesian_displ(i+1) = cartesian_displ(i) +                 &
                cartesian_npoint(i)
          end do
        end if
        call mpi_scatterv(h2o1d, cartesian_npoint, cartesian_displ,     &
           mpi_real,h2osub1d,iypgb*jxpgb, mpi_real, 0, mycomm, mpierr)
        call mpi_scatterv(port1d, cartesian_npoint, cartesian_displ,    &
           mpi_real,portsub1d,iypgb*jxpgb, mpi_real, 0, mycomm, mpierr)
        call myunpack_real_grid(h2osub1d,h2o_sub,ide1gb,ide2gb,jde1gb,   &
            jde2gb)
        call myunpack_real_grid(portsub1d,port_sub,ide1gb,ide2gb,jde1gb, &
            jde2gb)
      endif
      call MPI_BARRIER(mycomm,mpierr)

      if (myid == 0) then
        time=sdate
        write(tsdate,'(i0.10)') sdate
        call gmafromindex(time,hour,day,month,year)
        write(filename,'(a,a)') trim(tdnsim)//'_',trim(tsdate)// &
            '.nc'
        write(filenamerst,'(a,a)') trim(tdnsim)//'_',trim(tsdate)// &
            '_rst.nc'
        write(filenameqmax,'(a,a)') trim(tdnsim)//'_',trim(tsdate)// &
            '_qmax.nc'
        print*,"Create output file"
        call createfile(trim(filename),time)
        print*,"End creation output file"

      end if
      do js=2,nbc-1
         do is=2,nlc-1
            do ks=1,4
            if (fmap(is,js)==ks .and. fmap(is+ir(ks),js+jr(ks))==ks+4) then
              fmap(is,js) = 0
              fmap(is+ir(ks),js+jr(ks)) = 0
            end if
            end do
         end do
      end do
    end subroutine chym_init


    subroutine chym_run
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: t,i,j
      integer :: iv,k
#ifdef RUNOFF
      integer :: ii, jj, idir
#endif
      character(len=80) :: string, now

      hourstep = 0
      chym_runoff = 0.0

      do while (time <= edate)

        hourstep = hourstep + dstep
        if (myid == 0 ) then

          call createnc1(hourstep)
          write (string,'(15x,a,i8,a,i10)')                               &
              'CHyM integration step number ',&
              hourstep,': ',time
          write (6,'(a)') string(1:len_trim(string))
          call gmafromindex(time,hour,day,month,year)
          call dataorafromday (hour,day,month,year,now)
        end if
        call read_runoff(hourstep/dstep+inirun)
        do j=1,nbc
          do i=1,nlc
            chym_runoff(i,j) = chym_runoff(i,j) * 0.001
            if (chym_runoff(i,j) < 0.0 ) then
              chym_runoff(i,j) = 0.0
            end if
            ! Assume this is missing value
            if ( chym_runoff(i,j) > 1.0 ) then
              chym_runoff(i,j) = 0.0
            end if
          enddo
        enddo
        if (myid == 0 ) then
          do t=1,dstep
            time=increasetime(time)
          end do
        end if
        call mpi_bcast(time,nproc,mpi_integer,0,mycomm,mpierr)
        call gmafromindex(time,hour,day,month,year)
        call dataorafromday (hour,day,month,year,now)

        deltat = (3600.0*dstep)/real(step)
        do i = 1, step
          wkm1_sub=0.0
          call chymmodel(chym_runoff,month,day)
        enddo
        iv = 1
        do i=ide1,ide2
          do j=jde1,jde2
            portsub1d_n(iv)=port_sub(j,i)
            h2osub1d_n(iv)=h2o_sub(j,i)
            iv = iv+1
          end do
        end do
        iv = cartesian_np(myid+1)
        call mpi_gatherv(portsub1d_n,iv,MPI_REAL,port1d_n,  &
                cartesian_np, cartesian_dis,                &
                MPI_REAL,0,cartesian_communicator, mpierr)
        call mpi_gatherv(h2osub1d_n,iv,MPI_REAL,h2o1d_n,    &
                cartesian_np, cartesian_dis,                &
                MPI_REAL,0,cartesian_communicator, mpierr)
        if (myid == 0) then
          iv = 1
          do k=1,nproc
            do i=ide1p(k),ide2p(k)
              do j=jde1p(k),jde2p(k)
                port(j,i) = port1d_n(iv)
                h2o(j,i) = h2o1d_n(iv)
                iv = iv + 1
              end do
            end do
          end do
        end if
!
!-----------------------------------------------------------------------
!     Write to restart file
!-----------------------------------------------------------------------
!
        if (myid == 0 ) then
          if (iorstfreq /= 0) then
            call createnc2(hourstep,iorstfreq)
          end if
        end if
!-----------------------------------------------------------------------
!     Write output to file
!-----------------------------------------------------------------------
!
        if (myid == 0 ) then
          port_out = int(port)
          call add_timestep(chymout%ncid,chymout%varid(3),iostep)
          call write_dynvar(chymout%ncid,chymout%varid(4),port_out, &
                            iostep)
#ifdef RUNOFF
          roff_out = port_out/(chym_area*1.0e6)
          do j=2,nbc-1
            do i=2,nlc-1
              if ( chym_lsm(i,j) > 0.0 ) then
                idir = fmap(i,j)
                if ( idir == 0 ) then
                  call find_nearest_ocean(i,j,ii,jj)
                  roff_out(ii,jj) = roff_out(i,j)
                else
                  roff_out(i+ir(idir),j+jr(idir)) = roff_out(i,j)
                end if
              end if
            end do
          end do
          call write_dynvar(chymout%ncid,chymout%varid(5),roff_out, &
                            iostep)
#endif
          do j=2,nbc-1
            do i=2,nlc-1
              if (port_qmaxs(i,j) < port_out(i,j)) then
                port_qmaxs(i,j) = port_out(i,j)
              end if
            end do
          end do
          call createnc3(hourstep)
        end if

      end do ! edate reached

    end subroutine chym_run

    subroutine chym_close
        implicit none
        if (myid == 0) then
          write (6,'(/12x,a)') 'Closing all files'
          call closefile(chymout%ncid)
          write (6,'(/12x,a)') 'Simulation completed!'
        end if
        call mpi_barrier(mycomm,mpierr)
        call mpi_finalize(mpierr)
    end subroutine

end module mod_iface
