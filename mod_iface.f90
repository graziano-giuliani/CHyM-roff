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
      integer :: i1,j1
      integer :: i,j
      integer :: displacem
      integer :: is,js,ks
      character(len=256) :: namelistfile
      
      call getarg(1, namelistfile)
      call read_config(trim(namelistfile))

      call read_init()
      call runoffspeed

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
      end if
      call mpi_reduce(nngb, nngbp, 1, mpi_integer, mpi_sum, 0, mycomm,  &
        mpierr)
      if (myid  == 0) then
        allocate(port1d(nngbp),bwet1d(nngbp),h2o1d(nngbp))
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
        write(filename,'(a,a)') trim(sim_name)//'_',trim(tsdate)// &
            '.nc'
        write(filenamerst,'(a,a)') trim(sim_name)//'_',trim(tsdate)// &
            '_rst.nc'
        write(filenameqmax,'(a,a)') trim(sim_name)//'_',trim(tsdate)// &
            '_qmax.nc'
        print*,"Create output file"
        call createfile(trim(filename),time)
        print*,"End creation output file"

        if (iswrit /= 0) then
        end if

      end if
      do js=1,nbc
         do is=1,nlc
            do ks=1,4
            if (fmap(is,js)==ks .and. fmap(is+ir(ks),js+jr(ks))==ks+4) then
              fmap(is,js) = 0
              fmap(is+ir(ks),js+jr(ks)) = 0
            end if     
            end do
         end do
      end do
      end subroutine chym_init


      subroutine chym_run(istart, iend)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      integer, intent(in) :: istart
      integer, intent(in) :: iend
!
!-----------------------------------------------------------------------
!     Local variable declarations  
!-----------------------------------------------------------------------
!

      integer :: t,i,j, istep, icount
      integer :: iv,k
      character(len=80) :: string, now
      if ( myid == 0 ) then
        if (isread /= 0 .and. iswrit /= 0) then
          icount = mod(istart, iswrit)
        else
          icount = 0
        end if
      end if
      hourstep = 0
      do while (time <= edate)
        hourstep = hourstep + dstep
        if (myid == 0 ) then
          
          call createnc1(hourstep)
          write (string,'(15x,a,i5,a,i10)')                               &
              'CHyM integration step number ',&
              hourstep,': ',time
          write (6,'(a)') string(1:len_trim(string))           
          call gmafromindex(time,hour,day,month,year)
          call dataorafromday (hour,day,month,year,now)
        end if
        if (istep == 1) then
          chym_runoff = 0.0
        end if
        call read_runoff(hourstep/dstep+inirun)
        do j=1,nbc
          do i=1,nlc
            chym_runoff(i,j) = chym_runoff(i,j) * 0.001
            if (chym_runoff(i,j)>1 .or. chym_runoff(i,j)<0 )              &
                chym_runoff(i,j) = 0
          enddo
        enddo
        if (myid == 0 ) then
          do t=1,dstep
            time=increasetime(time)
          end do
          call gmafromindex(time,hour,day,month,year)
          call dataorafromday (hour,day,month,year,now)
        end if
      call mpi_bcast(time,nproc,mpi_integer,0,mycomm,mpierr)

      deltat = 86400/step
      do i=1,step
        wkm1_sub=0.0
        call chymmodel(istep,chym_runoff)
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
      call mpi_gatherv(portsub1d_n,iv,MPI_REAL,port1d_n,                &
         cartesian_np, cartesian_dis,                                   &
         MPI_REAL,0,cartesian_communicator, mpierr)
      call mpi_gatherv(h2osub1d_n,iv,MPI_REAL,h2o1d_n,                  &
         cartesian_np, cartesian_dis,                                   &
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
      if (iswrit /= 0) then
          call createnc2(hourstep)
      end if
      end if
!-----------------------------------------------------------------------
!     Write output to file
!-----------------------------------------------------------------------
!
      if (myid == 0 ) then 
        port_outs = port
        where(chym_drai < 50) port_outs = 0
        
        port_out=int((port_outs-add_offset)/scale_factor)
        call add_timestep(chymout%ncid,chymout%varid(3),iostep)
        call write_dynvar(chymout%ncid,chymout%varid(4),port_out,      &
             iostep)
        do j=2,nbc-1
          do i=2,nlc-1
            if (port_qmaxs(i,j) < port_outs(i,j)) then
              port_qmaxs(i,j) = port_outs(i,j)
            end if
          end do
        end do
          call createnc3(hourstep)
      end if
      end do
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

      subroutine runoffspeed
      implicit none
      integer i,j,idir,land,intstep
      real mann,wk(nlc,nbc)
      real vmax,alfamin,enne,gamma,delta,tresh,hrad
      alfa=0.0
      gamma=0.33
      delta=4.5                                       !cpar(8) in CHyM
      tresh=100.0                                     !cpar(6) in CHyM
      alfamin=0.1
      do j=2,nbc-1
        do i=2,nlc-1
          if (chym_lat(i,j) > 50.50407 .and. chym_lat(i,j) < 51.79574   &
              .and. chym_lon(i,j) > 5.7957 .and. chym_lon(i,j) < 7.5957 &
              .and. luse(i,j) == mare ) then
              luse(i,j) = lago  
          end if
          idir=fmap(i,j) ; land=luse(i,j) ; mann=manning(luse(i,j))
          if (idir.ge.1.and.idir.le.8.and.land.ne.mare.and.land.gt.0)   &
            then
            if (land.gt.lntypes.or.land.le.0) then
              print*,"Error in line: 845   in file: mod_hd_io"
              call exit(0)
            end if
            chym_dx(i,j)=geodistance(chym_lat(i,j),chym_lon(i,j),       &
                 chym_lat(i+ir(idir),j+jr(idir)),                       &
                 chym_lon(i+ir(idir),j+jr(idir)))
            if (chym_drai(i,j).gt.tresh) then
               enne=mann/delta
            else
               enne=mann/(1+(delta-1)*(1+(chym_drai(i,j)-tresh)/tresh))
            endif
            hrad=0.0015+0.050*((chym_drai(i,j)*1.e00)**gamma)
            !In CHyM 0.0015 = cpar( 2) ---> Alpha coefficients for
            !hydraulic radius (0.0015)
            !In CHyM 0.050 = cpar( 3) ---> Beta coefficients for
            !hydraulic radius (0.050)
            alfa(i,j)=((hrad**0.6666*accl(i,j)**0.5)/(enne))
            if (chym_drai(i,j)>5000 .and. alfa(i,j)>0.5) alfa(i,j) = 0.5
            if (alfa(i,j).lt.alfamin) alfa(i,j)=alfamin
          endif
        enddo
      enddo
      end subroutine runoffspeed
!
      real function geodistance(latt1,lonn1,latt2,lonn2)
      implicit none
      real rad,dpi ; parameter(rad=6371000.0,dpi=6.2831855)
      real latt1,lonn1,latt2,lonn2,lt1,lt2,ln1,ln2,x,y
      lt1=latt1*dpi/360. ; lt2=latt2*dpi/360.
      ln1=lonn1*dpi/360. ; ln2=lonn2*dpi/360.
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

      end module mod_iface
