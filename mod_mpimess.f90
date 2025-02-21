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

module mod_mpimess

  use mpi
  use mod_varandtypes
  implicit none

  interface exchange_bdy
    module procedure exchange_bdy_real8 , &
                     exchange_bdy_integer4
  end interface exchange_bdy

  logical , parameter :: lreorder = .false.
  integer :: cartesian_communicator
  integer :: mpierr
  type(model_area) , public :: ma
  real, allocatable, dimension(:) :: r8vector1,r8vector2
  real :: r81,r82
  integer, allocatable, dimension(:) :: i4vector1,i4vector2
  integer :: i41,i42

  integer , parameter :: tag_bt = 1     ! FROM bottom TO top
  integer , parameter :: tag_tb = 2     ! FROM top TO bottom
  integer , parameter :: tag_lr = 3     ! FROM left TO right
  integer , parameter :: tag_rl = 4     ! FROM right TO left
  integer , parameter :: tag_brtl = 5   ! FROM bottomrigth TO topleft
  integer , parameter :: tag_tlbr = 6   ! FROM topleft TO bottomright
  integer , parameter :: tag_bltr = 7   ! FROM bottomleft TO topright
  integer , parameter :: tag_trbl = 8   ! FROM topright TO bottomleft

  contains

    subroutine setup_model_indexes
      implicit none
      ma%jbl1 = 1
      ma%jbl2 = 2
      ma%jbl4 = 4
      ma%jbl6 = 6
      ma%jbr1 = 1
      ma%jbr2 = 2
      ma%jbr4 = 4
      ma%jbr6 = 6
      ma%ibt1 = 1
      ma%ibt2 = 2
      ma%ibt4 = 4
      ma%ibt6 = 6
      ma%ibb1 = 1
      ma%ibb2 = 2
      ma%ibb4 = 4
      ma%ibb6 = 6
      if ( ma%has_bdyleft ) then
        ma%jbl1 = 0
        ma%jbl2 = 0
        ma%jbl4 = 0
        ma%jbl6 = 0
      end if
      if ( ma%has_bdyright ) then
        ma%jbr1 = 0
        ma%jbr2 = 0
        ma%jbr4 = 0
        ma%jbr6 = 0
      end if
      if ( ma%has_bdytop ) then
        ma%ibt1 = 0
        ma%ibt2 = 0
        ma%ibt4 = 0
        ma%ibt6 = 0
      end if
      if ( ma%has_bdybottom ) then
        ma%ibb1 = 0
        ma%ibb2 = 0
        ma%ibb4 = 0
        ma%ibb6 = 0
      end if
      jde1  = global_dot_jstart
      jdi1  = global_dot_jstart
      jdii1 = global_dot_jstart
      jde2  = global_dot_jend
      jdi2  = global_dot_jend
      jdii2 = global_dot_jend
      ide1  = global_dot_istart
      idi1  = global_dot_istart
      idii1 = global_dot_istart
      ide2  = global_dot_iend
      idi2  = global_dot_iend
      idii2 = global_dot_iend
      if ( ma%has_bdyleft ) then
        jdi1 = jde1 + 1
        jdii1 = jde1 + 2
      end if
      if ( ma%has_bdyright ) then
        jdi2 = jde2 - 1
        jdii2 = jde2 - 2
      end if
      if ( ma%has_bdybottom ) then
        idi1 = ide1 + 1
        idii1 = ide1 + 2
      end if
      if ( ma%has_bdytop ) then
        idi2 = ide2 - 1
        idii2 = ide2 - 2
      end if
      idi1ga = idi1 - ma%ibb1
      idi2ga = idi2 + ma%ibt1
      jdi1ga = jdi1 - ma%jbl1
      jdi2ga = jdi2 + ma%jbr1
      ide1ga = ide1 - ma%ibb1
      ide2ga = ide2 + ma%ibt1
      jde1ga = jde1 - ma%jbl1
      jde2ga = jde2 + ma%jbr1
      idi1gb = idi1 - ma%ibb2
      idi2gb = idi2 + ma%ibt2
      jdi1gb = jdi1 - ma%jbl2
      jdi2gb = jdi2 + ma%jbr2
      ide1gb = ide1 - ma%ibb2
      ide2gb = ide2 + ma%ibt2
      jde1gb = jde1 - ma%jbl2
      jde2gb = jde2 + ma%jbr2

      idi1sl = idi1gb - ma%ibb2
      idi2sl = idi2gb + ma%ibt2
      jdi1sl = jdi1gb - ma%jbl2
      jdi2sl = jdi2gb + ma%jbr2
      ide1sl = ide1gb - ma%ibb2
      ide2sl = ide2gb + ma%ibt2
      jde1sl = jde1gb - ma%jbl2
      jde2sl = jde2gb + ma%jbr2

      iypga = ((ide2ga-ide1ga)+1)
      jxpga = ((jde2ga-jde1ga)+1)
      nnga  = iypga * jxpga
      iypgb = ((ide2gb-ide1gb)+2)
      jxpgb = ((jde2gb-jde1gb)+2)
      nngb  = iypgb * jxpgb
!  For NorthEast and SouthEast directions
      jde1ne  = jde1
      jde2ne  = jde2
      ide1ne  = ide1
      ide2ne  = ide2
      jde2nega = jde2ga
      jde1nega = jde1ga
      ide1nega = ide1ga
      ide2nega = ide2ga
      if (.not.ma%has_bdytop) then
        ide2nega = ide2nega + 1
        ide2ne = ide2ne + 1
      end if
      if (.not.ma%has_bdyright) then
        jde2nega = jde2nega + 1
        jde2ne = jde2ne + 1
      end if

      iypnega = ((ide2nega-ide1nega)+1)
      jxpnega = ((jde2nega-jde1nega)+1)
      nnnega  = iypnega * jxpnega

!  For West-East directions
      jde1w  = jde1
      jde2w  = jde2
      ide1w  = ide1
      ide2w  = ide2
      jde2wga = jde2ga
      jde1wga = jde1ga
      ide1wga = ide1ga
      ide2wga = ide2ga
      ide2wga = ide2wga + 1
      ide2w = ide2w + 1

      iypwga = ((ide2wga-ide1wga)+1)
      jxpwga = ((jde2wga-jde1wga)+1)
      nnwga  = iypwga * jxpwga

!  For North-South directions
      jde1l  = jde1
      jde2l  = jde2
      ide1l  = ide1
      ide2l  = ide2
      jde2lga = jde2ga
      jde1lga = jde1ga
      ide1lga = ide1ga
      ide2lga = ide2ga
      jde2lga = jde2lga + 1
      jde2l = jde2l + 1
      iyplga = ((ide2lga-ide1lga)+1)
      jxplga = ((jde2lga-jde1lga)+1)
      nnlga  = iyplga * jxplga
    end subroutine setup_model_indexes

    subroutine set_nproc
      implicit none
      integer , dimension(2) :: cpus_per_dim
      logical , dimension(2) :: dim_period
        integer , dimension(2) :: isearch
      integer :: imaxcpus , imax1 , imax2 , imiss
      data dim_period /.false.,.false./
    ma%bandflag    = .false.
    ma%top         = mpi_proc_null
    ma%bottom      = mpi_proc_null
    ma%left        = mpi_proc_null
    ma%right       = mpi_proc_null
    ma%bottomleft  = mpi_proc_null
    ma%bottomright = mpi_proc_null
    ma%topleft     = mpi_proc_null
    ma%topright    = mpi_proc_null

    ma%has_bdyright       = .false.
    ma%has_bdyleft        = .false.
    ma%has_bdytop         = .false.
    ma%has_bdybottom      = .false.
    ma%has_bdytopleft     = .false.
    ma%has_bdytopright    = .false.
    ma%has_bdybottomleft  = .false.
    ma%has_bdybottomright = .false.
    if ( nproc < 4 ) then
      cpus_per_dim(1) = nproc
      cpus_per_dim(2) = 1
    else if ( nproc >= 4 ) then
      cpus_per_dim(1) = (nint(sqrt(dble(nproc)))/2)*2
      if ( iy > int(1.5*dble(jx)) ) then
        cpus_per_dim(1) = cpus_per_dim(1) - 1
        do while ( mod(nproc,cpus_per_dim(1)) /= 0 )
          cpus_per_dim(1) = cpus_per_dim(1) - 1
        end do
      else if ( jx > int(1.5*dble(iy)) ) then
        cpus_per_dim(1) = cpus_per_dim(1) + 1
        do while ( mod(nproc,cpus_per_dim(1)) /= 0 )
          cpus_per_dim(1) = cpus_per_dim(1) + 1
        end do
      else
        do while ( mod(nproc,cpus_per_dim(1)) /= 0 )
          cpus_per_dim(1) = cpus_per_dim(1) + 1
        end do
      end if
      cpus_per_dim(2) = nproc/cpus_per_dim(1)
      imaxcpus = cpus_per_dim(1)*cpus_per_dim(2)
      if ( mod(nproc,imaxcpus) /= 0 ) then
        write(*,*) 'Work does not evenly divide.'
        write(*,*) 'I have calculated : '
        write(*,*) 'CPUS DIM1 = ', cpus_per_dim(1)
        write(*,*) 'CPUS DIM2 = ', cpus_per_dim(2)
        imax1 = ((jx/3)/2)*2
        imax2 = ((iy/3)/2)*2
        write(*,*) 'Suggested maximum number of CPUS jx: ', imax1
        write(*,*) 'Suggested maximum number of CPUS iy: ', imax2
        write(*,*) 'Closest number : ' , imaxcpus
      end if
    end if
    call MPI_CART_CREATE(mycomm,2,cpus_per_dim,dim_period,lreorder, &
      cartesian_communicator,mpierr)
    call mpi_comm_rank(cartesian_communicator,myid,mpierr)
    call MPI_CART_COORDS(cartesian_communicator,myid,2,ma%location,mpierr)
    ! Set coordinates in the grid for the other processors
    isearch(1) = ma%location(1)
    isearch(2) = ma%location(2)+1
    if ( isearch(2) < cpus_per_dim(2) ) then
      call mpi_cart_rank(cartesian_communicator,isearch,ma%top,mpierr)
    end if
    isearch(1) = ma%location(1)
    isearch(2) = ma%location(2)-1
    if ( isearch(2) >= 0 ) then
      call mpi_cart_rank(cartesian_communicator,isearch,ma%bottom,mpierr)
    end if
    isearch(1) = ma%location(1)-1
    isearch(2) = ma%location(2)
    if ( ma%bandflag .or. ( isearch(1) >= 0 ) ) then
      call mpi_cart_rank(cartesian_communicator,isearch,ma%left,mpierr)
    end if
    isearch(1) = ma%location(1)+1
    isearch(2) = ma%location(2)
    if ( ma%bandflag .or. ( isearch(1) < cpus_per_dim(1) ) ) then
      call mpi_cart_rank(cartesian_communicator,isearch,ma%right,mpierr)
    end if
    isearch(1) = ma%location(1)+1
    isearch(2) = ma%location(2)+1
    if ( ( ma%bandflag .or. ( isearch(1) < cpus_per_dim(1) ) ) .and. &
           isearch(2) < cpus_per_dim(2) ) then
      call mpi_cart_rank(cartesian_communicator,isearch,ma%topright,mpierr)
    end if
    isearch(1) = ma%location(1)-1
    isearch(2) = ma%location(2)+1
    if ( ( ma%bandflag .or. ( isearch(1) >= 0 ) ) .and. &
           isearch(2) < cpus_per_dim(2) ) then
      call mpi_cart_rank(cartesian_communicator,isearch,ma%topleft,mpierr)
    end if
    isearch(1) = ma%location(1)+1
    isearch(2) = ma%location(2)-1
    if ( ( ma%bandflag .or. ( isearch(1) < cpus_per_dim(1) ) ) .and. &
           isearch(2) >= 0 ) then
      call mpi_cart_rank(cartesian_communicator,isearch,ma%bottomright,mpierr)
    end if
    isearch(1) = ma%location(1)-1
    isearch(2) = ma%location(2)-1
    if ( ( ma%bandflag .or. ( isearch(1) >= 0 ) ) .and. isearch(2) >= 0 ) then
      call mpi_cart_rank(cartesian_communicator,isearch,ma%bottomleft,mpierr)
    end if
    ma%has_bdytop    = (ma%top    == mpi_proc_null)
    ma%has_bdybottom = (ma%bottom == mpi_proc_null)
    ma%has_bdyright  = (ma%right  == mpi_proc_null)
    ma%has_bdyleft   = (ma%left   == mpi_proc_null)
    ma%has_bdy       = ( ma%has_bdytop .or. ma%has_bdybottom .or. &
                         ma%has_bdyright .or. ma%has_bdyleft )
    ma%has_bdytopleft     = (ma%has_bdytop .and. ma%has_bdyleft)
    ma%has_bdytopright    = (ma%has_bdytop .and. ma%has_bdyright)
    ma%has_bdybottomleft  = (ma%has_bdybottom .and. ma%has_bdyleft)
    ma%has_bdybottomright = (ma%has_bdybottom .and. ma%has_bdyright)

    jxp =  jx/cpus_per_dim(1)
    iyp =  iy/cpus_per_dim(2)

    global_dot_jstart = ma%location(1)*jxp+1
    global_dot_istart = ma%location(2)*iyp+1
    if ( jxp * cpus_per_dim(1) < jx ) then
      imiss = jx - jxp * cpus_per_dim(1)
      if ( ma%location(1) < imiss ) then
        global_dot_jstart = global_dot_jstart + ma%location(1)
        jxp = jxp + 1
      else
        global_dot_jstart = global_dot_jstart + imiss
      end if
    end if
    if ( iyp * cpus_per_dim(2) < iy ) then
      imiss = iy - iyp * cpus_per_dim(2)
      if ( ma%location(2) < imiss ) then
        global_dot_istart = global_dot_istart + ma%location(2)
        iyp = iyp + 1
      else
        global_dot_istart = global_dot_istart + imiss
      end if
    end if

    global_dot_jend = global_dot_jstart+jxp-1
    global_dot_iend = global_dot_istart+iyp-1
    if ( global_dot_iend > iy .or. global_dot_jend > jx ) then
      write(*,*) 'Cannot evenly divide!!!!'
      write(*,*) 'Processor ',myid,' has I : ', global_dot_istart, &
                                                     global_dot_iend
      write(*,*) 'Processor ',myid,' has J : ', global_dot_jstart, &
                                                     global_dot_jend
    end if

    if ( jxp < 3 .or. iyp < 3 ) then
      write(*,*) 'CPUS DIM1 = ', cpus_per_dim(1)
      write(*,*) 'CPUS DIM2 = ', cpus_per_dim(2)
      write(*,*) &
        'Cannot have one processor with less than 3x3 points.'
      write(*,*) &
        'Processor ',myid,' has ',jxp*iyp,' (',jxp,'x',iyp,')'
    end if
    end subroutine set_nproc

    subroutine mypack_real_grid(matrix,vector,i1,i2,j1,j2,disp)
      implicit none
      integer :: i, j, i1 , i2 , j1, j2, iv, disp
      real , dimension(:) , intent(inout) :: vector
      real , dimension(:,:) , intent(in) :: matrix
      iv = 1+disp
      do i = i1,i2
        do j = j1,j2
          vector(iv) = matrix(j,i)
          iv = iv + 1
        end do
      end do
    end subroutine mypack_real_grid

    subroutine mypack_integer_grid(matrix,vector,i1,i2,j1,j2,disp)
      implicit none
      integer :: i, j, i1 , i2 , j1, j2, iv, disp
      integer , dimension(:) , intent(inout) :: vector
      integer , dimension(:,:) , intent(in) :: matrix
      iv = 1+disp
      do i = i1,i2
        do j = j1,j2
          vector(iv) = matrix(j,i)
          iv = iv + 1
        end do
      end do
    end subroutine mypack_integer_grid

    subroutine myunpack_real_grid(vector,matrix,i1,i2,j1,j2)
      implicit none
      integer :: i , j , iv, i1,i2,j1,j2
      real , dimension(:) , intent(in) :: vector
      real , dimension(j1:j2,i1:i2) , intent(inout) :: matrix
      iv = 1
      do i = i1,i2
        do j = j1,j2
          matrix(j,i) = vector(iv)
          iv = iv + 1
        end do
      end do
    end subroutine myunpack_real_grid

    subroutine myunpack_integer_grid(vector,matrix,i1,i2,j1,j2)
      implicit none
      integer :: i , j , iv, i1,i2,j1,j2
      integer , dimension(:) , intent(in) :: vector
      integer , dimension(j1:j2,i1:i2) , intent(inout) :: matrix
      iv = 1
      do i = i1,i2
        do j = j1,j2
          matrix(j,i) = vector(iv)
          iv = iv + 1
        end do
      end do
    end subroutine myunpack_integer_grid

    subroutine recv_array_real4(rval,isize,icpu,itag)
      implicit none
      real , dimension(:) , intent(out) :: rval
      integer , intent(in) :: isize , icpu , itag
      call mpi_recv(rval,isize,mpi_real4,icpu,itag, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end subroutine recv_array_real4

    subroutine exchange_bdy_real8(ml,nex,j1,j2,i1,i2,jj1,jj2,ii1,ii2)
      integer :: i , j , isize , jsize , ssize, ib
      integer , intent(in) ::  ii1,ii2,jj1,jj2,nex
      integer , intent(in) ::  i1,i2,j1,j2
      real, dimension(jj1:jj2,ii1:ii2) , intent(inout) :: ml
      isize = i2-i1+1
      jsize = j2-j1+1
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        if (allocated(r8vector1)) deallocate(r8vector1)
        if (allocated(r8vector2)) deallocate(r8vector2)
        allocate(r8vector1(ssize))
        allocate(r8vector2(ssize))
        ib = 1
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j2-j+1,i)
            ib = ib + 1
          end do
        end do
        call exchange_array_real8(ssize,ma%right,tag_rl,tag_lr)
        ib = 1
        do i = i1 , i2
          do j = 1 , nex
            ml(j2+j,i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        if (allocated(r8vector1)) deallocate(r8vector1)
        if (allocated(r8vector2)) deallocate(r8vector2)
        allocate(r8vector1(ssize))
        allocate(r8vector2(ssize))
        ib = 1
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j1+j-1,i)
            ib = ib + 1
          end do
        end do
        call exchange_array_real8(ssize,ma%left,tag_lr,tag_rl)
        ib = 1
        do i = i1 , i2
          do j = 1 , nex
            ml(j1-j,i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end if
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize
        if (allocated(r8vector1)) deallocate(r8vector1)
        if (allocated(r8vector2)) deallocate(r8vector2)
        allocate(r8vector1(ssize))
        allocate(r8vector2(ssize))
        ib = 1
        do i = 1 , nex
          do j = j1 , j2
            r8vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do
        call exchange_array_real8(ssize,ma%top,tag_tb,tag_bt)
        ib = 1
        do i = 1 , nex
          do j = j1 , j2
            ml(j,i2+i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize
        if (allocated(r8vector1)) deallocate(r8vector1)
        if (allocated(r8vector2)) deallocate(r8vector2)
        allocate(r8vector1(ssize))
        allocate(r8vector2(ssize))
        ib = 1
        do i = 1 , nex
          do j = j1 , j2
            r8vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do
        call exchange_array_real8(ssize,ma%bottom,tag_bt,tag_tb)
        ib = 1
        do i = 1 , nex
          do j = j1 , j2
            ml(j,i1-i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end if
      if ( ma%topleft /= mpi_proc_null ) then
        ssize = nex*nex
        if (allocated(r8vector1)) deallocate(r8vector1)
        if (allocated(r8vector2)) deallocate(r8vector2)
        allocate(r8vector1(ssize))
        allocate(r8vector2(ssize))
       ib = 1
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(ib) = ml(j1+j-1,i2-i+1)
           ib = ib + 1
          end do
        end do
        call exchange_array_real8(ssize,ma%topleft,tag_tlbr,tag_brtl)
       ib = 1
       do i = 1 , nex
         do j = 1 , nex
           ml(j1-j,i2+i) = r8vector2(ib)
           ib = ib + 1
         end do
       end do
      end if
      if ( ma%bottomright /= mpi_proc_null ) then
        ssize = nex*nex
        if (allocated(r8vector1)) deallocate(r8vector1)
        if (allocated(r8vector2)) deallocate(r8vector2)
        allocate(r8vector1(ssize))
        allocate(r8vector2(ssize))
        ib = 1
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(ib) = ml(j2-j+1,i1+i-1)
            ib = ib + 1
          end do
        end do
        call exchange_array_real8(ssize,ma%bottomright,tag_brtl,tag_tlbr)
        ib = 1
        do i = 1 , nex
          do j = 1 , nex
            ml(j2+j,i1-i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex
        if (allocated(r8vector1)) deallocate(r8vector1)
        if (allocated(r8vector2)) deallocate(r8vector2)
        allocate(r8vector1(ssize))
        allocate(r8vector2(ssize))
        ib = 1
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(ib) = ml(j2-j+1,i2-i+1)
            ib = ib + 1
          end do
        end do
        call exchange_array_real8(ssize,ma%topright,tag_trbl,tag_bltr)
        ib = 1
        do i = 1 , nex
          do j = 1 , nex
            ml(j2+j,i2+i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex
        if (allocated(r8vector1)) deallocate(r8vector1)
        if (allocated(r8vector2)) deallocate(r8vector2)
        allocate(r8vector1(ssize))
        allocate(r8vector2(ssize))
        ib = 1
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(ib) = ml(j1+j-1,i1+i-1)
            ib = ib + 1
          end do
        end do
        call exchange_array_real8(ssize,ma%bottomleft,tag_bltr,tag_trbl)
        ib = 1
        do i = 1 , nex
          do j = 1 , nex
            ml(j1-j,i1-i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end if
    end subroutine exchange_bdy_real8

    subroutine exchange_bdy_integer4(ml,j1,j2,i1,i2,jj1,jj2,ii1,ii2)
      integer :: i , j , isize , jsize , ssize, ib
      integer , intent(in) ::  ii1,ii2,jj1,jj2
      integer , intent(in) ::  i1,i2,j1,j2
      integer, dimension(jj1:jj2,ii1:ii2) , intent(inout) :: ml
      isize = i2-i1+1
      jsize = j2-j1+1
      if ( ma%right /= mpi_proc_null) then
        ssize = isize
        if (allocated(i4vector1)) deallocate(i4vector1)
        if (allocated(i4vector2)) deallocate(i4vector2)
        allocate(i4vector1(ssize))
        allocate(i4vector2(ssize))
        ib = 1
        do i = i1 , i2
          i4vector1(ib) = ml(j2,i)
          ib = ib + 1
        end do
        call exchange_array_integer4(ssize,ma%right,tag_rl,tag_lr)
        ib = 1
        do i = i1 , i2
          ml(j2+1,i) = i4vector2(ib)
          ib = ib + 1
        end do
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = isize
        if (allocated(i4vector1)) deallocate(i4vector1)
        if (allocated(i4vector2)) deallocate(i4vector2)
        allocate(i4vector1(ssize))
        allocate(i4vector2(ssize))
        ib = 1
        do i = i1 , i2
          i4vector1(ib) = ml(j1,i)
          ib = ib + 1
        end do
        call exchange_array_integer4(ssize,ma%left,tag_lr,tag_rl)
        ib = 1
        do i = i1 , i2
          ml(j1-1,i) = i4vector2(ib)
          ib = ib + 1
        end do
      end if
      if ( ma%top /= mpi_proc_null) then
        ssize = jsize
        if (allocated(i4vector1)) deallocate(i4vector1)
        if (allocated(i4vector2)) deallocate(i4vector2)
        allocate(i4vector1(ssize))
        allocate(i4vector2(ssize))
        ib = 1
        do j = j1 , j2
          i4vector1(ib) = ml(j,i2)
          ib = ib + 1
        end do
        call exchange_array_integer4(ssize,ma%top,tag_tb,tag_bt)
        ib = 1
        do j = j1 , j2
          ml(j,i2+1) = i4vector2(ib)
          ib = ib + 1
        end do
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = jsize
        if (allocated(i4vector1)) deallocate(i4vector1)
        if (allocated(i4vector2)) deallocate(i4vector2)
        allocate(i4vector1(ssize))
        allocate(i4vector2(ssize))
        ib = 1
        do j = j1 , j2
          i4vector1(ib) = ml(j,i1)
          ib = ib + 1
        end do
        call exchange_array_integer4(ssize,ma%bottom,tag_bt,tag_tb)
        ib = 1
        do j = j1 , j2
          ml(j,i1-1) = i4vector2(ib)
          ib = ib + 1
        end do
      end if
      if ( ma%topleft /= mpi_proc_null ) then
        ssize = 1
        i41 = ml(j1,i2)
        call exchange_array_integer4(ssize,ma%topleft,tag_tlbr,tag_brtl)
        ml(j1-1,i2+1) = i42
      end if
      if ( ma%bottomright /= mpi_proc_null ) then
        ssize = 1
        i41 = ml(j2,i1)
        call exchange_array_integer4(ssize,ma%bottomright,tag_brtl,tag_tlbr)
        ml(j2+1,i1-1) = i42
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = 1
        i41 = ml(j2,i2)
        call exchange_array_integer4(ssize,ma%topright,tag_trbl,tag_bltr)
        ml(j2+1,i2+1) = i42
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = 1
        i41 = ml(j1,i1)
        call exchange_array_integer4(ssize,ma%bottomleft,tag_bltr,tag_trbl)
        ml(j1-1,i1-1) = i42
      end if

    end subroutine exchange_bdy_integer4

    subroutine exchange_array_real8(isize,icpu,tag1,tag2)
      implicit none
      integer , intent(in) :: isize , icpu , tag1 , tag2
      integer :: ireq
      call mpi_irecv(r8vector2,isize,mpi_real,icpu,tag1, &
                     cartesian_communicator,ireq,mpierr)
      call mpi_send(r8vector1,isize,mpi_real,icpu,tag2, &
                    cartesian_communicator,mpierr)
      call mpi_wait(ireq,mpi_status_ignore,mpierr)
    end subroutine exchange_array_real8

    subroutine exchange_array_integer4(isize,icpu,tag1,tag2)
      implicit none
      integer , intent(in) :: isize , icpu , tag1 , tag2
      integer :: ireq
      call mpi_irecv(i4vector2,isize,mpi_integer,icpu,tag1, &
                     cartesian_communicator,ireq,mpierr)
      call mpi_send(i4vector1,isize,mpi_integer,icpu,tag2, &
                    cartesian_communicator,mpierr)
      call mpi_wait(ireq,mpi_status_ignore,mpierr)
    end subroutine exchange_array_integer4

end module mod_mpimess
