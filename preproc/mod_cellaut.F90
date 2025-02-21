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

module mod_cellaut

  use, intrinsic :: iso_fortran_env
  use mod_vector

  private

  integer , dimension(512) :: irule

  real, dimension(8) :: rflg = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ]

  public :: d2cellcycle

  ! Unused in code ?
  public :: cellularrule , cellularcycle , cellularinit , cavectorfill

  contains

  subroutine cellularrule(dindex,v)
    implicit none
    integer :: dindex , v
    intent (in) dindex , v
    integer :: i , j , k
    integer , dimension(8) :: ip
    if ( dindex<0 ) then
      do i = 0 , 255
        irule(i+1) = 0
      end do
      if ( dindex==-1 ) then
        do i = 0 , 255
          call bitsfrominteger(i,8,ip)
          k = 0
          do j = 1 , 8
            k = k + ip(j)
          end do
          if ( k<=1 ) irule(i+1) = -1
          if ( k>=4 ) irule(i+1) = -1
          if ( k==3 ) irule(i+1) = 1
        end do
      else if ( dindex==-2 ) then
        irule(66+1) = 12
        irule(68+1) = 13
        irule(65+1) = 11
        do i = 0 , 255
          call bitsfrominteger(i,8,ip)
          if ( ip(1)+ip(2)+ip(3)+ip(4)+ip(5)==0 ) irule(i+1) = -1
        end do
        irule(0+1) = -1
        irule(41+1) = -1
        irule(148+1) = -1
      end if
      return
    end if
    if ( dindex>255 .or. v>2 .or. v<-1 ) then
      write(error_unit,'(7x,a,2i4)') &
          'Cellular rule: bad values ' , dindex , v
    else
      irule(dindex+1) = v
    end if
  end subroutine cellularrule

  subroutine cellularcycle(mat,m,n)
    implicit none
    integer , parameter :: cmax = 200
    integer :: m , n
    integer , dimension(m,n) :: mat
    intent (in) m
    intent (inout) mat , n
    integer :: i , j , k
    integer , dimension(cmax,cmax) :: lmat
    real :: x
    integer , dimension(9) :: xc , yc
    data xc/ -1 ,  0 , +1 , -1 , +1 , -1 ,  0 , +1 , 0 /
    data yc/ +1 , +1 , +1 ,  0 ,  0 , -1 , -1 , -1 , 0 /
    do i = 1 , m
      do j = 1 , n
        if ( mat(i,j)>0 ) then
          lmat(i,j) = 1
        else
          lmat(i,j) = 0
        end if
      end do
    end do
    do i = 2 , m - 1
      do j = 2 , n - 1
        k = lmat(i-1,j+1)*1 + lmat(i,j+1)*2 + lmat(i+1,j+1)*4 + &
            lmat(i-1,j)*8 + lmat(i+1,j)*16 + lmat(i-1,j-1)*32 + &
            lmat(i,j-1)*64 + lmat(i+1,j-1)*128
        if ( irule(k+1)==1 ) then
          mat(i,j) = 1
        else if ( irule(k+1)==-1 ) then
          mat(i,j) = 0
        else if ( irule(k+1)==2 ) then
          if ( lmat(i,j)==0 ) then
            mat(i,j) = 1
          else
            mat(i,j) = 0
          end if
        else if ( irule(k+1)>=11 .and. irule(k+1)<=18 ) then
          mat(i,j) = mat(i+xc(irule(k+1)-10),j+yc(irule(k+1)-10))
        else if ( irule(k+1)==20 ) then
          x = 0.0
          n = 0
          do k = 1 , 8
            if ( lmat(i+xc(k),j+yc(k))==1 ) then
              x = x + mat(i+xc(k),j+yc(k))
              n = n + 1
            end if
          end do
          if ( n/=0 ) mat(i,j) = nint(x/n)
        end if
      end do
    end do
  end subroutine cellularcycle

  subroutine cellularinit(mat,m,n,xtype,xc,yc)
    implicit none
    integer :: m , n , xtype , xc , yc
    integer , dimension(m,n) :: mat
    intent (in) m , n , xtype , xc , yc
    intent (out) mat
    integer :: i
    if ( xtype==1 ) then
      do i = 1 , 3
        mat(xc-1+i,yc) = 1
        mat(xc-1+i,yc-4) = 1
      end do
      mat(xc-1,yc-1) = 1
      mat(xc+3,yc-1) = 1
      mat(xc-1,yc-3) = 1
      mat(xc+3,yc-3) = 1
      mat(xc-2,yc-2) = 1
      mat(xc+4,yc-2) = 1
    else if ( xtype==2 ) then
      do i = xc - 3 , xc + 3
        mat(i,yc) = 1
      end do
      do i = yc - 3 , yc + 3
        mat(xc,i) = 1
      end do
    else if ( xtype==3 ) then
      mat(xc,yc-1) = 1
      mat(xc+1,yc) = 1
      mat(xc-1,yc+1) = 1
      mat(xc-0,yc+1) = 1
      mat(xc+1,yc+1) = 1
    else
      write(error_unit,'(7x,a,i10)') &
          'MVLib cellulardef: Uknown type: ' , xtype
      call exit(0)
    end if
  end subroutine cellularinit

  subroutine d2cellcycle(matr,rule,work,ix,jx,alpha)
    implicit none
    real :: alpha
    integer :: ix , jx
    real , dimension(ix,jx) :: matr , work
    integer , dimension(ix,jx) :: rule
    intent (inout) matr
    integer :: i , j
    logical :: ifirst
    real , dimension(8) , save :: w
    data ifirst/.true./
    if ( ifirst ) then
      do i = 1 , 8
        w(i) = rflg(i)
      end do
      ifirst = .false.
    end if
    do j = 1 , jx
      do i = 1 , ix
        call d2cellapplyrule(matr,rule,work,ix,jx,i,j,w,alpha)
      end do
    end do
    do j = 1 , jx
      do i = 1 , ix
        matr(i,j) = matr(i,j) + work(i,j)
      end do
    end do
  end subroutine d2cellcycle

  subroutine d2cellapplyrule(matr,rule,work,ix,jx,i,j,w,alpha)
    implicit none
    real :: alpha
    integer :: i , ix , j , jx
    real , dimension(ix,jx) :: matr , work
    integer , dimension(ix,jx) :: rule
    real , dimension(8) :: w
    intent (in) alpha , i , ix , j , jx , matr , rule , w
    intent (inout) work
    integer , dimension(8) :: di , dj
    integer :: ii , jj , k , ntot
    real :: wtot
    data di /  0 ,  1 ,  0 , -1 ,  1 ,  1 , -1 , -1 /
    data dj /  1 ,  0 , -1 ,  0 ,  1 , -1 , -1 ,  1 /
    ! rule = 4 o 8 ---> intorno di Von Neumann o Moore
    ! rule = 0     ---> l'automa viene usato ma non updatato
    ! rule < 0     ---> l'automa non viene usato ne updatato
    work(i,j) = 0.0
    if ( rule(i,j)==4 .or. rule(i,j)==8 ) then
      if ( i>1 .and. i<ix .and. j>1 .and. j<jx ) then
        ntot = 0
        wtot = 0.0
        do k = 1 , rule(i,j)
          if ( rule(i+di(k),j+dj(k))>=0 ) then
            ntot = ntot + 1
            wtot = wtot + w(k)
            work(i,j) = work(i,j) + matr(i+di(k),j+dj(k))*w(k)
          end if
        end do
        if ( ntot>0 ) then
          work(i,j) = work(i,j)/wtot
          work(i,j) = (work(i,j)-matr(i,j))*alpha
        end if
      else
        ntot = 0
        wtot = 0.0
        do k = 1 , rule(i,j)
          ii = i + di(k)
          jj = j + dj(k)
          if ( ii>=1 .and. ii<=ix .and. jj>=1 .and. jj<=jx ) then
            if ( rule(ii,jj)>=0 ) then
              ntot = ntot + 1
              wtot = wtot + w(k)
              work(i,j) = work(i,j) + matr(ii,jj)*w(k)
            end if
          end if
        end do
        if ( ntot>0 ) then
          work(i,j) = work(i,j)/wtot
          work(i,j) = (work(i,j)-matr(i,j))*alpha
        end if
      end if
    end if
  end subroutine d2cellapplyrule

  subroutine d1cellcycle(matr,rule,work,n)
    implicit none
    integer :: n
    real , dimension(n) :: matr , work
    integer , dimension(n) :: rule
    intent (in) n , rule
    intent (inout) matr , work
    integer :: i
    do i = 2 , n - 1
      if ( rule(i)==2 ) work(i) = 0.5*(matr(i+1)+matr(i-1))
    end do
    do i = 2 , n - 1
      if ( rule(i)==2 ) matr(i) = work(i)
    end do
  end subroutine d1cellcycle

  subroutine cavectorfill(v,w,z,ndat)
    implicit none
    integer :: ndat
    real , dimension(:) :: v , w
    integer , dimension(:) :: z
    intent (inout) v , z
    real :: aver
    integer :: i , idist , maxd , ncyc , nrep
    nrep = 0
    aver = 0.0
    maxd = 0
    idist = 0
    do i = 1 , ndat
      if ( nint(v(i))==-9999 ) then
        nrep = nrep + 1
        z(i) = 2
        idist = idist + 1
      else
        z(i) = 0
        aver = aver + v(i)
        if ( maxd<idist ) maxd = idist
        idist = 0
      end if
    end do
    if ( maxd<idist ) maxd = idist
    if ( nrep==0 ) then
      return
    else if ( nrep==ndat ) then
      write(error_unit,'(10x,a,i5,a,i5,a)') &
          'CAVectorFill routine: no good data '//  &
          'found, cannot repair the vector.'
    else
      aver = aver/(ndat-nrep)
      do i = 1 , ndat
        if ( z(i)==2 ) v(i) = aver
      end do
      ncyc = maxd*5
      if ( ncyc>5000 ) ncyc = 5000
      do i = 1 , ncyc
        call d1cellcycle(v,z,w,ndat)
      end do
    end if
  end subroutine cavectorfill

end module mod_cellaut
