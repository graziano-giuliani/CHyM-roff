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

module mod_vector

  public

  contains

  subroutine vec_ordinac(vector,i1,i2)
    implicit none
    integer, intent(in) :: i1, i2
    real, dimension(:), intent(inout) :: vector
    integer :: i , j
    real :: xxx
    do i = i1 , i2 - 1
      do j = i + 1 , i2
        if ( vector(j)<vector(i) ) then
          xxx = vector(i)
          vector(i) = vector(j)
          vector(j) = xxx
        end if
      end do
    end do
  end subroutine vec_ordinac

  subroutine vec_ordinad(vector,i1,i2)
    implicit none
    integer, intent(in) :: i1, i2
    real, dimension(:), intent(inout) :: vector
    integer :: i , j
    real :: xxx
    do i = i1 , i2 - 1
      do j = i + 1 , i2
        if ( vector(j)>vector(i) ) then
          xxx = vector(i)
          vector(i) = vector(j)
          vector(j) = xxx
        end if
      end do
    end do
  end subroutine vec_ordinad

  integer function ifromk(k,ixm)
    implicit none
    integer, intent(in) :: ixm, k
    if ( mod(k,ixm) /= 0 ) then
      ifromk = k - (k/ixm)*ixm
    else
      ifromk = ixm
    end if
  end function ifromk

  integer function jfromk(k,ixm)
    implicit none
    integer, intent(in) :: ixm, k
    if ( mod(k,ixm) /= 0 ) then
      jfromk = mod(k,ixm) + k/ixm - (ifromk(k,ixm)-1)
    else
      jfromk = mod(k,ixm) + k/ixm
    end if
  end function jfromk

  integer function kfromij(i,j,ixm)
    implicit none
    integer, intent(in) :: i, ixm, j
    kfromij = (j-1)*ixm + i
  end function kfromij

  subroutine bitsfrominteger(i,n,ip)
    implicit none
    integer, intent(in) :: i, n
    integer, dimension(n), intent(out) :: ip
    integer :: j , k
    k = i
    do j = n , 1 , -1
      if ( k>=2**(j-1) ) then
        k = k - 2**(j-1)
        ip(j) = 1
      else
        ip(j) = 0
      end if
    end do
  end subroutine bitsfrominteger

  subroutine matxprod(a,c,n,m,k,d)
    implicit none
    integer, intent(in) :: k , m , n
    real, dimension(n,m), intent(in) :: a
    real, dimension(m,k), intent(in) :: c
    real, dimension(n,k), intent(out) :: d
    integer :: i , j , l
    do i = 1 , n
      do l = 1 , k
        d(i,l) = 0.0
        do j = 1 , m
          d(i,l) = d(i,l) + (a(i,j)*c(j,l))
        end do
      end do
    end do
  end subroutine matxprod

  subroutine iexchange(i,j)
    implicit none
    integer, intent(inout) :: i, j
    integer :: l
    l = i
    i = j
    j = l
  end subroutine iexchange

  subroutine rexchange(i,j)
    implicit none
    real, intent(inout) :: i , j
    real :: l
    l = i
    i = j
    j = l
  end subroutine rexchange

end module mod_vector
