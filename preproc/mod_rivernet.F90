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

module mod_rivernet

  use, intrinsic :: iso_fortran_env
  use mod_common

  public :: carve_rivernet

  private

  contains

  subroutine partition(a, marker)
    implicit none
    real , dimension(:) , intent(inout) :: a
    integer, intent(out) :: marker
    integer :: np , left , right
    real :: temp , pivot

    np = size(a)
    pivot = (a(1) + a(np))/2.0
    left = 0
    right = np + 1

    do while (left < right)
      right = right - 1
      do while (a(right) > pivot)
        right = right-1
      end do
      left = left + 1
      do while (a(left) < pivot)
        left = left + 1
      end do
      if (left < right) then
        temp = a(left)
        a(left) = a(right)
        a(right) = temp
      end if
    end do
    if (left == right) then
      marker = left + 1
    else
      marker = left
    end if
  end subroutine partition

  recursive subroutine qsort(a)
    implicit none
    real, dimension(:) , intent(inout) :: a
    integer :: np , isplit

    np = size(a)
    if (np > 1) then
     call partition(a, isplit)
     call qsort(a(:isplit-1))
     call qsort(a(isplit:))
    end if
  end subroutine qsort

  subroutine carve_rivernet
    implicit none
    integer :: i, j

    write(output_unit, '(7x,a,a)') 'Carving Rivernet'

    ! Water should get here! Carve the dem to help it in.
    do j = 2 , nlat-1
      do i = 2 , nlon-1
        if ( luc(i,j) /= 15 ) then
          if ( network(i,j) >= 0 ) then
            dem(i,j) = dem(i,j) - 20.0
            if ( dem(i,j) < 0.0 ) dem(i,j) = 10.0
          end if
        end if
      end do
    end do
    write(output_unit, '(12x,a)') 'Done.'
  end subroutine carve_rivernet

end module mod_rivernet
