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

module mod_varandtypes

  type model_area
    logical :: bandflag
    logical :: has_bdy
    logical :: has_bdyleft , has_bdyright , has_bdytop , has_bdybottom
    logical :: has_bdytopleft , has_bdytopright
    logical :: has_bdybottomleft , has_bdybottomright
    integer, dimension(2) :: location
    integer :: left , right , top , bottom
    integer :: topleft , topright , bottomleft , bottomright
    integer :: ibt1 , ibt2 , ibt4 , ibt6 , ibb1 , ibb2 , ibb4 , ibb6
    integer :: jbl1 , jbl2 , jbl4 , jbl6 , jbr1 , jbr2 , jbr4 , jbr6
  end type model_area

  !####################### MPI parameters ################################

  integer, public :: mycomm
  integer :: nproc
  integer :: myid
  integer :: njxcpus , niycpus
  integer :: iyp , jxp
  integer :: iypsg , jxpsg

  integer , allocatable, dimension(:) :: cartesian_dis
  integer , allocatable, dimension(:) :: cartesian_np
  integer , allocatable, dimension(:) :: cartesian_displ
  integer , allocatable, dimension(:) :: cartesian_dispne
  integer , allocatable, dimension(:) :: cartesian_dispse
  integer , allocatable, dimension(:) :: cartesian_npoint
  integer , allocatable, dimension(:) :: cartesian_displl
  integer , allocatable, dimension(:) :: cartesian_displw
  integer , allocatable, dimension(:) :: cartesian_displne
  integer , allocatable, dimension(:) :: cartesian_displse
  integer , allocatable, dimension(:) :: cartesian_npointl
  integer , allocatable, dimension(:) :: cartesian_npointw
  integer , allocatable, dimension(:) :: cartesian_npointne
  integer , allocatable, dimension(:) :: cartesian_npointse
  integer , allocatable, dimension(:) :: ide1p, ide2p, jde1p, jde2p
  integer , allocatable, dimension(:) :: iypp , jxpp
  integer , allocatable, dimension(:) :: ide1gap, ide2gap, jde1gap, jde2gap
  integer , allocatable, dimension(:) :: ide1gbp, ide2gbp, jde1gbp, jde2gbp
  integer , allocatable, dimension(:) :: ide1lgap, ide2lgap, jde1lgap, jde2lgap
  integer , allocatable, dimension(:) :: ide1wgap, ide2wgap, jde1wgap, jde2wgap
  integer , allocatable, dimension(:) :: ide1negap, ide2negap, jde1negap, jde2negap
  integer , allocatable, dimension(:) :: ide1segap, ide2segap, jde1segap, jde2segap
  integer , allocatable, dimension(:) :: iypgap , jxpgap
  integer , allocatable, dimension(:) :: iypgbp , jxpgbp
  integer , allocatable, dimension(:) :: iyplgap , jxplgap
  integer , allocatable, dimension(:) :: iypwgap , jxpwgap
  integer , allocatable, dimension(:) :: iypnegap , jxpnegap
  integer , allocatable, dimension(:) :: iypsegap , jxpsegap
  !################### GRID DIMENSION ####################################
  !

  ! Point in Y (latitude) direction

  integer :: iy

  ! Point in X (longitude) direction

  integer :: jx



  !####################### MPI parameters ################################


  ! D stands for DOT
  integer :: ide1 , ide2 ! External i (included bdy) (latitude)
  integer :: jde1 , jde2 ! External j (included bdy) (longitude)
  integer :: ide1l , ide2l ! External i (included bdy) (latitude) for links
  integer :: jde1l , jde2l ! External j (included bdy) (longitude) for links
  integer :: ide1w , ide2w ! External i (included bdy) (latitude) for links
  integer :: jde1w , jde2w ! External j (included bdy) (longitude) for links
  integer :: jde1ne , jde2ne ! External j (included bdy) (longitude) for links
  integer :: ide1ne , ide2ne ! External i (included bdy) (latitude) for links
  integer :: jde1se , jde2se ! External j (included bdy) (longitude) for links
  integer :: idi1 , idi2 ! Internal (excluded first and last line) i
  integer :: jdi1 , jdi2 ! Internal (excluded first and last column) j
  integer :: idii1 , idii2 ! Internal (excluded 2 lines and cols) i
  integer :: jdii1 , jdii2 ! Internal (excluded 2 lines and cols) j

  ! Ghost points state A

  integer , public :: ide1ga , ide2ga , jde1ga , jde2ga
  integer , public :: ide1lga , ide2lga , jde1lga , jde2lga      ! needed for links
  integer , public :: ide1wga , ide2wga , jde1wga , jde2wga      ! needed for links
  integer , public :: ide1nega , ide2nega , jde1nega , jde2nega      ! needed for links
  integer , public :: ide1sega , ide2sega , jde1sega , jde2sega      ! needed for links
  integer , public :: idi1ga , idi2ga , jdi1ga , jdi2ga

  ! Ghost points state B

  integer , public :: ide1gb , ide2gb , jde1gb , jde2gb
  integer , public :: idi1gb , idi2gb , jdi1gb , jdi2gb

  ! SL stencil points

  integer , public :: ide1sl , ide2sl , jde1sl , jde2sl
  integer , public :: idi1sl , idi2sl , jdi1sl , jdi2sl

  ! Global reference in global grid jx*iy of dot points
  ! The CROSS grid is contained within
  integer :: global_dot_jstart
  integer :: global_dot_jend
  integer :: global_dot_istart
  integer :: global_dot_iend

  ! Number of grid points of each processor considering the ghost points

  integer :: iypga, jxpga, nnga
  integer :: iypgb, jxpgb, nngb
  integer :: iyplga, jxplga, nnlga    ! for links north-south
  integer :: iypwga, jxpwga, nnwga    ! for links west-east
  integer :: iypnega, jxpnega, nnnega    ! for links west-east
  integer :: iypsega, jxpsega, nnsega    ! for links west-east

  ! Total number of grid points considering the ghost points

  integer :: nngap
  integer :: nngbp
  integer :: nnlgap
  integer :: nnwgap
  integer :: nnnegap
  integer :: nnsegap

end module mod_varandtypes
