module types_mod
  use, intrinsic :: iso_fortran_env

  implicit none

  !! Everything is private unless otherwise stated
  private
  public :: SP

  integer, parameter :: SP = REAL32
  integer, parameter :: DP = REAL64
  integer, parameter :: SI = INT32
  integer, parameter :: DI = INT64

contains

end module types_mod
