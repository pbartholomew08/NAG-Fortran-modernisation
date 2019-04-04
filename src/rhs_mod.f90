module rhs_mod

  use types_mod

  implicit none

  !! Make everything private
  private
  public :: func

contains

  function func(j, x) result (d)
    implicit none

    !! Inputs
    integer, intent(in) :: j
    real (kind=dp), intent(in) :: x(:)

    !! Outputs
    real (kind=dp) :: d

    d = 0.0e+00_dp
  end function

end module rhs_mod
