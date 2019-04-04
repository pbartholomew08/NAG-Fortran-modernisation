module solver_mod

  use types_mod, only : dp
  use rhs_mod, only : func

  !! Private by default
  private
  public :: fd1d_heat_explicit

contains

  subroutine fd1d_heat_explicit(x, t, dt, cfl, h, h_new)
    implicit none

    !! Inputs
    real (kind=dp), intent(in) :: h(:)
    real (kind=dp), intent(in) :: dt
    real (kind=dp), intent(in) :: cfl
    real (kind=dp), intent(in) :: t

    !! Outputs
    real (kind=dp), intent(out) :: h_new(:)
    real (kind=dp), intent(out) :: x(:)

    !! Locals
    real (kind=dp) :: f(size(x))
    integer :: j

    do j = 1, size(x)
      f(j) = func(j, x)
    end do

    h_new(1) = 0.0e+00_dp

    do j = 2, size(h_new) - 1
      h_new(j) = h(j) + dt*f(j) + cfl*(h(j-1)-2.0e+00_dp*h(j)+h(j+1))
    end do

! set the boundary conditions again
    h_new(1) = 90.0e+00_dp
    h_new(size(h_new)) = 70.0e+00_dp
  end subroutine fd1d_heat_explicit

end module solver_mod
