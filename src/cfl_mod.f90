module cfl_mod

  use types_mod

  implicit none

  !! Private by default
  private
  public :: fd1d_heat_explicit_cfl

contains

  subroutine fd1d_heat_explicit_cfl(k, t_num, t_min, t_max, x_num, x_min, &
    x_max, cfl)

    implicit none

    !! Inputs
    real (kind=dp), intent(in) :: k
    integer, intent(in) :: t_num
    integer, intent(in) :: x_num
    real (kind=dp), intent(in) :: t_max
    real (kind=dp), intent(in) :: t_min
    real (kind=dp), intent(in) :: x_max
    real (kind=dp), intent(in) :: x_min

    !! Outputs
    real (kind=dp), intent(out) :: cfl

    !! Locals
    real (kind=dp) :: dx
    real (kind=dp) :: dt

    dx = (x_max-x_min)/real(x_num-1, kind=dp)
    dt = (t_max-t_min)/real(t_num-1, kind=dp)

    cfl = k*dt/dx/dx

    write (*, '(a)') ' '
    write (*, '(a,g14.6)') '  CFL stability criterion value = ', cfl

  end subroutine

end module cfl_mod
