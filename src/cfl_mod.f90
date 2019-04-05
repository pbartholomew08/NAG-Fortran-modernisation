!> This module calculates the CFL number

module cfl_mod

  use types_mod

  implicit none

  !! Private by default
  private
  public :: fd1d_heat_explicit_cfl

contains

  !> Calculates the CFL number
  !> \begin{equation}
  !>   \text{CFL} = \kappa \frac{\Delta t}{\Delta x^2}
  !> \end{equation}
  subroutine fd1d_heat_explicit_cfl(k, t_num, t_min, t_max, x_num, x_min, &
    x_max, cfl)

    implicit none

    real (kind=dp), intent(in) :: k     !! The heat constant $$\kappa$$
    integer, intent(in) :: t_num        !! Number of intervals in t-axis
    integer, intent(in) :: x_num        !! Number of intervals in x-axis
    real (kind=dp), intent(in) :: t_max !! Upper bound of t-axis
    real (kind=dp), intent(in) :: t_min !! Lower bound of t-axis
    real (kind=dp), intent(in) :: x_max !! Upper bound of x-axis
    real (kind=dp), intent(in) :: x_min !! Lower bound of x-axis

    real (kind=dp), intent(out) :: cfl  !! Calculated CFL number

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
