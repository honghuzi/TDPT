module Functions
  use Const, only:DOUBLE, SINGLE, IM
  implicit none
  private
  public::sigma, shapefunc, CheckDelete

contains

  function sigma(E, paras)! sigma0, E0, yw, ya, P, y0, y1)
    real(DOUBLE), intent(in) :: E!, sigma0, E0, yw, ya, P, y0, y1
    real(DOUBLE), dimension(7), intent(in) :: paras!, sigma0, E0, yw, ya, P, y0, y1
    real(DOUBLE) :: sigma, E0, sigma0,  ya, P, yw, y0, y1
    real(DOUBLE) :: x, y, F
    E0 = paras(1)
    sigma0 = paras(2)
    ya = paras(3)
    P = paras(4)
    yw = paras(5)
    y0 = paras(6)
    y1 = paras(7)

    x = E/E0 - y0
    y = sqrt(x**2 + y1**2)
    F = ((x - 1d0)**2 + yw**2)*y**(0.5*P - 5.5)*(1 + sqrt(y/ya))**(-P)

    sigma = sigma0*F
  end function sigma

  function shapefunc(E, Ei, Ef, t, ele, dt)
    real(DOUBLE), intent(in) :: E, Ei, Ef, dt
    real(DOUBLE), dimension(:), intent(in) :: t, ele
    complex(DOUBLE) :: shapefunc, mid
    real(DOUBLE) :: w1, w2
    integer(SINGLE) :: nt, i, j

    w1 = E - Ei
    w2 = Ef - E
    nt = size(t)
    shapefunc = 0.0d0

    do j = 1, nt
       mid = 0d0
       do i = 1, j
          mid = mid + ele(i)*exp(IM*w1*t(i))
       end do
       shapefunc = shapefunc + mid*exp(IM*w2*t(j))*ele(j)
    end do

    shapefunc = shapefunc*dt**2
  end function shapefunc

  subroutine CheckDelete(filename)
    character(*), intent(in) :: filename
    integer(SINGLE) stat
    open (unit=1234, iostat=stat, file=filename, status='old')
    if (stat == 0) close (1234, status='delete')
  end subroutine CheckDelete

end module Functions
