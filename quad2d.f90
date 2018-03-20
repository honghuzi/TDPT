module quad2d_qgaus_mod
  use Const, only:SINGLE, DOUBLE
  use Laser
  private
  public quad2d_qgaus
  real(DOUBLE) :: xsav

  interface
     function shapefun(x, y, w1, w2, f1, f2, f3)
       use Const, only:SINGLE, DOUBLE
       use Laser
       implicit none
       real(DOUBLE), intent(in) :: x, w1, w2
       type(cos2field), intent(in) :: f1, f2, f3
       real(DOUBLE), dimension(:), intent(in) :: y
       complex(DOUBLE), dimension(size(y)) :: shapefun
     end function shapefun
  end interface

contains
  function h(x, w1, w2, f1, f2, f3)
    real(DOUBLE), dimension(:), intent(in) :: x
    real(DOUBLE), intent(in) :: w1, w2
    type(cos2field), intent(in) :: f1, f2, f3
    complex(DOUBLE), dimension(size(x)) :: h
    integer(SINGLE) :: i
    do i = 1, size(x)
       xsav = x(i)
       h(i) = qgaus(g, 0d0, xsav, w1, w2, f1, f2, f3)
    end do
  end function h
  function g(y, w1, w2, f1, f2, f3)
    real(DOUBLE), dimension(:), intent(in) :: y
    real(DOUBLE), intent(in) :: w1, w2
    type(cos2field), intent(in) :: f1, f2, f3
    complex(DOUBLE), dimension(size(y)) :: g
    g = shapefun(xsav, y, w1, w2, f1, f2, f3)
  end function g

  recursive function qgaus(func, a, b, w1, w2, f1, f2, f3)
    real(DOUBLE), intent(in) :: a, b, w1, w2
    type(cos2field), intent(in) :: f1, f2, f3
    complex(DOUBLE) :: qgaus
    interface
       function func(x, w1, w2, f1, f2, f3)
         use Const, only:SINGLE, DOUBLE
         use Laser
         real(DOUBLE), dimension(:), intent(in) :: x
         real(DOUBLE), intent(in) :: w1, w2
         type(cos2field), intent(in) :: f1, f2, f3
         complex(DOUBLE), dimension(size(x)) :: func
       end function func
    end interface
    real(DOUBLE) :: xm, xr
    real(DOUBLE), dimension(5) :: dx, w = (/0.2955242247d0, 0.2692667193d0, &
         0.2190863625d0, 0.1494513491d0, 0.0666713443d0/), &
         x = (/0.1488743389d0, 0.4333953941d0, 0.6794095682d0, &
         0.8650633666d0, 0.9739065285d0/)
    xm = 0.5d0*(b + a)
    xr = 0.5d0*(b - a)
    dx(:) = xr*x(:)
    qgaus = xr*sum(w(:)*(func(xm + dx, w1, w2, f1, f2, f3) + func(xm - dx, w1, w2, f1, f2, f3)))
  end function qgaus

  subroutine quad2d_qgaus(x1, x2, ss, w1, w2, f1, f2, f3)
    real(DOUBLE), intent(in) :: x1, x2, w1, w2
    type(cos2field), intent(in) :: f1, f2, f3
    complex(DOUBLE), intent(out) :: ss
    ss = qgaus(h, x1, x2, w1, w2, f1, f2, f3)
  end subroutine quad2d_qgaus
end module quad2d_qgaus_mod

module quad2d_qromb_mod
  use Const, only:SINGLE, DOUBLE
  use nrtype; use nrutil
  use Laser
  private
  public quad2d_qromb
  real(DOUBLE) :: xsav

  interface
     function shapefun(x, y, w1, w2, f1, f2, f3)
       use Const, only:SINGLE, DOUBLE
       use Laser
       implicit none
       real(DOUBLE), intent(in) :: x, w1, w2
       type(cos2field), intent(in) :: f1, f2, f3
       real(DOUBLE), dimension(:), intent(in) :: y
       complex(DOUBLE), dimension(size(y)) :: shapefun
     end function shapefun
  end interface

contains
  function h(x, w1, w2, f1, f2, f3)
    real(DOUBLE), dimension(:), intent(in) :: x
    real(DOUBLE), intent(in) :: w1, w2
    type(cos2field), intent(in) :: f1, f2, f3
    complex(DOUBLE), dimension(size(x)) :: h
    integer(SINGLE) :: i
    do i = 1, size(x)
       xsav = x(i)
       h(i) = qromb(g, 0d0, xsav, w1, w2, f1, f2, f3)
       ! h(i) = qromb(g, 0d0, maxval(x), w1, w2, f1, f2, f3)
       ! h(i) = qromb(g, -sqrt(100**2 - xsav**2), sqrt(100**2-xsav**2), w1, w2, f1, f2, f3)
    end do
  end function h
  function g(y, w1, w2, f1, f2, f3)
    real(DOUBLE), dimension(:), intent(in) :: y
    real(DOUBLE), intent(in) :: w1, w2
    type(cos2field), intent(in) :: f1, f2, f3
    complex(DOUBLE), dimension(size(y)) :: g
    g = shapefun(xsav, y, w1, w2, f1, f2, f3)
  end function g

  recursive function qromb(func, a, b, w1, w2, f1, f2, f3)
    use nrtype; use nrutil, only:nrerror
    implicit none
    real(DOUBLE), intent(in) :: a, b, w1, w2
    type(cos2field), intent(in) :: f1, f2, f3
    complex(DOUBLE) :: qromb

    interface
       function func(x, w1, w2, f1, f2, f3)
         use Const, only:SINGLE, DOUBLE
         use Laser
         real(DOUBLE), dimension(:), intent(in) :: x
         real(DOUBLE), intent(in) :: w1, w2
         type(cos2field), intent(in) :: f1, f2, f3
         complex(DOUBLE), dimension(size(x)) :: func
       end function func
    end interface

    integer(SINGLE), parameter :: jmax = 40, jmaxp = jmax + 1, k = 8, km = k - 1
    real(DOUBLE), parameter :: eps = 1d-4
    real(DOUBLE), dimension(jmaxp) :: h
    complex(DOUBLE), dimension(jmaxp) :: s
    complex(DOUBLE) :: dqromb
    integer(SINGLE) :: j
    h(1) = 1d0
    s = 0d0
    do j = 1, jmax
       call trapzd(func, a, b, s(j), j, w1, w2, f1, f2, f3)
       if (j >= k) then
          call polint(h(j - km:j), s(j - km:j), 0d0, qromb, dqromb)
          if (abs(dqromb) <= eps*abs(qromb)) return
       end if
       s(j + 1) = s(j)
       h(j + 1) = 0.25d0*h(j)
    end do
    call nrerror('qromb: too many steps')
  end function qromb

  subroutine quad2d_qromb(x1, x2, ss, w1, w2, f1, f2, f3)
    real(DOUBLE), intent(in) :: x1, x2, w1, w2
    type(cos2field), intent(in) :: f1, f2, f3
    complex(DOUBLE), intent(out) :: ss
    ss = qromb(h, x1, x2, w1, w2, f1, f2, f3)
  end subroutine quad2d_qromb

  recursive subroutine trapzd(func, a, b, s, n, w1, w2, f1, f2, f3)
    use nrtype; use nrutil, only:arth
    implicit none
    real(DOUBLE), intent(in) :: a, b
    real(DOUBLE), intent(in) :: w1, w2
    type(cos2field), intent(in) :: f1, f2, f3
    complex(DOUBLE), intent(inout) :: s
    integer(SINGLE), intent(in) :: n
    interface
       function func(x, w1, w2, f1, f2, f3)
         use Const, only:SINGLE, DOUBLE
         use Laser
         real(DOUBLE), dimension(:), intent(in) :: x
         real(DOUBLE), intent(in) :: w1, w2
         type(cos2field), intent(in) :: f1, f2, f3
         complex(DOUBLE), dimension(size(x)) :: func
       end function func
    end interface

    real(DOUBLE) :: del
    complex(DOUBLE) :: fsum
    integer(SINGLE) :: it
    if (n == 1) then
       s = 0.5d0*(b - a)*sum(func((/a, b/), w1, w2, f1, f2, f3))
    else
       it = 2**(n - 2)
       del = (b - a)/it
       fsum = sum(func(arth(a + 0.5d0*del, del, it), w1, w2, f1, f2, f3))
       s = 0.5d0*(s + del*fsum)
    end if
  end subroutine trapzd

  subroutine polint(xa, ya, x, y, dy)
    use nrutil, only:assert_eq, iminloc, nrerror
    implicit none
    real(DOUBLE), dimension(:), intent(in) :: xa
    complex(DOUBLE), dimension(:), intent(in) :: ya
    real(DOUBLE), intent(in) :: x
    complex(DOUBLE), intent(out) :: y, dy
    integer(SINGLE) :: m, n, ns
    real(DOUBLE), dimension(size(xa)) :: den, ho
    complex(DOUBLE), dimension(size(xa)) :: c, d
    n = assert_eq(size(xa), size(ya), 'polint')
    c = ya
    d = ya
    ho = xa - x
    ns = iminloc(abs(x - xa))
    y = ya(ns)
    ns = ns - 1
    do m = 1, n - 1
       den(1:n - m) = ho(1:n - m) - ho(1 + m:n)
       if (any(den(1:n - m) == 0.0)) &
            call nrerror('polint: calculation failure')
       den(1:n - m) = REAL(c(2:n - m + 1) - d(1:n - m))/den(1:n - m)
       d(1:n - m) = ho(1 + m:n)*den(1:n - m)
       c(1:n - m) = ho(1:n - m)*den(1:n - m)
       if (2*ns < n - m) then
          dy = c(ns + 1)
       else
          dy = d(ns)
          ns = ns - 1
       end if
       y = y + dy
    end do
  end subroutine polint

end module quad2d_qromb_mod
