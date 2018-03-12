module Functions
    use Const, only:DOUBLE, SINGLE, IM
    implicit none
    private
    public::sigma, shapefunc, CheckDelete

contains

    function sigma(E, sigma0, E0, yw, ya, P, y0, y1)
        real(DOUBLE), intent(in) :: E, sigma0, E0, yw, ya, P, y0, y1
        real(DOUBLE) :: sigma
        real(DOUBLE) :: x, y, F

        x = E/E0 - y0
        y = sqrt(x**2 + y1**2)
        F = ((x - 1)**2 + yw**2)*y**(0.5*P - 5.5)*(1 + sqrt(y/ya))**(-P)

        sigma = sigma0*F
    end function sigma

    function shapefunc(E, Ei, Ef, t, ele, dt)
        real(DOUBLE), intent(in) :: E, Ei, Ef, dt
        real(DOUBLE), dimension(:), intent(in) :: t, ele
        complex(DOUBLE) :: shapefunc
        real(DOUBLE) :: w1, w2
        integer(SINGLE) :: nt, i, j

        w1 = E - Ei
        w2 = Ef - E
        nt = size(t)
        shapefunc = 0.0d0

        do j = 1, nt
            do i = 1, j
                shapefunc = shapefunc + ele(i)*ele(j)*exp(IM*w1*t(i))*exp(IM*w2*t(j))
            end do
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
