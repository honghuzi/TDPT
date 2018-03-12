program main
    use MKL_DFTI
    use Const
    use Laser
    use Functions
    implicit none

    real(DOUBLE), dimension(2) :: sigma0, E0, yw, ya, P, y0, y1
    real(DOUBLE), dimension(:), allocatable :: t, ele
    real(DOUBLE) :: omega, tau, dt
    real(DOUBLE) :: E1, E2, Ea, Eb, Ei, Ef
    complex(DOUBLE) :: Ka, Kb
    integer(SINGLE), parameter :: nE = 350
    real(DOUBLE) :: dE, Eg, sigma1, sigma2, sigma3, sigma4
    real(DOUBLE), dimension(nE) :: energy
    complex(DOUBLE), dimension(nE) :: KK
    real(DOUBLE), dimension(nE, nE) :: spectrum
    real(DOUBLE), dimension(nE*nE) :: Ex, Ey
    type(cos2field) :: f1, f2
    integer(SINGLE) :: iE1, iE2, i, j
    character(10) nowtime

    namelist/He/sigma0, E0, yw, ya, P, y0, y1

    call CheckDelete("data/ele.bin")
    call CheckDelete("data/je.bin")
    sigma0 = [9.492d2, 1.369d4]
    E0 = [1.361d1, 1.720d0]
    yw = [2.039d0, 0d0]
    ya = [1.469d0, 3.288d1]
    P = [3.188d0, 2.963d0]
    y0 = [4.434d-1, 0d0]
    y1 = [2.136d0, 0d0]

    open (100, file='input.nml', delim='apostrophe')
    read (unit=100, nml=He)
    close (100)

    omega = 65/27.2114
    tau = 3*2*PI/omega
    f1 = cos2field(ItoE0(1.0d15), omega, 3*PI/2, tau, tau/2)
    f2 = cos2field(ItoE0(1.0d15), omega, 3*PI/2, tau, tau/2 + 10.0)
    ! f1 = cos2field(ItoE0(1.0d15), 1.86, 0.0, 500.0, 500.0)
    ! f2 = cos2field(ItoE0(1.0d15)/2, 1.92, 0.0, 500.0, 500.0)
    dt = 0.02
    nt = int((max(f1%tau + f1%tm, f2%tau + f2%tm))/dt)

    allocate (t(nt), ele(nt))
    forall (i=1:nt) t(i) = i*dt
    do i = 1, nt
        ele(i) = GetEle(f1, t(i)) + GetEle(f2, t(i))
    end do

    open (100, file='data/ele.bin', access='stream')
    write (100) (t(i), ele(i), i=1, nt)
    close (100)

    dE = 0.2
    energy = [1:nE]*dE
    Ei = -79.0
    Eg = -54.4
    call time(nowtime)
    print'("time is:", T20, A)', nowtime
    spectrum = 0d0
    do iE2 = 1, nE
        do iE1 = 1, iE2
            E1 = energy(iE1)
            E2 = energy(iE2)
            Ea = E1 + Eg
            Eb = E2 + Eg
            Ef = E1 + E2
            Ka = shapefunc(Ea, Ei, Ef, t, ele, dt)
            Kb = shapefunc(Eb, Ei, Ef, t, ele, dt)
            sigma1 = sigma(E1, sigma0(1), E0(1), yw(1), ya(1), P(1), y0(1), y1(1))
            sigma2 = sigma(E2, sigma0(2), E0(2), yw(2), ya(2), P(2), y0(2), y1(2))
            sigma3 = sigma(E2, sigma0(1), E0(1), yw(1), ya(1), P(1), y0(1), y1(1))
            sigma4 = sigma(E1, sigma0(2), E0(2), yw(2), ya(2), P(2), y0(2), y1(2))
            spectrum(iE1, iE2) = abs(sqrt(sigma1*sigma2/(Ea - Ei)/(Ef - Ea))*Ka &
                                     + sqrt(sigma3*sigma4/(Eb - Ei)/(Ef - Eb))*Kb)**2
            spectrum(iE2, iE1) = spectrum(iE1, iE2)
            Ex(nE*(iE2 - 1) + iE1) = iE1*dE
            Ey(nE*(iE1 - 1) + iE2) = iE2*dE
            Ex(nE*(iE1 - 1) + iE2) = iE2*dE
            Ey(nE*(iE2 - 1) + iE1) = iE1*dE
        end do

        if (mod(iE2, 20) == 0) then
            print *, "iE2=", iE2
        end if
    end do

    open (100, file='data/je.bin', access='stream')
    write (100) ((Ex(nE*(j - 1) + i), Ey(nE*(i - 1) + j), spectrum(i, j), i=1, nE), j=1, nE)
    close (100)
    call time(nowtime)
    print'("time is:", T20, A)', nowtime
    deallocate (t, ele)
    call system('cd fig & gnuplot ele.gp & gnuplot je.gp')
end program main
