program main
  use MKL_DFTI
  use quad2d_qgaus_mod
  use quad2d_qromb_mod
  use Const
  use Laser
  use Functions
  implicit none

  real(DOUBLE), dimension(7) :: paras1, paras2, paras3, paras4
  real(DOUBLE), dimension(:), allocatable :: t, ele
  real(DOUBLE) :: omega, tau, dt
  real(DOUBLE) :: E1, E2, Ea, Eb, Ei, Ef
  complex(DOUBLE) :: Ka, Kb
  complex(DOUBLE), dimension(:), allocatable :: Kas, Kbs
  integer(SINGLE), parameter :: nE = 200
  real(DOUBLE) :: dE, Eg, sigmas(4)
  real(DOUBLE), dimension(nE) :: energy, ss(nE)
  complex(DOUBLE), dimension(nE) :: KK
  real(DOUBLE), dimension(nE, nE) :: spectrum
  real(DOUBLE), dimension(nE*nE) :: Ex, Ey
  type(cos2field) :: f1, f2
  integer(SINGLE) :: iE1, iE2, i, j
  character(10) nowtime

  ! namelist/He/sigma0, E0, yw, ya, P, y0, y1

  call CheckDelete("data/ele.bin")
  call CheckDelete("data/je.bin")
  call CheckDelete("data/sigma.bin")

  paras1 = [1.361d1, 9.492d2, 1.469d0, 3.188d0, 2.039d0, 4.434d-1, 2.136d0]
  paras2 = [1.720d0, 1.369d4, 3.288d1, 2.963d0, 0d0, 0d0, 0d0]
  paras3 = [1.709d1, 2.16d1, 2.645d2, 4.796d0, 4.185d-1, 1.688d0, 8.942d-1]
  paras4 = [2.556d0, 4.140d0, 1.337d1, 1.191d1, 1.570d0, 6.634d0, 1.272d-1]
  ! open (100, file='input.nml', delim='apostrophe')
  ! read (unit=100, nml=He)
  ! close (100)

  omega = 42d0/transE
  ! tau = 3*2*PI/omega
  tau = 4d-15/2.4189d-17
  f1 = cos2field(ItoE0(1d15), omega, 3*PI/2, tau, tau/2 + 0d0)
  f2 = cos2field(ItoE0(0d15), omega, 3*PI/2, tau, tau/2 + 10d0)
  ! f1 = cos2field(ItoE0(1.0d15), 1.86, 0.0, 500.0, 500.0)
  ! f2 = cos2field(ItoE0(1.0d15)/2, 1.92, 0.0, 500.0, 500.0)
  dt = 0.05
  nt = int((max(f1%tau + f1%tm, f2%tau + f2%tm)+0d0)/dt)

  allocate (t(nt), ele(nt))
  allocate (Kas(nt-1), Kbs(nt-1))
  forall (i=1:nt) t(i) = i*dt
  do i = 1, nt
     ele(i) = GetEle(f1, t(i)) + GetEle(f2, t(i))
  end do

  open (100, file='data/ele.bin', access='stream')
  write (100) (t(i), ele(i), i=1, nt)
  close (100)

  dE = 27.2d0/nE/transE
  ! dE = 2d0/nE
  energy = [1:nE]*dE
  ! energy = [-nE/2:nE/2 - 1]*dE
  Ei = -79.0 / transE
  Eg = -54.4 / transE
  call time(nowtime)
  print'("time is:", T20, A)', nowtime
  spectrum = 0d0
  do i = 1, nE
     ss(i) = sigma(i*dE, paras2)
  end do

  open (100, file="data/sigma.bin", access="stream")
  write (100) (i*dE, ss(i), i=1, nE)
  close (100)

  do iE2 = 1, nE
     do iE1 = 1, iE2
        ! E1 = energy(iE1)
        ! E2 = energy(iE2)
        E1 = (energy(iE1))**2/2
        E2 = (energy(iE2))**2/2
        Ea = E1 + Eg
        Eb = E2 + Eg
        Ef = E1 + E2
        ! Ka = shapefunc(Ea, Ei, Ef, t, ele, dt)
        ! Kb = shapefunc(Eb, Ei, Ef, t, ele, dt)

        ! call quad2d_qromb(t(1), t(nt), Ka, Ea-Ei, Ef-Ea, f1, f2)  ! energy in a. u.
        ! call quad2d_qromb(t(1), t(nt), Kb, Eb-Ei, Ef-Eb, f1, f2)


        do i = 1, nt-1
           call quad2d_qgaus(t(i), t(i+1), Kas(i), Ea-Ei, Ef-Ea, f1, f2)
           call quad2d_qgaus(t(i), t(i+1), Kbs(i), Eb-Ei, Ef-Eb, f1, f2)
        end do
        Ka = sum(Kas)
        Kb = sum(Kbs)

        ! Ka = Ka*(1+exp(IM*(E1+E2-Ei)*10.0-2*PI))
        ! Kb = Kb*(1+exp(IM*(E1+E2-Ei)*10.0-2*PI))
        sigmas(1) = sigma(Ea - Ei, paras1)
        sigmas(2) = sigma(Ef - Ea, paras2)
        sigmas(3) = sigma(Eb - Ei, paras1)
        sigmas(4) = sigma(Ef - Eb, paras2)

        ! spectrum(iE1, iE2) = abs(energy(iE1)*energy(iE2))*abs(sqrt(sigmas(1)*sigmas(2)/(Ea - Ei)/(Ef - Ea))*Ka &
        !      + sqrt(sigmas(3)*sigmas(4)/(Eb - Ei)/(Ef - Eb))*Kb)**2

        spectrum(iE1, iE2) = abs(sqrt(sigmas(1)*sigmas(2)/(Ea - Ei)/(Ef - Ea))*Ka &
             + sqrt(sigmas(3)*sigmas(4)/(Eb - Ei)/(Ef - Eb))*Kb)**2

        spectrum(iE2, iE1) = spectrum(iE1, iE2)
        ! Ex(nE*(iE2 - 1) + iE1) = (iE1 - nE/2)*dE
        ! Ey(nE*(iE1 - 1) + iE2) = (iE2 - nE/2)*dE
        ! Ex(nE*(iE1 - 1) + iE2) = (iE2 - nE/2)*dE
        ! Ey(nE*(iE2 - 1) + iE1) = (iE1 - nE/2)*dE
        Ex(nE*(iE2 - 1) + iE1) = iE1*dE
        Ey(nE*(iE1 - 1) + iE2) = iE2*dE
        Ex(nE*(iE1 - 1) + iE2) = iE2*dE
        Ey(nE*(iE2 - 1) + iE1) = iE1*dE
     end do

     if (mod(iE2, 20) == 0) then
        print *, "iE2=", iE2
     end if
  end do
  Ex = Ex * transE
  Ey = Ey * transE

  open (100, file='data/je.bin', access='stream')
  write (100) ((Ex(nE*(j - 1) + i), Ey(nE*(i - 1) + j), spectrum(i, j), i=1, nE), j=1, nE)
  close (100)

  call time(nowtime)
  print'("time is:", T20, A)', nowtime
  deallocate (t, ele)
  call system('cd fig & gnuplot ele.gp & gnuplot je.gp & gnuplot sigma.gp')
end program main

function shapefun(x, y, w1, w2, f1, f2)  !!! w1, w2 in atomic unit
  use Const
  use Laser
  implicit none
  real(DOUBLE), intent(in) :: x, w1, w2
  type(cos2field) :: f1, f2
  real(DOUBLE) :: Ex
  real(DOUBLE), dimension(:), intent(in) :: y
  real(DOUBLE), dimension(size(y)) :: Ey
  complex(DOUBLE), dimension(size(y)) :: shapefun
  integer(SINGLE) :: i

  Ex = GetEle(f1, x) + GetEle(f2, x)

  do i = 1, size(y)
     Ey(i) = GetEle(f1, y(i)) + GetEle(f2, y(i))
  end do
  shapefun = Ex*Ey*exp(IM*(w1*x + w2*y))
end function shapefun
