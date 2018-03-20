module Const
   use mkl_dfti
   implicit none
   private
   public :: SINGLE, DOUBLE, PI, C, IM, transE,&
             handle, status, nt

   integer, parameter :: SINGLE = KIND(0.0), DOUBLE = KIND(0.0d0)
   real(DOUBLE), parameter :: PI = 3.14159d0
   real(DOUBLE), parameter :: C = 3.0d8
   real(DOUBLE), parameter :: transE = 27.2114d0
   complex(DOUBLE), parameter :: IM = DCMPLX(0.0d0, 1.0d0)

   integer(SINGLE) :: nt, status
   real(DOUBLE) ::  dt
   type(dfti_descriptor), pointer :: handle

end module Const
