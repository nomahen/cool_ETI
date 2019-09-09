!!****f* source/physics/sourceTerms/Cool/Cool_init
!!
!! NAME
!!
!!  Cool_init
!!
!!
!! SYNOPSIS
!!
!!  Cool_init()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***



subroutine Cool_init()

  use Cool_data
  use Grid_data, ONLY: gr_globalMe
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use Driver_interface, ONLY: Driver_abortFlash
implicit none

#include "constants.h"
#include "Flash.h"

! NKO: Variables for parsing user-defined piecewise polynomial fit
real :: tmp_bound, tmp_norm, tmp_exp

! general
integer :: AllocateStat, i

! for normalization
real :: mass_unit, length_unit, time_unit, energy_unit

  call RuntimeParameters_get("coolMode", cl_mode)
  call RuntimeParameters_get("maxCoolFactor", cl_coolFactor)
  call RuntimeParameters_get("coolDelay", cl_coolDelay)
  call RuntimeParameters_get("useCool", cl_useCool)
  call RuntimeParameters_get("tfloor", cl_tfloor)

  call PhysicalConstants_get("proton mass", cl_mH)
  call PhysicalConstants_get("Boltzmann", cl_kB)
  call PhysicalConstants_get("Newton", cl_G)

  if (.not. cl_useCool) then
    print*,"[Cool_init.F90] Not using cooling!"
    return
  endif

  ! NKO: Do initialization based on cl_mode
  select case (cl_mode)
    case ("default")

      cl_length = 10

      allocate ( cl_tempBounds(cl_length+1), stat = AllocateStat)
        if (AllocateStat /= 0) stop '*** not enough memory***'
      allocate ( cl_tempNorms(cl_length), stat = AllocateStat)
        if (AllocateStat /= 0) stop '*** not enough memory***'
      allocate ( cl_tempExps(cl_length), stat = AllocateStat)
        if (AllocateStat /= 0) stop '*** not enough memory***'
      allocate ( cl_emBounds(cl_length+1), stat = AllocateStat)
        if (AllocateStat /= 0) stop '*** not enough memory***'
      allocate ( cl_tefBounds(cl_length+1), stat = AllocateStat)
        if (AllocateStat /= 0) stop '*** not enough memory***'
     
      ! NKO: Edges of temperature bins, in Kelvin
      cl_tempBounds(1) = 10**(3.65) 
      cl_tempBounds(2) = 10**(3.9)
      cl_tempBounds(3) = 10**(4.3)
      cl_tempBounds(4) = 10**(4.6)
      cl_tempBounds(5) = 10**(4.9)
      cl_tempBounds(6) = 10**(5.4)
      cl_tempBounds(7) = 10**(5.75)
      cl_tempBounds(8) = 10**(6.3)
      cl_tempBounds(9) = 10**(7.0)
      cl_tempBounds(10) = 10**(7.6)
      cl_tempBounds(11) = 10**(8.0)

      ! NKO: Normalizations for each temperature range
      ! &  emissivity (ergs cm^3 /s)
      cl_tempNorms(1) = (10**(-5.97))**(11.7)
      cl_tempNorms(2) = (10**(-7.85))**(6.15)
      cl_tempNorms(3) = 10**(-21.85)
      cl_tempNorms(4) = 10**(-31.0)
      cl_tempNorms(5) = 10**(-21.2)
      cl_tempNorms(6) = 10**(-10.4)
      cl_tempNorms(7) = 10**(-21.94)
      cl_tempNorms(8) = 10**(-17.73)
      cl_tempNorms(9) = 10**(-18.21)
      cl_tempNorms(10) = 10**(-26.57)

      ! NKO: temperature exponents for each temperature range
      cl_tempExps(1) = 11.7
      cl_tempExps(2) = 6.15
      cl_tempExps(3) = 0.0
      cl_tempExps(4) = 2.0
      cl_tempExps(5) = 0.0
      cl_tempExps(6) = -2.0
      cl_tempExps(7) = 0.0
      cl_tempExps(8) = -2.0/3.0
      cl_tempExps(9) = -0.6
      cl_tempExps(10) = 0.5

      ! NKO: emissivity bounds (corresponding to each temperature bound)
      cl_emBounds(1) = cl_tempNorms(1)*cl_tempBounds(1)**(cl_tempExps(1))  
      cl_emBounds(2) = cl_tempNorms(2)*cl_tempBounds(2)**(cl_tempExps(2))
      cl_emBounds(3) = cl_tempNorms(3)*cl_tempBounds(3)**(cl_tempExps(3))
      cl_emBounds(4) = cl_tempNorms(4)*cl_tempBounds(4)**(cl_tempExps(4))
      cl_emBounds(5) = cl_tempNorms(5)*cl_tempBounds(5)**(cl_tempExps(5))
      cl_emBounds(6) = cl_tempNorms(6)*cl_tempBounds(6)**(cl_tempExps(6))
      cl_emBounds(7) = cl_tempNorms(7)*cl_tempBounds(7)**(cl_tempExps(7)) 
      cl_emBounds(8) = cl_tempNorms(8)*cl_tempBounds(8)**(cl_tempExps(8))
      cl_emBounds(9) = cl_tempNorms(9)*cl_tempBounds(9)**(cl_tempExps(9))
      cl_emBounds(10) = cl_tempNorms(10)*cl_tempBounds(10)**(cl_tempExps(10))
      cl_emBounds(11) = cl_tempNorms(10)*cl_tempBounds(11)**(cl_tempExps(10))

      ! NKO: temporal evolution function bounds (corresponding to each temperature bound)
      ! & see appendix of Townsend 2009
      cl_tefBounds(cl_length+1) = 0.0
      do i = cl_length,1,-1
        if (cl_tempExps(i) .eq. 1.0) then
          cl_tefBounds(i) = cl_tefbounds(i+1) - (cl_emBounds(cl_length+1)/cl_emBounds(i))*&
                            (cl_tempBounds(i)/cl_tempBounds(cl_length+1))*log(cl_tempBounds(i)/(cl_tempBounds(i+1)))
        else
          cl_tefBounds(i) = cl_tefBounds(i+1) - (1.0/(1.0-cl_tempExps(i)))*(cl_emBounds(cl_length+1)/cl_emBounds(i))*&
                            (cl_tempBounds(i)/cl_tempBounds(cl_length+1))*(1.0-(cl_tempBounds(i)/cl_tempBounds(i+1))**&
                            (cl_tempExps(i)-1.0))
        endif
      enddo
    case default
      call Driver_abortFlash("[Cool] Error! Invalid cooling mode")
  end select 
  
  ! NKO: Prevent temperature from falling out of range of cooling function
  cl_tempMin = cl_tempBounds(1)
  cl_tempMax = cl_tempBounds(cl_length+1)
  cl_tfloor = cl_tfloor 
  cl_tfloor = max(cl_tfloor, cl_tempMin)

  if (gr_globalMe  == MASTER_PE) then
  
  
     print*,"######## Checking ETI parameters ##################"
     do i = 1,cl_length
       print*,"Tef bound ", i, ": ", cl_tefBounds(i)
     enddo
     print*,"cl_mode = ", cl_mode
     print*,"cl_coolDelay = ", cl_coolDelay
     print*,"cl_coolFactor = ", cl_coolFactor
     print*,"min temp bound (code unit) = ", cl_tempBounds(1)
     print*,"max temp bound (code unit) = ", cl_tempBounds(cl_length+1)
     print*,"temp floor (code unit) = ", cl_tfloor
     print*,"cl_mH = ", cl_mH
     print*,"cl_kB = ", cl_kB
     print*,"cl_G = ", cl_G
     print*,"##############################################"

  endif

  return
end subroutine Cool_init
