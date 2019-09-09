!!****ih* source/physics/sourceTerms/Cool/CoolMain/cie/Cool_data
!!
!! NAME
!!  Cool_data
!!
!! SYNOPSIS
!!
!!  use Cool_data
!!
!! DESCRIPTION 
!!  Cool_data is a fortran module that holds variables with
!!  the Cool Unit scope.  All variables located in Cool_data are
!!  accessible to subroutines in the Cool unit only.  (Note this 
!!  is a convention, there is nothing in the fortran language that
!!  prevents another routine from using Cool_data.  It is the FLASH3
!!  architecture that tries to enforce these rules.)
!!  
!!  All variables located in the Cool_data fortran module start with 
!! "cl_".  This is to indicate where they come from, and to make it easier
!! on the developer and user to see which variables are local variables and 
!! which belong to the Cool unit.
!!
!!
!!***

Module Cool_data

  implicit none

#include "constants.h"
#include "Flash.h"

! Runtime Parameters

! NKO: Piecwise polynomial fit
real, save, allocatable, dimension(:) :: cl_tempBounds
real, save, allocatable, dimension(:) :: cl_tempNorms
real, save, allocatable, dimension(:) :: cl_tempExps
real, save, allocatable, dimension(:) :: cl_emBounds  ! Emissivity corresponding to each temp. bound
real, save, allocatable, dimension(:) :: cl_tefBounds ! "temporal evolution function"
character(len=MAX_STRING_LENGTH), save :: cl_mode
integer, save :: cl_length

! NKO: minimum / maximum temperature
real, save :: cl_tempMin, cl_tempMax

! NKO: Define hydrogen mass, boltzmann constant, gravitational constant in CGS
real, save :: cl_mH, cl_kB, cl_G

! NKO: For timestep limiter
real, save :: cl_coolFactor

! NKO: delay time or turn off
real, save :: cl_coolDelay
logical, save :: cl_useCool

! NKO: minimum temperature to cool to
real, save :: cl_tfloor

end Module Cool_data
