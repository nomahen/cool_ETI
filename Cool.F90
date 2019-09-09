!!****f* source/physics/sourceTerms/Cool/CoolMain/cie/Cool
!!
!! NAME
!!
!!  Cool
!!
!! SYNOPSIS
!!
!!  Cool(integer(IN) :: blockCount
!!       integer(IN) :: blockList(blockCount),
!!          real(IN) :: dt)
!!
!!
!!
!! DESCRIPTION
!!  Apply a piecewise polynomial cooling operator on the list of blocks provided as input
!!
!! ARGUMENTS
!!
!!  blockCount : The number of blocks in the list
!!  blockList(:) : The list of blocks on which to apply the cooling operator
!!  dt : the current timestep
!!
!!***



subroutine Cool(blockCount,blockList,dt)

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getDeltas
  use Eos_interface, ONLY: Eos_wrapped, Eos
  use Driver_data, ONLY: dr_simTime
  use Cool_data

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer, intent(IN) :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  real,intent(IN) :: dt

  ! Variables for getting and putting block data
  integer :: block_no, xCoordSize, yCoordSize, zCoordSize
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, allocatable, dimension(:) :: xCoord, yCoord, zCoord
  real, pointer, dimension(:,:,:,:)            :: solnData

  ! Integers for do loops
  integer :: i,j,k,t

  ! related to cooling
  integer :: temp_index, tef_index
  real    :: cell_temp, cell_temp_new, cell_dens
  real    :: cell_emissivity, cell_tcool, cell_tef
  real, dimension(MDIM) :: delta
  real :: dx,dy,dz,dvol

  ! related to equation of state
  real, dimension(EOS_NUM) :: eosData

  !! !!
  if (.not. cl_useCool) return
  if (dr_simTime .lt. cl_coolDelay) return

  do block_no = 1, blockCount

     ! get dimensions/limits and coordinates
     call Grid_getBlkIndexLimits(blockList(block_no),blkLimits,blkLimitsGC)
     call Grid_getDeltas(blockList(block_no),delta)

     xCoordSize = blkLimitsGC(HIGH,IAXIS)
     yCoordSize = blkLimitsGC(HIGH,JAXIS)
     zCoordSize = blkLimitsGC(HIGH,KAXIS)
     !! allocate space for dimensions
     allocate(xCoord(xCoordSize))
     allocate(yCoord(yCoordSize))
     allocate(zCoord(zCoordSize))

     !! fill x/y/z coordinates for current block
     call Grid_getCellCoords(IAXIS,blockList(block_no),CENTER,.TRUE.,xCoord,xCoordSize)
     call Grid_getCellCoords(JAXIS,blockList(block_no),CENTER,.TRUE.,yCoord,yCoordSize)
     call Grid_getCellCoords(KAXIS,blockList(block_no),CENTER,.TRUE.,zCoord,zCoordSize)

     ! Get a pointer to solution data 
     call Grid_getBlkPtr(blockList(block_no),solnData)

     ! Begin looping trhough cells in block
     ! loop z
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        if (NDIM>2) then
          dz = delta(3)
        endif
        ! loop y
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           if (NDIM>1) then
             dy = delta(2)
           endif
           ! loop x
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
             dx = delta(1)
             dvol = dx*dy*dz

             cell_temp = solnData(TEMP_VAR,i,j,k)
             cell_dens = solnData(DENS_VAR,i,j,k)

             ! Force temperature to stay within bounds
             if (cell_temp .le. cl_tfloor) then
               solnData(TEMP_VAR,i,j,k) = cl_tfloor
               exit
             elseif (cell_temp .ge. cl_tempBounds(cl_length+1)) then
               solnData(TEMP_VAR,i,j,k) = cl_tempBounds(cl_length+1)
               exit
             endif

             temp_index = -1
             do t = 1, cl_length
               if (cell_temp .le. cl_tempBounds(t+1)) then
                 temp_index = t
                 exit
               endif 
             end do ! t loop

             ! Calculate emissivity of this cell
             cell_emissivity = cl_tempNorms(temp_index)*cell_temp**cl_tempExps(temp_index) 
             ! Calculate cooling time of this cell (assumes g = 5/3 and X = 1)
             cell_tcool = cl_kB*2.0*cl_mH*cell_temp/((2.0/3.0)*cell_dens*cell_emissivity) 
             ! Calculate temporal evolution function of cell
             if (cl_tempExps(temp_index) .eq. 1.0) then
                cell_tef = cl_tefBounds(temp_index) + (cl_emBounds(cl_length+1)/cl_emBounds(temp_index))*&
                           (cl_tempBounds(temp_index)/cl_tempBounds(cl_length+1))*log(cl_tempBounds(temp_index)&
                           /cell_temp)
             else
                cell_tef = cl_tefBounds(temp_index) + (1.0/(1.0-cl_tempExps(temp_index)))*(cl_emBounds(cl_length+1)&
                           /cl_emBounds(temp_index))*(cl_tempBounds(temp_index)/cl_tempBounds(cl_length+1))*&
                           (1.0 - (cl_tempBounds(temp_index)/cell_temp)**(cl_tempExps(temp_index)-1.0))
             endif

             cell_tef = cell_tef + (cell_temp/cl_tempBounds(cl_length+1))*&
                        (cl_emBounds(cl_length+1)/cell_emissivity)*(dt/cell_tcool)             

             ! if TEF is greater than the upper bound, the cell is trying to
             ! cool below tfloor; enforce tfloor here
             if (cell_tef .ge. cl_tefBounds(1)) then    
               solnData(TEMP_VAR,i,j,k) = cl_tfloor
               exit
             endif


             tef_index = -1
             do t = 1, cl_length
               if (cell_tef .ge. cl_tefBounds(t+1)) then
                 tef_index = t
                 exit
               endif 
             end do ! tef loop

             ! now get temp update
             if (cl_tempExps(tef_index) .eq. 1.0) then
               cell_temp_new = cl_tempBounds(tef_index)*exp(-(cl_emBounds(tef_index)/cl_emBounds(cl_length+1))*&
                               (cl_tempBounds(cl_length+1)/cl_tempBounds(tef_index))*(cell_tef - cl_tefBounds(tef_index)))
             else  
               cell_temp_new = cl_tempBounds(tef_index)*(1.0 - (1.0 - cl_tempExps(tef_index))*&
                               (cl_emBounds(tef_index)/cl_emBounds(cl_length+1))*(cl_tempBounds(cl_length+1)&
                               /cl_tempBounds(tef_index))*(cell_tef - cl_tefBounds(tef_index)))**(1.0/(1.0 - &
                               cl_tempExps(tef_index)))
             endif

             ! Update solution data
             solnData(TEMP_VAR,i,j,k) = max(cell_temp_new,cl_tfloor)

           end do ! i loop
        end do ! j loop
     end do ! k loop

     call Eos_wrapped(MODE_DENS_TEMP, blkLimitsGC, blockList(block_no))
     call Grid_releaseBlkPtr(blockList(block_no), solnData)

     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)
  end do ! block_no loop

  return

end subroutine Cool
