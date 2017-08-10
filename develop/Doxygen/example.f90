!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration & Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE: Module Name
!
!> @author
!> Module Author Name and Affiliation
!
! DESCRIPTION: 
!> Brief description of module.
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
 
module MyModule_mod
   
   use AnotherModule_mod
   
   implicit none
   
   public MyModule_type            ! TODO_brief_description
 
   public someFunction             ! TODO_brief_description
   
   ! NOTE_avoid_public_variables_if_possible
   
contains
 
   !---------------------------------------------------------------------------  
   !> @author 
   !> Routine Author Name and Affiliation.
   !
   ! DESCRIPTION: 
   !> Brief description of routine. 
   !> @brief
   !> Flow method (rate of change of position) used by integrator.
   !> Compute \f$ \frac{d\lambda}{dt} , \frac{d\phi}{dt},  \frac{dz}{dt} \f$
   !
   ! REVISION HISTORY:
   ! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
   !
   !> @param[in] inParam      
   !> @param[out] outParam      
   !> @return returnValue
   !---------------------------------------------------------------------------  
   
   function someFunction
      use AnotherModule  
      
      real, intent(in) :: inParam        
      real, intent(out) :: outParam       
      real, intent(inout) :: inOutParam   !TODO_description
      real :: returnValue                 
      
      !> @var Variable description
      real :: someVariable
 
      ! TODO_insert_code_here
 
   end function someFunction
   
end module MyModule_mod
