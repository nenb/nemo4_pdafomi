!>##Controlling Pre- and Post-Processing of the PDAF output
!>
!> - For global filters (e.g. SEIK), the routine is called
!>before the analysis and after the ensemble transformation.
!> - For local filters (e.g. LSEIK), the routine is called
!>before and after the loop over all local analysis
!>domains.
!>
!>The routine provides full access to the state
!>estimate and the state ensemble to the user.
!>Thus, user-controlled pre- and poststep
!>operations can be performed here. For example
!>the forecast and the analysis states and ensemble
!>covariance matrix can be analyzed, e.g. by
!>computing the estimated variances.
!>For the offline mode, this routine is the place
!>in which the writing of the analysis ensemble
!>can be performed.
!>
!>If a user considers to perform adjustments to the
!>estimates (e.g. for balances), this routine is
!>the right place for it.
!>
!>**Calling Sequence**
!>
!> - Called by: `PDAF_get_state` (as U_prepoststep) `PDAF_X_update` (as U_prepoststep)
!>
SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
                                state_p, Uinv, ens_p, flag)

   USE mod_kind_pdaf

   IMPLICIT NONE

   !> Current time step (negative for call after forecast)
   INTEGER, INTENT(in) :: step
   !> PE-local state dimension
   INTEGER, INTENT(in) :: dim_p
   !> Size of state ensemble
   INTEGER, INTENT(in) :: dim_ens
   !> PE-local size of ensemble
   INTEGER, INTENT(in) :: dim_ens_p
   !> PE-local dimension of observation vector
   INTEGER, INTENT(in) :: dim_obs_p
   !> PE-local forecast/analysis state
   !> The array 'state_p' is already initialised and can be used
   !> freely here (not for SEEK!)
   REAL(pwp), INTENT(inout) :: state_p(dim_p)
   !> Inverse of matrix U
   REAL(pwp), INTENT(inout) :: Uinv(dim_ens - 1, dim_ens - 1)
   !> PE-local state ensemble
   REAL(pwp), INTENT(inout) :: ens_p(dim_p, dim_ens)
   !> PDAF status flag
   INTEGER, INTENT(in) :: flag

   ! Deallocate observation arrays
   EXTERNAL :: deallocate_obs_pdafomi

   WRITE (*, '(/1x,a)') 'Insert prepoststep routines here. '

   ! Deallocate observation arrays - DO NOT REMOVE
   CALL deallocate_obs_pdafomi(step)

END SUBROUTINE prepoststep_ens_pdaf
