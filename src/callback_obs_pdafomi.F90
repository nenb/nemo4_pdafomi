!>##PDAFOMI interface routines
!>
!>This file provides interface routines between the call-back routines
!>of PDAF and the observation-specific routines in PDAFOMI. This structure
!>collects all calls to observation-specifc routines in this single file
!>to make it easier to find the routines that need to be adapted.
!>
!>The routines here are mainly pure pass-through routines. Thus they
!>simply call one of the routines from PDAF-OMI.
!>
SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)

   USE mod_kind_pdaf
   USE mod_obs_ssh_mgrid_pdafomi, &
      ONLY: assim_ssh_mgrid, init_dim_obs_ssh_mgrid

   IMPLICIT NONE

   !> Current time step
   INTEGER, INTENT(in)  :: step
   !> Dimension of full observation vector
   INTEGER, INTENT(out) :: dim_obs

   !> Observation dimensions
   INTEGER :: dim_obs_ssh_mgrid

   dim_obs_ssh_mgrid = 0

   IF (assim_ssh_mgrid) CALL init_dim_obs_ssh_mgrid(step, dim_obs_ssh_mgrid)

   dim_obs = dim_obs_ssh_mgrid

END SUBROUTINE init_dim_obs_pdafomi

SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

   USE mod_kind_pdaf
   USE mod_obs_ssh_mgrid_pdafomi, &
      ONLY: obs_op_ssh_mgrid

   IMPLICIT NONE

   !> Current time step
   INTEGER, INTENT(in) :: step
   !> PE-local state dimension
   INTEGER, INTENT(in) :: dim_p
   !> Dimension of full observed state
   INTEGER, INTENT(in) :: dim_obs
   !> PE-local model state
   REAL(pwp), INTENT(in)    :: state_p(dim_p)
   !> PE-local full observed state
   REAL(pwp), INTENT(inout) :: ostate(dim_obs)

   CALL obs_op_ssh_mgrid(dim_p, dim_obs, state_p, ostate)

END SUBROUTINE obs_op_pdafomi

!>This routine calls the routine `PDAFomi_init_dim_obs_l`
!>for each observation type.
!>
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

   USE mod_kind_pdaf
   USE mod_assimilation_pdaf, &
      ONLY: indx_dom_l
   USE mod_obs_ssh_mgrid_pdafomi, &
      ONLY: init_dim_obs_l_ssh_mgrid
   USE dom_oce, &
      ONLY: nldi, nldj, tmask

   IMPLICIT NONE

   !> Index of current local analysis domain
   INTEGER, INTENT(in) :: domain_p
   !> Current time step
   INTEGER, INTENT(in) :: step
   !> Full dimension of observation vector
   INTEGER, INTENT(in) :: dim_obs
   !> Local dimension of observation vector
   INTEGER, INTENT(out) :: dim_obs_l

   !> Halo offset for local PE
   INTEGER :: i0, j0
   !> Grid coordinates for local analysis domain
   INTEGER :: i, j
   ! Indicate whether surface gridpoint land
   LOGICAL :: land

   ! Hack for dealing with case when no valid local domains on PE.
   ! See init_n_domains for details.
   land = .FALSE.
   IF (domain_p == 1) THEN
      ! Compute halo offset
      i0 = nldi - 1
      j0 = nldj - 1
      ! Compute i,j indices
      i = indx_dom_l(1, domain_p)
      j = indx_dom_l(2, domain_p)
      ! Check whether the local domain is actually a land point
      ! (and hence not a valid local domain).
      IF (tmask(i + i0, j + j0, 1) == 0.0_pwp) land = .TRUE.
   END IF

   IF (land) THEN
      dim_obs_l = 0
   ELSE
      CALL init_dim_obs_l_ssh_mgrid(domain_p, step, dim_obs, dim_obs_l)
   END IF

END SUBROUTINE init_dim_obs_l_pdafomi

!>This routine calls the routine `PDAFomi_deallocate_obs`
!>for each observation type.
!>
SUBROUTINE deallocate_obs_pdafomi(step)

   USE mod_kind_pdaf
   USE PDAFomi, &
      ONLY: PDAFomi_deallocate_obs
   USE mod_obs_ssh_mgrid_pdafomi, &
      ONLY: obs_ssh_mgrid => thisobs
   USE mod_assimilation_pdaf, &
      ONLY: indx_dom_l

   IMPLICIT NONE

   !> Current time step
   INTEGER, INTENT(in) :: step

   CALL PDAFomi_deallocate_obs(obs_ssh_mgrid)

   ! Tidy-up from init_n_domains_pdaf
   IF (ALLOCATED(indx_dom_l)) DEALLOCATE (indx_dom_l)

END SUBROUTINE deallocate_obs_pdafomi
