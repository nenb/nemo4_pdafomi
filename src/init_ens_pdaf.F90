!>##Ensemble Initialisation
!>This routine calls the routines for initialising the ensemble.
!>
!>Separate calls are made for the 2D and 3D state variables to
!>allow for differences in how these variables are initialised.
!>
!>The routine is called when the filter is initialized in
!>`PDAF_filter_init`.
!>
!>The routine is called by all filter processes and
!>initializes the ensemble for the *PE-local domain*.
!>
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
                         ens_p, flag)

   USE mod_kind_pdaf
   USE mod_parallel_pdaf, &
      ONLY: mype_ens
   USE mod_assimilation_pdaf, &
      ONLY: istate_ssh, istate_s, istate_t, istate_u, istate_v, &
            screen
   USE mod_statevector_pdaf, &
      ONLY: fill2d_ensarray, fill3d_ensarray

   IMPLICIT NONE

   !> Type of filter to initialize
   INTEGER, INTENT(in) :: filtertype
   !> PE-local state dimension
   INTEGER, INTENT(in) :: dim_p
   !> Size of ensemble
   INTEGER, INTENT(in) :: dim_ens
   !> PE-local model state
   !> It is not necessary to initialize the array 'state_p' for SEIK.
   !> It is available here only for convenience and can be used freely.
   REAL(pwp), INTENT(inout) :: state_p(dim_p)
   !> Array not referenced for SEIK
   REAL(pwp), INTENT(inout) :: Uinv(dim_ens - 1, dim_ens - 1)
   !> PE-local state ensemble
   REAL(pwp), INTENT(out) :: ens_p(dim_p, dim_ens)
   !> PDAF status flag
   INTEGER, INTENT(inout) :: flag

   IF (screen > 0) THEN
      IF (mype_ens == 0) THEN
         WRITE (*, '(/1x,a)') '------- Reading Initial State --------'
         WRITE (*, '(/1x,a)') 'Calling fill2d_ensarray'
         WRITE (*, '(/9x, a, 3x, a,)') "Initial state file:", TRIM(istate_ssh)
      END IF
   END IF

   CALL fill2d_ensarray(istate_ssh, ens_p)

   IF (screen > 0) THEN
      IF (mype_ens == 0) THEN
         WRITE (*, '(/1x,a)') '------- Reading Initial State --------'
         WRITE (*, '(/1x,a)') 'Calling fill3d_ensarray'
         WRITE (*, '(/9x, a, 3x, a, 3x, a, 3x, a, 3x, a)') &
            "Initial state files:", TRIM(istate_t), TRIM(istate_s), &
            TRIM(istate_u), TRIM(istate_v)
      END IF
   END IF

   CALL fill3d_ensarray(istate_t, 'T', ens_p)
   CALL fill3d_ensarray(istate_s, 'S', ens_p)
   CALL fill3d_ensarray(istate_u, 'U', ens_p)
   CALL fill3d_ensarray(istate_v, 'V', ens_p)

END SUBROUTINE init_ens_pdaf
