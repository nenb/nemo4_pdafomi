!>##Determining the Next Analysis Step
!>
!>The subroutine is called before each forecast phase
!>by `PDAF_get_state`. It has to initialize the number
!>of time steps until the next available observation
!>(`nsteps`). It indicates if the data assimilation process
!>is completed such that the ensemble loop in the model
!>routine can be exited.
!>
!>The routine is called by all processes.
!>
!> **Calling Sequence**
!>
!> - Called from: `init_pdaf/PDAF_get_state` (as U_next_obs)
!>
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

   USE mod_kind_pdaf
   USE mod_assimilation_pdaf, &
      ONLY: delt_obs
   USE mod_parallel_pdaf, &
      ONLY: mype_ens
   USE in_out_manager, &
      ONLY: nitend

   IMPLICIT NONE

   !> Number of the current time step
   INTEGER, INTENT(in)  :: stepnow
   !> Number of time steps until next obs
   INTEGER, INTENT(out) :: nsteps
   !> Whether to exit forecasting (1 for exit)
   INTEGER, INTENT(out) :: doexit
   !> Current model (physical) time
   REAL(pwp), INTENT(out) :: time

   ! *******************************************************
   ! *** Set number of time steps until next observation ***
   ! *******************************************************

   ! Not used in this implementation
   time = 0.0

   IF (stepnow + delt_obs <= nitend) THEN
      ! *** During the assimilation process ***
      nsteps = delt_obs   ! This assumes a constant time step interval
      doexit = 0          ! Not used in this implementation

      IF (mype_ens == 0) WRITE (*, '(i7, 3x, a, i7)') &
         stepnow, 'Next observation at time step', stepnow + nsteps
   ELSE
      ! *** End of assimilation process ***
      nsteps = 0          ! No more steps
      doexit = 1          ! Not used in this implementation

      IF (mype_ens == 0) WRITE (*, '(i7, 3x, a)') &
         stepnow, 'No more observations - end assimilation'
   END IF

END SUBROUTINE next_observation_pdaf
