!>##Using the Incremental Analysis Method
!>The material in this module is **heavily** based on a similar
!>implementation for the NEMOVAR data assimilation system. Please
!>refer to the `ASM` subdirectory in the `OCE` source code for
!>precise details.
!>
MODULE mod_iau_pdaf

   USE mod_kind_pdaf
   USE par_oce, &
      ONLY: jpi, jpj, jpk, jpkm1, jp_tem, jp_sal

   IMPLICIT NONE
   SAVE

   !> Array to store the ssh IAU
   REAL(pwp), DIMENSION(:, :), ALLOCATABLE :: ssh_iau_pdaf
   !> Array to store the T, S IAU
   REAL(pwp), DIMENSION(:, :, :), ALLOCATABLE :: t_iau_pdaf, s_iau_pdaf
   !> Array to store the U, V IAU
   REAL(pwp), DIMENSION(:, :, :), ALLOCATABLE :: u_iau_pdaf, v_iau_pdaf

CONTAINS

   !>This routine initialises the arrays for the IAU. The arrays
   !>**must* be initialised to zero, as PDAF does not compute an
   !>analysis update for all gridpoints. (In this case, no increment
   !>should be added to the model at this gridpoint.)
   !>
   !>The arrays are deallocated in `finalise_pdaf`.
   SUBROUTINE asm_inc_init_pdaf()

      ! 2D variables
      ALLOCATE (ssh_iau_pdaf(jpi, jpj))
      ssh_iau_pdaf = 0._pwp

      ! 3D variables
      ALLOCATE (t_iau_pdaf(jpi, jpj, jpk), s_iau_pdaf(jpi, jpj, jpk))
      t_iau_pdaf = 0._pwp
      s_iau_pdaf = 0._pwp
      ALLOCATE (u_iau_pdaf(jpi, jpj, jpk), v_iau_pdaf(jpi, jpj, jpk))
      u_iau_pdaf = 0._pwp
      v_iau_pdaf = 0._pwp

   END SUBROUTINE asm_inc_init_pdaf

   !>This routine is almost identical to a similar routine
   !>from NEMOVAR. It applies the IAU to the dynamical fields.
   !>
   !>*Called from:* `step.F90`
   SUBROUTINE dyn_asm_inc_pdaf(kt)

      USE mod_assimilation_pdaf, &
         ONLY: delt_obs
      USE in_out_manager, &
         ONLY: nit000
      USE oce, &
         ONLY: ua, va

      !> Current time step
      INTEGER, INTENT(IN) :: kt

      !> Counter
      INTEGER :: jk

      ! Check whether to update the dynamic tendencies
      IF (MOD(kt - nit000, delt_obs) == 0 .AND. kt > nit000) THEN
         DO jk = 1, jpkm1
            ua(:, :, jk) = ua(:, :, jk) + u_iau_pdaf(:, :, jk)
            va(:, :, jk) = va(:, :, jk) + v_iau_pdaf(:, :, jk)
         END DO
      END IF

   END SUBROUTINE dyn_asm_inc_pdaf

   !>This routine is almost identical to a similar routine
   !>from NEMOVAR. It applies the IAU to the tracer fields.
   !>
   !>*Called from:* `step.F90`
   SUBROUTINE tra_asm_inc_pdaf(kt)

      USE mod_assimilation_pdaf, &
         ONLY: delt_obs, salfixmin
      USE eosbn2, &
         ONLY: eos_fzp
      USE dom_oce, &
         ONLY: gdept_n
      USE in_out_manager, &
         ONLY: nit000
      USE oce, &
         ONLY: tsn, tsa

      !> Current time step
      INTEGER, INTENT(IN) :: kt

      !> Counter
      INTEGER  :: jk
      !> ! 3d freezing point values
      !> Nick: Taken from NEMOVAR. Will this lead to stack overflow?
      REAL(pwp), DIMENSION(jpi, jpj, jpk) :: fzptnz

      ! Check whether to update the tracer tendencies
      IF (MOD(kt - nit000, delt_obs) == 0 .AND. kt > nit000) THEN
         ! Freezing point calculation taken from oc_fz_pt (but calculated for
         ! all depths). Used to prevent the applied increments taking the
         ! temperature below the local freezing point.
         DO jk = 1, jpkm1
            CALL eos_fzp(tsn(:, :, jk, jp_sal), fzptnz(:, :, jk), gdept_n(:, :, jk))
         END DO

         ! Do not apply nonnegative increments.
         ! Do not apply increments if the temperature will fall below freezing
         ! or if the salinity will fall below a specified minimum value.
         DO jk = 1, jpkm1
            WHERE (t_iau_pdaf(:, :, jk) > 0.0_pwp .OR. &
                   tsn(:, :, jk, jp_tem) + tsa(:, :, jk, jp_tem) + t_iau_pdaf(:, :, jk) &
                   > fzptnz(:, :, jk))
               tsa(:, :, jk, jp_tem) = tsa(:, :, jk, jp_tem) + t_iau_pdaf(:, :, jk)
            END WHERE
            WHERE (s_iau_pdaf(:, :, jk) > 0.0_pwp .OR. &
                   tsn(:, :, jk, jp_sal) + tsa(:, :, jk, jp_sal) + s_iau_pdaf(:, :, jk) &
                   > salfixmin)
               tsa(:, :, jk, jp_sal) = tsa(:, :, jk, jp_sal) + s_iau_pdaf(:, :, jk)
            END WHERE
         END DO
      END IF

   END SUBROUTINE tra_asm_inc_pdaf

   !>This routine is almost identical to a similar routine
   !>from NEMOVAR. It applies the IAU to the ssh divergence term.
   !>
   !>Currently, the implementation is only valid when using the
   !>linear free surface assumption.
   !>
   !>*Called from:* `divhor.F90`
   SUBROUTINE ssh_asm_div_pdaf(kt, phdivn)

      USE mod_assimilation_pdaf, &
         ONLY: delt_obs
      USE mod_parallel_pdaf, &
         ONLY: abort_parallel
      USE dom_oce, &
         ONLY: e3t_n, ln_linssh, tmask
      USE in_out_manager, &
         ONLY: nit000

      !> Current time step
      INTEGER, INTENT(IN) :: kt
      !> Horizontal divergence
      REAL(pwp), DIMENSION(:, :, :), INTENT(inout) :: phdivn

      ! Check whether to update the tracer tendencies
      IF (MOD(kt - nit000, delt_obs) == 0 .AND. kt > nit000) THEN
         ! NEMO-PDAF currently only implemented for linear free surface.
         IF (ln_linssh) THEN
            phdivn(:, :, 1) = phdivn(:, :, 1) - &
                              ssh_iau_pdaf(:, :)/e3t_n(:, :, 1)*tmask(:, :, 1)
         ELSE
            WRITE (*, '(/9x, a)') &
               'NEMO-PDAF has not yet been implemented for nonlinear free &
               &surface (see ssh_asm_div_pdaf).'
            CALL abort_parallel()
         END IF
      END IF

   END SUBROUTINE ssh_asm_div_pdaf

END MODULE mod_iau_pdaf
