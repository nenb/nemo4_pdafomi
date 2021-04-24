!>##Utility Routines
!>This module contains several routines useful for common
!>model tasks. The initial routines included output configuration
!>information about the PDAF library, and configuration information
!>about the assimilation parameters.
!>
MODULE mod_util_pdaf

   USE mod_kind_pdaf

   IMPLICIT NONE
   SAVE

CONTAINS

   !> This routine performs a model-sided screen output about
   !> the coniguration of the data assimilation system.
   !>
   !> **Calling Sequence**
   !>
   !> - Called from: `init_pdaf`
   SUBROUTINE init_info_pdaf()

      USE mod_assimilation_pdaf, & ! Variables for assimilation
         ONLY: filtertype, subtype, dim_ens, delt_obs, model_error, &
               model_err_amp, forget, rank_analysis_enkf, int_rediag

      ! *****************************
      ! *** Initial Screen output ***
      ! *****************************

      IF (filtertype == 0) THEN
         WRITE (*, '(/21x, a)') 'Filter: SEEK'
         IF (subtype == 2) THEN
            WRITE (*, '(6x, a)') '-- fixed basis filter with update of matrix U'
            WRITE (*, '(6x, a)') '-- no re-diagonalization of VUV^T'
         ELSE IF (subtype == 3) THEN
            WRITE (*, '(6x, a)') '-- fixed basis filter & no update of matrix U'
            WRITE (*, '(6x, a)') '-- no re-diagonalization of VUV^T'
         ELSE IF (subtype == 5) THEN
            WRITE (*, '(6x, a)') '-- Offline mode'
         END IF
         WRITE (*, '(13x, a, i5)') 'number of EOFs:', dim_ens
         IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
         WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
         IF (subtype /= 5) THEN
            IF ((int_rediag > 0) .AND. ((subtype /= 2) .OR. (subtype /= 3))) &
               WRITE (*, '(10x, a, i4, a)') &
               'Re-diag each ', int_rediag, '-th analysis step'
         ELSE
            IF (int_rediag == 1) THEN
               WRITE (*, '(10x, a)') 'Perform re-diagonalization'
            ELSE
               WRITE (*, '(10x, a)') 'No re-diagonalization'
            END IF
         END IF
      ELSE IF (filtertype == 1) THEN
         WRITE (*, '(21x, a)') 'Filter: SEIK'
         IF (subtype == 2) THEN
            WRITE (*, '(6x, a)') '-- fixed error-space basis'
         ELSE IF (subtype == 3) THEN
            WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
         ELSE IF (subtype == 4) THEN
            WRITE (*, '(6x, a)') '-- use ensemble transformation'
         ELSE IF (subtype == 5) THEN
            WRITE (*, '(6x, a)') '-- Offline mode'
         END IF
         WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
         IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
         WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
         IF (model_error) THEN
            WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
         END IF
      ELSE IF (filtertype == 2) THEN
         WRITE (*, '(21x, a)') 'Filter: EnKF'
         IF (subtype == 5) THEN
            WRITE (*, '(6x, a)') '-- Offline mode'
         END IF
         WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
         IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
         WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
         IF (model_error) THEN
            WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
         END IF
         IF (rank_analysis_enkf > 0) THEN
            WRITE (*, '(6x, a, i5)') &
               'analysis with pseudo-inverse of HPH, rank:', rank_analysis_enkf
         END IF
      ELSE IF (filtertype == 3) THEN
         WRITE (*, '(21x, a)') 'Filter: LSEIK'
         IF (subtype == 2) THEN
            WRITE (*, '(6x, a)') '-- fixed error-space basis'
         ELSE IF (subtype == 3) THEN
            WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
         ELSE IF (subtype == 4) THEN
            WRITE (*, '(6x, a)') '-- use ensemble transformation'
         ELSE IF (subtype == 5) THEN
            WRITE (*, '(6x, a)') '-- Offline mode'
         END IF
         WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
         IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
         WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
         IF (model_error) THEN
            WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
         END IF
      ELSE IF (filtertype == 4) THEN
         WRITE (*, '(21x, a)') 'Filter: ETKF'
         IF (subtype == 0) THEN
            WRITE (*, '(6x, a)') '-- Variant using T-matrix'
         ELSE IF (subtype == 1) THEN
            WRITE (*, '(6x, a)') '-- Variant following Hunt et al. (2007)'
         ELSE IF (subtype == 5) THEN
            WRITE (*, '(6x, a)') '-- Offline mode'
         END IF
         WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
         IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
         WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
         IF (model_error) THEN
            WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
         END IF
      ELSE IF (filtertype == 5) THEN
         WRITE (*, '(21x, a)') 'Filter: LETKF'
         IF (subtype == 0) THEN
            WRITE (*, '(6x, a)') '-- Variant using T-matrix'
         ELSE IF (subtype == 1) THEN
            WRITE (*, '(6x, a)') '-- Variant following Hunt et al. (2007)'
         ELSE IF (subtype == 5) THEN
            WRITE (*, '(6x, a)') '-- Offline mode'
         END IF
         WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
         IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
         WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
         IF (model_error) THEN
            WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
         END IF
      ELSE IF (filtertype == 6) THEN
         WRITE (*, '(21x, a)') 'Filter: ESTKF'
         IF (subtype == 0) THEN
            WRITE (*, '(6x, a)') '-- Standard mode'
         ELSE IF (subtype == 5) THEN
            WRITE (*, '(6x, a)') '-- Offline mode'
         END IF
         WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
         IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
         WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
         IF (model_error) THEN
            WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
         END IF
      ELSE IF (filtertype == 7) THEN
         WRITE (*, '(21x, a)') 'Filter: LESTKF'
         IF (subtype == 0) THEN
            WRITE (*, '(6x, a)') '-- Standard mode'
         ELSE IF (subtype == 5) THEN
            WRITE (*, '(6x, a)') '-- Offline mode'
         END IF
         WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
         IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
         WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
         IF (model_error) THEN
            WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
         END IF
      END IF

   END SUBROUTINE init_info_pdaf

   !> This routine reads the namelist file with parameters
   !> controlling data assimilation with PDAF and outputs to
   !> screen.
   !>
   !> **Calling Sequence**
   !>
   !> - Called from: `init_pdaf`
   SUBROUTINE read_config_pdaf()

      USE mod_parallel_pdaf, &
         ONLY: mype_ens
      USE mod_assimilation_pdaf, &
         ONLY: filtertype, subtype, dim_ens, delt_obs, &
               screen, forget, local_range, locweight, srange, istate_t, &
               istate_s, istate_u, istate_v, istate_ssh

      !> Namelist file
      CHARACTER(lc) :: nmlfile

      NAMELIST /pdaf_nml/ filtertype, subtype, dim_ens, &
         delt_obs, screen, forget, local_range, locweight, &
         srange, istate_s, istate_t, istate_u, istate_v, istate_ssh

      ! ****************************************************
      ! ***   Initialize PDAF parameters from namelist   ***
      ! ****************************************************

      nmlfile = 'namelist.pdaf'

      OPEN (20, file=nmlfile)
      READ (20, NML=pdaf_nml)
      CLOSE (20)

      ! Print PDAF parameters to screen
      showconf: IF (mype_ens == 0) THEN

         WRITE (*, '(/1x,a)') '-- Overview of PDAF configuration --'
         WRITE (*, '(3x,a)') 'PDAF [pdaf_nml]:'
         WRITE (*, '(5x,a,i10)') 'filtertype   ', filtertype
         WRITE (*, '(5x,a,i10)') 'subtype      ', subtype
         WRITE (*, '(5x,a,i10)') 'dim_ens      ', dim_ens
         WRITE (*, '(5x,a,i10)') 'delt_obs     ', delt_obs
         WRITE (*, '(5x,a,i10)') 'screen       ', screen
         WRITE (*, '(5x,a,f10.2)') 'forget       ', forget
         WRITE (*, '(5x,a,es10.2)') 'local_range  ', local_range
         WRITE (*, '(5x,a,i10)') 'locweight    ', locweight
         WRITE (*, '(5x,a,es10.2)') 'srange       ', srange
         WRITE (*, '(5x,a,a)') 'istate_t   ', istate_t
         WRITE (*, '(5x,a,a)') 'istate_s   ', istate_s
         WRITE (*, '(5x,a,a)') 'istate_u   ', istate_u
         WRITE (*, '(5x,a,a)') 'istate_v   ', istate_v
         WRITE (*, '(5x,a,a)') 'istate_ssh ', istate_ssh
         WRITE (*, '(1x,a)') '-- End of PDAF configuration overview --'

      END IF showconf

   END SUBROUTINE read_config_pdaf

   SUBROUTINE finalize_pdaf()

      !>Timing and clean-up of PDAF
      !>
      !> **Calling Sequence**
      !>
      !> - Called from: `nemogcm`
      !>
      !> - Calls: `PDAF_deallocate`
      USE mod_parallel_pdaf, &
         ONLY: mype_ens

      ! Show allocated memory for PDAF
      ! DOES NOT CURRENTLY WORK WITH XIOS CONFIGURATION - TBD WITH LARS
      !IF (mype_ens==0) CALL PDAF_print_info(2)

      ! Print PDAF timings onto screen
      ! DOES NOT CURRENTLY WORK WITH XIOS CONFIGURATION - TBD WITH LARS
      !IF (mype_ens==0) CALL PDAF_print_info(1)

      ! *** Deallocate PDAF arrays
      CALL PDAF_deallocate()

   END SUBROUTINE finalize_pdaf

END MODULE mod_util_pdaf
