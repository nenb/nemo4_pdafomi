!>##Initialise PDAF
!>This modules contains the initialisation routine for PDAF
!>`init_pdaf`. Here the ensemble is initialised and distributed
!>and the statevector and state variable information is computed.
!>
MODULE mod_init_pdaf

   USE mod_kind_pdaf

   IMPLICIT NONE
   SAVE

CONTAINS

   !> This routine collects the initialization of variables for PDAF.
   !> In addition, the initialization routine `PDAF_init` is called
   !> such that the internal initialization of PDAF is performed.
   !> This variant is for the online mode of PDAF.
   !>
   !> The ensemble is initialised in `init_ens_pdaf`, and is then
   !> distributed to the model in `distribute_state_pdaf`.
   !>
   !> The statevector dimension, and the offset and dimension of the
   !> statevector variables is calculated in `calc_statevector_dim`.
   !>
   !> Much of the initialisation is read from a PDAF-specific namelist.
   !> This is performed in `read_config_pdaf`.
   !>
   !> **Calling Sequence**
   !> *Called from:* `nemogcm.F90`
   !> *Calls:* `calc_statevector_dim`
   !> *Calls:* `read_config_pdaf`
   !> *Calls:* `init_pdaf_info`
   !> *Calls:* `PDAF_init`
   !> *Calls:* `PDAF_get_state`
   !>
   SUBROUTINE init_pdaf()

      USE mod_parallel_pdaf, &
         ONLY: n_modeltasks, task_id, COMM_model, COMM_filter, &
               COMM_couple, mype_ens, filterpe, abort_parallel
      USE mod_assimilation_pdaf, &
         ONLY: dim_state_p, screen, filtertype, subtype, dim_ens, &
               incremental, covartype, type_forget, forget, rank_analysis_enkf, &
               type_trans, type_sqrt, delt_obs, locweight, local_range, srange
      USE mod_statevector_pdaf, &
         ONLY: calc_statevector_dim
      USE mod_util_pdaf, &
         ONLY: init_info_pdaf, read_config_pdaf

      !> Integer parameter array for filter
      INTEGER :: filter_param_i(7)
      !> Real parameter array for filter
      REAL(pwp) :: filter_param_r(2)
      !> PDAF status flag
      INTEGER :: status_pdaf
      !> Not used in this implementation
      INTEGER :: doexit, steps
      !> Not used in this implementation
      REAL(pwp) :: timenow

      ! Ensemble initialization
      EXTERNAL :: init_ens_pdaf
      ! Determine how long until next observation
      EXTERNAL :: next_observation_pdaf
      ! Routine to distribute a state vector to model fields
      EXTERNAL :: distribute_state_pdaf
      ! User supplied pre/poststep routine
      EXTERNAL :: prepoststep_ens_pdaf

      ! ***************************
      ! ***   Initialize PDAF   ***
      ! ***************************

      IF (mype_ens == 0) THEN
         WRITE (*, '(/1x,a)') 'INITIALIZE PDAF - ONLINE MODE'
      END IF

      ! Compute dimension of local statevector. Also compute offset
      ! and dimension of local state variables
      CALL calc_statevector_dim(dim_state_p)

      ! **********************************************************
      ! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
      ! **********************************************************

      ! *** IO options ***
      screen = 2  ! Write screen output (1) for output, (2) add timings

      ! *** Filter specific variables
      filtertype = 5    ! Type of filter
      !   (1) SEIK
      !   (2) EnKF
      !   (3) LSEIK
      !   (4) ETKF
      !   (5) LETKF
      !   (6) ESTKF
      !   (7) LESTKF
      dim_ens = n_modeltasks  ! Size of ensemble for all ensemble filters
      !   We use n_modeltasks here, initialized in init_parallel_pdaf
      subtype = 0       ! subtype of filter:
      !   ESTKF:
      !     (0) Standard form of ESTKF
      !   LESTKF:
      !     (0) Standard form of LESTKF
      type_trans = 0    ! Type of ensemble transformation
      !   SEIK/LSEIK and ESTKF/LESTKF:
      !     (0) use deterministic omega
      !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
      !     (2) use product of (0) with random orthonormal matrix with
      !         eigenvector (1,...,1)^T
      !   ETKF/LETKF:
      !     (0) use deterministic symmetric transformation
      !     (2) use product of (0) with random orthonormal matrix with
      !         eigenvector (1,...,1)^T
      type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
      !   (0) fixed
      !   (1) global adaptive
      !   (2) local adaptive for LSEIK/LETKF/LESTKF
      forget = 1.0     ! Forgetting factor
      type_sqrt = 0     ! Type of transform matrix square-root
      !   (0) symmetric square root, (1) Cholesky decomposition
      incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
      covartype = 1     ! Definition of factor in covar. matrix used in SEIK
      !   (0) for dim_ens^-1 (old SEIK)
      !   (1) for (dim_ens-1)^-1 (real ensemble covariance matrix)
      !   This parameter has also to be set internally in PDAF_init.
      rank_analysis_enkf = 0   ! rank to be considered for inversion of HPH
      ! in analysis of EnKF; (0) for analysis w/o eigendecomposition

      ! *********************************************************************
      ! ***   Settings for analysis steps  - used in call-back routines   ***
      ! *********************************************************************

      ! *** Forecast length (time interval between analysis steps) ***
      delt_obs = 24      ! Number of time steps between analysis/assimilation steps

      ! *** Localization settings
      locweight = 0     ! Type of localizating weighting
      !   (0) constant weight of 1
      !   (1) exponentially decreasing with SRANGE
      !   (2) use 5th-order polynomial
      !   (3) regulated localization of R with mean error variance
      !   (4) regulated localization of R with single-point error variance
      local_range = 0  ! Range in grid points for observation domain in local filters
      srange = local_range  ! Support range for 5th-order polynomial
      ! or range for 1/e for exponential weighting

      ! **************************
      ! Namelist and screen output
      ! **************************

      ! Read namelist file for PDAF
      CALL read_config_pdaf()

      ! Screen output for PDAF parameters
      IF (mype_ens == 0) CALL init_info_pdaf()

      ! *****************************************************
      ! *** Call PDAF initialization routine on all PEs.  ***
      ! ***                                               ***
      ! *** Here, the full selection of filters is        ***
      ! *** implemented. In a real implementation, one    ***
      ! *** reduce this to selected filters.              ***
      ! ***                                               ***
      ! *** For all filters, first the arrays of integer  ***
      ! *** and real number parameters are initialized.   ***
      ! *** Subsequently, PDAF_init is called.            ***
      ! *****************************************************

      whichinit: IF (filtertype == 2) THEN
         ! *** EnKF with Monte Carlo init ***
         filter_param_i(1) = dim_state_p ! State dimension
         filter_param_i(2) = dim_ens     ! Size of ensemble
         filter_param_i(3) = rank_analysis_enkf ! Rank of speudo-inverse in analysis
         filter_param_i(4) = incremental ! Whether to perform incremental analysis
         filter_param_i(5) = 0           ! Smoother lag (not implemented here)
         filter_param_r(1) = forget      ! Forgetting factor

         CALL PDAF_init(filtertype, subtype, 0, &
                        filter_param_i, 6, &
                        filter_param_r, 2, &
                        COMM_model, COMM_filter, COMM_couple, &
                        task_id, n_modeltasks, filterpe, init_ens_pdaf, &
                        screen, status_pdaf)
      ELSE
         ! *** All other filters                       ***
         ! *** SEIK, LSEIK, ETKF, LETKF, ESTKF, LESTKF ***
         filter_param_i(1) = dim_state_p ! State dimension
         filter_param_i(2) = dim_ens     ! Size of ensemble
         filter_param_i(3) = 0           ! Smoother lag (not implemented here)
         filter_param_i(4) = incremental ! Whether to perform incremental analysis
         filter_param_i(5) = type_forget ! Type of forgetting factor
         filter_param_i(6) = type_trans  ! Type of ensemble transformation
         filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
         filter_param_r(1) = forget      ! Forgetting factor

         CALL PDAF_init(filtertype, subtype, 0, &
                        filter_param_i, 7, &
                        filter_param_r, 2, &
                        COMM_model, COMM_filter, COMM_couple, &
                        task_id, n_modeltasks, filterpe, init_ens_pdaf, &
                        screen, status_pdaf)
      END IF whichinit

      ! *** Check whether initialization of PDAF was successful ***
      IF (status_pdaf /= 0) THEN
         WRITE (*, '(/1x,a6,i3,a43,i4,a1/)') &
            'ERROR ', status_pdaf, &
            ' in initialization of PDAF - stopping! (PE ', mype_ens, ')'
         CALL abort_parallel()
      END IF

      ! **********************************
      ! *** Prepare ensemble forecasts ***
      ! **********************************

      CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
                          distribute_state_pdaf, prepoststep_ens_pdaf, status_pdaf)

   END SUBROUTINE init_pdaf

END MODULE mod_init_pdaf
