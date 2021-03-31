!>##Assimilation Parameters
!>This module provides variables needed for the
!>assimilation.
!>
!>See `mod_init_pdaf` for where many of these
!>variables are initialised.
!>
MODULE mod_assimilation_pdaf

   USE mod_kind_pdaf

   IMPLICIT NONE
   SAVE

   !> Global model state dimension
   INTEGER :: dim_state
   !> Model state dimension for PE-local domain
   INTEGER :: dim_state_p
   !> Process-local number of observations
   INTEGER :: dim_obs_p
   !> Vector holding process-local observations
   REAL(pwp), ALLOCATABLE :: obs_p(:)
   !> Vector holding state-vector indices of observations
   INTEGER, ALLOCATABLE :: obs_index_p(:)
   !> Vector holding full vector of observations
   REAL(pwp), ALLOCATABLE :: obs_f(:)
   !> Array for full observation coordinates
   REAL(pwp), ALLOCATABLE :: coords_obs_f(:, :)
   !> Vector holding local state-vector indices of observations
   INTEGER, ALLOCATABLE :: obs_index_l(:)
   !> Vector holding distances of local observations
   REAL(pwp), ALLOCATABLE :: distance_l(:)

   !> Control application of model error
   LOGICAL :: model_error
   !> Amplitude for model error
   REAL(pwp) :: model_err_amp

   !> time step interval between assimilation steps
   INTEGER   :: delt_obs
   !> Number of observations
   INTEGER   :: dim_obs

   !> Control verbosity of PDAF
   !> (0) no outputs, (1) progess info, (2) add timings
   !> (3) debugging output
   INTEGER :: screen
   !> Size of ensemble for SEIK/LSEIK/EnKF/ETKF
   !> Number of EOFs to be used for SEEK
   INTEGER :: dim_ens
   !> Select filter algorithm:
   !> SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
   !> LETKF (5), ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
   INTEGER :: filtertype
   !> Subtype of filter algorithm
   !>   SEEK:
   !>     (0) evolve normalized modes
   !>     (1) evolve scaled modes with unit U
   !>     (2) fixed basis (V); variable U matrix
   !>     (3) fixed covar matrix (V,U kept static)
   !>   SEIK:
   !>     (0) ensemble forecast; new formulation
   !>     (1) ensemble forecast; old formulation
   !>     (2) fixed error space basis
   !>     (3) fixed state covariance matrix
   !>     (4) SEIK with ensemble transformation
   !>   EnKF:
   !>     (0) analysis for large observation dimension
   !>     (1) analysis for small observation dimension
   !>   LSEIK:
   !>     (0) ensemble forecast;
   !>     (2) fixed error space basis
   !>     (3) fixed state covariance matrix
   !>     (4) LSEIK with ensemble transformation
   !>   ETKF:
   !>     (0) ETKF using T-matrix like SEIK
   !>     (1) ETKF following Hunt et al. (2007)
   !>       There are no fixed basis/covariance cases, as
   !>       these are equivalent to SEIK subtypes 2/3
   !>   LETKF:
   !>     (0) LETKF using T-matrix like SEIK
   !>     (1) LETKF following Hunt et al. (2007)
   !>       There are no fixed basis/covariance cases, as
   !>       these are equivalent to LSEIK subtypes 2/3
   !>   ESTKF:
   !>     (0) standard ESTKF
   !>       There are no fixed basis/covariance cases, as
   !>       these are equivalent to SEIK subtypes 2/3
   !>   LESTKF:
   !>     (0) standard LESTKF
   !>       There are no fixed basis/covariance cases, as
   !>       these are equivalent to LSEIK subtypes 2/3
   !>   NETF:
   !>     (0) standard NETF
   !>   LNETF:
   !>     (0) standard LNETF
   INTEGER :: subtype
   !> Perform incremental updating in LSEIK
   INTEGER :: incremental
   !> Number of time instances for smoother
   INTEGER :: dim_lag

   !> Type of forgetting factor
   INTEGER :: type_forget
   !> Forgetting factor for filter analysis
   REAL(pwp) :: forget
   !> dimension of bias vector
   INTEGER :: dim_bias

   !> Interval to perform re-diagonalization in SEEK
   INTEGER :: int_rediag
   !> Epsilon for gradient approx. in SEEK forecast
   REAL(pwp) :: epsilon

   !> Rank to be considered for inversion of HPH
   INTEGER :: rank_analysis_enkf

   !> Type of ensemble transformation
   !> SEIK/LSEIK:
   !> (0) use deterministic omega
   !> (1) use random orthonormal omega orthogonal to (1,...,1)^T
   !> (2) use product of (0) with random orthonormal matrix with
   !>     eigenvector (1,...,1)^T
   !> ETKF/LETKF with subtype=4:
   !> (0) use deterministic symmetric transformation
   !> (2) use product of (0) with random orthonormal matrix with
   !>     eigenvector (1,...,1)^T
   !> ESTKF/LESTKF:
   !> (0) use deterministic omega
   !> (1) use random orthonormal omega orthogonal to (1,...,1)^T
   !> (2) use product of (0) with random orthonormal matrix with
   !>     eigenvector (1,...,1)^T
   !> NETF/LNETF:
   !> (0) use random orthonormal transformation orthogonal to (1,...,1)^T
   !> (1) use identity transformation
   !>    ! LSEIK/LETKF/LESTKF
   INTEGER :: type_trans
   !> Range for local observation domain - NEMO grid
   REAL(pwp) :: local_range
   !> Type of localizing weighting of observations

   !>   (0) constant weight of 1
   !>   (1) exponentially decreasing with SRANGE
   !>   (2) use 5th-order polynomial
   !>   (3) regulated localization of R with mean error variance
   !>   (4) regulated localization of R with single-point error variance
   INTEGER :: locweight

   !> Support range for 5th order polynomial - NEMO grid
   !>   or radius for 1/e for exponential weighting
   !>    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
   REAL(pwp) :: srange

   !> Type of the transform matrix square-root
   !> (0) symmetric square root, (1) Cholesky decomposition
   INTEGER :: type_sqrt

   !> file for ssh initial state estimate
   CHARACTER(lc) :: istate_ssh
   !> file for t initial state estimate
   CHARACTER(lc) :: istate_t
   !> file for t initial state estimate
   CHARACTER(lc) :: istate_s
   !> file for u initial state estimate
   CHARACTER(lc) :: istate_u
   !> file for v initial state estimate
   CHARACTER(lc) :: istate_v

   !> For SEIK: Definition of ensemble covar matrix
   !> (0): Factor (r+1)^-1 (or N^-1)
   !> (1): Factor r^-1 (or (N-1)^-1) - real ensemble covar.
   !> This setting is only for the model part; The definition
   !> of P has also to be specified in PDAF_filter_init.
   !> Only for upward-compatibility of PDAF
   INTEGER :: covartype
   !> model time
   REAL(pwp)    :: time

   !> Coordinates of local analysis domain
   REAL(pwp) :: coords_l(2)
   !> Indices of local analysis domain
   REAL(pwp), ALLOCATABLE :: indx_dom_l(:, :)

!$OMP THREADPRIVATE(coords_l)

CONTAINS

   !>##Performing the Assimilation Step
   !>This routine is called during the model integrations at each timestep.
   !>It calls PDAF to check whether the forecast phase is completed and if
   !>so, PDAF will perform the analysis step.
   !>
   !>**Calling Sequence**
   !>*Called from:* `step.F90'
   !>*Calls:* `PDAFomi_assimilate_local`
   !>
   SUBROUTINE assimilate_pdaf()

      USE pdaf_interfaces_module, &
         ONLY: PDAFomi_assimilate_local, PDAF_get_localfilter
      USE mod_parallel_pdaf, &
         ONLY: mype_ens, abort_parallel

      !> PDAF status flag
      INTEGER :: status_pdaf
      !> Flag for domain-localized filter (1=true)
      INTEGER :: localfilter

      ! Collect a state vector from model fields
      EXTERNAL :: collect_state_pdaf
      ! Distribute a state vector to model fields
      EXTERNAL :: distribute_state_pdaf
      ! Provide time step of next observation
      EXTERNAL :: next_observation_pdaf
      ! User supplied pre/poststep routine
      EXTERNAL :: prepoststep_ens_pdaf
      ! Provide number of local analysis domains
      EXTERNAL :: init_n_domains_pdaf
      ! Initialize state dimension for local analysis domain
      EXTERNAL :: init_dim_l_pdaf
      ! Get state on local analysis domain from global state
      EXTERNAL :: g2l_state_pdaf
      ! Update global state from state on local analysis domain
      EXTERNAL :: l2g_state_pdaf
      ! Get dimension of full obs. vector for PE-local domain
      EXTERNAL :: init_dim_obs_pdafomi
      ! Obs. operator for full obs. vector for PE-local domain
      EXTERNAL ::  obs_op_pdafomi
      ! Get dimension of obs. vector for local analysis domain
      EXTERNAL ::  init_dim_obs_l_pdafomi

      ! *********************************
      ! *** Call assimilation routine ***
      ! *********************************

      ! Check  whether the filter is domain-localized
      CALL PDAF_get_localfilter(localfilter)

      IF (localfilter == 1) THEN
         CALL PDAFomi_assimilate_local(collect_state_pdaf, &
                                       distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
                                       prepoststep_ens_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
                                       init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, &
                                       next_observation_pdaf, status_pdaf)
      ELSE
         WRITE (*, '(a)') 'ERROR - global filter not implemented, stopping.'
         CALL abort_parallel()
      END IF

      ! Check for errors during execution of PDAF
      IF (status_pdaf /= 0) THEN
         WRITE (*, '(/1x,a6,i3,a43,i4,a1/)') &
            'ERROR ', status_pdaf, &
            ' in PDAF_put_state - stopping! (PE ', mype_ens, ')'
         CALL abort_parallel()
      END IF

   END SUBROUTINE assimilate_pdaf

END MODULE mod_assimilation_pdaf
