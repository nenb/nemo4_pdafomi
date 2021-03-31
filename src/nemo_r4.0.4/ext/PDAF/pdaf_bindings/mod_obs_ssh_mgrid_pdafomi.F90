!>##PDAF-OMI observation module for ssh observations (on model grid)
!>
!>The subroutines in this module are for the particular handling of
!>ssh observations available on the model grid.
!>
!>The routines are called by the different call-back routines of PDAF.
!>Most of the routines are generic so that in practice only 2 routines
!>need to be adapted for a particular data type. These are the routines
!>for the initialization of the observation information (`init_dim_obs`)
!>and for the observation operator (`obs_op`).
!>
!>
!>The module uses two derived data type (obs_f and obs_l), which contain
!>all information about the full and local observations. Only variables
!>of the type obs_f need to be initialized in this module. The variables
!>in the type obs_l are initialized by the generic routines from `PDAFomi`.
!>
MODULE mod_obs_ssh_mgrid_pdafomi

   USE mod_kind_pdaf
   USE mod_parallel_pdaf, &
      ONLY: mype_filter, abort_parallel
   USE PDAFomi, &
      ONLY: obs_f, obs_l
   USE netcdf

   IMPLICIT NONE
   SAVE

   !> Whether to assimilate this data type
   LOGICAL :: assim_ssh_mgrid = .TRUE.
   !> Observation error standard deviation (for constant errors)
   REAL(pwp) :: rms_ssh_mgrid = 1
   !> Whether to perform an identical twin experiment
   LOGICAL :: twin_exp_ssh_mgrid = .FALSE.
   !> Standard deviation for Gaussian noise in twin experiment
   REAL(pwp) :: noise_amp_ssh_mgrid = 1
   ! NetCDF file holding observations
   CHARACTER(lc) :: file_ssh_mgrid = 'my_nemo_ssh_file.nc'

   !> Instance of full observation data type - see `PDAFomi` for details.
   TYPE(obs_f), TARGET, PUBLIC :: thisobs
   !> Instance of local observation data type - see `PDAFomi` for details.
   TYPE(obs_l), TARGET, PUBLIC :: thisobs_l

!$OMP THREADPRIVATE(thisobs_l)

CONTAINS

   !>###Initialize information on the observation
   !>
   !> The routine is called by each filter process.
   !> at the beginning of the analysis step before
   !> the loop through all local analysis domains.
   !>
   !> It has to count the number of process-local and full
   !> observations, initialize the vector of observations
   !> and their inverse variances, initialize the coordinate
   !> array and index array for indices of observed elements
   !> of the state vector.
   !>
   !> The following four variables have to be initialized in this routine: <br/>
   !> **thisobs%doassim** - Whether to assimilate ssh <br/>
   !> **thisobs%disttype** - type of distance computation for localization
   !> with ssh <br/>
   !> **thisobs%ncoord** - number of coordinates used for distance
   !> computation <br/>
   !> **thisobs%id_obs_p** - index of module-type observation in PE-local state
   !> vector <br/>
   !>
   !> Optional is the use of: <br/>
   !> **thisobs%icoeff_p** - Interpolation coefficients for obs. operator
   !> (only if interpolation is used) <br/>
   !> **thisobs%domainsize** - Size of domain for periodicity for *disttype=1*
   !> (<0 for no periodicity) <br/>
   !> **thisobs%obs_err_type** - Type of observation errors for particle filter
   !> and NETF <br/>
   !> **thisobs%use_global_obs** - Whether to use global observations or
   !> restrict the observations to the relevant ones (default: *.true.* i.e use
   !> global full observations) <br/>
   !>
   !> The following variables are set in the routine gather_obs:
   !> **thisobs%dim_obs_p** - PE-local number of ssh observations <br/>
   !> **thisobs%dim_obs** - full number of ssh observations <br/>
   !> **thisobs%obs_f** - full vector of ssh observations <br/>
   !> **thisobs%ocoord_f** - coordinates of observations in OBS_MOD_F <br/>
   !> **thisobs%ivar_obs_f** - full vector of inverse obs. error variances of
   !> module-type <br/>
   !> **thisobs%dim_obs_g** - Number of global observations (only if
   !> *use_global_obs=.false*) <br/>
   !> **thisobs%id_obs_f_lim** - Ids of full observations in global observations
   !> (if *use_global_obs=.false*) <br/>
   !>
   SUBROUTINE init_dim_obs_ssh_mgrid(step, dim_obs)

      USE PDAFomi, &
         ONLY: PDAFomi_gather_obs
      USE mod_assimilation_pdaf, &
         ONLY: filtertype, local_range, delt_obs
      USE mod_statevector_pdaf, &
         ONLY: mpi_subd_lon, mpi_subd_lat
      USE mod_parallel_pdaf, &
         ONLY: COMM_filter
      USE par_oce, &
         ONLY: jpiglo, jpjglo
      USE dom_oce, &
         ONLY: nldj, nldi, glamt, gphit, nimpp, njmpp, ndastp

      !> Current time step
      INTEGER, INTENT(in)    :: step
      !> Dimension of full observation vector
      INTEGER, INTENT(inout) :: dim_obs

      !> Counters
      INTEGER :: i, j, s
      !> Step for observations in NetCDF file
      INTEGER :: nc_step = 0
      !> Status array for NetCDF operations
      INTEGER :: stat(50)
      !> ID for NetCDF file
      INTEGER :: ncid_in
      !> IDs for fields
      INTEGER :: id_var
      !> NetCDF position arrays for 3D field
      INTEGER :: pos(3), cnt(3)
      !> Global observation field
      REAL(pwp), ALLOCATABLE :: obs(:, :, :)

      !> Number of process-local observations
      INTEGER :: dim_obs_p
      !> Counters
      INTEGER :: cnt_p, cnt0_p
      !> Global gridbox coordinates of observations
      INTEGER :: i_obs, j_obs
      !> Halo offset for local PE
      INTEGER :: i0, j0
      !> PE-local observation vector
      REAL(pwp), ALLOCATABLE :: obs_p(:)
      !> PE-local inverse observation error variance
      REAL(pwp), ALLOCATABLE :: ivar_obs_p(:)
      !> PE-local observation coordinates
      REAL(pwp), ALLOCATABLE :: ocoord_p(:, :)
      !> Degree to radian conversion
      REAL(pwp) :: rad_conv = 3.141592653589793/180.0

      ! *****************************
      ! *** Global setting config ***
      ! *****************************

      IF (mype_filter == 0) &
         WRITE (*, '(8x,a)') 'Assimilate observations - obs_ssh_mgrid'

      ! Store whether to assimilate this observation type (used in routines
      ! below)
      IF (assim_ssh_mgrid) thisobs%doassim = 1

      ! Specify type of distance computation
      thisobs%disttype = 3   ! 3=Haversine

      ! Number of coordinates used for distance computation.
      ! The distance compution starts from the first row
      thisobs%ncoord = 2

      ! **********************************
      ! *** Read PE-local observations ***
      ! **********************************

      ! Format of ndastp is YYYYMMDD
      IF (mype_filter == 0) WRITE (*, '(/9x, a, i8)') &
         'obs_ssh_mgrid current date:', ndastp

      s = 1
      stat(s) = NF90_OPEN(file_ssh_mgrid, NF90_NOWRITE, ncid_in)
      s = s + 1

      stat(s) = NF90_INQ_VARID(ncid_in, 'sossheig', id_var)
      s = s + 1

      ALLOCATE (obs(jpiglo, jpjglo, 1))
      ! Increment time in NetCDF file so correct obs read
      nc_step = nc_step + delt_obs
      pos = (/1, 1, nc_step/)
      cnt = (/jpiglo, jpjglo, 1/)

      stat(s) = NF90_GET_VAR(ncid_in, id_var, obs, start=pos, count=cnt)
      s = s + 1

      stat(s) = NF90_CLOSE(ncid_in)
      s = s + 1

      DO j = 1, s - 1
         IF (stat(j) .NE. NF90_NOERR) THEN
            WRITE (*, '(/9x, a, 3x, a, 3x, a, i2)') &
               'NetCDF error in reading obs file:', file_ssh_mgrid, &
               'status array, j=', j
            CALL abort_parallel()
         END IF
      END DO

      ! ***********************************************************
      ! *** Count available observations for the process domain ***
      ! *** and initialize index and coordinate arrays.         ***
      ! ***********************************************************

      ! Compute halo offset
      i0 = nldi - 1
      j0 = nldj - 1

      cnt_p = 0

      DO j = 1, mpi_subd_lat
         DO i = 1, mpi_subd_lon
            ! Convert to global coordinates
            i_obs = nimpp + i0 + i - 1
            j_obs = njmpp + j0 + j - 1
            cnt_p = cnt_p + 1
         END DO
      END DO

      ! Set number of local observations
      dim_obs_p = cnt_p

      IF (cnt_p == 0) WRITE (*, '(/9x, a, i3, 3x, a, i4)') &
         'WARNING: No ssh_mgrid observations on PE:', mype_filter, &
         'NetCDF file step=', nc_step

      obs_nonzero: IF (dim_obs_p > 0) THEN
         ! Vector of observations on the process sub-domain
         ALLOCATE (obs_p(dim_obs_p))
         ! Coordinate array of observations on the process sub-domain
         ALLOCATE (ocoord_p(2, dim_obs_p))
         ! Coordinate array for observation operator
         ALLOCATE (thisobs%id_obs_p(1, dim_obs_p))
         ALLOCATE (ivar_obs_p(dim_obs_p))

         ! Compute halo offset
         i0 = nldi - 1
         j0 = nldj - 1

         cnt_p = 0
         cnt0_p = 0

         DO j = 1, mpi_subd_lat
            DO i = 1, mpi_subd_lon
               ! State vector index counter for observation operator.
               cnt0_p = cnt0_p + 1

               ! Convert to global coordinates.
               i_obs = nimpp + i0 + i - 1
               j_obs = njmpp + j0 + j - 1

               cnt_p = cnt_p + 1
               obs_p(cnt_p) = obs(i_obs, j_obs, 1)

               ! Observation coordinates - must be in radians for PDAFOMI
               ocoord_p(1, cnt_p) = glamt(i + i0, j + j0)*rad_conv
               ocoord_p(2, cnt_p) = gphit(i + i0, j + j0)*rad_conv

               ! Coordinates for observation operator (gridpoint)
               thisobs%id_obs_p(1, cnt_p) = cnt0_p
            END DO
         END DO
      ELSE
         ! No observations on PE, create dummy arrays to pass to PDAFOMI
         ALLOCATE (obs_p(1))
         ALLOCATE (ivar_obs_p(1))
         ALLOCATE (ocoord_p(2, 1))
         ALLOCATE (thisobs%id_obs_p(1, 1))
         obs_p = -999999.0
         ivar_obs_p = EPSILON(ivar_obs_p)
         ocoord_p = 0
         thisobs%id_obs_p = 1
      END IF obs_nonzero

      ! ****************************************************************
      ! *** Define observation errors for process-local observations ***
      ! ****************************************************************

      ! Set inverse observation error variances
      ivar_obs_p(:) = 1.0/(rms_ssh_mgrid*rms_ssh_mgrid)

      ! *********************************************************
      ! *** For twin experiment: Read synthetic observations  ***
      ! *********************************************************

      IF (twin_exp_ssh_mgrid) THEN
         IF (dim_obs_p > 0) CALL add_noise(dim_obs_p, obs_p)
      END IF

      ! ****************************************
      ! *** Gather global observation arrays ***
      ! ****************************************

      CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                              thisobs%ncoord, local_range, dim_obs)

      ! ********************
      ! *** Finishing up ***
      ! ********************

      ! Deallocate all local arrays
      DEALLOCATE (obs)
      DEALLOCATE (obs_p, ocoord_p, ivar_obs_p)

      ! Arrays in THISOBS have to be deallocated after the analysis step
      ! by a call to deallocate_obs() in prepoststep_pdaf.

   END SUBROUTINE init_dim_obs_ssh_mgrid

   !>###Implementation of observation operator
   !>
   !>This routine applies the full observation operator
   !>for the ssh observations.
   !>
   !>The routine is called by all filter processes.
   !>
   SUBROUTINE obs_op_ssh_mgrid(dim_p, dim_obs, state_p, ostate)

      USE PDAFomi, &
         ONLY: PDAFomi_obs_op_gridpoint

      !> PE-local state dimension
      INTEGER, INTENT(in) :: dim_p
      !> Dimension of full observed state (all observed fields)
      INTEGER, INTENT(in) :: dim_obs
      !> PE-local model state
      REAL(pwp), INTENT(in) :: state_p(dim_p)
      !> Full observed state
      REAL(pwp), INTENT(inout) :: ostate(dim_obs)

      ! ******************************************************
      ! *** Apply observation operator H on a state vector ***
      ! ******************************************************

      IF (thisobs%doassim == 1) THEN
         CALL PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
      END IF

   END SUBROUTINE obs_op_ssh_mgrid

   !>###Initialize local information on the module-type observation
   !>
   !>The routine is called during the loop over all local
   !>analysis domains. It has to initialize the information
   !>about local ssh observations.
   !>
   !>This routine calls the routine `PDAFomi_init_dim_obs_l`
   !>for each observation type. The call allows to specify a
   !>different localization radius and localization functions
   !>for each observation type and local analysis domain.
   !>
   SUBROUTINE init_dim_obs_l_ssh_mgrid(domain_p, step, dim_obs, dim_obs_l)

      USE PDAFomi, &
         ONLY: PDAFomi_init_dim_obs_l
      USE mod_assimilation_pdaf, &
         ONLY: coords_l, local_range, locweight, srange

      !> Index of current local analysis domain
      INTEGER, INTENT(in)  :: domain_p
      !> Current time step
      INTEGER, INTENT(in)  :: step
      !> Full dimension of observation vector
      INTEGER, INTENT(in)  :: dim_obs
      !> Local dimension of observation vector
      INTEGER, INTENT(out) :: dim_obs_l

      ! **********************************************
      ! *** Initialize local observation dimension ***
      ! **********************************************

      CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
                                  locweight, local_range, srange, dim_obs_l)

   END SUBROUTINE init_dim_obs_l_ssh_mgrid

   !>###Routine to add model error.
   !>
   SUBROUTINE add_noise(dim_obs_p, obs)

      !> Number of process-local observations
      INTEGER, INTENT(in) :: dim_obs_p
      !> Process-local observations
      REAL, INTENT(inout) :: obs(dim_obs_p)

      !> Random noise
      REAL, ALLOCATABLE :: noise(:)
      !> Seed for random number generator
      INTEGER, SAVE :: iseed(4)
      !> Flag for first call
      LOGICAL, SAVE :: firststep = .TRUE.

      ! Seeds taken from PDAF Lorenz96 routine
      IF (firststep) THEN
         WRITE (*, '(9x, a)') '--- Initialize seed for ssh_mgrid noise'
         iseed(1) = 2*220 + 1
         iseed(2) = 2*100 + 5
         iseed(3) = 2*10 + 7
         iseed(4) = 2*30 + 9
         firststep = .FALSE.
      END IF

      ! Generate random Gaussian noise
      ALLOCATE (noise(dim_obs_p))
      CALL dlarnv(3, iseed, dim_obs_p, noise)

      obs = obs + (noise_amp_ssh_mgrid*noise)

      DEALLOCATE (noise)

   END SUBROUTINE add_noise

END MODULE mod_obs_ssh_mgrid_pdafomi
