!>##Set number of local analysis domains
!>The routine is called in `PDAF_X_update`
!>at the beginning of the analysis step before
!>the loop through all local analysis domains.
!>It has to set the number of local analysis
!>domains for the PE-local domain.
!>
!> - Called from: `PDAFomi_assimilate_local`/`mod_assimilation_pdaf`
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)

   USE mod_kind_pdaf
   USE mod_assimilation_pdaf, &
      ONLY: indx_dom_l
   USE mod_statevector_pdaf, &
      ONLY: mpi_subd_lon, mpi_subd_lat
   USE mod_parallel_pdaf, &
      ONLY: mype_filter
   USE dom_oce, &
      ONLY: nldi, nldj, tmask

   IMPLICIT NONE

   !> Current time step
   INTEGER, INTENT(in)  :: step
   !> PE-local number of analysis domains
   INTEGER, INTENT(out) :: n_domains_p

   !> Counters
   INTEGER :: i, j, cnt
   !> Halo offset
   INTEGER :: i0, j0

   ! ************************************
   ! *** Initialize number of domains ***
   ! ************************************

   ! *******************************************
   !
   ! The number of local domains is defined as
   ! the number of grid points at the surface
   ! where tmask is 1 ie horizontal localization
   ! is used, and land points are ignored.
   !
   ! *******************************************

   n_domains_p = 0

   ! Compute halo offset
   i0 = nldi - 1
   j0 = nldj - 1

   DO j = 1, mpi_subd_lat
      DO i = 1, mpi_subd_lon
         IF (tmask(i + i0, j + j0, 1) == 1.0_pwp) THEN
            n_domains_p = n_domains_p + 1
         END IF
      END DO
   END DO

   ! *******************************************
   ! *** Store local domain i,j index values ***
   ! *******************************************

   ! ******************************************
   !
   ! Important to consider case where there
   ! are no valid local domains. In this
   ! case we use a hack. We set the number
   ! of local domains (n_domains_p) to 1,
   ! and set the dimension of the observations
   ! for this local domain to zero (see
   ! init_dim_obs_l). This hack means no update
   ! will be performed on this local domain.
   !
   ! ******************************************

   IF (n_domains_p > 0) THEN
      ! Deallocated in deallocate_obs_pdafomi
      IF (.NOT. ALLOCATED(indx_dom_l)) ALLOCATE (indx_dom_l(2, n_domains_p))
      cnt = 0
      DO j = 1, mpi_subd_lat
         DO i = 1, mpi_subd_lon
            IF (tmask(i + i0, j + j0, 1) == 1.0_pwp) THEN
               cnt = cnt + 1
               indx_dom_l(1, cnt) = i
               indx_dom_l(2, cnt) = j
            END IF
         END DO
      END DO
   ELSE
      WRITE (*, '(8x,a,i3)') 'WARNING: No valid local domains, PE=', mype_filter
      n_domains_p = 1
      ! Deallocated in deallocate_obs_pdafomi
      IF (.NOT. ALLOCATED(indx_dom_l)) ALLOCATE (indx_dom_l(2, 1))
      indx_dom_l(1, 1) = 1
      indx_dom_l(2, 1) = 1
   END IF

END SUBROUTINE init_n_domains_pdaf
