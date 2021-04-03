!>##Set dimension of local model state
!>The routine is called during analysis step
!>in `PDAF_X_update` in the loop over all local
!>analysis domains. It has to set the dimension
!>of the local model state on the current analysis
!>domain.
!>
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

   USE mod_kind_pdaf
   USE mod_assimilation_pdaf, &
      ONLY: coords_l, indx_dom_l
   USE mod_statevector_pdaf, &
      ONLY: mpi_subd_vert, var2d_p_offset, var3d_p_offset
   USE mod_parallel_pdaf, &
      ONLY: abort_parallel
   USE dom_oce, &
      ONLY: nldj, nldi, tmask, glamt, gphit

   IMPLICIT NONE

   !> Current time step
   INTEGER, INTENT(in)  :: step
   !> Current local analysis domain
   INTEGER, INTENT(in)  :: domain_p
   !> Local state dimension
   INTEGER, INTENT(out) :: dim_l

   !> Counters
   INTEGER :: idx, cnt
   !> Grid coordinates for local analysis domain
   INTEGER :: i, j
   !> Halo offset for local PE
   INTEGER :: i0, j0
   !> Longitude, latitude for local analysis domain
   REAL(pwp) :: lon, lat
   !> Degree to radian conversion
   REAL(pwp) :: rad_conv = 3.141592653589793/180.0

   ! Compute halo offset
   i0 = nldi - 1
   j0 = nldj - 1

   ! Compute coordinates
   i = indx_dom_l(1, domain_p)
   j = indx_dom_l(2, domain_p)

   ! **********************************************
   ! *** Initialize coordinates of local domain ***
   ! **********************************************

   ! Use T-values to get local coordinates
   lat = gphit(i + i0, j + j0)
   lon = glamt(i + i0, j + j0)

   ! Convert local domain coordinates to radians (as required by PDAFOMI)
   coords_l(1) = lon*rad_conv
   coords_l(2) = lat*rad_conv

   ! ****************************************
   ! *** Initialize local state dimension ***
   ! ****************************************

   ! *************************************************************
   ! dimension = (number of 2D state variables) + (number of 3D
   ! variables * number of ocean vertical points).
   !
   ! We need to calculate the number of ocean vertical points ie
   ! we need to determine how many points in the vertical are
   ! ocean and how may are land (we do not include land points in
   ! our local state vector).
   ! *************************************************************

   cnt = 0

   DO idx = 1, mpi_subd_vert
      IF (tmask(i + i0, j + j0, idx) == 1.0_pwp) cnt = cnt + 1
   END DO

   dim_l = SIZE(var2d_p_offset) + (SIZE(var3d_p_offset)*cnt)

END SUBROUTINE init_dim_l_pdaf
