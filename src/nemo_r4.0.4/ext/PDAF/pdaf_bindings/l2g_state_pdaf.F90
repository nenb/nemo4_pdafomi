!>##Initialize full state from local analysis
!>
!>The routine is called during the loop over all
!>local analysis domains in `PDAF_X_update`
!>after the analysis and ensemble transformation
!>on a single local analysis domain. It has to
!>initialize elements of the PE-local full state
!>vector from the provided analysis state vector
!>on the local analysis domain.
!>
SUBROUTINE l2g_state_pdaf(step, domain_p, dim_l, state_l, dim_p, state_p)

   USE mod_kind_pdaf
   USE mod_parallel_pdaf, &
      ONLY: abort_parallel
   USE mod_assimilation_pdaf, &
      ONLY: indx_dom_l
   USE mod_statevector_pdaf, &
      ONLY: mpi_subd_vert, mpi_subd_lon, mpi_subd_lat, &
            var2d_p_offset, var3d_p_offset
   USE dom_oce, &
      ONLY: nldj, nldi, tmask

   IMPLICIT NONE

   !> Current time step
   INTEGER, INTENT(in) :: step
   !> Current local analysis domain
   INTEGER, INTENT(in) :: domain_p
   !> PE-local full state dimension
   INTEGER, INTENT(in) :: dim_l
   !> Local state dimension
   INTEGER, INTENT(in) :: dim_p
   !> PE-local full state vector
   REAL(pwp), INTENT(in)    :: state_l(dim_l)
   !> State vector on local analysis domain
   REAL(pwp), INTENT(out)   :: state_p(dim_p)

   !> Counter
   INTEGER :: idx
   !> Halo offset for local PE
   INTEGER :: i0, j0
   !> Grid coordinates for local analysis domain
   INTEGER :: i, j
   !> Array for converting vertical coordinate in local statevector
   INTEGER :: vcoord_conv(100)
   !> Number of vertical ocean points in local statevector
   INTEGER :: dim_vert_l
   !> 2D state variable coordinate in statevector
   INTEGER :: loc_2dvar
   !> Variables for 3D state variable index
   INTEGER :: a, b, c

   ! *******************************************************
   ! *** Initialise vertical coordinate conversion array ***
   ! *******************************************************

   IF (SIZE(vcoord_conv) < mpi_subd_vert) THEN
      WRITE (*, '(/1x,a59/)') &
         'ERROR: automatic array v_coord in l2g_state_pdaf too small '
      CALL abort_parallel()
   END IF

   vcoord_conv = 0
   dim_vert_l = 0

   ! Compute halo offset
   i0 = nldi - 1
   j0 = nldj - 1

   ! Compute coordinates
   i = indx_dom_l(1, domain_p)
   j = indx_dom_l(2, domain_p)

   ! Count number of ocean vertical points
   DO idx = 1, mpi_subd_vert
      IF (tmask(i + i0, j + j0, idx) == 1.0_pwp) THEN
         dim_vert_l = dim_vert_l + 1
         vcoord_conv(dim_vert_l) = idx
      END IF
   END DO

   ! *************************************
   ! *** Initialize local state vector ***
   ! *************************************

   ! **********************************************************
   ! A local domain consists of all ocean points in a vertical
   ! column. Such a domain will have coordinates (x,y,:).
   ! The 2d state variables in the global statevector will be
   ! located at ( (y-1)*dim_longitude ) + x + 2d_offset.
   ! The 3d state variables in the global statevector will be
   ! located at ( (z-1)*dim_longitude*dim_latitude ) +
   ! ( (y-1)*dim_longitude ) + x + 3d_offset, where z can vary
   ! over all *ocean* points in the vertical column.
   ! **********************************************************

   ! Compute location of 2d variables in statevector
   loc_2dvar = (j - 1)*(mpi_subd_lon) + i

   DO idx = 1, dim_l
      ! Compute 2d state variables
      IF (idx <= SIZE(var2d_p_offset)) THEN
         state_p(loc_2dvar + var2d_p_offset(idx)) = state_l(idx)
      ELSE ! Compute 3d state variables
         a = idx - SIZE(var2d_p_offset)
         b = MOD(a - 1, dim_vert_l) + 1
         c = (a - b)/dim_vert_l + 1
         state_p(mpi_subd_lon*mpi_subd_lat*(vcoord_conv(b) - 1) &
                 + loc_2dvar + var3d_p_offset(c)) = state_l(idx)
      END IF
   END DO

END SUBROUTINE l2g_state_pdaf
