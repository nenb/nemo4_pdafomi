!>##Building the Statevector
!>This module provides variables & routines for
!>building the state vector.
!>
MODULE mod_statevector_pdaf

   USE mod_kind_pdaf

   IMPLICIT NONE
   SAVE

   !> 2d statevector variables - start index
   INTEGER :: ssh_p_offset

   !> Array holding 2d state variable offsets
   INTEGER :: var2d_p_offset(1)

   !> 2d statevector variables - dimension size
   INTEGER :: ssh_p_dim

   !> 3d statevector variables - start index for t
   INTEGER :: t_p_offset
   !> 3d statevector variable - start index for s
   INTEGER :: s_p_offset
   !> 3d statevector variable - start index for u
   INTEGER :: u_p_offset
   !> 3d statevector variable - start index for v
   INTEGER :: v_p_offset

   !> Array holding 3d state variable offsets
   !> index i is the i-th variable
   INTEGER :: var3d_p_offset(4)

   !> 3d statevector variables (t) - dimension size
   INTEGER :: t_p_dim
   !> 3d statevector variables (s) - dimension size
   INTEGER :: s_p_dim
   !> 3d statevector variables (u) - dimension size
   INTEGER :: u_p_dim
   !> 3d statevector variables (v) - dimension size
   INTEGER :: v_p_dim

   !> Dimensions for MPI subdomain that is included
   !> in local statevector. Necessary so that halo
   !> regions are not included in multiple local
   !> statevectors

   !> size of local lat domain excluding halo region
   INTEGER :: mpi_subd_lat
   !> size of local lon domain excluding halo region
   INTEGER :: mpi_subd_lon
   !> size of local vertical domain
   INTEGER :: mpi_subd_vert

CONTAINS

   !>##This routine calculates the dimensions of the MPI subdomain
   !> that is used to fill the local statevector.
   !>
   !> **Calling Sequence**
   !>
   !> - Called from: `calc_statevar_dim`
   SUBROUTINE calc_mpi_dim()

      USE par_oce, ONLY: jpk
      USE dom_oce, ONLY: nldi, nldj, nlei, nlej

      mpi_subd_lon = nlei - nldi + 1
      mpi_subd_lat = nlej - nldj + 1
      mpi_subd_vert = jpk

   END SUBROUTINE calc_mpi_dim

   !>##This routine calculates the dimension of each of the local
   !> statevector variables.
   !>
   !> **Calling Sequence**
   !>
   !> - Called from: `calc_offset`
   !>
   !> - Calls: `calc_mpi_dim`
   SUBROUTINE calc_statevar_dim()

      ! Compute MPI subdomain dimensions
      CALL calc_mpi_dim()

      ssh_p_dim = mpi_subd_lat*mpi_subd_lon
      t_p_dim = mpi_subd_lat*mpi_subd_lon*mpi_subd_vert
      s_p_dim = mpi_subd_lat*mpi_subd_lon*mpi_subd_vert
      u_p_dim = mpi_subd_lat*mpi_subd_lon*mpi_subd_vert
      v_p_dim = mpi_subd_lat*mpi_subd_lon*mpi_subd_vert

   END SUBROUTINE calc_statevar_dim

   !>##This routine calculates the offset values for each of the
   !> local statevector variables.
   !>
   !> It then stores the 2d/3d offset values in separate arrays.
   !>
   !> **Calling Sequence**
   !>
   !> - Called from: `calc_statevector_dim`
   !>
   !> - Calls: `calc_statevar_dim`
   SUBROUTINE calc_offset()

      ! Compute local statevector dimensions
      CALL calc_statevar_dim()

      ssh_p_offset = 0
      t_p_offset = ssh_p_offset + ssh_p_dim
      s_p_offset = t_p_offset + t_p_dim
      u_p_offset = s_p_offset + s_p_dim
      v_p_offset = u_p_offset + u_p_dim

      ! Fill array of 2D state variable offsets for local PE
      var2d_p_offset(1) = ssh_p_offset

      ! Fill array of 3D state variable offsets for local PE
      var3d_p_offset(1) = t_p_offset
      var3d_p_offset(2) = s_p_offset
      var3d_p_offset(3) = u_p_offset
      var3d_p_offset(4) = v_p_offset

   END SUBROUTINE calc_offset

   !>##This routine calculates the dimension of the local statevector.
   !>
   !> **Calling Sequence**
   !>
   !> - Called from: `init_pdaf`
   !>
   !> - Calls: `calc_offset`
   SUBROUTINE calc_statevector_dim(dim_p)

      !> Local statevector dimension
      INTEGER, INTENT(inout) :: dim_p

      ! Calculate statevector variable offset and dimension.
      ! *DO NOT REMOVE* as offset is not calculated anywhere else.
      CALL calc_offset()

      dim_p = ssh_p_dim + t_p_dim + s_p_dim + u_p_dim + v_p_dim

   END SUBROUTINE calc_statevector_dim

   !>##Fill local ensemble array with 2d state variables from
   !> initial state file.
   !>
   !> @todo
   !> User should give values to the statevector
   !> for the DA.
   !> @endtodo
   SUBROUTINE fill2d_ensarray(fname, ens_p)

      USE netcdf

      !> Name of netCDF file
      CHARACTER(lc), INTENT(in) :: fname
      !> PE-local state ensemble
      REAL(pwp), INTENT(inout) :: ens_p(:, :)

      WRITE (*, '(/1x,a,a)') '2D initial state file =', TRIM(fname)
      WRITE (*, '(/1x,a)') 'Insert routine for reading 2D initial state here.'

   END SUBROUTINE fill2d_ensarray

   !>##Fill local ensemble array with 3d state variables from
   !> initial state file.
   !>
   !> @todo
   !> User should give values to the statevector
   !> for the DA.
   !> @endtodo
   SUBROUTINE fill3d_ensarray(fname, statevar, ens_p)

      USE netcdf

      !> Name of netCDF file
      CHARACTER(lc), INTENT(in) :: fname
      !> Name of state variable
      CHARACTER(len=1), INTENT(in) :: statevar
      !> PE-local state ensemble
      REAL(pwp), INTENT(inout) :: ens_p(:, :)

      SELECT CASE (statevar)
      CASE ('T')
         WRITE (*, '(/1x,a,a)') 'T initial state file =', TRIM(fname)
         WRITE (*, '(/1x,a)') 'Insert routine for reading T initial state here.'
      CASE ('S')
         WRITE (*, '(/1x,a,a)') 'S initial state file =', TRIM(fname)
         WRITE (*, '(/1x,a)') 'Insert routine for reading S initial state here.'
      CASE ('U')
         WRITE (*, '(/1x,a,a)') 'U initial state file =', TRIM(fname)
         WRITE (*, '(/1x,a)') 'Insert routine for reading U initial state here.'
      CASE ('V')
         WRITE (*, '(/1x,a,a)') 'V initial state file =', TRIM(fname)
         WRITE (*, '(/1x,a)') 'Insert routine for reading V initial state here.'
      END SELECT

   END SUBROUTINE fill3d_ensarray

END MODULE mod_statevector_pdaf
