!>##Collecting the statevector variables
!>The routine has to initialize the statevector of PDAF
!>from the fields of the model.
!>
!>The routine is executed by each process that is
!>participating in the model integrations.
!>
!>**Calling Sequence**
!>
!> - Called from:* `PDAFomi_assimilate_local`/`mod_assimilation_pdaf` (as U_coll_state)
!>
SUBROUTINE collect_state_pdaf(dim_p, state_p)

   USE mod_kind_pdaf
   USE mod_statevector_pdaf, &
      ONLY: mpi_subd_lat, mpi_subd_lon, mpi_subd_vert, ssh_p_offset, &
            t_p_offset, s_p_offset, u_p_offset, v_p_offset
   USE dom_oce, &
      ONLY: nldj, nldi
   USE oce, &
      ONLY: sshb, tsb, ub, vb
   USE par_oce, &
      ONLY: jp_tem, jp_sal, jpi, jpj, jpk

   IMPLICIT NONE

   !> PE-local state dimension
   INTEGER, INTENT(in) :: dim_p
   !> PE-local state vector
   REAL(pwp), INTENT(inout) :: state_p(dim_p)

   !> Counters
   INTEGER :: i, j, k
   !> Start index for MPI subdomain
   INTEGER :: i0, j0

   ! Set the starting index after the halo region
   j0 = nldj - 1
   i0 = nldi - 1

   ! *********************************
   ! Collect state vector 2d variables
   ! *********************************
   ! SSH
   DO j = 1, mpi_subd_lat
      DO i = 1, mpi_subd_lon
         state_p(i + (j - 1)*mpi_subd_lon + ssh_p_offset) = &
            sshb(i + i0, j + j0)
      END DO
   END DO

   ! *********************************
   ! Collect state vector 3d variables
   ! *********************************
   ! T
   DO k = 1, mpi_subd_vert
      DO j = 1, mpi_subd_lat
         DO i = 1, mpi_subd_lon
            state_p(i + (j - 1)*mpi_subd_lon + &
                    (k - 1)*mpi_subd_lat*mpi_subd_lon + t_p_offset) = tsb(i + i0, &
                                                                          j + j0, k, jp_tem)
         END DO
      END DO
   END DO

   ! S
   DO k = 1, mpi_subd_vert
      DO j = 1, mpi_subd_lat
         DO i = 1, mpi_subd_lon
            state_p(i + (j - 1)*mpi_subd_lon + &
                    (k - 1)*mpi_subd_lat*mpi_subd_lon + s_p_offset) = tsb(i + i0, &
                                                                          j + j0, k, jp_sal)
         END DO
      END DO
   END DO

   ! U
   DO k = 1, mpi_subd_vert
      DO j = 1, mpi_subd_lat
         DO i = 1, mpi_subd_lon
            state_p(i + (j - 1)*mpi_subd_lon + &
                    (k - 1)*mpi_subd_lat*mpi_subd_lon + u_p_offset) = ub(i + i0, &
                                                                         j + j0, k)
         END DO
      END DO
   END DO

   ! V
   DO k = 1, mpi_subd_vert
      DO j = 1, mpi_subd_lat
         DO i = 1, mpi_subd_lon
            state_p(i + (j - 1)*mpi_subd_lon + &
                    (k - 1)*mpi_subd_lat*mpi_subd_lon + v_p_offset) = vb(i + i0, &
                                                                         j + j0, k)
         END DO
      END DO
   END DO

END SUBROUTINE collect_state_pdaf
