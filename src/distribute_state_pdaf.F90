!>##Distributing the statevector variables
!>The routine has to initialize the fields of the
!>model from the statevector of PDAF.
!>
!>The routine is executed by each process that is
!>participating in the model integrations.
!>
!>**Calling Sequence**
!>*Called from:* `PDAF_get_state` (as U_dist_state)
!>*Called from:* `PDAFomi_assimilate_local` (as U_dist_state)
!>
SUBROUTINE distribute_state_pdaf(dim_p, state_p)

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
   USE lbclnk, &
      ONLY: lbc_lnk, lbc_lnk_multi

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

   ! ************************************
   ! Distribute state vector 2d variables
   ! ************************************
   ! SSH
   DO j = 1, mpi_subd_lat
      DO i = 1, mpi_subd_lon
         sshb(i + i0, j + j0) = state_p(i + (j - 1)*mpi_subd_lon &
                                        + ssh_p_offset)
      END DO
   END DO

   ! Fill halo regions
   CALL lbc_lnk('distribute_state_pdaf', sshb, 'T', 1.)

   ! ************************************
   ! Distribute state vector 3d variables
   ! ************************************
   ! T
   DO k = 1, mpi_subd_vert
      DO j = 1, mpi_subd_lat
         DO i = 1, mpi_subd_lon
            tsb(i + i0, j + j0, k, jp_tem) = state_p(i + (j - 1)*mpi_subd_lon &
                                                     + (k - 1)*mpi_subd_lat*mpi_subd_lon + t_p_offset)
         END DO
      END DO
   END DO

   ! S
   DO k = 1, mpi_subd_vert
      DO j = 1, mpi_subd_lat
         DO i = 1, mpi_subd_lon
            tsb(i + i0, j + j0, k, jp_sal) = state_p(i + (j - 1)*mpi_subd_lon &
                                                     + (k - 1)*mpi_subd_lat*mpi_subd_lon + s_p_offset)
         END DO
      END DO
   END DO

   ! Fill halo regions
   CALL lbc_lnk_multi('distribute_state_pdaf', tsb(:, :, :, jp_tem), 'T', &
                      1., tsb(:, :, :, jp_sal), 'T', 1.)

   ! U
   DO k = 1, mpi_subd_vert
      DO j = 1, mpi_subd_lat
         DO i = 1, mpi_subd_lon
            ub(i + i0, j + j0, k) = state_p(i + (j - 1)*mpi_subd_lon + &
                                            (k - 1)*mpi_subd_lat*mpi_subd_lon + u_p_offset)
         END DO
      END DO
   END DO

   ! V
   DO k = 1, mpi_subd_vert
      DO j = 1, mpi_subd_lat
         DO i = 1, mpi_subd_lon
            vb(i + i0, j + j0, k) = state_p(i + (j - 1)*mpi_subd_lon + &
                                            (k - 1)*mpi_subd_lat*mpi_subd_lon + v_p_offset)
         END DO
      END DO
   END DO

   ! Fill halo regions
   CALL lbc_lnk_multi('distribute_state_pdaf', ub, 'U', -1., vb, 'V', -1.)

END SUBROUTINE distribute_state_pdaf
