!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This module contains the netCDF include file and a netcdf error
!     handling routine.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: netcdf.f,v 1.2 2000/04/19 21:56:25 pwjones Exp $
!
!     Copyright (c) 1997, 1998 the Regents of the University of 
!       California.
!
!     This software and ancillary information (herein called software) 
!     called SCRIP is made available under the terms described here.  
!     The software has been approved for release with associated 
!     LA-CC Number 98-45.
!
!     Unless otherwise indicated, this software has been authored
!     by an employee or employees of the University of California,
!     operator of the Los Alamos National Laboratory under Contract
!     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
!     Government has rights to use, reproduce, and distribute this
!     software.  The public may copy and use this software without
!     charge, provided that this Notice and any statement of authorship
!     are reproduced on all copies.  Neither the Government nor the
!     University makes any warranty, express or implied, or assumes
!     any liability or responsibility for the use of this software.
!
!     If software is modified to produce derivative works, such modified
!     software should be clearly marked, so as not to confuse it with 
!     the version available from Los Alamos National Laboratory.
!
!***********************************************************************

      module netcdf_mod

!-----------------------------------------------------------------------

      use kinds_mod
      use constants

      implicit none

      include 'netcdf.inc'

!***********************************************************************

      contains

!***********************************************************************

      subroutine netcdf_error_handler(istat, mess)

!-----------------------------------------------------------------------
!
!     This routine provides a simple interface to netCDF error message
!     routine.
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(in) ::  &
          istat   ! integer status returned by netCDF function call
      character (len=*), intent(in), optional :: mess

!-----------------------------------------------------------------------

      if (istat /= NF_NOERR) then
       if (present(mess)) then
         print *,'Error in netCDF: ',nf_strerror(istat), 'Message: ',mess
       else
         print *,'Error in netCDF: ',nf_strerror(istat)
       endif
        stop
      endif

!-----------------------------------------------------------------------

      end subroutine netcdf_error_handler

!***********************************************************************

      end module netcdf_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
