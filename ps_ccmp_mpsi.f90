!>============================================================================<
!           FILE: PS_CCMP_MPSI.f90
!  SUBPROGRAM(S): SUBROUTINE PS_CCMP_MPSI
!     ORIGINALLY: SUBR_CCMP_ZM_MPSI
!         AUTHOR: Paul W. Staten, IU Bloomington
!   ADAPTED FROM: David Stepaniak, /NCAR/CGD/CAS
! DATE INITIATED: 11 November 1999
!  LAST MODIFIED: 
!
!    DESCRIPTION: Adapted from the zonal mean stream function calculator
!                 by Stepaniak. Compuates the stream function for
!                 each longitude for the divergent wind.
!                 The remainder of this description is original:
!
!                 Computes a ZONAL MEAN MERIDIONAL STREAM FUNCTION using a 
!                 modified definition of the CCM Processor zonal mean
!                 meridional stream function. The modified definition used
!                 here is given by 
!
!                                                           / PS
!                                                          |
!                                          2 pi a cos(lat) |  _
!                    ZM_MPSI(lat,lev)  =   --------------- |  V(lat,lev) dp
!                                                g         |
!                                                          |
!                                                         / p
!
!                       _
!                 where V is the zonal mean meridional wind, p pressure, PS
!                 surface pressure, a the radius of the earth, and g the
!                 acceleration of gravity. The coordinate dimensions are given
!                 by lon, the longitude dimension, lat, the latitude dimension,
!                 and lev, the pressure level dimension.
!
!      REFERENCE: Buja, L. E. (1994) CCM Processor User's Guide (Unicos
!                 Version). NCAR Technical Note NCAR/TN-384+IA, pages B-17
!                 to B-18.
!
!      INVOKE AS:
!
!             CALL PS_CCMP_MPSI(nlon, nlat, nlev, V, lat, p, PS, msg, ZM_MPSI )
!
!     INPUT ARGS:
!
!           nlon: Number of longitudes of longitude dimension. Type INTEGER.
!           nlat: Number of latitudes of latitude dimension. Type INTEGER.
!           nlev: Number of levels of pressure level dimension. Type 
!                 INTEGER.
!              V: Three-dimensional (lon,lat,lev) array of meridional wind
!                 values in which THE PRESSURE LEVEL DIMENSION MUST BE ORDERED
!                 TOP TO BOTTOM. Units must be m s^-1. Type REAL.
!                                            _ 
!                 (The zonal mean of V, i.e. V, will be computed in this sub-
!                 routine in allocated memory.)
!            lat: One-dimensional array of latitude values of latitude dimension
!                 in degrees. Type REAL.
!              p: One-dimensional array of pressure level values of vertical
!                 dimension ORDERED TOP TO BOTTOM. Units must be Pa. Type REAL.
!                 The first value must be greater than 50 Pa (0.5mb), and the
!                 last value must be less than 100500 Pa (1005mb).
!             PS: Two-dimensional (lon,lat) array of surface pressures. Units
!                 must be Pa. Type REAL.
!                                                                      
!            msg: Missing value (fill value) of stream function where V occurs
!                 entirely below the earth's surface for all longitudes at fixed
!                 latitude and pressure level. Type REAL.
!
!    OUTPUT ARGS: 
!
!        ZM_MPSI: Three-dimensional(lon,lat,lev) array of meridional stream
!                 function for divergent wind.
!                 Original documentation follows:
!                 Two-dimensional (lat,lev) array of zonal mean meridional
!                 stream function values in which the first dimension is lati-
!                 tude and the second dimension is pressure level. THE PRESSURE
!                 LEVEL DIMENSION IS ORDERED TOP TO BOTTOM UPON RETURN. Missing
!                 values (msg) are assigned to ZM_MPSI where V occurs below
!                 ground for all longitudes at fixed latitude and pressure level.
!                 Units are kg s^-1. Type REAL.
!
!>============================================================================<
! SUBROUTINE CCMP_ZM_MPSI( nlon, nlat, nlev, V, lat, p, PS, msg, ZM_MPSI )
  SUBROUTINE PS_CCMP_MPSI( nlon, nlat, nlev, V, lat, p, PS, msg, MPSI )

  IMPLICIT NONE

  INTEGER, INTENT(IN)                               :: nlon
  !f2py intent(hide),depend(V) :: nlon=shape(V,0)
  INTEGER, INTENT(IN)                               :: nlat
  !f2py intent(hide),depend(V) :: nlat=shape(V,1)
  INTEGER, INTENT(IN)                               :: nlev
  !f2py intent(hide),depend(V) :: nlev=shape(V,3)
  REAL, DIMENSION(1:nlon,1:nlat,1:nlev), INTENT(IN) :: V
  !f2py intent(in),dimension(nlon,nlat,nlev) :: V
  REAL, DIMENSION(1:nlat), INTENT(IN)               :: lat
  !f2py intent(in),dimension(nlat) :: lat
  REAL, DIMENSION(1:nlev), INTENT(IN)               :: p
  !f2py intent(in),dimension(nlev) :: p
  REAL, DIMENSION(1:nlon,1:nlat), INTENT(IN)        :: PS
  !f2py intent(in),dimension(nlon,nlat) :: PS
  REAL, INTENT(IN)                                  :: msg
  !f2py intent(in) :: msg

! REAL, DIMENSION(1:nlat,1:nlev), INTENT(OUT)       :: ZM_MPSI
  REAL, DIMENSION(1:nlon,1:nlat,1:nlev), INTENT(OUT):: MPSI
  !f2py intent(out),dimension(nlon,nlat,nlev) :: MPSI

!>============================================================================<

  REAL, PARAMETER                                   :: pi = 3.141593
  REAL, PARAMETER                                   ::  g = 9.80616
  REAL, PARAMETER                                   ::  a = 6.37122E6

!>----------------------------------------------------------------------------<

! REAL, PARAMETER                                   :: p_t = 500.
! 5 hPa is too low a ceiling IMO. pwstaten, 2017.07.28
  REAL, PARAMETER                                   :: p_t = 50.
  REAL, PARAMETER                                   :: p_o = 100500.
  REAL, DIMENSION(:), ALLOCATABLE                   :: ptmp
  REAL, DIMENSION(:), ALLOCATABLE                   :: dp
  REAL, DIMENSION(:,:,:), ALLOCATABLE               :: VTMP
  REAL, DIMENSION(:), ALLOCATABLE                   :: psitmp
  REAL, DIMENSION(:), ALLOCATABLE                   :: vvprof
  REAL                                              :: c
  INTEGER                                           :: ilon
  INTEGER                                           :: jlat
  INTEGER                                           :: klev
  INTEGER                                           :: flag

!>----------------------------------------------------------------------------<
! Check for valid range of pressure levels:

  IF ( ( ANY( p < p_t ) ) .OR. ( ANY( p > p_o ) ) ) THEN
   WRITE(*,'(A49)') "One or more pressure levels of vertical dimension"
   WRITE(*,'(A34)') "outside of range 50 to 100500 Pa."
   WRITE(*,'(A34)') "Execution stopped in subroutine CCMP_ZM_MPSI."
   STOP
  END IF

!>----------------------------------------------------------------------------<
! Assign pressure values at data levels, intermediate half levels, and
! pressure thicknesses at each data level:

  ALLOCATE( ptmp(0:2*nlev) )

  ptmp(1:2*nlev-1:2) = p(1:nlev)
                       ! data levels

  ptmp(2:2*nlev-2:2) = ( ptmp(3:2*nlev-1:2) + ptmp(1:2*nlev-3:2) ) / 2.
                       ! intemediate half levels

             ptmp(0) = p_t
                       ! pressure at top 

        ptmp(2*nlev) = p_o
                       ! pressure at bottom

  ALLOCATE( dp(0:2*nlev) )

    dp(1:2*nlev-1:2) = ptmp(2:2*nlev:2) - ptmp(0:2*nlev-2:2)
                       ! pressure thickness at each data level

      dp(0:2*nlev:2) = 0.
                       ! pressure thickness at intemediate hal levels,
                       ! unused

!>----------------------------------------------------------------------------<
! Compute zonal mean of V, taking into account below-ground grid points:
  
  ALLOCATE( VTMP(1:nlon,1:nlat,1:nlev) )

  VTMP(:,:,:) = V(:,:,:)

  DO ilon = 1, nlon

   DO jlat = 1, nlat

    WHERE( p(:) > PS(ilon,jlat) ) VTMP(ilon,jlat,:) = msg

   END DO

  END DO

! no need to calculate zonal mean for the divergent wind
!  ALLOCATE ( VBAR(1:nlat,1:nlev) )
!
!  DO jlat = 1, nlat
!
!   DO klev = 1, nlev
!
!    n = COUNT( VTMP(:,jlat,klev) /= msg )
!
!    IF ( n == 0 ) THEN
!
!      VBAR(jlat,klev) = msg
!
!    ELSE
!
!      VBAR(jlat,klev) = SUM( VTMP(:,jlat,klev), DIM = 1, &
!                                    MASK = VTMP(:,jlat,klev) /= msg ) / REAL(n)
!
!    END IF
!
!   END DO
!
!  END DO
! further references to VBAR will be replaced with VTMP

!>----------------------------------------------------------------------------<

  ALLOCATE( psitmp(0:2*nlev) )

  ALLOCATE( vvprof(0:2*nlev) )


  DO ilon = 1, nlon
    DO jlat = 1, nlat
    ! loop over latitude

            psitmp(0) = 0.

       psitmp(1:2*nlev) = msg

       vvprof(0:2*nlev:2) = 0.

      !vvprof(1:2*nlev-1:2) = VBAR(jlat,1:nlev:1)
      vvprof(1:2*nlev-1:2) = VTMP(ilon,jlat,1:nlev:1)

                      c = 2. * pi * a * COS( lat(jlat) * pi / 180. ) / g

       DO klev = 1, 2*nlev - 1, 2
       ! integrate from TOA down for each level where VBAR (vvprof) is not msg

          flag = klev

          IF ( vvprof(klev) == msg ) EXIT

          psitmp(klev+1) = psitmp(klev-1) - c * vvprof(klev) * dp(klev)

       END DO

                    psitmp(flag+1) = -psitmp(flag-1)
                                     ! impose lower boundary condition to
                                     ! ensure that MPSI is 0 at bottom
                                     ! boundary

                  psitmp(1:flag:2) = ( psitmp(2:flag+1:2) + &
                                                        psitmp(0:flag-1:2) ) / 2.
                                     ! stream function is obtained as average
                                     ! from its values at the intermediate half
                                     ! levels

          MPSI(ilon,jlat,1:nlev:1) = psitmp(1:2*nlev-1:2)

    END DO
  END DO

!>----------------------------------------------------------------------------<

  DEALLOCATE( ptmp )
  DEALLOCATE( dp )
  DEALLOCATE( VTMP )
  DEALLOCATE( psitmp )
  DEALLOCATE( vvprof )

!>============================================================================<

  END SUBROUTINE PS_CCMP_MPSI

!>============================================================================<
