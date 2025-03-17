Module FWHDefinitions  !

    use SMConstants
    Implicit None
    private
    public OB_BUFFER_SIZE_DEFAULT, STR_LEN_OBSERVER, OBS_LENGTH
    public OB_BUFFER_SIZE
    public rho0, P0, c0, U0, M0, fwGamma2
    public getMeanStreamValues

    real(kind=RP)                                           :: rho0, P0, c0, fwGamma2
    real(kind=RP), dimension(NDIM)                          :: U0, M0
    integer, parameter         :: OB_BUFFER_SIZE_DEFAULT = 200
    integer                    :: OB_BUFFER_SIZE = OB_BUFFER_SIZE_DEFAULT
    integer, parameter         :: STR_LEN_OBSERVER = 128
    integer, parameter         :: OBS_LENGTH = 10

    contains

    Subroutine  getMeanStreamValues()

        use fluiddata

        ! local variables
        real(kind=RP)                                       :: theta, phi, U0Magnitud

        theta = refvalues % AOAtheta*(pi/180.0_RP)
        phi   = refvalues % AOAphi*(pi/180.0_RP)

        ! set 1 by default
        ! TODO use values of boundary conditions (inflow if exists or outflow, or set this defaults if not exists)
        U0Magnitud = 1.0_RP
        rho0 = 1.0_RP

        U0(1)  = U0Magnitud*cos(theta)*cos(phi)
        U0(2)  = U0Magnitud*sin(theta)*cos(phi)
        U0(3)  = U0Magnitud*sin(phi)

        M0 = U0 * dimensionless % Mach
        c0 = U0Magnitud / dimensionless % Mach
        fwGamma2 = 1.0_RP / (1.0_RP - dimensionless % Mach**2)

        ! default initial condition and outflow BC for energy without external pressure
        ! TODO include external pressure
        P0 = 1.0_RP / (dimensionless % gammaM2)
        ! rhoe0 = P0 / thermodynamics%gammaMinus1 + 0.5_RP*rho0*(U0Magnitud**2)

    End Subroutine getMeanStreamValues

End Module FWHDefinitions