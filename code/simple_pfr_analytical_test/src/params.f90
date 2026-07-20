! ==============================================================================
! params.f90 — simple A -> B PFR (Wehner-Wilhelm baseline for Ch 4 Test 2)
! ==============================================================================
! This is the sibling project to [1_isothermal_model]/. Same FVM machinery,
! same solver code, single first-order reaction A -> B with closed-form
! Wehner-Wilhelm steady-state solution. Errors measured here against the
! analytical reference give the relative-error vs N_cells data for Ch 4.
!
! Geometry, velocity and D_ax are inherited unchanged from the real model so
! the global Peclet (~1032) and residence time (~38.8 s) match the real-model
! regime. k_rxn is sized for Da ~ 1 to keep the analytical conversion in a
! sensitive range.
!
! solvers.f90 / output.f90 / main.f90 are byte-identical to the [1_isothermal_model]/
! versions. To make them compile and link, this module retains the legacy
! species-index symbols (iCH4, iC6H6, iC9H10, iO2, iCO, iH2, iH2O, iCO2, iC,
! iN2) and the eight Arrhenius pair symbols (A1..A8, E1..E8) as zero-valued
! placeholders. The species-index symbols are mapped to the two real species
! (A=1, B=2) so all y(idx + i*) lookups stay in range; the resulting
! stdout/parameters.txt labels are misleading by design but the numerical
! data is correct.
!
! T_react and P_react are also inherited from the real model but the simple
! A -> B kinetics have no Arrhenius / ideal-gas dependence, so they only
! influence echoed banner text and not the integration.
! ==============================================================================
module params_mod
    use iso_fortran_env, only: DP => real64
    implicit none

    real(DP), parameter :: R_gas = 8.314_DP
    real(DP), parameter :: pi    = 3.141592653589793_DP

    ! Reactor geometry / operating point — inherited from Table 3.1 of the
    ! real model so Pe_L = u L / D_ax is identical (~1032).
    real(DP), save :: L_reactor = 20.0_DP
    real(DP), save :: d_reactor = 0.60_DP
    real(DP), save :: A_cross   = 0.282743339_DP
    real(DP), save :: u_vel     = 0.516_DP
    real(DP), save :: T_react   = 1173.15_DP
    real(DP), save :: P_react   = 101325.0_DP
    real(DP), save :: D_ax      = 1.0e-2_DP

    integer,  save      :: N_cells   = 100
    integer,  parameter :: N_species = 2

    integer,  save :: N_eqs = 0
    real(DP), save :: dx    = 0.0_DP

    ! --------------------------------------------------------------------------
    ! Simple-PFR kinetics: A -> B, rate = k_rxn * C_A.
    ! Da = k_rxn * L / u = 0.026 * 20 / 0.516 ~ 1.008.
    ! --------------------------------------------------------------------------
    real(DP), parameter :: k_rxn  = 0.026_DP    ! [1/s], first-order rate constant
    real(DP), parameter :: C_A_in = 10.4_DP     ! [mol/m^3], inlet conc. of A

    ! Real species indices (used inside model.f90's RHS).
    integer, parameter :: iA = 1, iB = 2

    ! --------------------------------------------------------------------------
    ! Legacy species-index aliases. Kept so the byte-identical solvers.f90 and
    ! output.f90 compile. Mapped onto the two real species so every
    ! y(idx + i*) / C_inlet(i*) lookup stays in [1, N_species]:
    !   iCH4   -> A   (the only inlet species; "y_CH4 ~ 1.0" prints in
    !                   parameters.txt as a stand-in for "feed is 100% A")
    !   others -> B   (zero at inlet, so the matching y_XXX lines in
    !                   parameters.txt are guarded out by C_inlet > 0)
    ! --------------------------------------------------------------------------
    integer, parameter :: iCH4 = 1
    integer, parameter :: iC6H6 = 2, iC9H10 = 2, iO2 = 2, iCO = 2
    integer, parameter :: iH2 = 2, iH2O = 2, iCO2 = 2, iC = 2, iN2 = 2

    ! Legacy Arrhenius placeholders. output.f90 prints these as "R1: A1 = ...
    ! E1 = ..."; zero values are honest about there being no Arrhenius
    ! reactions in the simple model.
    real(DP), parameter :: A1 = 0.0_DP, E1 = 0.0_DP
    real(DP), parameter :: A2 = 0.0_DP, E2 = 0.0_DP
    real(DP), parameter :: A3 = 0.0_DP, E3 = 0.0_DP
    real(DP), parameter :: A4 = 0.0_DP, E4 = 0.0_DP
    real(DP), parameter :: A5 = 0.0_DP, E5 = 0.0_DP
    real(DP), parameter :: A6 = 0.0_DP, E6 = 0.0_DP
    real(DP), parameter :: A7 = 0.0_DP, E7 = 0.0_DP
    real(DP), parameter :: A8 = 0.0_DP, E8 = 0.0_DP

    real(DP), save :: t_start = 0.0_DP
    real(DP), save :: t_end   = 200.0_DP    ! ~5 tau; needed to reach steady state
                                            ! from the vacuum IC for validation.
    real(DP), save :: dt_out  = 1.0_DP

    ! Injector arrays kept allocatable + always-false. main.f90 calls
    ! setup_injectors() (now a no-op) and output.f90 walks inj_active to
    ! list active injectors (none, so the section in parameters.txt is empty).
    logical,  save, allocatable :: inj_active(:, :)
    real(DP), save, allocatable :: inj_rate  (:, :)

    ! Cached inlet composition (set once by main; read by model RHS).
    real(DP), save :: C_inlet_const(N_species) = 0.0_DP

    namelist /reactor/      L_reactor, d_reactor, u_vel, T_react, P_react, &
                            D_ax, N_cells
    namelist /time_window/  t_start, t_end, dt_out

contains

    subroutine read_config(path)
        character(len=*), intent(in) :: path
        integer :: u, ios

        open (newunit=u, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (*, '(A,A)') '  ERROR: cannot open config file: ', trim(path)
            stop 1
        end if

        read (u, nml=reactor, iostat=ios)
        if (ios /= 0 .and. ios /= -1) then
            write (*, '(A,I0)') '  WARNING: &reactor namelist read failed (iostat=', ios
            write (*, '(A)')    '           keeping defaults for that group.'
        end if

        rewind (u)
        read (u, nml=time_window, iostat=ios)
        if (ios /= 0 .and. ios /= -1) then
            write (*, '(A,I0)') '  WARNING: &time_window namelist read failed (iostat=', ios
            write (*, '(A)')    '           keeping defaults for that group.'
        end if

        close (u)
        write (*, '(A,A)') '  Loaded config from: ', trim(path)
    end subroutine read_config

    subroutine finalize_params()
        A_cross = 0.25_DP * pi * d_reactor**2

        N_eqs = N_species * N_cells
        dx    = L_reactor / real(N_cells, DP)

        if (allocated(inj_active)) deallocate (inj_active)
        if (allocated(inj_rate))   deallocate (inj_rate)
        allocate (inj_active(N_species, N_cells))
        allocate (inj_rate  (N_species, N_cells))
        inj_active = .false.
        inj_rate   = 0.0_DP
    end subroutine finalize_params

end module params_mod
