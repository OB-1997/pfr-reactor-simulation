module params_mod
    use iso_fortran_env, only: DP => real64
    implicit none

    real(DP), parameter :: R_gas = 8.314_DP
    real(DP), parameter :: pi = 3.141592653589793_DP

    ! Reactor geometry and operating point — defaults match Saša's MATLAB
    ! reference. These can be overridden at runtime by passing
    ! `--config=<path>` to the binary; see read_config() below for the
    ! NAMELIST format. Without --config, defaults are used.

    real(DP), save :: L_reactor = 20.0_DP
    real(DP), save :: d_reactor = 0.60_DP
    real(DP), save :: A_cross = 0.282743339_DP   ! recomputed from d_reactor in finalize_params
    real(DP), save :: u_vel = 0.516_DP
    real(DP), save :: T_react = 1173.15_DP
    real(DP), save :: P_react = 101325.0_DP

    ! Axial dispersion coefficient. Choosing D_ax for this reactor is a
    ! physics-and-numerics balancing act:
    !
    !   - MATLAB reference uses D_ax = 0 (cold inlet kinetics
    !     dominate the upstream zone, axial mixing is irrelevant).
    !   - The classical laminar Aris-Taylor estimate
    !       D_ax = D_AB + u^2 R^2 / (48 D_AB)
    !     gives ~2.5 m^2/s with D_AB ~ 2e-4 m^2/s, but I tested that
    !     value and it breaks the iso-thermal chemistry — the cell-1
    !     H2O boost smears across the whole reactor in less than one
    !     residence time, the local steam concentration goes negative
    !     in 99/100 cells, max(C,0) zeros R4/R5/R8 everywhere, and
    !     benzene/CH4 leak out the outlet without being reformed.
    !
    ! The Aris-Taylor formula assumes laminar flow, but the actual flow
    ! is transitional/early turbulent. With T = 1173 K, P = 1 atm and
    ! a fuel-gas mean molar mass M ~ 30 g/mol the density is
    ! rho = P M / (R T) ~ 0.31 kg/m^3; the dynamic viscosity of hot
    ! combustion gases is mu ~ 4.5e-5 Pa s, so
    !     Re = rho u d / mu ~ 0.31 * 0.516 * 0.60 / 4.5e-5 ~ 2150,
    ! which is at the laminar-turbulent transition for empty pipes.
    ! In that regime the right physical model is turbulent eddy
    ! dispersion (Levenspiel, "Chemical Reaction Engineering" Ch 13),
    ! where empirical correlations give a Bodenstein number
    ! Pe_d = u d / D_ax ~ 1-10 across the Re = 2000-10000 range. The
    ! conservative end (Pe_d ~ 30, slightly more dispersion-suppressed
    ! than the empirical band) yields D_ax ~ 0.01 m^2/s.
    !
    ! This value gives a global Peclet Pe_L = u L / D_ax ~ 1030
    ! (still strongly convection-dominated, upwind FVM is well
    ! justified), a cell Peclet Pe_cell = u dx / D_ax ~ 10 (firmly
    ! above the central-difference threshold of 2), and a smearing
    ! length sqrt(2 D_ax tau) ~ 0.88 m over the 38.8 s residence
    ! (about 4 cells, well below the 13.9 m that broke things at
    ! D_ax = 2.5). The non-isothermal extension is the proper fix:
    ! once cold-inlet kinetics suppress R4/R5/R8 in cells 1-15, much
    ! larger D_ax should be tractable.

    real(DP), save :: D_ax = 1.0e-2_DP

    ! N_cells is the main runtime knob — the wall-time-vs-N_cells benchmark
    ! sweeps it. N_species and the kinetic model below stay compile-time
    ! constants since they're tied to the chosen 8-reaction scheme.
    integer, save      :: N_cells = 100
    integer, parameter :: N_species = 10

    ! Derived from N_cells / L_reactor; set by finalize_params() at startup.
    integer, save :: N_eqs = 0
    real(DP), save :: dx = 0.0_DP

    integer, parameter :: iCH4 = 1, iC6H6 = 2, iC9H10 = 3, iO2 = 4, iCO = 5
    integer, parameter :: iH2 = 6, iH2O = 7, iCO2 = 8, iC = 9, iN2 = 10

    ! R8's rate law is written in terms of the mass concentration of soot,
    ! so I need the carbon molar mass here to do the conversion.
    real(DP), parameter :: M_carbon = 12.011e-3_DP

    ! Arrhenius pairs (A_j, E_j), with activation energies in J/mol.
    ! Sources I used:
    !   R1-R3 : Smoot & Pratt (1979)
    !   R4    : Jones & Lindstedt (1988)
    !   R5    : Jess (1996) — I reduced E from 493 to 443 kJ/mol during
    !           calibration. Without this reduction R5 is essentially
    !           frozen at 900 C and the aromatic profile doesn't match
    !           Saša's reference.
    !   R6-R7 : Taralas (2003)
    !   R8    : Jess (1996)

    real(DP), parameter :: A1 = 59.8_DP, E1 = 101.4e3_DP
    real(DP), parameter :: A2 = 2.07e4_DP, E2 = 80.2e3_DP
    real(DP), parameter :: A3 = 2.07e4_DP, E3 = 80.2e3_DP
    real(DP), parameter :: A4 = 3.0e8_DP, E4 = 125.5e3_DP
    real(DP), parameter :: A5 = 2.0e16_DP, E5 = 443.0e3_DP
    real(DP), parameter :: A6 = 2.3e15_DP, E6 = 356.0e3_DP
    real(DP), parameter :: A7 = 3.3e10_DP, E7 = 250.0e3_DP
    real(DP), parameter :: A8 = 3.0e11_DP, E8 = 310.0e3_DP

    real(DP), save :: t_start = 0.0_DP
    real(DP), save :: t_end = 50.0_DP
    real(DP), save :: dt_out = 0.5_DP

    ! Stoichiometry matrix nu(species, reaction): negative entries mean
    ! the species are consumed, positive means produced. I've written the
    ! literal values out one reaction per row because that's easier to
    ! read and check. Fortran stores arrays in column-major order, so the
    ! reshape below ends up placing each reaction into its own column
    ! automatically — no transpose needed.

    real(DP), parameter :: nu(N_species, 8) = reshape([ &
                                                      !  CH4     C6H6     C9H10     O2     CO        H2        H2O     CO2       C      N2
                                            -1.0_DP, 0.0_DP, 0.0_DP, -0.5_DP, 1.0_DP, 2.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, & ! R1
                                            0.0_DP, -1.0_DP, 0.0_DP, -3.0_DP, 6.0_DP, 3.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, & ! R2
                                            0.0_DP, 0.0_DP, -1.0_DP, -4.5_DP, 9.0_DP, 5.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, & ! R3
                                            -1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 3.0_DP, -1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, & ! R4
                                            2.5_DP, -1.0_DP, 0.0_DP, 0.0_DP, 2.0_DP, 0.0_DP, -2.0_DP, 0.0_DP, 1.5_DP, 0.0_DP, & ! R5
                                          0.0_DP, 0.0_DP, -1.0_DP, 0.0_DP, 6.5_DP, 16.5_DP, -11.5_DP, 2.5_DP, 0.0_DP, 0.0_DP, & ! R6
                                            3.0_DP, 1.0_DP, -1.0_DP, 0.0_DP, 0.0_DP, -4.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, & ! R7
                                             0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 1.0_DP, -1.0_DP, 0.0_DP, -1.0_DP, 0.0_DP & ! R8
                                                      ], [N_species, 8])

    ! Injector state: which (species, cell) pairs have an active injector
    ! and at what molar rate. These arrays are filled in at startup by
    ! setup_injectors() and then just read from in the RHS routine. They
    ! are allocatable because their second dimension is N_cells, which is
    ! a runtime knob now.
    logical, save, allocatable :: inj_active(:, :)
    real(DP), save, allocatable :: inj_rate(:, :)

    ! Cached inlet composition. The Dirichlet ghost cell at x = 0 is a
    ! constant — y_mol(:) * P/(R T) — so there is no point recomputing it
    ! inside model() on every RHS evaluation. main fills this once after
    ! finalize_params() (and after any --config= overrides have updated
    ! P_react / T_react), and model() reads it as a plain module variable.
    real(DP), save :: C_inlet_const(N_species) = 0.0_DP

    ! Implicit-Euler time step. Default 10 ms matches what was hardcoded prior
    ! to 2026-05-18; expose it as a namelist/CLI knob so the explicit-vs-
    ! implicit comparison of thesis §4.5 can sweep it without recompiling.
    real(DP), save :: dt_step_impl = 1.0e-2_DP

    ! Global side-injection mass-flow multiplier. setup_injectors() scales every
    ! per-injector rate (O2 and H2O alike) by F_scale, so the H2O/O2 split ratio
    ! is preserved. Thesis §4.2.4 sweeps this knob at base T and u to probe
    ! side-feed intensity as a process lever distinct from main-stream throughput.
    real(DP), save :: F_scale = 1.0_DP

    ! NAMELIST groups for runtime overrides. Geometry + numerics in one
    ! group, time window in another, solver knobs in a third. See
    ! runs/baseline.nml for the format.
    namelist /reactor/ L_reactor, d_reactor, u_vel, T_react, P_react, &
        D_ax, N_cells, F_scale
    namelist /time_window/ t_start, t_end, dt_out
    namelist /solver/ dt_step_impl

contains

    ! --------------------------------------------------------------------------
    ! read_config — load &reactor and &time_window NAMELISTs from a text file.
    ! Either (or both) groups may be omitted; missing values fall back to the
    ! module defaults. Call finalize_params() afterwards to refresh derived
    ! quantities and (re)allocate injector arrays.
    ! --------------------------------------------------------------------------
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
            write (*, '(A)') '           — keeping defaults for that group.'
        end if

        rewind (u)
        read (u, nml=time_window, iostat=ios)
        if (ios /= 0 .and. ios /= -1) then
            write (*, '(A,I0)') '  WARNING: &time_window namelist read failed (iostat=', ios
            write (*, '(A)') '           — keeping defaults for that group.'
        end if

        rewind (u)
        read (u, nml=solver, iostat=ios)
        if (ios /= 0 .and. ios /= -1) then
            write (*, '(A,I0)') '  WARNING: &solver namelist read failed (iostat=', ios
            write (*, '(A)') '           — keeping defaults for that group.'
        end if

        close (u)
        write (*, '(A,A)') '  Loaded config from: ', trim(path)
    end subroutine read_config

    ! --------------------------------------------------------------------------
    ! finalize_params — recompute derived quantities and (re)allocate injector
    ! arrays. Always call this after either reading a config or accepting the
    ! defaults; it is safe to call repeatedly.
    ! --------------------------------------------------------------------------
    subroutine finalize_params()
        ! Cross-section is always computed from the diameter; users only need
        ! to set d_reactor in the config, A_cross follows.
        A_cross = 0.25_DP * pi * d_reactor**2

        N_eqs = N_species * N_cells
        dx = L_reactor / real(N_cells, DP)

        if (allocated(inj_active)) deallocate (inj_active)
        if (allocated(inj_rate)) deallocate (inj_rate)
        allocate (inj_active(N_species, N_cells))
        allocate (inj_rate(N_species, N_cells))
        inj_active = .false.
        inj_rate = 0.0_DP
    end subroutine finalize_params

end module params_mod
