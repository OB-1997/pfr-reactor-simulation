! ==============================================================================
! main.f90 — Top-level driver
! ==============================================================================
! Responsibilities:
!   - parse the --solver=<name> CLI argument (no default; user must specify)
!   - print the run banner and echo model parameters
!   - assemble the initial condition and the side-stream injectors
!   - call into solvers_mod and output_mod
!
! All file I/O lives in output.f90; all integrator code lives in solvers.f90;
! the model RHS and kinetics live in model.f90; reactor and kinetic constants
! live in params.f90.
!
! Usage:
!     ./prototype --solver=odepack
! ==============================================================================
program prototype_reactor
    use iso_fortran_env, only: DP => real64
    use params_mod
    use model_helpers_mod, only: compute_inlet, setup_injectors
    use solvers_mod
    use output_mod
    implicit none

    integer  :: neq, Nstep, i, j
    real(DP) :: tbeg, tend_local, tstep
    real(DP) :: C_inlet(N_species)
    real(DP) :: cell_volume, peclet
    real(DP) :: tspan(2)

    real(DP), allocatable :: y(:), yy(:, :), tt(:)

    character(len=64)     :: solver_choice
    character(len=512)    :: run_dir, config_path
    type(solver_config_t) :: config
    type(solver_result_t) :: result

    ! --- Parse CLI args -----------------------------------------------------
    call parse_args(solver_choice, config_path)

    ! --- Load runtime config (or accept defaults) ---------------------------
    if (len_trim(config_path) > 0) then
        call read_config(trim(config_path))
    end if
    call finalize_params()

    ! --- Banner -------------------------------------------------------------
    write (*, '(A)') ''
    write (*, '(A)') '  =================================================='
    write (*, '(A,A,A)') '   TUBULAR REACTOR PROTOTYPE  (solver=', trim(solver_choice), ')'
    write (*, '(A)') '  =================================================='
    write (*, '(A,F12.4,A)')  '    L                           = ', L_reactor, ' m'
    write (*, '(A,F12.4,A)')  '    u                           = ', u_vel,     ' m/s'
    write (*, '(A,F12.2,A)')  '    T                           = ', T_react,   ' K'
    write (*, '(A,ES12.3,A)') '    D_ax                        = ', D_ax,      ' m^2/s'
    write (*, '(A,I12)')      '    N_cells                     = ', N_cells
    write (*, '(A,I12)')      '    N_species                   = ', N_species
    write (*, '(A,I12)')      '    NEQ                         = ', N_eqs
    write (*, '(A,F12.4,A)')  '    dx                          = ', dx, ' m'
    if (D_ax > 0.0_DP) then
        ! Use a runtime variable so gfortran does not const-fold the divide
        ! and emit a divide-by-zero warning under -Wall when D_ax = 0 elsewhere.
        block
          real(DP) :: d_safe
          d_safe = D_ax
          peclet = u_vel * L_reactor / d_safe
        end block
        write (*, '(A,F12.1)') '    Pe                          = ', peclet
    else
        write (*, '(A)')       '    Pe                          =     infinity  (pure PFR, D_ax=0)'
    end if
    write (*, '(A)') ''

    ! --- Initial condition + injectors --------------------------------------
    neq = N_eqs
    allocate (y(neq))

    ! Vacuum IC: every cell starts at C = 0. The uniform-feed alternative
    ! (\bar C_{i,j}(0) = C_{i,in}) places C9H10 in every cell at t = 0, so
    ! reaction R6 immediately consumes "phantom" H2O downstream before the
    ! cell-1 wavefront arrives. The vacuum start avoids that artefact.
    y(:) = 0.0_DP
    call compute_inlet(C_inlet)
    C_inlet_const(:) = C_inlet(:)   ! cache for model() — see params_mod
    call setup_injectors()
    cell_volume = A_cross * dx

    write (*, '(A)') '  Injectors (active cells):'
    write (*, '(A)') '      species     cell        rate [mol/s]'
    write (*, '(A)') '      -------     ----        ------------'
    do j = 1, N_cells
        do i = 1, N_species
            if (inj_active(i, j)) then
                write (*, '(6X, I7, 5X, I4, 8X, ES12.4)') i, j, inj_rate(i, j)
            end if
        end do
    end do
    write (*, '(A)') ''

    ! --- Time integration window -------------------------------------------
    tbeg       = t_start
    tend_local = t_end
    tstep      = dt_out
    Nstep      = int((tend_local - tbeg) / tstep) + 1
    tspan(1)   = tbeg
    tspan(2)   = tend_local
    allocate (tt(Nstep))
    allocate (yy(Nstep, neq))

    ! --- Pick solver and initialise its config ------------------------------
    select case (trim(solver_choice))
    case ('odepack')
        config = init_odepack_config(neq)
    case ('explicit_euler')
        config = init_euler_explicit_config()
    case ('implicit_euler')
        config = init_euler_implicit_config()
    case default
        write (*, '(A,A)') '  ERROR: unknown solver: ', trim(solver_choice)
        write (*, '(A)') '  Available solvers: odepack, explicit_euler, implicit_euler'
        stop 1
    end select

    write (*, '(A,I12)') '    Integration steps           = ', Nstep - 1
    if (trim(config%name) == 'odepack') then
        write (*, '(A,ES12.3)') '    rtol                        = ', config%rtol
        write (*, '(A,ES12.3)') '    atol                        = ', config%atol
        write (*, '(A,I12)')    '    LRW                         = ', config%lrw
        write (*, '(A,I12)')    '    LIW                         = ', config%liw
    else if (trim(config%name) == 'explicit_euler') then
        write (*, '(A,ES12.3,A)') '    dt_step                     = ', config%dt_step, ' s'
    else if (trim(config%name) == 'implicit_euler') then
        write (*, '(A,ES12.3,A)') '    dt_step                     = ', config%dt_step, ' s'
        write (*, '(A,ES12.3)')   '    Newton tol (rtol)           = ', config%rtol
        write (*, '(A,I12)')      '    Max Newton iters            = ', config%mxstep
        write (*, '(A,I12)')      '    Banded ML / MU              = ', config%ml
    end if
    write (*, '(A)') ''

    ! --- Set up run dir and lay down parameters.txt -------------------------
    call setup_run_dir(config%name, run_dir)
    call write_parameters(run_dir, config, tbeg, tend_local, tstep, Nstep, C_inlet)

    ! --- Run the solver -----------------------------------------------------
    select case (trim(solver_choice))
    case ('odepack')
        call solve_odepack(neq, y, tspan, tstep, tt, yy, Nstep, config, result)
    case ('explicit_euler')
        call solve_euler_explicit(neq, y, tspan, tstep, tt, yy, Nstep, config, result)
    case ('implicit_euler')
        call solve_euler_implicit(neq, y, tspan, tstep, tt, yy, Nstep, config, result)
    end select

    ! --- Save results -------------------------------------------------------
    call write_output(run_dir, tt, yy)
    call write_performance(run_dir, config, result, Nstep, neq, tbeg, tend_local, tstep)

    deallocate (y, yy, tt)

contains

    ! --------------------------------------------------------------------------
    ! parse_args — read --solver=<name> (required) and --config=<path>
    ! (optional). Also handles -h / --help.
    ! --------------------------------------------------------------------------
    subroutine parse_args(solver, cfg_path)
        character(len=*), intent(out) :: solver
        character(len=*), intent(out) :: cfg_path
        integer :: k, nargs
        character(len=512) :: arg

        solver   = ''
        cfg_path = ''
        nargs = command_argument_count()
        do k = 1, nargs
            call get_command_argument(k, arg)
            if (index(arg, '--solver=') == 1) then
                solver = value_after(arg, '--solver=')
            else if (index(arg, '--config=') == 1) then
                cfg_path = value_after(arg, '--config=')
            else if (trim(arg) == '-h' .or. trim(arg) == '--help') then
                call print_usage()
                stop 0
            end if
        end do

        if (len_trim(solver) == 0) then
            write (*, '(A)') ''
            write (*, '(A)') '  ERROR: --solver=<name> is required.'
            call print_usage()
            stop 1
        end if
    end subroutine parse_args

    ! --------------------------------------------------------------------------
    ! value_after — return the substring of `arg` that follows `prefix`.
    ! Used by parse_args so that adding a flag with a different prefix length
    ! (e.g. --max-newton-iters=) doesn't silently mis-slice the value.
    ! --------------------------------------------------------------------------
    function value_after(arg, prefix) result(v)
        character(len=*), intent(in) :: arg, prefix
        character(len=len(arg))      :: v
        v = trim(arg(len(prefix) + 1:))
    end function value_after

    subroutine print_usage()
        write (*, '(A)') ''
        write (*, '(A)') '  Usage: ./prototype --solver=<name> [--config=<path>]'
        write (*, '(A)') ''
        write (*, '(A)') '  Available solvers:'
        write (*, '(A)') '    odepack         DLSODE BDF, banded FD Jacobian (mf=25)'
        write (*, '(A)') '    explicit_euler  forward (explicit) Euler, fixed dt'
        write (*, '(A)') '    implicit_euler  backward (implicit) Euler with banded-Jacobian Newton'
        write (*, '(A)') ''
        write (*, '(A)') '  --config=<path>   load a NAMELIST config file with &reactor and/or'
        write (*, '(A)') '                    &time_window groups (see runs/baseline.nml).'
        write (*, '(A)') '                    Without --config the built-in defaults are used.'
        write (*, '(A)') ''
    end subroutine print_usage

end program prototype_reactor
