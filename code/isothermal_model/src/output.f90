! ==============================================================================
! output.f90 — All file I/O and run-directory bookkeeping
! ==============================================================================
! Public API:
!   setup_run_dir(name, run_dir)         - find next free output/<name>_<N>
!                                          and mkdir it
!   write_parameters(run_dir, config,    - one-shot record of geometry, kinetics,
!                    tbeg, tend, tstep,    inlet, injectors, time-window, solver
!                    Nstep, C_inlet)
!   write_output(run_dir, tt, yy)        - the simulation table (ES15.6E3)
!   write_performance(run_dir, config,   - performance.txt + stdout summary
!                     result, Nstep, neq,
!                     tbeg, tend, tstep)
!
! Algorithm name is taken from config%name throughout — output_mod has no
! solver-specific knowledge. To support a new solver, just extend solvers.f90
! and main; nothing here needs to change.
! ==============================================================================
module output_mod
    use iso_fortran_env, only: DP => real64
    use params_mod
    use solvers_mod, only: solver_config_t, solver_result_t
    use model_helpers_mod, only: rates
    implicit none

    public :: setup_run_dir, write_parameters, write_output, write_performance, &
              write_integrated_rates

contains

    ! --------------------------------------------------------------------------
    ! setup_run_dir — find next free output/<name>_<N> directory and create it.
    ! --------------------------------------------------------------------------
    subroutine setup_run_dir(name, run_dir)
        character(len=*), intent(in)  :: name
        character(len=*), intent(out) :: run_dir

        character(len=128) :: temp
        integer :: test_num, mkdir_status
        logical :: dir_exists

        test_num = 1
        do
            write (temp, '(A,A,A,I0)') 'output/', trim(name), '_', test_num
            inquire (file=trim(temp)//'/.', exist=dir_exists)
            if (.not. dir_exists) exit
            test_num = test_num + 1
            if (test_num > 9999) then
                write (*, '(A)') '  ERROR: too many existing run dirs in output/, clean up first'
                stop 1
            end if
        end do
        run_dir = trim(temp)
        call execute_command_line('mkdir "'//trim(run_dir)//'"', &
                                  wait=.true., exitstat=mkdir_status)
        if (mkdir_status /= 0) then
            write (*, '(A)') '  WARNING: mkdir returned non-zero; run files may not be saved'
        end if
        write (*, '(A,A)') '  Run directory                 = ', trim(run_dir)
        write (*, '(A)') ''
    end subroutine setup_run_dir

    ! --------------------------------------------------------------------------
    ! write_parameters — full record of model + solver settings for this run.
    ! Written before integration so we have the file even if the solver fails.
    ! --------------------------------------------------------------------------

    subroutine write_parameters(run_dir, config, tbeg, tend, tstep, Nstep, C_inlet)
        character(len=*), intent(in) :: run_dir
        type(solver_config_t), intent(in) :: config
        real(DP), intent(in) :: tbeg, tend, tstep
        integer, intent(in) :: Nstep
        real(DP), intent(in) :: C_inlet(N_species)

        integer  :: param_unit, i, j
        real(DP) :: C_total

        C_total = P_react / (R_gas * T_react)
        open (newunit=param_unit, file=trim(run_dir)//'/parameters.txt', &
              status='replace', action='write')
        write (param_unit, '(A)') '======================================================'
        write (param_unit, '(A,A)') '  ', trim(run_dir)
        write (param_unit, '(A,A,A)') '  TUBULAR REACTOR PROTOTYPE  (', trim(config%name), ')'
        write (param_unit, '(A)') '  Run parameters'
        write (param_unit, '(A)') '======================================================'
        write (param_unit, '(A)') ''
        write (param_unit, '(A)') '  Reactor geometry and operating conditions'
        write (param_unit, '(A,F12.4,A)') '    L_reactor                   = ', L_reactor, ' m'
        write (param_unit, '(A,F12.4,A)') '    d_reactor                   = ', d_reactor, ' m'
        write (param_unit, '(A,F12.6,A)') '    A_cross                     = ', A_cross, ' m^2'
        write (param_unit, '(A,F12.4,A)') '    u_vel                       = ', u_vel, ' m/s'
        write (param_unit, '(A,F12.2,A)') '    T_react                     = ', T_react, ' K'
        write (param_unit, '(A,F12.2,A)') '    P_react                     = ', P_react, ' Pa'
        write (param_unit, '(A,ES12.3,A)') '    D_ax                        = ', D_ax, ' m^2/s'
        write (param_unit, '(A,F12.4)')    '    F_scale                     = ', F_scale
        write (param_unit, '(A)') ''
        write (param_unit, '(A)') '  Grid'
        write (param_unit, '(A,I12)') '    N_cells                     = ', N_cells
        write (param_unit, '(A,I12)') '    N_species                   = ', N_species
        write (param_unit, '(A,I12)') '    NEQ                         = ', N_eqs
        write (param_unit, '(A,F12.4,A)') '    dx                          = ', dx, ' m'
        write (param_unit, '(A)') ''
        write (param_unit, '(A)') '  Kinetics  (Arrhenius pre-exp A and activation energy E in J/mol)'
        write (param_unit, '(A,ES12.4,A,ES12.4)') '    R1: A1 = ', A1, '   E1 = ', E1
        write (param_unit, '(A,ES12.4,A,ES12.4)') '    R2: A2 = ', A2, '   E2 = ', E2
        write (param_unit, '(A,ES12.4,A,ES12.4)') '    R3: A3 = ', A3, '   E3 = ', E3
        write (param_unit, '(A,ES12.4,A,ES12.4)') '    R4: A4 = ', A4, '   E4 = ', E4
        write (param_unit, '(A,ES12.4,A,ES12.4)') '    R5: A5 = ', A5, '   E5 = ', E5
        write (param_unit, '(A,ES12.4,A,ES12.4)') '    R6: A6 = ', A6, '   E6 = ', E6
        write (param_unit, '(A,ES12.4,A,ES12.4)') '    R7: A7 = ', A7, '   E7 = ', E7
        write (param_unit, '(A,ES12.4,A,ES12.4)') '    R8: A8 = ', A8, '   E8 = ', E8
        write (param_unit, '(A)') ''
        write (param_unit, '(A,ES12.4,A)') '  Total inlet concentration C_tot = ', C_total, ' mol/m^3'
        write (param_unit, '(A)') '  Inlet mole fractions (species not listed are zero)'
        if (C_inlet(iCH4) > 0.0_DP) write (param_unit, '(A,F12.5)') &
            '    y_CH4                       = ', C_inlet(iCH4) / C_total
        if (C_inlet(iC9H10) > 0.0_DP) write (param_unit, '(A,F12.5)') &
            '    y_C9H10                     = ', C_inlet(iC9H10) / C_total
        if (C_inlet(iH2) > 0.0_DP) write (param_unit, '(A,F12.5)') &
            '    y_H2                        = ', C_inlet(iH2) / C_total
        if (C_inlet(iCO) > 0.0_DP) write (param_unit, '(A,F12.5)') &
            '    y_CO                        = ', C_inlet(iCO) / C_total
        if (C_inlet(iCO2) > 0.0_DP) write (param_unit, '(A,F12.5)') &
            '    y_CO2                       = ', C_inlet(iCO2) / C_total
        write (param_unit, '(A)') ''
        write (param_unit, '(A)') '  Initial condition: vacuum  (all cells C = 0)'
        write (param_unit, '(A)') ''
        write (param_unit, '(A)') '  Active injectors   (species index, cell, rate [mol/s])'
        do j = 1, N_cells
            do i = 1, N_species
                if (inj_active(i, j)) then
                    write (param_unit, '(A,I3,A,I4,A,ES12.4)') &
                        '    species ', i, '   cell ', j, '   rate = ', inj_rate(i, j)
                end if
            end do
        end do
        write (param_unit, '(A)') ''
        write (param_unit, '(A)') '  Time integration'
        write (param_unit, '(A,F12.3,A)') '    t_start                     = ', tbeg, ' s'
        write (param_unit, '(A,F12.3,A)') '    t_end                       = ', tend, ' s'
        write (param_unit, '(A,F12.4,A)') '    dt_out                      = ', tstep, ' s'
        write (param_unit, '(A,I12)') '    Output time steps           = ', Nstep - 1
        write (param_unit, '(A)') ''
        write (param_unit, '(A,A,A)') '  Solver  (', trim(config%name), ')'
        if (trim(config%name) == 'odepack') then
            write (param_unit, '(A,I12)') '    method_flag (mf)            = ', config%mf
            write (param_unit, '(A,I12)') '    ML (lower bandwidth)        = ', config%ml
            write (param_unit, '(A,I12)') '    MU (upper bandwidth)        = ', config%mu
            write (param_unit, '(A,I12)') '    MXSTEP                      = ', config%mxstep
            write (param_unit, '(A,ES12.3)') '    rtol                        = ', config%rtol
            write (param_unit, '(A,ES12.3)') '    atol                        = ', config%atol
            write (param_unit, '(A,I12)') '    LRW                         = ', config%lrw
            write (param_unit, '(A,I12)') '    LIW                         = ', config%liw
        else if (trim(config%name) == 'explicit_euler') then
            write (param_unit, '(A,ES12.4,A)') '    dt_step                     = ', config%dt_step, ' s'
        else if (trim(config%name) == 'implicit_euler') then
            write (param_unit, '(A,ES12.4,A)') '    dt_step                     = ', config%dt_step, ' s'
            write (param_unit, '(A,ES12.3)') '    Newton tol (rtol)           = ', config%rtol
            write (param_unit, '(A,I12)') '    Max Newton iters (mxstep)   = ', config%mxstep
            write (param_unit, '(A,I12)') '    Banded ML (lower)           = ', config%ml
            write (param_unit, '(A,I12)') '    Banded MU (upper)           = ', config%mu
        end if
        write (param_unit, '(A)') ''
        write (param_unit, '(A)') '======================================================'
        close (param_unit)
        write (*, '(A,A)') '  Parameters written to ', trim(run_dir)//'/parameters.txt'
        write (*, '(A)') ''
    end subroutine write_parameters

    ! --------------------------------------------------------------------------
    ! write_output — the simulation data table.
    ! ES15.6E3 forces a 3-digit exponent so numpy.loadtxt works for any range.
    ! --------------------------------------------------------------------------

    subroutine write_output(run_dir, tt, yy)
        character(len=*), intent(in) :: run_dir
        real(DP), intent(in) :: tt(:), yy(:, :)
        integer :: i, n

        n = size(tt)
        open (33, file=trim(run_dir)//'/output.txt', status='replace', action='write')
        do i = 1, n
            ! Unlimited-repeat format (F2008+): writes one ES15.6E3 token per
            ! value automatically, so this scales to any NEQ without the old
            ! 10000-column ceiling silently truncating large grids.
            write (33, '(*(ES15.6E3))') tt(i), yy(i, :)
        end do
        close (33)
        write (*, '(A,A)') '  Output written to ', trim(run_dir)//'/output.txt'
    end subroutine write_output

    ! --------------------------------------------------------------------------
    ! write_performance — performance.txt + stdout summary.
    ! Same numbers in both, formatted the same way; stdout first so the user
    ! sees them while the file is still being written.
    ! --------------------------------------------------------------------------

    subroutine write_performance(run_dir, config, result, Nstep, neq, tbeg, tend, tstep)
        character(len=*), intent(in) :: run_dir
        type(solver_config_t), intent(in) :: config
        type(solver_result_t), intent(in) :: result
        integer, intent(in) :: Nstep, neq
        real(DP), intent(in) :: tbeg, tend, tstep

        integer  :: perf_unit
        real(DP) :: rhs_per_step, time_per_step

        if (result%nst > 0) then
            rhs_per_step = real(result%nfe, DP) / real(result%nst, DP)
            time_per_step = result%wall_time / real(result%nst, DP)
        else
            rhs_per_step = 0.0_DP
            time_per_step = 0.0_DP
        end if

        ! ----- stdout summary -----
        write (*, '(A)') ''
        write (*, '(A)') '  =================================================='
        write (*, '(A,A,A)') '   ', trim(config%name), '  PERFORMANCE SUMMARY'
        write (*, '(A)') '  =================================================='
        write (*, '(A,ES12.4,A)') '    Wall-clock time             = ', result%wall_time, ' s'
        write (*, '(A,ES12.4,A)') '    CPU time                    = ', result%cpu_time, ' s'
        write (*, '(A,ES12.4,A)') '    CPU/Wall ratio              = ', &
            result%cpu_time / max(result%wall_time, 1e-9_DP), '     (1.0 = single-threaded)'
        write (*, '(A)') ''
        write (*, '(A,I12)') '    Output time steps           = ', Nstep - 1
        write (*, '(A,I12)') '    Internal steps (NST)        = ', result%nst
        write (*, '(A,I12)') '    RHS evaluations (NFE)       = ', result%nfe
        write (*, '(A,I12)') '    Jacobian evaluations (NJE)  = ', result%nje
        write (*, '(A,I12)') '    Last BDF order              = ', result%bdf_order_last
        write (*, '(A)') ''
        write (*, '(A,F12.2)') '    RHS evals per step          = ', rhs_per_step
        write (*, '(A,ES12.4,A)') '    Wall time per step          = ', time_per_step, ' s'
        write (*, '(A)') '  =================================================='
        write (*, '(A)') ''

        ! ----- performance.txt -----
        open (newunit=perf_unit, file=trim(run_dir)//'/performance.txt', &
              status='replace', action='write')
        write (perf_unit, '(A)') '======================================================'
        write (perf_unit, '(A,A,A)') '  TUBULAR REACTOR PROTOTYPE  (', trim(config%name), ')'
        write (perf_unit, '(A)') '  PERFORMANCE SUMMARY'
        write (perf_unit, '(A)') '======================================================'
        write (perf_unit, '(A)') ''
        write (perf_unit, '(A)') '  Problem'
        write (perf_unit, '(A,I12)') '    N_cells                     = ', N_cells
        write (perf_unit, '(A,I12)') '    N_species                   = ', N_species
        write (perf_unit, '(A,I12)') '    NEQ                         = ', neq
        write (perf_unit, '(A,F12.3,A)') '    t_start                     = ', tbeg, ' s'
        write (perf_unit, '(A,F12.3,A)') '    t_end                       = ', tend, ' s'
        write (perf_unit, '(A,F12.4,A)') '    dt_out                      = ', tstep, ' s'
        write (perf_unit, '(A,I12)') '    Output time steps           = ', Nstep - 1
        write (perf_unit, '(A)') ''
        write (perf_unit, '(A,A,A)') '  Solver  (', trim(config%name), ')'
        if (trim(config%name) == 'odepack') then
            write (perf_unit, '(A,I12)') '    method_flag (mf)            = ', config%mf
            write (perf_unit, '(A,I12)') '    ML (lower bandwidth)        = ', config%ml
            write (perf_unit, '(A,I12)') '    MU (upper bandwidth)        = ', config%mu
            write (perf_unit, '(A,I12)') '    MXSTEP                      = ', config%mxstep
            write (perf_unit, '(A,ES12.3)') '    rtol                        = ', config%rtol
            write (perf_unit, '(A,ES12.3)') '    atol                        = ', config%atol
            write (perf_unit, '(A,I12)') '    LRW                         = ', config%lrw
            write (perf_unit, '(A,I12)') '    LIW                         = ', config%liw
        else if (trim(config%name) == 'explicit_euler') then
            write (perf_unit, '(A,ES12.4,A)') '    dt_step                     = ', config%dt_step, ' s'
        else if (trim(config%name) == 'implicit_euler') then
            write (perf_unit, '(A,ES12.4,A)') '    dt_step                     = ', config%dt_step, ' s'
            write (perf_unit, '(A,ES12.3)') '    Newton tol (rtol)           = ', config%rtol
            write (perf_unit, '(A,I12)') '    Max Newton iters (mxstep)   = ', config%mxstep
            write (perf_unit, '(A,I12)') '    Banded ML (lower)           = ', config%ml
            write (perf_unit, '(A,I12)') '    Banded MU (upper)           = ', config%mu
        end if
        write (perf_unit, '(A)') ''
        write (perf_unit, '(A)') '  Timing'
        write (perf_unit, '(A,ES12.4,A)') '    Wall-clock time             = ', result%wall_time, ' s'
        write (perf_unit, '(A,ES12.4,A)') '    CPU time                    = ', result%cpu_time, ' s'
        write (perf_unit, '(A,ES12.4,A)') '    CPU/Wall ratio              = ', &
            result%cpu_time / max(result%wall_time, 1e-9_DP), '     (1.0 = single-threaded)'
        write (perf_unit, '(A)') ''
        write (perf_unit, '(A)') '  Solver counters'
        write (perf_unit, '(A,I12)') '    Internal steps (NST)        = ', result%nst
        write (perf_unit, '(A,I12)') '    RHS evaluations (NFE)       = ', result%nfe
        write (perf_unit, '(A,I12)') '    Jacobian evaluations (NJE)  = ', result%nje
        write (perf_unit, '(A,I12)') '    Last BDF order              = ', result%bdf_order_last
        write (perf_unit, '(A)') ''
        write (perf_unit, '(A)') '  Derived'
        write (perf_unit, '(A,F12.2)') '    RHS evals per step          = ', rhs_per_step
        write (perf_unit, '(A,ES12.4,A)') '    Wall time per step          = ', time_per_step, ' s'
        write (perf_unit, '(A)') ''
        write (perf_unit, '(A)') '======================================================'
        close (perf_unit)
        write (*, '(A,A)') '  Performance summary written to ', trim(run_dir)//'/performance.txt'
        write (*, '(A)') ''
    end subroutine write_performance

    ! --------------------------------------------------------------------------
    ! write_integrated_rates — spatial-integrated steady-state reaction rates.
    !
    ! Takes the last time row of yy as the (approximate) steady state, computes
    ! r_k(C(x_j)) per cell for k = 1..8, multiplies by cell volume, sums over
    ! cells, and writes a CSV. Two files:
    !   integrated_rates.csv    — one row, columns R1..R8 [mol/s] total reactor
    !   rates_per_cell.csv      — one row per cell, columns x, R1..R8 [mol/m^3/s]
    !
    ! Hypothesis being tested (defence prep, May 2026): does the R6 fraction of
    ! total aromatic-carbon conversion grow with T in a way that explains the
    ! non-monotonic CO outlet? This routine produces the data to find out.
    ! --------------------------------------------------------------------------
    subroutine write_integrated_rates(run_dir, tt, yy)
        character(len=*), intent(in) :: run_dir
        real(DP), intent(in) :: tt(:), yy(:, :)

        integer  :: rates_unit, cell_unit, n_t, j, idx, k
        real(DP) :: Cj(N_species), r(8), cell_volume
        real(DP) :: R_total(8)
        real(DP) :: x_cell

        n_t = size(tt)
        cell_volume = A_cross * dx
        R_total(:) = 0.0_DP

        ! Per-cell file: x, R1..R8 (mol/m^3/s, local rate).
        open (newunit=cell_unit, file=trim(run_dir)//'/rates_per_cell.csv', &
              status='replace', action='write')
        write (cell_unit, '(A)') 'x_m,R1,R2,R3,R4,R5,R6,R7,R8'

        do j = 1, N_cells
            idx = (j - 1) * N_species
            Cj(:) = yy(n_t, idx + 1 : idx + N_species)
            call rates(Cj, r)
            x_cell = (real(j, DP) - 0.5_DP) * dx
            write (cell_unit, '(F12.6, *(",", ES15.6E3))') x_cell, (r(k), k = 1, 8)
            do k = 1, 8
                R_total(k) = R_total(k) + r(k) * cell_volume
            end do
        end do
        close (cell_unit)

        ! Reactor-total file: one row, total molar conversion rate [mol/s] per reaction.
        open (newunit=rates_unit, file=trim(run_dir)//'/integrated_rates.csv', &
              status='replace', action='write')
        write (rates_unit, '(A)') 'R1,R2,R3,R4,R5,R6,R7,R8'
        write (rates_unit, '(*(ES15.6E3, :, ","))') (R_total(k), k = 1, 8)
        close (rates_unit)

        write (*, '(A,A)') '  Integrated rates written to ', trim(run_dir)//'/integrated_rates.csv'
        write (*, '(A,A)') '  Per-cell rates written to   ', trim(run_dir)//'/rates_per_cell.csv'
    end subroutine write_integrated_rates

end module output_mod
