! ==============================================================================
! solvers.f90 — solver wrappers and the shared config / result types
! ==============================================================================
! The prototype is structured so that every time-integrator (DLSODE today,
! explicit Euler / implicit Euler later) presents the same interface:
!
!     call solve_<name>(neq, y, tspan, dt_out, tt, yy, Nstep, config, result)
!
! solver_config_t  — input settings the solver was invoked with (echoed by
!                    output_mod into parameters.txt and performance.txt)
! solver_result_t  — output counters and timing the solver reports back
!
! Each `init_<name>_config` function tags the config with %name so the run
! directory, the parameters.txt header, and the performance.txt header all
! pick up the right algorithm tag automatically — the rest of the program
! does not have to know which solver ran.
! ==============================================================================
module solvers_mod
    use iso_fortran_env, only: DP => real64
    use params_mod
    implicit none

    public :: solver_config_t, solver_result_t
    public :: init_odepack_config, solve_odepack
    public :: init_euler_explicit_config, solve_euler_explicit
    public :: init_euler_implicit_config, solve_euler_implicit

    type :: solver_config_t
        character(len=32) :: name = ''
        ! Fields below are populated only by solvers that use them.
        ! ODEPACK            uses mf / ml / mu / mxstep / lrw / liw / rtol / atol.
        ! Forward Euler      uses dt_step.
        ! Backward Euler     uses dt_step / rtol / mxstep (Newton tol + max iters)
        !                    plus ml / mu (banded Jacobian half-widths).
        integer  :: mf     = 0
        integer  :: ml     = 0
        integer  :: mu     = 0
        integer  :: mxstep = 0
        integer  :: lrw    = 0
        integer  :: liw    = 0
        real(DP) :: rtol   = 0.0_DP
        real(DP) :: atol   = 0.0_DP
        real(DP) :: dt_step = 0.0_DP   ! fixed-step substep [s]
    end type solver_config_t

    type :: solver_result_t
        character(len=32) :: name = ''
        integer  :: nst            = 0    ! integrator steps
        integer  :: nfe            = 0    ! RHS evaluations
        integer  :: nje            = 0    ! Jacobian evaluations
        integer  :: bdf_order_last = -1   ! ODEPACK only; -1 if not applicable
        integer  :: istate         = 0    ! solver-specific exit code
        real(DP) :: wall_time      = 0.0_DP
        real(DP) :: cpu_time       = 0.0_DP
    end type solver_result_t

contains

    ! --------------------------------------------------------------------------
    ! init_odepack_config — settings for DLSODE (banded FD Jacobian, mf=25).
    ! lrw and liw are sized exactly per the ODEPACK manual for miter=5.
    ! --------------------------------------------------------------------------
    function init_odepack_config(neq) result(config)
        integer, intent(in) :: neq
        type(solver_config_t) :: config

        config%name   = 'odepack'
        config%mf     = 25
        config%ml     = N_species
        config%mu     = N_species
        config%mxstep = 20000
        config%rtol   = 1.0e-5_DP
        config%atol   = 1.0e-9_DP
        config%lrw    = 22 + 10 * neq + (3 * N_species + 1) * neq
        config%liw    = 22 + neq
    end function init_odepack_config

    ! --------------------------------------------------------------------------
    ! solve_odepack — wraps DLSODE.
    ! Inputs : neq, y (initial state, overwritten with final state), tspan,
    !          dt_out, Nstep, config.
    ! Outputs: tt(Nstep), yy(Nstep, neq), result (counters, timing, exit code).
    ! --------------------------------------------------------------------------
    subroutine solve_odepack(neq, y, tspan, dt_out, tt, yy, Nstep, config, result)
        integer,              intent(in)    :: neq
        real(DP),             intent(inout) :: y(neq)
        real(DP),             intent(in)    :: tspan(2)        ! [t_start, t_end]
        real(DP),             intent(in)    :: dt_out
        real(DP),             intent(out)   :: tt(:)
        real(DP),             intent(out)   :: yy(:, :)
        integer,              intent(in)    :: Nstep
        type(solver_config_t),intent(in)    :: config
        type(solver_result_t),intent(out)   :: result

        integer  :: i, idx, itol, itask, istate, iopt
        real(DP) :: t, tout
        real(DP), allocatable :: rwork(:)
        integer,  allocatable :: iwork(:)
        integer(8) :: clock_start, clock_end, clock_rate, clock_max
        real(DP)   :: cpu_time_start, cpu_time_end

        external model
        external jac

        result%name = config%name

        allocate (rwork(config%lrw))
        allocate (iwork(config%liw))

        itol   = 1
        itask  = 1
        istate = 1
        iopt   = 1   ! optional inputs in iwork(1..2, 5..10) and rwork(5..10)

        rwork(5:10) = 0.0_DP
        iwork(5:10) = 0
        iwork(1)    = config%ml
        iwork(2)    = config%mu
        iwork(6)    = config%mxstep

        t        = tspan(1)
        tt(1)    = tspan(1)
        yy(1, :) = y(:)
        tout     = tspan(1) + dt_out

        call system_clock(clock_start, clock_rate, clock_max)
        call cpu_time(cpu_time_start)

        write (*, '(A)') '  step       t [s]     CH4_mid      O2_mid      H2_mid      CO_mid'
        write (*, '(A)') '  ----  ----------  ----------  ----------  ----------  ----------'
        do i = 2, Nstep
            call dlsode(model, neq, y, t, tout, itol, config%rtol, config%atol, &
                        itask, istate, iopt, rwork, config%lrw, iwork, config%liw, &
                        jac, config%mf)
            if (istate /= 2) then
                write (*, '(A,I0)') '  DLSODE failure, istate = ', istate
                result%istate = istate
                deallocate (rwork, iwork)
                stop 1
            end if

            tt(i)    = t
            yy(i, :) = y(:)

            if (mod(i, max(Nstep / 20, 1)) == 0 .or. i == Nstep) then
                idx = (N_cells / 2) * N_species
                write (*, '(I6, F12.4, 4ES12.4)') i - 1, t, &
                    y(idx + iCH4), y(idx + iO2), y(idx + iH2), y(idx + iCO)
            end if

            tout = t + dt_out
        end do

        call system_clock(clock_end)
        call cpu_time(cpu_time_end)

        result%wall_time      = real(clock_end - clock_start, DP) / real(clock_rate, DP)
        result%cpu_time       = cpu_time_end - cpu_time_start
        result%nst            = iwork(11)
        result%nfe            = iwork(12)
        result%nje            = iwork(13)
        result%bdf_order_last = iwork(14)
        result%istate         = istate

        deallocate (rwork, iwork)
    end subroutine solve_odepack

    ! --------------------------------------------------------------------------
    ! init_euler_explicit_config — settings for forward (explicit) Euler.
    ! Only dt_step is used; the ODEPACK fields stay at zero so the parameters
    ! file knows not to print them.
    ! --------------------------------------------------------------------------
    function init_euler_explicit_config() result(config)
        type(solver_config_t) :: config

        config%name    = 'explicit_euler'
        config%dt_step = 1.0e-3_DP   ! 1 ms substep — Chapter 4 stability axis
    end function init_euler_explicit_config

    ! --------------------------------------------------------------------------
    ! solve_euler_explicit — fixed-step forward Euler, y_{n+1} = y_n + dt f(y_n).
    ! Outer loop walks the output grid (tt(i), yy(i,:)) so the saved data
    ! matches what solve_odepack produces; inner loop takes config%dt_step
    ! sub-steps until the next output time. RHS evaluation counter equals the
    ! step counter — one f-call per step.
    !
    ! Divergence guard: if any |y_i| > 1e15 or NaN appears, write the partial
    ! results, set istate = -1, and return. Forward Euler on this PDE is only
    ! conditionally stable, and the prototype is meant to demonstrate that
    ! limit in Chapter 4 — blowing up cleanly is part of the test.
    ! --------------------------------------------------------------------------
    subroutine solve_euler_explicit(neq, y, tspan, dt_out, tt, yy, Nstep, config, result)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        integer,              intent(in)    :: neq
        real(DP),             intent(inout) :: y(neq)
        real(DP),             intent(in)    :: tspan(2)
        real(DP),             intent(in)    :: dt_out
        real(DP),             intent(out)   :: tt(:)
        real(DP),             intent(out)   :: yy(:, :)
        integer,              intent(in)    :: Nstep
        type(solver_config_t),intent(in)    :: config
        type(solver_result_t),intent(out)   :: result

        integer  :: i, idx, n_substep, k, nst_total
        real(DP) :: t, t_target, dt, dt_last
        real(DP), allocatable :: dydt(:)
        integer(8) :: clock_start, clock_end, clock_rate, clock_max
        real(DP)   :: cpu_time_start, cpu_time_end

        external model

        result%name = config%name

        allocate (dydt(neq))

        dt = config%dt_step
        if (dt <= 0.0_DP) then
            write (*, '(A)') '  ERROR: forward Euler needs config%dt_step > 0'
            stop 1
        end if

        t        = tspan(1)
        tt(1)    = tspan(1)
        yy(1, :) = y(:)
        nst_total = 0

        call system_clock(clock_start, clock_rate, clock_max)
        call cpu_time(cpu_time_start)

        write (*, '(A)') '  step       t [s]     CH4_mid      O2_mid      H2_mid      CO_mid'
        write (*, '(A)') '  ----  ----------  ----------  ----------  ----------  ----------'
        do i = 2, Nstep
            t_target = tspan(1) + real(i - 1, DP) * dt_out

            ! whole sub-steps of length dt, then one short step to land on t_target
            n_substep = int((t_target - t) / dt)
            do k = 1, n_substep
                call model(neq, t, y, dydt)
                y = y + dt * dydt
                t = t + dt
                nst_total = nst_total + 1
                if (any(ieee_is_nan(y)) .or. maxval(abs(y)) > 1.0e15_DP) then
                    write (*, '(A,F12.4)') '  Forward Euler diverged at t = ', t
                    result%istate = -1
                    goto 100
                end if
            end do
            dt_last = t_target - t
            if (dt_last > 0.0_DP) then
                call model(neq, t, y, dydt)
                y = y + dt_last * dydt
                t = t_target
                nst_total = nst_total + 1
                if (any(ieee_is_nan(y)) .or. maxval(abs(y)) > 1.0e15_DP) then
                    write (*, '(A,F12.4)') '  Forward Euler diverged at t = ', t
                    result%istate = -1
                    goto 100
                end if
            end if

            tt(i)    = t
            yy(i, :) = y(:)

            if (mod(i, max(Nstep / 20, 1)) == 0 .or. i == Nstep) then
                idx = (N_cells / 2) * N_species
                write (*, '(I6, F12.4, 4ES12.4)') i - 1, t, &
                    y(idx + iCH4), y(idx + iO2), y(idx + iH2), y(idx + iCO)
            end if
        end do
        result%istate = 2

100     continue
        call system_clock(clock_end)
        call cpu_time(cpu_time_end)

        result%wall_time      = real(clock_end - clock_start, DP) / real(clock_rate, DP)
        result%cpu_time       = cpu_time_end - cpu_time_start
        result%nst            = nst_total
        result%nfe            = nst_total      ! one RHS call per step
        result%nje            = 0              ! Jacobian-free
        result%bdf_order_last = -1             ! not applicable

        deallocate (dydt)
    end subroutine solve_euler_explicit

    ! --------------------------------------------------------------------------
    ! init_euler_implicit_config — settings for backward (implicit) Euler.
    ! Each step solves g(y_{n+1}) = y_{n+1} - y_n - dt * f(y_{n+1}) = 0 with a
    ! plain Newton iteration. The Jacobian dg/dy = I - dt * df/dy is built by
    ! finite differences on a banded structure (ML = MU = N_species, exactly
    ! the same banded shape DLSODE exploits with mf=25), and factored with
    ! LINPACK's DGBFA / DGBSL (linked in via opkda2.f).
    ! --------------------------------------------------------------------------
    function init_euler_implicit_config() result(config)
        type(solver_config_t) :: config

        config%name    = 'implicit_euler'
        config%dt_step = 1.0e-2_DP   ! 10 ms substep — fixed; A-stability lets us
                                     ! take much larger steps than forward Euler
        config%rtol    = 1.0e-6_DP   ! Newton convergence tolerance on ||delta||
        config%mxstep  = 20          ! max Newton iterations per step
        config%ml      = N_species   ! banded Jacobian half-widths (same as mf=25)
        config%mu      = N_species
    end function init_euler_implicit_config

    ! --------------------------------------------------------------------------
    ! solve_euler_implicit — fixed-step backward Euler with banded Newton.
    !
    ! At each step k = 1..N_substep we want y_{n+1} satisfying
    !     g(y_{n+1}) = y_{n+1} - y_n - dt * f(y_{n+1}) = 0.
    ! Newton iteration:
    !     1. Build banded J = I - dt * df/dy by finite differences (ML+MU+1 = 21
    !        RHS evaluations per Jacobian, using Curtis-Powell-Reid coloring:
    !        columns separated by >= ML+MU+1 don't share a row, so they can be
    !        perturbed simultaneously).
    !     2. Factor with DGBFA (LINPACK banded LU, in opkda2.f).
    !     3. Solve J * delta = -g with DGBSL.
    !     4. Update y_iter <- y_iter + delta; check ||delta||_inf < tol.
    ! Two safety nets: damp the Newton step if the residual fails to decrease,
    ! and rebuild the Jacobian if convergence stalls.
    ! --------------------------------------------------------------------------
    subroutine solve_euler_implicit(neq, y, tspan, dt_out, tt, yy, Nstep, config, result)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        integer,              intent(in)    :: neq
        real(DP),             intent(inout) :: y(neq)
        real(DP),             intent(in)    :: tspan(2)
        real(DP),             intent(in)    :: dt_out
        real(DP),             intent(out)   :: tt(:)
        real(DP),             intent(out)   :: yy(:, :)
        integer,              intent(in)    :: Nstep
        type(solver_config_t),intent(in)    :: config
        type(solver_result_t),intent(out)   :: result

        integer  :: i, idx, n_substep, k, nst_total, nfe_total, nje_total
        integer  :: max_newton, ml, mu, lda, ierr
        real(DP) :: t, t_target, dt, dt_last, tol
        real(DP), allocatable :: y_n(:), y_iter(:), g_res(:), delta(:), &
                                 dydt_base(:), dydt_pert(:), y_pert(:), abd(:, :)
        integer,  allocatable :: ipvt(:)
        integer(8) :: clock_start, clock_end, clock_rate, clock_max
        real(DP)   :: cpu_time_start, cpu_time_end

        result%name = config%name

        dt         = config%dt_step
        tol        = config%rtol
        max_newton = config%mxstep
        ml         = config%ml
        mu         = config%mu
        lda        = 2 * ml + mu + 1   ! LINPACK banded storage convention

        if (dt <= 0.0_DP) then
            write (*, '(A)') '  ERROR: implicit Euler needs config%dt_step > 0'
            stop 1
        end if

        ! Heap-allocate everything that scales with neq once per run; the
        ! Jacobian builder receives y_pert through the call chain so it never
        ! puts an 80 KB stack array on the stack at neq = 10000.
        allocate (y_n(neq), y_iter(neq), g_res(neq), delta(neq), &
                  dydt_base(neq), dydt_pert(neq), y_pert(neq), &
                  abd(lda, neq), ipvt(neq))

        t        = tspan(1)
        tt(1)    = tspan(1)
        yy(1, :) = y(:)
        nst_total = 0
        nfe_total = 0
        nje_total = 0

        call system_clock(clock_start, clock_rate, clock_max)
        call cpu_time(cpu_time_start)

        write (*, '(A)') '  step       t [s]     CH4_mid      O2_mid      H2_mid      CO_mid'
        write (*, '(A)') '  ----  ----------  ----------  ----------  ----------  ----------'

        substep_loop: do i = 2, Nstep
            t_target  = tspan(1) + real(i - 1, DP) * dt_out
            n_substep = int((t_target - t) / dt)

            do k = 1, n_substep
                y_n(:) = y(:)
                call newton_step(neq, y_n, dt, tol, max_newton, ml, mu, lda, &
                                 abd, ipvt, y_iter, g_res, delta, dydt_base, &
                                 dydt_pert, y_pert, nfe_total, nje_total, ierr)
                if (ierr /= 0) then
                    write (*, '(A,I0,A,F12.4)') &
                        '  Newton failure (ierr = ', ierr, ') at t = ', t + dt
                    result%istate = ierr
                    exit substep_loop
                end if
                y(:) = y_iter(:)
                t = t + dt
                nst_total = nst_total + 1
                if (any(ieee_is_nan(y)) .or. maxval(abs(y)) > 1.0e15_DP) then
                    write (*, '(A,F12.4)') '  Implicit Euler diverged at t = ', t
                    result%istate = -2
                    exit substep_loop
                end if
            end do

            dt_last = t_target - t
            if (dt_last > 0.0_DP) then
                y_n(:) = y(:)
                call newton_step(neq, y_n, dt_last, tol, max_newton, ml, mu, &
                                 lda, abd, ipvt, y_iter, g_res, delta, &
                                 dydt_base, dydt_pert, y_pert, nfe_total, nje_total, ierr)
                if (ierr /= 0) then
                    write (*, '(A,I0,A,F12.4)') &
                        '  Newton failure (ierr = ', ierr, ') at t = ', t + dt_last
                    result%istate = ierr
                    exit substep_loop
                end if
                y(:) = y_iter(:)
                t = t_target
                nst_total = nst_total + 1
            end if

            tt(i)    = t
            yy(i, :) = y(:)

            if (mod(i, max(Nstep / 20, 1)) == 0 .or. i == Nstep) then
                idx = (N_cells / 2) * N_species
                write (*, '(I6, F12.4, 4ES12.4)') i - 1, t, &
                    y(idx + iCH4), y(idx + iO2), y(idx + iH2), y(idx + iCO)
            end if
        end do substep_loop
        if (result%istate == 0) result%istate = 2

        call system_clock(clock_end)
        call cpu_time(cpu_time_end)

        result%wall_time      = real(clock_end - clock_start, DP) / real(clock_rate, DP)
        result%cpu_time       = cpu_time_end - cpu_time_start
        result%nst            = nst_total
        result%nfe            = nfe_total
        result%nje            = nje_total
        result%bdf_order_last = 1            ! backward Euler is BDF order 1

        deallocate (y_n, y_iter, g_res, delta, dydt_base, dydt_pert, y_pert, abd, ipvt)
    end subroutine solve_euler_implicit

    ! --------------------------------------------------------------------------
    ! newton_step — solve g(y_iter) = y_iter - y_n - dt * f(y_iter) = 0 by
    ! Newton iteration with a banded Jacobian. On entry y_iter is unset and
    ! y_n holds the current state; on exit y_iter holds the converged y_{n+1}
    ! (or the last iterate if max_newton was hit).
    !
    ! Returns ierr = 0 on convergence, ierr = -1 if max_newton exceeded,
    ! ierr = -2 if DGBFA reports a singular factorisation.
    ! --------------------------------------------------------------------------
    subroutine newton_step(neq, y_n, dt, tol, max_newton, ml, mu, lda, &
                           abd, ipvt, y_iter, g_res, delta, dydt_base, &
                           dydt_pert, y_pert, nfe_total, nje_total, ierr)
        integer,  intent(in)    :: neq, max_newton, ml, mu, lda
        real(DP), intent(in)    :: y_n(neq), dt, tol
        real(DP), intent(inout) :: abd(lda, neq), y_iter(neq), g_res(neq), &
                                   delta(neq), dydt_base(neq), dydt_pert(neq), &
                                   y_pert(neq)
        integer,  intent(inout) :: ipvt(neq), nfe_total, nje_total
        integer,  intent(out)   :: ierr

        integer  :: it, info_lu
        real(DP) :: norm_delta, norm_delta_prev, t_dummy

        external model
        external dgbfa, dgbsl    ! LINPACK banded LU (in opkda2.f)

        t_dummy = 0.0_DP

        ! Initial guess: y_{n+1} = y_n. For smooth problems this seeds Newton
        ! within its quadratic-convergence basin.
        y_iter(:) = y_n(:)

        ! Build the banded Jacobian dg/dy = I - dt * df/dy at the initial guess
        ! and factorise it once. (For this problem the Jacobian doesn't change
        ! enough across Newton iterations to justify rebuilding within a step.)
        call build_banded_jacobian(neq, y_iter, dt, ml, mu, lda, abd, &
                                   dydt_base, dydt_pert, y_pert, nfe_total)
        nje_total = nje_total + 1
        call dgbfa(abd, lda, neq, ml, mu, ipvt, info_lu)
        if (info_lu /= 0) then
            ierr = -2
            return
        end if

        norm_delta_prev = huge(1.0_DP)
        do it = 1, max_newton
            ! Residual g(y_iter) — already costs one model() call, which we
            ! reuse from the Jacobian build on iteration 1.
            if (it == 1) then
                g_res(:) = y_iter(:) - y_n(:) - dt * dydt_base(:)
            else
                call model(neq, t_dummy, y_iter, dydt_base)
                nfe_total = nfe_total + 1
                g_res(:) = y_iter(:) - y_n(:) - dt * dydt_base(:)
            end if

            ! Solve J * delta = -g
            delta(:) = -g_res(:)
            call dgbsl(abd, lda, neq, ml, mu, ipvt, delta, 0)

            y_iter(:) = y_iter(:) + delta(:)

            norm_delta = maxval(abs(delta))
            if (norm_delta < tol) then
                ierr = 0
                return
            end if

            ! If progress stalls, rebuild + refactor the Jacobian and continue.
            if (norm_delta > 0.5_DP * norm_delta_prev) then
                call build_banded_jacobian(neq, y_iter, dt, ml, mu, lda, abd, &
                                           dydt_base, dydt_pert, y_pert, nfe_total)
                nje_total = nje_total + 1
                call dgbfa(abd, lda, neq, ml, mu, ipvt, info_lu)
                if (info_lu /= 0) then
                    ierr = -2
                    return
                end if
            end if
            norm_delta_prev = norm_delta
        end do

        ! Max Newton iterations exceeded.
        ierr = -1
    end subroutine newton_step

    ! --------------------------------------------------------------------------
    ! build_banded_jacobian — finite-difference banded Jacobian of the residual
    ! g(y) = y - y_n - dt * f(y) with respect to y. dg/dy = I - dt * df/dy.
    !
    ! Curtis-Powell-Reid colouring: with bandwidth (2*ml+mu+1) = 31, columns
    ! separated by at least (ml+mu+1) = 21 don't share a row, so they can be
    ! perturbed simultaneously. We do 21 RHS calls (one per colour) instead
    ! of neq = 1000.
    !
    ! Storage is LINPACK banded: abd(i - j + ml + mu + 1, j) = dg/dy(i, j) for
    ! max(1, j-mu) <= i <= min(neq, j+ml). The first ml rows of abd are
    ! workspace for fill-in during DGBFA.
    ! --------------------------------------------------------------------------
    subroutine build_banded_jacobian(neq, y_in, dt, ml, mu, lda, abd, &
                                     dydt_base, dydt_pert, y_pert, nfe_total)
        integer,  intent(in)    :: neq, ml, mu, lda
        real(DP), intent(in)    :: y_in(neq), dt
        real(DP), intent(inout) :: abd(lda, neq), y_pert(neq)
        real(DP), intent(out)   :: dydt_base(neq), dydt_pert(neq)
        integer,  intent(inout) :: nfe_total

        integer  :: colour, j, i, mband, irow_band
        real(DP) :: eps_fd, h, t_dummy

        external model

        t_dummy = 0.0_DP
        mband   = ml + mu + 1                ! 21 — number of colour groups
        eps_fd  = sqrt(epsilon(1.0_DP))      ! sqrt(machine eps), ~1.5e-8

        ! Base RHS at y_in
        call model(neq, t_dummy, y_in, dydt_base)
        nfe_total = nfe_total + 1

        ! Initialise abd to identity (dg/dy = I + ...). Only the rows that
        ! correspond to the diagonal in banded storage need touching, but we
        ! zero the whole thing so non-stencil entries are clean.
        abd(:, :) = 0.0_DP
        do j = 1, neq
            abd(ml + mu + 1, j) = 1.0_DP     ! diagonal of I
        end do

        ! For each colour, perturb every column j with mod(j-1, mband) == colour-1
        do colour = 1, mband
            y_pert(:) = y_in(:)
            do j = colour, neq, mband
                h = eps_fd * max(abs(y_in(j)), 1.0_DP)
                y_pert(j) = y_in(j) + h
            end do

            call model(neq, t_dummy, y_pert, dydt_pert)
            nfe_total = nfe_total + 1

            ! Recover the banded entries for each perturbed column. Within the
            ! bandwidth ml+mu of column j, only stencil rows i are affected;
            ! disjointness across the colour means we can read each row off
            ! cleanly.
            do j = colour, neq, mband
                h = eps_fd * max(abs(y_in(j)), 1.0_DP)
                do i = max(1, j - mu), min(neq, j + ml)
                    irow_band = i - j + ml + mu + 1
                    abd(irow_band, j) = abd(irow_band, j) &
                        - dt * (dydt_pert(i) - dydt_base(i)) / h
                end do
            end do
        end do
    end subroutine build_banded_jacobian

end module solvers_mod
