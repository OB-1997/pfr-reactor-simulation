! ==============================================================================
! model.f90 — physics (kinetics + FVM RHS) and ODEPACK callbacks
! ==============================================================================
! Layout:
!   * model_helpers_mod (module) — setup_injectors, rates, source_terms,
!     compute_inlet. Anything called from main.f90 or from inside the RHS.
!   * model(neq, t, y, dydt) — free subroutine, the DLSODE callback. Must be a
!     free symbol because ODEPACK is F77 and links the RHS by name with an
!     implicit interface.
!   * JAC(...) — free subroutine, dummy Jacobian. Same F77-linkage reason as
!     model(). Never actually called when mf=25 (DLSODE builds its own
!     banded FD Jacobian) but the symbol must exist for linking.
! ==============================================================================
module model_helpers_mod
    use iso_fortran_env, only: DP => real64
    use params_mod
    implicit none

    public :: setup_injectors, rates, source_terms, compute_inlet

contains
    ! Injector positions and flow rates. Positions are kept in metres so that
    ! changing N_cells at runtime doesn't move the injectors physically — the
    ! cell index is computed from the position and the current dx.
    !
    ! The downstream injectors keep their original MATLAB-reference flow rates;
    ! the inlet H2O injector has been bumped from 1.33/8 = 0.166 mol/s to
    ! 6.0 mol/s to close the iso-thermal steam deficit.
    !
    ! Reasoning: at T = 1173 K everywhere R5/R6 run at full Arrhenius rate
    ! from the inlet, so the MATLAB-reference distribution (calibrated for
    ! a cool inlet where the reforming reactions are Arrhenius-frozen)
    ! cannot meet the local H2O demand near the inlet. Pre-loading H2O
    ! at cell 1 supplies the inlet zone via convection; the downstream
    ! injectors top up the steam budget further along the reactor. This
    ! preserves the original injection layout while closing the mass
    ! balance for steam.
    !
    !   pos [m]  = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 19.0]
    !   nO2      = [1.64/8, 0.24, 0.20, 0.21, 0.17, 0.18, 0.16, 0.15, 0.16]  mol/s
    !   nH2O     = [6.0,    0.59, 0.56, 0.43, 0.39, 0.24, 0.06, 0.01, 1e-6]  mol/s
    subroutine setup_injectors()
        integer, parameter :: N_inj = 9
        real(DP), parameter :: inj_pos_m(N_inj) = &
                               [0.0_DP, 0.5_DP, 1.0_DP, &
                                2.0_DP, 3.0_DP, 5.0_DP, &
                                8.0_DP, 13.0_DP, 19.0_DP]
        real(DP), parameter :: rate_O2(N_inj) = [1.64_DP / 8.0_DP, 0.24_DP, 0.20_DP, &
                                                 0.21_DP, 0.17_DP, 0.18_DP, &
                                                 0.16_DP, 0.15_DP, 0.16_DP]
        real(DP), parameter :: rate_H2O(N_inj) = &
                               [6.0_DP, 0.59_DP, 0.56_DP, &
                                0.43_DP, 0.39_DP, 0.24_DP, &
                                0.06_DP, 0.01_DP, 1.0e-6_DP]
        integer :: k, j

        do k = 1, N_inj
            ! Map axial position to the cell whose centre is closest. Cell j has
            ! centre x_j = (j - 0.5)*dx, so j = floor(x/dx) + 1; clamp to [1, N_cells].
            j = max(1, min(N_cells, int(inj_pos_m(k) / dx) + 1))
            inj_active(iO2, j) = .true.; inj_rate(iO2, j) = rate_O2(k)
            inj_active(iH2O, j) = .true.; inj_rate(iH2O, j) = rate_H2O(k)
        end do

        ! Scale every active injector by the global F_scale knob (§4.2.4).
        ! Default F_scale = 1.0 leaves the calibrated rates unchanged.
        inj_rate = inj_rate * F_scale
    end subroutine setup_injectors

    subroutine rates(C, r)
        real(DP), intent(in)  :: C(N_species)
        real(DP), intent(out) :: r(8)

        real(DP) :: k1c, k2c, k3c, k4c, k5c, k6c, k7c, k8c
        real(DP) :: rho_soot

        k1c = A1 * exp(-E1 / (R_gas * T_react))
        k2c = A2 * exp(-E2 / (R_gas * T_react))
        k3c = A3 * exp(-E3 / (R_gas * T_react))
        k4c = A4 * exp(-E4 / (R_gas * T_react))
        k5c = A5 * exp(-E5 / (R_gas * T_react))
        k6c = A6 * exp(-E6 / (R_gas * T_react))
        k7c = A7 * exp(-E7 / (R_gas * T_react))
        k8c = A8 * exp(-E8 / (R_gas * T_react))

        r(1) = k1c * max(C(iCH4), 0.0_DP) * max(C(iO2), 0.0_DP)
        r(2) = k2c * max(C(iC6H6), 0.0_DP) * max(C(iO2), 0.0_DP)
        r(3) = k3c * max(C(iC9H10), 0.0_DP) * max(C(iO2), 0.0_DP)
        r(4) = k4c * max(C(iCH4), 0.0_DP) * max(C(iH2O), 0.0_DP)

        ! R5 has fractional reaction orders, including a [H2]^(-0.4) term.
        ! That blows up when H2 is near zero at the inlet, so I only evaluate
        ! the rate when all three reactants are safely positive; otherwise I
        ! just set it to zero.

        if (C(iC6H6) > 0.0_DP .and. C(iH2) > 1.0e-10_DP .and. C(iH2O) > 0.0_DP) then
            r(5) = k5c * C(iC6H6)**1.3_DP * C(iH2)**(-0.4_DP) * C(iH2O)**0.2_DP
        else
            r(5) = 0.0_DP
        end if

        ! R6: pseudo-first-order in C9H10 per Taralas (2003), matching the
        ! thesis kinetic model (eq:rate6 in Ch 1). The H2O dependence is
        ! folded into A6 under a large-excess assumption — that breaks
        ! down at the inlet of the iso-thermal model, but we test here
        ! whether the simulation can still complete without compensating
        ! tricks like f_H2O scaling.

        r(6) = k6c * max(C(iC9H10), 0.0_DP)

        if (C(iH2) > 0.0_DP) then
            r(7) = k7c * max(C(iC9H10), 0.0_DP) * max(C(iH2), 0.0_DP)**0.5_DP
        else
            r(7) = 0.0_DP
        end if

        ! R8 is written in terms of soot mass concentration, so I convert
        ! the molar concentration of C here before plugging it into the rate.

        rho_soot = max(C(iC), 0.0_DP) * M_carbon
        r(8) = k8c * rho_soot * max(C(iH2O), 0.0_DP)
    end subroutine rates

    subroutine source_terms(C, Rs)
        real(DP), intent(in)  :: C(N_species)
        real(DP), intent(out) :: Rs(N_species)

        real(DP) :: r(8)
        integer  :: i, j

        call rates(C, r)
        do i = 1, N_species
            Rs(i) = 0.0_DP
            do j = 1, 8
                Rs(i) = Rs(i) + nu(i, j) * r(j)
            end do
        end do
    end subroutine source_terms

    ! Inlet composition. The mole fractions below come from MATLAB ref.
    ! read_pars.m (normalised from his molar flow rates, which sum to
    ! 2.279 mol/s). Note that O2 and H2O don't appear in the main feed —
    ! they enter the reactor only through the side injectors — so their
    ! inlet mole fractions are zero here.

    subroutine compute_inlet(C_inlet)
        real(DP), intent(out) :: C_inlet(N_species)

        real(DP) :: C_total, y_mol(N_species)

        C_total = P_react / (R_gas * T_react)

        y_mol(iCH4) = 0.51382_DP
        y_mol(iC6H6) = 0.0_DP
        y_mol(iC9H10) = 0.39711_DP
        y_mol(iO2) = 0.0_DP
        y_mol(iCO) = 0.01141_DP
        y_mol(iH2) = 0.06714_DP
        y_mol(iH2O) = 0.0_DP
        y_mol(iCO2) = 0.01053_DP
        y_mol(iC) = 0.0_DP
        y_mol(iN2) = 0.0_DP

        C_inlet(:) = y_mol(:) * C_total
    end subroutine compute_inlet

end module model_helpers_mod

! Right-hand side for DLSODE. I'm packing the state vector so that
! y( (j-1)*N_species + i ) holds the concentration of species i in cell
! j — species index varies fastest, cell index slowest. This keeps the
! per-cell chunks contiguous, which is convenient for the finite-volume
! stencil below.

subroutine model(neq, t, y, dydt)
    use iso_fortran_env, only: DP => real64
    use params_mod
    use model_helpers_mod, only: source_terms
    implicit none

    integer, intent(in)  :: neq
    real(DP), intent(in)  :: t, y(neq)
    real(DP), intent(out) :: dydt(neq)

    real(DP) :: Cj(N_species), Cjm(N_species), Cjp(N_species)
    real(DP) :: Rs(N_species)

    real(DP) :: conv, diff, inj_src, cell_volume
    integer  :: i, j, idx, idx_jm, idx_jp

    ! My RHS doesn't actually depend on t (the system is autonomous),
    ! but DLSODE still requires t in the signature. The empty associate
    ! block lexically references t so -Wall doesn't warn about an unused
    ! argument, without leaving an unreachable assignment in the body.

    associate (dummy_t => t)
    end associate

    cell_volume = A_cross * dx

    do j = 1, N_cells
        idx = (j - 1) * N_species
        do i = 1, N_species
            Cj(i) = y(idx + i)
        end do

        if (j == 1) then
            ! Dirichlet ghost cell at the inlet — the "cell to the left" of
            ! the first cell is the cached feed composition (set by main once
            ! at startup; see params_mod%C_inlet_const).
            Cjm(:) = C_inlet_const(:)
        else
            idx_jm = (j - 2) * N_species
            do i = 1, N_species
                Cjm(i) = y(idx_jm + i)
            end do
        end if

        if (j == N_cells) then
            ! Neumann (zero-gradient) ghost cell at the outlet — I mirror
            ! the last cell's composition so dC/dx = 0 at the boundary.
            Cjp(:) = Cj(:)
        else
            idx_jp = j * N_species
            do i = 1, N_species
                Cjp(i) = y(idx_jp + i)
            end do
        end if

        call source_terms(Cj, Rs)

        do i = 1, N_species
            conv = -u_vel * (Cj(i) - Cjm(i)) / dx
            diff = D_ax * (Cjp(i) - 2.0_DP * Cj(i) + Cjm(i)) / (dx * dx)
            if (inj_active(i, j)) then
                inj_src = inj_rate(i, j) / cell_volume
            else
                inj_src = 0.0_DP
            end if
            dydt(idx + i) = conv + diff + Rs(i) + inj_src
        end do
    end do
end subroutine model

! Dummy Jacobian routine. I'm using mf = 25 (DLSODE computes its own
! banded finite-difference Jacobian with ML = MU = N_species), so this
! routine is never actually called. But ODEPACK is F77 and doesn't
! support optional procedure arguments — the linker still wants a JAC
! symbol, so I provide an empty one. The assignment inside the
! `if (.false.)` block is just to silence the unused-argument warnings.
! Option C (associate) doesn't fit cleanly here because PD is intent-out
! and seven F77-style dummy args would need separate handling; the
! one-liner is simpler for an F77-shaped routine.

SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD)
    INTEGER          :: NEQ, ML, MU, NROWPD
    DOUBLE PRECISION :: T, Y(*), PD(NROWPD, *)
    if (.false.) PD(1, 1) = T + Y(1) + real(ML + MU + NROWPD + NEQ, kind=kind(T))
END SUBROUTINE
