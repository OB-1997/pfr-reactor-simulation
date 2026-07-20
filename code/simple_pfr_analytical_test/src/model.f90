! ==============================================================================
! model.f90 — simple A -> B PFR (kinetics + FVM RHS) and ODEPACK callbacks
! ==============================================================================
! Same overall layout as [1_isothermal_model]/src/model.f90:
!   * model_helpers_mod       — setup_injectors (no-op), source_terms,
!                               compute_inlet. Keeps the same public surface
!                               main.f90 expects.
!   * model(neq, t, y, dydt)  — free subroutine, the DLSODE callback. F77
!                               linkage requires it to be a free symbol.
!   * JAC(...)                — free dummy Jacobian (mf=25 builds its own).
!
! Physics is deliberately trivial: convection (first-order upwind) +
! diffusion (central differences, identical to real model) + a single
! first-order reaction A -> B at rate k_rxn * C_A. No clamps, no
! Arrhenius, no injectors. Linear in y, so DLSODE's Newton converges in
! one step and isolates spatial discretisation error from nonlinear-
! iteration error.
! ==============================================================================
module model_helpers_mod
    use iso_fortran_env, only: DP => real64
    use params_mod
    implicit none

    public :: setup_injectors, source_terms, compute_inlet

contains

    ! Simple model has no side-stream injectors. Kept as a no-op so main.f90's
    ! `call setup_injectors()` still resolves; inj_active stays .false. and
    ! inj_rate stays 0 from finalize_params.
    subroutine setup_injectors()
    end subroutine setup_injectors


    ! Source term for the simple A -> B reaction:
    !   r   =  k_rxn * C_A
    !   R_A = -r,   R_B = +r
    subroutine source_terms(C, Rs)
        real(DP), intent(in)  :: C(N_species)
        real(DP), intent(out) :: Rs(N_species)

        real(DP) :: rate

        rate  = k_rxn * C(iA)
        Rs(iA) = -rate
        Rs(iB) = +rate
    end subroutine source_terms


    ! Inlet composition: pure A at C_A_in, B at zero.
    subroutine compute_inlet(C_inlet)
        real(DP), intent(out) :: C_inlet(N_species)

        C_inlet(iA) = C_A_in
        C_inlet(iB) = 0.0_DP
    end subroutine compute_inlet

end module model_helpers_mod


! Right-hand side for DLSODE. State packing matches [1_isothermal_model]/:
! y((j-1)*N_species + i) holds species i in cell j. Boundary treatment is
! also identical: Dirichlet ghost (cached inlet) at j=1, Neumann mirror at
! j=N_cells.
subroutine model(neq, t, y, dydt)
    use iso_fortran_env, only: DP => real64
    use params_mod
    use model_helpers_mod, only: source_terms
    implicit none

    integer,  intent(in)  :: neq
    real(DP), intent(in)  :: t, y(neq)
    real(DP), intent(out) :: dydt(neq)

    real(DP) :: Cj(N_species), Cjm(N_species), Cjp(N_species)
    real(DP) :: Rs(N_species)

    real(DP) :: conv, diff
    integer  :: i, j, idx, idx_jm, idx_jp

    ! Suppress unused-argument warning on t (system is autonomous).
    associate(dummy_t => t)
    end associate

    do j = 1, N_cells
      idx = (j - 1) * N_species
      do i = 1, N_species
        Cj(i) = y(idx + i)
      end do

      ! Inlet BC: Dirichlet ghost. C at the (virtual) cell j=0 is pinned to
      ! the feed composition. This is the BC the analytical D/N reference in
      ! validate.py is built against, so the FVM converges to that reference
      ! as dx -> 0.
      if (j == 1) then
        Cjm(:) = C_inlet_const(:)
      else
        idx_jm = (j - 2) * N_species
        do i = 1, N_species
          Cjm(i) = y(idx_jm + i)
        end do
      end if

      ! Outlet BC: Neumann zero-gradient. The ghost cell j=N+1 mirrors cell
      ! N, so (Cjp - Cj)/dx = 0 in the diffusion stencil below.
      if (j == N_cells) then
        Cjp(:) = Cj(:)
      else
        idx_jp = j * N_species
        do i = 1, N_species
          Cjp(i) = y(idx_jp + i)
        end do
      end if

      call source_terms(Cj, Rs)

      do i = 1, N_species
        ! conv: first-order upwind (u_vel > 0 in this project, so backward
        !       difference). diff: standard central-difference second derivative.
        conv = -u_vel * (Cj(i) - Cjm(i)) / dx
        diff =  D_ax  * (Cjp(i) - 2.0_DP * Cj(i) + Cjm(i)) / (dx * dx)
        dydt(idx + i) = conv + diff + Rs(i)
      end do
    end do
end subroutine model


! Dummy Jacobian: mf = 25 builds its own banded FD Jacobian, so this is
! never called. Symbol must exist for F77 linkage.
SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD)
    INTEGER          :: NEQ, ML, MU, NROWPD
    DOUBLE PRECISION :: T, Y(*), PD(NROWPD, *)
    if (.false.) PD(1,1) = T + Y(1) + real(ML+MU+NROWPD+NEQ, kind=kind(T))
END SUBROUTINE
