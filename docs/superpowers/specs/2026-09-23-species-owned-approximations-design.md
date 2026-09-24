# Species-owned approximations, and Cold NEDE

**Status (2026-09-23).** Shipped with Cold NEDE (`species/nede_species.*`), its first
and only user. The module's own approximations (tca, rsa, ufa, ncdmfa, the idm ones)
are unchanged and still owned by `PerturbationsModule` by name.

## Problem

Cold New Early Dark Energy (H0 Olympics, arXiv:2107.10291 sec. 2.4.5; Niedermann &
Sloth, arXiv:2006.06686) is a cosmological constant with no perturbations until a
phase transition at a_*, and a fluid with delta and theta after it. Before a_*, the
perturbation of a trigger field is evolved instead, and at a_* the fluid is seeded from
it by junction conditions. So the set of evolved variables changes at a time that
belongs to the species, and the state is transformed across it. That is exactly what
the module's approximation machinery does for tca or ncdmfa, but there it is owned by
the module and names each approximation.

## Approach

A species declares how many approximations it owns and answers which regime it is in:

    virtual int ApproximationCount() const;                            // default 0
    virtual int ApproximationRegimeAt(int which, double k, const double* pvecback) const;
    int ApproximationRegime(const perturb_workspace* ppw, int which = 0) const;

The module needed two touch points **[C]**. `perturb_workspace_init` reserves the
slots after its own and records each species' first slot in
`ppw->species_approx_index`, a vector parallel to `all_species_` like
`ppw->species_scratch`. `perturb_approximations` asks each species for its regime at
(k, tau). Everything downstream was already generic over `ap_size` and needed nothing
**[C]**: counting the intervals, bisecting for the switch times, rebuilding the
perturbation vector, and `CopyPerturbationsAcrossSwitch`, which is where the junction
conditions live. `RegisterPerturbationIndices` reads the regime through
`ApproximationRegime(ppw)`.

Regimes are numbered chronologically, as the switch bisection assumes.

## Alternatives rejected

- **The species stores its slot index** (the prototype on `h0-olympics-plugins`,
  `8f96e8b9`). `perturb_workspace_init` runs inside the concurrent per-k tasks, so
  every task wrote the same member of a shared species: a data race, benign in practice
  but undefined. The index belongs to the workspace.
- **A context struct** for `ApproximationRegimeAt`, carrying k, tau, a, aH, pvecback,
  pvecthermo and ppr. NEDE reads a single number, `a` from pvecback. The struct can
  come when an approximation needs more.
- **A regime count per approximation**, used only to check that it is at least two.

## Limits before moving the module's approximations across

- **Scalars only.** Slots are reserved in scalar workspaces, and species register
  perturbations only for scalars.
- **NEDE is the easy case.** Its switch is k-independent and depends on `a` alone.
  tca reads kappa' and the baryon-photon ratio, rsa depends on tca having finished,
  and the module's approximations are entangled with each other and with cross-species
  state. This interface does not express that yet.

## A second hook NEDE needed

f_NEDE is NEDE's share of rho_tot at a_*, where H = H_* = H_over_m m_NEDE, so a_* is
where everything but NEDE, curvature included, contributes H_*^2 (1 - f_NEDE) to H^2. That has to be known at
construction: NEDE's density today enters the budget closure. So
`BaseSpecies::BackgroundDensityOverH0Sq(a, H0)` gives a species' density without a
background table. The default scales `GetOmega0()` by the energy type, which is exact
for a constant equation of state, and refuses `EnergyType::Other`. NCDM overrides it
from its momentum quadrature. That override matters: a 0.06 eV neutrino is still
relativistic at a_*, and counting it as matter costs 8% on rho_tot there **[P]** (the
prototype's measurement). NEDE is built last, so it sees every non-closure species. The
closure species, dark energy, is negligible at a_*.

## Verification

- LambdaCDM is bit-identical to master in lensed C_l and P(k) **[M]**.
- rho_NEDE = f_NEDE H_*^2 before the transition to 1e-12, and the other species reach
  H_*^2 (1 - f_NEDE) at z_* to 8e-6, with and without a massive neutrino **[M]**.
  A python regression test pins both.
- **Against TriggerCLASS** (flo1984/TriggerCLASS, built from the local checkout),
  f_NEDE = 0.1, w = c_s^2 = 2/3, the same explicit cosmology in both codes **[M]**:
  - Shooting on TriggerCLASS's z_decay_NEDE = 4589.54, the z_* CLASS++ gives for
    m_NEDE = 316.2278/Mpc, returns 316.2535/Mpc (8e-5 apart).
  - NEDE/LambdaCDM agrees to 2.1e-3 at most in TT (rms 2.5e-4) and 2.3e-3 in EE
    (rms 5.7e-4), on an effect of 11% in TT and 23% in EE.
  - A python regression test pins the TriggerCLASS ratios at l = 500, 1000, 1500 to
    3e-3; the prototype misses them by up to 1.2e-2.
- **The prototype missed that bar.** Its `CopyPerturbationsAcrossSwitch` carried NEDE's
  variables only across its own switch, so every later switch (tca, rsa, ...) zeroed
  them. Against TriggerCLASS it was off by 1.3% in TT and 12% in EE **[M]**. With that
  copy disabled, the port reproduces the prototype to 1e-4 **[M]**, which pins the
  difference on this bug. The prototype also lacked NEDE's p', and did not refuse
  Newtonian gauge, where its trigger equation misses the metric terms.
- f_NEDE = 1e-8 reproduces LambdaCDM to 1e-4 in C_l **[M]**. The residual is the
  integrator restarting at the extra switch.

## Not done

TriggerCLASS shoots on m_NEDE to hit a target z_*; here z_* is derived
(`NEDE.z_star`). The trigger is left out of the budget, which is the subdominant limit
the H0 Olympics assumes; TriggerCLASS can also carry it as dark matter. NEDE publishes
no transfer source. TriggerCLASS's optional shear viscosity defaults to zero, and so
the shear is omitted here.
