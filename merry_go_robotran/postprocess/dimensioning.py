"""
postprocess/dimensioning.py
----------------------------
Translates the internal forces from the NERi simulation into stresses,
then sizes the pendulum cross-section with a safety factor of 10.

Beam theory applied (Euler-Bernoulli)
──────────────────────────────────────
The pendulum is modelled as a hollow circular cross-section tube
(outer radius r_o, inner radius r_i = r_o * k_hollow).

At the critical cross-section s* the dominant stress is:

    σ_max = M_b * r_o / I_b  +  N / A

where
    M_b = |Mb_cut|        bending moment magnitude [N·m]
    N   = axial component of F_cut along pendulum axis [N]
    r_o = outer radius    [m]
    I_b = π/4 * (r_o⁴ - r_i⁴)   area moment of inertia [m⁴]
    A   = π * (r_o² - r_i²)      cross-section area [m²]

Sizing criterion (safety factor SF = 10):
    σ_max  ≤  σ_yield / SF

Given σ_yield (material property), SF, and the worst-case (M_b, N),
we solve for the minimum r_o satisfying the criterion.

Iteration strategy
──────────────────
1. Solve analytically for a solid section (r_i = 0) to get a lower bound.
2. Iterate r_o from that lower bound upward until σ_max ≤ σ_yield / SF,
   with the chosen hollow ratio k_hollow.

Public API
──────────
    compute_stress(Mb_mag, N_axial, r_o, k_hollow) → σ_bending, σ_axial, σ_total
    size_section(Mb_max, N_max, sigma_yield, SF, k_hollow) → r_o_min, results
    print_dimensioning_report(results)
"""

import numpy as np


# ═══════════════════════════════════════════════════════════════════════════ #
#  Material and safety constants                                              #
# ═══════════════════════════════════════════════════════════════════════════ #

SAFETY_FACTOR    = 10            # as required by the project brief
SIGMA_YIELD_PA   = 250e6         # [Pa]  structural steel S235 yield strength
K_HOLLOW_DEFAULT = 0.8           # r_i / r_o — standard structural tube ratio


# ═══════════════════════════════════════════════════════════════════════════ #
#  Section geometry helpers                                                   #
# ═══════════════════════════════════════════════════════════════════════════ #

def section_properties(r_o: float, k_hollow: float) -> tuple:
    """
    Compute cross-section area A and second moment of area I for a
    hollow circular tube.

    Parameters
    ----------
    r_o      : outer radius [m]
    k_hollow : r_i / r_o  (0 = solid, 1 = infinitely thin shell)

    Returns
    -------
    A   : cross-section area [m²]
    I_b : second moment of area (bending) [m⁴]
    r_i : inner radius [m]
    """
    r_i = k_hollow * r_o
    A   = np.pi * (r_o**2 - r_i**2)
    I_b = (np.pi / 4.0) * (r_o**4 - r_i**4)
    return A, I_b, r_i


# ═══════════════════════════════════════════════════════════════════════════ #
#  Stress calculation                                                         #
# ═══════════════════════════════════════════════════════════════════════════ #

def compute_stress(Mb_mag:   float,
                   N_axial:  float,
                   r_o:      float,
                   k_hollow: float = K_HOLLOW_DEFAULT
                   ) -> tuple:
    """
    Compute the maximum combined stress at the outer fibre of the section.

    Bending stress (maximum at outer fibre):
        σ_b = M_b * r_o / I_b

    Axial (normal) stress:
        σ_n = N / A   (positive = tension, negative = compression)

    Combined (worst case, same sign as bending):
        σ_total = |σ_b| + |σ_n|   (conservative superposition)

    Parameters
    ----------
    Mb_mag   : bending moment magnitude [N·m]
    N_axial  : axial force (along beam axis) [N]
    r_o      : outer radius [m]
    k_hollow : hollow ratio r_i / r_o

    Returns
    -------
    sigma_bending : float [Pa]
    sigma_axial   : float [Pa]
    sigma_total   : float [Pa]  (conservative sum)
    """
    A, I_b, _ = section_properties(r_o, k_hollow)

    sigma_bending = Mb_mag * r_o / I_b       # [Pa]
    sigma_axial   = abs(N_axial) / A         # [Pa] magnitude only (conservative)
    sigma_total   = sigma_bending + sigma_axial

    return sigma_bending, sigma_axial, sigma_total


# ═══════════════════════════════════════════════════════════════════════════ #
#  Section sizing                                                             #
# ═══════════════════════════════════════════════════════════════════════════ #

def size_section(Mb_max:      float,
                 N_max:       float,
                 sigma_yield: float  = SIGMA_YIELD_PA,
                 SF:          float  = SAFETY_FACTOR,
                 k_hollow:    float  = K_HOLLOW_DEFAULT,
                 r_o_min_m:   float  = 0.005,
                 r_o_max_m:   float  = 0.500,
                 n_iter:      int    = 10000
                 ) -> tuple:
    """
    Find the minimum outer radius r_o such that σ_total ≤ σ_yield / SF.

    Uses a linear sweep from r_o_min_m to r_o_max_m with n_iter steps.

    Parameters
    ----------
    Mb_max      : maximum bending moment [N·m]
    N_max       : maximum axial force [N]
    sigma_yield : material yield stress [Pa]
    SF          : safety factor (default 10)
    k_hollow    : hollow ratio r_i / r_o
    r_o_min_m   : lower bound for search [m]
    r_o_max_m   : upper bound for search [m]
    n_iter      : number of candidate radii to test

    Returns
    -------
    r_o_min  : minimum admissible outer radius [m]
    results  : dict with full sizing details
    """

    sigma_allow = sigma_yield / SF      # [Pa] allowable stress

    # Sweep outer radius from small to large; stop at first compliant r_o
    r_candidates = np.linspace(r_o_min_m, r_o_max_m, n_iter)
    r_o_solution = None

    for r_o in r_candidates:
        _, _, sigma_total = compute_stress(Mb_max, N_max, r_o, k_hollow)
        if sigma_total <= sigma_allow:
            r_o_solution = r_o
            break

    if r_o_solution is None:
        raise ValueError(
            f"No admissible section found up to r_o = {r_o_max_m*1000:.0f} mm. "
            f"Increase r_o_max_m or check load inputs."
        )

    # Recompute detailed stresses at solution
    A_sol, I_sol, r_i_sol = section_properties(r_o_solution, k_hollow)
    sigma_b, sigma_n, sigma_tot = compute_stress(Mb_max, N_max, r_o_solution, k_hollow)

    results = {
        # Sizing inputs
        'Mb_max':        Mb_max,
        'N_max':         N_max,
        'sigma_yield':   sigma_yield,
        'SF':            SF,
        'sigma_allow':   sigma_allow,
        'k_hollow':      k_hollow,
        # Solution geometry
        'r_o':           r_o_solution,
        'r_i':           r_i_sol,
        'd_outer_mm':    r_o_solution * 2000.0,     # [mm] outer diameter
        'd_inner_mm':    r_i_sol      * 2000.0,     # [mm] inner diameter
        'wall_mm':       (r_o_solution - r_i_sol) * 1000.0,  # [mm]
        'A_m2':          A_sol,
        'I_m4':          I_sol,
        # Stress breakdown
        'sigma_bending': sigma_b,
        'sigma_axial':   sigma_n,
        'sigma_total':   sigma_tot,
        'utilisation':   sigma_tot / sigma_allow,   # should be ≤ 1.0
    }

    return r_o_solution, results


# ═══════════════════════════════════════════════════════════════════════════ #
#  Reporting                                                                  #
# ═══════════════════════════════════════════════════════════════════════════ #

def print_dimensioning_report(results: dict, s_crit: float = None) -> None:
    """
    Print a formatted dimensioning summary to stdout.

    Parameters
    ----------
    results : dict returned by size_section()
    s_crit  : critical abscissa [m], if known
    """
    print()
    print("=" * 60)
    print("  PENDULUM DIMENSIONING REPORT")
    print("=" * 60)

    if s_crit is not None:
        print(f"  Critical section abscissa  s* = {s_crit*1000:.1f} mm "
              f"from nacelle joint")

    print(f"\n  Loading (worst-case over simulation)")
    print(f"    Bending moment  Mb_max = {results['Mb_max']:.1f} N·m")
    print(f"    Axial force     N_max  = {results['N_max']:.1f} N")

    print(f"\n  Material & Safety")
    print(f"    Yield stress  σ_y  = {results['sigma_yield']/1e6:.0f} MPa")
    print(f"    Safety factor SF   = {results['SF']}")
    print(f"    Allowable     σ_a  = {results['sigma_allow']/1e6:.2f} MPa")

    print(f"\n  Required hollow circular section  (k = r_i/r_o = {results['k_hollow']})")
    print(f"    Outer diameter  d_o = {results['d_outer_mm']:.1f} mm")
    print(f"    Inner diameter  d_i = {results['d_inner_mm']:.1f} mm")
    print(f"    Wall thickness  t   = {results['wall_mm']:.1f} mm")

    print(f"\n  Stress check at critical section")
    print(f"    σ_bending = {results['sigma_bending']/1e6:.2f} MPa")
    print(f"    σ_axial   = {results['sigma_axial']/1e6:.2f} MPa")
    print(f"    σ_total   = {results['sigma_total']/1e6:.2f} MPa")
    print(f"    Utilisation = {results['utilisation']*100:.1f}%  "
          f"({'OK' if results['utilisation'] <= 1.0 else 'FAIL'})")
    print("=" * 60)
    print()