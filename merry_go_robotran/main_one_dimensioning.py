"""
main_one_dimensioning.py
------------------------
Run (or reload) a single dimensioning simulation at one cut position along
the arm, then plot the full time-domain internal loads at that section.

Unlike main_dimensioning.py (which sweeps many positions and keeps only the
peak force at each), this file keeps every timestep so you can see how the
load evolves over time.

Usage
-----
    # Run simulation at the default section (half arm length), save cache:
    python main_one_dimensioning.py

    # Run at a specific position [m]:
    python main_one_dimensioning.py --section 1.2

    # Re-use cached result (skips integration):
    python main_one_dimensioning.py --section 1.2
    (cache is found automatically; delete the .npz file to force re-run)

    # Force re-simulation even if cache exists:
    python main_one_dimensioning.py --section 1.2 --rerun
"""

import numpy as np
import sys
import os
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from data.carousel_data    import build_bodies_split, q0, qd0, L_pendulum
from neri.joint            import make_all_joints
from neri.assembly         import get_dof_indices
from simulation.integrator import run_simulation

from main_dimensioning import (
    compute_split_reactions, _assemble_initial_conditions,
    save_dimensioning_data, load_dimensioning_data, _DATA_FILE,
)

from merry_go_robotran.postprocess.plots import plot_force_vs_time

# ─── defaults ────────────────────────────────────────────────────────────────

RTOL         = 1e-6
ATOL         = 1e-9
T_END        = 3.0
DEFAULT_SECTION = 1e-6   # effectively 0 cm — safe epsilon for split geometry

_RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")


def _cache_path(len_section: float) -> str:
    return os.path.join(_RESULTS_DIR, f"one_dim_cache_{len_section:.4f}.npz")


# ─── simulation / cache ───────────────────────────────────────────────────────

def run_or_load(len_section: float, rerun: bool = False) -> tuple:
    """
    Return (t, Q_react_c, F_cut, L_cut) either from cache or by running the
    full NERi simulation + split-reaction post-processing.

    Parameters
    ----------
    len_section : float   cut position along the arm [m]
    rerun       : bool    if True, ignore any existing cache

    Returns
    -------
    t           : (n_steps,)
    Q_react_c   : (n_steps, 3)  generalised cardan reactions [N·m]
    F_cut       : (n_steps, 3)  physical forces at cut [N]
    L_cut       : (n_steps, 3)  physical moments at cut [N·m]
    """
    cache = _cache_path(len_section)

    if not rerun and os.path.exists(cache):
        print(f"[one_dim] Loading cached results from {cache}")
        data = np.load(cache)
        return data['t'], data['Q_react_c'], data['F_cut'], data['L_cut']

    print(f"[one_dim] Building system for section at {len_section:.4f} m...")
    bodies = build_bodies_split(len_section)
    joints = make_all_joints(bodies)

    q0_full, qd0_full = _assemble_initial_conditions(bodies)

    print(f"[one_dim] Running simulation (t_end={T_END} s)...")
    result = run_simulation(
        bodies, joints,
        q0_full, qd0_full,
        t_end=T_END,
        rtol=RTOL,
        atol=ATOL,
        case='dimensioning',
    )

    print("[one_dim] Computing split reactions...")
    Q_react_c, F_cut, L_cut = compute_split_reactions(result, bodies, joints)

    os.makedirs(_RESULTS_DIR, exist_ok=True)
    np.savez(cache,
             t=result['t'],
             Q_react_c=Q_react_c,
             F_cut=F_cut,
             L_cut=L_cut)
    print(f"[one_dim] Cache saved → {cache}")

    return result['t'], Q_react_c, F_cut, L_cut


# ─── entry point ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Time-domain internal loads at one cut position along the arm."
    )
    parser.add_argument(
        "--section", type=float, default=DEFAULT_SECTION, metavar="LENGTH",
        help=f"Cut position along the arm [m]  (default: {DEFAULT_SECTION:.3f} m)"
    )
    parser.add_argument(
        "--rerun", action="store_true",
        help="Ignore cached results and re-run the simulation."
    )
    args = parser.parse_args()

    t, Q_react_c, F_cut, L_cut = run_or_load(args.section, rerun=args.rerun)

    # ── Compute peak force (same logic as main_dimensioning.py loop) ──────── #
    max_force = 0
    hap_at_i = 0
    for ii in range(len(Q_react_c)):
        norm_force = np.sqrt(Q_react_c[ii, 0]**2 + Q_react_c[ii, 1]**2 + Q_react_c[ii, 2]**2)
        if norm_force > max_force:
            max_force = norm_force
            hap_at_i = ii

    fx, fy, fz = Q_react_c[hap_at_i]
    t_peak = t[hap_at_i]

    print(f"\n[one_dim] Peak force at section {args.section:.6f} m:")
    print(f"  max_force = {max_force:.4f} N")
    print(f"  fx = {fx:.6f} N,  fy = {fy:.6f} N,  fz = {fz:.6f} N")
    print(f"  at time   = {t_peak:.4f} s")

    # ── Prepend row to dimensioning_data.txt ─────────────────────────────── #
    new_row = {
        'length':    [args.section],
        'max_force': [max_force],
        'fx':        [fx],
        'fy':        [fy],
        'fz':        [fz],
        'time':      [t_peak],
    }
    if os.path.exists(_DATA_FILE):
        existing = load_dimensioning_data(_DATA_FILE)
        for key in new_row:
            new_row[key] = new_row[key] + existing[key]
    save_dimensioning_data(new_row)
    print(f"[one_dim] Row prepended to {_DATA_FILE}")

    print(f"[one_dim] Plotting time-domain loads at {args.section:.4f} m...")
    plot_force_vs_time(t, Q_react_c, args.section)
    print("[one_dim] Done.")
