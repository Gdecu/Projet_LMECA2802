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

from main_dimensioning import compute_split_reactions, _assemble_initial_conditions

from merry_go_robotran.postprocess.plots import plot_force_vs_time

# ─── defaults ────────────────────────────────────────────────────────────────

RTOL         = 1e-6
ATOL         = 1e-9
T_END        = 3.0
DEFAULT_SECTION = L_pendulum / 2.0   # mid-arm by default

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

    print(f"[one_dim] Plotting time-domain loads at {args.section:.4f} m...")
    plot_force_vs_time(t, Q_react_c, args.section)
    print("[one_dim] Done.")
