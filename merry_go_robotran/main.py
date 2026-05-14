"""
simulation/main.py
------------------
Main entry point for the Merry-Go-Round NERi simulation.

Execution pipeline
──────────────────
  1. Load system data           (carousel_data.build_bodies)
  2. Build joints               (joint.make_all_joints)
  3. Set initial conditions     (carousel_data.q0, qd0)
  4. Integrate ODE              (integrator.run_simulation)
  5. Post-process               (internal_forces, dimensioning, plots)
  6. Go/No-Go check             (gonogo_comparison)

Run from the project root:
    python -m simulation.main
or:
    python simulation/main.py
"""

import numpy as np
import sys
import os

# Allow running directly from the project root
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from data.carousel_data    import build_bodies, q0, qd0, t_end
from neri.joint            import make_all_joints
from neri.assembly         import get_dof_indices
from simulation.integrator import run_simulation

CACHE_FILE = os.path.join(os.path.dirname(__file__), "results", "simulation_cache.npz")

from merry_go_robotran.postprocess.plots import (
    plot_pole_omega,
    plot_pendulum_positions,
    plot_all_nacelle_angles,
    plot_pole,
    plot_nacelle_pendulum_grid,
    load_reference_data
)

# ═══════════════════════════════════════════════════════════════════════════ #
#  Configuration                                                              #
# ═══════════════════════════════════════════════════════════════════════════ #

RTOL             = 1e-6       # ODE relative tolerance
ATOL             = 1e-9       # ODE absolute tolerance

# ═══════════════════════════════════════════════════════════════════════════ #
#  Helper: build global initial condition vector                              #
# ═══════════════════════════════════════════════════════════════════════════ #

def _assemble_initial_conditions(bodies: list) -> tuple:
    """
    Map the named initial conditions (q0, qd0 from carousel_data) onto
    the global q-vector using each body's q_indices.

    Returns
    -------
    q0_full  : (n_dof,) initial positions
    qd0_full : (n_dof,) initial velocities
    """
    # Count total DOFs
    n_dof = sum(b.n_dof for b in bodies)
    q0_full  = np.zeros(n_dof)
    qd0_full = np.zeros(n_dof)

    for body in bodies:
        for k, (var_type, q_idx) in enumerate(zip(body.var_types, body.q_indices)):
            name = body.name
            if name in q0 and k < len(body.q_indices):
                q0_full[q_idx]  = q0.get(name, 0.0)
                qd0_full[q_idx] = qd0.get(name, 0.0)

    return q0_full, qd0_full


# ═══════════════════════════════════════════════════════════════════════════ #
#  Main simulation loop (with optional design iteration)                     #
# ═══════════════════════════════════════════════════════════════════════════ #

def main():
    print("\n" + "=" * 60)
    print("  MERRY-GO-ROUND — NERi Multibody Simulation")
    print("=" * 60 + "\n")

    # ── Step 1: Load system ──────────────────────────────────────────── #
    print("[main] Building bodies and joints...")
    bodies = build_bodies()
    joints = make_all_joints(bodies)

    idx_u, idx_c = get_dof_indices(bodies)
    n_dof = sum(b.n_dof for b in bodies)
    print(f"[main] System: {len(bodies)} bodies, {n_dof} DOFs "
          f"({len(idx_u)} independent, {len(idx_c)} driven)")

    # ── Step 2: Initial conditions ──────────────────────────────── #
    q0_full, qd0_full = _assemble_initial_conditions(bodies)

    # ── Step 3: Integrate ODE (or load cached result) ───────────── #
    if os.path.exists(CACHE_FILE):
        print(f"[main] Loading cached simulation from {CACHE_FILE}")
        data = np.load(CACHE_FILE)
        result = {
            't':      data['t'],
            'q':      data['q'],
            'qd':     data['qd'],
            'qdd':    data['qdd'],
            'Q':      data['Q'],
            'idx_u':  list(data['idx_u']),
            'idx_c':  list(data['idx_c']),
            'bodies': bodies,
            'joints': joints,
            'ctx':    {'t0_motor': float(data['t0_motor'])},
        }
    else:
        result = run_simulation(
            bodies, joints,
            q0_full, qd0_full,
            t_end=10,  # t_end, # shorter time for faster debugging
            rtol=RTOL,
            atol=ATOL,
        )
        os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)
        np.savez(CACHE_FILE,
                 t=result['t'],
                 q=result['q'],
                 qd=result['qd'],
                 qdd=result['qdd'],
                 Q=result['Q'],
                 idx_u=np.array(result['idx_u']),
                 idx_c=np.array(result['idx_c']),
                 t0_motor=np.array(result['ctx']['t0_motor']))
        print(f"[main] Simulation cached to {CACHE_FILE}")
    plot_pole_omega(result)
    plot_pendulum_positions(result)
    plot_all_nacelle_angles(result)  # or plot_nacelle_angles(result, nacelle_id=1) for a single one
    ref = load_reference_data('Robotran/lmerry_go_final/plot/reference_data.txt')
    plot_pole(result, ref=ref)
    plot_nacelle_pendulum_grid(result, ref=ref)
    #plot_nacelle_pendulum_grid(result)
    #plot_pole(result)

    print("\n[main] Simulation complete. Results saved to results/")

if __name__ == "__main__":
    main()