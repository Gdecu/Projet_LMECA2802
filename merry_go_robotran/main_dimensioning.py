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

from merry_go_robotran.data.carousel_data import L_pendulum

# Allow running directly from the project root
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from data.carousel_data    import build_bodies_split, q0, qd0, L_pendulum
from neri.joint            import make_all_joints
from neri.assembly         import get_dof_indices, assemble_Mc, solve_independent_accelerations
from simulation.integrator import run_simulation,_damping_torques, _wind_force
from neri.state            import MBState
from neri.forward          import forward_pass
from neri.backward         import backward_pass

from merry_go_robotran.postprocess.plots import (
    plot_pole_omega,
    plot_pendulum_positions,
    plot_all_nacelle_angles,
    plot_pole,
    plot_nacelle_pendulum_grid,
    load_reference_data,
    plot_force,
    plot_force_vs_length,
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
    step = 0.01 * L_pendulum
    split_sections = np.arange(step, 3, step)
    total = len(split_sections)
    max_force_at = {
        'time': [],
        'max_force': [],
        'length' : [],
        'fx': [],
        'fy': [],
        'fz': [],
    }
    for i, len_section in enumerate(split_sections):
        bodies = build_bodies_split(len_section)
        joints = make_all_joints(bodies)

        idx_u, idx_c = get_dof_indices(bodies)
        n_dof = sum(b.n_dof for b in bodies)
        print(f"[main] Simulation [{i+1}/{total}] running...")

        # ── Step 2: Initial conditions ──────────────────────────────── #
        q0_full, qd0_full = _assemble_initial_conditions(bodies)
        # ── Step 3: Integrate ODE ───────────────────────────────────── #
        result = run_simulation(
            bodies, joints,
            q0_full, qd0_full,
            t_end=3,  # t_end, # shorter time for faster debugging
            rtol=RTOL,
            atol=ATOL,
            case='dimensioning'
        )
        Q_react_c, F_cut, L_cut = compute_split_reactions(result, bodies, joints)

        #plot_force(result, Q_react_c, F_cut, L_cut)
        max_force = 0
        for ii in range(len(Q_react_c[:,0])):
            norm_force = np.sqrt(Q_react_c[ii,0]**2 + Q_react_c[ii,1]**2 + Q_react_c[ii,2]**2)
            if norm_force > max_force:
                max_force = norm_force
                hap_at_i = ii
        max_force_at['time'].append(result['t'][hap_at_i])
        max_force_at['max_force'].append(max_force)
        max_force_at['length'].append(len_section)
        max_force_at['fx'].append(Q_react_c[hap_at_i,0])
        max_force_at['fy'].append(Q_react_c[hap_at_i,1])
        max_force_at['fz'].append(Q_react_c[hap_at_i,2])
        print(f"[main] Simulation [{i+1}/{total}] complete. Results saved to results/\n")
    plot_force_vs_length(max_force_at)
    print(f"\n[main] All simulations complete. Results saved to results/")

    #plot_pole_omega(result)
    #plot_pendulum_positions(result)
    #plot_all_nacelle_angles(result)  # or plot_nacelle_angles(result, nacelle_id=1) for a single one
    #ref = load_reference_data('Robotran/lmerry_go_final/plot/reference_data.txt')
    #plot_pole(result, ref=ref)
    #plot_nacelle_pendulum_grid(result, ref=ref)
    #plot_nacelle_pendulum_grid(result)
    #plot_pole(result)

def compute_split_reactions(result, bodies, joints,
                            wind_fn=_wind_force, damping_fn=_damping_torques,  # same callables as in integrator
                            g=np.array([0,0,-9.81])):
    """
    Post-process sweep: at every saved timestep, rebuild M, c, Q_ext, solve
    for qdd_u, then evaluate the 3 driven-joint reactions at lower_pend.

    Returns
    -------
    Q_react_c : (n_steps, 3)  internal moment about the 3 cardan axes [N·m]
    F_cut    : (n_steps, 3)  full internal force at O^lower_pend [N]   (bonus)
    L_cut    : (n_steps, 3)  full internal moment at O^lower_pend [N·m] (bonus)
    """
    t_arr      = result['t']
    q_hist     = result['q']
    qd_hist    = result['qd']
    idx_u      = result['idx_u']
    idx_c      = result['idx_c']

    # Locate the driven body once
    lower_id   = next(b.body_id for b in bodies if b.name == 'lower_pend')
    n_bodies   = len(bodies) + 1          # +1 for base
    n_dof      = sum(b.n_dof for b in bodies)

    n_steps    = len(t_arr)
    Q_react_c  = np.zeros((n_steps, 3))   # 3 cardan axes
    F_cut      = np.zeros((n_steps, 3))   # bonus: full 3D reaction force
    L_cut      = np.zeros((n_steps, 3))   # bonus: full 3D reaction moment

    # Scratch state reused across all timesteps
    state = MBState(n_bodies=n_bodies, n_dof=n_dof)
    body_name_to_id = {b.name: b.body_id for b in bodies}

    for k, t in enumerate(t_arr):
        # --- 1. Load current kinematic state (driven DOFs already at 0) ----
        state.q[:]  = q_hist[k]
        state.qd[:] = qd_hist[k]
        # --- 2. External loads at this time (same callables as integrator) -
        F_ext = wind_fn(t, n_bodies, body_name_to_id)       # (n_bodies, 3)
        L_ext = np.zeros((n_bodies, 3))    # no external torques here
        Q_ext = damping_fn(bodies, state, body_name_to_id)  # (n_dof,) — motor + damping

        # --- 3. Assemble M and c via the last trick ------------------------
        M, c_vec = assemble_Mc(bodies, joints, state, F_ext, L_ext, g)

        # --- 4. Solve independent accelerations; driven stay at zero -------
        qdd_c_vec = np.zeros(len(idx_c))
        qdd_u = solve_independent_accelerations(
            M, c_vec, idx_u, idx_c, qdd_c_vec, Q_ext
        )
        qdd = np.zeros(n_dof)
        qdd[idx_u] = qdd_u                 # qdd[idx_c] stays 0 by construction

        # --- 5. Generalised reaction on the 3 driven cardan axes -----------
        # Q_react_c = (M·qdd + c - Q_ext)[idx_c]    ← the dimensioning number
        Q_full = M @ qdd + c_vec - Q_ext
        Q_react_c[k] = Q_full[idx_c[-3:]]

        # --- 6. (Bonus) full 6-component reaction from the backward pass ---
        # After assemble_Mc, run one clean pass with qdd injected to read F,L
        state.qdd[:] = qdd
        forward_pass(bodies, joints, state, g)
        backward_pass(bodies, joints, state, F_ext, L_ext)
        F_cut[k] = state.F[lower_id]       # force transmitted parent→child
        L_cut[k] = state.L[lower_id]       # moment at O^lower_pend

    return Q_react_c, F_cut, L_cut

if __name__ == "__main__":
    main()