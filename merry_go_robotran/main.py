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
from postprocess.internal_forces import (scan_all_sections,
                                          find_critical_section)
from postprocess.dimensioning    import size_section, print_dimensioning_report
from postprocess.plots           import generate_all_plots
from postprocess.gonogo_comparison import compute_passenger_acceleration

# ═══════════════════════════════════════════════════════════════════════════ #
#  Configuration                                                              #
# ═══════════════════════════════════════════════════════════════════════════ #

N_SECTIONS       = 50         # number of cross-sections along pendulum
RTOL             = 1e-6       # ODE relative tolerance
ATOL             = 1e-9       # ODE absolute tolerance
GONOGO_CSV       = None       # path to MOOC reference CSV (None = skip overlay)
CONVERGENCE_ITER = 3          # max design-loop iterations (mass/inertia update)


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

    # ── Design iteration loop ───────────────────────────────────────── #
    for iteration in range(1, CONVERGENCE_ITER + 1):
        print(f"\n[main] ── Design iteration {iteration}/{CONVERGENCE_ITER} ──")

        # ── Step 2: Initial conditions ──────────────────────────────── #
        q0_full, qd0_full = _assemble_initial_conditions(bodies)

        # ── Step 3: Integrate ODE ───────────────────────────────────── #
        result = run_simulation(
            bodies, joints,
            q0_full, qd0_full,
            t_end=t_end,
            rtol=RTOL,
            atol=ATOL,
        )

        # ── Step 4: Internal forces scan ────────────────────────────── #
        print("[main] Scanning internal forces along pendulum_1...")
        scan = scan_all_sections(N_SECTIONS, result, bodies, joints)

        s_crit, idx_t_crit, Mb_max, idx_s_crit, _ = find_critical_section(scan)
        t_crit = result['t'][idx_t_crit]

        # Axial force at critical section and time
        N_max = np.linalg.norm(scan['F_history'][idx_t_crit, idx_s_crit])

        print(f"[main] Critical section: s* = {s_crit*1000:.1f} mm, "
              f"t* = {t_crit:.2f} s, Mb_max = {Mb_max:.1f} N·m")

        # ── Step 5: Dimensioning ─────────────────────────────────────── #
        print("[main] Sizing pendulum cross-section (SF = 10)...")
        r_o_min, dim_results = size_section(Mb_max, N_max)
        print_dimensioning_report(dim_results, s_crit)

        # ── Convergence check ────────────────────────────────────────── #
        # If this is not the last iteration, update the pendulum mass/inertia
        # in the body list to reflect the newly sized section, then re-run.
        if iteration < CONVERGENCE_ITER:
            _update_pendulum_section(bodies, dim_results)
            joints = make_all_joints(bodies)   # rebuild joints (unchanged but safe)
            print(f"[main] Updated pendulum section — re-running simulation...")
        else:
            print("[main] Convergence loop complete.")

    # ── Step 6: Go/No-Go validation ─────────────────────────────────── #
    print("\n[main] Computing passenger acceleration for Go/No-Go check...")
    a_pass = compute_passenger_acceleration(result)
    a_max  = a_pass.max()
    print(f"[main] Max passenger acceleration = {a_max:.2f} m/s²  "
          f"({a_max/9.81:.2f} g)")

    # ── Step 7: Generate plots ───────────────────────────────────────── #
    print("\n[main] Generating plots...")
    generate_all_plots(
        result,
        scan=scan,
        s_crit=s_crit,
        idx_s_crit=idx_s_crit,
        reference_csv=GONOGO_CSV,
    )

    print("\n[main] Simulation complete. Results saved to results/")


# ═══════════════════════════════════════════════════════════════════════════ #
#  Helper: update body properties after re-sizing                             #
# ═══════════════════════════════════════════════════════════════════════════ #

def _update_pendulum_section(bodies: list, dim_results: dict) -> None:
    """
    Update the mass and inertia of all 4 pendulum bodies to match the
    newly computed cross-section from dimensioning.

    Uses a solid-equivalent mass distribution for the hollow tube:
        m_new = ρ_steel * A * L_pendulum  (uniform distribution)
        I_tra_new = m_new * (3*r_o² + 3*r_i² + L²) / 12  (hollow cylinder)
        I_rev_new = m_new * (r_o² + r_i²) / 2
    """
    import data.carousel_data as cd

    rho_steel = 7850.0   # [kg/m³] structural steel density
    r_o  = dim_results['r_o']
    r_i  = dim_results['r_i']
    A    = dim_results['A_m2']
    L    = cd.L_pendulum

    m_new   = rho_steel * A * L
    I_tra   = m_new * (3*r_o**2 + 3*r_i**2 + L**2) / 12.0
    I_rev   = m_new * (r_o**2 + r_i**2) / 2.0

    for body in bodies:
        if body.name.startswith("pendulum"):
            body.mass = m_new
            body.inertia_com_local = np.diag([I_tra, I_tra, I_rev])
            print(f"[main] {body.name}: m = {m_new:.1f} kg, "
                  f"I_tra = {I_tra:.2f} kg·m²")


if __name__ == "__main__":
    main()






'''from merry_go_robotran.data.carousel_data import build_bodies
from merry_go_robotran.neri.joint import make_all_joints
from merry_go_robotran.neri.state import MBState
from merry_go_robotran.neri.forward import forward_pass
from merry_go_robotran.neri.backward import backward_pass
from merry_go_robotran.neri.joint import make_all_joints

import numpy as np

T_END = 30
TIMESTEP = 0.5
Q_0 = np.zeros(15)
QD_0 = np.zeros_like(Q_0)

def main():
    bodies = build_bodies()
    joints = make_all_joints(bodies)
    state = MBState(n_bodies=len(bodies), n_dof=15)  # allocate once
    q = Q_0
    qd = QD_0
    for t_act in np.arange(0,T_END, TIMESTEP):
        state.set_q(q, qd, t_act)  # load current state
        forward_pass(bodies, joints, state)  # fill ω, α, R ...
        backward_pass(bodies, joints, state)  # fill W, F, L, Q

if __name__ == "__main__":
    main()'''