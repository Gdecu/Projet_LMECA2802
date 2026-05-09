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

from data.carousel_data    import build_bodies_split, q0, qd0, t_end
from neri.joint            import make_all_joints
from neri.assembly         import get_dof_indices
from simulation.integrator import run_simulation

from merry_go_robotran.postprocess.plots import (
    plot_pole_omega,
    plot_pendulum_positions,
    plot_all_nacelle_angles,
    plot_pole,
    plot_nacelle_pendulum_grid,
    load_reference_data,
    plot_force,
    plot_comparison_sections
)

# ═══════════════════════════════════════════════════════════════════════════ #
#  Configuration                                                              #
# ═══════════════════════════════════════════════════════════════════════════ #

RTOL             = 1e-6       # ODE relative tolerance
ATOL             = 1e-9       # ODE absolute tolerance

SIGMA_Y = 235e6  # Limite élastique de l'acier [Pa] (S235)
K_SAFETY = 10     # Facteur de sécurité imposé

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
    
    all_results = []

    step = 0.5
    split_sections = np.arange(3*step, 3, step)
    total = len(split_sections)

    M_max_global = 0
    section_critique = 0


    for i, len_section in enumerate(split_sections):
        bodies = build_bodies_split(len_section)
        joints = make_all_joints(bodies)

        idx_u, idx_c = get_dof_indices(bodies)
        n_dof = sum(b.n_dof for b in bodies)
        print(f"[main] System: {len(bodies)} bodies, {n_dof} DOFs "
          f"({len(idx_u)} independent, {len(idx_c)} driven)")

        # ── Step 2: Initial conditions ──────────────────────────────── #
        q0_full, qd0_full = _assemble_initial_conditions(bodies)
        # ── Step 3: Integrate ODE ───────────────────────────────────── #
        result = run_simulation(
            bodies, joints,
            q0_full, qd0_full,
            t_end=10,  # t_end, # shorter time for faster debugging
            rtol=RTOL,
            atol=ATOL,
            case='dimensioning'
        )

        # ── Step 4: Post-process results for dimensioning ────────────────── #
        all_results.append(result)

        # Identification du moment maximum sur l'axe 1 pour cette section
        lower_body = next(b for b in result['bodies'] if b.name == 'lower_pend')
        idx_1 = lower_body.q_indices[0]
        M_max_local = np.max(np.abs(result['Q'][:, idx_1]))
        
        if M_max_local > M_max_global:
            M_max_global = M_max_local
            section_critique = len_section

        fig = plot_force(result, len_section)
        # Sauvegarde avec la longueur dans le nom du fichier
        filename = f"./results/forces_internes_{len_section}m.pdf"
        fig.savefig(filename)
        print(f"Graphique généré : {filename}")
        print(f"\n[main] Simulation [{i+1}/{total}] complete. Results saved to results/\n")

    print(f"\n[main] All simulations complete. Results saved to results/")

# --- Étape de calcul finale ---
    print("\n" + "="*40)
    print("      RÉSULTATS DE DIMENSIONNEMENT")
    print("="*40)
    print(f"Moment fléchissant max détecté : {M_max_global:.2f} N.m")
    print(f"Localisation de la section critique : {section_critique} m")
    
    # Calcul du rayon minimal : R >= ( (4 * M) / (pi * (sigma_y / k)) )^(1/3)
    sigma_admissible = SIGMA_Y / K_SAFETY
    R_min = ((4 * M_max_global) / (np.pi * sigma_admissible))**(1/3)
    
    print(f"Contrainte admissible : {sigma_admissible/1e6:.1f} MPa")
    print(f"Rayon minimal requis : {R_min*1000:.2f} mm")
    print(f"Diamètre minimal conseillé : {R_min*2000:.2f} mm")
    
    # Tracer la comparaison
    plot_comparison_sections(all_results, split_sections)


if __name__ == "__main__":
    main()