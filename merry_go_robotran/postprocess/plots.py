"""
postprocess/plots.py
--------------------
All visualisation routines for the merry-go-round simulation.

Plots produced
──────────────
1. pole_kinematics  — pole angular position and velocity vs. time
2. pendulum_angles  — all 4 pendulum angles vs. time
3. bending_moment_envelope — |Mb(s,t)| as a 2D colour map
4. critical_section_history — Mb(t) at the critical abscissa s*
5. gonogo_comparison  — overlaid against the MOOC reference curve

All functions accept the `result` dict from integrator.run_simulation()
and the optional `scan` dict from internal_forces.scan_all_sections().
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")           # non-interactive backend (safe for all envs)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

import merry_go_robotran.data.carousel_data as cd

# Output folder for saved figures
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(RESULTS_DIR, exist_ok=True)


### New plots

# ═══════════════════════════════════════════════════════════════════════════ #
#  Plot A — Pole angular velocity (ω_pole about I3)                          #
# ═══════════════════════════════════════════════════════════════════════════ #

def plot_pole_omega(result: dict) -> plt.Figure:
    """
    Angular velocity of the main pole about I3 vs. time [rad/s].
    Highlights the motor shut-down threshold with a dashed line.
    """
    t       = result['t']
    bodies  = result['bodies']
    q_idx   = _find_q_index(bodies, 'pole', dof=0)   # pole spin = DOF 0

    omega_pole = result['qd'][:, q_idx]              # [rad/s]

    fig, ax = plt.subplots(figsize=(10, 4))
    fig.suptitle("Main Pole — Angular Velocity about I₃", fontweight='bold')

    ax.plot(t, omega_pole, color='steelblue', linewidth=1.8, label='ω_pole')
    ax.axhline(cd.omega_pole_threshold, ls='--', color='crimson', linewidth=1.2,
               label=f'Motor threshold  ({cd.omega_pole_threshold} rad/s)')

    ax.set_xlabel("Time  [s]")
    ax.set_ylabel("ω_pole  [rad/s]")
    ax.legend()
    ax.grid(True, alpha=0.35)
    fig.tight_layout()

    _save(fig, "pole_omega.pdf")
    return fig


# ═══════════════════════════════════════════════════════════════════════════ #
#  Plot B — All pendulum angular positions                                    #
# ═══════════════════════════════════════════════════════════════════════════ #

# Colour palette: one distinct colour per sub-system
_PEND_COLORS  = ['steelblue', 'crimson', 'seagreen', 'darkorange']
_PEND_LABELS  = ['Pendulum 1', 'Pendulum 2 (rusted)', 'Pendulum 3', 'Pendulum 4']

def plot_pendulum_positions(result: dict) -> plt.Figure:
    """
    Angular position of all 4 pendulum hinges vs. time [deg].
    Pendulum 2 (rusted hinge, high damping) is rendered with a dashed line
    to make its distinct behaviour immediately apparent.
    """
    t      = result['t']
    bodies = result['bodies']

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle("Pendulums — Angular Position", fontweight='bold')

    for k in range(1, 5):
        q_idx = _find_q_index(bodies, f'pendulum_{k}', dof=0)
        theta = np.degrees(result['q'][:, q_idx])   # convert rad → deg

        linestyle = '--' if k == 2 else '-'          # dashed for rusted hinge
        ax.plot(t, theta,
                color=_PEND_COLORS[k - 1],
                linestyle=linestyle,
                linewidth=1.6,
                label=_PEND_LABELS[k - 1])

    ax.set_xlabel("Time  [s]")
    ax.set_ylabel("θ_pendulum  [deg]")
    ax.legend()
    ax.grid(True, alpha=0.35)
    fig.tight_layout()

    _save(fig, "pendulum_positions.pdf")
    return fig


# ═══════════════════════════════════════════════════════════════════════════ #
#  Plot C — Nacelle Cardan angles (both DOFs on one plot, per nacelle)       #
# ═══════════════════════════════════════════════════════════════════════════ #

def plot_nacelle_angles(result: dict, nacelle_id: int = 1) -> plt.Figure:
    """
    Both Cardan angles (DOF 0 and DOF 1) of a single nacelle vs. time [deg].
    Call once per nacelle, or loop over nacelle_id ∈ {1,2,3,4}.

    Parameters
    ----------
    nacelle_id : int
        Sub-system index (1–4).
    """
    t      = result['t']
    bodies = result['bodies']

    # Retrieve both q-indices for the 2-DOF Cardan joint
    nacelle_body = next(b for b in bodies if b.name == f'nacelle_{nacelle_id}')
    q_idx_0, q_idx_1 = nacelle_body.q_indices[0], nacelle_body.q_indices[1]

    theta_0 = np.degrees(result['q'][:, q_idx_0])   # 1st Cardan angle [deg]
    theta_1 = np.degrees(result['q'][:, q_idx_1])   # 2nd Cardan angle [deg]

    fig, ax = plt.subplots(figsize=(10, 4))
    fig.suptitle(f"Nacelle {nacelle_id} — Cardan Angles", fontweight='bold')

    ax.plot(t, theta_0, color='steelblue', linewidth=1.6, label='Cardan angle 1 (⊥ arm & pole axes)')
    ax.plot(t, theta_1, color='darkorange', linewidth=1.6, linestyle='--',
            label='Cardan angle 2 (arm–pendulum plane)')

    ax.set_xlabel("Time  [s]")
    ax.set_ylabel("θ_nacelle  [deg]")
    ax.legend()
    ax.grid(True, alpha=0.35)
    fig.tight_layout()

    _save(fig, f"nacelle_{nacelle_id}_angles.pdf")
    return fig


def plot_all_nacelle_angles(result: dict) -> plt.Figure:
    """
    Both Cardan angles for all 4 nacelles on a single 4-panel figure.
    Gives a direct side-by-side comparison across sub-systems.
    """
    t      = result['t']
    bodies = result['bodies']

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True, sharey=True)
    fig.suptitle("All Nacelles — Cardan Angles", fontweight='bold')

    for k, ax in enumerate(axes.flat, start=1):
        nacelle_body = next(b for b in bodies if b.name == f'nacelle_{k}')
        q0, q1 = nacelle_body.q_indices[0], nacelle_body.q_indices[1]

        ax.plot(t, np.degrees(result['q'][:, q0]),
                color='steelblue', linewidth=1.4, label='Cardan 1')
        ax.plot(t, np.degrees(result['q'][:, q1]),
                color='darkorange', linewidth=1.4, linestyle='--', label='Cardan 2')

        ax.set_title(f'Nacelle {k}', fontsize=10)
        ax.grid(True, alpha=0.35)
        ax.set_ylabel("θ  [deg]")

    # Shared x-label and legend on bottom row only
    for ax in axes[1]:
        ax.set_xlabel("Time  [s]")
    axes[0, 0].legend(fontsize=8)

    fig.tight_layout()
    _save(fig, "all_nacelles_angles.pdf")
    return fig

# ═══════════════════════════════════════════════════════════════════════════ #
#  Utility                                                                    #
# ═══════════════════════════════════════════════════════════════════════════ #

def _save(fig: plt.Figure, filename: str) -> str:
    """Save figure to results/ and return the full path."""
    path = os.path.join(RESULTS_DIR, filename)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    print(f"[plots] Saved: {path}")
    return path


def _find_q_index(bodies: list, body_name: str, dof: int = 0) -> int:
    """Return the global q-index for a named body's k-th DOF."""
    body = next(b for b in bodies if b.name == body_name)
    return body.q_indices[dof]


# ═══════════════════════════════════════════════════════════════════════════ #
#  Plot 1 — Pole kinematics                                                   #
# ═══════════════════════════════════════════════════════════════════════════ #

def plot_pole_kinematics(result: dict) -> plt.Figure:
    """
    Angular position θ_pole and angular velocity ω_pole vs. time.
    Includes a horizontal dashed line at the motor threshold (0.8 rad/s).
    """
    t      = result['t']
    bodies = result['bodies']
    q_idx  = _find_q_index(bodies, 'pole')

    theta = result['q'][:,  q_idx]   # [rad]
    omega = result['qd'][:, q_idx]   # [rad/s]

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.suptitle("Pole — Main rotation kinematics", fontweight='bold')

    axes[0].plot(t, np.degrees(theta), color='steelblue')
    axes[0].set_ylabel("θ_pole  [deg]")
    axes[0].grid(True, alpha=0.4)

    axes[1].plot(t, omega, color='darkorange')
    axes[1].axhline(cd.omega_pole_threshold, ls='--', color='red',
                    label=f'Motor threshold ({cd.omega_pole_threshold} rad/s)')
    axes[1].set_ylabel("ω_pole  [rad/s]")
    axes[1].set_xlabel("Time  [s]")
    axes[1].legend()
    axes[1].grid(True, alpha=0.4)

    fig.tight_layout()
    _save(fig, "pole_kinematics.png")
    return fig


# ═══════════════════════════════════════════════════════════════════════════ #
#  Plot 2 — Pendulum angles                                                   #
# ═══════════════════════════════════════════════════════════════════════════ #

def plot_pendulum_angles(result: dict) -> plt.Figure:
    """
    All 4 pendulum hinge angles vs. time (degrees).
    Pendulum 2 (rusted hinge) is highlighted in red.
    """
    t      = result['t']
    bodies = result['bodies']

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle("Pendulum hinge angles", fontweight='bold')

    colors  = ['steelblue', 'red', 'forestgreen', 'mediumpurple']
    styles  = ['-', '-', '--', '--']

    for k, pname in enumerate(['pendulum_1','pendulum_2','pendulum_3','pendulum_4']):
        try:
            q_idx = _find_q_index(bodies, pname)
        except StopIteration:
            continue
        theta_deg = np.degrees(result['q'][:, q_idx])
        label = pname.replace('_', ' ').capitalize()
        if pname == 'pendulum_2':
            label += ' (rusted)'
        ax.plot(t, theta_deg, color=colors[k], ls=styles[k], label=label)

    ax.set_xlabel("Time  [s]")
    ax.set_ylabel("Pendulum angle  [deg]")
    ax.legend()
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    _save(fig, "pendulum_angles.png")
    return fig


# ═══════════════════════════════════════════════════════════════════════════ #
#  Plot 3 — Bending moment envelope (2D colour map)                           #
# ═══════════════════════════════════════════════════════════════════════════ #

def plot_bending_moment_envelope(scan: dict, result: dict) -> plt.Figure:
    """
    Colour map of  |Mb(s, t)|  over the full simulation.
    X-axis = time [s],  Y-axis = section abscissa s [m].
    """
    Mb_norm = scan['Mb_norm']           # (n_steps, n_sections)
    s_array = scan['s_array']
    t       = result['t']

    fig, ax = plt.subplots(figsize=(11, 5))
    fig.suptitle("Bending moment envelope — pendulum_1", fontweight='bold')

    cmap   = cm.plasma
    levels = np.linspace(0, Mb_norm.max(), 64)

    cf = ax.contourf(t, s_array, Mb_norm.T, levels=levels, cmap=cmap)
    fig.colorbar(cf, ax=ax, label="|Mb|  [N·m]")

    ax.set_xlabel("Time  [s]")
    ax.set_ylabel("Section abscissa  s  [m]  (0 = nacelle end)")
    ax.grid(True, alpha=0.2, color='white')

    fig.tight_layout()
    _save(fig, "bending_envelope.png")
    return fig


# ═══════════════════════════════════════════════════════════════════════════ #
#  Plot 4 — Bending moment history at critical section                        #
# ═══════════════════════════════════════════════════════════════════════════ #

def plot_critical_section_history(scan: dict,
                                  result: dict,
                                  idx_s_crit: int,
                                  s_crit: float) -> plt.Figure:
    """
    Time history of  |Mb|  at the dynamically critical cross-section s*.
    """
    Mb_norm = scan['Mb_norm'][:, idx_s_crit]   # (n_steps,) at s*
    t       = result['t']

    fig, ax = plt.subplots(figsize=(10, 4))
    fig.suptitle(f"Bending moment at critical section  s* = {s_crit*1000:.0f} mm",
                 fontweight='bold')

    ax.plot(t, Mb_norm, color='firebrick')
    ax.axhline(Mb_norm.max(), ls='--', color='black', alpha=0.6,
               label=f'Max = {Mb_norm.max():.0f} N·m')
    ax.set_xlabel("Time  [s]")
    ax.set_ylabel("|Mb(s*)|  [N·m]")
    ax.legend()
    ax.grid(True, alpha=0.4)

    fig.tight_layout()
    _save(fig, "critical_section_Mb.png")
    return fig


# ═══════════════════════════════════════════════════════════════════════════ #
#  Plot 5 — Go/No-Go validation                                               #
# ═══════════════════════════════════════════════════════════════════════════ #

def plot_gonogo_comparison(result: dict,
                           reference_csv: str = None) -> plt.Figure:
    """
    Overlay simulation results against the MOOC Go/No-Go reference curve.

    The MOOC validation metric is the magnitude of the absolute acceleration
    experienced by the passenger in nacelle_1, plotted vs. time.

    If `reference_csv` is provided it must be a two-column file: t, |a|.
    Otherwise only the simulation curve is plotted.

    Parameters
    ----------
    result        : integrator output dict
    reference_csv : path to MOOC reference data file (optional)
    """
    from merry_go_robotran.postprocess.gonogo_comparison import compute_passenger_acceleration

    t    = result['t']
    a_pg = compute_passenger_acceleration(result)   # (n_steps,) [m/s²]

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle("Go/No-Go validation — passenger acceleration", fontweight='bold')

    ax.plot(t, a_pg, color='steelblue', label='Simulation (this work)')

    if reference_csv is not None and os.path.isfile(reference_csv):
        ref = np.loadtxt(reference_csv, delimiter=',')
        ax.plot(ref[:, 0], ref[:, 1], 'r--', label='MOOC reference')

    ax.set_xlabel("Time  [s]")
    ax.set_ylabel("|a_passenger|  [m/s²]")
    ax.legend()
    ax.grid(True, alpha=0.4)

    fig.tight_layout()
    _save(fig, "gonogo_comparison.png")
    return fig


# ═══════════════════════════════════════════════════════════════════════════ #
#  Convenience: generate all plots at once                                    #
# ═══════════════════════════════════════════════════════════════════════════ #

def generate_all_plots(result: dict,
                       scan:   dict = None,
                       s_crit: float = None,
                       idx_s_crit: int = None,
                       reference_csv: str = None) -> list:
    """
    Generate and save all available plots.  Returns list of Figure objects.
    """
    figs = []

    figs.append(plot_pole_kinematics(result))
    figs.append(plot_pendulum_angles(result))

    if scan is not None:
        figs.append(plot_bending_moment_envelope(scan, result))
        if idx_s_crit is not None and s_crit is not None:
            figs.append(plot_critical_section_history(scan, result,
                                                       idx_s_crit, s_crit))

    figs.append(plot_gonogo_comparison(result, reference_csv))

    plt.close('all')
    return figs

# ═══════════════════════════════════════════════════════════════════════════ #
#  Plot P — Main pole kinematics: θ, ω, α on a single 1×3 figure             #
# ═══════════════════════════════════════════════════════════════════════════ #

def plot_pole(result: dict, ref: dict = None) -> plt.Figure:
    """
    Three-panel figure (1 row × 3 columns) for the main pole.
    If `ref` is provided (output of load_reference_data), the Robotran
    reference curves are overlaid in dashed red on each panel.
    """
    t      = result['t']
    bodies = result['bodies']
    q_idx  = _find_q_index(bodies, 'pole', dof=0)

    theta_deg = np.degrees(result['q']  [:, q_idx])
    omega_deg = np.degrees(result['qd'] [:, q_idx])
    alpha_deg = np.degrees(result['qdd'][:, q_idx])

    omega_thresh_deg = np.degrees(cd.omega_pole_threshold)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharex=True)
    fig.suptitle("Main Pole — Kinematics", fontweight='bold')

    signals   = [theta_deg,        omega_deg,         alpha_deg       ]
    ylabels   = ["θ_pole  [deg]",  "ω_pole  [deg/s]", "α_pole  [deg/s²]"]
    titles    = ["Angular Position","Angular Velocity","Angular Acceleration"]
    ref_keys  = ['pole_theta',     'pole_omega',      'pole_alpha'    ]
    colors    = ['steelblue',      'darkorange',      'seagreen'      ]

    for ax, sig, ylabel, title, rkey, color in zip(
            axes, signals, ylabels, titles, ref_keys, colors):

        ax.plot(t, sig, color=color, linewidth=1.6, label='Simulation')

        # Overlay reference curve if provided
        if ref is not None:
            ax.plot(ref['t'], ref[rkey], color='crimson', linewidth=1.2,
                    linestyle='--', label='Robotran ref.')

        # Motor threshold on velocity panel only
        if rkey == 'pole_omega':
            ax.axhline(omega_thresh_deg, ls=':', color='black', linewidth=1.0,
                       label=f'Motor threshold  ({omega_thresh_deg:.1f} deg/s)')

        ax.set_title(title)
        ax.set_xlabel("Time  [s]")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.35)

    fig.tight_layout()
    _save(fig, "pole_kinematics_full.pdf")
    return fig

# ═══════════════════════════════════════════════════════════════════════════ #
#  Plot Q — 4×3 grid: hinge + both Cardan angles for every nacelle           #
# ═══════════════════════════════════════════════════════════════════════════ #

# Column layout for the grid
_GRID_COL_LABELS  = ['Hinge Angle  [deg]',
                     'Cardan Angle 1  [deg]',
                     'Cardan Angle 2  [deg]']
_GRID_COL_TITLES  = ['Pendulum Hinge',
                     'Nacelle Cardan 1',
                     'Nacelle Cardan 2']
_GRID_COL_COLORS  = ['steelblue', 'darkorange', 'seagreen']

def plot_nacelle_pendulum_grid(result: dict, ref: dict = None) -> plt.Figure:
    """
    Large 4×3 figure. If `ref` is provided, Robotran reference curves are
    overlaid in dashed red on every cell.
    """
    t      = result['t']
    bodies = result['bodies']

    N_ROWS, N_COLS  = 4, 3
    ref_col_keys    = ['hinge', 'cadran1', 'cadran2']   # keys inside ref dict

    fig, axes = plt.subplots(N_ROWS, N_COLS,
                             figsize=(16, 12),
                             sharex=True,
                             sharey='col')
    fig.suptitle("Nacelle & Pendulum Angles — All Sub-Systems", fontweight='bold',
                 fontsize=13)

    for row in range(N_ROWS):
        nacelle_idx  = row + 1

        # Retrieve simulation signals
        hinge_idx    = _find_q_index(bodies, f'pendulum_{nacelle_idx}', dof=0)
        theta_hinge  = np.degrees(result['q'][:, hinge_idx])

        nacelle_body = next(b for b in bodies if b.name == f'nacelle_{nacelle_idx}')
        theta_c1     = np.degrees(result['q'][:, nacelle_body.q_indices[0]])
        theta_c2     = np.degrees(result['q'][:, nacelle_body.q_indices[1]])

        sim_signals  = [theta_hinge, theta_c1, theta_c2]

        for col, (signal, color, ylabel, rkey) in enumerate(zip(
                sim_signals, _GRID_COL_COLORS, _GRID_COL_LABELS, ref_col_keys)):

            ax = axes[row, col]
            ax.plot(t, signal, color=color, linewidth=1.4, label='Simulation')

            # Overlay reference curve (column `row` of the (n,4) ref array)
            if ref is not None:
                ref_signal = np.degrees(ref[rkey][:, row])   # convert rad→deg
                ax.plot(ref['t'], ref_signal, color='crimson', linewidth=1.0,
                        linestyle='--', label='Robotran ref.')

            ax.grid(True, alpha=0.35)

            if col == 0:
                ax.set_ylabel(f"Sub-system {nacelle_idx}\n{ylabel}", fontsize=8)
            else:
                ax.set_ylabel(ylabel, fontsize=8)

            if row == 0:
                ax.set_title(_GRID_COL_TITLES[col], fontsize=10, fontweight='bold')
                ax.legend(fontsize=7)   # legend once per column, at top

            if row == N_ROWS - 1:
                ax.set_xlabel("Time  [s]", fontsize=8)

    fig.tight_layout()
    _save(fig, "nacelle_pendulum_grid.pdf")
    return fig


# ═══════════════════════════════════════════════════════════════════════════ #
#  Reference data I/O                                                         #
# ═══════════════════════════════════════════════════════════════════════════ #

# Column indices inside the saved .txt file (matches save_reference_data order)
_REF_COL = {
    't':           0,
    'pole_theta':  1,  'pole_omega':  2,  'pole_alpha':  3,
    'hinge':      [4,  5,  6,  7],        # indices for sub-systems 1-4
    'cadran1':    [8,  9,  10, 11],
    'cadran2':    [12, 13, 14, 15],
}

def load_reference_data(filepath: str) -> dict:
    """
    Load the reference arrays saved by save_reference_data().

    Returns
    -------
    ref : dict with keys:
        't'          : (n,)    time vector [s]
        'pole_theta' : (n,)    pole angular position  [deg]
        'pole_omega' : (n,)    pole angular velocity  [deg/s]
        'pole_alpha' : (n,)    pole angular accel.    [deg/s²]
        'hinge'      : (n, 4)  pendulum hinge angles  [rad]
        'cadran1'    : (n, 4)  nacelle Cardan-1 angles [rad]
        'cadran2'    : (n, 4)  nacelle Cardan-2 angles [rad]
    """
    data = np.loadtxt(filepath, comments='#')   # skip header lines

    ref = {
        't':          data[:, _REF_COL['t']],
        'pole_theta': data[:, _REF_COL['pole_theta']],
        'pole_omega': data[:, _REF_COL['pole_omega']],
        'pole_alpha': data[:, _REF_COL['pole_alpha']],
        'hinge':      data[:, _REF_COL['hinge']],      # shape (n, 4)
        'cadran1':    data[:, _REF_COL['cadran1']],
        'cadran2':    data[:, _REF_COL['cadran2']],
    }
    print(f"[reference] Loaded reference data ← {filepath}  ({data.shape[0]} steps)")
    return ref

def plot_force(result):
    """
        Large 4×3 figure. If `ref` is provided, Robotran reference curves are
        overlaid in dashed red on every cell.
        """
    t = result['t']
    bodies = result['bodies']

    fig = plt.figure(figsize=(8, 5))

    # Retrieve simulation signals
    hinge_idx = _find_q_index(bodies, 'upper_pend', dof=0)
    theta_hinge = np.degrees(result['q'][:, hinge_idx])

    lower_body = next(b for b in bodies if b.name == 'lower_pend')
    acc1 = np.degrees(result['qdd'][:, lower_body.q_indices[0]])
    acc2 = np.degrees(result['qdd'][:, lower_body.q_indices[1]])
    acc3 = np.degrees(result['qdd'][:, lower_body.q_indices[2]])
    plt.plot(t, acc1, label='acc1')
    plt.plot(t, acc2, label='acc2')
    plt.plot(t, acc3, label='acc3')

    fig.tight_layout()
    _save(fig, "force_as_acc.pdf")
    return fig