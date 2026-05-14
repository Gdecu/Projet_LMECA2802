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

# ═══════════════════════════════════════════════════════════════════════════ #
#  Plot D — Peak reaction force vs. split position along arm                  #
# ═══════════════════════════════════════════════════════════════════════════ #

def plot_force_vs_length(max_force_at: dict) -> plt.Figure:
    """
    Two-panel figure summarising the dimensioning sweep.

    Top panel   — peak generalised-force magnitude |F| [N or N·m] at each
                  split position, i.e. the envelope used for sizing.
    Bottom panel — the three components (fx, fy, fz) at the instant when
                  |F| was maximum, so the dominant loading axis is visible.

    Parameters
    ----------
    max_force_at : dict with keys
        'length'    : list[float]  split positions tested [m]
        'max_force' : list[float]  peak |F| at each position
        'fx','fy','fz' : list[float] components at peak instant
        'time'      : list[float]  time at which the peak occurred [s]
    """
    length    = np.array(max_force_at['length'])
    max_force = np.array(max_force_at['max_force'])
    fx        = np.array(max_force_at['fx'])
    fy        = np.array(max_force_at['fy'])
    fz        = np.array(max_force_at['fz'])

    # Index of the worst-case section
    i_crit = int(np.argmax(max_force))

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Dimensioning — Peak Reaction Force vs. Split Position",
                 fontweight='bold')

    # — Top: magnitude envelope ——————————————————————————————————————————
    ax0 = axes[0]
    ax0.plot(length, max_force, color='steelblue', linewidth=1.8, label='|F| peak')
    ax0.axvline(length[i_crit], color='crimson', linewidth=1.0, linestyle='--',
                label=f'Worst section  ({length[i_crit]:.3f} m)')
    ax0.set_ylabel("Peak Force Magnitude  [N·m]")
    ax0.legend()
    ax0.grid(True, alpha=0.35)

    # — Bottom: force components ————————————————————————————————————————
    ax1 = axes[1]
    ax1.plot(length, fx, color='crimson',    linewidth=1.4, label='F₁')
    ax1.plot(length, fy, color='seagreen',   linewidth=1.4, label='F₂')
    ax1.plot(length, fz, color='darkorange', linewidth=1.4, label='F₃')
    ax1.axhline(0, color='black', linewidth=0.6, linestyle=':')
    ax1.set_xlabel("Split Position along Arm  [m]")
    ax1.set_ylabel("Force Components  [N·m]")
    ax1.legend()
    ax1.grid(True, alpha=0.35)

    fig.tight_layout()
    _save(fig, "force_vs_length.pdf")
    return fig


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


# ═══════════════════════════════════════════════════════════════════════════ #
#  Plot E — Internal reaction forces vs. time at one cross-section            #
# ═══════════════════════════════════════════════════════════════════════════ #

def plot_force_vs_time(t: np.ndarray,
                       Q_react_c: np.ndarray,
                       len_section: float) -> plt.Figure:
    """
    Two-panel time-domain plot for the generalised reactions at a single cross-section.

    Panel 1 — |Q| magnitude with peak marker
    Panel 2 — three cardan-axis components Q1, Q2, Q3

    Parameters
    ----------
    t            : (n_steps,)    time vector [s]
    Q_react_c    : (n_steps, 3)  generalised cardan reactions [N·m]
    len_section  : float         position of the cut along the arm [m]
    """
    mag = np.linalg.norm(Q_react_c, axis=1)
    i_peak = int(np.argmax(mag))

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    fig.suptitle(f"Internal Loads vs. Time — Cut at {len_section:.3f} m along arm",
                 fontweight='bold')

    # — Panel 1: magnitude ————————————————————————————————————————————————
    ax = axes[0]
    ax.plot(t, mag, color='steelblue', linewidth=1.6)
    ax.axvline(t[i_peak], color='crimson', linewidth=1.0, linestyle='--',
               label=f'Peak  {mag[i_peak]:.1f} N·m  @ t={t[i_peak]:.2f} s')
    ax.set_ylabel("|Q|  [N·m]")
    ax.set_title("Generalised Reaction Magnitude")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.35)

    # — Panel 2: cardan components ————————————————————————————————————————
    ax = axes[1]
    for j, (color, label) in enumerate(zip(
            ['crimson', 'seagreen', 'darkorange'], ['Q₁', 'Q₂', 'Q₃'])):
        ax.plot(t, Q_react_c[:, j], color=color, linewidth=1.4, label=label)
    ax.axhline(0, color='black', linewidth=0.5, linestyle=':')
    ax.set_xlabel("Time  [s]")
    ax.set_ylabel("Q  [N·m]")
    ax.set_title("Generalised Cardan-Axis Reactions")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.35)

    fig.tight_layout()
    _save(fig, f"force_vs_time_{len_section:.3f}m.pdf")
    return fig