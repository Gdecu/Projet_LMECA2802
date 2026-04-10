"""
simulation/driven_vars.py
-------------------------
Prescribed motion laws for all driven ('c') generalised coordinates.

Every driven DOF has an associated law that returns (qc, qdc, qddc) at
time t. The integrator calls `evaluate_driven(bodies, t)` once per
timestep to inject these values into the global q/qd/qdd vectors before
the NERi passes run.

────────────────────────────────────────────────────────────
  CAROUSEL DRIVEN DOFs  (from MOOC document & carousel_data.py)
  ─────────────────────────────────────────────────────────
  1. pole_tilt1  (q-index assigned in build_bodies)
       θ1(t) = A * (1 − cos(ω*t + φ1))
       φ1 = π/2,  A = 2.5 * 2π/360 rad,  ω = 0.4π rad/s

  2. pole_tilt2  (q-index assigned in build_bodies)
       θ2(t) = A * (1 − cos(ω*t + φ2))
       φ2 = 0

  3. nacelle_k  (2 DOF Cardan, 4 nacelles = 8 DOFs)
       Both Cardan angles are locked to zero — the nacelle is kept
       aligned with the pendulum (gravity keeps it vertical in practice;
       the constraint is that the gondola cannot spin about the pendulum
       axis, enforced by the Cardan joint topology).

────────────────────────────────────────────────────────────
  Generic design
  ──────────────
  Each law is a callable registered in DRIVEN_LAW_REGISTRY under the
  body name (or a pattern). Adding a new driven DOF only requires adding
  one entry to the registry — the integrator loop is unchanged.

  A law callable has signature:
      law(t: float) -> (qc: float, qdc: float, qddc: float)
"""

import numpy as np
from math import cos, sin, pi


# ═══════════════════════════════════════════════════════════════════════════ #
#  Physical constants (matching carousel_data.py)                            #
# ═══════════════════════════════════════════════════════════════════════════ #

_TILT_A     = 2.5 * (2.0 * pi / 360.0)   # [rad]  amplitude ≈ 2.5°
_TILT_OMEGA = 0.4 * pi                    # [rad/s] frequency
_TILT_PHI1  = pi / 2.0                   # [rad]  phase for tilt-1 (I1 axis)
_TILT_PHI2  = 0.0                        # [rad]  phase for tilt-2 (I2 axis)


# ═══════════════════════════════════════════════════════════════════════════ #
#  Individual law functions                                                   #
# ═══════════════════════════════════════════════════════════════════════════ #

def _law_tilt1(t: float) -> tuple:
    """
    Pole tilt about I1 — 1st driven DOF.
    θ = A*(1 - cos(ω*t + φ1))
    Returns (position [rad], velocity [rad/s], acceleration [rad/s²]).
    """
    phase   = _TILT_OMEGA * t + _TILT_PHI1
    qc      =  _TILT_A * (1.0 - cos(phase))
    qdc     =  _TILT_A * _TILT_OMEGA * sin(phase)
    qddc    =  _TILT_A * _TILT_OMEGA**2 * cos(phase)
    return qc, qdc, qddc


def _law_tilt2(t: float) -> tuple:
    """
    Pole tilt about I2 — 2nd driven DOF.
    θ = A*(1 - cos(ω*t + φ2))
    Returns (position [rad], velocity [rad/s], acceleration [rad/s²]).
    """
    phase   = _TILT_OMEGA * t + _TILT_PHI2
    qc      =  _TILT_A * (1.0 - cos(phase))
    qdc     =  _TILT_A * _TILT_OMEGA * sin(phase)
    qddc    =  _TILT_A * _TILT_OMEGA**2 * cos(phase)
    return qc, qdc, qddc


def _law_locked(t: float) -> tuple:
    """
    Zero-motion law for any Cardan DOF that is kept at rest.
    Used for nacelle joints: both Cardan angles remain 0 at all times.
    Returns (0, 0, 0).
    """
    return 0.0, 0.0, 0.0


# ═══════════════════════════════════════════════════════════════════════════ #
#  Registry: body name → list of law callables (one per driven DOF)          #
# ═══════════════════════════════════════════════════════════════════════════ #
#
# Keys must match body.name as defined in carousel_data.py / build_bodies().
# Multi-DOF bodies (Cardan2) get a list with one entry per DOF, in q-order.
#
DRIVEN_LAW_REGISTRY: dict = {
    # Pole tilt — each is a single 1-DOF revolute body
    "pole_tilt1":  [_law_tilt1],
    "pole_tilt2":  [_law_tilt2],

    # Nacelle Cardan joints — 2 DOF each, both locked to zero
    "nacelle_1":   [_law_locked, _law_locked],
    "nacelle_2":   [_law_locked, _law_locked],
    "nacelle_3":   [_law_locked, _law_locked],
    "nacelle_4":   [_law_locked, _law_locked],
}


# ═══════════════════════════════════════════════════════════════════════════ #
#  Main evaluation routine                                                    #
# ═══════════════════════════════════════════════════════════════════════════ #

def evaluate_driven(bodies: list, t: float) -> dict:
    """
    Evaluate all driven laws at time t and return a mapping from
    global q-index to (qc, qdc, qddc).

    Parameters
    ----------
    bodies : list of Body, same order as used in main.py
    t      : current simulation time [s]

    Returns
    -------
    driven_values : dict  { q_index: (qc, qdc, qddc) }
        Contains one entry per driven DOF.  Only 'c' entries are included.
    """
    driven_values = {}

    for body in bodies:
        # Skip bodies that have no driven DOFs
        if 'c' not in body.var_types:
            continue

        laws = DRIVEN_LAW_REGISTRY.get(body.name)
        if laws is None:
            raise KeyError(
                f"Body '{body.name}' has driven DOF(s) but no entry in "
                f"DRIVEN_LAW_REGISTRY. Add a law for this body."
            )

        # Iterate over each DOF of this body
        for var_type, q_idx, law in zip(body.var_types, body.q_indices, laws):
            if var_type == 'c':
                driven_values[q_idx] = law(t)   # (qc, qdc, qddc)

    return driven_values


def inject_driven_into_state(state, driven_values: dict) -> None:
    """
    Write the driven (qc, qdc, qddc) values from `driven_values` into
    the global state vectors (q, qd, qdd) at the correct indices.

    Call this at the START of every yd=f(y,t) evaluation, BEFORE any
    NERi forward/backward pass.

    Parameters
    ----------
    state         : MBState — modified in-place
    driven_values : dict returned by evaluate_driven()
    """
    for q_idx, (qc, qdc, qddc) in driven_values.items():
        state.q[q_idx]   = qc     # overwrite position with prescribed value
        state.qd[q_idx]  = qdc    # overwrite velocity
        state.qdd[q_idx] = qddc   # overwrite acceleration (used in forward pass)


def get_qddc_vector(bodies: list, t: float, n_dof: int) -> np.ndarray:
    """
    Build the full driven-acceleration vector qdd_c in global q-index order.

    This is consumed by assembly.partition_system() as the `qdd_c` argument.

    Parameters
    ----------
    bodies : list of Body
    t      : current time [s]
    n_dof  : total number of DOFs (size of global q)

    Returns
    -------
    qdd_c_full : (n_dof,) array — zeros at independent indices,
                 prescribed qddc at driven indices
    """
    driven_values = evaluate_driven(bodies, t)

    qdd_c_full = np.zeros(n_dof)
    for q_idx, (_qc, _qdc, qddc) in driven_values.items():
        qdd_c_full[q_idx] = qddc

    return qdd_c_full