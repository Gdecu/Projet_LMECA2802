"""
neri/state.py
-------------
MBState: the time-varying working state of one NERi evaluation.

Allocated once at the start of simulation and overwritten at every
timestep. Nothing here is permanent — it reflects the system state
at a single instant (q, qd, t).

Shape convention: axis 0 = body index (0 = base/ground, 1..N = bodies)
All vectors expressed in the INERTIAL frame unless noted otherwise.
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


@dataclass
class MBState:
    """
    Working state of the multibody system at one instant in time.

    All arrays are pre-allocated with zeros and overwritten by the
    forward/backward passes. Index 0 = base (ground), always zero/identity.

    Parameters
    ----------
    n_bodies : int
        Total number of bodies INCLUDING the base (so n_bodies = N+1
        where N is the number of moving bodies).
    n_dof : int
        Total number of generalised coordinates (size of q vector).
    """

    n_bodies: int
    n_dof:    int

    # ------------------------------------------------------------------ #
    # Current generalised state                                            #
    # ------------------------------------------------------------------ #

    t:   float = 0.0                          # current time [s]

    # Allocated in __post_init__ because shape depends on n_dof
    q:   np.ndarray = field(default=None)     # (n_dof,)  generalised positions
    qd:  np.ndarray = field(default=None)     # (n_dof,)  generalised velocities
    qdd: np.ndarray = field(default=None)     # (n_dof,)  generalised accelerations

    # ------------------------------------------------------------------ #
    # Forward pass quantities  (per body, inertial frame)                 #
    # ------------------------------------------------------------------ #

    # Rotation matrix: R[i] = R^{0i}, maps body-i frame → inertial frame
    R:    np.ndarray = field(default=None)    # (n_bodies, 3, 3)

    # Angular velocity and acceleration of each body
    omega:  np.ndarray = field(default=None)  # (n_bodies, 3)  ω^i
    omegad: np.ndarray = field(default=None)  # (n_bodies, 3)  ω̇^i

    # β^i = ω̃̇^i + ω̃^i·ω̃^i  (centripetal + angular acceleration matrix)
    beta:  np.ndarray = field(default=None)   # (n_bodies, 3, 3)

    # α^i = robotician's acceleration (absorbs gravity, Coriolis, centripetal)
    alpha: np.ndarray = field(default=None)   # (n_bodies, 3)

    # Geometry vectors rotated to inertial frame (recomputed each step)
    d_hi: np.ndarray = field(default=None)    # (n_bodies, 3)  O^h → O^i
    d_ii: np.ndarray = field(default=None)    # (n_bodies, 3)  O^i → G^i

    # Joint axis and linear velocity jacobian in inertial frame
    phi: np.ndarray = field(default=None)     # (n_bodies, 3)  φ^i  (or (n_bodies,3,k) for multi-dof)
    psi: np.ndarray = field(default=None)     # (n_bodies, 3)  ψ^i = φ^i × d^{ii}

    # ------------------------------------------------------------------ #
    # Backward pass quantities  (per body, inertial frame)                #
    # ------------------------------------------------------------------ #

    W: np.ndarray = field(default=None)       # (n_bodies, 3)  W^i  effective inertia force
    F: np.ndarray = field(default=None)       # (n_bodies, 3)  F^i  joint reaction force
    L: np.ndarray = field(default=None)       # (n_bodies, 3)  L^i  joint reaction moment

    # Projected generalised force Q^i = F^i·ψ^i + L^i·φ^i
    Q: np.ndarray = field(default=None)       # (n_dof,)

    # ------------------------------------------------------------------ #
    # Post-init: allocate all arrays                                       #
    # ------------------------------------------------------------------ #

    def __post_init__(self):
        N = self.n_bodies

        # Generalised coordinates
        self.q   = np.zeros(self.n_dof)
        self.qd  = np.zeros(self.n_dof)
        self.qdd = np.zeros(self.n_dof)

        # Forward pass — body 0 (base) stays at identity/zero throughout
        self.R      = np.zeros((N, 3, 3))
        self.R[0]   = np.eye(3)              # base frame = inertial frame

        self.omega  = np.zeros((N, 3))       # ω^0 = 0 (fixed base)
        self.omegad = np.zeros((N, 3))       # ω̇^0 = 0
        self.beta   = np.zeros((N, 3, 3))    # β^0 = 0
        self.alpha  = np.zeros((N, 3))       # α^0 = -g set by forward pass init

        self.d_hi   = np.zeros((N, 3))
        self.d_ii   = np.zeros((N, 3))
        self.phi    = np.zeros((N, 3))
        self.psi    = np.zeros((N, 3))

        # Backward pass
        self.W = np.zeros((N, 3))
        self.F = np.zeros((N, 3))
        self.L = np.zeros((N, 3))
        self.Q = np.zeros(self.n_dof)

    # ------------------------------------------------------------------ #
    # Convenience methods                                                  #
    # ------------------------------------------------------------------ #

    def set_q(self, q: np.ndarray, qd: np.ndarray, t: float):
        """Load a new (q, qd, t) into state and zero all computed quantities."""
        self.t  = t
        self.q  = q.copy()
        self.qd = qd.copy()

        # Reset computed arrays — forces forward/backward to recompute cleanly
        self.omega[:]  = 0.0
        self.omegad[:] = 0.0
        self.beta[:]   = 0.0
        self.alpha[:]  = 0.0
        self.W[:]      = 0.0
        self.F[:]      = 0.0
        self.L[:]      = 0.0
        self.Q[:]      = 0.0

        # Base rotation always stays identity
        self.R[0] = np.eye(3)

    def get_body_state(self, body_id: int) -> dict:
        """Return a readable snapshot of one body's state — useful for debugging."""
        return {
            'R':      self.R[body_id],
            'omega':  self.omega[body_id],
            'omegad': self.omegad[body_id],
            'alpha':  self.alpha[body_id],
            'F':      self.F[body_id],
            'L':      self.L[body_id],
        }

    def __repr__(self) -> str:
        return (f"MBState(t={self.t:.4f}s, n_bodies={self.n_bodies}, "
                f"n_dof={self.n_dof})")