"""
neri/assembly.py
----------------
Assembly of the mass matrix M(q) and bias vector c(q, qd) from the NERi
formalism, using the "last trick" described in NERi_V1.pdf slide 11.

────────────────────────────────────────────────────────────
  LAST TRICK (NERi slide 11)
  ──────────────────────────
  The full equation of motion is:  M(q) * qdd + c(q, qd) = Q_ext

  To extract M column-by-column:
    → run forward + backward with  qd = 0,  F_ext = 0,  L_ext = 0
      and  qdd = e_k  (k-th canonical unit vector)
    → the resulting Q vector is the k-th column of M

  To extract c:
    → run forward + backward with  qdd = 0  (qd and forces as usual)
    → the resulting Q vector equals c(q, qd, F_ext, L_ext)

  Both passes work on a TEMPORARY state so that the live simulation
  state is never corrupted.

────────────────────────────────────────────────────────────
  PARTITIONING (independent 'u' vs driven 'c' coordinates)
  ─────────────────────────────────────────────────────────
  Let  n  = total DOF,   nu = #independent,   nc = #driven.
  The global q-vector mixes 'u' and 'c' entries.  We partition:

      M = | Muu  Muc |      c = | cu |
          | Mcu  Mcc |          | cc |

  The integrator only needs the independent sub-system:
      Muu * qdd_u = Q_u - Muc * qdd_c - cu

  The driven accelerations are prescribed externally (driven_vars.py).

Public API
──────────
    assemble_Mc(bodies, joints, state, F_ext, L_ext)
        → M  (n×n),  c  (n,)

    partition_system(M, c, q_idx_u, q_idx_c, qdd_c)
        → A  (nu×nu),  rhs  (nu,)   ready for np.linalg.solve
"""

import numpy as np
import copy

from neri.forward  import forward_pass
from neri.backward import backward_pass
from neri.state    import MBState


# ═══════════════════════════════════════════════════════════════════════════ #
#  Internal helper                                                            #
# ═══════════════════════════════════════════════════════════════════════════ #

def _make_tmp_state(state: MBState) -> MBState:
    """
    Allocate a fresh MBState of the same size as `state`.
    Only q, qd, t are copied — all kinematic/dynamic arrays start at zero.
    This avoids polluting the live simulation state during column sweeps.
    """
    tmp = MBState(n_bodies=state.n_bodies, n_dof=state.n_dof)
    tmp.set_q(state.q.copy(), state.qd.copy(), state.t)
    return tmp


# ═══════════════════════════════════════════════════════════════════════════ #
#  Core assembly routine                                                      #
# ═══════════════════════════════════════════════════════════════════════════ #

def assemble_Mc(bodies:  list,
                joints:  list,
                state:   MBState,
                F_ext:   np.ndarray = None,
                L_ext:   np.ndarray = None,
                g:       np.ndarray = np.array([0.0, 0.0, -9.81])
                ) -> tuple:
    """
    Build the full mass matrix M(q) and bias vector c(q, qd).

    Parameters
    ----------
    bodies  : list of Body, topological order (parent before child)
    joints  : list of Joint, same order
    state   : MBState with q and qd already loaded for this timestep
    F_ext   : (n_bodies, 3) external forces at each CoM [N], or None
    L_ext   : (n_bodies, 3) external torques on each body [N·m], or None
    g       : gravity vector [m/s²]

    Returns
    -------
    M : np.ndarray, shape (n_dof, n_dof)
        Symmetric positive-definite mass (inertia) matrix.
    c : np.ndarray, shape (n_dof,)
        Bias vector (Coriolis + centripetal + gravity + external loads).
    """

    n_dof = state.n_dof
    n_bodies = state.n_bodies

    # Default external loads to zero if not provided
    if F_ext is None:
        F_ext = np.zeros((n_bodies, 3))
    if L_ext is None:
        L_ext = np.zeros((n_bodies, 3))

    # ── 1. Compute c: run with qdd = 0, keep qd and F_ext ──────────────── #
    tmp_c = _make_tmp_state(state)
    tmp_c.qdd[:] = 0.0                          # qdd = 0 → only bias terms
    forward_pass(bodies, joints, tmp_c, g)
    backward_pass(bodies, joints, tmp_c, F_ext, L_ext)
    c_vec = tmp_c.Q.copy()                       # c(q, qd, F_ext, L_ext)

    # ── 2. Compute M column by column (last trick) ───────────────────────── #
    # For each DOF k: set qdd = e_k, zero velocities and forces,
    # run both passes, read Q → that is M[:, k].
    M_mat = np.zeros((n_dof, n_dof))

    for k in range(n_dof):
        tmp_k = _make_tmp_state(state)
        tmp_k.qd[:]  = 0.0                       # zero velocities (no Coriolis)
        tmp_k.qdd[:] = 0.0
        tmp_k.qdd[k] = 1.0                       # unit impulse on DOF k

        # No external forces for mass column sweep
        forward_pass(bodies, joints, tmp_k, g)
        backward_pass(bodies, joints, tmp_k,
                      F_ext=np.zeros((n_bodies, 3)),
                      L_ext=np.zeros((n_bodies, 3)))

        M_mat[:, k] = tmp_k.Q.copy()            # k-th column of M

    return M_mat, c_vec


# ═══════════════════════════════════════════════════════════════════════════ #
#  System partitioning                                                        #
# ═══════════════════════════════════════════════════════════════════════════ #

def partition_system(M:       np.ndarray,
                     c:       np.ndarray,
                     idx_u:   list,
                     idx_c:   list,
                     qdd_c:   np.ndarray,
                     Q_ext:   np.ndarray = None
                     ) -> tuple:
    """
    Partition the full system  M*qdd + c = Q_ext  into:
      - the independent sub-system (to be solved for qdd_u)
      - the driven sub-system      (for recovering internal forces)

    The independent equations read:
        Muu * qdd_u = (Q_u - cu) - Muc * qdd_c
                    =  rhs

    Parameters
    ----------
    M      : (n, n) full mass matrix
    c      : (n,)   full bias vector
    idx_u  : list of int — global indices of independent ('u') DOFs
    idx_c  : list of int — global indices of driven ('c') DOFs
    qdd_c  : (nc,) prescribed driven accelerations
    Q_ext  : (n,) generalised external forces, or None (→ zero)

    Returns
    -------
    A   : (nu, nu) — Muu, the independent mass sub-matrix
    rhs : (nu,)    — right-hand side for  A * qdd_u = rhs
    """

    n = M.shape[0]
    if Q_ext is None:
        Q_ext = np.zeros(n)

    # Extract sub-blocks
    Muu = M[np.ix_(idx_u, idx_u)]   # (nu × nu) — independent inertia
    Muc = M[np.ix_(idx_u, idx_c)]   # (nu × nc) — coupling

    c_u = c[idx_u]                   # bias on independent equations

    # rhs = Q_u - c_u - Muc * qdd_c
    rhs = Q_ext[idx_u] - c_u - Muc @ qdd_c

    return Muu, rhs


def solve_independent_accelerations(M:     np.ndarray,
                                     c:     np.ndarray,
                                     idx_u: list,
                                     idx_c: list,
                                     qdd_c: np.ndarray,
                                     Q_ext: np.ndarray = None
                                     ) -> np.ndarray:
    """
    Convenience wrapper: partition and solve for qdd_u in one call.

    Returns
    -------
    qdd_u : (nu,) independent generalised accelerations
    """

    Muu, rhs = partition_system(M, c, idx_u, idx_c, qdd_c, Q_ext)

    # Solve Muu * qdd_u = rhs  (Muu is symmetric positive-definite)
    qdd_u = np.linalg.solve(Muu, rhs)

    return qdd_u


# ═══════════════════════════════════════════════════════════════════════════ #
#  Index extraction helpers                                                   #
# ═══════════════════════════════════════════════════════════════════════════ #

def get_dof_indices(bodies: list) -> tuple:
    """
    Scan the body list and return sorted index lists for independent ('u')
    and driven ('c') DOFs.

    Returns
    -------
    idx_u : list of int — global q-indices of independent DOFs
    idx_c : list of int — global q-indices of driven DOFs
    """

    idx_u = []
    idx_c = []

    for body in bodies:
        for var_type, q_idx in zip(body.var_types, body.q_indices):
            if var_type == 'u':
                idx_u.append(q_idx)
            elif var_type == 'c':
                idx_c.append(q_idx)

    return sorted(idx_u), sorted(idx_c)