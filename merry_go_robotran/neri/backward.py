"""
neri/backward.py
----------------
Phases 2 and 3 of NERi: backward dynamic recursion + projection.

For each body i (leaves → root):

    W[i]  = W^i = m^i·(α^i + β^i·d^{ii}) − F^i_ext          (4.73)
    F[i]  = F^i = Σ_{j∈children} F^j + W^i                   (4.74)
    L[i]  = L^i = Σ_{j∈children}(L^j + d̃^{ij}·F^j)
                 + d̃^{ii}·W^i − L^i_ext
                 + I^i·ω̇^i + ω̃^i·I^i·ω^i                    (4.75)

Projection onto each joint DOF k:

    Q[q_idx_k] = L^i · φ_k

    Note: the formula Q = F·ψ + L·φ found in the NERi slides uses L
    computed at the CoM G^i. Our L is at the joint origin O^i, which
    gives the equivalent (and simpler) form Q = L·φ.

Entry point
-----------
    backward_pass(bodies, joints, state, F_ext, L_ext)
"""

import numpy as np


def skew(v: np.ndarray) -> np.ndarray:
    """3×3 skew-symmetric: skew(v) @ u ≡ v × u."""
    return np.array([
        [ 0.0,  -v[2],  v[1]],
        [ v[2],  0.0,  -v[0]],
        [-v[1],  v[0],  0.0 ]
    ])


def backward_pass(bodies:  list,
                  joints:  list,
                  state,
                  F_ext:   np.ndarray = None,
                  L_ext:   np.ndarray = None) -> None:
    """
    Backward dynamic recursion — fills state.W, F, L, Q in-place.
    Must be called after forward_pass.

    Parameters
    ----------
    bodies : list of Body, TOPOLOGICAL order (reversed internally)
    joints : list of Joint, same order as bodies
    state  : MBState filled by forward_pass
    F_ext  : (n_bodies, 3) external forces at each CoM [N], or None
    L_ext  : (n_bodies, 3) external torques on each body [N·m], or None
    """

    n = state.R.shape[0]
    if F_ext is None: F_ext = np.zeros((n, 3))
    if L_ext is None: L_ext = np.zeros((n, 3))

    state.W[:] = 0.0
    state.F[:] = 0.0
    state.L[:] = 0.0
    state.Q[:] = 0.0

    # Build children lookup once
    children = {b.body_id: [] for b in bodies}
    children[0] = []
    for b in bodies:
        children[b.parent_id].append(b.body_id)

    # Traverse leaves → root
    for body, joint in zip(reversed(bodies), reversed(joints)):
        i = body.body_id
        m = body.mass
        I_loc = body.inertia_com_local

        # Rotate inertia to inertial frame: I_inertial = R · I_local · Rᵀ
        R0i      = state.R[i]
        I_world  = R0i @ I_loc @ R0i.T

        Om_i = skew(state.omega[i])

        # --- W^i = m^i·(α^i + β^i·d^{ii}) − F^i_ext  (4.73) ---
        # α^i is at O^i (pure transport); β^i·d^{ii} transports to G^i.
        # Together they give m·(ẍ_G − g), the inertia demand on the parent.
        state.W[i] = (m * (state.alpha[i] + state.beta[i] @ state.d_ii[i])
                      - F_ext[i])

        # --- F^i = Σ F^j + W^i  (4.74) ---
        state.F[i] = state.W[i].copy()
        for j in children[i]:
            state.F[i] += state.F[j]

        # --- L^i  (4.75) ---
        # Euler rotational inertia term at G^i (Euler's equation):
        euler = I_world @ state.omegad[i] + Om_i @ (I_world @ state.omega[i])

        # Moment from this body's own inertia force (lever arm from O^i to G^i)
        state.L[i] = skew(state.d_ii[i]) @ state.W[i] + euler - L_ext[i]

        # Add transported contributions from each child j
        # d_hi[j] = vector from O^i to O^j (already in inertial frame)
        for j in children[i]:
            state.L[i] += state.L[j] + skew(state.d_hi[j]) @ state.F[j]

        # --- Phase 3: project onto joint DOFs ---
        # Q_k = L^i · φ_k   (virtual work through angular displacement only;
        #                     O^i is fixed during virtual δq_k)
        if joint.n_dof == 0:
            continue

        q_i   = state.q[body.q_indices]
        phi_i = joint.axes_inertial(q_i, state.R[body.parent_id])  # (3, n_dof)

        for k, q_idx in enumerate(body.q_indices):
            state.Q[q_idx] = np.dot(state.L[i], phi_i[:, k])