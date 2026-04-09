"""
neri/forward.py
---------------
Phase 1 of the NERi formalism: forward kinematic recursion.

For each body i (root → leaves), computes:

    R[i]      = R^{0i}    rotation matrix (body frame → inertial)
    omega[i]  = ω^i       absolute angular velocity          (3.58)
    omegad[i] = ω̇^i       absolute angular acceleration      (3.59)
    beta[i]   = β^i       β = skew(ω̇) + skew(ω)²             (3.60)
    alpha[i]  = α^i       acceleration of joint origin O^i   (*)
    d_hi[i]   = d^{hi}    arm vector O^h→O^i in inertial frame
    d_ii[i]   = d^{ii}    CoM offset O^i→G^i in inertial frame
    phi[i]    = φ^i       joint axis in inertial frame
    psi[i]    = ψ^i       linear velocity jacobian φ^i × d^{ii}

(*) α^i is the "acceleration at the joint origin O^i minus g":
        α^i = p̈_{O^i} − g
    propagated via:   α^j = α^i + β^i · d^{hi}
    WITHOUT Coriolis (2ω̃·ψ·q̇) and WITHOUT direct (ψ·q̈) terms.
    Those effects enter through ω̇^i → β^i → W^i in the backward pass.
    α^0 = −g  (base: fixed, p̈=0).

Entry point
-----------
    forward_pass(bodies, joints, state, g)
"""

import numpy as np


def skew(v: np.ndarray) -> np.ndarray:
    """3×3 skew-symmetric matrix: skew(v) @ u ≡ v × u."""
    return np.array([
        [ 0.0,  -v[2],  v[1]],
        [ v[2],  0.0,  -v[0]],
        [-v[1],  v[0],  0.0 ]
    ])


def forward_pass(bodies: list,
                 joints: list,
                 state,
                 g: np.ndarray = np.array([0.0, 0.0, -9.81])) -> None:
    """
    Forward kinematic recursion — fills state in-place.

    Parameters
    ----------
    bodies : list of Body, topological order (parent before child)
    joints : list of Joint, same order as bodies
    state  : MBState with q, qd, qdd already set
    g      : gravity vector [m/s²], default [0, 0, -9.81]
    """

    # --- Base (body 0): fixed, no velocity or acceleration ---
    # α^0 = -g encodes that the base "sees" gravity; this propagates
    # down the tree so every body's W will automatically include gravity.
    state.R[0]      = np.eye(3)
    state.omega[0]  = np.zeros(3)
    state.omegad[0] = np.zeros(3)
    state.beta[0]   = np.zeros((3, 3))
    state.alpha[0]  = -g

    for body, joint in zip(bodies, joints):
        i = body.body_id
        h = body.parent_id

        q_i   = state.q[body.q_indices]
        qd_i  = state.qd[body.q_indices]
        qdd_i = state.qdd[body.q_indices]

        # --- A: rotation matrix R^{0i} = R^{0h} · R^{hi}(q) ---
        state.R[i] = state.R[h] @ joint.rotation_matrix(q_i)

        # --- B: geometry vectors in inertial frame ---
        state.d_hi[i] = state.R[h] @ body.d_parent_to_joint_in_parent
        state.d_ii[i] = state.R[i] @ body.d_joint_to_com_in_local

        # --- C: joint axes and ψ = φ × d^{ii} in inertial frame ---
        phi_i = joint.axes_inertial(q_i, state.R[h])  # (3, n_dof)
        psi_i = joint.psi(phi_i, state.d_ii[i])        # (3, n_dof)

        if joint.n_dof == 1:
            state.phi[i] = phi_i[:, 0]
            state.psi[i] = psi_i[:, 0]

        # --- D: angular velocity  ω^i = ω^h + Σ_k φ_k q̇_k  (3.58) ---
        state.omega[i] = state.omega[h].copy()
        for k in range(joint.n_dof):
            state.omega[i] += phi_i[:, k] * qd_i[k]

        # --- E: angular acceleration  (3.59) ---
        # ω̇^i = ω̇^h + ω̃^i·φ_k·q̇_k  (Coriolis on axis)
        #              + φ_k·q̈_k       (direct)
        # NOTE: ω̇ includes q̈, so β (used in W and child propagation)
        # automatically captures all centripetal and Coriolis effects.
        Om_i = skew(state.omega[i])
        state.omegad[i] = state.omegad[h].copy()
        for k in range(joint.n_dof):
            state.omegad[i] += Om_i @ phi_i[:, k] * qd_i[k]   # axis Coriolis
            state.omegad[i] += phi_i[:, k] * qdd_i[k]          # direct

        # --- F: β^i = skew(ω̇^i) + skew(ω^i)²  (3.60) ---
        state.beta[i] = skew(state.omegad[i]) + Om_i @ Om_i

        # --- G: α^i = acceleration of joint origin O^i minus g ---
        # α^i = α^h + β^h · d^{hi}
        # This is a PURE TRANSPORT: no Coriolis or q̈ terms.
        # Coriolis/centripetal effects enter through β^i in W (backward).
        state.alpha[i] = state.alpha[h] + state.beta[h] @ state.d_hi[i]