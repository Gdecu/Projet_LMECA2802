"""
simulation/integrator.py
------------------------
Defines the ODE right-hand side  yd = f(y, t)  and wraps scipy's
solve_ivp (RK45, which matches the RK4 spirit requested by the project).

State vector layout
───────────────────
  y = [q_u | qd_u]      shape: (2*nu,)
  where  q_u  are the nu INDEPENDENT generalised coordinates and
         qd_u  are their time derivatives.

  Driven coordinates (qc, qdc) are NOT part of y — they are
  re-injected at every function evaluation from driven_vars.py.

Pipeline inside f(y, t)
───────────────────────
  1. Unpack y  →  q_u, qd_u
  2. Evaluate driven laws  →  qc, qdc, qddc  (from driven_vars)
  3. Assemble full q, qd, qdd into state
  4. Compute external forces / torques  (wind, motor, spring-damper, damping)
  5. Assemble M(q) and c(q,qd)          (assembly.py  last trick)
  6. Partition and solve for qdd_u      (assembly.partition_system)
  7. Return yd = [qd_u | qdd_u]
"""

import numpy as np
from math import cos, pi
from scipy.integrate import solve_ivp

from merry_go_robotran.neri.state    import MBState
from merry_go_robotran.neri.forward  import forward_pass
from merry_go_robotran.neri.backward import backward_pass
from merry_go_robotran.neri.assembly import (assemble_Mc,
                            solve_independent_accelerations,
                            get_dof_indices)
from merry_go_robotran.simulation.driven_vars import (evaluate_driven,
                                     inject_driven_into_state,
                                     get_qddc_vector)
import merry_go_robotran.data.carousel_data as cd


# ═══════════════════════════════════════════════════════════════════════════ #
#  External force / torque builders                                           #
# ═══════════════════════════════════════════════════════════════════════════ #

def _motor_torque(t: float, omega_pole: float, t0_motor: float) -> float:
    """
    Three-phase motor torque law on the main pole rotation (I3 axis).

    Phase 1 — ω_pole has never reached threshold:  T = T_phase1
    Phase 2 — within 0.5 s after first crossing t0: T = A*(1+cos(ω*t + φ))
    Phase 3 — beyond 0.5 s after crossing:          T = 0

    Parameters
    ----------
    t           : current time [s]
    omega_pole  : current angular velocity of the pole about I3 [rad/s]
    t0_motor    : time of first threshold crossing (-1 if never reached)

    Returns
    -------
    torque [N·m]  (positive = accelerating spin direction)
    """
    if t0_motor < 0:
        # Phase 1: motor at full torque until threshold is crossed
        return cd.T_motor_phase1
    elif t < t0_motor + 0.5:
        # Phase 2: smooth shut-down ramp over 0.5 s
        phi_motor = -2.0 * pi * t0_motor
        return cd.T_motor_A * (1.0 + cos(cd.omega_ctrl_motor * t + phi_motor))
    else:
        # Phase 3: motor off
        return 0.0


def _wind_force(t: float, n_bodies: int, body_name_to_id: dict) -> np.ndarray:
    """
    Wind force F = A*(1-cos(ω*t)) applied along I1 to each nacelle CoM.

    Returns
    -------
    F_ext : (n_bodies, 3) array — zero everywhere except nacelle bodies
    """
    F_ext = np.zeros((n_bodies, 3))
    amplitude = cd.F_wind_A * (1.0 - cos(cd.omega_wind * t))   # [N]

    for name in ["nacelle_1", "nacelle_2", "nacelle_3", "nacelle_4"]:
        bid = body_name_to_id.get(name)
        if bid is not None:
            F_ext[bid, 0] = amplitude   # I1 direction = global x-axis
    return F_ext


def _damping_torques(bodies: list, state: MBState,
                     body_name_to_id: dict) -> np.ndarray:
    """
    Viscous damping generalised torques  Q_damp = -b * qd  for:
      - pendulum hinges   (per-pendulum damping coefficient)
      - nacelle Cardan joints (uniform 6 Ns/m)

    Returns
    -------
    Q_damp : (n_dof,) array — generalised damping forces (negative = dissipative)
    """
    Q_damp = np.zeros(state.n_dof)

    # Pendulum hinge damping
    for pname, b_coeff in cd.damping_pendulum_hinge.items():
        bid = body_name_to_id.get(pname)
        if bid is None:
            continue
        body = bodies[bid - 1]   # body_id 1-indexed; bodies list 0-indexed
        for var_type, q_idx in zip(body.var_types, body.q_indices):
            if var_type == 'u':
                Q_damp[q_idx] -= b_coeff * state.qd[q_idx]

    # Nacelle Cardan damping (both DOFs, driven but damping still acts)
    for nname in ["nacelle_1", "nacelle_2", "nacelle_3", "nacelle_4"]:
        bid = body_name_to_id.get(nname)
        if bid is None:
            continue
        body = bodies[bid - 1]
        for q_idx in body.q_indices:
            Q_damp[q_idx] -= cd.damping_nacelle_cardan * state.qd[q_idx]

    return Q_damp


def _spring_damper_forces(bodies: list, state: MBState,
                           body_name_to_id: dict) -> tuple:
    """
    Spring-damper force between each arm attachment point and the
    corresponding pendulum attachment point (4 elements).

    Attachment geometry (from carousel_data.py):
      - Arm side:     0.5 m from pole axis along arm direction
      - Pendulum side: 1.5 m below hinge, 0.1 m inward

    The force is aligned with the connector vector; it contributes to
    both F_ext (applied to pendulum CoM via virtual work) and, by
    reaction, to the arm/pole body.

    Returns
    -------
    F_ext  : (n_bodies, 3) spring-damper contributions to external forces
    L_ext  : (n_bodies, 3) spring-damper contributions to external torques
    """
    n = state.n_bodies
    F_ext = np.zeros((n, 3))
    L_ext = np.zeros((n, 3))

    arm_names     = ["arm_1",     "arm_2",     "arm_3",     "arm_4"]
    pendulum_names= ["pendulum_1","pendulum_2","pendulum_3","pendulum_4"]

    for arm_name, pend_name in zip(arm_names, pendulum_names):
        arm_bid  = body_name_to_id.get(arm_name)
        pend_bid = body_name_to_id.get(pend_name)
        if arm_bid is None or pend_bid is None:
            continue

        # Attachment point on arm in inertial frame
        # p_arm = O^arm + R^arm * [r_attach, 0, 0]
        p_arm_attach_local = np.array([cd.r_spring_arm_attach, 0.0, 0.0])
        p_arm = (state.alpha[arm_bid]           # reuse joint origin from forward pass?
                 # NOTE: we reconstruct from geometry for correctness
                 )
        # Reconstruct arm joint origin position in inertial frame
        # p_O^arm = p_O^pole + R^pole * d_arm_in_pole
        # For a clean computation we use the d_hi chain:
        # p_O^i = Σ_{k=0}^{i} d_hi[k]  (cumulative sum of joint vectors)
        # This is approximate for general trees; we use R[i] @ local_offset.
        p_O_arm   = _joint_origin_position(arm_bid,  bodies, state)
        p_O_pend  = _joint_origin_position(pend_bid, bodies, state)

        # Arm attachment in inertial frame
        p_A = p_O_arm + state.R[arm_bid] @ p_arm_attach_local

        # Pendulum attachment: 1.5 m below hinge, 0.1 m inward (local z-down, x-inward)
        p_pend_attach_local = np.array([-cd.r_spring_pendulum_attach_inward,
                                         0.0,
                                        -cd.r_spring_pendulum_attach_below])
        p_B = p_O_pend + state.R[pend_bid] @ p_pend_attach_local

        # Connector vector and length
        d_vec    = p_B - p_A                    # from arm attach to pendulum attach
        dist     = np.linalg.norm(d_vec)
        if dist < 1e-9:
            continue
        unit_vec = d_vec / dist

        # Relative velocity along connector (for damping)
        # Approximate: use CoM velocities — good enough for generalised forces
        # v_A ≈ 0 (arm is fixed to pole, slow); computed via ω × r
        v_B_approx = np.cross(state.omega[pend_bid],
                               p_B - (p_O_pend + state.d_ii[pend_bid]))
        v_rel = np.dot(v_B_approx, unit_vec)   # scalar relative velocity

        # Spring-damper force magnitude (F > 0 = tension)
        F_mag = cd.k_spring * (dist - cd.L0_spring) + cd.c_damper * v_rel

        # Force vector on pendulum (along connector, toward arm)
        F_on_pend = -F_mag * unit_vec

        # Apply as external force at pendulum CoM (via lever arm torque)
        r_com_from_B = state.d_ii[pend_bid] - (p_B - p_O_pend)
        F_ext[pend_bid] += F_on_pend
        L_ext[pend_bid] += np.cross(r_com_from_B, F_on_pend)

    return F_ext, L_ext


def _joint_origin_position(body_id: int, bodies: list, state: MBState) -> np.ndarray:
    """
    Reconstruct the inertial position of joint origin O^i by walking
    up the kinematic tree and summing d_hi vectors.

    This is only valid AFTER forward_pass has populated state.d_hi.
    """
    pos = np.zeros(3)
    bid = body_id
    # Walk from this body up to the base, accumulating joint-to-joint vectors
    visited = []
    while bid != 0:
        visited.append(bid)
        # Find parent
        body = next(b for b in bodies if b.body_id == bid)
        bid  = body.parent_id

    # Sum d_hi from root to this body
    for b in reversed(visited):
        pos += state.d_hi[b]

    return pos


# ═══════════════════════════════════════════════════════════════════════════ #
#  ODE right-hand side factory                                                #
# ═══════════════════════════════════════════════════════════════════════════ #

def make_ode(bodies: list, joints: list, state: MBState) -> callable:
    """
    Build and return the closure  yd = f(t, y)  compatible with solve_ivp.

    The closure captures `bodies`, `joints`, and `state` (pre-allocated).
    A mutable dict `ctx` carries persistent state across calls (motor t0).

    Parameters
    ----------
    bodies : list of Body
    joints : list of Joint
    state  : MBState, allocated once and reused every call

    Returns
    -------
    f_ode  : callable(t, y) → yd   (shape (2*nu,))
    ctx    : dict — runtime context (inspect t0_motor, etc.)
    """

    # Pre-compute index lists once — they never change during simulation
    idx_u, idx_c = get_dof_indices(bodies)
    nu = len(idx_u)

    # Body name → body_id lookup for force routines
    body_name_to_id = {b.name: b.body_id for b in bodies}

    # Mutable simulation context (persists across ODE evaluations)
    ctx = {"t0_motor": -1.0}    # -1 = threshold never yet crossed

    def f_ode(t: float, y: np.ndarray) -> np.ndarray:
        """
        ODE right-hand side: yd = [qd_u, qdd_u].

        Parameters
        ----------
        t : float — current time [s]
        y : (2*nu,) — [q_u | qd_u] for independent DOFs only

        Returns
        -------
        yd : (2*nu,) — [qd_u | qdd_u]
        """

        # ── 1. Unpack independent state ──────────────────────────────── #
        q_u  = y[:nu]
        qd_u = y[nu:]

        # ── 2. Evaluate driven laws at time t ────────────────────────── #
        driven_vals = evaluate_driven(bodies, t)

        # ── 3. Build full q, qd, qdd in state ───────────────────────── #
        state.set_q(np.zeros(state.n_dof), np.zeros(state.n_dof), t)
        state.q[idx_u]  = q_u
        state.qd[idx_u] = qd_u
        inject_driven_into_state(state, driven_vals)   # writes qc, qdc, qddc

        # ── 4. Check motor threshold for first time ───────────────────── #
        pole_body    = next(b for b in bodies if b.name == "pole")
        pole_q_idx   = pole_body.q_indices[0]
        omega_pole   = state.qd[pole_q_idx]            # pole angular velocity
        if ctx["t0_motor"] < 0 and omega_pole >= cd.omega_pole_threshold:
            ctx["t0_motor"] = t                         # record first crossing

        # ── 5. Build external force/torque arrays ────────────────────── #
        # We need one forward pass (without qdd) to get geometry for forces
        state.qdd[:] = 0.0
        forward_pass(bodies, joints, state)

        # Wind on nacelles
        F_ext, L_ext = _wind_force(t, state.n_bodies, body_name_to_id), \
                       np.zeros((state.n_bodies, 3))

        # Spring-damper between arms and pendulums
        F_sd, L_sd = _spring_damper_forces(bodies, state, body_name_to_id)
        F_ext += F_sd
        L_ext += L_sd

        # Motor torque on pole (projected onto I3 = pole spin axis)
        T_motor     = _motor_torque(t, omega_pole, ctx["t0_motor"])
        Q_ext       = np.zeros(state.n_dof)
        Q_ext[pole_q_idx] = T_motor

        # Viscous damping at pendulum and nacelle joints
        Q_ext += _damping_torques(bodies, state, body_name_to_id)

        # ── 6. Assemble M(q) and c(q,qd)  (last trick) ───────────────── #
        M, c_vec = assemble_Mc(bodies, joints, state, F_ext, L_ext)

        # ── 7. Solve for independent accelerations ───────────────────── #
        qdd_c = get_qddc_vector(bodies, t, state.n_dof)   # prescribed qddc

        qdd_u = solve_independent_accelerations(
            M, c_vec, idx_u, idx_c, qdd_c[idx_c], Q_ext
        )

        # ── 8. Return yd = [qd_u | qdd_u] ───────────────────────────── #
        yd = np.empty(2 * nu)
        yd[:nu] = qd_u
        yd[nu:] = qdd_u
        return yd

    return f_ode, ctx


# ═══════════════════════════════════════════════════════════════════════════ #
#  Integration runner                                                         #
# ═══════════════════════════════════════════════════════════════════════════ #

def run_simulation(bodies: list,
                   joints: list,
                   q0_full: np.ndarray,
                   qd0_full: np.ndarray,
                   t_end:    float = 30.0,
                   rtol:     float = 1e-6,
                   atol:     float = 1e-9
                   ) -> dict:
    """
    Integrate the multibody system from t=0 to t=t_end.

    Parameters
    ----------
    bodies     : list of Body
    joints     : list of Joint
    q0_full    : (n_dof,) initial generalised positions (full vector)
    qd0_full   : (n_dof,) initial generalised velocities (full vector)
    t_end      : simulation end time [s]
    rtol, atol : solver tolerances

    Returns
    -------
    result : dict with keys:
        't'        : (n_steps,)   time vector
        'q'        : (n_steps, n_dof)  full position history
        'qd'       : (n_steps, n_dof)  full velocity history
        'idx_u'    : list of int — independent DOF indices
        'idx_c'    : list of int — driven DOF indices
        'bodies'   : body list (for post-processing)
        'joints'   : joint list
        'ctx'      : simulation context (t0_motor, etc.)
    """

    idx_u, idx_c = get_dof_indices(bodies)
    nu = len(idx_u)
    n_dof = len(q0_full)

    # Allocate state once (reused inside ODE closure)
    state = MBState(n_bodies=len(bodies) + 1, n_dof=n_dof)

    # Build ODE closure
    f_ode, ctx = make_ode(bodies, joints, state)

    # Initial condition for independent DOFs only
    y0 = np.concatenate([q0_full[idx_u], qd0_full[idx_u]])
    print(f_ode(10.0, y0))
    print(f"[integrator] Starting RK45 integration — t_end={t_end}s, "
          f"nu={nu}, nc={len(idx_c)}")

    # Integrate with RK45 (dense output for post-processing)
    sol = solve_ivp(f_ode,
                    t_span=(0.0, t_end),
                    y0=y0,
                    method='RK45',
                    rtol=rtol,
                    atol=atol,
                    dense_output=False)

    if not sol.success:
        raise RuntimeError(f"[integrator] solve_ivp failed: {sol.message}")

    print(f"[integrator] Done — {sol.t.shape[0]} steps, "
          f"status: {sol.message}")

    # Reconstruct full q, qd history by re-injecting driven values
    n_steps = sol.t.shape[0]
    q_hist  = np.zeros((n_steps, n_dof))
    qd_hist = np.zeros((n_steps, n_dof))

    from merry_go_robotran.simulation.driven_vars import evaluate_driven  # local import to avoid circular

    for k, t_k in enumerate(sol.t):
        q_hist[k,  idx_u] = sol.y[:nu, k]
        qd_hist[k, idx_u] = sol.y[nu:, k]

        driven_vals = evaluate_driven(bodies, t_k)
        for q_idx, (qc, qdc, _qddc) in driven_vals.items():
            q_hist[k,  q_idx] = qc
            qd_hist[k, q_idx] = qdc

    return {
        't':      sol.t,
        'q':      q_hist,
        'qd':     qd_hist,
        'idx_u':  idx_u,
        'idx_c':  idx_c,
        'bodies': bodies,
        'joints': joints,
        'ctx':    ctx,
    }