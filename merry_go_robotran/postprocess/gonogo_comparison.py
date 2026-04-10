"""
validation/gonogo_comparison.py
--------------------------------
Computes the absolute acceleration at the passenger reference point in
nacelle_1, as required by the MOOC Go/No-Go validation criterion.

Passenger location (from carousel_data.py / MOOC document)
───────────────────────────────────────────────────────────
  50 cm below the nacelle Cardan joint, 37.5 cm radially inward,
  expressed in the nacelle local frame:
      r_pass_local = [-0.375, 0.0, -0.50]  m

Acceleration formula (rigid-body kinematics)
────────────────────────────────────────────
  p_pass = p_O_nacelle + R_nacelle @ r_pass_local

  ẍ_pass = α^nacelle + β^nacelle @ (p_pass - p_O_nacelle)
         + g                        (re-add gravity, stripped by NERi convention)

  |a_pass| = ||ẍ_pass||

  (The NERi α already absorbs gravity in its recursion, so we add g back
   when reporting the physical absolute acceleration felt by the passenger.)

Public API
──────────
    compute_passenger_acceleration(result) → a_mag_array  (n_steps,)
"""

import numpy as np
from neri.state   import MBState
from neri.forward import forward_pass
from neri.backward import backward_pass
import data.carousel_data as cd

# Gravity vector (inertial frame)
_G_VEC = np.array([0.0, 0.0, -9.81])


def _joint_origin_chain(body_id: int, bodies: list, state: MBState) -> np.ndarray:
    """
    Compute the inertial position of joint origin O^i by summing d_hi
    vectors from the base to body_id.
    Requires that forward_pass has already been called on `state`.
    """
    pos = np.zeros(3)
    bid = body_id
    chain = []
    # Walk up to base collecting body ids
    while bid != 0:
        chain.append(bid)
        body = next(b for b in bodies if b.body_id == bid)
        bid  = body.parent_id
    # Sum d_hi from root downward
    for b in reversed(chain):
        pos += state.d_hi[b]
    return pos


def compute_passenger_acceleration(result: dict) -> np.ndarray:
    """
    Evaluate the magnitude of the absolute (inertial-frame) acceleration
    at the passenger reference point in nacelle_1, for every timestep.

    Parameters
    ----------
    result : dict from integrator.run_simulation()

    Returns
    -------
    a_mag : (n_steps,) array of |a_passenger|  [m/s²]
    """

    bodies = result['bodies']
    joints = result['joints']
    t_arr  = result['t']
    q_hist = result['q']
    qd_hist= result['qd']
    n_steps= len(t_arr)
    n_dof  = q_hist.shape[1]

    # Locate nacelle_1 body
    nacelle_body  = next(b for b in bodies if b.name == "nacelle_1")
    nacelle_bid   = nacelle_body.body_id

    # Passenger position in nacelle local frame
    r_pass_local = cd.r_passenger_from_nacelle_cardan  # (3,)

    a_mag = np.zeros(n_steps)

    # Allocate one reusable state
    tmp_state = MBState(n_bodies=len(bodies) + 1, n_dof=n_dof)

    for k, t_k in enumerate(t_arr):
        # Reload and run forward pass for this timestep
        tmp_state.set_q(q_hist[k].copy(), qd_hist[k].copy(), t_k)
        forward_pass(bodies, joints, tmp_state)

        # Inertial position of nacelle joint origin O^nacelle
        p_O_nac = _joint_origin_chain(nacelle_bid, bodies, tmp_state)

        # Passenger position in inertial frame
        p_pass = p_O_nac + tmp_state.R[nacelle_bid] @ r_pass_local

        # Vector from joint origin to passenger (inertial frame)
        d_pass = p_pass - p_O_nac   # = R_nac @ r_pass_local

        # Absolute acceleration of the passenger (rigid-body formula)
        # ẍ_pass = α^nac + β^nac @ d_pass  +  g  (re-add gravity)
        # The NERi α encodes (ẍ_joint - g), so α + β@d gives (ẍ_pass - g).
        a_inertial_minus_g = (tmp_state.alpha[nacelle_bid]
                               + tmp_state.beta[nacelle_bid] @ d_pass)
        a_absolute = a_inertial_minus_g + _G_VEC

        a_mag[k] = np.linalg.norm(a_absolute)

    return a_mag


def load_mooc_reference(csv_path: str) -> tuple:
    """
    Load the MOOC Go/No-Go reference data from a two-column CSV.

    Parameters
    ----------
    csv_path : path to CSV file  (columns: time [s], |a| [m/s²])

    Returns
    -------
    t_ref   : (n,) time array [s]
    a_ref   : (n,) acceleration array [m/s²]
    """
    data = np.loadtxt(csv_path, delimiter=',', comments='#')
    return data[:, 0], data[:, 1]


def gonogo_check(a_sim: np.ndarray,
                 a_ref: np.ndarray,
                 t_sim: np.ndarray,
                 t_ref: np.ndarray,
                 tol:   float = 0.10) -> bool:
    """
    Compare the simulated passenger acceleration against the MOOC reference.

    Interpolates the reference onto the simulation time grid and checks
    that the relative RMS deviation is within `tol` (default 10%).

    Parameters
    ----------
    a_sim  : (n_sim,) simulated |a| [m/s²]
    a_ref  : (n_ref,) reference |a| [m/s²]
    t_sim  : (n_sim,) simulation time [s]
    t_ref  : (n_ref,) reference time [s]
    tol    : maximum relative RMS deviation (0.10 = 10%)

    Returns
    -------
    passed : bool  — True if within tolerance (GO), False otherwise (NO-GO)
    """
    # Interpolate reference onto simulation time grid
    a_ref_interp = np.interp(t_sim, t_ref, a_ref)

    # Relative RMS deviation
    rms_ref = np.sqrt(np.mean(a_ref_interp**2))
    rms_err = np.sqrt(np.mean((a_sim - a_ref_interp)**2))
    rel_err = rms_err / (rms_ref + 1e-12)

    print(f"[gonogo] RMS error = {rel_err*100:.2f}%  "
          f"({'GO ✓' if rel_err <= tol else 'NO-GO ✗'})")

    return rel_err <= tol