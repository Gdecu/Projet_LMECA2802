"""
postprocess/internal_forces.py
-------------------------------
Computes internal forces and bending moments along the pendulum beam
at any cross-section, for all timesteps of a completed simulation.

Theory (Euler-Bernoulli beam, free-body diagram from the free end)
──────────────────────────────────────────────────────────────────
For a cross-section at curvilinear abscissa  s  (measured from the
nacelle Cardan joint downward = 0, toward the arm hinge = L_pend):

  Cut the pendulum at s.  The "free" portion (length = L_pend - s,
  carrying the nacelle + passengers) is isolated.

  Internal forces at cut section:
    N(s,t)  = normal force  (along pendulum axis)
    V(s,t)  = shear force   (perpendicular to axis)
    Mb(s,t) = bending moment

  These are obtained from Newton's 2nd law applied to the free portion,
  using the inertial-frame forces and moments already computed by the
  NERi backward pass (state.F, state.L).

Public API
──────────
    compute_section_loads(s, pendulum_id, nacelle_id, bodies, state)
        → N, V_vec, Mb_vec  at a single section s and a single instant

    scan_all_sections(n_sections, result_dict, bodies, joints, state)
        → loads_history  shape (n_steps, n_sections, 3_components)
          Envelope of bending moment along the beam over the full simulation.

    find_critical_section(Mb_history)
        → s_crit, t_crit, Mb_max  (location and time of worst bending)
"""

import numpy as np
from merry_go_robotran.neri.state   import MBState
from merry_go_robotran.neri.forward import forward_pass
from merry_go_robotran.neri.backward import backward_pass
import merry_go_robotran.data.carousel_data as cd


# ═══════════════════════════════════════════════════════════════════════════ #
#  Single cross-section loader                                                #
# ═══════════════════════════════════════════════════════════════════════════ #

def compute_section_loads(s:           float,
                          pend_body,
                          nacelle_body,
                          state:       MBState
                          ) -> tuple:
    """
    Internal loads at a cross-section  s  of the pendulum.

    s = 0  at the nacelle Cardan joint (free end of the pendulum).
    s = L  at the arm hinge            (fixed end of the pendulum).

    The free-body diagram isolates the portion [0, s]:
      - Body forces: weight + inertial loads from nacelle (already in state.W)
      - Constraint: the reaction at the cut section = internal loads

    From equilibrium of the free portion:
        F_cut = state.F[nacelle_id]           (reaction from parent chain)
        M_cut = state.L[nacelle_id] + r × F_cut
    where r is from the nacelle Cardan joint to the cut section.

    Parameters
    ----------
    s            : abscissa from nacelle joint [m], 0 ≤ s ≤ L_pendulum
    pend_body    : Body object for the pendulum
    nacelle_body : Body object for the nacelle
    state        : MBState after a complete forward+backward pass

    Returns
    -------
    F_cut  : (3,) internal force vector at cut [N]   (inertial frame)
    Mb_cut : (3,) internal moment vector at cut [N·m] (inertial frame)
    """

    i_pend   = pend_body.body_id
    i_nacelle = nacelle_body.body_id

    # The force transmitted across the section is the sum of all forces
    # acting on the free portion (nacelle side).
    # The backward pass already accumulated these in state.F[nacelle]:
    #   state.F[nacelle] = total resultant (inertia + ext) of nacelle subtree.
    F_cut = state.F[i_nacelle].copy()   # [N], inertial frame

    # Moment at the Cardan joint (nacelle attachment) from the backward pass:
    #   state.L[nacelle] already includes nacelle inertia and any children.
    M_at_cardan = state.L[i_nacelle].copy()   # [N·m], inertial frame

    # Transport moment to cut at abscissa s along the pendulum axis.
    # The pendulum local z-axis points from hinge downward:
    #   pendulum_axis_inertial = R[pend] @ [0,0,-1]  (down when vertical)
    pend_axis = state.R[i_pend] @ np.array([0.0, 0.0, -1.0])

    # Vector from nacelle Cardan joint to cut section
    # Nacelle joint = bottom of pendulum (s=0); cut is s metres up the pendulum
    r_to_cut = s * pend_axis   # pointing up, along pendulum axis

    # Moment at cut = moment at nacelle joint + r_cut × F_cut
    Mb_cut = M_at_cardan + np.cross(r_to_cut, F_cut)

    return F_cut, Mb_cut


# ═══════════════════════════════════════════════════════════════════════════ #
#  Full simulation scan                                                       #
# ═══════════════════════════════════════════════════════════════════════════ #

def scan_all_sections(n_sections:  int,
                      result:      dict,
                      bodies:      list,
                      joints:      list
                      ) -> dict:
    """
    Evaluate internal loads at  n_sections  cross-sections, for every
    timestep stored in `result`, for pendulum_1 (the reference pendulum).

    Parameters
    ----------
    n_sections : number of cross-sections (uniformly spaced along L_pend)
    result     : output dict from integrator.run_simulation()
    bodies     : list of Body
    joints     : list of Joint

    Returns
    -------
    scan : dict with keys:
        's_array'   : (n_sections,) abscissa values [m]
        'F_history' : (n_steps, n_sections, 3) internal force [N]
        'Mb_history': (n_steps, n_sections, 3) bending moment [N·m]
        'Mb_norm'   : (n_steps, n_sections)   |Mb| envelope [N·m]
    """

    # Locate bodies of interest for pendulum_1
    pend_body    = next(b for b in bodies if b.name == "pendulum_1")
    nacelle_body = next(b for b in bodies if b.name == "nacelle_1")

    t_arr   = result['t']
    q_hist  = result['q']
    qd_hist = result['qd']
    n_steps = len(t_arr)
    n_dof   = q_hist.shape[1]

    # Section abscissa array: 0 (nacelle end) → L_pend (arm hinge)
    s_array = np.linspace(0.0, cd.L_pendulum, n_sections)

    F_history  = np.zeros((n_steps, n_sections, 3))
    Mb_history = np.zeros((n_steps, n_sections, 3))

    # Re-use a single MBState for all re-evaluations
    tmp_state = MBState(n_bodies=len(bodies) + 1, n_dof=n_dof)

    for k, t_k in enumerate(t_arr):
        # Reload state for this timestep
        tmp_state.set_q(q_hist[k].copy(), qd_hist[k].copy(), t_k)
        forward_pass(bodies, joints, tmp_state)
        backward_pass(bodies, joints, tmp_state)

        for j, s_j in enumerate(s_array):
            F_cut, Mb_cut = compute_section_loads(
                s_j, pend_body, nacelle_body, tmp_state
            )
            F_history[k,  j] = F_cut
            Mb_history[k, j] = Mb_cut

    Mb_norm = np.linalg.norm(Mb_history, axis=2)   # (n_steps, n_sections)

    return {
        's_array':    s_array,
        'F_history':  F_history,
        'Mb_history': Mb_history,
        'Mb_norm':    Mb_norm,
    }


# ═══════════════════════════════════════════════════════════════════════════ #
#  Critical section finder                                                    #
# ═══════════════════════════════════════════════════════════════════════════ #

def find_critical_section(scan: dict) -> tuple:
    """
    Identify the dynamically critical cross-section: the location s* that
    experiences the maximum bending moment magnitude over the entire
    simulation, and the time t* at which it occurs.

    Parameters
    ----------
    scan : dict returned by scan_all_sections()

    Returns
    -------
    s_crit  : float  — critical abscissa [m]
    t_crit  : float  — time of maximum [s]  (requires result['t'])
    Mb_max  : float  — maximum bending moment magnitude [N·m]
    idx_s   : int    — section index
    idx_t   : int    — timestep index
    """
    Mb_norm = scan['Mb_norm']           # (n_steps, n_sections)
    s_array = scan['s_array']

    # Global maximum over all time and all sections
    idx_flat = np.argmax(Mb_norm)
    idx_t, idx_s = np.unravel_index(idx_flat, Mb_norm.shape)

    Mb_max = Mb_norm[idx_t, idx_s]
    s_crit = s_array[idx_s]

    return s_crit, idx_t, Mb_max, idx_s, idx_t