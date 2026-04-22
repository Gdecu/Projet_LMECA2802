# -*- coding: utf-8 -*-
"""Module for the definition of driven joints."""
# Author: Robotran Team
# (c) Universite catholique de Louvain, 2020

import numpy as np

def user_DrivenJoints(mbs_data, tsim):
    """Set the values of the driven joints directly in the MbsData structure.

    The position, velocity and acceleration of the driven joints must be set in
    the attributes mbs_data.q, mbs_data.qd and mbs_data.qdd .

    Parameters
    ----------
    mbs_data : MBsysPy.MbsData
        The multibody system associated to this computation.
    tsim : float
        The current time of the simulation.

    Returns
    -------
    None
    """

    # Example: joint 5 under constant acceleration with non-zero initial
    #          coordinate (mbs_data.q0) and velocity (mbs_data.qd0).
    # mbs_data.qdd[5] = 2
    # mbs_data.qd[5]  = mbs_data.qd0[5] + mbs_data.qdd[5]*tsim
    # mbs_data.q[5]   = mbs_data.q0[5]  + mbs_data.qd0[5]*tsim + 0.5 * mbs_data.qdd[5]*tsim*tsim

    # Récupération des ID 
    id_tilt1 = mbs_data.joint_id["Tilt_angle_1"]
    id_tilt2 = mbs_data.joint_id["Tilt_angle_2"]

    # Paramètres du MOOC
    A = 2.5 * (2.0 * np.pi / 360.0) # Amplitude en radians
    omega = 0.4 * np.pi
    phi1 = np.pi / 2.0
    phi2 = 0.0

    # Tilt 1 (axe I1)
    mbs_data.q[id_tilt1]   = A * (1.0 - np.cos(omega * tsim + phi1))
    mbs_data.qd[id_tilt1]  = A * omega * np.sin(omega * tsim + phi1)
    mbs_data.qdd[id_tilt1] = A * (omega**2) * np.cos(omega * tsim + phi1)

    # Tilt 2 (axe I2)
    mbs_data.q[id_tilt2]   = A * (1.0 - np.cos(omega * tsim + phi2))
    mbs_data.qd[id_tilt2]  = A * omega * np.sin(omega * tsim + phi2)
    mbs_data.qdd[id_tilt2] = A * (omega**2) * np.cos(omega * tsim + phi2)

    return
