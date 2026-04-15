# -*- coding: utf-8 -*-
"""Module for the definition of joint forces."""
# Author: Robotran Team
# (c) Universite catholique de Louvain, 2020

import numpy as np


def user_JointForces(mbs_data, tsim):
    """Compute the force and torques in the joint.

    It fills the MBsysPy.MbsData.Qq array.

    Parameters
    ----------
    mbs_data : MBsysPy.MbsData
        The multibody system associated to this computation.
    tsim : float
        The current time of the simulation.

    Notes
    -----
    The numpy.ndarray MBsysPy.MbsData.Qq is 1D array with index starting at 1.
    The first index (array[0]) must not be modified. The first index to be
    filled is array[1].

    Returns
    -------
    None
    """
    # cleaning previous forces value
    mbs_data.Qq[1:] = 0.


    # Example: damping in joint number 5
    # D = 0.5 # N/(m/s)
    # mbs_data.Qq[5] = -D * mbs_data.qd[5]

    # Nettoyage obligatoire des forces précédentes
    mbs_data.Qq[1:] = 0.

    # MOTEUR DU POLE PRINCIPAL 
    id_pole = mbs_data.joint_id["Pole_rotation"]
    omega_pole = mbs_data.qd[id_pole]
    
    # Astuce Python : Créer une variable dynamique dans mbs_data pour retenir t0
    if not hasattr(mbs_data, 't0_motor'):
        mbs_data.t0_motor = -1.0 # -1 signifie que le seuil n'a pas encore été atteint
        
    # Vérification du seuil de vitesse (0.8 rad/s)
    if mbs_data.t0_motor < 0.0 and omega_pole >= 0.8:
        mbs_data.t0_motor = tsim # On enregistre l'instant t0
        
    # Lois des 3 phases
    if mbs_data.t0_motor < 0.0:
        # Phase 1 : Vitesse < 0.8 rad/s
        mbs_data.Qq[id_pole] = 1000.0
    elif tsim < mbs_data.t0_motor + 0.5:
        # Phase 2 : Transition de 0.5 secondes
        mbs_data.Qq[id_pole] = 500.0 * (1.0 + np.cos(2.0 * np.pi * tsim - 2.0 * np.pi * mbs_data.t0_motor))
    else:
        # Phase 3 : Moteur coupé
        mbs_data.Qq[id_pole] = 0.0


    # AMORTISSEMENT DES PENDULES
    id_p1 = mbs_data.joint_id["Hinge_Pend_1"]
    id_p2 = mbs_data.joint_id["Hinge_Pend_2"]
    id_p3 = mbs_data.joint_id["Hinge_Pend_3"]
    id_p4 = mbs_data.joint_id["Hinge_Pend_4"]

    mbs_data.Qq[id_p1] = -100.0 * mbs_data.qd[id_p1]
    mbs_data.Qq[id_p2] = -20000.0 * mbs_data.qd[id_p2]  # Charnière rouillée !
    mbs_data.Qq[id_p3] = -100.0 * mbs_data.qd[id_p3]
    mbs_data.Qq[id_p4] = -100.0 * mbs_data.qd[id_p4]


    # AMORTISSEMENT DES CARDANS DES NACELLES
    nacelle_joints = [
        "R1_nacelle_1", "R2_nacelle_1",
        "R1_nacelle_2", "R2_nacelle_2",
        "R1_nacelle_3", "R2_nacelle_3",
        "R1_nacelle_4", "R2_nacelle_4"
    ]
    for j_name in nacelle_joints:
        j_id = mbs_data.joint_id[j_name]
        mbs_data.Qq[j_id] = -6.0 * mbs_data.qd[j_id]

    return
