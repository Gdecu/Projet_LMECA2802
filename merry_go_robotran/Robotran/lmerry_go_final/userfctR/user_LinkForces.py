# -*- coding: utf-8 -*-
"""Module for the definition of user links forces."""
# Author: Robotran Team
# (c) Universite catholique de Louvain, 2020


def user_LinkForces(Z, Zd, mbs_data, tsim, identity):
    """Compute the force in the given link.

    Parameters
    ----------
    Z : float
        The distance between the two anchor points of the link.
    Zd : float
        The relative velocity between the two anchor points of the link.
    mbs_data : MBsysPy.MbsData
        The multibody system associated to this computation.
    tsim : float
        The current time of the simulation.
    identity : int
        The identity of the computed link.

    Returns
    -------
    Flink : float
        The force in the current link.

    """

    Flink = 0.0

    # Example: linear spring
    # k = 1000 #N/m
    # Z0= 0.1  #m
    # Flink = k*(Z-Z0)

    '''Flink = 0.0

    # Regrouper les ID de tes 4 liens créés dans MBsysPad
    spring_ids = [
        mbs_data.link_id["Link_0"],
        mbs_data.link_id["Link_1"],
        mbs_data.link_id["Link_2"],
        mbs_data.link_id["Link_3"]
    ]

    # Si le lien évalué par Robotran est un de nos ressorts
    if identity in spring_ids:
        K = 500.0  # Raideur N/m
        Z0 = 10.0  # Longueur à vide en mètres
        C = 700.0  # Coefficient de damping N/(m/s)
        
        # Uniquement la force élastique (traction si Z > Z0)
        Flink = K * (Z - Z0)
        Flink = K * (Z - Z0) + C * Zd # Force élastique + force de damping
        '''
    return Flink
