# =============================================================================
# LMECA2802 - Merry-Go-Robotran
# carousel_data.py
#
# Données complètes du carrousel du MOOC (UCLouvain / Polytechnique Montréal).
# Source : Mooc_MerrygoRound_V1.pdf
#
# Toutes les grandeurs sont exprimées dans le repère inertiel I = (I1, I2, I3)
# avec I3 aligné sur l'axe de rotation du pôle (vertical, vers le haut).
# La gravité est dans la direction -I3, g = 9.81 m/s².
#
# Convention d'arbre (topologie) :
#   0  = base (inertielle, fixe)
#   1  = pôle          (joint R1 : rotation principale autour de I3)
#   2  = pôle tilt 1   (joint R2 : joint Cardan axe I1 - commandé)
#   3  = pôle tilt 2   (joint R3 : joint Cardan axe I2 - commandé)
#   4  = bras 1        (corps rigide, solidaire du pôle — pas de ddl propre)
#   5  = bras 2        (corps rigide, solidaire du pôle — pas de ddl propre)
#   6  = bras 3        (corps rigide, solidaire du pôle — pas de ddl propre)
#   7  = bras 4        (corps rigide, solidaire du pôle — pas de ddl propre)
#   8  = pendule 1     (joint R4 : charnière par rapport au bras 1)
#   9  = pendule 2     (joint R5 : charnière par rapport au bras 2)
#   10 = pendule 3     (joint R6 : charnière par rapport au bras 3)
#   11 = pendule 4     (joint R7 : charnière par rapport au bras 4)
#   12 = nacelle 1     (Cardan 2 ddl par rapport au pendule 1 - commandé)
#   13 = nacelle 2     (Cardan 2 ddl par rapport au pendule 2 - commandé)
#   14 = nacelle 3     (Cardan 2 ddl par rapport au pendule 3 - commandé)
#   15 = nacelle 4     (Cardan 2 ddl par rapport au pendule 4 - commandé)
#
# Remarque : dans une première implémentation simplifiée, les bras sont
# intégrés au corps "pôle" (masse et inertie additionnées). Les nacelles
# sont supposées alignées sur leur pendule (Cardan verrouillé au départ).
# =============================================================================

import numpy as np
from math import pi, cos, sin
from merry_go_robotran.neri.body import Body


# =============================================================================
# 1. CONSTANTES GÉNÉRALES
# =============================================================================

g = 9.81            # [m/s²] accélération gravitationnelle (direction -I3)
t_end = 30.0        # [s]    durée totale de simulation (consigne MOOC)
N_sub = 4           # nombre de sous-systèmes (bras + pendule + nacelle)


# =============================================================================
# 2. TOPOLOGIE DE L'ARBRE MULTICORPS
#
# parent[i] = indice du corps parent de i (0 = base)
# joint_type[i] = type de joint entre i et parent[i]
#   'R' = rotation (1 ddl)
#   'C2' = cardan (2 ddl)
#   'fixed' = corps rigidement lié au parent (0 ddl)
#   'none' = corps de base (pas de joint)
# joint_axis[i] = axe(s) de rotation du joint, exprimé(s) dans le repère parent
# =============================================================================

# Corps :   0       1       2      3       4-7       8-11     12-15
body_names = [
    "base",
    "pole_tilt1",   # Cardan axe I1 (commandé)
    "pole_tilt2",   # Cardan axe I2 (commandé)
    "pole",
    "arm_1", "arm_2", "arm_3", "arm_4",
    "pendulum_1", "pendulum_2", "pendulum_3", "pendulum_4",
    "nacelle_1", "nacelle_2", "nacelle_3", "nacelle_4",
]

# Indice parent de chaque corps
parent = {
    "base": None,
    "pole_tilt1": "base",
    "pole_tilt2": "pole_tilt1",
    "pole": "pole_tilt2",
    "arm_1": "pole",
    "arm_2": "pole",
    "arm_3": "pole",
    "arm_4": "pole",
    "pendulum_1": "arm_1",
    "pendulum_2": "arm_2",
    "pendulum_3": "arm_3",
    "pendulum_4": "arm_4",
    "nacelle_1": "pendulum_1",
    "nacelle_2": "pendulum_2",
    "nacelle_3": "pendulum_3",
    "nacelle_4": "pendulum_4",
}

# Type de joint (côté enfant)
joint_type = {
    "base": "none",
    "pole": "R",          # rotation principale autour de I3 — variable INDÉPENDANTE
    "pole_tilt1": "R",    # tilt Cardan 1 autour de I1 — COMMANDÉ
    "pole_tilt2": "R",    # tilt Cardan 2 autour de I2 — COMMANDÉ
    "arm_1": "fixed",     # bras solidaire du pôle (intégré dans inertie pôle)
    "arm_2": "fixed",
    "arm_3": "fixed",
    "arm_4": "fixed",
    "pendulum_1": "R",    # charnière pendule/bras — variable INDÉPENDANTE
    "pendulum_2": "R",    # charnière pendule/bras — variable INDÉPENDANTE
    "pendulum_3": "R",    # charnière pendule/bras — variable INDÉPENDANTE
    "pendulum_4": "R",    # charnière pendule/bras — variable INDÉPENDANTE
    "nacelle_1": "C2",    # Cardan nacelle/pendule — COMMANDÉ (verrouillé sur pendule)
    "nacelle_2": "C2",
    "nacelle_3": "C2",
    "nacelle_4": "C2",
}

# Variables commandées vs indépendantes
# 'u' = indépendante (intégrée par RK4)
# 'c' = commandée (calculée par loi de commande, non intégrée)
var_type = {
    "pole_tilt1": "c",
    "pole_tilt2": "c",
    "pole": "u",        # rotating at given speed no? not driven, rotating because external torque applied
    "pendulum_1": "u",
    "pendulum_2": "u",
    "pendulum_3": "u",
    "pendulum_4": "u",
    "nacelle_1": "u",   # no driving motions for the nacelle, just hanging there with 2DOF
    "nacelle_2": "u",
    "nacelle_3": "u",
    "nacelle_4": "u",
}


# =============================================================================
# 3. GÉOMÉTRIE
#
# Toutes les distances sont en mètres.
# Les vecteurs sont exprimés dans le repère LOCAL du corps parent.
# "r_parent_to_joint" = vecteur du centre de masse du parent vers l'articulation
# "r_joint_to_com"    = vecteur de l'articulation vers le centre de masse du corps
# =============================================================================

# --- Pôle ---
L_pole = 4.5        # [m] hauteur totale du pôle
r_pole = 0.25       # [m] rayon du cylindre du pôle (diamètre 0.5 m)

# Vecteur base -> joint pied du pôle (Cardan au sol)
# Le pôle est ancré au sol ; le joint est à la base du pôle
r_base_to_pole_joint = np.array([0.0, 0.0, 0.0])

# Centre de masse du pôle = mi-hauteur
r_pole_joint_to_com = np.array([0.0, 0.0, 2])#  2m COM given in doc L_pole / 2.0])   # [0, 0, 2.25] m

# --- Bras (identiques, orientés à 90° l'un de l'autre) ---
# Les 4 bras partent du sommet du pôle (z = L_pole = 4.5 m)
# Chaque bras est dans le plan horizontal ; longueur = 1 m (jusqu'à l'axe de la charnière pendule)
L_arm = 1.0#-0.25         # [m] longueur du bras (axe pôle -> charnière pendule)
# with this arm starts outside of pole

# Angle de chaque bras dans le plan horizontal (quand pôle vertical, repère inertiel)
# Bras 1 : +I1 direction, Bras 2 : +I2, Bras 3 : -I1, Bras 4 : -I2
arm_angles_inertial = {
    "arm_1": 0.0,
    "arm_2": pi / 2.0,
    "arm_3": pi,
    "arm_4": 3.0 * pi / 2.0,
}

# Vecteur sommet pôle -> charnière bras k (dans repère PÔLE après tilt)
# Pour le bras k d'angle phi_k (dans repère pôle) :
#   r = L_arm * [cos(phi_k), sin(phi_k), 0]  +  [0, 0, L_pole/2] pour aller au sommet
# Le joint bras-pôle est au sommet du pôle, donc depuis le CDM pôle :
#   r_com_pole_to_arm_hinge_k = [L_arm*cos(phi_k), L_arm*sin(phi_k), L_pole/2]
# Note : dans le code NERi, r_joint_to_com sera utilisé directement par corps.
# Ici on stocke la géométrie de référence pour chaque bras.

def r_pole_com_to_arm_hinge(arm_name):
    """Vecteur CDM pôle -> charnière bras (dans repère pôle, pôle vertical)."""
    phi = arm_angles_inertial[arm_name]
    return np.array([L_arm * cos(phi), L_arm * sin(phi), 2.5]) # because pole COM at 2m, so not 4.5/2
## Function above is never used??

# Centre de masse des bras : les bras sont intégrés dans le pôle (corps rigide).
# Si on les sépare plus tard, le CDM du bras est à mi-longueur :
r_arm_joint_to_com = np.array([L_arm / 2.0, 0.0, 0.0])  # dans repère local du bras

# --- Pendules ---
L_pendulum = 3.0    # [m] longueur du pendule (charnière -> Cardan nacelle)
# Le centre de masse du pendule est à mi-hauteur (corps cylindrique uniforme)
r_pendulum_hinge_to_com = np.array([0.0, 0.0, -L_pendulum / 2.0])  # vers le bas
r_pendulum_hinge_to_nacelle_joint = np.array([0.0, 0.0, -L_pendulum])

# Axe de rotation de la charnière pendule/bras :
# Perpendiculaire à l'axe du pôle ET à l'axe du bras → direction tangentielle
# Dans le repère local du bras, l'axe de rotation est I2_local = perpendiculaire au bras

# --- Nacelles ---
# Centre de masse de la nacelle : 1 m en dessous du joint Cardan
r_nacelle_cardan_to_com = np.array([0.0, 0.0, -1.0])   # [m] (Fig. 2 : 1 m)

# --- Ancrage du ressort-amortisseur sur le bras ---
# Point d'attache sur le bras : 0.5 m de l'axe du pôle, dans la direction du bras
'''r_spring_arm_attach = 0.5           # [m] depuis axe pôle, dans direction bras'''
## ^^ Not working with linear springs - Fisette

# Point d'attache sur le pendule : 1.5 m sous la charnière et 0.1 m vers l'intérieur
'''r_spring_pendulum_attach_below  = 1.5   # [m] sous la charnière le long du pendule
r_spring_pendulum_attach_inward = 0.1   # [m] vers le centre (inward)'''


# =============================================================================
# 4. MASSES ET INERTIES
#
# Tenseurs d'inertie au CENTRE DE MASSE, exprimés dans les axes PRINCIPAUX
# de chaque corps. Matrice diagonale : I = diag(Ixx, Iyy, Izz)
#
# Convention :
#   "rev"  = autour de l'axe de révolution principal du corps
#   "tra"  = autour des deux axes transverses (supposés égaux par symétrie)
# =============================================================================

# --- Pôle (masse + bras intégrés dans un seul corps pour simplification) ---
# Attention : les 4 bras ne sont pas explicitement décrits dans le MOOC pour
# leur masse propre. Le document donne directement I_pole incluant les bras.
m_pole = 600.0                      # [kg]
I_pole_rev = 16.0                   # [kg·m²] autour de I3 (axe vertical)
I_pole_tra = 458.8                  # [kg·m²] autour de I1 ou I2 (axes transverses)
I_pole = np.diag([I_pole_tra, I_pole_tra, I_pole_rev])   # tenseur 3×3 au CDM

# --- Pendules (tous identiques) ---
m_pendulum = 44.0                   # [kg]
I_pendulum_main = 0.091             # [kg·m²] autour de l'axe principal (axe du pendule)
I_pendulum_tra  = 15.9              # [kg·m²] autour des axes transverses
# Axe principal du pendule = I3_local (axe du cylindre) → Izz = I_main, Ixx=Iyy=I_tra
I_pendulum = np.diag([I_pendulum_tra, I_pendulum_tra, I_pendulum_main])

# Coefficient d'amortissement visqueux dans les charnières pendule/bras
# Pendules 1, 3, 4 : normal ; Pendule 2 : rouillé (valeur ×200)
damping_pendulum_hinge = {
    "pendulum_1": 100.0,     # [N·s/m]
    "pendulum_2": 20000,#20000.0,   # [N·s/m] charnière rouillée
    "pendulum_3": 100.0,
    "pendulum_4": 100.0,
}

# --- Nacelles (toutes identiques) ---
m_nacelle = 123.3                   # [kg]
I_nacelle_sym = 61.7                # [kg·m²] autour de l'axe de symétrie (I_sym)
I_nacelle_tra = 15.4                # [kg·m²] autour des axes transverses
I_nacelle = np.diag([I_nacelle_tra, I_nacelle_tra, I_nacelle_sym])

# Amortissement visqueux dans les deux joints du Cardan nacelle/pendule
damping_nacelle_cardan = 6.0        # [N·s/m] (idem pour les deux axes du Cardan)


# =============================================================================
# 5. RESSORTS-AMORTISSEURS (spring-damper elements)
#
# Un élément ressort-amortisseur agit entre chaque bras et son pendule.
# La force est alignée entre les deux points d'attache.
# =============================================================================

'''k_spring    = 500.0     # [N/m]   raideur du ressort
L0_spring   = 10.0      # [m]     longueur à vide du ressort
c_damper    = 700.0     # [N·s/m] coefficient d'amortissement'''


# =============================================================================
# 6. FORCES ET COUPLES EXTERNES
# =============================================================================

# --- Couple moteur sur la rotation principale du pôle ---
# Loi de commande (voir MOOC p. 6) :
#   Phase 1 : T = 1000 Nm  si ω_pole < 0.8 rad/s (jamais atteint)
#   Phase 2 : T = A(1+cos(ω_ctrl*t + φ)) Nm  pendant 0.5 s après t0
#   Phase 3 : T = 0  ensuite
# où t0 = premier instant où ω_pole atteint 0.8 rad/s

T_motor_phase1   = 1000.0           # [N·m]
T_motor_A        = 500.0            # [N·m]  amplitude phase 2
omega_ctrl_motor = 2.0 * pi         # [rad/s] fréquence phase 2
omega_pole_threshold = 0.8          # [rad/s] seuil de vitesse

def motor_torque(t, t_cross):
    if t < 0.5 + t_cross:
        torque = T_motor_A * (1 + np.cos(omega_ctrl_motor * t - 2 * np.pi * t_cross))
    else:
        torque = 0
    return torque

# --- Force de vent sur chaque nacelle ---
# F = A*(1 - cos(ω*t)) dans la direction I1 (force sinusoïdale)
F_wind_A     = 300.0                # [N]    amplitude
omega_wind   = 6.0 * pi             # [rad/s] fréquence # in direction I1

def ext_force(t):
    F = F_wind_A * (1 - cos(omega_wind * t)) ## in direction I1
    return F

# =============================================================================
# 7. VARIABLES COMMANDÉES : lois de commande
#
# Les angles de tilt du pôle (Cardan) sont des variables commandées.
# θ1,2(t) = A * (1 - cos(ω*t + φ1,2))
#
# A      = 2.5 * (2π/360) rad  ≈ 2.5 degrés en amplitude maximale
# ω      = 0.4π rad/s
# φ1     = π/2
# φ2     = 0
# =============================================================================

tilt_A       = 2.5 * (2.0 * pi / 360.0)    # [rad]    ≈ 0.04363 rad (≈ 2.5°)
tilt_omega   = 0.4 * pi                     # [rad/s]
tilt_phi1    = pi / 2.0                     # [rad]    pour tilt axe I1
tilt_phi2    = 0.0                          # [rad]    pour tilt axe I2


def driven_tilt1(t):
    """Angle commandé du tilt Cardan 1 (autour de I1) [rad]."""
    return tilt_A * (1.0 - cos(tilt_omega * t + tilt_phi1))


def driven_tilt1_dot(t):
    """Vitesse commandée du tilt Cardan 1 [rad/s]."""
    return tilt_A * tilt_omega * sin(tilt_omega * t + tilt_phi1)


def driven_tilt2(t):
    """Angle commandé du tilt Cardan 2 (autour de I2) [rad]."""
    return tilt_A * (1.0 - cos(tilt_omega * t + tilt_phi2))


def driven_tilt2_dot(t):
    """Vitesse commandée du tilt Cardan 2 [rad/s]."""
    return tilt_A * tilt_omega * sin(tilt_omega * t + tilt_phi2)


# =============================================================================
# 8. CONDITIONS INITIALES
#
# Toutes les vitesses initiales sont nulles.
# Angles initiaux :
#   - Rotation principale du pôle : 0 rad
#   - Tilts Cardan : θ1(0) = driven_tilt1(0), θ2(0) = driven_tilt2(0)
#   - Pendule k : 75° = 75 * pi/180 rad (vers l'extérieur, voir Fig. 3)
#   - Nacelles  : alignées sur le pendule (Cardan à 0)
# =============================================================================

theta_pendulum_init = -15.0 * pi / 180.0    # [rad]  ≈ 1.3090 rad

q0 = {
    "pole_tilt1":  driven_tilt1(0.0),       # = 0 car cos(π/2) = 0 → tilt_A*(1-0)
    "pole_tilt2":  driven_tilt2(0.0),       # = 0 car cos(0) = 1   → 0
    "pole": 0.0,
    "pendulum_1":  theta_pendulum_init, #would be negative in my notes
    "pendulum_2":  theta_pendulum_init,
    "pendulum_3":  theta_pendulum_init,
    "pendulum_4":  theta_pendulum_init,
    "upper_pend": theta_pendulum_init,
}

qd0 = {k: 0.0 for k in q0}                # toutes les vitesses initiales nulles


# =============================================================================
# 9. POINT DE RÉFÉRENCE PASSAGER (pour calcul d'accélération - MOOC)
#
# Passager situé dans la nacelle 1 :
#   - 50 cm sous le joint Cardan nacelle-pendule
#   - 37.5 cm vers le centre du pôle (radially inward)
# =============================================================================

#r_passenger_from_nacelle_cardan = np.array([-0.375, 0.0, -0.50])   # [m] dans repère nacelle


# =============================================================================
# 10. RÉSUMÉ SYNTHÉTIQUE (affichage)
# =============================================================================

def print_summary():
    print("=" * 60)
    print("  MOOC Merry-Go-Round — Résumé des paramètres")
    print("=" * 60)
    print(f"\n  Géométrie")
    print(f"    Hauteur pôle          : {L_pole} m")
    print(f"    Rayon pôle            : {r_pole} m")
    print(f"    Longueur bras         : {L_arm} m")
    print(f"    Longueur pendule      : {L_pendulum} m")
    print(f"\n  Masses")
    print(f"    Pôle (+ bras)         : {m_pole} kg")
    print(f"    Pendule               : {m_pendulum} kg  (× {N_sub})")
    print(f"    Nacelle               : {m_nacelle} kg  (× {N_sub})")
    print(f"    Masse totale          : {m_pole + N_sub*(m_pendulum + m_nacelle):.1f} kg")
    print(f"\n  Inerties du pôle")
    print(f"    I_rev (I3)            : {I_pole_rev} kg·m²")
    print(f"    I_tra (I1,I2)         : {I_pole_tra} kg·m²")
    print(f"\n  Ressort-amortisseur")
    #print(f"    k = {k_spring} N/m,  L0 = {L0_spring} m,  c = {c_damper} N·s/m")
    print(f"\n  Couple moteur")
    print(f"    Phase 1 : T = {T_motor_phase1} N·m  (jusqu'à ω = {omega_pole_threshold} rad/s)")
    print(f"    Phase 2 : T = {T_motor_A}*(1+cos(...)) N·m  (0.5 s de transition)")
    print(f"    Phase 3 : T = 0")
    print(f"\n  Conditions initiales")
    print(f"    θ_pendule (init)      : {theta_pendulum_init*180/pi:.1f}°")
    print(f"    Toutes vitesses       : 0 rad/s")
    print(f"\n  Simulation")
    print(f"    Durée totale          : {t_end} s")
    print("=" * 60)


def build_bodies() -> list[Body]:
    """
    Instantiate all 9 bodies of the carousel using the parameters
    defined above in this file, and return them in topological order
    (parent always before child).
    """
    bodies = []

    # --- Body 1: Pole ---
    bodies.append(Body(
        body_id   = 1,
        name      = 'pole',
        parent_id = 0,
        mass              = 0,
        inertia_com_local = np.zeros((3,3)),
        d_parent_to_joint_in_parent = np.zeros(3),      # [0,0,0]
        d_joint_to_com_in_local     = np.zeros(3),   #r_pole_joint_to_com not this one because no mass # [0,0,2]
        joint_type           = 'revolute',
        joint_axis_in_parent = np.array([0, 0, 1]),
        n_dof     = 1,
        var_types = ['u'],
        q_indices = [0],
    ))
    bodies.append(Body(
        body_id = 2,
        name = 'pole_tilt1',
        parent_id = 1,
        mass = 0,
        inertia_com_local = np.zeros((3,3)),
        d_parent_to_joint_in_parent=np.zeros(3),
        d_joint_to_com_in_local=np.zeros(3),
        joint_type='revolute',
        joint_axis_in_parent=np.array([1, 0, 0]),
        n_dof = 1,
        var_types = ['c'],
        q_indices = [1],
    ))
    bodies.append(Body(
        body_id=3,
        name='pole_tilt2',
        parent_id=2,
        mass=m_pole,
        inertia_com_local=I_pole,
        d_parent_to_joint_in_parent=r_base_to_pole_joint,
        d_joint_to_com_in_local=r_pole_joint_to_com,
        joint_type='revolute',
        joint_axis_in_parent=np.array([0, 1, 0]),
        n_dof=1,
        var_types=['c'],
        q_indices=[2],
    ))

    # --- Bodies 2–9: 4× (pendulum + nacelle) ---
    q_idx = 3   # global q index counter, starts after pole's 3 DOFs

    for k in range(N_sub):
        pend_id   = 4 + 2*k     # 2, 4, 6, 8  all +2
        nacelle_id = 5 + 2*k    # 3, 5, 7, 9  all +2
        label = k + 1           # 1, 2, 3, 4

        # Arm tip location in pole frame (varies per sub-system)
        phi   = arm_angles_inertial[f'arm_{label}']
        arm_tip_in_pole = np.array([L_arm*cos(phi), L_arm*sin(phi), L_pole])

        # Hinge axis: perpendicular to both pole axis (I3) and arm direction
        # = tangential direction in the horizontal plane
        arm_dir   = np.array([cos(phi), sin(phi), 0.0])
        hinge_axis = np.array([-sin(phi), cos(phi), 0.0])  # I3 × arm_dir

        # Damping value differs for pendulum 2 (rusted hinge)
        pend_name = f'pendulum_{label}'
        #no arm (included in pole), straight to pendulum)
        bodies.append(Body.revolute(
            body_id   = pend_id,
            name      = pend_name,
            parent_id = 3,                              # all pendulums attach to pole
            mass         = m_pendulum,
            inertia_diag = [I_pendulum_tra, I_pendulum_tra, I_pendulum_main],
            d_parent_to_joint = arm_tip_in_pole,          # skipping arms as they are attached, no joint no motion
            d_joint_to_com    = r_pendulum_hinge_to_com,  # [0,0,-1.5]
            axis_in_parent    = hinge_axis,
            var_type = 'u',
            q_index  = q_idx,
        ))
        q_idx += 1

        bodies.append(Body.cardan2(
            body_id   = nacelle_id,
            name      = f'nacelle_{label}',
            parent_id = pend_id,
            mass         = m_nacelle,
            inertia_diag = [I_nacelle_tra, I_nacelle_tra, I_nacelle_sym],
            d_parent_to_joint = r_pendulum_hinge_to_nacelle_joint,  # [0,0,-3]
            d_joint_to_com    = r_nacelle_cardan_to_com,             # [0,0,-1]
            axes_in_parent    = [hinge_axis,
                                 arm_dir],    # axis 1 = hinge axis, axis 2 = arm direction
            var_types = ['u', 'u'],
            q_indices = [q_idx, q_idx+1],
        ))
        q_idx += 2

    return bodies

def build_bodies_split(len_section: float) -> list[Body]:
    """
        Same as build_bodies() but for the last pendulum arm, splits it in 2
        and add driven condition of no motion so forces can be measured and
        dimensioning can be done.
        """
    bodies = []

    # --- Body 1: Pole ---
    bodies.append(Body(
        body_id=1,
        name='pole',
        parent_id=0,
        mass=0,
        inertia_com_local=np.zeros((3, 3)),
        d_parent_to_joint_in_parent=np.zeros(3),  # [0,0,0]
        d_joint_to_com_in_local=np.zeros(3),  # r_pole_joint_to_com not this one because no mass # [0,0,2]
        joint_type='revolute',
        joint_axis_in_parent=np.array([0, 0, 1]),
        n_dof=1,
        var_types=['u'],
        q_indices=[0],
    ))
    bodies.append(Body(
        body_id=2,
        name='pole_tilt1',
        parent_id=1,
        mass=0,
        inertia_com_local=np.zeros((3, 3)),
        d_parent_to_joint_in_parent=np.zeros(3),
        d_joint_to_com_in_local=np.zeros(3),
        joint_type='revolute',
        joint_axis_in_parent=np.array([1, 0, 0]),
        n_dof=1,
        var_types=['c'],
        q_indices=[1],
    ))
    bodies.append(Body(
        body_id=3,
        name='pole_tilt2',
        parent_id=2,
        mass=m_pole,
        inertia_com_local=I_pole,
        d_parent_to_joint_in_parent=r_base_to_pole_joint,
        d_joint_to_com_in_local=r_pole_joint_to_com,
        joint_type='revolute',
        joint_axis_in_parent=np.array([0, 1, 0]),
        n_dof=1,
        var_types=['c'],
        q_indices=[2],
    ))

    # --- Bodies 2–9: 4× (pendulum + nacelle) ---
    q_idx = 3  # global q index counter, starts after pole's 3 DOFs

    for k in range(N_sub-1):
        pend_id = 4 + 2 * k  #  4, 6, 8
        nacelle_id = 5 + 2 * k  # 5, 7, 9
        label = k + 1  # 1, 2, 3

        # Arm tip location in pole frame (varies per sub-system)
        phi = arm_angles_inertial[f'arm_{label}']
        arm_tip_in_pole = np.array([L_arm * cos(phi), L_arm * sin(phi), L_pole])

        # Hinge axis: perpendicular to both pole axis (I3) and arm direction
        # = tangential direction in the horizontal plane
        arm_dir = np.array([cos(phi), sin(phi), 0.0])
        hinge_axis = np.array([-sin(phi), cos(phi), 0.0])  # I3 × arm_dir

        # Damping value differs for pendulum 2 (rusted hinge)
        pend_name = f'pendulum_{label}'
        # no arm (included in pole), straight to pendulum)
        bodies.append(Body.revolute(
            body_id=pend_id,
            name=pend_name,
            parent_id=3,  # all pendulums attach to pole
            mass=m_pendulum,
            inertia_diag=[I_pendulum_tra, I_pendulum_tra, I_pendulum_main],
            d_parent_to_joint=arm_tip_in_pole,  # skipping arms as they are attached, no joint no motion
            d_joint_to_com=r_pendulum_hinge_to_com,  # [0,0,-1.5]
            axis_in_parent=hinge_axis,
            var_type='u',
            q_index=q_idx,
        ))
        q_idx += 1

        bodies.append(Body.cardan2(
            body_id=nacelle_id,
            name=f'nacelle_{label}',
            parent_id=pend_id,
            mass=m_nacelle,
            inertia_diag=[I_nacelle_tra, I_nacelle_tra, I_nacelle_sym],
            d_parent_to_joint=r_pendulum_hinge_to_nacelle_joint,  # [0,0,-3]
            d_joint_to_com=r_nacelle_cardan_to_com,  # [0,0,-1]
            axes_in_parent=[hinge_axis,
                            arm_dir],  # axis 1 = hinge axis, axis 2 = arm direction
            var_types=['u', 'u'],
            q_indices=[q_idx, q_idx + 1],
        ))
        q_idx += 2
    # Arm tip location in pole frame (varies per sub-system)
    phi = arm_angles_inertial[f'arm_{label}']
    arm_tip_in_pole = np.array([L_arm * cos(phi), L_arm * sin(phi), L_pole])
        # Hinge axis: perpendicular to both pole axis (I3) and arm direction
        # = tangential direction in the horizontal plane
    arm_dir = np.array([cos(phi), sin(phi), 0.0])
    hinge_axis = np.array([-sin(phi), cos(phi), 0.0])  # I3 × arm_dir
    # calculating the new variables for upper and lower part of the pendulum arm
    m_upper = len_section * m_pendulum / L_pendulum
    m_lower = m_pendulum - m_upper
    I_tra_upper = I_pendulum_tra * len_section / L_pendulum
    I_tra_lower = I_pendulum_tra - I_tra_upper
    # for inertia of a slender bar at its end: I = 1/3 ML**2
    # here M is linear w/ len_section, and L too so len_section**3
    I_main_upper = I_pendulum_main * (len_section / L_pendulum)**3
    I_main_lower = I_pendulum_main * (1 - len_section/L_pendulum)**3



        # no arm (included in pole), straight to pendulum)
    bodies.append(Body.revolute(
        body_id=10,
        name="upper_pend",
        parent_id=3,  # all pendulums attach to pole
        mass=m_upper,
        inertia_diag=[I_tra_upper, I_tra_upper, I_main_upper],
        d_parent_to_joint=arm_tip_in_pole,  # skipping arms as they are attached, no joint no motion
        d_joint_to_com=[0, 0, len_section/2],  # [0,0,-1.5]
        axis_in_parent=hinge_axis,
        var_type='u',
        q_index=q_idx,
    ))
    q_idx += 1
    bodies.append(Body.cardan3(
        body_id=11,
        name="lower_pend",
        parent_id=10,  # all pendulums attach to pole
        mass=m_lower,
        inertia_diag=[I_tra_lower, I_tra_lower, I_main_lower],
        d_parent_to_joint=[0, 0, -len_section],  # straight down with length of split section
        d_joint_to_com=[0, 0, (L_pendulum-len_section)/2],  # [0,0,-1.5]
        axes_in_parent=[[1, 0, 0], [0, 1, 0], [0, 0 ,1]],
        var_types=['c','c','c'],
        q_indices=[q_idx, q_idx+1, q_idx+2],
    ))
    q_idx += 3

    bodies.append(Body.cardan2(
        body_id=12,
        name='nacelle_4',
        parent_id=11,
        mass=m_nacelle,
        inertia_diag=[I_nacelle_tra, I_nacelle_tra, I_nacelle_sym],
        d_parent_to_joint=r_pendulum_hinge_to_nacelle_joint,  # [0,0,-3]
        d_joint_to_com=r_nacelle_cardan_to_com,  # [0,0,-1]
        axes_in_parent=[hinge_axis, arm_dir],  # axis 1 = hinge axis, axis 2 = arm direction
        var_types=['u', 'u'],
        q_indices=[q_idx, q_idx + 1],
    ))
    q_idx += 2
    return bodies
#if __name__ == "__main__":
#    print_summary()