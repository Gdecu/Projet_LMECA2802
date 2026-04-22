#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script to run a direct dynamic analysis on a multibody system.

Summary
-------
This template loads the data file *.mbs and execute:
 - the coordinate partitioning module
 - the direct dynamic module (time integration of equations of motion).
 - if available, plot the time evolution of the first generalized coordinate.

It may have to be adapted and completed by the user.


Universite catholique de Louvain
CEREM : Centre for research in mechatronics

http://www.robotran.eu
Contact : info@robotran.be

(c) Universite catholique de Louvain
"""

# %%============================================================================
# Packages loading
# =============================================================================
try:
    import MBsysPy as Robotran
except:
    raise ImportError("MBsysPy not found/installed."
                      "See: https://www.robotran.eu/download/how-to-install/"
                      )
import numpy as np

# %%===========================================================================
# Project loading
# =============================================================================
print("Loading data file...")
mbs_data = Robotran.MbsData('../dataR/merry_go.mbs')
print("Data file loaded successfully.\n")
# %%===========================================================================
# Partitionning
# =============================================================================
print("Running partitioning...")
mbs_data.process = 1
mbs_part = Robotran.MbsPart(mbs_data)
mbs_part.set_options(rowperm=1, verbose=1)
mbs_part.run()
print("Partitioning completed successfully.\n")

# %%===========================================================================
# Direct Dynamics
# =============================================================================
print("Running direct dynamics...")
mbs_data.process = 3
mbs_dirdyn = Robotran.MbsDirdyn(mbs_data)
mbs_dirdyn.set_options(dt0=1e-3, tf=10.0, save2file=1)
results = mbs_dirdyn.run()
print("Direct dynamics completed successfully.\n")

# %%===========================================================================
# Plotting results
# =============================================================================
try:
    import matplotlib.pyplot as plt
except Exception:
    raise RuntimeError('Unable to load matplotlib, plotting results unavailable.')

print("Plotting main pole angular position, velocity and acceleration...")
# Récupération de l'ID du joint du pôle principal (celui qui tourne)
id_pole = mbs_data.joint_id["Pole_rotation"]

# Dans les résultats Robotran, la colonne 0 est le temps.
t = results.q[:, 0]
theta_rad = results.q[:, id_pole]
omega_rad_s = results.qd[:, id_pole]
alpha_rad_s2 = results.qdd[:, id_pole]

# Conversion en degrés
theta_deg = theta_rad * (180.0 / np.pi)
omega_deg_s = omega_rad_s * (180.0 / np.pi)
alpha_deg_s2 = alpha_rad_s2 * (180.0 / np.pi)

# Création de 3 sous-graphes horizontaux: position, vitesse, accélération
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)

axes[0].plot(t, theta_deg, color='darkblue')
axes[0].grid(True, linestyle='--', alpha=0.7)
axes[0].set_xlabel('Time [s]', fontsize=11)
axes[0].set_ylabel('Angular position [deg]', fontsize=11)
axes[0].set_title('Main pole angular position', fontsize=12)
axes[0].set_xlim(0, 3)

axes[1].plot(t, omega_deg_s, color='darkblue')
axes[1].grid(True, linestyle='--', alpha=0.7)
axes[1].set_xlabel('Time [s]', fontsize=11)
axes[1].set_ylabel('Angular velocity [deg/s]', fontsize=11)
axes[1].set_title('Main pole angular velocity', fontsize=12)
axes[1].set_xlim(0, 3)

axes[2].plot(t, alpha_deg_s2, color='darkblue')
axes[2].grid(True, linestyle='--', alpha=0.7)
axes[2].set_xlabel('Time [s]', fontsize=11)
axes[2].set_ylabel('Angular acceleration [deg/s^2]', fontsize=11)
axes[2].set_title('Main pole angular acceleration', fontsize=12)
axes[2].set_xlim(0, 3)

fig.suptitle('Main pole angular kinematics around the vertical axis', fontsize=13)
plt.tight_layout()

plt.savefig('../plot/pole_rot.png', dpi=300, bbox_inches='tight')
#plt.show()


print("Main pole plotting completed successfully.")

print("Plotting Hinge_Pend_i, Cadran1_nacelle_i and Cadran2_nacelle_i for i = 1, 2, 3, 4.")
fig, axes = plt.subplots(4, 3, figsize=(15, 12))

nacelle_names = ["Hinge_Pend_i", "Cadran1_nacelle_i", "Cadran2_nacelle_i"]

for i in range(1, 5):
    for j, joint_name in enumerate(nacelle_names):
        actual_name = joint_name.replace("_i", f"_{i}")
        joint_id = mbs_data.joint_id[actual_name]
    
        t = results.q[:, 0]
        q = results.q[:, joint_id]
        
        ax = axes[i-1, j]
        ax.plot(t, q, color='darkblue')
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.set_xlabel('Time [s]', fontsize=10)
        ax.set_ylabel('Coordinate', fontsize=10)
        ax.set_title(f'{actual_name}', fontsize=11)
        ax.set_ylim(bottom=-0.9, top=0.9)

plt.tight_layout()
plt.savefig('../plot/nacelles.png', dpi=300, bbox_inches='tight')
#plt.show()

print("Nacelles plotting completed successfully.")

# %%
