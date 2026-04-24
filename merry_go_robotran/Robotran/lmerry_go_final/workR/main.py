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

def save_reference_data(results, mbs_data, filepath: str) -> None:
    """
    Extract and save all arrays needed for plot_pole and plot_nacelle_pendulum_grid
    into a single .txt file with a named header.

    Columns
    -------
    t                          : time vector [s]
    pole_theta, pole_omega,
    pole_alpha                 : pole kinematics [deg, deg/s, deg/s²]
    hinge_1..4                 : pendulum hinge angles [rad]  (one per sub-system)
    cadran1_1..4               : nacelle Cardan-1 angles [rad]
    cadran2_1..4               : nacelle Cardan-2 angles [rad]
    """
    t = results.q[:, 0]

    # --- Pole kinematics (converted to degrees) ---
    id_pole     = mbs_data.joint_id["Pole_rotation"]
    pole_theta  = np.degrees(results.q  [:, id_pole])
    pole_omega  = np.degrees(results.qd [:, id_pole])
    pole_alpha  = np.degrees(results.qdd[:, id_pole])

    # --- Nacelle / pendulum angles (kept in radians, plots.py converts) ---
    hinge_cols   = []
    cadran1_cols = []
    cadran2_cols = []
    for i in range(1, 5):
        if i == 2 or i ==3:
            hinge_cols.append(-results.q[:, mbs_data.joint_id[f"Hinge_Pend_{i}"  ]])
            cadran1_cols.append(-results.q[:, mbs_data.joint_id[f"Cadran1_nacelle_{i}"]])

        else:
            hinge_cols  .append(results.q[:, mbs_data.joint_id[f"Hinge_Pend_{i}"  ]])
            cadran1_cols.append(results.q[:, mbs_data.joint_id[f"Cadran1_nacelle_{i}"]])
        if i == 3 or i ==4:
            cadran2_cols.append(-results.q[:, mbs_data.joint_id[f"Cadran2_nacelle_{i}"]])
        else:
            cadran2_cols.append(results.q[:, mbs_data.joint_id[f"Cadran2_nacelle_{i}"]])

    # --- Assemble into one (n_steps × 13) array ---
    data = np.column_stack([
        t,
        pole_theta, pole_omega, pole_alpha,
        *hinge_cols,
        *cadran1_cols,
        *cadran2_cols,
    ])

    header = (
        "t  "
        "pole_theta_deg  pole_omega_deg_s  pole_alpha_deg_s2  "
        "hinge_1  hinge_2  hinge_3  hinge_4  "
        "cadran1_1  cadran1_2  cadran1_3  cadran1_4  "
        "cadran2_1  cadran2_2  cadran2_3  cadran2_4"
    )
    np.savetxt(filepath, data, header=header, comments='# ', fmt='%.8e')
    print(f"[reference] Saved reference data → {filepath}")


# Call it right after the existing plots, pointing to your preferred path
save_reference_data(results, mbs_data, '../plot/reference_data.txt')

# %%
