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

print("Plotting results...")
# Récupération de l'ID du joint du pôle principal (celui qui tourne)
# /!\ Modifie "R3_pole_main" par le nom exact que tu as donné à ce joint dans MBsysPad
id_pole = mbs_data.joint_id["Pole_rotation"]

# Dans les résultats de Robotran (results.qd), la colonne 0 est le temps,
# et la colonne correspondant à l'ID du joint contient sa vitesse.
t = results.qd[:, 0]
omega_rad_s = results.qd[:, id_pole]

# Conversion de rad/s en deg/s
omega_deg_s = omega_rad_s * (180.0 / np.pi)

# Création de la figure (Graphe Go/No-Go)
fig = plt.figure(figsize=(10, 6))
axis = fig.gca()

axis.plot(t, omega_deg_s, label='Main pole angular velocity', color='darkblue')

# Mise en forme pour ressembler au graphe du MOOC
axis.grid(True, linestyle='--', alpha=0.7)
axis.set_xlim(left=0.0, right=3.0)
axis.set_ylim(bottom=-50.0, top=125.0)
axis.set_xlabel('Time [s]', fontsize=12)
axis.set_ylabel('Angular velocity [deg/s]', fontsize=12)
axis.set_title('Angular velocity around the vertical axis of the main pole as a function of time', fontsize=13)
axis.legend(loc='upper left', fontsize=11)

plt.show()
print("Plotting completed successfully.")