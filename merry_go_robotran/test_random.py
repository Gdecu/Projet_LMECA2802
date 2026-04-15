import numpy as np
from data.carousel_data import build_bodies
from neri.joint import make_all_joints
from neri.state import MBState
from neri.forward import forward_pass
from neri.backward import backward_pass
from neri.assembly import assemble_Mc
from simulation.driven_vars import inject_driven_into_state, evaluate_driven

# --- 1. Build the system ---
bodies = build_bodies()
joints = make_all_joints(bodies)
n_dof  = sum(b.n_dof for b in bodies)
state  = MBState(n_bodies=len(bodies) + 1, n_dof=n_dof)

# --- 2. Set a specific configuration ---
q  = np.zeros(n_dof)
qd = np.zeros(n_dof)
qdd = np.zeros(n_dof)
# Example: pole spinning at 0.5 rad/s
tilt1 = next(b for b in bodies if b.name == "pole_tilt1")
tilt2 = next(b for b in bodies if b.name == "pole_tilt2")
pendulum_idx = [4, 6, 8, 10]
pend_q = -75 * np.pi / 180
q[pendulum_idx] = pend_q
#qd[pole.q_indices[0]] = 0.5  # ω_pole = 0.5 rad/s
qd1 = np.pi**2 / 180 * np.sin(np.pi / 2)
qd[tilt1.q_indices[0]] = qd1
qd2 = np.pi**2 / 180 * np.sin(0)
qd[tilt2.q_indices[0]] = qd2
qdd1 = np.pi**3 / 450 * np.cos(np.pi / 2)
qdd[tilt1.q_indices[0]] = qdd1
qdd2 = np.pi**3 / 450 * np.cos(0)
qdd[tilt2.q_indices[0]] = qdd2
state.set_q(q, qd, t=0.0)
t = 0.0
driven_vals = evaluate_driven(bodies, t)
print(driven_vals)
inject_driven_into_state(state, driven_vals)
# --- 3. Run forward pass ---
forward_pass(bodies, joints, state)

# --- 4. Inspect per-body quantities ---
for body in bodies:
    continue
    i = body.body_id
    print(f"\n{'='*50}")
    print(f"Body {i}: {body.name}")
    print(f"  R[{i}]      =\n{state.R[i]}")
    print(f"  omega[{i}]  = {state.omega[i]}")
    print(f"  omegad[{i}] = {state.omegad[i]}")
    print(f"  beta[{i}]   =\n{state.beta[i]}")
    print(f"  alpha[{i}]  = {state.alpha[i]}")
    print(f"  d_hi[{i}]   = {state.d_hi[i]}")
    print(f"  d_ii[{i}]   = {state.d_ii[i]}")

# --- 5. Inspect the mass matrix ---
backward_pass(bodies, joints, state)
#M, c = assemble_Mc(bodies, joints, state)
np.set_printoptions(precision=1, suppress=True)
#print(f"\nM (mass matrix):\n{M}")
#print(f"\nc (bias vector):\n{c}")
#print(f"\nM symmetry check: {np.allclose(M, M.T)}")
#print(f"M eigenvalues: {np.linalg.eigvalsh(M)}")