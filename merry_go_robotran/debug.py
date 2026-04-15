"""
debug_backward_pass.py
----------------------
Rigorous debugging of the backward pass and assembly pipeline.

STRATEGY:
  We use a STATIC configuration (all velocities zero, t=0) so that many
  terms vanish analytically and every intermediate quantity can be checked
  against hand calculations.

  Configuration at t=0:
    - Pole spin q[0] = 0
    - Tilt angles q[1]=q[2] = 0  (driven, but zero at t=0)
    - Pendulums at 75° from vertical → q_pend = 75° = 5π/12 rad
    - Nacelle Cardan angles = 0
    - ALL velocities = 0
    - ALL accelerations = 0 (we check bias only, no qdd)
    - Wind force at t=0: F = A*(1-cos(0)) = 0  → no external forces
    - No motor torque yet (ω_pole = 0 < threshold)

  Under these conditions:
    - ω^i  = 0 for all bodies     (no velocities)
    - ω̇^i  = 0 for all bodies     (no accelerations, no Coriolis)
    - β^i  = 0 for all bodies     (both skew(ω̇) and ω̃² vanish)
    - α^i  = α^0 for all bodies   (β^h · d_hi = 0 when β = 0)
    - α^0  = -g = [0, 0, 9.81]    (NERi convention: α^0 = -g)

  Therefore:
    - W^i  = m^i * α^0 - F_ext^i = m^i * [0, 0, 9.81]
      This is the gravity pseudo-force: mass × (upward g).
      For zero-mass bodies (pole body 1, pole_tilt1 body 2): W = 0.
    - Euler term = I·ω̇ + ω̃·I·ω = 0  (both ω and ω̇ are zero)
    - F^i  = W^i + Σ F^j           (pure accumulation of weight)
    - L^i  = skew(d_ii)·W^i + Σ(L^j + skew(d_hi^j)·F^j)
      Pure static moment balance (no Euler, no external torques).

WHAT TO LOOK FOR:
  Any deviation from the above means a bug in:
    - The sign of α^0 in forward.py  (g vs -g, critically important)
    - The children map construction in backward.py
    - The traversal order (must be leaves → root)
    - The transport term d_hi[j] in the moment equation
    - The projection Q = L · φ
    - The mass matrix assembly (last trick)
"""

import numpy as np
import sys
import os

# ──────────────────────────────────────────────────────────────────────── #
#  Adjust import path — adapt if your project layout differs              #
# ──────────────────────────────────────────────────────────────────────── #
# sys.path.insert(0, os.path.abspath(".."))  # uncomment if needed

from merry_go_robotran.data.carousel_data import build_bodies
from merry_go_robotran.neri.joint import make_all_joints
from merry_go_robotran.neri.state import MBState
from merry_go_robotran.neri.forward import forward_pass
from merry_go_robotran.neri.backward import backward_pass
from merry_go_robotran.neri.assembly import assemble_Mc, get_dof_indices, \
    solve_independent_accelerations
from merry_go_robotran.simulation.driven_vars import evaluate_driven, \
    inject_driven_into_state, get_qddc_vector

np.set_printoptions(precision=6, suppress=True, linewidth=120)

# ──────────────────────────────────────────────────────────────────────── #
#  Helper: separator for readability                                      #
# ──────────────────────────────────────────────────────────────────────── #
def section(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")

def check(name, value, expected, tol=1e-6):
    """Compare a value to its expected result and flag mismatches."""
    val = np.asarray(value)
    exp = np.asarray(expected)
    ok  = np.allclose(val, exp, atol=tol)
    tag = "  OK  " if ok else "**FAIL**"
    print(f"  [{tag}] {name}")
    print(f"           got:      {val}")
    print(f"           expected: {exp}")
    if not ok:
        print(f"           diff:     {val - exp}")
    return ok

# ====================================================================== #
#  1. BUILD SYSTEM                                                        #
# ====================================================================== #
section("1. SYSTEM ASSEMBLY")

bodies = build_bodies()
joints = make_all_joints(bodies)
n_dof  = sum(b.n_dof for b in bodies)
n_bodies = len(bodies) + 1          # +1 for base body 0

state = MBState(n_bodies=n_bodies, n_dof=n_dof)

print(f"  Number of bodies (incl. base): {n_bodies}")
print(f"  Number of DOFs:                {n_dof}")

# Print body table for reference
print(f"\n  {'ID':>3} {'Name':<16} {'Parent':>6} {'Mass':>8} {'Joint':<10} "
      f"{'n_dof':>5} {'q_idx'}")
for b in bodies:
    print(f"  {b.body_id:>3} {b.name:<16} {b.parent_id:>6} {b.mass:>8.1f} "
          f"{b.joint_type:<10} {b.n_dof:>5} {b.q_indices}")

# ====================================================================== #
#  2. SET UP STATIC CONFIGURATION (t=0, all velocities zero)              #
# ====================================================================== #
section("2. STATIC CONFIGURATION SETUP")

t = 0.0
q  = np.zeros(n_dof)
qd = np.zeros(n_dof)

# Inject driven values at t=0 (tilt angles should be zero at t=0)
#driven_vals = evaluate_driven(bodies, t)
print(f"\n  Driven values at t={t}:")
#for q_idx, (qc, qdc, qddc) in driven_vals.items():
 #   print(f"    q[{q_idx}]:  q={qc:.6f}  qd={qdc:.6f}  qdd={qddc:.6f}")

state.set_q(q, qd, t)
#inject_driven_into_state(state, driven_vals)

# Set pendulum angles to 75° (initial condition from MOOC)
# EXPECT: the pendulum DOF indices (one per pendulum).
#         75° = 5π/12 ≈ 1.3090 rad
theta_pend = - 15.0 * np.pi / 180.0
for b in bodies:
    if b.name.startswith("pendulum"):
        for qi in b.q_indices:
            state.q[qi] = theta_pend
            print(f"  Set q[{qi}] ({b.name}) = {theta_pend:.6f} rad ({75.0}°)")

print(f"\n  Full q  = {state.q}")
print(f"  Full qd = {state.qd}")

# ====================================================================== #
#  3. FORWARD PASS — VERIFY STATIC CONDITIONS                            #
# ====================================================================== #
section("3. FORWARD PASS CHECKS")

g_vec = np.array([0.0, 0.0, -9.81])
forward_pass(bodies, joints, state, g_vec)

# ── CHECK 3a: α^0 should be -g = [0, 0, +9.81] ────────────────────── #
# The NERi convention is α^0 = -g (see slide 6-7 of NERi_V1.pdf).
# If g = [0,0,-9.81], then α^0 = [0,0,+9.81].
# THIS IS THE MOST COMMON SIGN BUG. If α^0 = g = [0,0,-9.81],
# gravity will act UPWARD in the equations of motion.
print("\n  --- Check 3a: Sign of alpha[0] ---")
check("alpha[0] = -g (must be [0, 0, +9.81])",
      state.alpha[0], [0.0, 0.0, 9.81])

# ── CHECK 3b: All ω should be zero (static case) ───────────────────── #
print("\n  --- Check 3b: All omega = 0 (static) ---")
all_omega_zero = True
for b in bodies:
    i = b.body_id
    if not np.allclose(state.omega[i], 0.0, atol=1e-12):
        print(f"  **FAIL** omega[{i}] ({b.name}) = {state.omega[i]}  (should be 0)")
        all_omega_zero = False
if all_omega_zero:
    print("  [  OK  ] All omega are zero")

# ── CHECK 3c: All β should be zero (no velocity → no centripetal) ──── #
print("\n  --- Check 3c: All beta = 0 (static) ---")
all_beta_zero = True
for b in bodies:
    i = b.body_id
    if not np.allclose(state.beta[i], 0.0, atol=1e-12):
        print(f"  **FAIL** beta[{i}] ({b.name}) =\n{state.beta[i]}")
        all_beta_zero = False
if all_beta_zero:
    print("  [  OK  ] All beta are zero")

# ── CHECK 3d: All α should equal α^0 (β=0 → no transport effect) ──── #
# Since α^i = α^h + β^h · d_hi and β^h = 0, we get α^i = α^0 everywhere.
print("\n  --- Check 3d: All alpha = alpha[0] (no transport when beta=0) ---")
for b in bodies:
    i = b.body_id
    check(f"alpha[{i}] ({b.name})", state.alpha[i], state.alpha[0])

# ── CHECK 3e: Geometry vectors d_hi and d_ii ─────────────────────── #
print("\n  --- Check 3e: Geometry vectors ---")
for b in bodies:
    i = b.body_id
    print(f"  Body {i:>2} ({b.name:<16}): "
          f"d_hi={state.d_hi[i]}  d_ii={state.d_ii[i]}")
    # EXPECT for zero-mass virtual bodies (pole, tilt1, tilt2 if mass=0):
    #   d_ii = 0 (no CoM offset when d_joint_to_com = 0)
    # EXPECT for pendulums: d_ii should point ~downward (along -I3
    #   when vertical), but rotated by 75° around hinge axis.
    # EXPECT for nacelles: d_ii ≈ R^0i · [0, 0, -1]

# ====================================================================== #
#  4. BACKWARD PASS — DETAILED CHECK                                      #
# ====================================================================== #
section("4. BACKWARD PASS CHECKS")

# No external forces at t=0 (wind = 0)
F_ext = np.zeros((n_bodies, 3))
L_ext = np.zeros((n_bodies, 3))

backward_pass(bodies, joints, state, F_ext, L_ext)

# ── CHECK 4a: W^i = m^i * alpha (gravity pseudo-force) ──────────── #
# Since α^i = α^0 and β^i = 0, we have:
#   W^i = m^i * (α + β·d_ii) - F_ext = m^i * α^0
# For zero-mass bodies: W = 0.
# For massive bodies: W should point in +I3 direction (upward, m*9.81).
print("\n  --- Check 4a: W^i = m * alpha[0] (static gravity pseudo-force) ---")
for b in bodies:
    i = b.body_id
    expected_W = b.mass * state.alpha[0]
    check(f"W[{i}] ({b.name}, m={b.mass:.1f}kg)", state.W[i], expected_W)

# ── CHECK 4b: F^i accumulation (leaves → root) ──────────────────── #
# Build children map to manually verify accumulation
print("\n  --- Check 4b: F^i = W^i + Σ F^j (force accumulation) ---")
children = {b.body_id: [] for b in bodies}
children[0] = []
for b in bodies:
    children[b.parent_id].append(b.body_id)
print(f"  Children map: {children}")

# For leaf bodies (nacelles): F = W (no children)
for b in bodies:
    i = b.body_id
    if len(children.get(i, [])) == 0:
        # Leaf body: F must equal W
        check(f"F[{i}] ({b.name}, leaf) = W[{i}]", state.F[i], state.W[i])

# For all bodies: verify F[i] = W[i] + sum of F[j] for j in children
for b in reversed(bodies):
    i = b.body_id
    expected_F = state.W[i].copy()
    for j in children.get(i, []):
        expected_F += state.F[j]
    check(f"F[{i}] ({b.name}) = W + ΣF_children", state.F[i], expected_F)

# ── CHECK 4c: Total force at root should equal total weight ─────── #
# The root body (body with parent=0) carries the total system weight.
# F_root = total_mass * g_up = total_mass * [0, 0, 9.81]
print("\n  --- Check 4c: Root force = total weight ---")
total_mass = sum(b.mass for b in bodies)
print(f"  Total system mass: {total_mass:.1f} kg")
# Find root bodies (parent_id=0)
root_F_sum = np.zeros(3)
for b in bodies:
    if b.parent_id == 0:
        root_F_sum += state.F[b.body_id]
expected_total_F = total_mass * state.alpha[0]
check("Sum of root F = total_mass * alpha[0]", root_F_sum, expected_total_F)

# ── CHECK 4d: L^i moment equation ───────────────────────────────── #
# In the static case (ω=ω̇=0), the Euler term vanishes:
#   euler = I·ω̇ + ω̃·I·ω = 0
# So L^i = skew(d_ii)·W^i + Σ(L^j + skew(d_hi^j)·F^j) - L_ext
# This is a pure static moment balance.
print("\n  --- Check 4d: L^i moment (static, no Euler term) ---")
def skew(v):
    return np.array([[ 0, -v[2], v[1]],
                     [ v[2], 0, -v[0]],
                     [-v[1], v[0], 0 ]])

for b in reversed(bodies):
    i = b.body_id
    # Manually compute expected L
    expected_L = skew(state.d_ii[i]) @ state.W[i]  # moment from own weight
    # Euler term should be zero — verify
    R0i = state.R[i]
    I_world = R0i @ b.inertia_com_local @ R0i.T
    euler = I_world @ state.omegad[i] + skew(state.omega[i]) @ (I_world @ state.omega[i])
    expected_L += euler   # should add zero in static case
    expected_L -= L_ext[i]
    # Children contributions
    for j in children.get(i, []):
        expected_L += state.L[j] + skew(state.d_hi[j]) @ state.F[j]
    check(f"L[{i}] ({b.name})", state.L[i], expected_L)

# ── CHECK 4e: Q projection ──────────────────────────────────────── #
# Q[q_idx] = L^i · φ_k
# This is the generalized torque at each joint DOF.
# For pendulums in static equilibrium at 75°, Q should represent the
# gravitational torque trying to swing them back to vertical.
print("\n  --- Check 4e: Q = L · phi (generalized force projection) ---")
for b in bodies:
    i = b.body_id
    joint = joints[bodies.index(b)]
    if joint.n_dof == 0:
        continue
    q_i   = state.q[b.q_indices]
    phi_i = joint.axes_inertial(q_i, state.R[b.parent_id])
    for k, q_idx in enumerate(b.q_indices):
        expected_Q = np.dot(state.L[i], phi_i[:, k])
        check(f"Q[{q_idx}] ({b.name}, dof {k})", state.Q[q_idx], expected_Q)

# ====================================================================== #
#  5. MASS MATRIX ASSEMBLY — STRUCTURAL CHECKS                           #
# ====================================================================== #
section("5. MASS MATRIX CHECKS")

M, c_vec = assemble_Mc(bodies, joints, state, F_ext, L_ext, g_vec)

print(f"\n  M shape: {M.shape}")
print(f"  c shape: {c_vec.shape}")

# ── CHECK 5a: M must be symmetric ─────────────────────────────────── #
sym_err = np.max(np.abs(M - M.T))
check("M is symmetric (max|M - M^T|)", sym_err, 0.0, tol=1e-10)

# ── CHECK 5b: M must be positive semi-definite ────────────────────── #
# Eigenvalues of M should all be ≥ 0. Zero eigenvalues correspond
# to DOFs of zero-mass bodies (pole, tilt1 if mass=0).
eigvals = np.linalg.eigvalsh(M)
print(f"\n  M eigenvalues: {eigvals}")
n_negative = np.sum(eigvals < -1e-10)
n_zero     = np.sum(np.abs(eigvals) < 1e-10)
print(f"  Negative eigenvalues: {n_negative}  (EXPECT: 0)")
print(f"  Zero eigenvalues:     {n_zero}  (EXPECT: = number of zero-mass DOFs)")
if n_negative > 0:
    print("  **FAIL** Negative eigenvalues → M is not positive semi-definite!")
    print("           This means the backward pass or assembly is wrong.")

# ── CHECK 5c: Diagonal of M (rough sanity) ─────────────────────── #
# M[k,k] ≥ 0 always. For pendulum hinges: M_kk should be on the order
# of m_pend * L² + I_pend ≈ 44 * 1.5² + 15.9 ≈ 115 kg·m² (very rough).
print(f"\n  M diagonal: {np.diag(M)}")
print("  EXPECT: all ≥ 0. Pendulum DOFs ~ O(100) kg·m².")
print("          Zero-mass DOFs (tilt1, tilt2 if massless) may be small.")

# ── CHECK 5d: Print full M for inspection ──────────────────────── #
print(f"\n  Full mass matrix M:")
print(M)

# ── CHECK 5e: c vector (bias = gravity + Coriolis at qd=0) ──────── #
# With qd=0 and F_ext=0: c should contain ONLY gravitational torques.
# For pendulums at 75°, c should show a restoring torque from gravity.
print(f"\n  Bias vector c: {c_vec}")
print("  EXPECT: nonzero only where gravity creates torque on a DOF.")
print("          Pendulum DOFs should have significant gravity torque.")
print("          DOFs of zero-mass bodies should have zero bias.")

# ====================================================================== #
#  6. DOF PARTITIONING CHECK                                              #
# ====================================================================== #
section("6. DOF PARTITION (u vs c)")

idx_u, idx_c = get_dof_indices(bodies)
print(f"  Independent (u) DOF indices: {idx_u}")
print(f"  Driven (c) DOF indices:      {idx_c}")

# Verify var_types match
for b in bodies:
    for vt, qi in zip(b.var_types, b.q_indices):
        if vt == 'u' and qi not in idx_u:
            print(f"  **FAIL** {b.name} DOF {qi} is 'u' but not in idx_u!")
        if vt == 'c' and qi not in idx_c:
            print(f"  **FAIL** {b.name} DOF {qi} is 'c' but not in idx_c!")

# ── CHECK 6a: Muu sub-matrix ──────────────────────────────────────── #
Muu = M[np.ix_(idx_u, idx_u)]
Muc = M[np.ix_(idx_u, idx_c)]
print(f"\n  Muu ({len(idx_u)}×{len(idx_u)}):")
print(Muu)
print(f"\n  Muc ({len(idx_u)}×{len(idx_c)}):")
print(Muc)

# Muu must be invertible (positive definite) for the solver to work
eigvals_uu = np.linalg.eigvalsh(Muu)
print(f"\n  Muu eigenvalues: {eigvals_uu}")
if np.any(eigvals_uu < 1e-12):
    print("  **FAIL** Muu is singular or near-singular!")
    print("           The ODE solver CANNOT invert this → divergence.")
else:
    print("  [  OK  ] Muu is positive definite (invertible)")

# ====================================================================== #
#  7. SOLVE FOR ACCELERATIONS (static case)                               #
# ====================================================================== #
section("7. SOLVE FOR qdd_u (static, no external forces)")

qdd_c = get_qddc_vector(bodies, t, n_dof)
Q_ext = np.zeros(n_dof)

qdd_u = solve_independent_accelerations(M, c_vec, idx_u, idx_c,
                                         qdd_c[idx_c], Q_ext)

print(f"  qdd_u = {qdd_u}")
print("\n  EXPECT for pendulums:")
print("    Gravity pulls pendulums toward vertical → negative angular")
print("    acceleration (restoring torque). Magnitude should be physically")
print("    reasonable: |qdd| ~ g*sin(75°)/L ~ 9.81*0.966/1.5 ≈ 6.3 rad/s²")
print("    (rough single-pendulum estimate, actual value differs due to")
print("    coupled inertias and nacelle mass).")

for i, qi in enumerate(idx_u):
    body = next(b for b in bodies if qi in b.q_indices)
    print(f"    qdd[{qi}] ({body.name}): {qdd_u[i]:.6f} rad/s²")

# ── Sanity: accelerations should be finite and reasonable ──────── #
if np.any(np.isnan(qdd_u)) or np.any(np.isinf(qdd_u)):
    print("\n  **FAIL** NaN or Inf in accelerations → solver will diverge!")
if np.any(np.abs(qdd_u) > 1000):
    print("\n  **WARNING** Very large accelerations detected (>1000 rad/s²).")
    print("             This suggests a mass matrix or bias vector problem.")

# ====================================================================== #
#  8. FULL ODE EVALUATION (mimics what solve_ivp calls)                   #
# ====================================================================== #
section("8. FULL ODE EVALUATION at t=0")

from merry_go_robotran.simulation.integrator import make_ode

state2 = MBState(n_bodies=n_bodies, n_dof=n_dof)
f_ode, ctx = make_ode(bodies, joints, state2)

# Build y0: [q_u | qd_u] at t=0
q0_full  = np.zeros(n_dof)
qd0_full = np.zeros(n_dof)
for b in bodies:
    if b.name.startswith("pendulum"):
        for qi in b.q_indices:
            q0_full[qi] = theta_pend

y0 = np.concatenate([q0_full[idx_u], qd0_full[idx_u]])

print(f"  y0 = {y0}")
print(f"  Calling f_ode(t=0, y0)...")

try:
    yd = f_ode(0.0, y0)
    print(f"  yd = {yd}")
    print(f"\n  yd[:nu]  = qd_u  = {yd[:len(idx_u)]}")
    print(f"  yd[nu:]  = qdd_u = {yd[len(idx_u):]}")

    if np.any(np.isnan(yd)) or np.any(np.isinf(yd)):
        print("\n  **FAIL** ODE returns NaN/Inf → solver WILL diverge!")
    elif np.any(np.abs(yd) > 1e6):
        print("\n  **FAIL** ODE returns extremely large values → solver will struggle!")
    else:
        print("\n  [  OK  ] ODE evaluation produced finite, reasonable values.")

except Exception as e:
    print(f"\n  **FAIL** ODE evaluation crashed: {e}")
    import traceback
    traceback.print_exc()

# ====================================================================== #
#  9. SECOND EVALUATION at t=10 (nonzero driven angles + wind)            #
# ====================================================================== #
section("9. ODE EVALUATION at t=10 (with driven motion + wind)")

try:
    yd10 = f_ode(10.0, y0)
    print(f"  yd at t=10: {yd10}")
    if np.any(np.isnan(yd10)) or np.any(np.isinf(yd10)):
        print("  **FAIL** NaN/Inf at t=10!")
    elif np.any(np.abs(yd10) > 1e6):
        print("  **WARNING** Very large derivatives at t=10.")
    else:
        print("  [  OK  ] t=10 evaluation looks reasonable.")
except Exception as e:
    print(f"  **FAIL** Crashed at t=10: {e}")
    import traceback
    traceback.print_exc()

# ====================================================================== #
#  SUMMARY                                                                #
# ====================================================================== #
section("DEBUGGING COMPLETE")
print("""
  Key things to investigate if failures appeared above:

  1. SIGN OF alpha[0]:
     forward.py must set  state.alpha[0] = -g  (NOT +g).
     If g = [0,0,-9.81], then alpha[0] = [0,0,+9.81].
     Getting this wrong flips gravity direction in the entire simulation.

  2. ZERO-MASS BODIES:
     Bodies 1 and 2 (pole, pole_tilt1) have mass=0 and I=0.
     Their W must be exactly zero. If not, something is injecting
     spurious forces.

  3. CHILDREN MAP:
     The backward pass builds a children dict from body.parent_id.
     If bodies are not in topological order, or parent IDs are wrong,
     force accumulation will be incorrect.

  4. d_hi[j] vs d_hi[i] IN MOMENT TRANSPORT:
     The line  L[i] += L[j] + skew(d_hi[j]) @ F[j]  uses d_hi[j],
     which is the vector from O^i to O^j. Using d_hi[i] here would
     be wrong (that's the vector from O^h to O^i).

  5. Muu SINGULARITY:
     If Muu has zero or negative eigenvalues, the linear solve in the
     ODE will fail or produce garbage. This typically means a body in
     the independent partition has zero mass AND zero children inertia.

  6. NACELLE CARDAN AXES:
     The axes_in_parent for nacelle Cardan2 joints must match the
     physical configuration. First axis = hinge axis (tangential),
     second axis = arm direction (radial). If swapped, the nacelle
     rotation is about the wrong axes.
""")