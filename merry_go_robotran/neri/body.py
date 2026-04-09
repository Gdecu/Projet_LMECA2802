"""
neri/body.py
------------
Defines the Body dataclass representing one rigid body in the multibody tree.

Each Body stores:
  - its identity and position in the tree (parent link)
  - its inertial properties (mass, inertia tensor at CoM)
  - its constant geometry vectors (in local body frames)
  - its joint metadata (type, DOF count, variable classification)

These are the FIXED properties of a body — nothing here changes during
simulation. All time-varying quantities (ω, α, F, L, R) live in MBState.

Notation (following NERi slides):
  O^i   : origin of joint i (attachment point of body i to its parent)
  G^i   : center of mass of body i
  d^{hi}: vector from O^h (parent joint) to O^i (this joint), in PARENT frame
  d^{ii}: vector from O^i (this joint)  to G^i (CoM),        in THIS body frame
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


@dataclass
class Body:
    """
    One rigid body in an open-chain (tree) multibody system.

    Parameters
    ----------
    body_id : int
        Unique body index. 0 is reserved for the fixed base (ground).
        Bodies are numbered 1 … N_body in topological order (parent before child).
    name : str
        Human-readable label (e.g. 'pole', 'pendulum_1', 'nacelle_3').
    parent_id : int
        Index of the parent body. 0 means the body connects directly to ground.

    mass : float
        Total mass of the body [kg].
    inertia_com_local : np.ndarray, shape (3, 3)
        Inertia tensor at the CoM, expressed in the body-fixed frame [kg·m²].
        Must be symmetric and positive semi-definite.
        For bodies with axial symmetry (pole, pendulum, nacelle):
            diag([I_transverse, I_transverse, I_axial])

    d_parent_to_joint_in_parent : np.ndarray, shape (3,)
        Vector from the PARENT joint origin O^h to THIS joint origin O^i,
        expressed in the PARENT body frame [m].
        This is the constant "arm" vector that locates where this body
        attaches to its parent. It gets rotated to the inertial frame
        during the forward kinematic pass via:  d^{hi} = R^{0h} @ d_parent_to_joint_in_parent
    d_joint_to_com_in_local : np.ndarray, shape (3,)
        Vector from THIS joint origin O^i to the CoM G^i,
        expressed in THIS body's local frame [m].
        Gets rotated to inertial frame during forward pass via:
            d^{ii} = R^{0i} @ d_joint_to_com_in_local

    joint_type : str
        One of: 'revolute'  — 1-DOF rotation about a fixed axis
                'cardan2'   — 2-DOF Cardan/universal joint (sequence of 2 rotations)
                'cardan3'   — 3-DOF spherical-like (main rot + 2 tilt), used for pole
                'fixed'     — 0-DOF rigid attachment (no relative motion)
    joint_axis_in_parent : np.ndarray, shape (3,) or (3, 2) for cardan2
        For 'revolute'  : unit vector of the rotation axis in the PARENT frame.
        For 'cardan2'   : shape (3, 2) — columns are the two successive rotation axes,
                          each expressed in the frame BEFORE that rotation is applied.
        For 'cardan3'   : shape (3, 3) — three successive axes (main rot, tilt1, tilt2).
        For 'fixed'     : not used, set to np.zeros(3).

    n_dof : int
        Number of degrees of freedom contributed by the joint (0, 1, 2, or 3).
    var_types : list of str, length n_dof
        Classification of each DOF:
            'u' — independent variable, integrated by Runge-Kutta
            'c' — driven/commanded variable, prescribed by control law
    q_indices : list of int, length n_dof
        Indices in the global generalised coordinate vector q for each DOF.
        Assigned once when the full system is assembled.
    """

    # --- Identity ---
    body_id:  int
    name:     str
    parent_id: int

    # --- Inertial properties ---
    mass:               float
    inertia_com_local:  np.ndarray   # (3,3), in body frame, at CoM

    # --- Geometry (constant, in local frames) ---
    d_parent_to_joint_in_parent: np.ndarray  # (3,) O^h -> O^i, in parent frame
    d_joint_to_com_in_local:     np.ndarray  # (3,) O^i -> G^i, in this body frame

    # --- Joint metadata ---
    joint_type:            str
    joint_axis_in_parent:  np.ndarray   # (3,) or (3,k) depending on joint_type
    n_dof:                 int
    var_types:             list          # ['u'], ['c'], ['u','c'], etc.
    q_indices:             list          # filled by system assembler, default empty

    # ------------------------------------------------------------------ #
    #  Post-init validation                                                #
    # ------------------------------------------------------------------ #

    def __post_init__(self):
        # Enforce numpy arrays (allows passing plain lists from data file)
        self.inertia_com_local             = np.asarray(self.inertia_com_local,            dtype=float)
        self.d_parent_to_joint_in_parent   = np.asarray(self.d_parent_to_joint_in_parent,  dtype=float)
        self.d_joint_to_com_in_local       = np.asarray(self.d_joint_to_com_in_local,      dtype=float)
        self.joint_axis_in_parent          = np.asarray(self.joint_axis_in_parent,         dtype=float)

        # Shape checks
        assert self.inertia_com_local.shape == (3, 3),  \
            f"Body '{self.name}': inertia_com_local must be (3,3), got {self.inertia_com_local.shape}"
        assert self.d_parent_to_joint_in_parent.shape == (3,), \
            f"Body '{self.name}': d_parent_to_joint_in_parent must be (3,), got {self.d_parent_to_joint_in_parent.shape}"
        assert self.d_joint_to_com_in_local.shape == (3,), \
            f"Body '{self.name}': d_joint_to_com_in_local must be (3,), got {self.d_joint_to_com_in_local.shape}"

        # Inertia symmetry check (numerical tolerance)
        assert np.allclose(self.inertia_com_local, self.inertia_com_local.T, atol=1e-10), \
            f"Body '{self.name}': inertia tensor must be symmetric."

        # var_types / n_dof consistency
        assert len(self.var_types) == self.n_dof, \
            f"Body '{self.name}': len(var_types)={len(self.var_types)} != n_dof={self.n_dof}"
        assert all(v in ('u', 'c') for v in self.var_types), \
            f"Body '{self.name}': var_types entries must be 'u' or 'c'."

    # ------------------------------------------------------------------ #
    #  Convenience constructors                                            #
    # ------------------------------------------------------------------ #

    @classmethod
    def revolute(cls,
                 body_id:   int,
                 name:      str,
                 parent_id: int,
                 mass:      float,
                 inertia_diag: list,          # [Ixx, Iyy, Izz] in body frame
                 d_parent_to_joint: list,     # 3-vector, in parent frame
                 d_joint_to_com:    list,     # 3-vector, in local frame
                 axis_in_parent:    list,     # 3-vector, rotation axis
                 var_type: str = 'u',
                 q_index:  int = -1) -> Body:
        """Convenience constructor for the common 1-DOF revolute joint case."""
        return cls(
            body_id=body_id,
            name=name,
            parent_id=parent_id,
            mass=mass,
            inertia_com_local=np.diag(inertia_diag),
            d_parent_to_joint_in_parent=np.array(d_parent_to_joint, dtype=float),
            d_joint_to_com_in_local=np.array(d_joint_to_com, dtype=float),
            joint_type='revolute',
            joint_axis_in_parent=np.array(axis_in_parent, dtype=float),
            n_dof=1,
            var_types=[var_type],
            q_indices=[q_index],
        )

    @classmethod
    def cardan2(cls,
                body_id:   int,
                name:      str,
                parent_id: int,
                mass:      float,
                inertia_diag: list,
                d_parent_to_joint: list,
                d_joint_to_com:    list,
                axes_in_parent:    list,    # (3,2) — two successive axes
                var_types: list = None,
                q_indices: list = None) -> Body:
        """Convenience constructor for a 2-DOF Cardan joint."""
        axes = np.array(axes_in_parent, dtype=float)      # shape (2,3) → transpose to (3,2)
        if axes.shape == (2, 3):
            axes = axes.T
        assert axes.shape == (3, 2), "axes_in_parent must be (2,3) or (3,2)"
        return cls(
            body_id=body_id,
            name=name,
            parent_id=parent_id,
            mass=mass,
            inertia_com_local=np.diag(inertia_diag),
            d_parent_to_joint_in_parent=np.array(d_parent_to_joint, dtype=float),
            d_joint_to_com_in_local=np.array(d_joint_to_com, dtype=float),
            joint_type='cardan2',
            joint_axis_in_parent=axes,
            n_dof=2,
            var_types=var_types or ['u', 'u'],
            q_indices=q_indices or [-1, -1],
        )

    # ------------------------------------------------------------------ #
    #  Helpers                                                             #
    # ------------------------------------------------------------------ #

    @property
    def is_leaf(self) -> bool:
        """True if this body has no children (determined externally by MBSystem)."""
        # Actual leaf status is set by the system assembler — this is a placeholder.
        return False

    def __repr__(self) -> str:
        return (f"Body(id={self.body_id}, name='{self.name}', "
                f"parent={self.parent_id}, joint='{self.joint_type}', "
                f"n_dof={self.n_dof}, vars={self.var_types})")