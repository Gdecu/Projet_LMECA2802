"""
neri/joint.py
-------------
Joint classes for the NERi formalism.

Each joint knows how to:
  1. Build its relative rotation matrix R^{hi}(q)
  2. Return its joint axis φ^i in the inertial frame
  3. Return ψ^i = φ^i × d^{ii} (linear velocity jacobian)

All joints operate on numpy arrays. No state is stored here —
joints are pure functions of (q, R0h).
"""

import numpy as np
from abc import ABC, abstractmethod


# ------------------------------------------------------------------ #
#  Utility: skew-symmetric (cross-product) matrix                     #
# ------------------------------------------------------------------ #

def skew(v: np.ndarray) -> np.ndarray:
    """Return the 3×3 skew-symmetric matrix of vector v such that
    skew(v) @ u = v × u."""
    return np.array([
        [ 0.0,  -v[2],  v[1]],
        [ v[2],  0.0,  -v[0]],
        [-v[1],  v[0],  0.0 ]
    ])


def rodrigues(axis: np.ndarray, angle: float) -> np.ndarray:
    """
    Rotation matrix for a rotation of `angle` [rad] about unit `axis`,
    using Rodrigues' formula:  R = I + sin(θ)ñ + (1-cos(θ))ñ²
    """
    n  = axis / np.linalg.norm(axis)   # ensure unit vector
    N  = skew(n)
    return np.eye(3) + np.sin(angle)*N + (1 - np.cos(angle))*(N @ N)


# ------------------------------------------------------------------ #
#  Abstract base                                                       #
# ------------------------------------------------------------------ #

class Joint(ABC):
    """Base class for all joint types."""

    @property
    @abstractmethod
    def n_dof(self) -> int:
        """Number of degrees of freedom."""

    @abstractmethod
    def rotation_matrix(self, q: np.ndarray) -> np.ndarray:
        """
        Relative rotation matrix R^{hi}(q).
        Maps vectors from THIS body's frame to PARENT's frame.

        Parameters
        ----------
        q : np.ndarray, shape (n_dof,)
            Current joint coordinates.

        Returns
        -------
        R : np.ndarray, shape (3, 3)
        """

    @abstractmethod
    def axes_inertial(self, q: np.ndarray,
                      R0h: np.ndarray) -> np.ndarray:
        """
        Joint axes expressed in the INERTIAL frame.

        Parameters
        ----------
        q   : current joint coordinates, shape (n_dof,)
        R0h : rotation matrix of the PARENT body (inertial ← parent),
              shape (3, 3)

        Returns
        -------
        phi : np.ndarray, shape (3, n_dof)
              Each column is one joint axis in inertial frame.
              For a 1-DOF joint this is shape (3, 1).
        """

    def psi(self, phi: np.ndarray, d_ii: np.ndarray) -> np.ndarray:
        """
        Linear velocity jacobian ψ^i = φ^i × d^{ii}.

        For multi-DOF joints: each column k → φ_k × d_ii.

        Parameters
        ----------
        phi  : joint axes in inertial frame, shape (3, n_dof)
        d_ii : vector O^i → G^i in inertial frame, shape (3,)

        Returns
        -------
        psi : np.ndarray, shape (3, n_dof)
        """
        # phi is (3, n_dof); compute cross product column-wise
        psi = np.zeros_like(phi)
        for k in range(self.n_dof):
            psi[:, k] = np.cross(phi[:, k], d_ii)
        return psi
        #return np.zeros([3,3])#[0, 0, 0]


# ------------------------------------------------------------------ #
#  Fixed joint  (0 DOF)                                               #
# ------------------------------------------------------------------ #

class FixedJoint(Joint):
    """Rigid attachment — no relative motion between bodies."""

    @property
    def n_dof(self) -> int:
        return 0

    def rotation_matrix(self, q: np.ndarray) -> np.ndarray:
        return np.eye(3)

    def axes_inertial(self, q: np.ndarray,
                      R0h: np.ndarray) -> np.ndarray:
        return np.zeros((3, 0))   # no axes


# ------------------------------------------------------------------ #
#  Revolute joint  (1 DOF)                                            #
# ------------------------------------------------------------------ #

class RevoluteJoint(Joint):
    """
    1-DOF rotation about a fixed axis.
    Used for the 4 pendulum hinges.

    Parameters
    ----------
    axis_in_parent : np.ndarray, shape (3,)
        Rotation axis expressed in the PARENT body frame.
    """

    def __init__(self, axis_in_parent: np.ndarray):
        ax = np.asarray(axis_in_parent, dtype=float)
        self._axis = ax / np.linalg.norm(ax)    # normalise once

    @property
    def n_dof(self) -> int:
        return 1

    def rotation_matrix(self, q: np.ndarray) -> np.ndarray:
        """R^{hi}(q[0]) — rotation about the stored axis by angle q[0]."""
        return rodrigues(self._axis, q[0])

    def axes_inertial(self, q: np.ndarray,
                      R0h: np.ndarray) -> np.ndarray:
        """
        φ^i = R^{0h} @ axis_in_parent
        The joint axis expressed in the inertial frame.
        Returns shape (3, 1).
        """
        phi = R0h @ self._axis          # rotate axis to inertial frame
        return phi.reshape(3, 1)


# ------------------------------------------------------------------ #
#  Cardan2 joint  (2 DOF)                                             #
# ------------------------------------------------------------------ #

class Cardan2Joint(Joint):
    """
    2-DOF Cardan (universal) joint — two successive rotations.
    Used for the 4 nacelle joints.

    Rotation sequence: first about axis_1 (in parent frame) by q[0],
    then about axis_2 (in the intermediate frame after first rotation)
    by q[1].

    Parameters
    ----------
    axes_in_parent : np.ndarray, shape (3, 2)
        Column 0 = first rotation axis in parent frame.
        Column 1 = second rotation axis in parent frame (at q=0).
    """

    def __init__(self, axes_in_parent: np.ndarray):
        ax = np.asarray(axes_in_parent, dtype=float)
        assert ax.shape == (3, 2), "axes_in_parent must be shape (3, 2)"
        # normalise each column
        self._axes = ax / np.linalg.norm(ax, axis=0)

    @property
    def n_dof(self) -> int:
        return 2

    def rotation_matrix(self, q: np.ndarray) -> np.ndarray:
        """R^{hi} = R1(q[0]) @ R2(q[1])."""
        R1 = rodrigues(self._axes[:, 0], q[0])
        # second axis lives in the frame AFTER the first rotation
        axis2_rotated = R1 @ self._axes[:, 1]
        R2 = rodrigues(axis2_rotated, q[1])
        return R1 @ R2

    def axes_inertial(self, q: np.ndarray,
                      R0h: np.ndarray) -> np.ndarray:
        """
        φ_1 = R^{0h} @ axis_1            (first axis, before any rotation)
        φ_2 = R^{0h} @ R1(q[0]) @ axis_2 (second axis, after first rotation)
        Returns shape (3, 2).
        """
        R1   = rodrigues(self._axes[:, 0], q[0])
        phi1 = R0h @ self._axes[:, 0]
        phi2 = R0h @ R1 @ self._axes[:, 1]
        return np.column_stack([phi1, phi2])


# ------------------------------------------------------------------ #
#  Cardan3 joint  (3 DOF)                                             #
# ------------------------------------------------------------------ #

class Cardan3Joint(Joint):
    """
    3-DOF joint: main spin + 2-axis Cardan tilt.
    Used for the pole (main rotation independent, 2 tilts driven).

    Rotation sequence:
      q[0] about axis_0 (main spin, typically I3)
      q[1] about axis_1 (tilt 1, typically I1)
      q[2] about axis_2 (tilt 2, typically I2)

    Parameters
    ----------
    axes_in_parent : np.ndarray, shape (3, 3)
        Columns are the three successive rotation axes in the parent frame.
    """

    def __init__(self, axes_in_parent: np.ndarray):
        ax = np.asarray(axes_in_parent, dtype=float)
        assert ax.shape == (3, 3), "axes_in_parent must be shape (3, 3)"
        self._axes = ax / np.linalg.norm(ax, axis=0)

    @property
    def n_dof(self) -> int:
        return 3

    def rotation_matrix(self, q: np.ndarray) -> np.ndarray:
        """R^{hi} = R0(q[0]) @ R1(q[1]) @ R2(q[2])."""
        R0 = rodrigues(self._axes[:, 0], q[0])
        axis1_r = R0 @ self._axes[:, 1]
        R1 = rodrigues(axis1_r, q[1])
        axis2_r = R0 @ R1 @ self._axes[:, 2]
        R2 = rodrigues(axis2_r, q[2])
        return R0 @ R1 @ R2

    def axes_inertial(self, q: np.ndarray,
                      R0h: np.ndarray) -> np.ndarray:
        """
        φ_0 = R^{0h} @ axis_0
        φ_1 = R^{0h} @ R0(q[0]) @ axis_1
        φ_2 = R^{0h} @ R0(q[0]) @ R1(q[1]) @ axis_2
        Returns shape (3, 3).
        """
        R0   = rodrigues(self._axes[:, 0], q[0])
        axis1_r = R0 @ self._axes[:, 1]
        R1   = rodrigues(axis1_r, q[1])

        phi0 = R0h @ self._axes[:, 0]
        phi1 = R0h @ R0 @ self._axes[:, 1]
        phi2 = R0h @ R0 @ R1 @ self._axes[:, 2]
        return np.column_stack([phi0, phi1, phi2])


# ------------------------------------------------------------------ #
#  Factory function                                                    #
# ------------------------------------------------------------------ #

def make_joint(body) -> Joint:
    """
    Instantiate the correct Joint subclass from a Body object.
    Call this once during system assembly, store result alongside body.
    """
    jt = body.joint_type
    ax = body.joint_axis_in_parent

    if jt == 'fixed':
        return FixedJoint()
    elif jt == 'revolute':
        return RevoluteJoint(ax)
    elif jt == 'cardan2':
        return Cardan2Joint(ax)
    elif jt == 'cardan3':
        return Cardan3Joint(ax)
    else:
        raise ValueError(f"Unknown joint type '{jt}' for body '{body.name}'")


def make_all_joints(bodies):
    joints = []
    for body in bodies:
        joint = make_joint(body)
        joints.append(joint)
    return joints