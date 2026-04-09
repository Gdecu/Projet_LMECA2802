from merry_go_robotran.data.carousel_data import build_bodies
from merry_go_robotran.neri.joint import make_all_joints
from merry_go_robotran.neri.state import MBState
from merry_go_robotran.neri.forward import forward_pass
from merry_go_robotran.neri.backward import backward_pass
from merry_go_robotran.neri.joint import make_all_joints

import numpy as np

T_END = 30
TIMESTEP = 0.5
Q_0 = np.zeros(15)
QD_0 = np.zeros_like(Q_0)

def main():
    bodies = build_bodies()
    joints = make_all_joints(bodies)
    state = MBState(n_bodies=len(bodies), n_dof=15)  # allocate once
    q = Q_0
    qd = QD_0
    for t_act in np.arange(0,T_END, TIMESTEP):
        state.set_q(q, qd, t_act)  # load current state
        forward_pass(bodies, joints, state)  # fill ω, α, R ...
        backward_pass(bodies, joints, state)  # fill W, F, L, Q

if __name__ == "__main__":
    main()