#
# LMECA2802 - Small example of project
#
# @date 2026
# @author 
# 
# Universite catholique de Louvain

# To run this code, you should have installed MBsysPy (see installation instructions in the Robotran tutorials (robotran.be))
# You should also have downloaded cartpendulum.mbs and placed it in the same folder


# import useful modules and functions
from math import sin, cos, pi
import numpy as np
from scipy.integrate import solve_ivp                # To do time integration
import MBsysPy as Robotran                           # To read your .mbs file
import matplotlib.pyplot as plt                      # To plot nice graphs of your results
from MBsysPy import algebra                          # To have access to useful mathematical functions

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# Loading data from your .mbs file (in the current folder)
mbs_data = Robotran.MbsData("Cartpendulum.mbs")      

# General parameters
g = mbs_data.g[3]                                    # Recovering gravity along z axis

# Masses
m1 = mbs_data.m[mbs_data.body_id["Cart"]]            # Recovering the mass of the body "Cart" 
m2 = mbs_data.m[mbs_data.body_id["Pendulum_Mass"]]   # Recovering the mass of the body "Pendulum_Mass"

# Pendulum length
Lp_id = mbs_data.points_id["Pendulum_String"]["Lp"]  # Recovering the "id" that identifies an anchor point in MBSysPad : mbs_data.points_id["body_name"]["anchor_point_name"]
Lp_array = mbs_data.dpt[1:4, Lp_id]                  # Recovering the [X, Y, Z] data related to the anchor point through its id. 
Lp = Lp_array[2]                                     # Here, only interested in the Z component for the length of the pendulum

# Initial conditions
q1 = mbs_data.q0[mbs_data.joint_id["Cart_T2"]]       # Initial position coordinate of the cart
q2 = mbs_data.q0[mbs_data.joint_id["Pendulum_R1"]]   # Initial position coordinate of the pendulum
qd1 = mbs_data.qd0[mbs_data.joint_id["Cart_T2"]]     # Initial velocity coordinate of the cart
qd2 = mbs_data.qd0[mbs_data.joint_id["Pendulum_R1"]] # Initial velocity coordinate of the pendulum

# Maximum force applied to cart
Fmax = mbs_data.user_model["Force"]["Fmax"]          # Recovering the maximum force applied to the cart through the User Model
t0 = mbs_data.user_model["Time"]["t0"]               # Recovering the starting time through the User Model
t1 = mbs_data.user_model["Time"]["t1"]               # Recovering the ending time through the User Model

# To know what other information is contained inside mbs_data check the Robotran website
# robotran.be > Docs > Documentation > MBsysPy > Search "mbs_data" in the search bar (top left)
# An interesting module (called algebra) for mathematical operations on vectors starting at index 1 is available with Robotran
# robotran.be > Docs > Documentation > MBsysPy > Search "algebra" in the search bar (top left)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *    
# Writing the equations of movement of my system "Cart + pendulum"
# Your code should be a generalized implementation of a multibody formalism
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
def compute_derivatives(t, y, data):
    """ Compute the derivatives yd for a given state y of the system.
    The derivatives are computed at the given time t with
    the parameters values in the given data structure.
    
    It is assumed that the state vector y contains the following states:
      y = [q1, q2, qd1, qd2] with:
         - q1: the cart position
         - q2: the pendulum angular position 
         - qd1: the cart velocity
         - qd2: the pendulum angular velocity 

    :param  t: the time instant when to compute the derivatives.
    :param  y: the numpy array containing the states 
    :return: yd a numpy array containing the states derivatives  yd = [qd1, qd2, qdd1, qdd2]
    :param data: the mbs_data object containing the parameters of the model
    """
    # Recovering parameters -> redundant from earlier, do it once only, organize your code as you feel best
    m1 = mbs_data.m[mbs_data.body_id["Cart"]]
    m2 = mbs_data.m[mbs_data.body_id["Pendulum_Mass"]]
    Lp_id = mbs_data.points_id["Pendulum_String"]["Lp"]  
    Lp_array = mbs_data.dpt[1:4, Lp_id]                   
    Lp = Lp_array[2]
    Fmax = mbs_data.user_model["Force"]["Fmax"]                                      

    #### Initialization of matrices and vectors ####
    M = np.zeros((2,2))
    c = np.zeros((2,1))
    Q = np.zeros((2,1))  
    yd = np.zeros(4)
    F = np.zeros((2,1)) 
    
    ### Values of matrices et vectors ###
    #M
    M[0][0] = m1 + m2
    M[0][1] = m2*Lp*np.cos(y[1])/2
    M[1][0] = m2*Lp*np.cos(y[1])/2
    M[1][1] = m2*Lp**2/3
    
    #Minv
    invM = np.zeros((2,2))
    detM = M[0][0]*M[1][1] - M[0][1]*M[1][0]
    invM[0][0] = 1/detM*M[1][1]
    invM[0][1] = -1/detM*M[0][1]
    invM[1][0] = -1/detM*M[1][0] 
    invM[1][1] = 1/detM*M[0][0]   
    
    # Apply constant force on cart
    F_cart = Fmax 
        
    #c
    c[0][0] = -(m2*Lp/2)*y[3]**2*np.sin(y[1])
    c[1][0] = -g*np.sin(y[1])*m2*Lp/2
    
    ### Calculate yd ###
    yd[0] = y[2]
    yd[1] = y[3]
    
    Q[0][0] = F_cart
    Q[1][0] = 0.0
    
    F[0][0] = Q[0][0] - c[0][0]
    F[1][0] = Q[1][0] - c[1][0]
    
    yd[2] = invM[0][0]*F[0][0] + invM[0][1]*F[1][0]
    yd[3] = invM[1][0]*F[0][0] + invM[1][1]*F[1][0]
    
    return yd


#  * * * * * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * *
# Calling the time integration of the equations to recover the dynamics
# of the system
#  * * * * * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * *
def compute_dynamic_response(data):
    """  Compute the time evolution of the dynamic response of the cart-pendulum system
         for the given data.
    Initial and final time are determined
    by the t0 and t1 parameter of the user models.
    Results are saved to three text files named dirdyn_q.res, dirdyn_qd.res and dirdyn_qdd.res
 
    Time evolution is computed using a time integrator (typically Runge-Kutta).
 
    :param data: the mbs_data object containing the parameters of the model
    """
    # ### Runge Kutta ###   should be called via solve_ivp()
    # to pass the mbs_data object to compute_derivative function in solve_ivp, you may use lambda mechanism:
    #
    #    fprime = lambda t,y: compute_derivatives(t, y, data)
    #
    # fprime can be viewed as a function that takes two arguments: t, y
    # this fprime function can be provided to solve_ivp
    # Note that you can change the tolerances with rtol and atol options (see online solve_iv doc)
    #
    
    ### Data ### -> redundant from earlier, do it once only, organize your code as you feel best
    t0 = mbs_data.user_model["Time"]["t0"]               # Recovering the starting time through the User Model
    t1 = mbs_data.user_model["Time"]["t1"]               # Recovering the ending time through the User Model

    tspan=(t0, t1)
    n_state = 2 #number of states  of q and qd
    
    #### Matrices initialization ####
    q = np.zeros(n_state)
    qd = np.zeros(n_state)

    
    ### Initial values ###
    q1 = mbs_data.q0[mbs_data.joint_id["Cart_T2"]]       # Initial position coordinate of the cart
    q2 = mbs_data.q0[mbs_data.joint_id["Pendulum_R1"]]   # Initial position coordinate of the pendulum
    qd1 = mbs_data.qd0[mbs_data.joint_id["Cart_T2"]]     # Initial velocity coordinate of the cart
    qd2 = mbs_data.qd0[mbs_data.joint_id["Pendulum_R1"]] # Initial velocity coordinate of the pendulum
    #q
    q0 = np.array([q1, q2])
    qd0 = np.array([qd1, qd2])
    
    q[:] = q0[:]
    qd[:] = qd0[:]
    
    #y
    y0 = np.zeros(4)
    y0[:n_state] = q[:]
    y0[n_state:2*n_state] = qd[:]

    ### Runge Kutta ###   
    sol = solve_ivp(fun=lambda t,y: compute_derivatives(t, y, data), t_span=tspan, y0=y0, rtol=1e-5, atol=1e-8)
    
    # Adjust initial step and max step to find Robotran precision
    # sol = solve_ivp(fun=lambda t,y: compute_derivatives(t, y, data), t_span=tspan, y0=y0, first_step=1e-4, max_step=1e-4)
    
    
    ### Display balance values ###
    print('q1 = ',sol.y[0][-1],'q2 = ',sol.y[1][-1])
    

    ### Calculate qdd in post process ###
    ydmatrix = np.zeros([2*n_state, sol.t.__len__()])
    for i in range(sol.t.__len__()):
        np.array([sol.y[0][i],sol.y[1][i],sol.y[2][i],sol.y[3][i]])
        a = compute_derivatives(sol.t[i], np.array([sol.y[0][i],sol.y[1][i],sol.y[2][i],sol.y[3][i]]),data)
        ydmatrix[0][i] = a[0]
        ydmatrix[1][i] = a[1] 
        ydmatrix[2][i] = a[2] 
        ydmatrix[3][i] = a[3]

    
    ### Files creation and results record ###
    np.savetxt('dirdyn_q.res', np.transpose(np.array([sol.t,sol.y[0],sol.y[1]])))
    np.savetxt('dirdyn_qd.res', np.transpose(np.array([sol.t,sol.y[2],sol.y[3]])))
    np.savetxt('dirdyn_qdd.res', np.transpose(np.array([sol.t,ydmatrix[2],ydmatrix[3]])))


#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# Plotting the results of the time simulation
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
def plots():
    q = np.loadtxt('dirdyn_q.res')
    qd = np.loadtxt('dirdyn_qd.res')
    qdd = np.loadtxt('dirdyn_qdd.res')

    # Creation of 1 figure with 6 plots
    fig, axs = plt.subplots(3, 2, sharex=True, figsize=(16,9))

    fig.suptitle("Cart with pendulum kinematics")

    axs[0][0].plot(q[:, 0], q[:, 1])
    axs[0][0].set_ylabel("Position (m)", fontsize=12)
    axs[0][0].grid()
    axs[0][0].set_title("Kinematics of the cart", fontsize=12)

    axs[0][1].plot(q[:, 0], q[:, 2] * (180/np.pi))
    axs[0][1].set_ylabel("Position (deg)", fontsize=12)
    axs[0][1].grid()
    axs[0][1].set_title("Kinematics of the pendulum", fontsize=12)

    axs[1][0].plot(qd[:, 0], qd[:, 1])
    axs[1][0].set_ylabel("Velocity (m/s)", fontsize=12)
    axs[1][0].grid()

    axs[1][1].plot(qd[:, 0], qd[:, 2] * (180/np.pi))
    axs[1][1].set_ylabel("Velocity (deg/s)", fontsize=12)
    axs[1][1].grid()

    axs[2][0].plot(qdd[:, 0], qdd[:, 1])
    axs[2][0].set_ylabel("Acceleration (m/s2)", fontsize=12)
    axs[2][0].set_xlabel("Time (s)", fontsize=12)
    axs[2][0].grid()

    axs[2][1].plot(qdd[:, 0], qdd[:, 2] * (180/np.pi))
    axs[2][1].set_ylabel("Acceleration (deg/s2)", fontsize=12)
    axs[2][1].set_xlabel("Time (s)", fontsize=12)
    axs[2][1].grid()

    plt.show()



# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# Main function

if __name__ == '__main__':
    mbs_data = Robotran.MbsData("Cartpendulum.mbs")
    compute_dynamic_response(mbs_data) 
    plots()
    
