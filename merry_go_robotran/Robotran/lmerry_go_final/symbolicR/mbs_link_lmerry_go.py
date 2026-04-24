#
#	MBsysTran - Release 8.1
#
#	Copyright 
#	Universite catholique de Louvain (UCLouvain) 
#	Mechatronic, Electrical Energy, and Dynamic systems (MEED Division) 
#	2, Place du Levant
#	1348 Louvain-la-Neuve 
#	Belgium 
#
#	http://www.robotran.be 
#
#	==> Generation Date: Thu Apr 23 18:19:45 2026
#	==> using automatic loading with extension .mbs 
#
#	==> Project name: lmerry_go
#
#	==> Number of joints: 15
#
#	==> Function: F7 - Link Forces (1D)
#
#	==> Git hash: 0e4e6a608eeee06956095d2f2ad315abdb092777
#
##

from math import sin, cos, sqrt

def link(frc, trq, Flink, Z, Zd, s, tsim):
    q = s.q
    qd = s.qd
 
# Trigonometric functions

    S4 = sin(q[4])
    C4 = cos(q[4])
    S7 = sin(q[7])
    C7 = cos(q[7])
    S10 = sin(q[10])
    C10 = cos(q[10])
    S13 = sin(q[13])
    C13 = cos(q[13])
 
# Augmented Joint Position Vectors

 
# Link anchor points Kinematics

    RLlnk2_12 = s.dpt[1,10]*C4+s.dpt[3,10]*S4
    RLlnk2_32 = -s.dpt[1,10]*S4+s.dpt[3,10]*C4
    POlnk2_12 = RLlnk2_12+s.dpt[1,1]
    POlnk2_32 = RLlnk2_32+s.dpt[3,1]
    ORlnk2_12 = qd[4]*RLlnk2_32
    ORlnk2_32 = -qd[4]*RLlnk2_12
    Plnk11 = POlnk2_12-s.dpt[1,5]
    Plnk31 = POlnk2_32-s.dpt[3,5]
    PPlnk1 = Plnk11*Plnk11+Plnk31*Plnk31
    Z1 = sqrt(PPlnk1)
    e11 = Plnk11/Z1
    e31 = Plnk31/Z1
    Zd1 = ORlnk2_12*e11+ORlnk2_32*e31
    RLlnk4_22 = s.dpt[2,13]*C7-s.dpt[3,13]*S7
    RLlnk4_32 = s.dpt[2,13]*S7+s.dpt[3,13]*C7
    POlnk4_22 = RLlnk4_22+s.dpt[2,2]
    POlnk4_32 = RLlnk4_32+s.dpt[3,2]
    ORlnk4_22 = -qd[7]*RLlnk4_32
    ORlnk4_32 = qd[7]*RLlnk4_22
    Plnk22 = POlnk4_22-s.dpt[2,6]
    Plnk32 = POlnk4_32-s.dpt[3,6]
    PPlnk2 = Plnk22*Plnk22+Plnk32*Plnk32
    Z2 = sqrt(PPlnk2)
    e22 = Plnk22/Z2
    e32 = Plnk32/Z2
    Zd2 = ORlnk4_22*e22+ORlnk4_32*e32
    RLlnk6_12 = s.dpt[1,16]*C10+s.dpt[3,16]*S10
    RLlnk6_32 = -s.dpt[1,16]*S10+s.dpt[3,16]*C10
    POlnk6_12 = RLlnk6_12+s.dpt[1,3]
    POlnk6_32 = RLlnk6_32+s.dpt[3,3]
    ORlnk6_12 = qd[10]*RLlnk6_32
    ORlnk6_32 = -qd[10]*RLlnk6_12
    Plnk13 = POlnk6_12-s.dpt[1,7]
    Plnk33 = POlnk6_32-s.dpt[3,7]
    PPlnk3 = Plnk13*Plnk13+Plnk33*Plnk33
    Z3 = sqrt(PPlnk3)
    e13 = Plnk13/Z3
    e33 = Plnk33/Z3
    Zd3 = ORlnk6_12*e13+ORlnk6_32*e33
    RLlnk8_22 = s.dpt[2,19]*C13-s.dpt[3,19]*S13
    RLlnk8_32 = s.dpt[2,19]*S13+s.dpt[3,19]*C13
    POlnk8_22 = RLlnk8_22+s.dpt[2,4]
    POlnk8_32 = RLlnk8_32+s.dpt[3,4]
    ORlnk8_22 = -qd[13]*RLlnk8_32
    ORlnk8_32 = qd[13]*RLlnk8_22
    Plnk24 = POlnk8_22-s.dpt[2,8]
    Plnk34 = POlnk8_32-s.dpt[3,8]
    PPlnk4 = Plnk24*Plnk24+Plnk34*Plnk34
    Z4 = sqrt(PPlnk4)
    e24 = Plnk24/Z4
    e34 = Plnk34/Z4
    Zd4 = ORlnk8_22*e24+ORlnk8_32*e34

# Link Forces 

    Flink1 = s.user_LinkForces(Z1,Zd1,s,tsim,1)
    Flink2 = s.user_LinkForces(Z2,Zd2,s,tsim,2)
    Flink3 = s.user_LinkForces(Z3,Zd3,s,tsim,3)
    Flink4 = s.user_LinkForces(Z4,Zd4,s,tsim,4)
 
# Link Dynamics: forces projection on body-fixed frames

    fPlnk11 = Flink1*e11
    fPlnk31 = Flink1*e31
    trqlnk3_1_2 = fPlnk11*(s.dpt[3,5]-s.l[3,3])-fPlnk31*s.dpt[1,5]
    fSlnk11 = Flink1*(e11*C4-e31*S4)
    fSlnk31 = Flink1*(e11*S4+e31*C4)
    trqlnk4_1_2 = -fSlnk11*(s.dpt[3,10]-s.l[3,4])+fSlnk31*s.dpt[1,10]
    fPlnk22 = Flink2*e22
    fPlnk32 = Flink2*e32
    frclnk3_2_3 = fPlnk31+fPlnk32
    trqlnk3_2_1 = -fPlnk22*(s.dpt[3,6]-s.l[3,3])+fPlnk32*s.dpt[2,6]
    fSlnk22 = Flink2*(e22*C7+e32*S7)
    fSlnk32 = Flink2*(-e22*S7+e32*C7)
    trqlnk7_2_1 = fSlnk22*(s.dpt[3,13]-s.l[3,7])-fSlnk32*s.dpt[2,13]
    fPlnk13 = Flink3*e13
    fPlnk33 = Flink3*e33
    frclnk3_3_1 = fPlnk11+fPlnk13
    frclnk3_3_3 = fPlnk33+frclnk3_2_3
    trqlnk3_3_2 = trqlnk3_1_2+fPlnk13*(s.dpt[3,7]-s.l[3,3])-fPlnk33*s.dpt[1,7]
    fSlnk13 = Flink3*(e13*C10-e33*S10)
    fSlnk33 = Flink3*(e13*S10+e33*C10)
    trqlnk10_3_2 = -fSlnk13*(s.dpt[3,16]-s.l[3,10])+fSlnk33*s.dpt[1,16]
    fPlnk24 = Flink4*e24
    fPlnk34 = Flink4*e34
    frclnk3_4_2 = fPlnk22+fPlnk24
    frclnk3_4_3 = fPlnk34+frclnk3_3_3
    trqlnk3_4_1 = trqlnk3_2_1-fPlnk24*(s.dpt[3,8]-s.l[3,3])+fPlnk34*s.dpt[2,8]
    fSlnk24 = Flink4*(e24*C13+e34*S13)
    fSlnk34 = Flink4*(-e24*S13+e34*C13)
    trqlnk13_4_1 = fSlnk24*(s.dpt[3,19]-s.l[3,13])-fSlnk34*s.dpt[2,19]
 
# Symbolic model output

    frc[1,3] = s.frc[1,3]+frclnk3_3_1
    frc[2,3] = s.frc[2,3]+frclnk3_4_2
    frc[3,3] = s.frc[3,3]+frclnk3_4_3
    trq[1,3] = s.trq[1,3]+trqlnk3_4_1
    trq[2,3] = s.trq[2,3]+trqlnk3_3_2
    frc[1,4] = s.frc[1,4]-fSlnk11
    frc[3,4] = s.frc[3,4]-fSlnk31
    trq[2,4] = s.trq[2,4]+trqlnk4_1_2
    frc[2,7] = s.frc[2,7]-fSlnk22
    frc[3,7] = s.frc[3,7]-fSlnk32
    trq[1,7] = s.trq[1,7]+trqlnk7_2_1
    frc[1,10] = s.frc[1,10]-fSlnk13
    frc[3,10] = s.frc[3,10]-fSlnk33
    trq[2,10] = s.trq[2,10]+trqlnk10_3_2
    frc[2,13] = s.frc[2,13]-fSlnk24
    frc[3,13] = s.frc[3,13]-fSlnk34
    trq[1,13] = s.trq[1,13]+trqlnk13_4_1
 
# Symbolic model output

    Z[1] = Z1
    Zd[1] = Zd1
    Flink[1] = Flink1
    Z[2] = Z2
    Zd[2] = Zd2
    Flink[2] = Flink2
    Z[3] = Z3
    Zd[3] = Zd3
    Flink[3] = Flink3
    Z[4] = Z4
    Zd[4] = Zd4
    Flink[4] = Flink4

# Number of continuation lines = 0


