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
#	==> Generation Date: Wed Apr 15 02:34:06 2026
#	==> using automatic loading with extension .mbs 
#
#	==> Project name: lmerry_go
#
#	==> Number of joints: 15
#
#	==> Function: F6 - Sensors Kinematics
#
#	==> Git hash: 0e4e6a608eeee06956095d2f2ad315abdb092777
#
##

from math import sin, cos, sqrt

def sensor(sens, s, isens):
  q = s.q
  qd = s.qd
  qdd = s.qdd

  dpt = s.dpt
 
# Trigonometric functions

  S1 = sin(q[1])
  C1 = cos(q[1])
  S2 = sin(q[2])
  C2 = cos(q[2])
  S3 = sin(q[3])
  C3 = cos(q[3])
  S4 = sin(q[4])
  C4 = cos(q[4])
  S5 = sin(q[5])
  C5 = cos(q[5])
  S6 = sin(q[6])
  C6 = cos(q[6])
  S7 = sin(q[7])
  C7 = cos(q[7])
  S8 = sin(q[8])
  C8 = cos(q[8])
  S9 = sin(q[9])
  C9 = cos(q[9])
  S10 = sin(q[10])
  C10 = cos(q[10])
  S11 = sin(q[11])
  C11 = cos(q[11])
  S12 = sin(q[12])
  C12 = cos(q[12])
  S13 = sin(q[13])
  C13 = cos(q[13])
  S14 = sin(q[14])
  C14 = cos(q[14])
  S15 = sin(q[15])
  C15 = cos(q[15])
 
# Augmented Joint Position Vectors

 
# Sensor Kinematics


  if (isens == 1): 

    ROcp1_42 = -S1*C2
    ROcp1_52 = C1*C2
    ROcp1_72 = S1*S2
    ROcp1_82 = -C1*S2
    ROcp1_13 = -ROcp1_72*S3+C1*C3
    ROcp1_23 = -ROcp1_82*S3+S1*C3
    ROcp1_33 = -C2*S3
    ROcp1_73 = ROcp1_72*C3+C1*S3
    ROcp1_83 = ROcp1_82*C3+S1*S3
    ROcp1_93 = C2*C3
    ROcp1_14 = ROcp1_13*C4-ROcp1_73*S4
    ROcp1_24 = ROcp1_23*C4-ROcp1_83*S4
    ROcp1_34 = ROcp1_33*C4-ROcp1_93*S4
    ROcp1_74 = ROcp1_13*S4+ROcp1_73*C4
    ROcp1_84 = ROcp1_23*S4+ROcp1_83*C4
    ROcp1_94 = ROcp1_33*S4+ROcp1_93*C4
    ROcp1_45 = ROcp1_42*C5+ROcp1_74*S5
    ROcp1_55 = ROcp1_52*C5+ROcp1_84*S5
    ROcp1_65 = ROcp1_94*S5+S2*C5
    ROcp1_75 = -ROcp1_42*S5+ROcp1_74*C5
    ROcp1_85 = -ROcp1_52*S5+ROcp1_84*C5
    ROcp1_95 = ROcp1_94*C5-S2*S5
    ROcp1_16 = ROcp1_14*C6-ROcp1_75*S6
    ROcp1_26 = ROcp1_24*C6-ROcp1_85*S6
    ROcp1_36 = ROcp1_34*C6-ROcp1_95*S6
    ROcp1_76 = ROcp1_14*S6+ROcp1_75*C6
    ROcp1_86 = ROcp1_24*S6+ROcp1_85*C6
    ROcp1_96 = ROcp1_34*S6+ROcp1_95*C6
    OMcp1_12 = qd[2]*C1
    OMcp1_22 = qd[2]*S1
    OPcp1_12 = qdd[2]*C1-qd[1]*qd[2]*S1
    OPcp1_22 = qdd[2]*S1+qd[1]*qd[2]*C1
    OMcp1_13 = OMcp1_12+ROcp1_42*qd[3]
    OMcp1_23 = OMcp1_22+ROcp1_52*qd[3]
    OMcp1_33 = qd[1]+qd[3]*S2
    OPcp1_13 = OPcp1_12+ROcp1_42*qdd[3]+qd[3]*(OMcp1_22*S2-ROcp1_52*qd[1])
    OPcp1_23 = OPcp1_22+ROcp1_52*qdd[3]+qd[3]*(-OMcp1_12*S2+ROcp1_42*qd[1])
    OPcp1_33 = qdd[1]+qdd[3]*S2+qd[3]*(OMcp1_12*ROcp1_52-OMcp1_22*ROcp1_42)
    RLcp1_14 = ROcp1_13*s.dpt[1,1]+ROcp1_73*s.dpt[3,1]
    RLcp1_24 = ROcp1_23*s.dpt[1,1]+ROcp1_83*s.dpt[3,1]
    RLcp1_34 = ROcp1_33*s.dpt[1,1]+ROcp1_93*s.dpt[3,1]
    OMcp1_14 = OMcp1_13+ROcp1_42*qd[4]
    OMcp1_24 = OMcp1_23+ROcp1_52*qd[4]
    OMcp1_34 = OMcp1_33+qd[4]*S2
    ORcp1_14 = OMcp1_23*RLcp1_34-OMcp1_33*RLcp1_24
    ORcp1_24 = -OMcp1_13*RLcp1_34+OMcp1_33*RLcp1_14
    ORcp1_34 = OMcp1_13*RLcp1_24-OMcp1_23*RLcp1_14
    OPcp1_14 = OPcp1_13+ROcp1_42*qdd[4]+qd[4]*(OMcp1_23*S2-OMcp1_33*ROcp1_52)
    OPcp1_24 = OPcp1_23+ROcp1_52*qdd[4]+qd[4]*(-OMcp1_13*S2+OMcp1_33*ROcp1_42)
    OPcp1_34 = OPcp1_33+qdd[4]*S2+qd[4]*(OMcp1_13*ROcp1_52-OMcp1_23*ROcp1_42)
    ACcp1_14 = OMcp1_23*ORcp1_34-OMcp1_33*ORcp1_24+OPcp1_23*RLcp1_34-OPcp1_33*RLcp1_24
    ACcp1_24 = -OMcp1_13*ORcp1_34+OMcp1_33*ORcp1_14-OPcp1_13*RLcp1_34+OPcp1_33*RLcp1_14
    ACcp1_34 = OMcp1_13*ORcp1_24-OMcp1_23*ORcp1_14+OPcp1_13*RLcp1_24-OPcp1_23*RLcp1_14
    RLcp1_15 = ROcp1_74*s.dpt[3,9]
    RLcp1_25 = ROcp1_84*s.dpt[3,9]
    RLcp1_35 = ROcp1_94*s.dpt[3,9]
    POcp1_15 = RLcp1_14+RLcp1_15
    POcp1_25 = RLcp1_24+RLcp1_25
    POcp1_35 = RLcp1_34+RLcp1_35
    OMcp1_15 = OMcp1_14+ROcp1_14*qd[5]
    OMcp1_25 = OMcp1_24+ROcp1_24*qd[5]
    OMcp1_35 = OMcp1_34+ROcp1_34*qd[5]
    ORcp1_15 = OMcp1_24*RLcp1_35-OMcp1_34*RLcp1_25
    ORcp1_25 = -OMcp1_14*RLcp1_35+OMcp1_34*RLcp1_15
    ORcp1_35 = OMcp1_14*RLcp1_25-OMcp1_24*RLcp1_15
    VIcp1_15 = ORcp1_14+ORcp1_15
    VIcp1_25 = ORcp1_24+ORcp1_25
    VIcp1_35 = ORcp1_34+ORcp1_35
    OPcp1_15 = OPcp1_14+ROcp1_14*qdd[5]+qd[5]*(OMcp1_24*ROcp1_34-OMcp1_34*ROcp1_24)
    OPcp1_25 = OPcp1_24+ROcp1_24*qdd[5]+qd[5]*(-OMcp1_14*ROcp1_34+OMcp1_34*ROcp1_14)
    OPcp1_35 = OPcp1_34+ROcp1_34*qdd[5]+qd[5]*(OMcp1_14*ROcp1_24-OMcp1_24*ROcp1_14)
    ACcp1_15 = ACcp1_14+OMcp1_24*ORcp1_35-OMcp1_34*ORcp1_25+OPcp1_24*RLcp1_35-OPcp1_34*RLcp1_25
    ACcp1_25 = ACcp1_24-OMcp1_14*ORcp1_35+OMcp1_34*ORcp1_15-OPcp1_14*RLcp1_35+OPcp1_34*RLcp1_15
    ACcp1_35 = ACcp1_34+OMcp1_14*ORcp1_25-OMcp1_24*ORcp1_15+OPcp1_14*RLcp1_25-OPcp1_24*RLcp1_15
    OMcp1_16 = OMcp1_15+ROcp1_45*qd[6]
    OMcp1_26 = OMcp1_25+ROcp1_55*qd[6]
    OMcp1_36 = OMcp1_35+ROcp1_65*qd[6]
    OPcp1_16 = OPcp1_15+ROcp1_45*qdd[6]+qd[6]*(OMcp1_25*ROcp1_65-OMcp1_35*ROcp1_55)
    OPcp1_26 = OPcp1_25+ROcp1_55*qdd[6]+qd[6]*(-OMcp1_15*ROcp1_65+OMcp1_35*ROcp1_45)
    OPcp1_36 = OPcp1_35+ROcp1_65*qdd[6]+qd[6]*(OMcp1_15*ROcp1_55-OMcp1_25*ROcp1_45)
    sens.P[1] = POcp1_15
    sens.P[2] = POcp1_25
    sens.P[3] = POcp1_35
    sens.R[1,1] = ROcp1_16
    sens.R[1,2] = ROcp1_26
    sens.R[1,3] = ROcp1_36
    sens.R[2,1] = ROcp1_45
    sens.R[2,2] = ROcp1_55
    sens.R[2,3] = ROcp1_65
    sens.R[3,1] = ROcp1_76
    sens.R[3,2] = ROcp1_86
    sens.R[3,3] = ROcp1_96
    sens.V[1] = VIcp1_15
    sens.V[2] = VIcp1_25
    sens.V[3] = VIcp1_35
    sens.OM[1] = OMcp1_16
    sens.OM[2] = OMcp1_26
    sens.OM[3] = OMcp1_36
    sens.A[1] = ACcp1_15
    sens.A[2] = ACcp1_25
    sens.A[3] = ACcp1_35
    sens.OMP[1] = OPcp1_16
    sens.OMP[2] = OPcp1_26
    sens.OMP[3] = OPcp1_36

  if (isens == 2): 

    ROcp2_42 = -S1*C2
    ROcp2_52 = C1*C2
    ROcp2_72 = S1*S2
    ROcp2_82 = -C1*S2
    ROcp2_13 = -ROcp2_72*S3+C1*C3
    ROcp2_23 = -ROcp2_82*S3+S1*C3
    ROcp2_33 = -C2*S3
    ROcp2_73 = ROcp2_72*C3+C1*S3
    ROcp2_83 = ROcp2_82*C3+S1*S3
    ROcp2_93 = C2*C3
    ROcp2_17 = ROcp2_13*C7-ROcp2_73*S7
    ROcp2_27 = ROcp2_23*C7-ROcp2_83*S7
    ROcp2_37 = ROcp2_33*C7-ROcp2_93*S7
    ROcp2_77 = ROcp2_13*S7+ROcp2_73*C7
    ROcp2_87 = ROcp2_23*S7+ROcp2_83*C7
    ROcp2_97 = ROcp2_33*S7+ROcp2_93*C7
    ROcp2_48 = ROcp2_42*C8+ROcp2_77*S8
    ROcp2_58 = ROcp2_52*C8+ROcp2_87*S8
    ROcp2_68 = ROcp2_97*S8+S2*C8
    ROcp2_78 = -ROcp2_42*S8+ROcp2_77*C8
    ROcp2_88 = -ROcp2_52*S8+ROcp2_87*C8
    ROcp2_98 = ROcp2_97*C8-S2*S8
    ROcp2_19 = ROcp2_17*C9-ROcp2_78*S9
    ROcp2_29 = ROcp2_27*C9-ROcp2_88*S9
    ROcp2_39 = ROcp2_37*C9-ROcp2_98*S9
    ROcp2_79 = ROcp2_17*S9+ROcp2_78*C9
    ROcp2_89 = ROcp2_27*S9+ROcp2_88*C9
    ROcp2_99 = ROcp2_37*S9+ROcp2_98*C9
    OMcp2_12 = qd[2]*C1
    OMcp2_22 = qd[2]*S1
    OPcp2_12 = qdd[2]*C1-qd[1]*qd[2]*S1
    OPcp2_22 = qdd[2]*S1+qd[1]*qd[2]*C1
    OMcp2_13 = OMcp2_12+ROcp2_42*qd[3]
    OMcp2_23 = OMcp2_22+ROcp2_52*qd[3]
    OMcp2_33 = qd[1]+qd[3]*S2
    OPcp2_13 = OPcp2_12+ROcp2_42*qdd[3]+qd[3]*(OMcp2_22*S2-ROcp2_52*qd[1])
    OPcp2_23 = OPcp2_22+ROcp2_52*qdd[3]+qd[3]*(-OMcp2_12*S2+ROcp2_42*qd[1])
    OPcp2_33 = qdd[1]+qdd[3]*S2+qd[3]*(OMcp2_12*ROcp2_52-OMcp2_22*ROcp2_42)
    RLcp2_14 = ROcp2_42*s.dpt[2,2]+ROcp2_73*s.dpt[3,2]
    RLcp2_24 = ROcp2_52*s.dpt[2,2]+ROcp2_83*s.dpt[3,2]
    RLcp2_34 = ROcp2_93*s.dpt[3,2]+s.dpt[2,2]*S2
    OMcp2_14 = OMcp2_13+ROcp2_42*qd[7]
    OMcp2_24 = OMcp2_23+ROcp2_52*qd[7]
    OMcp2_34 = OMcp2_33+qd[7]*S2
    ORcp2_14 = OMcp2_23*RLcp2_34-OMcp2_33*RLcp2_24
    ORcp2_24 = -OMcp2_13*RLcp2_34+OMcp2_33*RLcp2_14
    ORcp2_34 = OMcp2_13*RLcp2_24-OMcp2_23*RLcp2_14
    OPcp2_14 = OPcp2_13+ROcp2_42*qdd[7]+qd[7]*(OMcp2_23*S2-OMcp2_33*ROcp2_52)
    OPcp2_24 = OPcp2_23+ROcp2_52*qdd[7]+qd[7]*(-OMcp2_13*S2+OMcp2_33*ROcp2_42)
    OPcp2_34 = OPcp2_33+qdd[7]*S2+qd[7]*(OMcp2_13*ROcp2_52-OMcp2_23*ROcp2_42)
    ACcp2_14 = OMcp2_23*ORcp2_34-OMcp2_33*ORcp2_24+OPcp2_23*RLcp2_34-OPcp2_33*RLcp2_24
    ACcp2_24 = -OMcp2_13*ORcp2_34+OMcp2_33*ORcp2_14-OPcp2_13*RLcp2_34+OPcp2_33*RLcp2_14
    ACcp2_34 = OMcp2_13*ORcp2_24-OMcp2_23*ORcp2_14+OPcp2_13*RLcp2_24-OPcp2_23*RLcp2_14
    RLcp2_15 = ROcp2_77*s.dpt[3,12]
    RLcp2_25 = ROcp2_87*s.dpt[3,12]
    RLcp2_35 = ROcp2_97*s.dpt[3,12]
    POcp2_15 = RLcp2_14+RLcp2_15
    POcp2_25 = RLcp2_24+RLcp2_25
    POcp2_35 = RLcp2_34+RLcp2_35
    OMcp2_15 = OMcp2_14+ROcp2_17*qd[8]
    OMcp2_25 = OMcp2_24+ROcp2_27*qd[8]
    OMcp2_35 = OMcp2_34+ROcp2_37*qd[8]
    ORcp2_15 = OMcp2_24*RLcp2_35-OMcp2_34*RLcp2_25
    ORcp2_25 = -OMcp2_14*RLcp2_35+OMcp2_34*RLcp2_15
    ORcp2_35 = OMcp2_14*RLcp2_25-OMcp2_24*RLcp2_15
    VIcp2_15 = ORcp2_14+ORcp2_15
    VIcp2_25 = ORcp2_24+ORcp2_25
    VIcp2_35 = ORcp2_34+ORcp2_35
    OPcp2_15 = OPcp2_14+ROcp2_17*qdd[8]+qd[8]*(OMcp2_24*ROcp2_37-OMcp2_34*ROcp2_27)
    OPcp2_25 = OPcp2_24+ROcp2_27*qdd[8]+qd[8]*(-OMcp2_14*ROcp2_37+OMcp2_34*ROcp2_17)
    OPcp2_35 = OPcp2_34+ROcp2_37*qdd[8]+qd[8]*(OMcp2_14*ROcp2_27-OMcp2_24*ROcp2_17)
    ACcp2_15 = ACcp2_14+OMcp2_24*ORcp2_35-OMcp2_34*ORcp2_25+OPcp2_24*RLcp2_35-OPcp2_34*RLcp2_25
    ACcp2_25 = ACcp2_24-OMcp2_14*ORcp2_35+OMcp2_34*ORcp2_15-OPcp2_14*RLcp2_35+OPcp2_34*RLcp2_15
    ACcp2_35 = ACcp2_34+OMcp2_14*ORcp2_25-OMcp2_24*ORcp2_15+OPcp2_14*RLcp2_25-OPcp2_24*RLcp2_15
    OMcp2_16 = OMcp2_15+ROcp2_48*qd[9]
    OMcp2_26 = OMcp2_25+ROcp2_58*qd[9]
    OMcp2_36 = OMcp2_35+ROcp2_68*qd[9]
    OPcp2_16 = OPcp2_15+ROcp2_48*qdd[9]+qd[9]*(OMcp2_25*ROcp2_68-OMcp2_35*ROcp2_58)
    OPcp2_26 = OPcp2_25+ROcp2_58*qdd[9]+qd[9]*(-OMcp2_15*ROcp2_68+OMcp2_35*ROcp2_48)
    OPcp2_36 = OPcp2_35+ROcp2_68*qdd[9]+qd[9]*(OMcp2_15*ROcp2_58-OMcp2_25*ROcp2_48)
    RLcp2_17 = ROcp2_79*s.dpt[3,14]
    RLcp2_27 = ROcp2_89*s.dpt[3,14]
    RLcp2_37 = ROcp2_99*s.dpt[3,14]
    POcp2_17 = POcp2_15+RLcp2_17
    POcp2_27 = POcp2_25+RLcp2_27
    POcp2_37 = POcp2_35+RLcp2_37
    ORcp2_17 = OMcp2_26*RLcp2_37-OMcp2_36*RLcp2_27
    ORcp2_27 = -OMcp2_16*RLcp2_37+OMcp2_36*RLcp2_17
    ORcp2_37 = OMcp2_16*RLcp2_27-OMcp2_26*RLcp2_17
    VIcp2_17 = ORcp2_17+VIcp2_15
    VIcp2_27 = ORcp2_27+VIcp2_25
    VIcp2_37 = ORcp2_37+VIcp2_35
    ACcp2_17 = ACcp2_15+OMcp2_26*ORcp2_37-OMcp2_36*ORcp2_27+OPcp2_26*RLcp2_37-OPcp2_36*RLcp2_27
    ACcp2_27 = ACcp2_25-OMcp2_16*ORcp2_37+OMcp2_36*ORcp2_17-OPcp2_16*RLcp2_37+OPcp2_36*RLcp2_17
    ACcp2_37 = ACcp2_35+OMcp2_16*ORcp2_27-OMcp2_26*ORcp2_17+OPcp2_16*RLcp2_27-OPcp2_26*RLcp2_17
    sens.P[1] = POcp2_17
    sens.P[2] = POcp2_27
    sens.P[3] = POcp2_37
    sens.R[1,1] = ROcp2_19
    sens.R[1,2] = ROcp2_29
    sens.R[1,3] = ROcp2_39
    sens.R[2,1] = ROcp2_48
    sens.R[2,2] = ROcp2_58
    sens.R[2,3] = ROcp2_68
    sens.R[3,1] = ROcp2_79
    sens.R[3,2] = ROcp2_89
    sens.R[3,3] = ROcp2_99
    sens.V[1] = VIcp2_17
    sens.V[2] = VIcp2_27
    sens.V[3] = VIcp2_37
    sens.OM[1] = OMcp2_16
    sens.OM[2] = OMcp2_26
    sens.OM[3] = OMcp2_36
    sens.A[1] = ACcp2_17
    sens.A[2] = ACcp2_27
    sens.A[3] = ACcp2_37
    sens.OMP[1] = OPcp2_16
    sens.OMP[2] = OPcp2_26
    sens.OMP[3] = OPcp2_36

  if (isens == 3): 

    ROcp3_42 = -S1*C2
    ROcp3_52 = C1*C2
    ROcp3_72 = S1*S2
    ROcp3_82 = -C1*S2
    ROcp3_13 = -ROcp3_72*S3+C1*C3
    ROcp3_23 = -ROcp3_82*S3+S1*C3
    ROcp3_33 = -C2*S3
    ROcp3_73 = ROcp3_72*C3+C1*S3
    ROcp3_83 = ROcp3_82*C3+S1*S3
    ROcp3_93 = C2*C3
    ROcp3_110 = ROcp3_13*C10-ROcp3_73*S10
    ROcp3_210 = ROcp3_23*C10-ROcp3_83*S10
    ROcp3_310 = ROcp3_33*C10-ROcp3_93*S10
    ROcp3_710 = ROcp3_13*S10+ROcp3_73*C10
    ROcp3_810 = ROcp3_23*S10+ROcp3_83*C10
    ROcp3_910 = ROcp3_33*S10+ROcp3_93*C10
    ROcp3_411 = ROcp3_42*C11+ROcp3_710*S11
    ROcp3_511 = ROcp3_52*C11+ROcp3_810*S11
    ROcp3_611 = ROcp3_910*S11+C11*S2
    ROcp3_711 = -ROcp3_42*S11+ROcp3_710*C11
    ROcp3_811 = -ROcp3_52*S11+ROcp3_810*C11
    ROcp3_911 = ROcp3_910*C11-S11*S2
    ROcp3_112 = ROcp3_110*C12-ROcp3_711*S12
    ROcp3_212 = ROcp3_210*C12-ROcp3_811*S12
    ROcp3_312 = ROcp3_310*C12-ROcp3_911*S12
    ROcp3_712 = ROcp3_110*S12+ROcp3_711*C12
    ROcp3_812 = ROcp3_210*S12+ROcp3_811*C12
    ROcp3_912 = ROcp3_310*S12+ROcp3_911*C12
    OMcp3_12 = qd[2]*C1
    OMcp3_22 = qd[2]*S1
    OPcp3_12 = qdd[2]*C1-qd[1]*qd[2]*S1
    OPcp3_22 = qdd[2]*S1+qd[1]*qd[2]*C1
    OMcp3_13 = OMcp3_12+ROcp3_42*qd[3]
    OMcp3_23 = OMcp3_22+ROcp3_52*qd[3]
    OMcp3_33 = qd[1]+qd[3]*S2
    OPcp3_13 = OPcp3_12+ROcp3_42*qdd[3]+qd[3]*(OMcp3_22*S2-ROcp3_52*qd[1])
    OPcp3_23 = OPcp3_22+ROcp3_52*qdd[3]+qd[3]*(-OMcp3_12*S2+ROcp3_42*qd[1])
    OPcp3_33 = qdd[1]+qdd[3]*S2+qd[3]*(OMcp3_12*ROcp3_52-OMcp3_22*ROcp3_42)
    RLcp3_14 = ROcp3_13*s.dpt[1,3]
    RLcp3_24 = ROcp3_23*s.dpt[1,3]
    RLcp3_34 = ROcp3_33*s.dpt[1,3]
    OMcp3_14 = OMcp3_13+ROcp3_42*qd[10]
    OMcp3_24 = OMcp3_23+ROcp3_52*qd[10]
    OMcp3_34 = OMcp3_33+qd[10]*S2
    ORcp3_14 = OMcp3_23*RLcp3_34-OMcp3_33*RLcp3_24
    ORcp3_24 = -OMcp3_13*RLcp3_34+OMcp3_33*RLcp3_14
    ORcp3_34 = OMcp3_13*RLcp3_24-OMcp3_23*RLcp3_14
    OPcp3_14 = OPcp3_13+ROcp3_42*qdd[10]+qd[10]*(OMcp3_23*S2-OMcp3_33*ROcp3_52)
    OPcp3_24 = OPcp3_23+ROcp3_52*qdd[10]+qd[10]*(-OMcp3_13*S2+OMcp3_33*ROcp3_42)
    OPcp3_34 = OPcp3_33+qdd[10]*S2+qd[10]*(OMcp3_13*ROcp3_52-OMcp3_23*ROcp3_42)
    ACcp3_14 = OMcp3_23*ORcp3_34-OMcp3_33*ORcp3_24+OPcp3_23*RLcp3_34-OPcp3_33*RLcp3_24
    ACcp3_24 = -OMcp3_13*ORcp3_34+OMcp3_33*ORcp3_14-OPcp3_13*RLcp3_34+OPcp3_33*RLcp3_14
    ACcp3_34 = OMcp3_13*ORcp3_24-OMcp3_23*ORcp3_14+OPcp3_13*RLcp3_24-OPcp3_23*RLcp3_14
    RLcp3_15 = ROcp3_710*s.dpt[3,15]
    RLcp3_25 = ROcp3_810*s.dpt[3,15]
    RLcp3_35 = ROcp3_910*s.dpt[3,15]
    POcp3_15 = RLcp3_14+RLcp3_15
    POcp3_25 = RLcp3_24+RLcp3_25
    POcp3_35 = RLcp3_34+RLcp3_35
    OMcp3_15 = OMcp3_14+ROcp3_110*qd[11]
    OMcp3_25 = OMcp3_24+ROcp3_210*qd[11]
    OMcp3_35 = OMcp3_34+ROcp3_310*qd[11]
    ORcp3_15 = OMcp3_24*RLcp3_35-OMcp3_34*RLcp3_25
    ORcp3_25 = -OMcp3_14*RLcp3_35+OMcp3_34*RLcp3_15
    ORcp3_35 = OMcp3_14*RLcp3_25-OMcp3_24*RLcp3_15
    VIcp3_15 = ORcp3_14+ORcp3_15
    VIcp3_25 = ORcp3_24+ORcp3_25
    VIcp3_35 = ORcp3_34+ORcp3_35
    OPcp3_15 = OPcp3_14+ROcp3_110*qdd[11]+qd[11]*(OMcp3_24*ROcp3_310-OMcp3_34*ROcp3_210)
    OPcp3_25 = OPcp3_24+ROcp3_210*qdd[11]+qd[11]*(-OMcp3_14*ROcp3_310+OMcp3_34*ROcp3_110)
    OPcp3_35 = OPcp3_34+ROcp3_310*qdd[11]+qd[11]*(OMcp3_14*ROcp3_210-OMcp3_24*ROcp3_110)
    ACcp3_15 = ACcp3_14+OMcp3_24*ORcp3_35-OMcp3_34*ORcp3_25+OPcp3_24*RLcp3_35-OPcp3_34*RLcp3_25
    ACcp3_25 = ACcp3_24-OMcp3_14*ORcp3_35+OMcp3_34*ORcp3_15-OPcp3_14*RLcp3_35+OPcp3_34*RLcp3_15
    ACcp3_35 = ACcp3_34+OMcp3_14*ORcp3_25-OMcp3_24*ORcp3_15+OPcp3_14*RLcp3_25-OPcp3_24*RLcp3_15
    OMcp3_16 = OMcp3_15+ROcp3_411*qd[12]
    OMcp3_26 = OMcp3_25+ROcp3_511*qd[12]
    OMcp3_36 = OMcp3_35+ROcp3_611*qd[12]
    OPcp3_16 = OPcp3_15+ROcp3_411*qdd[12]+qd[12]*(OMcp3_25*ROcp3_611-OMcp3_35*ROcp3_511)
    OPcp3_26 = OPcp3_25+ROcp3_511*qdd[12]+qd[12]*(-OMcp3_15*ROcp3_611+OMcp3_35*ROcp3_411)
    OPcp3_36 = OPcp3_35+ROcp3_611*qdd[12]+qd[12]*(OMcp3_15*ROcp3_511-OMcp3_25*ROcp3_411)
    sens.P[1] = POcp3_15
    sens.P[2] = POcp3_25
    sens.P[3] = POcp3_35
    sens.R[1,1] = ROcp3_112
    sens.R[1,2] = ROcp3_212
    sens.R[1,3] = ROcp3_312
    sens.R[2,1] = ROcp3_411
    sens.R[2,2] = ROcp3_511
    sens.R[2,3] = ROcp3_611
    sens.R[3,1] = ROcp3_712
    sens.R[3,2] = ROcp3_812
    sens.R[3,3] = ROcp3_912
    sens.V[1] = VIcp3_15
    sens.V[2] = VIcp3_25
    sens.V[3] = VIcp3_35
    sens.OM[1] = OMcp3_16
    sens.OM[2] = OMcp3_26
    sens.OM[3] = OMcp3_36
    sens.A[1] = ACcp3_15
    sens.A[2] = ACcp3_25
    sens.A[3] = ACcp3_35
    sens.OMP[1] = OPcp3_16
    sens.OMP[2] = OPcp3_26
    sens.OMP[3] = OPcp3_36

  if (isens == 4): 

    ROcp4_42 = -S1*C2
    ROcp4_52 = C1*C2
    ROcp4_72 = S1*S2
    ROcp4_82 = -C1*S2
    ROcp4_13 = -ROcp4_72*S3+C1*C3
    ROcp4_23 = -ROcp4_82*S3+S1*C3
    ROcp4_33 = -C2*S3
    ROcp4_73 = ROcp4_72*C3+C1*S3
    ROcp4_83 = ROcp4_82*C3+S1*S3
    ROcp4_93 = C2*C3
    ROcp4_113 = ROcp4_13*C13-ROcp4_73*S13
    ROcp4_213 = ROcp4_23*C13-ROcp4_83*S13
    ROcp4_313 = ROcp4_33*C13-ROcp4_93*S13
    ROcp4_713 = ROcp4_13*S13+ROcp4_73*C13
    ROcp4_813 = ROcp4_23*S13+ROcp4_83*C13
    ROcp4_913 = ROcp4_33*S13+ROcp4_93*C13
    ROcp4_414 = ROcp4_42*C14+ROcp4_713*S14
    ROcp4_514 = ROcp4_52*C14+ROcp4_813*S14
    ROcp4_614 = ROcp4_913*S14+C14*S2
    ROcp4_714 = -ROcp4_42*S14+ROcp4_713*C14
    ROcp4_814 = -ROcp4_52*S14+ROcp4_813*C14
    ROcp4_914 = ROcp4_913*C14-S14*S2
    ROcp4_115 = ROcp4_113*C15-ROcp4_714*S15
    ROcp4_215 = ROcp4_213*C15-ROcp4_814*S15
    ROcp4_315 = ROcp4_313*C15-ROcp4_914*S15
    ROcp4_715 = ROcp4_113*S15+ROcp4_714*C15
    ROcp4_815 = ROcp4_213*S15+ROcp4_814*C15
    ROcp4_915 = ROcp4_313*S15+ROcp4_914*C15
    OMcp4_12 = qd[2]*C1
    OMcp4_22 = qd[2]*S1
    OPcp4_12 = qdd[2]*C1-qd[1]*qd[2]*S1
    OPcp4_22 = qdd[2]*S1+qd[1]*qd[2]*C1
    OMcp4_13 = OMcp4_12+ROcp4_42*qd[3]
    OMcp4_23 = OMcp4_22+ROcp4_52*qd[3]
    OMcp4_33 = qd[1]+qd[3]*S2
    OPcp4_13 = OPcp4_12+ROcp4_42*qdd[3]+qd[3]*(OMcp4_22*S2-ROcp4_52*qd[1])
    OPcp4_23 = OPcp4_22+ROcp4_52*qdd[3]+qd[3]*(-OMcp4_12*S2+ROcp4_42*qd[1])
    OPcp4_33 = qdd[1]+qdd[3]*S2+qd[3]*(OMcp4_12*ROcp4_52-OMcp4_22*ROcp4_42)
    RLcp4_14 = ROcp4_42*s.dpt[2,4]+ROcp4_73*s.dpt[3,4]
    RLcp4_24 = ROcp4_52*s.dpt[2,4]+ROcp4_83*s.dpt[3,4]
    RLcp4_34 = ROcp4_93*s.dpt[3,4]+s.dpt[2,4]*S2
    OMcp4_14 = OMcp4_13+ROcp4_42*qd[13]
    OMcp4_24 = OMcp4_23+ROcp4_52*qd[13]
    OMcp4_34 = OMcp4_33+qd[13]*S2
    ORcp4_14 = OMcp4_23*RLcp4_34-OMcp4_33*RLcp4_24
    ORcp4_24 = -OMcp4_13*RLcp4_34+OMcp4_33*RLcp4_14
    ORcp4_34 = OMcp4_13*RLcp4_24-OMcp4_23*RLcp4_14
    OPcp4_14 = OPcp4_13+ROcp4_42*qdd[13]+qd[13]*(OMcp4_23*S2-OMcp4_33*ROcp4_52)
    OPcp4_24 = OPcp4_23+ROcp4_52*qdd[13]+qd[13]*(-OMcp4_13*S2+OMcp4_33*ROcp4_42)
    OPcp4_34 = OPcp4_33+qdd[13]*S2+qd[13]*(OMcp4_13*ROcp4_52-OMcp4_23*ROcp4_42)
    ACcp4_14 = OMcp4_23*ORcp4_34-OMcp4_33*ORcp4_24+OPcp4_23*RLcp4_34-OPcp4_33*RLcp4_24
    ACcp4_24 = -OMcp4_13*ORcp4_34+OMcp4_33*ORcp4_14-OPcp4_13*RLcp4_34+OPcp4_33*RLcp4_14
    ACcp4_34 = OMcp4_13*ORcp4_24-OMcp4_23*ORcp4_14+OPcp4_13*RLcp4_24-OPcp4_23*RLcp4_14
    RLcp4_15 = ROcp4_713*s.dpt[3,18]
    RLcp4_25 = ROcp4_813*s.dpt[3,18]
    RLcp4_35 = ROcp4_913*s.dpt[3,18]
    POcp4_15 = RLcp4_14+RLcp4_15
    POcp4_25 = RLcp4_24+RLcp4_25
    POcp4_35 = RLcp4_34+RLcp4_35
    OMcp4_15 = OMcp4_14+ROcp4_113*qd[14]
    OMcp4_25 = OMcp4_24+ROcp4_213*qd[14]
    OMcp4_35 = OMcp4_34+ROcp4_313*qd[14]
    ORcp4_15 = OMcp4_24*RLcp4_35-OMcp4_34*RLcp4_25
    ORcp4_25 = -OMcp4_14*RLcp4_35+OMcp4_34*RLcp4_15
    ORcp4_35 = OMcp4_14*RLcp4_25-OMcp4_24*RLcp4_15
    VIcp4_15 = ORcp4_14+ORcp4_15
    VIcp4_25 = ORcp4_24+ORcp4_25
    VIcp4_35 = ORcp4_34+ORcp4_35
    OPcp4_15 = OPcp4_14+ROcp4_113*qdd[14]+qd[14]*(OMcp4_24*ROcp4_313-OMcp4_34*ROcp4_213)
    OPcp4_25 = OPcp4_24+ROcp4_213*qdd[14]+qd[14]*(-OMcp4_14*ROcp4_313+OMcp4_34*ROcp4_113)
    OPcp4_35 = OPcp4_34+ROcp4_313*qdd[14]+qd[14]*(OMcp4_14*ROcp4_213-OMcp4_24*ROcp4_113)
    ACcp4_15 = ACcp4_14+OMcp4_24*ORcp4_35-OMcp4_34*ORcp4_25+OPcp4_24*RLcp4_35-OPcp4_34*RLcp4_25
    ACcp4_25 = ACcp4_24-OMcp4_14*ORcp4_35+OMcp4_34*ORcp4_15-OPcp4_14*RLcp4_35+OPcp4_34*RLcp4_15
    ACcp4_35 = ACcp4_34+OMcp4_14*ORcp4_25-OMcp4_24*ORcp4_15+OPcp4_14*RLcp4_25-OPcp4_24*RLcp4_15
    OMcp4_16 = OMcp4_15+ROcp4_414*qd[15]
    OMcp4_26 = OMcp4_25+ROcp4_514*qd[15]
    OMcp4_36 = OMcp4_35+ROcp4_614*qd[15]
    OPcp4_16 = OPcp4_15+ROcp4_414*qdd[15]+qd[15]*(OMcp4_25*ROcp4_614-OMcp4_35*ROcp4_514)
    OPcp4_26 = OPcp4_25+ROcp4_514*qdd[15]+qd[15]*(-OMcp4_15*ROcp4_614+OMcp4_35*ROcp4_414)
    OPcp4_36 = OPcp4_35+ROcp4_614*qdd[15]+qd[15]*(OMcp4_15*ROcp4_514-OMcp4_25*ROcp4_414)
    sens.P[1] = POcp4_15
    sens.P[2] = POcp4_25
    sens.P[3] = POcp4_35
    sens.R[1,1] = ROcp4_115
    sens.R[1,2] = ROcp4_215
    sens.R[1,3] = ROcp4_315
    sens.R[2,1] = ROcp4_414
    sens.R[2,2] = ROcp4_514
    sens.R[2,3] = ROcp4_614
    sens.R[3,1] = ROcp4_715
    sens.R[3,2] = ROcp4_815
    sens.R[3,3] = ROcp4_915
    sens.V[1] = VIcp4_15
    sens.V[2] = VIcp4_25
    sens.V[3] = VIcp4_35
    sens.OM[1] = OMcp4_16
    sens.OM[2] = OMcp4_26
    sens.OM[3] = OMcp4_36
    sens.A[1] = ACcp4_15
    sens.A[2] = ACcp4_25
    sens.A[3] = ACcp4_35
    sens.OMP[1] = OPcp4_16
    sens.OMP[2] = OPcp4_26
    sens.OMP[3] = OPcp4_36

 


# Number of continuation lines = 0


