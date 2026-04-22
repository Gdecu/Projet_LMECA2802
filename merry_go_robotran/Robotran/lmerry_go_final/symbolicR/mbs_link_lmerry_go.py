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
#	==> Generation Date: Wed Apr 22 23:26:05 2026
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

# Number of continuation lines = 0

    print("ERROR : Your symbolic files seem obsolete, i.e. not up-to-date with your MBsysPad model. ")
    print("        Please regenerate your symbolic files (MBsysPad->Tools->Generate Symbolic Files). Exiting. ")
    print("        Error raised in mbs_link.")
    s.flag_stop = 1

