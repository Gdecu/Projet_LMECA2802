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
#	==> Function: F2 - Recursive Inverse Dynamics of tree-like MBS
#
#	==> Git hash: 0e4e6a608eeee06956095d2f2ad315abdb092777
#
##

from math import sin, cos

def invdyna(phi,s,tsim):
    Qq = phi  # compatibility with symbolic generation
    q = s.q
    qd = s.qd
    qdd = s.qdd
 
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

 
# Forward Kinematics

    OM22 = qd[1]*S2
    OM32 = qd[1]*C2
    OMp22 = qd[1]*qd[2]*C2+qdd[1]*S2
    OMp32 = -qd[1]*qd[2]*S2+qdd[1]*C2
    ALPHA22 = -s.g[3]*S2
    ALPHA32 = -s.g[3]*C2
    OM13 = qd[2]*C3-OM32*S3
    OM23 = qd[3]+OM22
    OM33 = qd[2]*S3+OM32*C3
    OMp13 = C3*(qdd[2]-qd[3]*OM32)-S3*(OMp32+qd[2]*qd[3])
    OMp23 = qdd[3]+OMp22
    OMp33 = C3*(OMp32+qd[2]*qd[3])+S3*(qdd[2]-qd[3]*OM32)
    BS13 = -OM23*OM23-OM33*OM33
    BS23 = OM13*OM23
    BS33 = OM13*OM33
    BS53 = -OM13*OM13-OM33*OM33
    BS63 = OM23*OM33
    BS93 = -OM13*OM13-OM23*OM23
    BETA23 = BS23-OMp33
    BETA33 = BS33+OMp23
    BETA43 = BS23+OMp33
    BETA63 = BS63-OMp13
    BETA73 = BS33-OMp23
    BETA83 = BS63+OMp13
    ALPHA13 = -ALPHA32*S3
    ALPHA33 = ALPHA32*C3
    OM14 = OM13*C4-OM33*S4
    OM24 = qd[4]+OM23
    OM34 = OM13*S4+OM33*C4
    OMp14 = C4*(OMp13-qd[4]*OM33)-S4*(OMp33+qd[4]*OM13)
    OMp24 = qdd[4]+OMp23
    OMp34 = C4*(OMp33+qd[4]*OM13)+S4*(OMp13-qd[4]*OM33)
    BS34 = OM14*OM34
    BS64 = OM24*OM34
    BS94 = -OM14*OM14-OM24*OM24
    BETA34 = BS34+OMp24
    BETA64 = BS64-OMp14
    ALPHA14 = C4*(ALPHA13+BETA33*s.dpt[3,1]+BS13*s.dpt[1,1])-S4*(ALPHA33+BETA73*s.dpt[1,1]+BS93*s.dpt[3,1])
    ALPHA24 = ALPHA22+BETA43*s.dpt[1,1]+BETA63*s.dpt[3,1]
    ALPHA34 = C4*(ALPHA33+BETA73*s.dpt[1,1]+BS93*s.dpt[3,1])+S4*(ALPHA13+BETA33*s.dpt[3,1]+BS13*s.dpt[1,1])
    OM15 = OM14*C5-OM34*S5
    OM25 = qd[5]+OM24
    OM35 = OM14*S5+OM34*C5
    OMp15 = C5*(OMp14-qd[5]*OM34)-S5*(OMp34+qd[5]*OM14)
    OMp25 = qdd[5]+OMp24
    OMp35 = C5*(OMp34+qd[5]*OM14)+S5*(OMp14-qd[5]*OM34)
    ALPHA15 = C5*(ALPHA14+BETA34*s.dpt[3,9])-S5*(ALPHA34+BS94*s.dpt[3,9])
    ALPHA25 = ALPHA24+BETA64*s.dpt[3,9]
    ALPHA35 = C5*(ALPHA34+BS94*s.dpt[3,9])+S5*(ALPHA14+BETA34*s.dpt[3,9])
    OM16 = qd[6]+OM15
    OM26 = OM25*C6+OM35*S6
    OM36 = -OM25*S6+OM35*C6
    OMp16 = qdd[6]+OMp15
    OMp26 = C6*(OMp25+qd[6]*OM35)+S6*(OMp35-qd[6]*OM25)
    OMp36 = C6*(OMp35-qd[6]*OM25)-S6*(OMp25+qd[6]*OM35)
    BS36 = OM16*OM36
    BS66 = OM26*OM36
    BS96 = -OM16*OM16-OM26*OM26
    BETA36 = BS36+OMp26
    BETA66 = BS66-OMp16
    ALPHA26 = ALPHA25*C6+ALPHA35*S6
    ALPHA36 = -ALPHA25*S6+ALPHA35*C6
    OM17 = qd[7]+OM13
    OM27 = OM23*C7+OM33*S7
    OM37 = -OM23*S7+OM33*C7
    OMp17 = qdd[7]+OMp13
    OMp27 = C7*(OMp23+qd[7]*OM33)+S7*(OMp33-qd[7]*OM23)
    OMp37 = C7*(OMp33-qd[7]*OM23)-S7*(OMp23+qd[7]*OM33)
    BS37 = OM17*OM37
    BS67 = OM27*OM37
    BS97 = -OM17*OM17-OM27*OM27
    BETA37 = BS37+OMp27
    BETA67 = BS67-OMp17
    ALPHA17 = ALPHA13+BETA23*s.dpt[2,2]+BETA33*s.dpt[3,2]
    ALPHA27 = C7*(ALPHA22+BETA63*s.dpt[3,2]+BS53*s.dpt[2,2])+S7*(ALPHA33+BETA83*s.dpt[2,2]+BS93*s.dpt[3,2])
    ALPHA37 = C7*(ALPHA33+BETA83*s.dpt[2,2]+BS93*s.dpt[3,2])-S7*(ALPHA22+BETA63*s.dpt[3,2]+BS53*s.dpt[2,2])
    OM18 = qd[8]+OM17
    OM28 = OM27*C8+OM37*S8
    OM38 = -OM27*S8+OM37*C8
    OMp18 = qdd[8]+OMp17
    OMp28 = C8*(OMp27+qd[8]*OM37)+S8*(OMp37-qd[8]*OM27)
    OMp38 = C8*(OMp37-qd[8]*OM27)-S8*(OMp27+qd[8]*OM37)
    ALPHA18 = ALPHA17+BETA37*s.dpt[3,12]
    ALPHA28 = C8*(ALPHA27+BETA67*s.dpt[3,12])+S8*(ALPHA37+BS97*s.dpt[3,12])
    ALPHA38 = C8*(ALPHA37+BS97*s.dpt[3,12])-S8*(ALPHA27+BETA67*s.dpt[3,12])
    OM19 = OM18*C9-OM38*S9
    OM29 = qd[9]+OM28
    OM39 = OM18*S9+OM38*C9
    OMp19 = C9*(OMp18-qd[9]*OM38)-S9*(OMp38+qd[9]*OM18)
    OMp29 = qdd[9]+OMp28
    OMp39 = C9*(OMp38+qd[9]*OM18)+S9*(OMp18-qd[9]*OM38)
    BS39 = OM19*OM39
    BS69 = OM29*OM39
    BS99 = -OM19*OM19-OM29*OM29
    BETA39 = BS39+OMp29
    BETA69 = BS69-OMp19
    ALPHA19 = ALPHA18*C9-ALPHA38*S9
    ALPHA39 = ALPHA18*S9+ALPHA38*C9
    OM110 = OM13*C10-OM33*S10
    OM210 = qd[10]+OM23
    OM310 = OM13*S10+OM33*C10
    OMp110 = C10*(OMp13-qd[10]*OM33)-S10*(OMp33+qd[10]*OM13)
    OMp210 = qdd[10]+OMp23
    OMp310 = C10*(OMp33+qd[10]*OM13)+S10*(OMp13-qd[10]*OM33)
    BS310 = OM110*OM310
    BS610 = OM210*OM310
    BS910 = -OM110*OM110-OM210*OM210
    BETA310 = BS310+OMp210
    BETA610 = BS610-OMp110
    ALPHA110 = C10*(ALPHA13+BETA33*s.dpt[3,3]+BS13*s.dpt[1,3])-S10*(ALPHA33+BETA73*s.dpt[1,3]+BS93*s.dpt[3,3])
    ALPHA210 = ALPHA22+BETA43*s.dpt[1,3]+BETA63*s.dpt[3,3]
    ALPHA310 = C10*(ALPHA33+BETA73*s.dpt[1,3]+BS93*s.dpt[3,3])+S10*(ALPHA13+BETA33*s.dpt[3,3]+BS13*s.dpt[1,3])
    OM111 = OM110*C11-OM310*S11
    OM211 = qd[11]+OM210
    OM311 = OM110*S11+OM310*C11
    OMp111 = C11*(OMp110-qd[11]*OM310)-S11*(OMp310+qd[11]*OM110)
    OMp211 = qdd[11]+OMp210
    OMp311 = C11*(OMp310+qd[11]*OM110)+S11*(OMp110-qd[11]*OM310)
    ALPHA111 = C11*(ALPHA110+BETA310*s.dpt[3,15])-S11*(ALPHA310+BS910*s.dpt[3,15])
    ALPHA211 = ALPHA210+BETA610*s.dpt[3,15]
    ALPHA311 = C11*(ALPHA310+BS910*s.dpt[3,15])+S11*(ALPHA110+BETA310*s.dpt[3,15])
    OM112 = qd[12]+OM111
    OM212 = OM211*C12+OM311*S12
    OM312 = -OM211*S12+OM311*C12
    OMp112 = qdd[12]+OMp111
    OMp212 = C12*(OMp211+qd[12]*OM311)+S12*(OMp311-qd[12]*OM211)
    OMp312 = C12*(OMp311-qd[12]*OM211)-S12*(OMp211+qd[12]*OM311)
    BS312 = OM112*OM312
    BS612 = OM212*OM312
    BS912 = -OM112*OM112-OM212*OM212
    BETA312 = BS312+OMp212
    BETA612 = BS612-OMp112
    ALPHA212 = ALPHA211*C12+ALPHA311*S12
    ALPHA312 = -ALPHA211*S12+ALPHA311*C12
    OM113 = qd[13]+OM13
    OM213 = OM23*C13+OM33*S13
    OM313 = -OM23*S13+OM33*C13
    OMp113 = qdd[13]+OMp13
    OMp213 = C13*(OMp23+qd[13]*OM33)+S13*(OMp33-qd[13]*OM23)
    OMp313 = C13*(OMp33-qd[13]*OM23)-S13*(OMp23+qd[13]*OM33)
    BS313 = OM113*OM313
    BS613 = OM213*OM313
    BS913 = -OM113*OM113-OM213*OM213
    BETA313 = BS313+OMp213
    BETA613 = BS613-OMp113
    ALPHA113 = ALPHA13+BETA23*s.dpt[2,4]+BETA33*s.dpt[3,4]
    ALPHA213 = C13*(ALPHA22+BETA63*s.dpt[3,4]+BS53*s.dpt[2,4])+S13*(ALPHA33+BETA83*s.dpt[2,4]+BS93*s.dpt[3,4])
    ALPHA313 = C13*(ALPHA33+BETA83*s.dpt[2,4]+BS93*s.dpt[3,4])-S13*(ALPHA22+BETA63*s.dpt[3,4]+BS53*s.dpt[2,4])
    OM114 = qd[14]+OM113
    OM214 = OM213*C14+OM313*S14
    OM314 = -OM213*S14+OM313*C14
    OMp114 = qdd[14]+OMp113
    OMp214 = C14*(OMp213+qd[14]*OM313)+S14*(OMp313-qd[14]*OM213)
    OMp314 = C14*(OMp313-qd[14]*OM213)-S14*(OMp213+qd[14]*OM313)
    ALPHA114 = ALPHA113+BETA313*s.dpt[3,18]
    ALPHA214 = C14*(ALPHA213+BETA613*s.dpt[3,18])+S14*(ALPHA313+BS913*s.dpt[3,18])
    ALPHA314 = C14*(ALPHA313+BS913*s.dpt[3,18])-S14*(ALPHA213+BETA613*s.dpt[3,18])
    OM115 = OM114*C15-OM314*S15
    OM215 = qd[15]+OM214
    OM315 = OM114*S15+OM314*C15
    OMp115 = C15*(OMp114-qd[15]*OM314)-S15*(OMp314+qd[15]*OM114)
    OMp215 = qdd[15]+OMp214
    OMp315 = C15*(OMp314+qd[15]*OM114)+S15*(OMp114-qd[15]*OM314)
    BS315 = OM115*OM315
    BS615 = OM215*OM315
    BS915 = -OM115*OM115-OM215*OM215
    BETA315 = BS315+OMp215
    BETA615 = BS615-OMp115
    ALPHA115 = ALPHA114*C15-ALPHA314*S15
    ALPHA315 = ALPHA114*S15+ALPHA314*C15
 
# Backward Dynamics

    Fs115 = -s.frc[1,15]+s.m[15]*(ALPHA115+BETA315*s.l[3,15])
    Fs215 = -s.frc[2,15]+s.m[15]*(ALPHA214+BETA615*s.l[3,15])
    Fs315 = -s.frc[3,15]+s.m[15]*(ALPHA315+BS915*s.l[3,15])
    Cq115 = -s.trq[1,15]+s.In[1,15]*OMp115-s.In[5,15]*OM215*OM315+s.In[9,15]*OM215*OM315-Fs215*s.l[3,15]
    Cq215 = -s.trq[2,15]+s.In[1,15]*OM115*OM315+s.In[5,15]*OMp215-s.In[9,15]*OM115*OM315+Fs115*s.l[3,15]
    Cq315 = -s.trq[3,15]-s.In[1,15]*OM115*OM215+s.In[5,15]*OM115*OM215+s.In[9,15]*OMp315
    Fq114 = Fs115*C15+Fs315*S15
    Fq314 = -Fs115*S15+Fs315*C15
    Cq114 = Cq115*C15+Cq315*S15
    Cq314 = -Cq115*S15+Cq315*C15
    Fs113 = -s.frc[1,13]+s.m[13]*(ALPHA113+BETA313*s.l[3,13])
    Fs213 = -s.frc[2,13]+s.m[13]*(ALPHA213+BETA613*s.l[3,13])
    Fs313 = -s.frc[3,13]+s.m[13]*(ALPHA313+BS913*s.l[3,13])
    Fq113 = Fq114+Fs113
    Fq213 = Fs213-Fq314*S14+Fs215*C14
    Fq313 = Fs313+Fq314*C14+Fs215*S14
    Cq113 = -s.trq[1,13]+Cq114+s.In[1,13]*OMp113-s.In[5,13]*OM213*OM313+s.In[9,13]*OM213*OM313-Fs213*s.l[3,13]- \
 	  s.dpt[3,18]*(-Fq314*S14+Fs215*C14)
    Cq213 = -s.trq[2,13]+s.In[1,13]*OM113*OM313+s.In[5,13]*OMp213-s.In[9,13]*OM113*OM313+Cq215*C14-Cq314*S14+Fq114* \
 	  s.dpt[3,18]+Fs113*s.l[3,13]
    Cq313 = -s.trq[3,13]-s.In[1,13]*OM113*OM213+s.In[5,13]*OM113*OM213+s.In[9,13]*OMp313+Cq215*S14+Cq314*C14
    Fs112 = -s.frc[1,12]+s.m[12]*(ALPHA111+BETA312*s.l[3,12])
    Fs212 = -s.frc[2,12]+s.m[12]*(ALPHA212+BETA612*s.l[3,12])
    Fs312 = -s.frc[3,12]+s.m[12]*(ALPHA312+BS912*s.l[3,12])
    Cq112 = -s.trq[1,12]+s.In[1,12]*OMp112-s.In[5,12]*OM212*OM312+s.In[9,12]*OM212*OM312-Fs212*s.l[3,12]
    Cq212 = -s.trq[2,12]+s.In[1,12]*OM112*OM312+s.In[5,12]*OMp212-s.In[9,12]*OM112*OM312+Fs112*s.l[3,12]
    Cq312 = -s.trq[3,12]-s.In[1,12]*OM112*OM212+s.In[5,12]*OM112*OM212+s.In[9,12]*OMp312
    Fq211 = Fs212*C12-Fs312*S12
    Fq311 = Fs212*S12+Fs312*C12
    Cq211 = Cq212*C12-Cq312*S12
    Cq311 = Cq212*S12+Cq312*C12
    Fs110 = -s.frc[1,10]+s.m[10]*(ALPHA110+BETA310*s.l[3,10])
    Fs210 = -s.frc[2,10]+s.m[10]*(ALPHA210+BETA610*s.l[3,10])
    Fs310 = -s.frc[3,10]+s.m[10]*(ALPHA310+BS910*s.l[3,10])
    Fq110 = Fs110+Fq311*S11+Fs112*C11
    Fq210 = Fq211+Fs210
    Fq310 = Fs310+Fq311*C11-Fs112*S11
    Cq110 = -s.trq[1,10]+s.In[1,10]*OMp110-s.In[5,10]*OM210*OM310+s.In[9,10]*OM210*OM310+Cq112*C11+Cq311*S11-Fq211* \
 	  s.dpt[3,15]-Fs210*s.l[3,10]
    Cq210 = -s.trq[2,10]+Cq211+s.In[1,10]*OM110*OM310+s.In[5,10]*OMp210-s.In[9,10]*OM110*OM310+Fs110*s.l[3,10]+ \
 	  s.dpt[3,15]*(Fq311*S11+Fs112*C11)
    Cq310 = -s.trq[3,10]-s.In[1,10]*OM110*OM210+s.In[5,10]*OM110*OM210+s.In[9,10]*OMp310-Cq112*S11+Cq311*C11
    Fs19 = -s.frc[1,9]+s.m[9]*(ALPHA19+BETA39*s.l[3,9])
    Fs29 = -s.frc[2,9]+s.m[9]*(ALPHA28+BETA69*s.l[3,9])
    Fs39 = -s.frc[3,9]+s.m[9]*(ALPHA39+BS99*s.l[3,9])
    Cq19 = -s.trq[1,9]+s.In[1,9]*OMp19-s.In[5,9]*OM29*OM39+s.In[9,9]*OM29*OM39-Fs29*s.l[3,9]
    Cq29 = -s.trq[2,9]+s.In[1,9]*OM19*OM39+s.In[5,9]*OMp29-s.In[9,9]*OM19*OM39+Fs19*s.l[3,9]
    Cq39 = -s.trq[3,9]-s.In[1,9]*OM19*OM29+s.In[5,9]*OM19*OM29+s.In[9,9]*OMp39
    Fq18 = Fs19*C9+Fs39*S9
    Fq38 = -Fs19*S9+Fs39*C9
    Cq18 = Cq19*C9+Cq39*S9
    Cq38 = -Cq19*S9+Cq39*C9
    Fs17 = -s.frc[1,7]+s.m[7]*(ALPHA17+BETA37*s.l[3,7])
    Fs27 = -s.frc[2,7]+s.m[7]*(ALPHA27+BETA67*s.l[3,7])
    Fs37 = -s.frc[3,7]+s.m[7]*(ALPHA37+BS97*s.l[3,7])
    Fq17 = Fq18+Fs17
    Fq27 = Fs27-Fq38*S8+Fs29*C8
    Fq37 = Fs37+Fq38*C8+Fs29*S8
    Cq17 = -s.trq[1,7]+Cq18+s.In[1,7]*OMp17-s.In[5,7]*OM27*OM37+s.In[9,7]*OM27*OM37-Fs27*s.l[3,7]-s.dpt[3,12]*(-Fq38* \
 	  S8+Fs29*C8)
    Cq27 = -s.trq[2,7]+s.In[1,7]*OM17*OM37+s.In[5,7]*OMp27-s.In[9,7]*OM17*OM37+Cq29*C8-Cq38*S8+Fq18*s.dpt[3,12]+Fs17* \
 	  s.l[3,7]
    Cq37 = -s.trq[3,7]-s.In[1,7]*OM17*OM27+s.In[5,7]*OM17*OM27+s.In[9,7]*OMp37+Cq29*S8+Cq38*C8
    Fs16 = -s.frc[1,6]+s.m[6]*(ALPHA15+BETA36*s.l[3,6])
    Fs26 = -s.frc[2,6]+s.m[6]*(ALPHA26+BETA66*s.l[3,6])
    Fs36 = -s.frc[3,6]+s.m[6]*(ALPHA36+BS96*s.l[3,6])
    Cq16 = -s.trq[1,6]+s.In[1,6]*OMp16-s.In[5,6]*OM26*OM36+s.In[9,6]*OM26*OM36-Fs26*s.l[3,6]
    Cq26 = -s.trq[2,6]+s.In[1,6]*OM16*OM36+s.In[5,6]*OMp26-s.In[9,6]*OM16*OM36+Fs16*s.l[3,6]
    Cq36 = -s.trq[3,6]-s.In[1,6]*OM16*OM26+s.In[5,6]*OM16*OM26+s.In[9,6]*OMp36
    Fq25 = Fs26*C6-Fs36*S6
    Fq35 = Fs26*S6+Fs36*C6
    Cq25 = Cq26*C6-Cq36*S6
    Cq35 = Cq26*S6+Cq36*C6
    Fs14 = -s.frc[1,4]+s.m[4]*(ALPHA14+BETA34*s.l[3,4])
    Fs24 = -s.frc[2,4]+s.m[4]*(ALPHA24+BETA64*s.l[3,4])
    Fs34 = -s.frc[3,4]+s.m[4]*(ALPHA34+BS94*s.l[3,4])
    Fq14 = Fs14+Fq35*S5+Fs16*C5
    Fq24 = Fq25+Fs24
    Fq34 = Fs34+Fq35*C5-Fs16*S5
    Cq14 = -s.trq[1,4]+s.In[1,4]*OMp14-s.In[5,4]*OM24*OM34+s.In[9,4]*OM24*OM34+Cq16*C5+Cq35*S5-Fq25*s.dpt[3,9]-Fs24* \
 	  s.l[3,4]
    Cq24 = -s.trq[2,4]+Cq25+s.In[1,4]*OM14*OM34+s.In[5,4]*OMp24-s.In[9,4]*OM14*OM34+Fs14*s.l[3,4]+s.dpt[3,9]*(Fq35*S5 \
 	  +Fs16*C5)
    Cq34 = -s.trq[3,4]-s.In[1,4]*OM14*OM24+s.In[5,4]*OM14*OM24+s.In[9,4]*OMp34-Cq16*S5+Cq35*C5
    Fs13 = -s.frc[1,3]+s.m[3]*(ALPHA13+BETA33*s.l[3,3])
    Fs23 = -s.frc[2,3]+s.m[3]*(ALPHA22+BETA63*s.l[3,3])
    Cq13 = -s.trq[1,3]+Cq113+Cq17+s.In[1,3]*OMp13-s.In[5,3]*OM23*OM33+s.In[9,3]*OM23*OM33+Cq110*C10+Cq14*C4+Cq310*S10 \
 	  +Cq34*S4-Fq210*s.dpt[3,3]-Fq24*s.dpt[3,1]-Fs23*s.l[3,3]+s.dpt[2,2]*(Fq27*S7+Fq37*C7)+s.dpt[2,4]*(Fq213*S13+Fq313*C13)- \
 	  s.dpt[3,2]*(Fq27*C7-Fq37*S7)-s.dpt[3,4]*(Fq213*C13-Fq313*S13)
    Cq23 = -s.trq[2,3]+Cq210+Cq24+s.In[1,3]*OM13*OM33+s.In[5,3]*OMp23-s.In[9,3]*OM13*OM33+Cq213*C13+Cq27*C7-Cq313*S13 \
 	  -Cq37*S7+Fq113*s.dpt[3,4]+Fq17*s.dpt[3,2]+Fs13*s.l[3,3]-s.dpt[1,1]*(-Fq14*S4+Fq34*C4)-s.dpt[1,3]*(-Fq110*S10+Fq310*C10 \
 	  )+s.dpt[3,1]*(Fq14*C4+Fq34*S4)+s.dpt[3,3]*(Fq110*C10+Fq310*S10)
    Cq33 = -s.trq[3,3]-s.In[1,3]*OM13*OM23+s.In[5,3]*OM13*OM23+s.In[9,3]*OMp33-Cq110*S10-Cq14*S4+Cq213*S13+Cq27*S7+ \
 	  Cq310*C10+Cq313*C13+Cq34*C4+Cq37*C7-Fq113*s.dpt[2,4]-Fq17*s.dpt[2,2]+Fq210*s.dpt[1,3]+Fq24*s.dpt[1,1]
    Cq12 = Cq13*C3+Cq33*S3
    Cq32 = -Cq13*S3+Cq33*C3
    Cq31 = Cq23*S2+Cq32*C2
 
# Symbolic model output

    Qq[1] = Cq31
    Qq[2] = Cq12
    Qq[3] = Cq23
    Qq[4] = Cq24
    Qq[5] = Cq25
    Qq[6] = Cq16
    Qq[7] = Cq17
    Qq[8] = Cq18
    Qq[9] = Cq29
    Qq[10] = Cq210
    Qq[11] = Cq211
    Qq[12] = Cq112
    Qq[13] = Cq113
    Qq[14] = Cq114
    Qq[15] = Cq215

# Number of continuation lines = 2


