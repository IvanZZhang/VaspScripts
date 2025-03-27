#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 11:55:30 2019

@author: Yaqiong Su

Rewritten on Mar 21 2025, with Python3.7, By Yifan Zhang
"""

import math
import datetime
import time
from scipy import constants as cont

######   timing   ######
start = time.time()
print('********** thermodynamic analysis from Yaqiong Su Eindhoven **********')
print('is performing thermodynamic analysis')
### current time ###
start_time = datetime.datetime.now()
print("Start time:\t" + start_time.strftime('%Y.%m.%d-%H:%M:%S'))   #strftime可以自定义时间的输出格式

with open('OUTCAR','rb') as f1:
    lines = f1.readlines()

with open('Frequency','wb') as f2:
    keyword = b'cm-1'
    for idx, line in enumerate(lines):
        if keyword in line:
            print(line)
            f2.write(line)

#Zero-Point Energy Correction (ZPE)
Evib = []
with open ('Frequency','rb') as Fre:
    for line in Fre.readlines():
        if b'f/i' not in line:
            E_vib = float(line[65:75])
            #print(E_vib)
            Evib.append(E_vib)
print('Vibrational Energy of each mode (meV):\t{}'.format(Evib))
ZPE = 0.5 * sum(Evib) / 1000   # unit: eV
print('ZPE Value: %12.10f eV' % ZPE)


#Calcultion of vibrational entropy
#Ref: Atkin's Physical Chemistry,10th edition.(15E.13a)
#Ref: http://210.45.168.34:8080/elite/wlhx/jiaocai/C_05.htm
#Ref: https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html
#Ref: https://www.bigbrosci.com/2018/11/07/ex69/
# script, nu = sys.argv # in this script, we do not use the outside input

# grep frequency
nu = []
with open ('Frequency','rb') as Fre:
    for line in Fre.readlines():
        if b'f/i' not in line:
            nu_f = float(line[46:57])
            #print(nu_f)
            nu.append(nu_f)
print('Vibrational Frequency of each mode (cm-1):\t{}'.format(nu))
#nu = float(nu)*100 # convert cm-1 to m-1 for frequency of each vibrational mode
h_p = cont.h    # 6.62606957E-34 (J*S) Plank constant
k_b = cont.k    # 1.38064852E-23 (m2*kg*s-2*K-1) Boltzman constant
R_gas = cont.R  # 8.3144598 (J*mol-1*k-1) Gas constant
c_ls = cont.c   # 299792458 (m*s-1) # light speed
Tem = 300       # Temperature 300 K
beta = 1 / (k_b * Tem)

def get_svib(nui):  # to get vibrational entropy
    x_i = h_p * nui * 100 * c_ls * beta     # convert cm-1 to m-1
    vpf_l = x_i / (math.exp(x_i) - 1)   # left part in brace in 15E.13a
    vpf_r = math.log(1 - math.exp(-x_i))  # right part in brace in 15E.13a
    vpf = vpf_l - vpf_r   # vibrational partition function
    S_vib = R_gas * vpf   # vibrational entropy of each mode
    return S_vib

#S_vib = get_vpf(nui)   # J*k-1*mol-1
#TS_vib = S_vib*Tem/1000/96.483   # eV  vibrational entropy to energy
#print ('%0.4f\t %.4f') %(S_vib, TS_vib)
S_vibn, TS_vibn = [], []

for i in range(len(nu)):
    S_vib = get_svib(nu[i])
    TS_vib = S_vib * Tem / 1000 / 96.483
    S_vibn.append(S_vib), TS_vibn.append(TS_vib)
    #print(S_vib, TS_vib)
    #print(S_vibn)   # group of vibrational entropy of each mode
print('Contribution of Vibrational Entropy of each mode (eV):\t{}'.format(TS_vibn)) # group of vibrational entropy energy
Total_TS_vibn = sum(TS_vibn)
print('Total Contribution of Vibrational Entropy (eV):\t{}'.format(Total_TS_vibn))  # Total Contribution of Vibrational Entropy (eV)


#calculation of translation entropy of adsorbates 2D-approximation
#Ref: Computational Catalysis (RSC), P32, E1.61
pi = cont.pi  # 3.141592653589793    pai contant
N0 = cont.Avogadro # 6.02214076×1E23 mol-1   Avogadro constant
Cs = 1E19   # sites*m-2     active sites per unit area
m_CO = 28.0101  # g*mol-1   molar mass of CO
m = (m_CO / 1000) / N0   # kg   molecular mass
S_trans = k_b * (math.log((2 * pi * m * k_b * Tem) / (Cs * (h_p ** 2)) + 2)) # J*k-1*mol-1
TS_trans = S_trans * Tem / 1000 / 96.483   # transitional entropty energy per molecule
TS_trans_mol = TS_trans * N0   # eV*mol-1   transitional contribution per mol
#print('Transitonal Entropy Contribution of surface adosrbates (eV):\t{}'.format(TS_trans_mol))


#calculation of heat capacity (vibrational)
#Ref: C.J. Cramer. Essentials of Computational Chemistry, 2nd Edition. Wiley, 2004.
#Ref: https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html
ivhc0 = []
for i in Evib:
    ivhc_mode = (i / 1000) / (math.exp((i * 96.483) / (R_gas * Tem)) - 1)   # integral of heat capacity of each vibrational mode
    ivhc0.append(ivhc_mode)
print('integral of heat capacity of each vibrational mode (eV):\t{}'.format(ivhc0))
ivhc = sum(ivhc0)   # sum of integral of heat capacity of all vibrational modes
print('Contribution of vibrational heat capacity (eV):\t{}'.format(ivhc))


# calculation of vibrational partition function of adsorbates
# Ref: Physical Chemistry, Fu Xiancai
# Ref: http://chem.xmu.edu.cn/teach/chemistry-net-teaching/wuhua/chapter7/part5/5-4.html
qv0 = []
for i in Evib:
    qvi = 1 / (1 - math.exp((-i * 96.483)/(R_gas * Tem)))   # vibrational PF of each mode
    qv0.append(qvi)
#print('Vibrational Partition Function of each mode:\t{}'.format(qv0))
#qv = reduce(operator.mul,qv0)
#print('Vibrational Partition Function:\t' + qv)


# calculation of transitional partition function gaseous molecule
# Ref: Physical Chemistry, Fu Xiancai
# Ref: http://chem.xmu.edu.cn/teach/chemistry-net-teaching/wuhua/chapter7/part5/5-2.html
P_molecule = 1  # unit atm
P0 = 1  # unit atom   one atmospheric (barometric) pressure
m_CH4 = 16.04   # g*mol-1   molar mass of CH4
m_CO2 = 44.0095 # g*mol-1   molar mass of CO2
m_H2O = 18.01524    # g*mol-1   molar mass of H2O
m_H2 = 2.01588  # g*mol-1   molar mass of H2O
V_molecule = (k_b * Tem) / (P_molecule / P0)
qt_CH4 = V_molecule * (2 * pi * (m_CH4/1000/N0) * k_b * Tem / (h_p ** 2)) ** (float(3) / float(2))
qt_CO = V_molecule * (2 * pi * (m_CO/1000/N0) * k_b * Tem / (h_p ** 2)) ** (float(3) / float(2))
qt_CO2 = V_molecule * (2 * pi * (m_CO2/1000/N0) * k_b * Tem / (h_p ** 2)) ** (float(3) / float(2))
qt_H2O = V_molecule * (2 * pi * (m_H2O/1000/N0) * k_b * Tem / (h_p ** 2)) ** (float(3) / float(2))
qt_H2 = V_molecule * (2 * pi * (m_H2/1000/N0) * k_b * Tem / (h_p ** 2)) ** (float(3) / float(2))
#print('Transitional Partition Function of gaseous CH4 at 1 atm:\t'+qt_CH4)
#print('Transitional Partition Function of gaseous CO at 1 atm:\t'+qt_CO)
#print('Transitional Partition Function of gaseous CO2 at 1 atm:\t'+qt_CO2)
#print('Transitional Partition Function of gaseous H2O at 1 atm:\t'+qt_H2O)
#print('Transitional Partition Function of gaseous H2 at 1 atm:\t'+qt_H2)


# calculation of transitional entropy of gaseous molecule
# Ref: Physical Chemistry, Fu Xiancai
# Ref: http://chem.xmu.edu.cn/teach/chemistry-net-teaching/wuhua/chapter7/part5/5-2.html
# Ref: Sackur-Tetrode formula
St_CH4 = R_gas * (math.log(qt_CH4 / N0) + float(5) / float(2)) # molar transitional entropy
#print('Molar Transitional Entropy of gaseous CH4 at 1 atm:\t'+St_CH4)
TSt_CH4 = St_CH4*Tem/1000/96.483   # eV
#print('Molar Transitional Entropy contribution of gaseous CH4 at 1 atm (eV):\t'+TSt_CH4)


# calculation of rotational entropy of gaseous molecule
# Ref: Physical Chemistry, Fu Xiancai
# Ref: http://chem.xmu.edu.cn/teach/chemistry-net-teaching/wuhua/chapter7/part5/5-2.html
# Ref: https://pubs.acs.org/doi/suppl/10.1021/acscatal.7b03295/suppl_file/cs7b03295_si_001.pdf
delta_CH4 = 0.655 * 1E-3   # eV    rotational constant of gaseous CH4
qr_CH4 = (32 * pi ** 2 / 3) * ((R_gas * Tem / 1000 / 96.483) / (4 * pi * delta_CH4)) ** (float(3) / float(2))
#print('Rotational Partition Function of gaseous CH4:\t'+qr_CH4)


##########   timing   #############
stop = time.time()
print("running time:\t" + str(stop-start) + " seconds")
terminal_time = datetime.datetime.now()
print("Terminal time:\t" + terminal_time.strftime('%Y.%m.%d-%H:%M:%S'))   #strftime可以自定义时间的输出格式
print('thermodynamic analysis has been done well')
print('********** thermodynamic analysis from Yaqiong Su Eindhoven **********')
