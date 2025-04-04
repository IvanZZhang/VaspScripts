#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
For python3 version
Ported from python2 to python3 by Yanyang Qin
Copyright: Yanyang Qin. XJTU
@author: Yanyang Qin, XJTU
Rewritten: Yifan Zhang, XJTU
"""
#  Copyright (c) 2025. By ZYF

######   split the DOSCAR to each atom   ######
import numpy as np
import datetime
import time
import sys
import linecache as lce
import pandas as pd
from scipy.integrate import simps
import os
import glob

######   get input   ######
np.set_printoptions(suppress=True)
try:
    atom_ini, atom_final = sys.argv[1], sys.argv[2]   # pass parameters  atom range
    argv_input = True
except IndexError:
    argv_input = False

######   timing   ######
start = time.time()
print('\n********** splitted dos and d-band center from Yaqiong Su Eindhoven **********\n')
print('is getting splitted dos and d-band center')
### current time ###
start_time = datetime.datetime.now()
print("Start time:\t" + start_time.strftime('%Y.%m.%d-%H:%M:%S'))   #strftime可以自定义时间的输出格式

######   reading DOSCAR   ######
line1 = lce.getline('DOSCAR', 1)
atoms_number = int(line1.split()[0])
line6 = lce.getline('DOSCAR', 6)
E_fermi = float(line6.split()[3])
epoints = int(line6.split()[2])
linec = lce.getline('DOSCAR', 6 + epoints + 1 + 1)
column_num = int(len(linec.split()))
print("atoms number:\t", atoms_number)
print("Fermi energy:\t", E_fermi)
print("Energy points:\t", epoints)
print("column number:\t", column_num)

if not argv_input:
    atom_ini, atom_final = 1, atoms_number

f2 = open('total_dos.dat', 'wb')
f3 = open('atoms_dos0.dat', 'wb')
f4 = open('atoms_dos.dat', 'wb')
with open('DOSCAR', 'rb') as f:
    i = 0
    while True:
        i += 1
        line = f.readline()
        if 6 < i < 6+epoints+1:
            f2.write(line)
        if 6+epoints < i < 1e5:
            f3.write(line)
            if b' 301  ' in line:
                continue
            f4.write(line)
        if i > 1e5 or line == b'':
            print('DOSCAR reading finished, {} lines in total.'.format(i))
            break
f2.close()
f3.close()
f4.close()

E_corrected = np.array([float(l.split()[0]) for l in open('total_dos.dat', 'rb')]) - E_fermi
E_corrected = E_corrected.reshape(epoints, 1)
rows = np.zeros([epoints,column_num-1])
E_f = np.zeros([epoints, 1]) - E_fermi
E_f0 = np.hstack((E_f, rows))
ones = np.ones([epoints,column_num])
for i in range(2,column_num,2):
    ones[:,i] = -1   # construct 1 -1 matrix

names = [ 's', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2-r2', 'd_xz', 'd_x2-y2' ]
all_names = []
for n in names:
    all_names.extend( [ '{}_up'.format(n), '{}_down'.format(n) ] )
all_names.insert( 0, 'energy(eV)' )

names0 = [ 's', 'p', 'd' ]
group_names = []
for n in names0:
    group_names.extend( [ '{}_up'.format(n), '{}_down'.format(n) ] )
group_names.insert( 0, 'energy(eV)' )
#print (all_names)

atoms_dos = np.loadtxt('atoms_dos.dat')
dos = []
path = 'DOS_data'
os.makedirs(path, exist_ok=True)
for i in range(1,atoms_number+1):
    dos.append('DOS' + str(i))
    dos[i-1] = np.vsplit(atoms_dos,atoms_number)[i-1].astype(float) # vertical split
    dos[i-1] = np.add(dos[i-1], E_f0) # corrected energy by E_fermi
    dos[i-1] = np.multiply(dos[i-1],ones)
    locals()['DOS' + str(i)] = dos[i-1]
    file = path + '/DOS' + str(i) + '.dat'
    np.savetxt(file, dos[i-1])   # projected to each atomic orbital orientation
    #fi = 'f' + str(i)
    fi = np.loadtxt(file)
    #dfi = 'df' + str(i)
    dfi = pd.DataFrame(fi, columns = all_names)
    dfi.set_index('energy(eV)',inplace=True)
    file_splitted = path + '/DOS' + str(i) + '_splitted' + '.dat'
    dfi.to_csv(file_splitted, sep='\t', float_format='%12.10f') # projected to each atomic orbital orientation

######   sum of orbital dos   ######
    #Eni = 'En' + str(i)
    Eni = dos[i-1][:,0].reshape(epoints,1)
    locals()['En' + str(i)] = Eni
    #s_upi = 's_up' + str(i)
    s_upi = dos[i-1][:,1].reshape(epoints,1)
    locals()['s_up' + str(i)] = s_upi
    #s_downi = 's_down' + str(i)
    s_downi = dos[i-1][:,2].reshape(epoints,1)
    locals()['s_down' + str(i)] = s_downi
    #p_upi = 'p_up' + str(i)
    p_upi = (dos[i-1][:,3]+dos[i-1][:,5]+dos[i-1][:,7]).reshape(epoints,1)
    locals()['p_up' + str(i)] = p_upi
    #p_downi = 'p_down' + str(i)
    p_downi = (dos[i-1][:,4]+dos[i-1][:,6]+dos[i-1][:,8]).reshape(epoints,1)
    locals()['p_down' + str(i)] = p_downi
    #d_upi = 'd_up' + str(i)
    d_upi = (dos[i-1][:,9]+dos[i-1][:,11]+dos[i-1][:,13]+dos[i-1][:,15]+dos[i-1][:,17]).reshape(epoints,1)
    locals()['d_up' + str(i)] = d_upi
    #d_downi = 'd_down' + str(i)
    d_downi = (dos[i-1][:,10]+dos[i-1][:,12]+dos[i-1][:,14]+dos[i-1][:,16]+dos[i-1][:,18]).reshape(epoints,1)
    locals()['d_down' + str(i)] = d_downi
    #groupi = 'group' + str(i)
    groupi = np.hstack((Eni,s_upi,s_downi,p_upi,p_downi,d_upi,d_downi))
    locals()['group' + str(i)] = groupi
    np.savetxt(file,groupi)   # projected to each atomic orbital
    #fi = 'f' + str(i)
    fi = np.loadtxt(file)
    #dfi = 'df' + str(i)
    dfi = pd.DataFrame(fi, columns = group_names)
    dfi.set_index('energy(eV)',inplace=True)
    dfi.to_csv(file,sep='\t',float_format='%12.10f') # projected to each atomic orbital

    ######   sum of dos   ######
    #dos_groupi = 'dos_group' + str(i)
    dos_groupi = np.hstack((s_upi,s_downi,p_upi,p_downi,d_upi,d_downi))
    locals()['dos_group' + str(i)] = dos_groupi
    filedos = path + '/dos' + str(i)
    np.savetxt(filedos,dos_groupi)

dos_sum = np.zeros([epoints,6])
#print dos_sum.shape
for j in range(int(atom_ini),int(atom_final)+1):
    dos_groupj = path + '/dos' + str(j)
    locals()['dos' + str(j)] = dos_groupj
    #foj = 'fo' + str(j)
    foj= np.loadtxt(dos_groupj)
#    print foj.shape
    dos_sum = np.add(dos_sum,foj)
dos_sum = np.hstack((E_corrected, dos_sum))
ddos = pd.DataFrame(dos_sum, columns = group_names)
ddos.set_index('energy(eV)',inplace=True)
ddos.to_csv('sum_dos.dat',sep=' ',float_format='%12.10f') # pdos sum of selected atoms


###### calculaton of p-band center ######
#file = sys.argv[1] # pass paramter
file = 'sum_dos.dat'  # the selected atom
list_energy = [l.split()[0] for l in open(file,'rb')]
list_energy.pop(0)
energy = np.array([float(i) for i in list_energy])  # extract energy

emin, emax = energy[0], 0.05   # integral energy range
erange = (energy[0],energy[-1])
emask = (energy >= emin) & (energy <= emax) # bool to make a mapping between energy and dos

list_pup = [l.split()[3] for l in open(file,'rb')]
list_pup.pop(0)
p_up = np.array([float(i) for i in list_pup])   # extract p_up
list_pdown = [l.split()[4] for l in open(file,'rb')]
list_pdown.pop(0)
p_down = np.array([float(i) for i in list_pdown])   # extract p_down

x = energy[emask]
y1 = p_up[emask]
y2 = p_down[emask]

pbc_up   = simps(y1*x, x) / simps(y1, x)
pbc_down = simps(y2*x, x) / simps(y2, x)
pbc = []
pbc.append(pbc_up)
pbc.append(pbc_down)

###### calculaton of d-band center ######

## same energy set ##

list_dup = [l.split()[5] for l in open(file,'rb')]
list_dup.pop(0)
d_up = np.array([float(i) for i in list_dup])   # extract d_up
list_ddown = [l.split()[6] for l in open(file,'rb')]
list_ddown.pop(0)
d_down = np.array([float(i) for i in list_ddown])   # extract d_down

y1_d = d_up[emask]
y2_d = d_down[emask]

if np.all(y1_d == 0) and np.all(y2_d == 0):
	dbc = 'No d-band electrons !!!'
else:
	dbc_up   = simps(y1_d*x, x) / simps(y1_d, x)
	dbc_down = simps(y2_d*x, x) / simps(y2_d, x)
	dbc = []
	dbc.append(dbc_up)
	dbc.append(dbc_down)

###### calculation of p_z-band center #######

## same energy set ##

dos_splitted = np.loadtxt('atoms_dos.dat')
dos_s = np.zeros([epoints,column_num])
for j in range(int(atom_ini),int(atom_final)+1):
    dos_s_add = np.vsplit(atoms_dos,atoms_number)[j-1].astype(float)
    dos_s_add = np.add(dos_s_add, E_f0) # corrected energy by efermi
    dos_s_add = np.multiply(dos_s_add,ones)
    dos_s = np.add(dos_s,dos_s_add)    # vertical stack

dos_s[:,0] = energy
df_s = pd.DataFrame(dos_s, columns=all_names)
df_s.set_index('energy(eV)', inplace=True)
df_s.to_csv('sum_dos_split.dat', sep=' ', float_format='%-12.10f')  # projected to each atomic orbital orientation

list_pzup = [l.split()[5] for l in open('sum_dos_split.dat','rb')]
list_pzup.pop(0)
pz_up = np.array([float(i) for i in list_pzup])   # extract p_up
list_pzdown = [l.split()[6] for l in open('sum_dos_split.dat','rb')]
list_pzdown.pop(0)
pz_down = np.array([float(i) for i in list_pzdown])   # extract p_down

y1_pz = pz_up[emask]
y2_pz = pz_down[emask]

pzbc_up   = simps(y1_pz*x, x) / simps(y1_pz, x)
pzbc_down = simps(y2_pz*x, x) / simps(y2_pz, x)
pzbc = []
pzbc.append(pzbc_up)
pzbc.append(pzbc_down)

psa = range(int(atom_ini), int(atom_final)+1)

print('the selected atoms:', list(range(int(atom_ini), int(atom_final)+1)))
print('dbc_up(eV), dbc_down(eV)')
print(dbc)
print('p-BC_up(eV), p-BC_down(eV)')
print(pbc)
print('p_z-BC_up(eV), p_z-BC_down(eV)')
print(pzbc)

if os.name == 'posix':
    os.system('rm atoms_dos* DOS_data/dos*')
elif os.name == 'nt':   # Windows系统
    for f in glob.glob('atoms_dos*'):
        try:
            os.remove(f)
        except Exception as e:
            print(f"删除失败: {f} - {e}")

        # 严格删除 DOS_data 目录下以 "dos" 开头（小写）的文件
    for f in glob.glob('DOS_data/dos*'):
        if os.path.basename(f).startswith('dos'):  # 精确检查文件名前缀
            try:
                os.remove(f)
            except Exception as e:
                print(f"删除失败: {f} - {e}")

##########   timing   #############
stop=time.time()
print("running time:\t" + str(stop-start) + " seconds")
terminal_time = datetime.datetime.now()
print("Terminal time:\t" + terminal_time.strftime('%Y.%m.%d-%H:%M:%S'))   #strftime可以自定义时间的输出格式
print('splitted dos and d/p/pz-band center have been obtained')
print('\n********** splitted dos and d&p-band center from Yaqiong Su & Yanyang Qin XJTU & Yifan Zhang XJTU **********\n')

