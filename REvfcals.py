#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 2020
for vabrational frequency calculating the correction terms of G
@author: Yanyang Qin, XJTU

Rewritten on Feb 25 2025, with Python3.7, By Yifan Zhang
"""
# Ref:https://www.bigbrosci.com/2018/11/07/ex69/
# email: 750881345@qq.com

from scipy import constants as con
import numpy as np
import linecache as lce

###Create the result file with raw frequency data###
outfile = open('tsresult', mode= 'w')
keyword = '2PiTHz'
with open('OUTCAR', mode= 'r') as getdatafile:
	for line in getdatafile:
		try:
			line_str = line.strip( ).split()
			if '=' in line_str[1]:
				if keyword in line_str[5]:
					outfile.write(line)
			else:
				if keyword in line_str[6]:
					outfile.write(line)
		except:
			continue

outfile.close()
print('OUTCAR read!')

###Calculater###
with open('tsresult', mode= 'r') as count:
	num = len(count.readlines())

with open('tsresult', mode= 'a') as rcal:
	rcal.write('\n\n######\n')
	rcal.write('%-15s\t%-15s\t%-15s\t%-15s\t%-15s\t\n' %('DOFnumber','S_J/K/mol','TS','Cpi','EmeV'))

class Sum:
	def __init__(self):
		self.S = float(0)
		self.TS = float(0)
		self.meV = float(0)
		self.Cp = float(0)

class OUTCAR:
	h = con.h  # get Plank Constant
	k = con.k  # get Boltzman Constant
	R = con.R  # get Gas Constant
	c = con.c  # get Lightspeed Constant
	T = 300  # Set Temperature as 300K
	beta = 1 / (k * T)

	def __init__(self, num):
		self.WaveNumber = np.zeros(num, dtype=float)
		self.EmeV = np.zeros(num, dtype=float)
		self.nu = np.zeros(num, dtype=float)

	def assign(self, wavenum, emev, row):
		self.WaveNumber[row-1] = float(wavenum)
		self.EmeV[row-1] = float(emev)
		if self.WaveNumber.all() != 0:
			self.nu = self.WaveNumber * 100  # convert unit from cm-1 to m-1

	def get_pf(self):
		x = self.beta * self.h * self.c * self.nu  # Numerator of the first term of the equation
		pf1 = x / (np.exp(x) - 1)  # First term of the equation
		pf2 = np.log(1 - np.exp(-x))  # Second term *
		pf = pf1 - pf2
		self.entropy = self.R * pf # unit is J*K-1*mol-1

	def calculate(self):
		self.get_pf()
		self.TS = self.entropy * self.T / 1000 / 96.485  # unit is eV
		epsilon_i = self.h * self.c * self.nu
		expot_i = epsilon_i * self.beta
		self.PCapacity_i = epsilon_i / con.e / (np.exp(expot_i) - 1)

OUTCAR = OUTCAR(num)
Sum = Sum()

for i in range(1, num+1):
	line_str = lce.getline('tsresult', i)
	str_spl = line_str.split()
	if '=' in str_spl[1]:
		OUTCAR.assign(str_spl[6], str_spl[8], i)
	else:
		OUTCAR.assign(str_spl[7], str_spl[9], i)

OUTCAR.calculate()

for i in range(0, num):
	line_str = lce.getline('tsresult', i+1)
	str_spl = line_str.split()
	if '=' not in str_spl[1]:
		Sum.S += OUTCAR.entropy[i]
		Sum.TS += OUTCAR.TS[i]
		Sum.meV += OUTCAR.EmeV[i]
		Sum.Cp += OUTCAR.PCapacity_i[i]
	with open('tsresult', mode='a') as rcal:  # write into the result file
		rcal.write('%-15s\t%-15.4f\t%-15.4f\t%-15.4f\t%-15.4f\t\n' % (str_spl[0], OUTCAR.entropy[i], OUTCAR.TS[i], OUTCAR.PCapacity_i[i], OUTCAR.EmeV[i]))
	print(str_spl)
	print('%-15s\t%-15.4f\t%-15.4f\t%-15.4f\t%-15.4f\t' % (str_spl[0], OUTCAR.entropy[i], OUTCAR.TS[i], OUTCAR.PCapacity_i[i], OUTCAR.EmeV[i]))

ZPEeV = Sum.meV / 2000
Dif = ZPEeV + Sum.Cp - Sum.TS

with open('tsresult',mode='a') as rcal:                    #write the sumup value
	rcal.write('######\n%-15s\t%-15s\t%-15s\t%-15s\t%-15s\t\n%-15s\t%-15.4f\t%-15.4f\t%-15.4f\t%-15.4f\t\n%-15s\t%-15.4f\n' %('/','Sum_S(J/K/mol)','Sum_TS(eV)','CvT(eV)','ZPE(eV)','Sum/Final',Sum.S,Sum.TS,Sum.Cp,ZPEeV,'TotDiff is',Dif))
