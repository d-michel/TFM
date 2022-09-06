#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 17:32:43 2022

@author: michel
"""

import numpy as np
import matplotlib.pyplot as p

str_txt = "Ec_150_spectrum_16.txt"
f = open(str_txt, 'r')
lines = f.readlines()
data_aux = []
eigvals = []
for x in lines:
    data_aux.append(x.split("\t"))
for i in range(3, len(data_aux)):
    eigvals.append(float(data_aux[i][0]))

def density_of_states(eigvals, nwindows):
    window_width = (eigvals[len(eigvals)-1] - eigvals[0])/nwindows
    density = np.zeros((2,nwindows))
    
    j = 0
    for i in range(nwindows):
        density[0][i] = eigvals[0] + (i+1)*window_width
        while ((j < len(eigvals)) and (eigvals[j] < density[0][i])):
            density[1][i] += 1
            j += 1
    
    return density

density = density_of_states(eigvals, 20)
p.plot(density[0]/16, density[1])

f = open("density_16.txt", 'w')
for i,j in zip(density[0], density[1]):
    f.write("{}\t{}\n".format(i, j))
f.close()
