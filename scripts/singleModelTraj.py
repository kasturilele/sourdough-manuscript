#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 12:30:40 2024

@author: kasturilele
file just to print a small subset of randomly generated trajectories - as a check of how we are exploring the parameter space
 
"""

import numpy as np
import csv
import sys
from scipy.integrate import solve_ivp
from collections import defaultdict

fileNum = sys.argv[1]
species = sys.argv[2]
t_int = np.linspace(0,68,69) #defines a vector of timepoints to solve over

#function that converts a list to a tuple - for the output of the next function
def convert(list):
    return tuple(list)

#final function - constructs the equations as it goes
def testFun5(t, y, r, a_s, N):
    outputs = []
    for j in range(0,N):
        outputs_calc = y[j]*(r[j] + a_s[j]*y[j])
        outputs.append(outputs_calc)
    outputs_tuple = convert(outputs)
    return(outputs_tuple)

# grab N0 vals for the species of interest
# Strain,Time,Replicate,TOTAL_CFUS_well
fh = open(sys.argv[3],'r')
N_init = []
for line in fh:
    data = line.strip().split(',')
    if data[0] == species and data[1] == '0':
        N_init.append(int(data[3]))
fh.close()

# account for 20x dilution and randomness
for i in range(len(N_init)):
    N_init[i] = 20* np.random.poisson(lam=(N_init[i]/20.))

# 1-species equation solver - sampling from gamma distribution for r and exponential distribution for aii and aij
y_endpoint = []
for olN in range(100):
    
    num_spec = 1

    #set the parameters that vary - these are now randomized for both single species as well as pairwise parameters
    #g_rate = gamma.rvs(a = 1, scale = 1, size = num_species)
    g_rate = []
    for spec in range(num_spec):
        g_rate.append(3*np.random.random())

    #alpha_single = -expon.rvs(scale = a_mean, size = num_species)
    alpha_single = []
    for spec in range(num_spec):
        alpha_single.append(-1*10**(-3-7*np.random.random()))
        
    sols = []
    for N in N_init:
        sol = solve_ivp(testFun5, [0, 68], [N], args = (g_rate, alpha_single, num_spec), dense_output=True)
        sols.append(sol)

    y_inters = []
    for sol in sols:
        y_inter = sol.sol(t_int) #stores values for all timepoints being solved over
        y_inters.append(y_inter)
    

    times  = [0,6,18,24,30,42,48,54,68]
    for y_inter in y_inters:
        for time in times:
            out = [species,olN,time]
            out.extend(y_inter[:,time])
            y_endpoint.append(out)
     
fN = '../singAllTraj.csv' 

with open(fN, 'a') as f:
      
    # using csv.writer method from CSV package
    write = csv.writer(f, dialect = 'excel')

    #outNames = ["g0","a00","sum_stats"]
 
    #write.writerow(outNames)
    write.writerows(y_endpoint)
    
