#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 12:30:40 2024

@author: kasturilele

change from 5.2 - change the parameter sampling process so that they are sampled from distributions 
 instead of being sampled from the limited subset of observed parameters
 approximate baysian computation
 - got rid of all the functions used to store/compute parameters
 
"""

import pandas as pd
import numpy as np
import csv
import random
import sys
from scipy.stats import truncnorm
from scipy.integrate import solve_ivp
from scipy.spatial import distance
from collections import defaultdict
from datetime import datetime

fileNum = sys.argv[1]

spec0 = sys.argv[2]
spec1 = sys.argv[3]

t_int = np.linspace(0,68,69) #defines a vector of timepoints to solve over

#function that converts a list to a tuple - for the output of the next function
def convert(list):
    return (*list, )

#final function - constructs the equations as it goes
def testFun5(t, y, r, a_s, a_p, N):
    outputs = []
    for j in range(0,N):
        pairs = []
        for k in range(0,N):
            if j == k:
                pairs.append(0)
            else:
                pairs.append(a_p[j][k]*y[k])
        outputs_calc = y[j]*(r[j] + a_s[j]*y[j] + sum(pairs))
        outputs.append(outputs_calc)
    outputs_tuple = convert(outputs)
    return(outputs_tuple)

# get each single-species growth rates and alphas

par_spec0 = [] 
fh = open('../singSpecEstimates/est.'+spec0+'.txt','r')
for line in fh:
    break
for line in fh:
    par_spec0.append([float(x) for x in line.strip().split(',')[3:]])

par_spec1 = [] 
fh = open('../singSpecEstimates/est.'+spec1+'.txt','r')
for line in fh:
    break
for line in fh:
    par_spec1.append([float(x) for x in line.strip().split(',')[3:]])

# 9-species equation solver - sampling from gamma distribution for r and exponential distribution for aii and aij
y_endpoint = []
for olN in range(10000):
    
    num_spec = 2
    #N_init = appender(all_strains, [], i, 0, num_species)
    N_init = [np.random.randint(25,5000),np.random.randint(25,5000)]

    #set the parameters that vary - these are now randomized for both single species as well as pairwise parameters
    #g_rate = gamma.rvs(a = 1, scale = 1, size = num_species)
    pos0 = np.random.randint(0,len(par_spec0))
    pos1 = np.random.randint(0,len(par_spec1))

    g_rate = [par_spec0[pos0][0],par_spec1[pos1][0]]

    #alpha_single = -expon.rvs(scale = a_mean, size = num_species)
    alpha_single = [par_spec0[pos0][1],par_spec1[pos1][1]]

    alpha_pair = []
    for spec in range(num_spec):
        pair = []
        for spec2 in range(num_spec):
            if spec == spec2:
                pair.append(0)
            else:
                #flip = np.random.random()
                #if flip < 0.5:
                #    pair.append(1*10**(-5-4*np.random.random()))
                #else:
                pair.append(-1*10**(-5-4*np.random.random()))
        alpha_pair.append(pair)
    
    sol = solve_ivp(testFun5, [0, 68], N_init, args = (g_rate, alpha_single, alpha_pair, num_spec), dense_output=True)

    y_inter = sol.sol(t_int) #stores values for all timepoints being solved over
    
    out = N_init+g_rate+alpha_single
    for pair in alpha_pair:
        for item in pair:
            if item != 0.:
                out += [item]
    times  = [0,6,18,24,30,42,48,54,68]
    for time in times:
        out.extend(y_inter[:,time])
    y_endpoint.append(out)
     
fN = '../rejDataPaired/paramsAndDens.'+spec0+'.'+spec1+'.'+fileNum+'.csv' 

with open(fN, 'w') as f:
      
    # using csv.writer method from CSV package
    write = csv.writer(f, dialect = 'excel')

    outNames = ["N0","N1","g0","g1","a00","a11","a01","a10","sum_stats"]
 
    write.writerow(outNames)
    write.writerows(y_endpoint)
    
