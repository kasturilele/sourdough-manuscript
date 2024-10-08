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

spec0 = sys.argv[1]
spec1 = sys.argv[2]

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

fh = open(sys.argv[3], 'r')
params = []
for line in fh:
    data = line.strip().split(",")
    if not ((spec0 == data[0] and spec1 == data[1]) or (spec1 == data[0] and spec0 == data[1])): 
        continue
    params.append([float(x) for x in data[3:]])

# 9-species equation solver - sampling from gamma distribution for r and exponential distribution for aii and aij
y_endpoint = []
print("time,spec1,spec2,rep")
for olN in range(len(params)):
    
    num_spec = 2
    #N_init = appender(all_strains, [], i, 0, num_species)
    N_init = [params[olN][0],params[olN][1]]

    #set the parameters that vary - these are now randomized for both single species as well as pairwise parameters
    #g_rate = gamma.rvs(a = 1, scale = 1, size = num_species)

    g_rate = [params[olN][2],params[olN][3]]

    #alpha_single = -expon.rvs(scale = a_mean, size = num_species)
    alpha_single = [params[olN][4],params[olN][5]]

    alpha_pair = [[0,params[olN][6]],[params[olN][7],0]]
    
    sol = solve_ivp(testFun5, [0, 68], N_init, args = (g_rate, alpha_single, alpha_pair, num_spec), dense_output=True)

    y_inter = sol.sol(t_int) #stores values for all timepoints being solved over
    
    out = N_init+g_rate+alpha_single
    for pair in alpha_pair:
        for item in pair:
            if item != 0.:
                out += [item]
    times  = [0,6,18,24,30,42,48,54,68]
  
    for time in times:
        print(time, end=",")
        for thing in y_inter[:,time]:
            print(thing, end=",")
        print(olN)
 
"""
fN = 'test.csv'    
with open(fN, 'w') as f:
      
    # using csv.writer method from CSV package
    write = csv.writer(f, dialect = 'excel')

    outNames = ["N0","N1","g0","g1","a00","a11","a01","a10","sum_stats"]
 
    write.writerow(outNames)
    write.writerows(y_endpoint)
""" 
