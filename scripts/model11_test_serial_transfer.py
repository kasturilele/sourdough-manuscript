#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 09:21:13 2024

@author: kasturilele
file just to create toy model for testing how serial transfer model affects coexistence
everything is the same until the actual model, there we use a two-species model and make graphs for both original vs serial transfer to see the difference

"""

import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import csv
import random

#read the data from text files
single_data = pd.read_csv("est_single_all_6-3.txt")
single_data = single_data.to_dict('records')

paired_data = pd.read_csv("est_all_pairs_6-3.txt")
paired_data = paired_data.to_dict('records')

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
              
t_int = np.linspace(0,48,100) #defines a vector of timepoints to solve over
t_int1 = np.linspace(0,288,600) #defines a vector of timepoints to solve over

strain_names = ['17B2','0092a','232','550','460','253','163','228','177']

#all_strains = ['550','163'] #use these to get results for the other pair of species that we explored
all_strains = ['232','460'] 

#strain_color = ['#aa4499','#88ccee']
strain_color = ['#ddcc77','#cc6677']

#making empty lists to save the model results and replicates picked for each run of the model for later
y_endpoint1 = []
y_endpoint2 = []

#initialize single parameter dictionary
single_dict = {x: [] for x in strain_names}

#populate single parameter dictionary
for i in range(0,len(single_data)):
     tempStrain = single_data[i]
     tempKey = tempStrain['Strains']
     tempPars = [tempStrain['N0'],tempStrain['r'],tempStrain['a11']]
     single_dict[tempKey].append(tempPars)

#initialize paired parameter dictionary
paired_dict = {x: {} for x in strain_names}

for i in range(0,len(strain_names)):
    outerVar = strain_names[i]
    delStrain = [x for x in strain_names if x != outerVar]
    paired_dict[outerVar] = {x: [] for x in delStrain}

#populate paired parameter dictionary
for i in range(0,len(paired_data)):
    tempStrain = paired_data[i]
    tempKey1 = tempStrain['Strain_1']
    tempKey2 = tempStrain['Strain_2']
    paired_dict[tempKey1][tempKey2].append(tempStrain['a12'])
    paired_dict[tempKey2][tempKey1].append(tempStrain['a21'])
    
#define a function to make the input list per replicate for single parameters
def appender(strains,array,iList,j,n):
    for a in range(0,n):
        S = strains[a]
        iTemp = iList[a]
        array.append(single_dict[S][iTemp][j])
    return array

# define a similar function for paired parameters
def appenderPairsNew(strains,iList,n):
    pairMatrix = [[0] * n for _ in range(n)]
    idx = 0
    for j in range(0,n):
        for k in range(j+1,n):
            S1 = strains[j]
            S2 = strains[k]
            i = iList[idx]
            pairMatrix[j][k] = paired_dict[S1][S2][i]
            pairMatrix[k][j] = paired_dict[S2][S1][i]
            idx += 1
    return pairMatrix

#function that counts how many replicates are there in a given pair of strains
def varCounter(strains,x,y):
    vS1 = strains[x]
    vS2 = strains[y]
    varLen = len(paired_dict[vS1][vS2])
    return varLen

# 9-species equation solver - with more randomized sampling from parameter distributions
seed_start = 1998

for olN in range(100):
    
    #seed_start = reps[olN]
    #the fixed parameters are set here for each run
    num_species = 2

    #N_init = [19,121] 
    N_init = [490,149] 
    
    random.seed(seed_start)
    testList1 = []
    testList2 = []
    
    for lN in range (num_species):
        testList1.append(random.randint(0,19))
    
    for xlN in range(0,num_species):
        for ylN in range(xlN+1,num_species):
            #uB = varCounter(all_strains, xlN, ylN) - 1
            testList2.append(random.randint(0,19))

    #set the parameters that vary - these are now randomized for both single species as well as pairwise parameters
    g_rate = appender(all_strains, [], testList1, 1, num_species)
    alpha_single = appender(all_strains, [], testList1, 2, num_species)
    
    alpha_pair = appenderPairsNew(all_strains, testList2, num_species)
    
    #model without bottlenecks
    sol = solve_ivp(testFun5, [0, 288], N_init, args = (g_rate, alpha_single, alpha_pair, num_species),dense_output=True)

    y_inter = sol.sol(t_int1) #stores values for all timepoints being solved over
    
    for m in range(600):
        repname = "M"+str(olN+1)
        y_appendee = [repname]
        y_appendee.append(t_int1[m])
        y_appendee.extend(y_inter[:,m])
        y_endpoint1.append(y_appendee)

    for pltnum in range(0,num_species):
        plt.plot(t_int1, y_inter[pltnum], color = strain_color[pltnum], linewidth = 1, label = all_strains[pltnum])
    
    plt.title(olN)
    #plt.legend(loc = "upper left")
    plt.show()
    
    y_all = [[],[]]
    N_temp = N_init #inital species abundances are set here - to the original values
    
    #model with bottlenecks
    for i in range(6): 
        #print(N_temp)
        sol1 = solve_ivp(testFun5, [0, 48], N_temp, args = (g_rate, alpha_single, alpha_pair, num_species),dense_output=True)
        y_inter1 = sol1.sol(t_int)
        for j in range(num_species):
            y_all[j].extend(y_inter1[j])
        #here's the new code - dilutes the N at the end for the start of the new growth cycle
        N_end = y_inter1[:,99]
        N_temp = np.round(N_end/20) #to remove overflow error, rounding the very small numbers

    y_inter1 = sol1.sol(t_int) #stores values for all timepoints being solved over

    t_ext = np.linspace(0,288,600)
    
    for m in range(600):
        repname = "M"+str(olN+1)
        y_appendee = [repname]
        y_appendee.append(t_ext[m])
        y_appendee.append(y_all[0][m])
        y_appendee.append(y_all[1][m])
        y_endpoint2.append(y_appendee)

    for pltnum in range(0,num_species):
        plt.plot(t_ext, y_all[pltnum], color = strain_color[pltnum], linewidth = 1, label = all_strains[pltnum])
    
    plt.title(olN)
    #plt.legend(loc = "upper left")
    plt.show()
    
    seed_start = seed_start + 1 
    print(seed_start)

#uncomment these to save growth curves to file
# with open('egs_no_transfer.csv', 'w') as f:
      
#     # using csv.writer method from CSV package
#     write = csv.writer(f, dialect = 'excel')
#     rowNames = ['ModelID','Time','C.paralimentarius','W.anomalus']
      
#     write.writerow(rowNames)
#     write.writerows(y_endpoint1)

# with open('egs_serial_transfer.csv', 'w') as f:
      
#     # using csv.writer method from CSV package
#     write = csv.writer(f, dialect = 'excel')
#     rowNames = ['ModelID','Time','C.paralimentarius','W.anomalus']
      
#     write.writerow(rowNames)
#     write.writerows(y_endpoint2)
