import sys
import csv
import numpy as np
from collections import defaultdict

species1 = sys.argv[3]
species2 = sys.argv[4]

def computeDist(obs,sim,keys,spec1,spec2):
    tots = []
    num = 0
    for tup in keys:
        tot = 0
        for rep in obs[spec1]:
            pos = 0
            for pair in sim[tup]:
                if obs[spec1][rep][pos] > 10**-6 and pair[0] > 10**-6:
                    tot += (np.log(pair[0])-np.log(obs[spec1][rep][pos]))**2
                else:
                    tot += 10
                pos += 1
        for rep in obs[spec2]:
            pos = 0
            for pair in sim[tup]:
                if obs[spec2][rep][pos] > 10**-6 and pair[1] > 10**-6:
                    tot += (np.log(pair[1])-np.log(obs[spec2][rep][pos]))**2
                else:
                    tot += 10
                pos += 1

        tot = tot**0.5
        tots.append([tot,num])
        num += 1
       
    tots = sorted(tots, key=lambda x: x[0])
    return tots[:50]  

traj = defaultdict(dict) 

#read file of obs sum stats
#Time,Strain_1,Replicate,CFU_Strain_1,Strain_2,replicate,CFU_Strain_2
fh = open(sys.argv[1],'r')

for line in fh:
    break
for line in fh:
    data = line.strip().split(',')

    spec1 = data[1]
    spec2 = data[4]
    if not ((species1 == spec1 and species2 == spec2) or (species1 == spec2 and species2 == spec1)):
        continue

    rep = int(data[2])
     
    if spec1 not in traj:
        traj[spec1][rep] = []
    elif rep not in traj[spec1]:
        traj[spec1][rep] = []
    traj[spec1][rep].append(float(data[3]))
    
    if spec2 not in traj:
        traj[spec2][rep] = []
    elif rep not in traj[spec2]:
        traj[spec2][rep] = []
    traj[spec2][rep].append(float(data[6]))

fh.close()

     
# read in file of simulated summary stats for this species pair
fh = open(sys.argv[2], 'r')

sumStats = {}
allData = []
for line in fh:
    break
for line in fh:
    data = line.strip().split(',')
    sumStats[tuple(data[:8])] = []
    allData.append([float(x) for x in data[8:]])
    for i in range(8,8+len(data[8:]),2):
        sumStats[tuple(data[:8])].append([max(float(data[i]),1e-5),max(float(data[i+1]),1e-5)])

"""
allDataT = np.transpose(allData)
means = []
stds = []
for row in allDataT:
    means.append(np.mean(row))
    stds.append(np.std(row)) 

for tup in sumStats:
    i = 0
    for pair in sumStats[tup]:
        pair[0] = pair[0] - means[i]
        pair[0]/=stds[i]
        pair[1] = pair[1] - means[i+1]
        pair[1]/=stds[i+1]
       
        i += 2

for tup in traj:
    for rep in traj[tup]:
        i = 0
        for pair in traj[tup][rep]:
            pair[0] = pair[0] - means[i]
            pair[0]/=stds[i]
            pair[1] = pair[1] - means[i+1]
            pair[1]/=stds[i+1]
            i += 2
"""
tVals = [0,6,18,24,30,42,48,54,68]
y_endpoint = []
# compute distance and sort
myKeys = list(sumStats.keys())

mins = computeDist(traj,sumStats,myKeys,species1,species2)

print("Strain_1,Strain_2,Replicate,N10,N20,r1,r2,a11,a22,a12,a21")

k = 0    
for thing in mins:
    print(species1+","+species2+","+str(k),end=",")
    for t in range(0,9):
        y_endpoint.append([species1,species2,k,tVals[t],sumStats[myKeys[thing[1]]][t][0], sumStats[myKeys[thing[1]]][t][1]])
    j = 0
    for item in myKeys[thing[1]]:
        if j < len( myKeys[thing[1]]) -1:
            print(item, end=",")
        else:
            print(item)
        j += 1
    k += 1


fN = "../pairTraj/Traj_"+species1+"_"+species2+".csv"
with open (fN, "w") as f:
    write = csv.writer(f, dialect = 'excel')
    rowNames = ["Strain_1","Strain_2","Replicate","Time","CFU_Strain_1","CFU_Strain_2"]
    
    write.writerow(rowNames)
    write.writerows(y_endpoint)

