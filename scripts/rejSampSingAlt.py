import sys
import numpy as np
import math
from collections import defaultdict
import csv

species =  sys.argv[1]

def computeDist(obs,sim,keys):
    tots = []
    num = 0
    for tup in keys:
        tot = 10000
        i = 0
        repStart = 0
        for rep in obs:
            j = 0
            temptot = 0
            for item in obs[rep]:
                if not math.isnan(item):
                    if obs[rep][j] > 0 and sim[tup][i] > 0:
                        temptot += (np.log(obs[rep][j])-np.log(sim[tup][i]))**2
                    else:
                        temptot += 10
                j += 1
                i += 1
            if temptot < tot:
                tot = temptot
                repStart = i-9
        tot = tot**0.5
        tots.append([tot,num,repStart])
        num += 1
       
    tots = sorted(tots, key=lambda x: x[0])
    return tots[:50]  

traj = defaultdict(dict) 

# grab N0 vals for the species of interest
# Strain,Time,Replicate,TOTAL_CFUS_well,scaled

fh = open(sys.argv[2],'r')
for line in fh:
    break
for line in fh:
    data = line.strip().split(',')

    spec = data[0]
    if spec != species:
        continue 
   
    rep = int(data[2])
    if spec not in traj:
        traj[spec][rep] = []    
    if rep not in traj[spec]:
        traj[spec][rep] = []
    traj[spec][rep].append(float(data[3]))
fh.close()

nreps = len(traj[species])    

# read in file of simulated summary stats
fh = open(sys.argv[3], 'r')

sumStats = {}
allData = []
for line in fh:
    break
for line in fh:
    data = line.strip().split(',')
    sumStats[tuple(data[:2])] = [float(x) for x in data[2:]]
    allData.append([float(x) for x in data[2:]])

# z-transform all sumStats
#allDataT = np.transpose(allData)
#means = []
#stds = []
#for row in allDataT:
#    means.append(np.mean(row))
#    stds.append(np.std(row)) 

#for tup in sumStats:
#    i = 0
#    for num in sumStats[tup]:
#        sumStats[tup][i] -= means[i]
#        if stds[i] > 0:
#            sumStats[tup][i] /= stds[i]
#        i += 1

#reps = list(traj[species].keys())
#reps = sorted(reps)
#i = 0
#for rep in reps:
#    j = 0
#    for point in traj[species][rep]:
#        traj[species][rep][j] -= means[i]
#        if stds[i] > 0:
#            traj[species][rep][j] /= stds[i]
#        i += 1
#        j += 1

# compute distance and sort
myKeys = list(sumStats.keys())
mins = computeDist(traj[species],sumStats,myKeys)
times = [0,6,18,24,30,42,48,54,68]
y_endpoint = []
print("Strains,Replicate,N0,r,a11")

rep = 0
for thing in mins:
    print (species+","+str(rep)+","+str(int(sumStats[myKeys[thing[1]]][0])),end=',')
    repStart = thing[2]
    for t in range(0,9):
        y_endpoint.append([species,rep,times[t],sumStats[myKeys[thing[1]]][repStart+t]])
    nu = 0
    for item in myKeys[thing[1]]:
        if nu < len(myKeys[thing[1]]) -1 : 
            print(item,end=',')
        else:
            print(item,end='')
        nu += 1
    print()
    rep += 1

fN = "../singTrajAlt/Traj_"+species+".csv"
with open (fN, "w") as f:
    write = csv.writer(f, dialect = 'excel')
    rowNames = ["Species","Replicate","Time","Population"]
    write.writerow(rowNames) 
    write.writerows(y_endpoint)
