#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 14:32:38 2021

@author: shkim
"""

import subprocess
import os
import sys
import rotate
import csv
import numpy as np
from util import calculate_stat

queueDir = '/home/shkim/Libraries/macro/queue'
queueTempDir = '/home/shkim/Libraries/macro/queue_temp'
resultDir = '/home/shkim/Libraries/macro/result'
resultTempDir = '/home/shkim/Libraries/macro/result_temp'

if __name__ == "__main__":
    innofitDir = sys.argv[1]
    targetDir = sys.argv[2]
    csvDir = sys.argv[3]
    isMesh = int(sys.argv[4])
    print(f"{innofitDir} {targetDir} {csvDir} {isMesh}")
    
    if isMesh == 0:
        raise NotImplementedError
        ''' If queue folder exists and is not empty,
              first, rename the existing queue folder as queue_temp,
              second, move innofitDir into queue folder and run iter.sh (and targetDir, consecutively)
              lastly, restore the renamed queue_temp folder.
            Else if queue folder dosen't exist or exists but empty,
              skip the first process '''
        if os.path.isdir(queueDir):
            assert not os.path.isdir(queueTempDir), "both the queue, queue_temp folders already exist."
            os.rename(queueDir, queueTempDir)
        if os.path.isdir(resultDir):
            assert not os.path.isdir(resultTempDir), "both the result, result_temp folders already exist."
            os.rename(resultDir, resultTempDir)
        
        os.mkdir(resultDir)
        os.rename(innofitDir, queueDir)
        subprocess.run(['source', 'iter.sh'])
        os.rename(queueDir, innofitDir)
        if innofitDir[-1] == '/':
            innofitDir = innofitDir[:-1] + '_mesh'
        else:
            innofitDir = innofitDir + '_mesh'
        os.rename(resultDir, innofitDir)
        
        os.mkdir(resultDir)
        os.rename(targetDir, queueDir)
        subprocess.run(['source', 'iter_mask.sh'])
        os.rename(queueDir, targetDir)
        if targetDir[-1] == '/':
            targetDir = targetDir[:-1] + '_mesh'
        else:
            targetDir = targetDir + '_mesh'
        os.rename(resultDir, targetDir)
        
        if os.path.isdir(queueTempDir):
            os.rename(queueTempDir, queueDir)
        else:
            os.mkdir(queueDir)
        if os.path.isdir(resultTempDir):
            os.rename(resultTempDir, resultDir)
        else:
            os.mkdir(resultDir)
            

    elif isMesh == 1:
        patientList = os.listdir(innofitDir)
        
        # check all the patients have their ground truth
        for p in patientList:
            assert os.path.isdir(os.path.join(targetDir, p)), f"Coudn't find ground truth of {p}"
    
        if not os.path.isdir(csvDir):
            os.mkdir(csvDir)
        
    # evaluate each patient and record it in csv directory
        csvList = []
        for p in patientList:
            
            rotate.run(os.path.join(innofitDir, p, "result.off"),
                       os.path.join(innofitDir, p, "result_rotated.off"),
                       os.path.join(innofitDir, p, "result_trimmed-Align.txt"))
            ''' rotate.run(os.path.join(targetDir, p, "result_trimmed.off"),
                       os.path.join(targetDir, p, "result_trimmed_rotated.off"),
                       os.path.join(targetDir, p, "rotation.txt")) '''
            
            subprocess.run(["./eval_two_mesh", os.path.join(innofitDir, p, "result_rotated.off"),
                            os.path.join(targetDir, p, p + ".off"), os.path.join(csvDir, p + ".csv")])
            csvList.append(os.path.join(csvDir, p + ".csv"))
        
    # calculate statistics using csv files
        means = []
        vars = []
        maximums = []
        ports = []
        for c in csvList:
            with open(c, 'r') as f:
                rdr = csv.reader(f)
                dist = []
                for l in rdr:
                    dist.append(float(l[1]))
            dist = np.array(dist)
            mean, var, maximum, port = calculate_stat(dist, 1.)
            means.append(mean)
            vars.append(var)
            maximums.append(maximum)
            ports.append(port)
        
        with open(os.path.join(csvDir, "summary.csv"), "w") as f:
            wr = csv.writer(f)
            wr.writerow(["", "Mean", "Variance", "Maximum", "Portion"])
            for i in range(len(csvList)):
                wr.writerow([patientList[i], means[i], vars[i], maximums[i], ports[i]])
                print(f"{patientList[i]} - Mean: {means[i]}, Variation: {vars[i]}, Max: {maximums[i]}, Port: {ports[i]}")
            
            print(f"Avg - Mean: {np.mean(np.array(means))}, Variation: {np.mean(np.array(vars))}, Max: {np.mean(np.array(maximums))}, Port: {np.mean(np.array(ports))}")
            wr.writerow(["Avg", np.mean(np.array(means)), np.mean(np.array(vars)),np.mean(np.array(maximums)), np.mean(np.array(ports))])
