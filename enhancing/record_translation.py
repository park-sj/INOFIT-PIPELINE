#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 10:43:45 2021

@author: shkim
"""

import os
import os.path as osp
from reader import loadDicom

patients_dir = os.getcwd() + '/../io/test'
# patients_dir = '/home/shkim/Libraries/pytorch-3dunet/datasets/JW/test_masks'
# patients = ['LeeGangHu']
patients = os.listdir(patients_dir)

for patient in patients:
    pathLoad = osp.join(patients_dir, patient)
    series_reader, _, _ = loadDicom(pathLoad)
    
    origin = series_reader.GetMetaData(0, "0020|0032").split('\\')
    with open(os.getcwd() + '/../io/landmark/result_trimmed-Align.txt', 'w') as f:
        f.write(f"Translation : [{float(origin[0])},{float(origin[1])},{float(origin[2])}]\n")
