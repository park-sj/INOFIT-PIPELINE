#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 11:02:18 2020

@author: shkim
"""
import os
import os.path as osp
import sys
import numpy as np
from scipy.spatial.transform import Rotation as R
from remesh_scanner import readOff, saveOff

patients_dir = os.getcwd() + '/../io/mesh/raw/'

if __name__ == "__main__":
    # patients_dir = sys.argv[1]
    # register target to the template
    patient = 'target.off'
    pathLoad = osp.join(patients_dir, patient)
    vertex, faces = readOff(pathLoad)
    
    with open(os.path.join(os.getcwd(),'..','mesh_alignment','scanner_offset.txt'), 'r') as f:
        offset = list(map(float, f.read().split(' ')))
    
    rot = R.from_euler('xyz', offset, degrees=True)
    vertex_rotated = rot.apply(vertex)
    saveOff(vertex_rotated, faces, osp.join(patients_dir, 'target.off'))