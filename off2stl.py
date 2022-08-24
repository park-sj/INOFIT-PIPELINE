#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:48:02 2020

@author: shkim
"""

import numpy as np
import sys
from stl import mesh

def readOff(in_file):
    with open(in_file, 'r') as f:
        src = f.read().split('\n')[:-1]
    if src[0] == 'OFF' or src[0] == 'COFF':
        src = src[1:]
    else:
        src[0] = src[0][3:]

    num_nodes, num_faces = [int(item) for item in src[0].split()[:2]]
    vert = np.zeros((num_nodes, 3))
    for i in range(num_nodes):
        l = src[i+1].split(' ')
        for j in range(3):
            vert[i, j] = float(l[j])
        # for j in range(4):
        #     color[i, j] = int(l[j+3])
    
    faces = np.zeros((num_faces, 3), dtype=np.int32)
    color = np.zeros((num_faces, 3), dtype=np.int16)
    for i in range(num_faces):
        l = src[1+num_nodes+i].split(' ')
        for j in range(3):
            faces[i, j] = int(l[j+1])
        for j in range(3):
            color[i, j] = int(l[j+4])
    return vert, faces, color

            
def saveStl(verts, faces, filepath):
    stlMesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            stlMesh.vectors[i][j] = verts[f[j],:]
    stlMesh.save(filepath, update_normals=False)

if __name__ == "__main__":
    if (len(sys.argv) != 3):
        raise Exception("Wrong number of argv")
    print(f"Convert {sys.argv[1]} to {sys.argv[2]}")
    vert, faces, _ = readOff(sys.argv[1])
    saveStl(vert, faces, sys.argv[2])