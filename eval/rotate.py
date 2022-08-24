#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 15:44:32 2021

@author: shkim
"""

import sys
import numpy as np
from scipy.spatial.transform import Rotation as R


def read_off_with_color(in_file):
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
    
    faces = np.zeros((num_faces, 3), dtype=np.int32)
    color = np.zeros((num_faces, 3), dtype=np.int16)
    for i in range(num_faces):
        l = src[1+num_nodes+i].split(' ')
        for j in range(3):
            faces[i, j] = int(l[j+1])
        for j in range(3):
            color[i, j] = int(l[j+4])
    return vert, faces, color

def save_off_with_color(verts, faces, color, filepath):
    with open(filepath, "w") as off_file:
        off_file.write("OFF\n")
        off_file.write("{} {} {}\n".format(len(verts), len(faces), 0))
        for i in verts:
            # off_file.write("{} {} {} {} {} {} {}\n".format(i[0], i[1], i[2], i[3], i[4], i[5], i[6]))
            off_file.write("{} {} {}\n".format(i[0], i[1], i[2]))
        for n, j in enumerate(faces):
            off_file.write("3 {} {} {} {} {} {}\n".format(j[0], j[1], j[2], color[n, 0], color[n, 1], color[n,2]))
            
def rotate(verts, x, y, z):
    rot = R.from_euler('zyx', [-z,-y,-x], degrees=True)
    verts = rot.apply(verts)
    return verts

def run(inputMesh, outputMesh, rotationTxt):
    verts, faces, colors = read_off_with_color(inputMesh)
    
    f = open(rotationTxt)
    lines = f.read().split('\n')
    for line in lines:
        line = line.split(' ')
        if line[0] == 'Rotation':
            x, y, z = map(float, line[2].strip('[').strip(']').split(','))

    print(x,y,z)
    
    verts = rotate(verts, x, y, z)
    save_off_with_color(verts, faces, colors, outputMesh)
    
if __name__ == "__main__":
    inputMesh = sys.argv[1]
    outputMesh = sys.argv[2]
    rotationTxt = sys.argv[3]
    
    run(inputMesh, outputMesh, rotationTxt)
    
