#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:48:02 2020

@author: shkim
"""

import numpy as np
import sys

def readOffColor(in_file):
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
    for i in range(num_faces):
        l = src[1+num_nodes+i].split(' ')
        for j in range(3):
            faces[i, j] = int(l[j+1])
    return vert, faces


def saveObj(verts, faces, filepath):
    with open(filepath, "w") as obj_file:
        for i in verts:
            obj_file.write("v {} {} {}\n".format(i[0], i[1], i[2]))
            # obj_file.write("v {} {} {}\n".format(-i[2]+80, i[1]-115, i[0]+660))
        for j in faces:
            obj_file.write("f {} {} {}\n".format(j[0]+1, j[1]+1, j[2]+1))

if __name__ == "__main__":
    if (len(sys.argv) != 3):
        raise Exception("Wrong number of argv")
    print(f"Convert {sys.argv[1]} to {sys.argv[2]}")
    #vert, faces, _ = readOffColor(sys.argv[1])
    vert, faces = readOff(sys.argv[1])
    saveObj(vert, faces, sys.argv[2])
