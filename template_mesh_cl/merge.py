#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 13:02:36 2021

@author: shkim
"""

import numpy as np
import sys


def read_off(in_file):
    with open(in_file, 'r') as f:
        src = f.read().split('\n')[:-1]
    if src[0] == 'OFF':
        src = src[1:]
    else:
        src[0] = src[0][3:]

    num_nodes, num_faces = [int(item) for item in src[0].split()[:2]]
    vertices = np.zeros((num_nodes, 3))
    for i in range(num_nodes):
        l = src[i+1].split(' ')
        vertices[i, 0] = float(l[0])
        vertices[i, 1] = float(l[1])
        vertices[i, 2] = float(l[2])
    # vertices = parse_txt_array(src[1:1 + num_nodes])
    
    faces = np.zeros((num_faces, 4), dtype=np.int32)
    for i in range(num_faces):
        l = src[1+num_nodes+i].split(' ')
        faces[i, 0] = int(l[1])
        faces[i, 1] = int(l[2])
        faces[i, 2] = int(l[3])
        faces[i, 3] = int(l[4])
    return vertices, faces

def save_off(verts, faces, filepath):
    with open(filepath, "w") as off_file:
        off_file.write("OFF\n")
        off_file.write("{} {} {}\n".format(len(verts), len(faces), 0))
        for i in verts:
            off_file.write("{} {} {}\n".format(i[0], i[1], i[2]))
            # off_file.write("{} {} {}\n".format(-i[2]+80, i[1]-115, i[0]+660)) # JW
            # off_file.write("{} {} {}\n".format(-i[2]+115, i[1]-125, i[0]+625)) # H
            # off_file.write("{} {} {}\n".format(-i[2]/4+80, i[1]/4-100, i[0]/4+660))
        for j in faces:
            off_file.write("4 {} {} {} {}\n".format(j[0], j[1], j[2], j[3]))
            
def read_indices(filepath):
    with open(filepath, "r") as txt_file:
        indices = txt_file.read().split('\n')[:-1]
    indices = list(map(int, indices))
    return indices

# def merge():
    
    

if __name__ == "__main__":
    # if len(sys.argv) != 4:
    #     raise ValueError
    
    right_path = sys.argv[1]
    left_path = sys.argv[2]
    indices_path = "./right_indices.txt"
    result_path = sys.argv[3]
    
    left_verts, left_faces = read_off(left_path)
    right_verts, right_faces = read_off(right_path)
    right_indices = read_indices(indices_path)
    
    left_verts[right_indices] = right_verts[right_indices]
    
    save_off(left_verts, left_faces, result_path)
