#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:00:40 2020

This function transfers vertex color of one mesh to another.

All the topology of two meshes should be identical.

@author: shkim
"""

#import os.path as osp
import os
import numpy as np

def read_off_without_color(in_file):
    with open(in_file, 'r') as f:
        src = f.read().split('\n')[:-1]
    if src[0] == 'OFF':
        src = src[1:]
    else:
        src[0] = src[0][3:]

    num_nodes, num_faces = [int(item) for item in src[0].split()[:2]]
    vert = np.zeros((num_nodes, 3))
    for i in range(num_nodes):
        l = src[i+1].split(' ')
        vert[i, 0] = float(l[0])
        vert[i, 1] = float(l[1])
        vert[i, 2] = float(l[2])
    # vert = parse_txt_array(src[1:1 + num_nodes])
    
    faces = np.zeros((num_faces, 3), dtype=np.int32)
    for i in range(num_faces):
        l = src[1+num_nodes+i].split(' ')
        faces[i, 0] = int(l[1])
        faces[i, 1] = int(l[2])
        faces[i, 2] = int(l[3])
    return vert, faces

def read_off_with_color(in_file):
    with open(in_file, 'r') as f:
        src = f.read().split('\n')[:-1]
    if src[0] == 'OFF':
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

def save_off_with_color(verts, faces, color, filepath):
    with open(filepath, "w") as off_file:
        off_file.write("OFF\n")
        off_file.write("{} {} {}\n".format(len(verts), len(faces), 0))
        for i in verts:
            # off_file.write("{} {} {} {} {} {} {}\n".format(i[0], i[1], i[2], i[3], i[4], i[5], i[6]))
            off_file.write("{} {} {}\n".format(i[0], i[1], i[2]))
        for n, j in enumerate(faces):
            off_file.write("3 {} {} {} {} {} {}\n".format(j[0], j[1], j[2], color[n, 0], color[n, 1], color[n,2]))

# def transfer_color(gvert, color):
#     gvert_colored = np.zeros((gvert.shape[0], 7))
    
#     return gvert
        
if __name__ == "__main__":
    path = "../io/landmark/"
    _, _, color = read_off_with_color('template.off')
    gvert, gfaces = read_off_without_color(os.path.join(path, 'result.off'))
    save_off_with_color(gvert, gfaces, color, os.path.join(path, 'result.off'))



    #path = "../io/landmark/"
    #patients = os.listdir(path)
     
    #pathLoad = osp.join(path, patients[0])
    #print(pathLoad)
    #_, _, color = read_off_with_color('template.off')
    #gvert, gfaces = read_off_without_color(os.path.join(pathLoad, 'target_0_aligned_D_final.off'))
    #save_off_with_color(gvert, gfaces, color, os.path.join(pathLoad, 'target_template.off'))
