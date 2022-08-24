#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 15:33:29 2020

@author: shkim
"""

# import open3d as o3d
import numpy as np
import os
import os.path as osp
import sys
from scipy.spatial.transform import Rotation as R
            
def read_off(in_file):
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

def save_off(verts, faces, filepath):
    with open(filepath, "w") as off_file:
        off_file.write("OFF\n")
        off_file.write("{} {} {}\n".format(len(verts), len(faces), 0))
        for i in verts:
            # off_file.write("{} {} {} {} {} {} {}\n".format(i[0], i[1], i[2], i[3], i[4], i[5], i[6]))
            off_file.write("{} {} {}\n".format(i[0], i[1], i[2]))
        for n, j in enumerate(faces):
            off_file.write("3 {} {} {}\n".format(j[0], j[1], j[2]))

def save_txt(content, filepath):
    with open(filepath, 'w') as f:
        for c in content:
            f.write(str(c[0]))
            f.write(' ')
            f.write(str(c[1]))
            f.write(' ')
            f.write(str(c[2]))
            f.write('\n')        
            
def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix


patients_dir = os.getcwd() + '/../io/landmark'

if __name__ == "__main__":
    # patients_dir = sys.argv[1]
    # register target to the template
    patient = 'target.off'
    pathLoad = osp.join(patients_dir, patient)
    vertex_target, faces_target = read_off(pathLoad)
    fish_target = np.zeros((24,3)).astype(np.float64)
    with open(osp.join(patients_dir, 'target.txt'), 'r') as f:
        for i in range(24):
            l = f.readline()
            fish_target[i] = np.array(l.split(' '))
    
    rotation_vec = (fish_target[14] - fish_target[5]) + (fish_target[11] - fish_target[0])
    rotation = rotation_matrix_from_vectors(rotation_vec, np.array([0, 0, 1], dtype=np.float64))
    x = R.from_matrix(rotation).as_euler('xyz', degrees=True)[0]
    x -= 5
    rot = R.from_euler('x', [x], degrees=True)
    fish_target = rot.apply(fish_target)

    canting_vec = (fish_target[11] - fish_target[14]) + (fish_target[0] - fish_target[5]) \
                 + (fish_target[1] - fish_target[4])
    canting = rotation_matrix_from_vectors(canting_vec, np.array([1, 0, 0], dtype=np.float64))
    y = R.from_matrix(canting).as_euler('xyz', degrees=True)[1]
    rot = R.from_euler('y', [y], degrees=True)
    fish_target = rot.apply(fish_target)
    
    yawing_vec = (fish_target[11] - fish_target[14]) + (fish_target[0] - fish_target[5]) \
                 + (fish_target[1] - fish_target[4])
    yawing = rotation_matrix_from_vectors(yawing_vec, np.array([1, 0, 0], dtype=np.float64))
    z = R.from_matrix(yawing).as_euler('xyz', degrees=True)[2]
    rot = R.from_euler('z', [z], degrees=True)
    fish_target = rot.apply(fish_target)
    
    rot = R.from_euler('xyz', [x,y,z], degrees=True)
    vertex_rotated = rot.apply(vertex_target)

    save_off(vertex_rotated, faces_target, osp.join(patients_dir, 'target.off'))
    save_txt(fish_target, osp.join(patients_dir,'target.txt'))
    
    with open(osp.join(patients_dir, 'result_trimmed-Align.txt'), 'a') as f:
        f.write(f"Rotation : [{x},{y},{z}]")
