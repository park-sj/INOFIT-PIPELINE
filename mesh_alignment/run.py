#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 10:41:15 2020

@author: shkim
"""
import numpy as np
import os
import os.path as osp

def yield_file(in_file):
    f = open(in_file)
    buf = f.read()
    f.close()
    for b in buf.split('\n'):
        if b.startswith('v '):
            yield ['v', [float(x) for x in b.split()[1:]]]
        elif b.startswith('f '):
            triangle = b.split(' ')[1:]
            # -1 as .obj is base 1 but the Data class expects base 0 indices
            yield ['f', [[int(t.split("/")[0]) - 1 for t in triangle]]]
            # yield ['f', [[int(i) - 1 for i in t.split("/")] for t in triangle]]
        else:
            yield ['', ""]

def readObj(in_file):
    verts = []
    faces = []

    for k, v in yield_file(in_file):
        if k == 'v':
            verts.append(v)
        elif k == 'f':
            for i in v:
                faces.append(i)

    if not len(faces) or not len(verts):
        return None, None

    return verts, faces

def saveOff(verts, faces, filepath):
    with open(filepath, "w") as off_file:
        off_file.write("OFF\n")
        off_file.write("{} {} {}\n".format(len(verts), len(faces), 0))
        for i in verts:
            # off_file.write("{} {} {}\n".format(i[0], i[1], i[2]))
            # off_file.write("{} {} {}\n".format(-i[2]+80, i[1]-115, i[0]+660)) # JW
            off_file.write("{} {} {}\n".format( (i[0]-80), (i[2]-660), -(i[1]+115))) # working # H
            # off_file.write("{} {} {}\n".format( (i[0]), (i[1]-2), (i[2])))
            # off_file.write("{} {} {}\n".format(-i[2]/4+80, i[1]/4-100, i[0]/4+660))
        for j in faces:
            off_file.write("3 {} {} {}\n".format(j[0], j[1], j[2]))

def readOff(in_file):
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
    
    faces = np.zeros((num_faces, 3), dtype=np.int32)
    for i in range(num_faces):
        l = src[1+num_nodes+i].split(' ')
        faces[i, 0] = int(l[1])
        faces[i, 1] = int(l[2])
        faces[i, 2] = int(l[3])
    return vertices, faces

def readTxt(in_file):
    t = np.zeros((24,3)).astype(np.float64)
    with open(in_file, 'r') as f:
        for i in range(24):
            l = f.readline()
            t[i] = np.array(l.split(' '))
    return t

def saveTxt(content, filepath):
    with open(filepath, 'w') as f:
        for c in content:
            f.write(str(c[0]-80))
            f.write(' ')
            f.write(str(c[2]-660))
            f.write(' ')
            f.write(str(-(c[1]+115)))
            f.write('\n')        

patients_dir = os.getcwd() + '/../io/landmark'
# patients_dir = '/home/shkim/Libraries/pytorch-3dunet/datasets/JW/test_masks'
# patients = ['LeeGangHu']
patients = os.listdir(patients_dir)

for patient in patients:
    pathLoad = osp.join(patients_dir, patient)
    if pathLoad[-3:] == 'off':
    	v, f = readOff(pathLoad)
    	saveOff(v, f, pathLoad)
    if pathLoad[-3:] == 'obj':
        v, f = readObj(pathLoad)
        saveOff(v, f, pathLoad[:-3] + 'off')
    if pathLoad[-3:] == 'txt':
        t = readTxt(pathLoad)
        saveTxt(t, pathLoad)
