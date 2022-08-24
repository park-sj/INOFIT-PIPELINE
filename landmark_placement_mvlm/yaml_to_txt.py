#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:33:37 2020

@author: shkim
"""
import numpy as np
import yaml
import os
import glob

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

def read_obj(in_file):
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
    
    faces = np.zeros((num_faces, 3), dtype=np.int32)
    for i in range(num_faces):
        l = src[1+num_nodes+i].split(' ')
        faces[i, 0] = int(l[1])
        faces[i, 1] = int(l[2])
        faces[i, 2] = int(l[3])
    return vertices, faces



def read_yaml(filepath):
    list_lm = []
    lm = yaml.safe_load(open(filepath, 'r'))
    lm = lm['landmark']
    for i in lm.keys():
        list_lm.append(lm[i])
    return list_lm

def save_txt(content, filepath):
    with open(filepath, 'w') as f:
        for c in content:
            for i, j in enumerate(c):
                f.write(str(j))
                if i != 2:
                    f.write(' ')
            f.write('\n')

if __name__ == "__main__":
    basedir = '../io/landmark/'
    #basedir = 'assets/CTMASK/'
    files = os.listdir(basedir)
    # files = ['assets/CTMASK/LeeGangHu.yaml']
    for f in files:
        f = basedir + f
        if f[-4:] != 'yaml':
            continue
        lm = read_yaml(f)
        #v, _ = read_obj(f[:-4] + 'obj')
        v, _ = read_off(f[:-4] + 'off')
        c = []
        for l in lm:
            c.append(v[l])
        save_txt(c, f[:-4] + 'txt')
    

    
