#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 10:41:15 2020

@author: shkim
"""
import numpy as np
import os
import os.path as osp

def readOff(in_file):
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
    
    faces = np.zeros((num_faces, 4), dtype=np.int32)
    color = np.zeros((num_faces, 3), dtype=np.int16)
    for i in range(num_faces):
        l = src[1+num_nodes+i].split(' ')
        for j in range(4):
            faces[i, j] = int(l[j+1])
        for j in range(3):
            color[i, j] = int(l[j+5])
    return vert, faces, color

def saveOff(verts, faces, color, filepath):
    with open(filepath, "w") as off_file:
        off_file.write("OFF\n")
        off_file.write("{} {} {}\n".format(len(verts), len(faces), 0))
        for i in verts:
            # off_file.write("{} {} {} {} {} {} {}\n".format(i[0], i[1], i[2], i[3], i[4], i[5], i[6]))
            off_file.write("{} {} {}\n".format(i[0], i[1], i[2]))
        for n, j in enumerate(faces):
            off_file.write("4 {} {} {} {} {} {} {}\n".format(j[0], j[1], j[2], j[3], color[n, 0], color[n, 1], color[n,2]))

def append_seq(*args):
    for arg in args:
        assert isinstance(arg, np.ndarray)
    result = args[0]
    for i in range(1, len(args)):
        result = np.append(result, args[i])
    return result

def _load_yaml(file):
    return yaml.safe_load(open(file, 'r'))

if __name__ == "__main__":
    # pathLoad = os.getcwd() + '/../io/mesh/raw/' + 'target' + '.off'
    # pathSave = os.getcwd() + '/../io/mesh/raw/' + patient + '.off'
    # pathSave = os.getcwd() + '/../io/mesh/raw/' + 'target' + '.off'
    pathLoad = 'template_subdiv2.off'
    pathSave = 'template.off'

    verts, faces, color = readOff(pathLoad)

    num2exchange = 500
    
    # v = np.zeros(verts.shape)
    # f = np.zeros(faces.shape, dtype = np.int32)
    
    # v[:num2exchange] = verts[-num2exchange:]
    # v[-num2exchange:] = verts[:num2exchange]
    # v[num2exchange:-num2exchange] = verts[num2exchange:-num2exchange]
    
    # # f = faces
    # for i, j in np.ndindex(f.shape):
    #     if faces[i,j] < num2exchange:
    #         f[i,j] = faces[i,j] + verts.shape[0] - num2exchange
    #     elif faces[i,j] >= verts.shape[0] - num2exchange:
    #         f[i,j] = faces[i,j] - (verts.shape[0] - num2exchange)
    #     else:
    #         f[i,j] = faces[i,j]
    
    # f2 = np.zeros(faces.shape, dtype = np.int32)
    # f2[:10000] = f[-10000:]
    # f2[-10000:] = f[:10000]
    # f2[10000:-10000] = f[10000:-10000]
    
    # c = np.zeros(color.shape, dtype = np.int16)
    # c[:10000] = color[-10000:]
    # c[-10000:] = color[:10000]
    # c[10000:-10000] = color[10000:-10000]
            
    # print(f"verts: {v.shape[0]}, faces: {f.shape[0]}")
    # saveOff(verts, np.flip(faces, axis=0), np.flip(color, axis=0), pathSave)

    # p = np.random.permutation(faces.shape[0])

    # a = np.random.permutation(np.arange(2550,2800))
    # p = append_seq(np.arange(3200,faces.shape[0]),np.arange(1500), a, np.arange(1500,2550), np.arange(2800,3200))


    step = 400
    bins = faces.shape[0] // step
    ps = []
    for i in range(bins):
        ps.append(np.arange(i*step,(i+1)*step))
    ps.append(np.arange(bins*step,faces.shape[0]))
    import random
    random.shuffle(ps)
    p = append_seq(*ps)
    saveOff(verts, faces[p,:], color[p,:], pathSave)

'''
    import yaml
    ldmkLoad = 'template.yaml'
    ldmkSave = 'template_reindexed.yaml'
    y = _load_yaml(ldmkLoad)
    ldmks = y['landmark']
    print(ldmks)
    ldmks_flipped = {}
    for key, value in ldmks.items():
        ldmks_flipped[key] = verts.shape[0] - value
    y = {}
    y['landmark'] = ldmks_flipped
    with open(ldmkSave, 'w') as f:
        yaml.dump(y, f)
'''
