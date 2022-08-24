#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 10:41:15 2020

@author: shkim
"""
import numpy as np
import os
import os.path as osp
import trimesh
import openmesh as om
from skimage import measure, filters

from save_obj import saveObj, saveOff
from dcmIO import loadDcm

def tri_to_om(mesh):
    om_mesh = om.TriMesh(np.array(mesh.vertices))
    om_mesh.add_faces(np.array(mesh.faces))
    return om_mesh

def om_to_tri(om_mesh):
    mesh = trimesh.Trimesh()
    mesh.vertices = om_mesh.points()
    mesh.faces = om_mesh.face_vertex_indices()
    return mesh

def mesh_simplication(mesh):
    om_mesh = tri_to_om(mesh)
    mod = om.TriMeshModQuadricHandle()
    decimator = om.TriMeshDecimater(om_mesh)
    decimator.add(mod)
    mod = decimator.module(mod)
    mod.set_binary(False)
    mod.set_max_err(1e-3)
    decimator.initialize()
    decimator.decimate()
    om_mesh.garbage_collection()
    return om_to_tri(om_mesh)

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


if __name__ == '__main__':
    verts, faces = readOff('../target.off')
    um = None
    idx = -1
    for i, v in enumerate(verts):
        if um is None or um > v[1]:
            um = v[1]
            idx = i
    print(idx)
    print(verts[idx])
    
    v_to_del = []
    for i, v in enumerate(verts):
        if v[2] > verts[idx][2] + 100:
            v_to_del.append(i)
    
    f_to_del = []
    for i, f in enumerate(faces):
        for v in v_to_del:
            if v in f:
                f_to_del.append(v)