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

def saveOff(verts, faces, filepath):
    with open(filepath, "w") as off_file:
        off_file.write("OFF\n")
        off_file.write("{} {} {}\n".format(len(verts), len(faces), 0))
        for i in verts:
            off_file.write("{} {} {}\n".format(i[0], i[1], i[2]))
            # off_file.write("{} {} {}\n".format(-i[2]+80, i[1]-115, i[0]+660)) # JW
            # off_file.write("{} {} {}\n".format(-i[2]+115, i[1]-125, i[0]+625)) # H
            # off_file.write("{} {} {}\n".format(-i[2]/4+80, i[1]/4-100, i[0]/4+660))
        for j in faces:
            off_file.write("3 {} {} {}\n".format(j[0], j[1], j[2]))


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

if __name__ == "__main__":
    pathLoad = os.getcwd() + '/../io/mesh/raw/' + 'target' + '.off'
    # pathSave = os.getcwd() + '/../io/mesh/raw/' + patient + '.off'
    pathSave = os.getcwd() + '/../io/mesh/raw/' + 'target' + '.off'

    verts, faces = readOff(pathLoad)

    # saveOff(verts, faces, pathSave)
    # mesh = trimesh.Trimesh(verts, faces)
    mesh = trimesh.load(pathLoad, process=False)
    # trimesh.smoothing.filter_laplacian(mesh, iterations=1)
    v, f = trimesh.sample.sample_surface(mesh, 100000)
    mesh = mesh_simplication(mesh)
    f = mesh.faces
    v = mesh.vertices
    print(f"verts: {v.shape[0]}, faces: {f.shape[0]}")

    saveOff(v, f, pathSave)
