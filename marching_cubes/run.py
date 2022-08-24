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

patients_dir = os.getcwd() + '/../io/mask'
# patients_dir = '/home/shkim/Libraries/pytorch-3dunet/datasets/JW/test_masks'
# patients = ['LeeGangHu']
patients = os.listdir(patients_dir)

for patient in patients:
    pathLoad = osp.join(patients_dir, patient)
    try:
        os.mkdir(os.getcwd() + '/../io/mesh/raw/mesh/')
    except FileExistsError:
        pass
    # pathSave = os.getcwd() + '/../io/mesh/raw/' + patient + '.off'
    pathSave = os.getcwd() + '/../io/mesh/raw/' + 'target' + '.off'

    mask, image, series_reader = loadDcm(pathLoad)
    sp_x, sp_y = series_reader.GetMetaData(0, "0028|0030").split('\\')
    sp_z = series_reader.GetMetaData(0, "0018|0050")
    spacing = tuple(map(float, [sp_z, sp_x, sp_y]))
    # print(spacing)
    # mask = filters.laplace(mask)
    # mask = mask[:,:,200:-200]
    verts, faces, _, _ = measure.marching_cubes_lewiner(mask,
                                                        spacing=spacing,
                                                        step_size=1,
                                                        gradient_direction='ascent',
                                                        allow_degenerate=False)
    # saveOff(verts, faces, pathSave)
    mesh = trimesh.Trimesh(verts, faces)
    trimesh.smoothing.filter_laplacian(mesh)
    v, f = trimesh.sample.sample_surface(mesh, 10000)
    mesh = mesh_simplication(mesh)
    f = mesh.faces
    v = mesh.vertices
    print(f"{patient} - verts: {v.shape[0]}, faces: {f.shape[0]}, spacing: {spacing}")
    saveOff(v, f, pathSave)