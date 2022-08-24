# -*- coding: utf-8 -*-
"""
Created on Wed May 13 12:06:26 2020

@author: user
"""


from functools import partial
import matplotlib.pyplot as plt
from pycpd import RigidRegistration
# import open3d as o3d
import copy
import numpy as np
import os
import os.path as osp
import yaml
            
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

def save_off(verts, faces, filepath):
    with open(filepath, "w") as off_file:
        off_file.write("OFF\n")
        off_file.write("{} {} {}\n".format(len(verts), len(faces), 0))
        for i in verts:
            off_file.write("{} {} {}\n".format(i[0], i[1], i[2]))
        for j in faces:
            off_file.write("3 {} {} {}\n".format(j[0], j[1], j[2]))

def visualize(iteration, error, X, Y, ax):
#    plt.clf()   
#    plt.close()    
    plt.cla()
    ax.scatter(X[:, 0],  X[:, 1], X[:, 2], color='red', label='Target')
    ax.scatter(Y[:, 0],  Y[:, 1], Y[:, 2], color='blue', label='Source')
    ax.text2D(0.87, 0.92, 'Iteration: {:d}\nQ: {:06.4f}'.format(
        iteration, error), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize='x-large')
    ax.legend(loc='upper left', fontsize='x-large')
    plt.draw()
    plt.pause(0.001)

def draw_registration_result(source, target, transformation):
    source_temp = copy.deepcopy(source)
    target_temp = copy.deepcopy(target)
    source_temp.paint_uniform_color([1, 0.706, 0])
    target_temp.paint_uniform_color([0, 0.651, 0.929])
    source_temp.transform(transformation)
    # o3d.visualization.draw_geometries([source_temp, target_temp])

def _load_yaml(file):
    return yaml.safe_load(open(file, 'r'))


def save_txt(content, filepath):
    with open(filepath, 'w') as f:
        for c in content:
            f.write(str(c[0]))
            f.write(' ')
            f.write(str(c[1]))
            f.write(' ')
            f.write(str(c[2]))
            f.write('\n')        

# cancer = np.asarray([342, 153, 134, 0.5], dtype = np.float64) # 375, 135, 147
patients_dir = os.getcwd() + '/../io/landmark'
patients = os.listdir(patients_dir)

# register template to the target
for patient in patients:
    pathLoad = osp.join(patients_dir, patient)
    patientOFF = 'target.off'
    pathLoadOFF = osp.join(patients_dir, patientOFF)
    vertex_target, faces_target = read_off(pathLoadOFF)

    if pathLoad[-3:] == 'off':
        fish_target = np.zeros((19,3)).astype(np.float64)
        with open(pathLoad[:-3] + 'txt', 'r') as f:
            for i in range(19):
                l = f.readline()
                fish_target[i] = np.array(l.split(' '))
            
        # vertex_target, faces_target = read_off(pathLoad)
        # y = _load_yaml(pathLoad[:-3] + 'yaml')
        # fish_target = np.zeros((24,3)).astype(np.int32)
        # landmark = y['landmark']
        # for i in range(24):
        #     fish_target[i] = vertex_target[landmark[i+1]]
        
        vertex_source, faces_source = read_off("templatePCAWide.off")
        y = _load_yaml('templatePCAWide.yaml')
        fish_source = np.zeros((19,3)).astype(np.float64)
        landmark = y['landmark']
        for i in range(19):
            fish_source[i] = vertex_source[landmark[i+1]]
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        callback = partial(visualize, ax=ax)
        
        reg = RigidRegistration(**{'X': fish_source, 'Y': fish_target})
        reg.register(callback)
        vertex_result = reg.transform_point_cloud(vertex_target)
        save_off(vertex_result, faces_target, osp.join(patients_dir, patient[:-4] + 'ToPCA.off'))

# register target to the template
# for patient in patients:
#     pathLoad = osp.join(patients_dir, patient)
#     if pathLoad[-3:] == 'off':
#         vertex_target, faces_target = read_off(pathLoad)    
#         fish_target = np.zeros((24,3)).astype(np.float64)
#         with open(pathLoad[:-3] + 'txt', 'r') as f:
#             for i in range(24):
#                 l = f.readline()
#                 fish_target[i] = np.array(l.split(' '))
        
#         vertex_source, faces_source = read_off("template.off")
#         y = _load_yaml('template.yaml')
#         fish_source = np.zeros((24,3)).astype(np.float64)
#         landmark = y['landmark']
#         for i in range(24):
#             fish_source[i] = vertex_source[landmark[i+1]]
        
#         fig = plt.figure()
#         ax = fig.add_subplot(111, projection='3d')
#         callback = partial(visualize, ax=ax)
        
#         reg = RigidRegistration(**{'X': fish_source, 'Y': fish_target})
#         reg.register(callback)
#         vertex_result = reg.transform_point_cloud(vertex_target)
#         landmark_result = reg.transform_point_cloud(fish_target)
#         save_txt(landmark_result, osp.join(patients_dir, patient[:-4] + '.txt'))
#         save_off(vertex_result, faces_target, osp.join(patients_dir, patient))
#         save_off(vertex_source, faces_source, osp.join(patients_dir, patient[:-4] + '_template.off'))
