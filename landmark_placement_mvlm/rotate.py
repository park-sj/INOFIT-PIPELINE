#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 14:48:56 2020

@author: shkim
"""

import argparse
import os

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
    vertices = []
    faces = []

    for k, v in yield_file(in_file):
        if k == 'v':
            vertices.append(v)
        elif k == 'f':
            for i in v:
                faces.append(i)

    if not len(faces) or not len(vertices):
        return None

    return vertices, faces

def save_obj(verts, faces, filepath):
    with open(filepath, "w") as obj_file:
        for i in verts:
            obj_file.write("v {} {} {}\n".format(i[0], i[1], i[2]))
        for j in faces:
            obj_file.write("f {} {} {}\n".format(j[0]+1, j[1]+1, j[2]+1))

parser = argparse.ArgumentParser(description='rotate')
# parser.add_argument('--l', type=str, default='/home/shkim/Libraries/pytorch_geometric/data/CTMASK/raw/CTMASK2/YuHoJeong.obj')
# parser.add_argument('--l', type=str, default='/home/shkim/Libraries/ACSCNN/landmark/build_dataset/YuHoJeong.obj')
parser.add_argument('--l', type=str, default=os.path.join(os.getcwd(), 'assets', 'BaeYuSeok.obj'))
parser.add_argument('--s', type=str, default=os.path.join(os.getcwd(), 'assets', 'rotated.obj'))
args = parser.parse_args()

origVerts, origFaces = read_obj(args.l)
# newVerts = origVerts
newVerts = []
for v in origVerts:
    newVerts.append([(-v[2]+600)/4, v[0]/4, (-v[1]+600)/4])
    # newVerts.append([(v[1]+600)/10, v[0]/10, (v[2]+600)/10])
save_obj(newVerts, origFaces, args.s)
