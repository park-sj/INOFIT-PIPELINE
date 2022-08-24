#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:48:02 2020

@author: shkim
"""

import numpy as np
import sys

def yield_file(in_file):
    f = open(in_file)
    buf = f.read()
    f.close()
    for b in buf.split('\n'):
        b = b.strip()
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
            off_file.write("{} {} {}\n".format(i[0], i[1], i[2]))
        for j in faces:
            off_file.write("3 {} {} {}\n".format(j[0], j[1], j[2]))

if __name__ == "__main__":
    if (len(sys.argv) != 3):
        raise Exception("Wrong number of argv")
    print(f"Convert {sys.argv[1]} to {sys.argv[2]}")
    vert, faces = readObj(sys.argv[1])
    saveOff(vert, faces, sys.argv[2])