#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 10:37:01 2020

@author: shkim
"""

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

