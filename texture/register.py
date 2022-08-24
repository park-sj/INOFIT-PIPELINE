# -*- coding: utf-8 -*-
"""
Created on Thu May 20 16:29:54 2021

@author: user
"""

import sys
import os
import numpy as np
from scipy.spatial.transform import Rotation as R
from objIO import MyObjIO

def align(verts, transform):
    def read(transform):
        f = open(transform)
        lines = f.read().split('\n')
        translation = None # scanner mesh doesn't have translation record
        for line in lines:
            line = line.split(' ')
            if line[0] == 'Rotation':
                rotation = list(map(float, line[2].strip('[').strip(']').split(',')))
            elif line[0] == 'Translation':
                translation = list(map(float, line[2].strip('[').strip(']').split(',')))
        assert 'rotation' in locals(), 'Could not read align txt'
        return rotation, translation
    
    def translate(verts, translation):
        verts += np.array(translation)
        return verts
    
    def rotate(verts, rotation):
        r = R.from_euler('xyz', rotation, degrees=True)
        verts = r.apply(verts)
        return verts
    
    def rotate_offset(verts):
        with open(os.path.join(os.getcwd(),'..','mesh_alignment','scanner_offset.txt'), 'r') as f:
            offset = list(map(float, f.read().split(' ')))
        verts = rotate(verts, offset)
        return verts
    
    rotation, _ = read(transform)
    verts = rotate_offset(verts)
    verts = rotate(verts, rotation)
    return verts


if __name__ == "__main__":
    if len(sys.argv) == 4:
        loadPath = sys.argv[1]
        savePath = sys.argv[2]
        alignPath = sys.argv[3]
    else:
        loadPath = '/home/shkim/Libraries/texture_transfer/data/sb.obj'
        savePath = '/home/shkim/Libraries/texture_transfer/data/sb_test.obj'
        alignPath = 'align.txt'
    verts, faces, texts, mat = MyObjIO.load(loadPath)
    verts = align(verts, alignPath)
    MyObjIO.save(verts, faces, texts, mat, savePath)
