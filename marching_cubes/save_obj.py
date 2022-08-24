#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 14:37:22 2020

@author: shkim
"""

import numpy as np

def saveObj(verts, faces, filepath):
    with open(filepath, "w") as obj_file:
        for i in verts:
            obj_file.write("v {} {} {}\n".format(i[2], i[1], i[0]))
            # obj_file.write("v {} {} {}\n".format(-i[2]+80, i[1]-115, i[0]+660))
        for j in faces:
            obj_file.write("f {} {} {}\n".format(j[0]+1, j[2]+1, j[1]+1))
            
def saveOff(verts, faces, filepath):
    with open(filepath, "w") as off_file:
        off_file.write("OFF\n")
        off_file.write("{} {} {}\n".format(len(verts), len(faces), 0))
        for i in verts:
            off_file.write("{} {} {}\n".format(i[2], i[1], i[0]))
            # off_file.write("{} {} {}\n".format(-i[2]+80, i[1]-115, i[0]+660)) # JW
            # off_file.write("{} {} {}\n".format(-i[2]+115, i[1]-125, i[0]+625)) # H
            # off_file.write("{} {} {}\n".format(-i[2]/4+80, i[1]/4-100, i[0]/4+660))
        for j in faces:
            off_file.write("3 {} {} {}\n".format(j[0], j[2], j[1]))

