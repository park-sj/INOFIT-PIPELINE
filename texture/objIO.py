#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 31 15:18:29 2021

@author: shkim
"""


class MyObjIO():
    ''' 
    Obj IO jobs here are somewhat peculiar and task-specific
    You must not reuse this code for general jobs
    To avoid possible error, I designed IO function as custumized class
    Some faces don't have texture property and we won't do anything with faces and vts
    So, we load only the verts as number and load others as just string
    '''
    @classmethod
    def load(cls, in_file):
        verts = []
        faces = []
        texts = []
        mat = []
    
        for k, v in cls._yield_file(in_file):
            if k == 'v':
                verts.append(v)
            elif k == 'f':
                faces.append(v)
            elif k == 'vt':
                texts.append(v)
            elif k == 'm':
                mat.append(v)
    
        if not len(faces) or not len(verts):
            return None, None
    
        return verts, faces, texts, mat
    
    @classmethod
    def save(cls, verts, faces, texts, mat, out_file):
        with open(out_file, "w") as obj_file:
            for m in mat:
                obj_file.write(m+'\n')
            for i in verts:
                obj_file.write("v {} {} {}\n".format(i[0], i[1], i[2]))
            obj_file.write('usemtl material_0\n')
            for j in texts:
                obj_file.write(j+'\n')
            for k in faces:
                obj_file.write(k+'\n')
    
    @staticmethod
    def _yield_file(in_file):
        f = open(in_file)
        buf = f.read()
        f.close()
        for b in buf.split('\n'):
            b = b.strip()
            if b.startswith('v '):
                yield ['v', [float(x) for x in b.split()[1:]]]
            elif b.startswith('vt '):
                yield ['vt', b]
            elif b.startswith('f '):
                yield ['f', b]
            elif b.startswith('mtllib '):
                yield ['m', b]
            else:
                yield ['', ""]
