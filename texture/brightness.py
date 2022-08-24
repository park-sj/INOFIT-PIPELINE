#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 14:08:51 2021

@author: shkim
"""

import matplotlib.pyplot as plt
import numpy as np

def gamma_correction(array, g=0.5):
    return ((array/255)**g)*255

def contrast_stretching(array):
    xp = [0, 36, 64, 128, 192, 255]
    fp = [0, 18, 56, 200, 240, 255]
    x = np.arange(256)
    table = np.interp(x, xp, fp).astype(np.uint8)
    return table[array]

import os
ps = os.listdir('data')
for p in ps:
    if not p.endswith('jpg') or p.endswith('1.jpg'):
        continue
    
    im = plt.imread(os.path.join('data', p), format='jpg')
    im = contrast_stretching(im)

    im = gamma_correction(im.astype(np.float32), 0.55).astype(np.uint8)
    plt.imsave(os.path.join('data', p[:-5]+'1.jpg'), im)