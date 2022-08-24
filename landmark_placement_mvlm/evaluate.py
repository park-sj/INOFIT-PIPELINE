#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 11:01:21 2020

@author: shkim
"""

import os
import glob
import numpy as np

def save_txt(content, filepath):
    with open(filepath, 'w') as f:
        for c in content:
            for i, j in enumerate(c):
                f.write(str(j))
                if i != 2:
                    f.write(' ')
            f.write('\n')

def read_txt(filepath):
    print(filepath)
    with open(filepath, 'r') as f:
        c = np.zeros((24,3))
        lines = f.readlines()
        for i in range(len(lines)):
            l = lines[i]
            # print(l)
            c[i] = list(map(float, l.split(' ')))
    return c

if __name__ == "__main__":
    basedir = 'assets/CTMASK/'
    files = ['ShinHeaJin', 'ShinYeRin', 'WonWuBin', 'WuBoRam', 'YangYeongSuk', 'YuHoJeong', 'YunJeongSuk']
    # files = os.listdir(basedir)
    # files = ['assets/CTMASK/LeeGangHu.yaml']
    dist = 0
    for f in files:
        f = basedir + f
        gt = read_txt(f + '.txt')
        ev = read_txt(f + '_landmarks.txt')
        d = np.power(gt - ev, 2)
        d = np.sum(d, axis=-1)
        d = np.sqrt(d)
        d = np.sum(d)
        dist += d / gt.shape[0]
    dist = dist / len(files)