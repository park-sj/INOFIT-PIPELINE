#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 17:15:44 2021

@author: shkim
"""

import numpy as np

def calculate_stat(dist, th):
    portion = len(dist[dist>th]) / len(dist)
    return np.mean(dist), np.var(dist), np.max(dist), portion