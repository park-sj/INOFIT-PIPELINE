#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 14:42:05 2020

@author: shkim
"""

import SimpleITK as sitk

def loadDicom(filepath):
    # dir = '/home/shkim/SpyderProjects/0611/'
    # dir = '/home/shkim/Libraries/pytorch-3dunet/datasets/JW/test/'
    # dir = '/home/shkim/Libraries/pytorch-3dunet/datasets/H/test/'
    # patient = os.listdir(dir)
    # dir = os.path.join(dir, patient[0])
    reader = sitk.ImageSeriesReader()
    dicomFiles = reader.GetGDCMSeriesFileNames(filepath)
    reader.SetFileNames(dicomFiles)
    reader.MetaDataDictionaryArrayUpdateOn()
    reader.LoadPrivateTagsOn()
    reader.Execute()
    imgShape=[len(reader.GetFileNames()),
              int(reader.GetMetaData(0, "0028|0010")),
              int(reader.GetMetaData(0, "0028|0011"))]
    image = reader.Execute()
    img3d = sitk.GetArrayFromImage(image)
    # img3d = img3d.transpose((1, 2, 0))
    return reader, img3d, imgShape
