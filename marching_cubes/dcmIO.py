#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 10:45:24 2020

@author: shkim
"""

import SimpleITK as sitk

def loadDcm(filepath):
    series_IDs = sitk.ImageSeriesReader.GetGDCMSeriesIDs(filepath)
    if not series_IDs:
        print("No DICOM here :(")
        return
    series_file_names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(filepath, series_IDs[0])
    series_reader = sitk.ImageSeriesReader()
    series_reader.SetFileNames(series_file_names)
    series_reader.MetaDataDictionaryArrayUpdateOn()
    series_reader.LoadPrivateTagsOn()
    image = series_reader.Execute()
    img3d = sitk.GetArrayFromImage(image)
#    img3d = img3d.transpose((1,2,0))
    return img3d, image, series_reader