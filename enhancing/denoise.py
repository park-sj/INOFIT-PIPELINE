#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 15:32:21 2020

@author: shkim
"""

import itk
import SimpleITK as sitk
import os
import os.path as osp
import sys
import numpy as np
import time
from reader import loadDicom



# def array2itkImage(array):
#     image = itk.PyBuffer[ImageType].GetImageFromArray(array.astype(np.float32))
#     return image

# def itkImage2array(image):
#     return itk.GetArrayFromImage(image)
    
# def denoise(image):
#     FilterType = itk.GradientAnisotropicDiffusionImageFilter[ImageType, ImageType]
#     filter = FilterType.New()
#     filter.SetInput(image)
#     filter.SetNumberOfIterations(10)
#     filter.SetTimeStep(1.0)
#     filter.SetConductanceParameter(1.0)
    
#     RescaleFilterType = itk.RescaleIntensityImageFilter[ImageType, ImageType]
#     rescaler = RescaleFilterType.New()
#     rescaler.SetInput(filter.GetOutput())
    
#     pixelTypeMinimum = itk.NumericTraits[PixelType].min()
#     pixelTypeMaximum = itk.NumericTraits[PixelType].max()
    
#     rescaler.SetoutputMinimum(pixelTypeMinimum)
#     rescaler.SetoutputMaximum(pixelTypeMaximum)
    
#     return filter.GetOutput()

def saveDicom(newArray, filepath, templatepath):
    def _load_template():
        dir = templatepath
        assert os.path.isdir(dir), 'Cannot find the template directory'
        reader = sitk.ImageSeriesReader()
        dicomFiles = reader.GetGDCMSeriesFileNames(dir)
        reader.SetFileNames(dicomFiles)
        reader.MetaDataDictionaryArrayUpdateOn()
        reader.LoadPrivateTagsOn()
        reader.Execute()
        imgShape=[len(reader.GetFileNames()),
                  int(reader.GetMetaData(0, "0028|0010")),
                  int(reader.GetMetaData(0, "0028|0011"))]
        # image = reader.Execute()
        # img3d = sitk.GetArrayFromImage(image)
        # img3d = img3d.transpose((1, 2, 0))
        return reader, imgShape

    if not os.path.isdir(filepath):
        os.mkdir(filepath)

    reader, imgShape = _load_template()
    
    newArray = np.squeeze(newArray)
    paddedArray = newArray.astype(np.int16)
    # paddedArray = newArray
    newImage = sitk.GetImageFromArray(paddedArray)
    newImage = sitk.Cast(newImage, sitk.sitkInt16)

    # newImage.CopyInformation(oldImage)
    writer = sitk.ImageFileWriter()
    writer.KeepOriginalImageUIDOn()
    
    sp_x, sp_y = reader.GetMetaData(0, "0028|0030").split('\\')
    # sp_z = reader.GetMetaData(0, "0018|0050")
    _, _, z_0 = reader.GetMetaData(0, "0020|0032").split('\\')
    _, _, z_1 = reader.GetMetaData(1, "0020|0032").split('\\')
    sp_z = abs(float(z_0) - float(z_1))

    modification_time = time.strftime("%H%M%S")
    modification_date = time.strftime("%Y%m%d")
    direction = newImage.GetDirection()
    series_tag_values = [(k, reader.GetMetaData(0, k)) for k in reader.GetMetaDataKeys(0)] + \
                         [("0008|0031", modification_time),
                         ("0008|0021", modification_date),
                         ("0028|0100", "16"),
                         ("0028|0101", "16"),
                         ("0028|0102", "15"),
                         ("0028|0103", "1"),
                         ("0028|0002", "1"),
                         ("0008|0008", "DERIVED\\SECONDARY"),
                         ("0020|000e", "1.2.826.0.1.3680043.2.1125." + modification_date + ".1" + modification_time),
                         ("0020|0037", '\\'.join(map(str, (direction[0], direction[3], direction[6], direction[1], direction[4], direction[7]))))]
#                         ("0008|103e", reader.GetMetaData(0, "0008|103e") + " Processed-SimpleITK")]
#    print(series_tag_values)
    
    tags_to_skip = ['0010|0010', '0028|0030', '7fe0|0010', '7fe0|0000', '0028|1052',
                    '0028|1053', '0028|1054', '0010|4000', '0008|1030', '0010|1001',
                    '0008|0080']
    for i in range(newImage.GetDepth()):
        image_slice = newImage[:, :, i]
        # image_slice.CopyInformation(oldImage[:, :, i])
        for tag, value in series_tag_values:
            if (tag in tags_to_skip):
                continue                
            if i == 0:
                try:
                    print(f'{tag} | {value}')
                except:
                    continue
            image_slice.SetMetaData(tag, value)
        image_slice.SetMetaData("0008|0012", time.strftime("%Y%m%d"))
        image_slice.SetMetaData("0008|0013", time.strftime("%H%M%S"))
        image_slice.SetMetaData('0020|0032', '\\'.join(map(str, [0, 0, i*sp_z])))
        image_slice.SetMetaData("0020|0013", str(i))
        image_slice.SetMetaData('0028|0030', '\\'.join(map(str, [sp_x, sp_y])))
        image_slice.SetSpacing([float(sp_x), float(sp_y)])
        image_slice.SetMetaData("0018|0050", str(sp_z))
        writer.SetFileName(os.path.join(filepath, str(i).zfill(3) + '.dcm'))
        writer.Execute(image_slice)


numberOfIterations = 10
timeStep = 0.0625
conductance = 5.0 # the lower, the more preserving edge values

InputPixelType = itk.F
OutputPixelType = itk.UC
Dimension = 3

InputImageType = itk.Image[InputPixelType, Dimension]
OutputImageType = itk.Image[OutputPixelType, Dimension]

if __name__ == '__main__':
    # _, array, _ = loadDicom('./data/JW_JeongSeon')
    # image = array2itkImage(array)
    # imageDenoised = denoise(image)
    # array2 = itkImage2array(image)
    # arrayDenoised = itkImage2array(imageDenoised)
    
    if len(sys.argv) < 3:
        patients_dir = os.getcwd() + '/../io/temp'
        patients = os.listdir(patients_dir)
        inputImage = osp.join(patients_dir, patients[0])
        outputImage = osp.join(os.getcwd(), '../io/test', patients[0])
    else:
        inputImage = sys.argv[1]
        outputImage = sys.argv[2]
    if not os.path.isdir(outputImage):
        os.mkdir(outputImage)

    
    namesGenerator = itk.GDCMSeriesFileNames.New()
    namesGenerator.SetUseSeriesDetails(True)
    namesGenerator.AddSeriesRestriction("0008|0021")
    namesGenerator.SetGlobalWarningDisplay(False)
    namesGenerator.SetDirectory(inputImage)
    
    seriesUID = namesGenerator.GetSeriesUIDs()
    
    if len(seriesUID) < 1:
        print('No DICOMs in: ' + inputImage)
        sys.exit(1)
    
    print('The directory: ' + inputImage)
    print('Contains the following DICOM Series: ')
    for uid in seriesUID:
        print(uid)
    
    seriesFound = False
    for uid in seriesUID:
        seriesIdentifier = uid
        if len(sys.argv) > 3:
            seriesIdentifier = sys.argv[3]
            seriesFound = True
        print('Reading: ' + seriesIdentifier)
        fileNames = namesGenerator.GetFileNames(seriesIdentifier)
    
        reader = itk.ImageSeriesReader[InputImageType].New()
        dicomIO = itk.GDCMImageIO.New()
        reader.SetImageIO(dicomIO)
        reader.SetFileNames(fileNames)
        reader.ForceOrthogonalDirectionOff()
    
        FilterType = itk.GradientAnisotropicDiffusionImageFilter[
            InputImageType, InputImageType]
        # FilterType = itk.CurvatureAnisotropicDiffusionImageFilter[
        #     InputImageType, InputImageType]
        AnisotropicDiffusionFilter = FilterType.New()
        
        AnisotropicDiffusionFilter.SetInput(reader.GetOutput())
        AnisotropicDiffusionFilter.SetNumberOfIterations(numberOfIterations)
        AnisotropicDiffusionFilter.SetTimeStep(timeStep)
        AnisotropicDiffusionFilter.SetConductanceParameter(conductance)
        
        RescaleFilterType = itk.RescaleIntensityImageFilter[
            InputImageType, OutputImageType]
        rescaler = RescaleFilterType.New()
        rescaler.SetInput(AnisotropicDiffusionFilter.GetOutput())
        
        outputPixelTypeMinimum = itk.NumericTraits[OutputPixelType].min()
        outputPixelTypeMaximum = itk.NumericTraits[OutputPixelType].max()
        
        rescaler.SetOutputMinimum(outputPixelTypeMinimum)
        rescaler.SetOutputMaximum(outputPixelTypeMaximum)
        
        # WriterType = itk.ImageSeriesWriter[OutputImageType, itk.Image[OutputPixelType,3]]
        # writer = WriterType.New()
        # writer.SetImageIO(dicomIO)
        # outputImages = []
        # for i in range(600):
        #     outputImages.append(outputImage +  str(i).zfill(3) + '.dcm')
        # writer.SetFileNames(outputImages)
        # writer.SetInput(rescaler.GetOutput())
        
        # writer.Update()
        
        output = rescaler.GetOutput()
        print("Saving...")
        saveDicom(itk.GetArrayFromImage(output), outputImage, inputImage)
        