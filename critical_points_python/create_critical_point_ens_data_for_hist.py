from matplotlib import colors, cm
import matplotlib.pyplot as plt

import numpy as np
import math
import time
import utpy.utils
import utpy.vis
import flatpy
import os

import vtk
from vtk.util import numpy_support
from vtkmodules.vtkCommonDataModel import vtkStructuredPoints

def writeStructuredDs(fname, ds):
    writer = vtk.vtkStructuredPointsWriter()
    writer.SetFileName(fname)
    writer.SetFileVersion(42)
    writer.SetInputData(ds)
    writer.Update()
    writer.Write() 

def writeVTKDataFromArray(xdim, ydim, zdim, file_name, inputArray):
    structured_dataset = vtkStructuredPoints()
    structured_dataset.SetDimensions(xdim, ydim, zdim)
    structured_dataset.SetOrigin(0, 0, 0)
    #print(np.array(g).shape)
    vtkArray = numpy_support.numpy_to_vtk(np.array(inputArray).flatten())
    vtkArray.SetNumberOfComponents(1)
    vtkArray.SetName("TestField")

    structured_dataset.GetPointData().AddArray(vtkArray)
    structured_dataset.GetPointData().SetActiveScalars("TestField")
    print("write file", file_name)
    writeStructuredDs(file_name,structured_dataset)

# define normalized 2D gaussian
def gaus2d(x=0, y=0, mx=0, my=0, sx=1, sy=1):
    return 1. / (2. * np.pi * sx * sy) * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.)))


def add_Gaussian_noise(field, noise_level, numMembers):
    h,w = field.shape
    ensemble = np.zeros((h,w,numMembers))
    
    for i in range(numMembers):
        epsilon = np.random.normal(0,noise_level,[h,w])
        ensemble[:,:,i] = field + epsilon
        
    return ensemble

x = np.linspace(-5, 5)
y = np.linspace(-5, 5)
x, y = np.meshgrid(x, y) # get 2D variables instead of 1D
z1 = gaus2d(x, y, 2.5, 2.5, 0.3,0.3)
z2 = gaus2d(x, y, -2.5, -2.5, 0.3,0.3)
truthDataset = -(z1+z2)
img = plt.imshow(truthDataset)
plt.colorbar(img)
plt.savefig("demo_truth_dataset.png")

plt.clf()

x = np.linspace(-5, 5)
y = np.linspace(-5, 5)
x, y = np.meshgrid(x, y) # get 2D variables instead of 1D
z1 = gaus2d(x, y, -2.5, 2.5, 0.3, 0.3)
z2 = gaus2d(x, y, 2.5, -2.5, 0.3, 0.3)
rotatedOutlierDataset = -(z1+z2)
img = plt.imshow(rotatedOutlierDataset)
plt.colorbar(img)
plt.savefig("demo_rotate_dataset.png")

plt.clf()

noise_level = 0.001

truthEnsembleDataset = add_Gaussian_noise(truthDataset,noise_level,40)
plt.figure()
img = plt.imshow(truthEnsembleDataset[:,:,0])
plt.colorbar(img)
print(truthEnsembleDataset.shape)
plt.savefig("demo_truth_dataset_gauss_noise.png")

rotatedOutlierEnsembleDataset = add_Gaussian_noise(rotatedOutlierDataset,noise_level,10)
plt.figure()
img = plt.imshow(rotatedOutlierEnsembleDataset[:,:,0])
plt.colorbar(img)
plt.savefig("demo_rotate_dataset_gauss_noise.png")

ensemble = np.concatenate((truthEnsembleDataset, rotatedOutlierEnsembleDataset), axis=2)
print("ensemble shape", ensemble.shape)

plt.clf()
img = plt.imshow(ensemble[:,:,1])
plt.colorbar(img)
plt.savefig("demo_critial_point_hist.png")

print("all ens shape", ensemble.shape)
print("each ens shape", ensemble[:,:,0].shape)

dir_name="critical_point_ens_data_hist_50by50"

isExist = os.path.exists(dir_name)
if not isExist:
    os.mkdir(dir_name)

#xdim=50
#ydim=50

xdim=50
ydim=50
zdim=1
count=50

for ensid in range(count):
    #print("create ensid", ensid)
    output_filename=dir_name+"/critical_point_ens_data_hist_50by50_"+str(ensid)+".vtk"
    print("ensid",ensid)
    #print(ensemble[:,:,ensid])
    writeVTKDataFromArray(xdim,ydim,zdim,output_filename,ensemble[:,:,ensid])
