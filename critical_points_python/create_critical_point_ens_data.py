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

# create input data
foo = flatpy.nD.available_functions["ackley"]
fractional_noise_level=20
persistence = 0.665
n_clusters = 9
#number of ensemble data
count = 3
noise_level = 0.01*persistence*fractional_noise_level
noise_model = "uniform"
ground_truth, ensemble = utpy.utils.generate_ensemble(foo, noise_level, count, noise_model)
print(np.min(ground_truth))
img = plt.imshow(ensemble[:,:,1])
plt.colorbar(img)
plt.savefig("demo_critial_point.png")

print("all ens shape", ensemble.shape)
print("each ens shape", ensemble[:,:,0].shape)

dir_name="critical_point_ens_data_5by5"

isExist = os.path.exists(dir_name)
if not isExist:
    os.mkdir(dir_name)

#xdim=50
#ydim=50

xdim=5
ydim=5
zdim=1

for ensid in range(count):
    print("create ensid", ensid)
    output_filename=dir_name+"/critical_point_ens_5by5_"+str(ensid)+".vtk"
    writeVTKDataFromArray(xdim,ydim,zdim,output_filename,ensemble[:,:,ensid])
