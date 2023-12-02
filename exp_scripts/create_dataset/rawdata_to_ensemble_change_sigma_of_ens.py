import numpy as np
import vtk
import os
import sys
from vtk.util import numpy_support
import matplotlib.pyplot as plt
from vtkmodules.vtkCommonDataModel import vtkStructuredPoints

def writeStructuredDs(fname, ds):
    writer = vtk.vtkStructuredPointsWriter()
    writer.SetFileName(fname)
    writer.SetFileVersion(42)
    writer.SetInputData(ds)
    writer.Update()
    writer.Write() 

def readDS(fname) :
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(fname)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    ds = reader.GetOutput()
    return ds

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

# load the point data
if __name__ == "__main__":
    filename = sys.argv[1]
    fieldarray= "TestField"
    numEnsembleMem=20

    ds = readDS(filename)
    # convert the data into the numpy format
    pointArray = ds.GetPointData().GetArray(fieldarray)
    xdim,ydim,zdim=ds.GetDimensions()
    pointArrayNp = numpy_support.vtk_to_numpy(pointArray)
    print("pointArrayNp.shape", pointArrayNp.shape)

    # go through the data, 
    # at each point, extract ensemble memebrs by adding some noises 
    ensemble_all=[]
    for p in pointArrayNp:
        # mu is current point, sigma is 1
        mu1=p
        sigma1=0.02

        mu2=p+0.5
        sigma2=0.03

        # TODO, maybe adding more gaussain for testing
        # mixture, how well
        # shuffling of members
        # clustering + outline

        # let half of values from distribution1 and another half values from distribution2
        ensembles_each_point_first_half = np.random.normal(mu1, sigma1, 10)
        ensembles_each_point_second_half = np.random.normal(mu2, sigma2, 10)

        ensembles_each_point= np.concatenate((ensembles_each_point_first_half,ensembles_each_point_second_half),axis=0)

        ensemble_all.append(ensembles_each_point)
        
    # write out the ensemble members into a separate files.
    # plt.hist(ensemble_all[0],bins=20)
    # plt.savefig("test_ensemble_data.png")
    print("ensemble_all num points", len(ensemble_all), "num of ensemble for each point", len(ensemble_all[0]))
    
    # for each ensemble member
    # create a list for it
    for i in range(len(ensemble_all[0])):
        print("output ensemble with",i)
        ensemble_list_local=[]
        output_filename = filename[:-4] + "_" + str(i) + ".vtk"
        for j in range(len(ensemble_all)):
            ensemble_list_local.append(ensemble_all[j][i])

        writeVTKDataFromArray(xdim,ydim,zdim,output_filename,ensemble_list_local)
        

