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
    print(xdim,ydim,zdim)
    pointArrayNp = numpy_support.vtk_to_numpy(pointArray)
    print("pointArrayNp.shape", pointArrayNp.shape)
    
    for ens in range(numEnsembleMem):
        # using deep copy here
        pointArrayNpEns=np.copy(pointArrayNp)
        # go through each row of pointArray
        # twist x direaction
        print(ens)
        if ens<numEnsembleMem/2:
            print("twist x") 
            for r in range(ydim):
                mu=0
                sigma=0.05
                moveV = np.random.normal(mu, sigma)
                for c in range (xdim):
                    index = r*xdim+c
                    pointArrayNpEns[index]=pointArrayNpEns[index]+moveV
        else:
            # twist y direaction
            print("twist y")
            for c in range (xdim):
                mu=0
                sigma=0.05
                moveV = np.random.normal(mu, sigma)
                for r in range (ydim):
                    index = r*xdim+c
                    pointArrayNpEns[index]=pointArrayNpEns[index]+moveV
    
        output_filename = filename[:-4] + "_" + str(ens) + ".vtk"
        writeVTKDataFromArray(xdim,ydim,zdim,output_filename,pointArrayNpEns)


        

