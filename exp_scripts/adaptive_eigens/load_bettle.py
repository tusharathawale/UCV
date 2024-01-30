
import vtk
from vtk.util import numpy_support
import os
import sys


def readDS(fname) :
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(fname)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    ds = reader.GetOutput()
    return ds

if __name__ == "__main__":
    filename = sys.argv[1]
    fieldarray= "ground_truth"
    ds = readDS(filename)
    #print(ds)

    #get vtkarray
    point_array = ds.GetPointData().GetArray(fieldarray)
    #print(point_array)

    #get numpy array
    nparray = numpy_support.vtk_to_numpy(point_array)

    print(nparray.shape)

    # TODO, operations based on nparray