import vtk
import os
import sys
from vtk.util import numpy_support
import numpy as np
import math
from vtkmodules.vtkCommonDataModel import vtkStructuredPoints
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def readDS(fname) :
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(fname)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    ds = reader.GetOutput()
    return ds

if __name__ == "__main__":

    folder_name = sys.argv[1]
    
    ds = readDS(file_name)

    # get two field
    point_array_nickel = ds.GetPointData().GetArray("Nickel")
    point_array_iron = ds.GetPointData().GetArray("Iron")

    nparray_n = numpy_support.vtk_to_numpy(point_array_nickel)
    nparray_i = numpy_support.vtk_to_numpy(point_array_iron)
    fig, ax = plt.subplots(figsize=(7,6))
    plt.scatter(nparray_i, nparray_n, s=0.5, alpha=0.1)
    #plt.hist2d(nparray_n,nparray_i,bins=100)
    ax.set_xlabel('Iron',  fontsize='large')
    ax.set_ylabel('Nickel', fontsize='large')
    plt.savefig("supernova_scatter.png",bbox_inches='tight')

