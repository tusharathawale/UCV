import vtk
import os
import sys
from vtk.util import numpy_support
import numpy as np
import math
from vtkmodules.vtkCommonDataModel import vtkStructuredPoints
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle


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


    #go through all ens and compute the average values
    ens_num=64
    nparray_summ_n=None
    nparray_summ_i=None
    for i in range(0,ens_num,1):
        file_name = folder_name+"/ens_"+str(i)+".vtk"
        print("load file",file_name)

        ds = readDS(file_name)

        # get two field
        point_array_nickel = ds.GetPointData().GetArray("Nickel")
        point_array_iron = ds.GetPointData().GetArray("Iron")

        nparray_n = numpy_support.vtk_to_numpy(point_array_nickel)
        nparray_i = numpy_support.vtk_to_numpy(point_array_iron)

        if i==0:
            nparray_summ_n=np.zeros(nparray_n.size)
            nparray_summ_i=np.zeros(nparray_i.size)
    
        nparray_summ_n=np.add(nparray_summ_n,nparray_n)
        nparray_summ_i=np.add(nparray_summ_i,nparray_i)

    fig, ax = plt.subplots(figsize=(6,5))
    plt.scatter(np.divide(nparray_summ_i, ens_num), np.divide(nparray_summ_n,ens_num), s=0.5, alpha=0.1)
    #plt.hist2d(nparray_n,nparray_i,bins=100)
    ax.set_xlabel('Iron',  fontsize=20)
    ax.set_ylabel('Nickel', fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

    # Create a Rectangle patch
    rect = Rectangle((0.1,0.1),0.4,0.2,linewidth=2,edgecolor='r',facecolor='none')

    # Add the patch to the Axes
    ax.add_patch(rect)

    plt.savefig("supernova_scatter.png",bbox_inches='tight')

