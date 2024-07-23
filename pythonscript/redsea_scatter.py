import vtk
import os
import sys
from vtk.util import numpy_support
import numpy as np
import math
from vtkmodules.vtkCommonDataModel import vtkStructuredPoints
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
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

    nparray_mean_curl_file="/Users/g1e/Desktop/vtk-mjay/data/dataForGautam/curlZ/meanVol/meanVorticity.vtk"
    nparray_mean_vort="/Users/g1e/Desktop/vtk-mjay/data/dataForGautam/vorticityMagnitude/meanVol/meanVorticity.vtk"
    
    ds_curl = readDS(nparray_mean_curl_file)
    vtkarray_curl=ds_curl.GetPointData().GetArray("meanVorticity")
    nparray_curl = numpy_support.vtk_to_numpy(vtkarray_curl)

    ds_vort = readDS(nparray_mean_vort)
    vtkarray_vort=ds_vort.GetPointData().GetArray("meanVorticity")
    nparray_vort = numpy_support.vtk_to_numpy(vtkarray_vort)

    fig, ax = plt.subplots(figsize=(6,5))
    plt.scatter(nparray_curl, nparray_vort, s=0.5, alpha=0.1)
    #plt.hist2d(nparray_n,nparray_i,bins=100)
    ax.set_xlabel('curlZ',  fontsize=20)
    ax.set_ylabel('vorticityMagnitude', fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

    # Create a Rectangle patch
    rect = Rectangle((-15,0.6),14.7,14.4,linewidth=2,edgecolor='r',facecolor='none')

    # Add the patch to the Axes
    ax.add_patch(rect)

    plt.savefig("redsea_scatter.png",bbox_inches='tight')

