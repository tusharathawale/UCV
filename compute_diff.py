import os
import sys
import vtk
from vtk.util import numpy_support
import numpy as np
import math
from vtkmodules.vtkCommonDataModel import vtkStructuredPoints

def readDS(fname) :
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(fname)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    ds = reader.GetOutput()
    return ds

# compute the RMSE of two input vtk data set
# assuming they have same field name
# the property of the field is point

if __name__ == "__main__":
    
    print("input parameter len", len(sys.argv))
    if (len(sys.argv)!=4):
        print("<script> <file1> <file2> <fieldName>")
        exit(0)
    
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    fieldname = sys.argv[3]

    ds1 = readDS(filename1)
    ds2 = readDS(filename2)

    pointArrayName= fieldname

    pointArray1 = ds1.GetPointData().GetArray(pointArrayName)
    pointArray2 = ds2.GetPointData().GetArray(pointArrayName)

    nparray1 = numpy_support.vtk_to_numpy(pointArray1)
    nparray2 = numpy_support.vtk_to_numpy(pointArray2)

    # compute diff 
    diff_list = []
    max_diff = 0.0
    max_idx = 0.0
    rmse=0.0
    for idx, _ in enumerate(nparray1):
        diff = math.fabs(nparray1[idx]-nparray2[idx])
        rmse = rmse+diff*diff
        diff_list.append(diff)
        # keep track of the id with max diff
        if diff>max_diff:
            max_diff = diff
            max_idx = idx

    #print(np.array(diff_list))
    print("max diff is", max_diff)
    #print("max idx", max_idx)
    print("rmse is", math.sqrt(rmse/(1.0*len(diff_list))))



