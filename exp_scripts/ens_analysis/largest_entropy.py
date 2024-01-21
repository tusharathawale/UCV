
# use ensemble_contour_output_with_no_ens as input
# create a new ens data
# for each ens member, show the ensemble id with the largest entropy value
# need to do the sorting in this step
# the largest entropy difference is caused by which ensemble member

import numpy as np
import vtk
import os
import sys
from vtk.util import numpy_support
import matplotlib.pyplot as plt
from vtkmodules.vtkCommonDataModel import vtkStructuredPoints

def readDS(fname) :
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(fname)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    ds = reader.GetOutput()
    return ds

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
    vtkArray.SetName("MaxEntropyEnsId")

    structured_dataset.GetPointData().AddArray(vtkArray)
    structured_dataset.GetPointData().SetActiveScalars("MaxEntropyEnsId")
    print("write file", file_name)
    writeStructuredDs(file_name,structured_dataset)

if __name__ == "__main__":
    filename = sys.argv[1]
    # load data set
    ds=readDS(filename)

    suffix = "entropy_diff_no_ens_"
    all_entroly_diff_list=[]
    for id in range(0,20,1):
        field_name=suffix+str(id)
        # load field
        field_array = ds.GetCellData().GetArray(field_name)
        all_entroly_diff_list.append(numpy_support.vtk_to_numpy(field_array))
        print("load diff no ens for", id)

    #print(len(all_entroly_diff_list), len(all_entroly_diff_list[0]))

    # for each row in the entropy list
    largest_contribute_num=[]
    elem_num = len(all_entroly_diff_list[0])
    for element_index in range(0,elem_num,1):
        colum = [row[element_index] for row in all_entroly_diff_list]
        #check the largest value for current row
        if max(colum)>0.01:
            index_max = np.argmax(colum)
            largest_contribute_num.append(index_max)
        else:
            largest_contribute_num.append(-1)

    #print(largest_contribute_num)
    writeVTKDataFromArray(60,80,1,"max_entropy_id.vtk",largest_contribute_num)



