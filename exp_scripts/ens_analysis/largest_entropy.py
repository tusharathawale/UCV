
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

def writeVTKDataFromArray(xdim, ydim, zdim, file_name, field_name, inputArray):
    structured_dataset = vtkStructuredPoints()
    structured_dataset.SetDimensions(xdim, ydim, zdim)
    structured_dataset.SetOrigin(0, 0, 0)
    #print(np.array(g).shape)
    vtkArray = numpy_support.numpy_to_vtk(np.array(inputArray).flatten())
    vtkArray.SetNumberOfComponents(1)
    vtkArray.SetName(field_name)

    structured_dataset.GetPointData().AddArray(vtkArray)
    structured_dataset.GetPointData().SetActiveScalars(field_name)
    print("write file", file_name)
    writeStructuredDs(file_name,structured_dataset)

def get_num_larger_than_3IQR(input_num_list):
    q75 = np.quantile(input_num_list,0.75)
    #print("q 0.75", actual_q75)
    q50 = np.quantile(input_num_list,0.5)
    #print("q 0.5", actual_q50)
    q25 = np.quantile(input_num_list,0.25)
    #print("q 0.25", actual_q25)

    iqr = q75-q25
    higher_bound = q50+1.5*iqr

    num_elem_larger_than_3iqr=0
    for v in input_num_list:
        if v>higher_bound:
            num_elem_larger_than_3iqr=num_elem_larger_than_3iqr+1

    return num_elem_larger_than_3iqr

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
    num_outlayer=[]
    elem_num = len(all_entroly_diff_list[0])
    for element_index in range(0,elem_num,1):
        colum = [row[element_index] for row in all_entroly_diff_list]
        #check the largest value for current row
        if max(colum)>0.01:
            index_max = np.argmax(colum)
            largest_contribute_num.append(index_max)
        else:
            largest_contribute_num.append(-1)

        #check the 3 iqr value
        num_outlayer_for_colm=get_num_larger_than_3IQR(colum)
        num_outlayer.append(num_outlayer_for_colm)

    #print(largest_contribute_num)
    writeVTKDataFromArray(60,80,1,"max_entropy_id.vtk","MaxEntropyDiffId",largest_contribute_num)
    writeVTKDataFromArray(60,80,1,"num_outlayer.vtk","NumOfOutlayer", num_outlayer)



