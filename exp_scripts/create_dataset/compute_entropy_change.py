# load files 
# test_syntheticdata_el_sequence_ig_ens1_iso0.80.vtk...test_syntheticdata_el_sequence_ig_ens10_iso0.80.vtk
# compare the change of entropy values
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

if __name__ == "__main__":
    file_dir = "/Users/zw1/Documents/cworkspace/src/UCV/install_scripts/mac/install/UCV"
    num_ensemble=10
    acc_diff=[]
    prev_array=[]
    for i in range(1,num_ensemble+1,1):
        file_name = file_dir+"/"+"test_syntheticdata_el_sequence_ig_ens"+str(i)+"_iso0.80.vtk"
        #print("process file", file_name)
        #load data
        field_name= "entropy0.80"
        ds = readDS(file_name)
        # convert the data into the numpy format
        pointArray = ds.GetCellData().GetArray(field_name)
        pointArrayNp = numpy_support.vtk_to_numpy(pointArray)
        #print("pointArrayNp.shape", pointArrayNp.shape)
        #extract field
        #compare accumulated difference between two files
        #put results into a list
        if len(prev_array)==0:
            # just replace the array and do nothing
            prev_array=pointArrayNp
            continue

            
        # when the prev is not none, comare difference between two array
        diffs = np.abs(pointArrayNp-prev_array)
        #print(diffs)
        diff_sum = np.sum(diffs)
        print("diff sum between first %d and first %d is %f"%(i,i-1,diff_sum))
        # update prev array
        prev_array=pointArrayNp