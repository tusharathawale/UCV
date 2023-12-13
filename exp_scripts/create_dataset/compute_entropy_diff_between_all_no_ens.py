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
from sklearn.cluster import KMeans

def readDS(fname) :
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(fname)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    ds = reader.GetOutput()
    return ds

if __name__ == "__main__":
    #file_dir = "/Users/zw1/Documents/cworkspace/src/UCV/install_scripts/mac/install/UCV"
    
    acc_diff=[]
    prev_array=[]

    if(len(sys.argv)!=4):
        print("<executable> <fieldSuffix> <iso> <num_ens>")
        exit(0)
        
    fieldSuffix = str(sys.argv[1])
    isovalue = str(sys.argv[2])
    num_ensemble=int(sys.argv[3])

    # get case that use all ens
    file_name_all_ens = fieldSuffix+"_using_all_ens_iso_"+isovalue+".vtk"
    ds_all_ens = readDS(file_name_all_ens)
    entropy_field_name= "entropy"+isovalue
    pointArray = ds_all_ens.GetCellData().GetArray(entropy_field_name)
    pointArrayAllEnsNp = numpy_support.vtk_to_numpy(pointArray)
    
    diff_sum_list=[]
    for i in range(0,num_ensemble,1):
        #file_name = file_dir+"/"+"test_syntheticdata_el_sequence_ig_ens"+str(i)+"_iso0.80.vtk"
        #file_name = file_dir+"/"+"test_2ddata_el_using_ens_"+str(i)+"_iso0.80.vtk"
        
        #file_name = file_dir+"/"+"test_2ddata_el_velocityMagnitude_using_ens_"+str(i)+"_iso0.10.vtk"
        file_name = fieldSuffix+"_no_ens_"+str(i)+"_iso_"+isovalue+".vtk"
        #print("process file", file_name)
        #load data
        #field_name= "entropy0.80"
        #for red sea
        field_name= "entropy"+isovalue
        ds = readDS(file_name)
        # convert the data into the numpy format
        pointArray = ds.GetCellData().GetArray(field_name)
        pointArrayNp = numpy_support.vtk_to_numpy(pointArray)

        # when the prev is not none, comare difference between two array
        diffs = np.abs(pointArrayNp-pointArrayAllEnsNp)
        #print(diffs)
        diff_sum = np.sum(diffs)
        diff_sum_list.append(diff_sum)
        print("diff sum between all and no ens %d is %f"%(i,diff_sum))

no_ens_id_list=list(range(0, num_ensemble,1))
plt.plot(diff_sum_list, np.zeros_like(diff_sum_list) + 0, 'x')
plt.savefig("no_ens_diff_1d.png")
plt.clf()
plt.scatter(no_ens_id_list, diff_sum_list)
plt.savefig("no_ens_diff_2d.png")

# TODO Using GMM to do the classification for the diff_sum_list
# using kmeans to get the classification
# km = KMeans()
# km.fit(diff_sum_list)  # -1 will be calculated to be 13876 here
# print(km)

# Create some data
# data = np.random.randn(10)
# make sure the column is 1 for the output data set
print("diff_sum_list before reshape: ", diff_sum_list)
data_reshaped = np.array(diff_sum_list).reshape(-1,1)
print(data_reshaped)

# Create a KMeans object
kmeans = KMeans(n_clusters=3)

# Fit the model to the data
kmeans.fit(data_reshaped)

# Predict the cluster labels for each data point
labels = kmeans.predict(data_reshaped)

# Print the cluster labels
print(labels)
