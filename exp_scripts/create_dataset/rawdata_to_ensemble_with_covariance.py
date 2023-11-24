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
    numEnsembleMem=80

    ds = readDS(filename)
    # convert the data into the numpy format
    pointArray = ds.GetPointData().GetArray(fieldarray)
    xdim,ydim,zdim=ds.GetDimensions()
    pointArrayNp = numpy_support.vtk_to_numpy(pointArray)
    print("pointArrayNp.shape", pointArrayNp.shape)

    # go through the data, 
    # at each point, extract ensemble memebrs by adding some noises 
    ensemble_all=[]
    for p in pointArrayNp:
        # mu is current point, sigma is 1
        mu=p
        sigma=0.02
        # use multimodal distribution to add noise
        # create numEnsembleMem values from uniform distribution
        uniform_samples = np.random.uniform(0,1,numEnsembleMem)
        # for each sample, if it is larger than 0.8 use first Gaussian
        ensembles_values = []
        for s in uniform_samples:
            if s<0.6:
                ensembles_values.append(np.random.normal(mu, sigma))
            else:
                # Try large gap to get some vis results
                # refer to the probabalisitc mc paper, Figure6, changing the correlation can change vis effects
                # ensembles_values.append(np.random.normal(mu+5, sigma))
                ensembles_values.append(np.random.normal(mu+5, sigma))

        ensemble_all.append(ensembles_values)
    # write out the ensemble members into a separate files.
    my_formatted_list = [ float('%.5f' % elem) for elem in ensemble_all[0] ]
    print ("ensemble_all[0]",my_formatted_list)

    plt.hist(ensemble_all[0],bins=50)
    print("ensemble_all num points", len(ensemble_all), "num of ensemble for each point", len(ensemble_all[0]))
    
    plt.savefig("test_ensemble_multi_peaks_data.png")
    
    # kde1d results from cpp backend
    kde1d=[0.104024,0.104784,0.105468,0.106075,0.106605,0.107059,0.107437,0.10774,0.107969,0.108125,0.108211,0.108228,0.108178,0.108065,0.10789,0.107657,0.107368,0.107028,0.106639,0.106206,0.105733,0.105222,0.104679,0.104108,0.103512,0.102896,0.102265,0.101623,0.100974,0.100322,0.0996714,0.0990265,0.0983911,0.0977692,0.0971643,0.0965799,0.0960195,0.0954861,0.0949826,0.0945117,0.0940757,0.0936768,0.0933169,0.0929975,0.0927199,0.0924851,0.0922938,0.0921463,0.0920427,0.0919827,0.0919658,0.0919912,0.0920576,0.0921637,0.0923078,0.0924878,0.0927016,0.0929466,0.0932201,0.0935193,0.093841,0.0941819,0.0945385,0.0949073,0.0952846,0.0956664,0.096049,0.0964284,0.0968006,0.0971616,0.0975074,0.097834,0.0981375,0.0984142,0.0986601,0.0988718,0.0990456,0.0991783,0.0992666,0.0993074]
    kde1d_x = []
    for i in range(80):
        # get kde value, the input is from 0 to 2.5
        kde1d_x.append(i*2.5/80.0)
    plt.clf()
    plt.plot(kde1d_x,kde1d)
    plt.savefig("test_ensemble_multi_peaks_data_kde_output.png")

    
    # for each ensemble member
    # create a list for it
    # for i in range(len(ensemble_all[0])):
    #     print("output ensemble with",i)
    #     ensemble_list_local=[]
    #     output_filename = filename[:-4] + "_" + str(i) + ".vtk"
    #     for j in range(len(ensemble_all)):
    #         ensemble_list_local.append(ensemble_all[j][i])

    #     writeVTKDataFromArray(xdim,ydim,zdim,output_filename,ensemble_list_local)
        

