# load the wind data set

# build the histogram for each point

# compute the capability of different kernel
from matplotlib import pyplot
import numpy as np 
import os
from numpy.random import default_rng 
import random 
import time
import matplotlib.pyplot as plt
# import multiprocessing as mp
import glob
import pandas as pd
import statistics
from numpy import exp
from sklearn.neighbors import KernelDensity
import seaborn as sns

def load_txt_data(dir):
    path = dir+"Lead_33_*.txt"
    print("data file path: ", path)
    data_ensembles=[]
    for filename in glob.glob(path):
        data = np.loadtxt(filename)
        data_ensembles.append(data)   
    return np.array(data_ensembles)

if __name__ == '__main__':
    # # print("Number of processors: ", mp.cpu_count())
    # # pool = mp.Pool(mp.cpu_count())
    # # data_dir = "../datasets/wind_pressure_200/Lead_33.npy"
    # isovalue = 0.2

    # #data = np.load(data_dir)
    # data_dir = "../dataset/txt_files/wind_pressure_200/"
    # data = load_txt_data(data_dir)
    # xdim=121
    # ydim=240
    # new_data = []

    # # each d is an ensemble snapshot
    # print("ensemble version:",len(data))

    # for d in data:
    #     new_data.append(d.flatten())
    
    # # each row of new_data is flatten ensemble for different point
    # # each colum of new_data is 15 ensemble values for same point
    # print(len(new_data))
    # print(len(new_data[0]))
    
    # data_matrix = np.array(new_data)
    # res_list = []
    # res_max=0
    # res_max_id=0
    # for i in range(xdim*ydim):
    #     colum = data_matrix[:,i]
    #     #pyplot.hist(colum0, bins=8, density=False)
    #     #pyplot.savefig('colm0_hist.png')
    #     res = statistics.variance(colum)
    #     res_list.append(res)
    #     if res>res_max:
    #         res_max=res
    #         res_max_id=i

    # np_array_res_list=np.asarray(res_list)
    # reshaped_array = np_array_res_list.reshape(xdim,ydim)
    # plt.imshow(reshaped_array, interpolation='nearest')
    # plt.savefig('reshaped_array.png')

    # print("res_max_id",res_max_id)
    # max_var_colm=data_matrix[:,res_max_id]
    # print(max_var_colm)
    
    var_colm=np.array([0.83298728,0.50671819,0.59565333,0.46039456,0.65027173,0.77037451
,0.5241643,0.49069771,0.3680336,0.84358042,0.82585091,0.6840288,0.60165993,0.71032339,0.61138151])

    sample_min=min(var_colm)
    sample_max=max(var_colm)
    print(sample_min,sample_min)

    plt.clf()
    sns.kdeplot(np.array(var_colm), bw_method=0.3)
    plt.savefig('max_var_colm_sns_kdeplot.png')
    
    #pyplot.clf()
    #there is a bug when the density is set as true
    #https://stackoverflow.com/questions/55555466/matplotlib-hist-function-argument-density-not-working
    var_colm_bin=pyplot.hist(var_colm, bins=6, density=False, range=[sample_min,sample_max])
    
    print("var_colm_bin")
    print(var_colm_bin[0])
    var_colm_bin_density=var_colm_bin[0]/15.0
    print("var_colm_bin_density")
    print(var_colm_bin_density)
    
    print("sample_range for bar")
    sample_range=np.linspace(sample_min,sample_max,6)
    print(sample_range)
    
    pyplot.clf()
    
    pyplot.bar(sample_range,var_colm_bin[0], width=0.1)

    model = KernelDensity(bandwidth=0.2, kernel='gaussian')
    sample = var_colm.reshape((len(var_colm), 1))
    model.fit(sample)

    x_values=np.linspace(0,sample_max+0.5,10)
    x_values = x_values.reshape((len(x_values), 1))
    probabilities = model.score_samples(x_values)
    probabilities = exp(probabilities)
    
    print("x_values",x_values)
    print("probabilities",probabilities)

    pyplot.plot(x_values, probabilities, color="red")
    pyplot.savefig('colm_max_hist_density.png')

    # pyplot.clf()
    # # # create a model and fit the observed data
    # # # update, the input here should be the histogramed value
    # # # instead of the observation
    # model = KernelDensity(bandwidth=0.1, kernel='epanechnikov')
    # sample = max_var_colm.reshape((len(max_var_colm), 1))
    # model.fit(sample)
    # print("sample")
    # print(sample)

    # values=np.linspace(0.0,1.0,20)
    # values = values.reshape((len(values), 1))
 
    # probabilities = model.score_samples(values)
    # # # # the value returned by score samples are log likelihood
    # # # # use exp to get the original probability density value
    # probabilities = exp(probabilities)

    # pyplot.hist(sample, bins=10, density=True)
    # #pyplot.plot(values[:], probabilities)
    # pyplot.savefig('clm_max_kde.png')