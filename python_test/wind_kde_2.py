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
from scipy.stats import norm
from scipy import stats, optimize

# there kernel are standard unit
def kernel(k: str):
    """Kernel Functions.
    Ref: https://en.wikipedia.org/wiki/Kernel_(statistics)

    Args:
        k (str): Kernel name. Can be one of ['gaussian', 'epanechnikov', 'cosine', 'linear'.]
    """

    if k not in ['gaussian', 'epanechnikov', 'cosine', 'linear']:
        raise ValueError('Unknown kernel.')

    def bounded(f):
        def _f(x):
            return f(x) if np.abs(x) <= 1 else 0
        return _f

    if k == 'gaussian':
        return lambda u: 1 / np.sqrt(2 * np.pi) * np.exp(-1 / 2 * u * u)
    elif k == 'epanechnikov':
        return bounded(lambda u: (3 / 4 * (1 - u * u)))
    elif k =='cosine':
        return bounded(lambda u: np.pi / 4 * np.cos(np.pi / 2 * u))
    elif k == 'linear':
        return bounded(lambda u: 1 - np.abs(u))

# three ways to find the suitable h
def bw_scott(data: np.ndarray):
    std_dev = np.std(data, axis=0, ddof=1)
    n = len(data)
    return 3.49 * std_dev * n ** (-0.333)

def bw_silverman(data: np.ndarray):
    def _select_sigma(x):
        normalizer = 1.349
        iqr = (stats.scoreatpercentile(x, 75) - stats.scoreatpercentile(x, 25)) / normalizer
        std_dev = np.std(x, axis=0, ddof=1)
        return np.minimum(std_dev, iqr) if iqr > 0 else std_dev
    sigma = _select_sigma(data)
    n = len(data)
    return 0.9 * sigma * n ** (-0.2)

def bw_mlcv(data: np.ndarray, k):
    """
    Ref: https://rdrr.io/cran/kedd/src/R/MLCV.R
    """
    n = len(data)
    x = np.linspace(np.min(data), np.max(data), n)
    def mlcv(h):
        fj = np.zeros(n)
        for j in range(n):
            for i in range(n):
                if i == j: continue
                fj[j] += k((x[j] - data[i]) / h)
            fj[j] /= (n - 1) * h
        return -np.mean(np.log(fj[fj > 0]))
    h = optimize.minimize(mlcv, 1)
    if np.abs(h.x[0]) > 10:
        return bw_scott(data)
    return h.x[0]

def kde(data, k=None, h=None, x=None):
    """Kernel Density Estimation.

    Args:
        data (np.ndarray): Data.
        k (function): Kernel function.
        h (float): Bandwidth.
        x (np.ndarray, optional): Grid. Defaults to None.

    Returns:
        np.ndarray: Kernel density estimation.
    """
    # line space is between the min and max of the input data
    if x is None:
        x = np.linspace(np.min(data), np.max(data), 100)
    if h is None:
        h = bw_silverman(data)
    if k is None:
        k = kernel('gaussian')
    n = len(data)
    kde = np.zeros_like(x)
    # using the kde function to accumulate the positions
    # at each obvervation points
    for j in range(len(x)):
        for i in range(n):
            kde[j] += k((x[j] - data[i]) / h)
        kde[j] /= n * h
    return kde


def load_txt_data(dir):
    path = dir+"Lead_33_*.txt"
    print("data file path: ", path)
    data_ensembles=[]
    for filename in glob.glob(path):
        data = np.loadtxt(filename)
        data_ensembles.append(data)   
    return np.array(data_ensembles)

if __name__ == '__main__':
    var_colm=np.array([0.83298728,0.50671819,0.59565333,0.46039456,0.65027173,0.77037451
,0.5241643,0.49069771,0.3680336,0.84358042,0.82585091,0.6840288,0.60165993,0.71032339,0.61138151])
    x_points = np.linspace(np.min(var_colm), np.max(var_colm), 100)
    pyplot.hist(var_colm, bins=8, density=True)
    

    print("dedicated h")
    print(bw_scott(var_colm))
    print(bw_silverman(var_colm))
    print(bw_mlcv(var_colm,kernel('gaussian')))



    # pay attention, this kde_result is actually the histogram 
    kde_result=kde(var_colm,h=0.05)
    # this is the histogram results
    pyplot.plot(x_points,kde_result)
    pyplot.savefig('wind_kde_2.png')


    # compute the mse of the kde?

