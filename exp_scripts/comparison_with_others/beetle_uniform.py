from matplotlib import colors, cm
import matplotlib.pyplot as plt
import scipy.io
import itertools
#from skimage import io, filters, morphology, restoration, feature, transform
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D
#import seaborn as sns
from itertools import cycle
#import sklearn.cluster
#from sklearn.preprocessing import MinMaxScaler
#import skimage.segmentation as seg

import numpy as np
import math
import numpy.ma as ma


#import nglpy as ngl
#import topopy
import time
import os
from functools import partial
from itertools import cycle
#from pyevtk.hl import imageToVTK
import random

#there are some issues loading these libraries
#try to see if they are necessary

#import nglpy_cuda as ngl
#import pdir
#import utpy.utils
#import utpy.vis
#import flatpy

# Load ground truth
myarray = np.fromfile('stagbeetle832x832x494.dat', dtype='H')

print(len(myarray),flush=True)

newarray = myarray[3: len(myarray)]

ground_truth = newarray.reshape((494, 832, 832))

#newarray.astype('H').tofile('beetle832x832x484.raw')

#print(832*832*494)

# Preprare ensemble from the ground truth
ensemble = np.zeros((123,208,208,64),dtype='H')

time1 = time.perf_counter()
print("start timer1",flush=True)
# why it is 492 not 494 here (divisable by 4)
# the value of i is 0 4 8 12 ...
for i in range(0,492,4):
    for j in range(0,832,4):
        for k in range(0,832,4):
            # l, m, n is the index for the 4*4*4 data block
            # 4*4*4 block
            l = int(i/4) 
            m = int(j/4)
            n = int(k/4)
            
            count = 0
            # go through the 4*4*4 data block
            for n1 in range(4):
                for n2 in range(4):
                    for n3 in range(4):
                        # ensemble datasize size is 123*208*208*64
                        ensemble[l,m,n,count] = ground_truth[i+n1,j+n2,k+n3]
                        count += 1

# Compute the probabillity of values being less than the isovalue and greater than the isovalue at each grid vertex
def computePositiveNegativeProbabilities(ensemble, isovalue, confidence, distributionModel):
    h, w, d, numMembers = ensemble.shape
    minimumGrid = np.zeros((h,w,d))
    maximumGrid = np.zeros((h,w,d))
    meanGrid = np.zeros((h,w,d))
    mostProbableGrid = np.zeros((h,w,d))
    
    meanGaussian = np.zeros((h,w,d))
    stdGaussian = np.zeros((h,w,d))
    for i in range(h):
        for j in range(w):
            for k in range(d):
                # compute min max for ensemble data
                minimum = math.inf
                maximum = -math.inf
                t = np.array([])
                # numMembers contains multiple values for each position
                # of the ensemble data
                for idx in range(numMembers):
                    t = np.append(t, [ensemble[i,j,k,idx]])
                    if(ensemble[i,j,k,idx] < minimum):
                        minimum = ensemble[i,j,k,idx]
                    if(ensemble[i,j,k,idx] > maximum):
                        maximum = ensemble[i,j,k,idx]    
                tMean = np.mean(t,0)
                tVar = np.var(t,0)
                tStd = math.sqrt(tVar)
            
                #Uniform
                minimumGrid[i,j,k] = minimum
                maximumGrid[i,j,k] = maximum   
            
                #Gaussian
                meanGaussian[i,j,k] = tMean
                stdGaussian[i,j,k] = tStd
            
    # Compute probabilitis based of confidence value (between 0 and 100)
    for i in range(h):
        for j in range(w):
            for k in range(d):
            
                minimum=minimumGrid[i,j,k]
                maximum=maximumGrid[i,j,k]
            
                minimumGrid[i,j,k]=minimum
                maximumGrid[i,j,k]=maximum
            
                meanVal = (maximum + minimum)/2
            
                confidenceFraction = 0.01*confidence*(meanVal - minimum)
            
                newMinimum = meanVal - confidenceFraction
                newMaximum = meanVal + confidenceFraction
            
                # Revise minimum and maximum grids based on confidence information
                minimumGrid[i,j,k]=newMinimum
                maximumGrid[i,j,k]=newMaximum
                meanGrid[i,j,k] = (newMinimum + newMaximum)/2
    
                #plt.figure()
                #plt.imshow(minimumGrid)
                #plt.contour(minimumGrid, levels=[isovalue])
    
                #plt.figure()
                #plt.imshow(maximumGrid)
                #plt.contour(maximumGrid, levels=[isovalue])
    
    # Compute negative and positive probabilities
    # Pr(value <= isovalue)
    negativeProbabilitiesGrid = np.zeros((h,w,d))
    # Pr(value > isovalue)
    positiveProbabilitiesGrid = np.zeros((h,w,d))
    # h is 128 w is 208 d is 208
    for i in range(h):
        for j in range(w):
            for k in range(d):
                # fetch minimal maximal and mean value for each voxel of the ensemble data
                minimum = minimumGrid[i,j,k]
                maximum = maximumGrid[i,j,k]
                meanVal = meanGrid[i,j,k]
                if (isovalue <= minimum):
                    positiveProbabilitiesGrid[i,j,k] = 1.0
                    negativeProbabilitiesGrid[i,j,k] = 0.0
                elif (isovalue >= maximum):    
                    positiveProbabilitiesGrid[i,j,k] = 0.0
                    negativeProbabilitiesGrid[i,j,k] = 1.0
                else:
                    if distributionModel == 'uniform':
                        positiveProbabilitiesGrid[i,j,k] = (maximum - isovalue)/(maximum - minimum)
                        negativeProbabilitiesGrid[i,j,k] = 1.0 - positiveProbabilitiesGrid[i,j,k]
                    elif distributionModel == 'Gaussian':
                        tempMean = meanGaussian[i,j,k]
                        tempStd = stdGaussian[i,j,k]
                        negativeProbabilitiesGrid[i,j,k] = 0.5*(1 + math.erf((isovalue - tempMean)/(math.sqrt(2)*tempStd)))
                        positiveProbabilitiesGrid[i,j,k] = 1.0 - negativeProbabilitiesGrid[i,j,k]
                    else:
                        print('No valid model')
            
                # Compute most probable grid
                if(positiveProbabilitiesGrid[i,j,k] > negativeProbabilitiesGrid[i,j,k]):
                    mostProbableGrid[i,j,k] = (meanVal + maximum)/2
                else:
                    mostProbableGrid[i,j,k] = (minimum + meanVal)/2  
            
    return positiveProbabilitiesGrid, negativeProbabilitiesGrid, maximumGrid, minimumGrid, mostProbableGrid

def computeMSTopologyCasesProbabilities(positiveProbabilitiesGrid, negativeProbabilitiesGrid):
    h, w, d = positiveProbabilitiesGrid.shape
    MStopologyCasesProbabilities = np.zeros((h-1,w-1,d-1, 256))
    # Record number of cases per cell as well
    numCases = np.zeros((h-1,w-1,d-1))

    tempNumCases = 1
    for i in range(h-1):
        for j in range(w-1):
            for k in range(d-1):
            
                x1Positive = positiveProbabilitiesGrid[i,j,k]
                x1Negative = negativeProbabilitiesGrid[i,j,k]
                x1Sign = np.array([x1Negative, x1Positive])
                if(x1Positive > 0) and (x1Positive < 1):
                    tempNumCases = tempNumCases*2
            
                x2Positive = positiveProbabilitiesGrid[i+1,j,k]
                x2Negative = negativeProbabilitiesGrid[i+1,j,k]
                x2Sign = np.array([x2Negative, x2Positive])
                if(x2Positive > 0) and (x2Positive < 1):
                    tempNumCases = tempNumCases*2
            
                x3Positive = positiveProbabilitiesGrid[i+1,j+1,k]
                x3Negative = negativeProbabilitiesGrid[i+1,j+1,k]
                x3Sign = np.array([x3Negative, x3Positive])
                if(x3Positive > 0) and (x3Positive < 1):
                    tempNumCases = tempNumCases*2
            
                x4Positive = positiveProbabilitiesGrid[i,j+1,k]
                x4Negative = negativeProbabilitiesGrid[i,j+1,k]
                x4Sign = np.array([x4Negative, x4Positive])
                if(x4Positive > 0) and (x4Positive < 1):
                    tempNumCases = tempNumCases*2
                    
                x5Positive = positiveProbabilitiesGrid[i,j,k+1]
                x5Negative = negativeProbabilitiesGrid[i,j,k+1]
                x5Sign = np.array([x5Negative, x5Positive])
                if(x5Positive > 0) and (x5Positive < 1):
                    tempNumCases = tempNumCases*2
            
                x6Positive = positiveProbabilitiesGrid[i+1,j,k+1]
                x6Negative = negativeProbabilitiesGrid[i+1,j,k+1]
                x6Sign = np.array([x6Negative, x6Positive])
                if(x6Positive > 0) and (x6Positive < 1):
                    tempNumCases = tempNumCases*2
            
                x7Positive = positiveProbabilitiesGrid[i+1,j+1,k+1]
                x7Negative = negativeProbabilitiesGrid[i+1,j+1,k+1]
                x7Sign = np.array([x7Negative, x7Positive])
                if(x7Positive > 0) and (x7Positive < 1):
                    tempNumCases = tempNumCases*2
            
                x8Positive = positiveProbabilitiesGrid[i,j+1,k+1]
                x8Negative = negativeProbabilitiesGrid[i,j+1,k+1]
                x8Sign = np.array([x8Negative, x8Positive])
                if(x8Positive > 0) and (x8Positive < 1):
                    tempNumCases = tempNumCases*2
                
                numCases[i,j,k] = tempNumCases
                # Reset tempNumCases for next iteration
                tempNumCases = 1
                
                count = 0
                # the final probability based on 8 vertices
                for s1 in range(2):
                    for s2 in range(2):
                        for s3 in range(2):
                            for s4 in range(2):
                                for s5 in range(2):
                                    for s6 in range(2):
                                        for s7 in range(2):
                                            for s8 in range(2):
                                                MStopologyCasesProbabilities[i,j,k,count] = x1Sign[s1]*x2Sign[s2]*x3Sign[s3]*x4Sign[s4]*x5Sign[s5]*x6Sign[s6]*x7Sign[s7]*x8Sign[s8]
                                                count+=1
              
            
    # Compute entropy of probability distribution
    entropy = np.zeros((h-1,w-1,d-1))
        
    #for i in range(256):
        
    #    probs = MStopologyCasesProbabilities[:,:,:,i]
    #    res = ma.log2(probs);
    #    entropy += np.multiply(-probs,res)
        
        
    for i in range(h-1):
        for j in range(w-1):
            for k in range(d-1):
                entropyValue = 0
                res = 0
                for num in range(256):
                    probs = MStopologyCasesProbabilities[i,j,k,num]
                    if probs > 0:
                        res = ma.log2(probs)
                    entropyValue += np.multiply(-probs,res)    
                entropy[i,j,k] = entropyValue
        
    return MStopologyCasesProbabilities, entropy, numCases

def computeMSTopologyCasesProbabilities(positiveProbabilitiesGrid, negativeProbabilitiesGrid):
    h, w, d = positiveProbabilitiesGrid.shape
    MStopologyCasesProbabilities = np.zeros((h-1,w-1,d-1, 256))
    # Record number of cases per cell as well
    numCases = np.zeros((h-1,w-1,d-1))

    tempNumCases = 1
    for i in range(h-1):
        for j in range(w-1):
            for k in range(d-1):
            
                x1Positive = positiveProbabilitiesGrid[i,j,k]
                x1Negative = negativeProbabilitiesGrid[i,j,k]
                x1Sign = np.array([x1Negative, x1Positive])
                if(x1Positive > 0) and (x1Positive < 1):
                    tempNumCases = tempNumCases*2
            
                x2Positive = positiveProbabilitiesGrid[i+1,j,k]
                x2Negative = negativeProbabilitiesGrid[i+1,j,k]
                x2Sign = np.array([x2Negative, x2Positive])
                if(x2Positive > 0) and (x2Positive < 1):
                    tempNumCases = tempNumCases*2
            
                x3Positive = positiveProbabilitiesGrid[i+1,j+1,k]
                x3Negative = negativeProbabilitiesGrid[i+1,j+1,k]
                x3Sign = np.array([x3Negative, x3Positive])
                if(x3Positive > 0) and (x3Positive < 1):
                    tempNumCases = tempNumCases*2
            
                x4Positive = positiveProbabilitiesGrid[i,j+1,k]
                x4Negative = negativeProbabilitiesGrid[i,j+1,k]
                x4Sign = np.array([x4Negative, x4Positive])
                if(x4Positive > 0) and (x4Positive < 1):
                    tempNumCases = tempNumCases*2
                    
                x5Positive = positiveProbabilitiesGrid[i,j,k+1]
                x5Negative = negativeProbabilitiesGrid[i,j,k+1]
                x5Sign = np.array([x5Negative, x5Positive])
                if(x5Positive > 0) and (x5Positive < 1):
                    tempNumCases = tempNumCases*2
            
                x6Positive = positiveProbabilitiesGrid[i+1,j,k+1]
                x6Negative = negativeProbabilitiesGrid[i+1,j,k+1]
                x6Sign = np.array([x6Negative, x6Positive])
                if(x6Positive > 0) and (x6Positive < 1):
                    tempNumCases = tempNumCases*2
            
                x7Positive = positiveProbabilitiesGrid[i+1,j+1,k+1]
                x7Negative = negativeProbabilitiesGrid[i+1,j+1,k+1]
                x7Sign = np.array([x7Negative, x7Positive])
                if(x7Positive > 0) and (x7Positive < 1):
                    tempNumCases = tempNumCases*2
            
                x8Positive = positiveProbabilitiesGrid[i,j+1,k+1]
                x8Negative = negativeProbabilitiesGrid[i,j+1,k+1]
                x8Sign = np.array([x8Negative, x8Positive])
                if(x8Positive > 0) and (x8Positive < 1):
                    tempNumCases = tempNumCases*2
                
                numCases[i,j,k] = tempNumCases
                # Reset tempNumCases for next iteration
                tempNumCases = 1
                
                count = 0
                for s1 in range(2):
                    for s2 in range(2):
                        for s3 in range(2):
                            for s4 in range(2):
                                for s5 in range(2):
                                    for s6 in range(2):
                                        for s7 in range(2):
                                            for s8 in range(2):
                                                MStopologyCasesProbabilities[i,j,k,count] = x1Sign[s1]*x2Sign[s2]*x3Sign[s3]*x4Sign[s4]*x5Sign[s5]*x6Sign[s6]*x7Sign[s7]*x8Sign[s8]
                                                count+=1
              
            
    # Compute entropy of probability distribution
    entropy = np.zeros((h-1,w-1,d-1))
        
    #for i in range(256):
        
    #    probs = MStopologyCasesProbabilities[:,:,:,i]
    #    res = ma.log2(probs);
    #    entropy += np.multiply(-probs,res)
        
        
    for i in range(h-1):
        for j in range(w-1):
            for k in range(d-1):
                entropyValue = 0
                res = 0
                for num in range(256):
                    probs = MStopologyCasesProbabilities[i,j,k,num]
                    if probs > 0:
                        res = ma.log2(probs)
                    entropyValue += np.multiply(-probs,res)    
                entropy[i,j,k] = entropyValue
        
    return MStopologyCasesProbabilities, entropy, numCases

def probabilisticMarchingSquares(positiveProbabilitiesGrid, negativeProbabilitiesGrid):
    
    h, w, d = positiveProbabilitiesGrid.shape
    # Record number of cases per cell as well
    crossingProbabilities = np.zeros((h-1,w-1,d-1))
    for i in range(h-1):
        for j in range(w-1):
            for k in range(d-1):
                nonCrossingProbability1 = positiveProbabilitiesGrid[i,j,k]*positiveProbabilitiesGrid[i+1,j,k]*positiveProbabilitiesGrid[i,j+1,k]*positiveProbabilitiesGrid[i+1,j+1,k]*positiveProbabilitiesGrid[i,j,k+1]*positiveProbabilitiesGrid[i+1,j,k+1]*positiveProbabilitiesGrid[i,j+1,k+1]*positiveProbabilitiesGrid[i+1,j+1,k+1]
                nonCrossingProbability2 = negativeProbabilitiesGrid[i,j,k]*negativeProbabilitiesGrid[i+1,j,k]*negativeProbabilitiesGrid[i,j+1,k]*negativeProbabilitiesGrid[i+1,j+1,k]*negativeProbabilitiesGrid[i,j,k+1]*negativeProbabilitiesGrid[i+1,j,k+1]*negativeProbabilitiesGrid[i,j+1,k+1]*negativeProbabilitiesGrid[i+1,j+1,k+1]
                nonCrossingProbability = nonCrossingProbability1 + nonCrossingProbability2
                crossingProbability = 1.0 - nonCrossingProbability
                crossingProbabilities[i,j,k] = crossingProbability
    return crossingProbabilities

time2 = time.perf_counter()
print("start timer2",flush=True)
isovalue = 900

#[positiveProbabilitiesGrid, negativeProbabilitiesGrid,  maximumGrid, minimumGrid, mostProbableGrid] = computePositiveNegativeProbabilities(ensemble, isovalue, 100, 'Gaussian')
[positiveProbabilitiesGrid, negativeProbabilitiesGrid,  maximumGrid, minimumGrid, mostProbableGrid] = computePositiveNegativeProbabilities(ensemble, isovalue, 100, 'uniform')

time3 = time.perf_counter()
print("start timer3",flush=True)
#plt.figure()
#plt.imshow(positiveProbabilitiesGrid)
#plt.contour(positiveProbabilitiesGrid, levels=[0.5])

MStopologyCasesProbabilities,entropy,numCases = computeMSTopologyCasesProbabilities(positiveProbabilitiesGrid, negativeProbabilitiesGrid)
time4 = time.perf_counter()

print(f"extracting ensembles in {time2 - time1:0.3f} seconds",flush=True)
print(f"computing properties of hixel blocks in {time3 - time2:0.3f} seconds",flush=True)
print(f"computing uncertainty metrics in {time4 - time3:0.3f} seconds",flush=True)

## do not do this for performance test
'''
imageToVTK("./mostProbableGrid", pointData = {"mostProbableGrid" : mostProbableGrid})
imageToVTK("./entropy", pointData = {"entropy" : entropy})
imageToVTK("./numCases", pointData = {"numCases" : numCases})

crossingProbabilities = probabilisticMarchingSquares(positiveProbabilitiesGrid, negativeProbabilitiesGrid)
imageToVTK("./crossingProbabilities", pointData = {"crossingProbabilities" : crossingProbabilities})


# Test single instance for isovalue 900
instance = np.zeros((123,208,208))
for i in range(123):
    for j in range(208):
        for k in range(208):
            instance[i,j,k] = ensemble[i,j,k,0]

print(instance.shape)
print(instance.flags)
print(ground_truth.flags)
print(ensemble.flags)
imageToVTK("./instance", pointData = {"instance" : instance})
imageToVTK("./ground_truth", pointData = {"ground_truth" : ground_truth})
'''
