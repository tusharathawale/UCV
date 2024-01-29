import math
import numpy as np
from numpy import linalg as LA


cov = np.array([[207399.6875000,-8085.0615234, -664.2767944, 76.0896912, 19.5749016, 37946.8789062, 0.0000000, 0.0000000],
                [-8085.0615234,368099.3750000, -997.2315674, -1403.7443848, -361.1280212 ,-20408.1562500, 0.0000000, 0.0000000],
                [-664.2767944,-997.2315674, 495.0614929, -10.6393852, -2.7371032 ,2588.3798828 ,0.0000000, 0.0000000],
                [76.0896912, -1403.7443848, -10.6393852, 907.5153809, 233.4692078, -680.5007324, 0.0000000, 0.0000000],
                [19.5749016, -361.1280212, -2.7371032, 233.4692078, 60.0625000, -175.0664978, 0.0000000 ,0.0000000],
                [37946.8789062, -20408.1562500, 2588.3798828, -680.5007324 ,-175.0664978, 346069.5000000, 0.0000000, 0.0000000],
                [0.0000000, 0.0000000, 0.0000000, 0.0000000 ,0.0000000 ,0.0000000, 0.0000000, 0.0000000],
                [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000 ,0.0000000, 0.0000000, 0.0000000]])


w, v = LA.eig(cov)

idx = w.argsort()[::-1]   
w = w[idx]
v = v[:,idx]

print(w)
print(v)

meanArray=[235.109375, 635.953125, 2.781250, 3.765625, 0.968750 ,501.890625 ,0.000000, 0.000000]

isovalue=900

numSamples=1000

indexImpEigenvalue = 0
threshold2 = 0.01

useall=True

        
# search important eigenvalue from end. Not checked for w[3] as it would have use independent model if
# w[3] > 0.2*w[0]
# Two threshold w7
if (w[7] > threshold2*w[0]):
    indexImpEigenvalue = 7
        
elif (w[6] > threshold2*w[0]):
    indexImpEigenvalue = 6
            
elif (w[5] > threshold2*w[0]):
    indexImpEigenvalue = 5
            
elif (w[4] > threshold2*w[0]):
    indexImpEigenvalue = 4
            
elif (w[3] > threshold2*w[0]):
    indexImpEigenvalue = 3
        
elif (w[2] > threshold2*w[0]):
    indexImpEigenvalue = 2
        
elif (w[1] > threshold2*w[0]):
    indexImpEigenvalue = 1

if useall:
    indexImpEigenvalue=6

print("indexImpEigenvalue", indexImpEigenvalue) 
        
# Transformed isovalue
transformedIso1 = (isovalue-meanArray[0]) 
transformedIso2 = (isovalue-meanArray[1])    
transformedIso3 = (isovalue-meanArray[2]) 
transformedIso4 = (isovalue-meanArray[3])
transformedIso5 = (isovalue-meanArray[4]) 
transformedIso6 = (isovalue-meanArray[5])    
transformedIso7 = (isovalue-meanArray[6]) 
transformedIso8 = (isovalue-meanArray[7]) 
    
xTest = np.zeros((numSamples))
yTest = np.zeros((numSamples))
zTest = np.zeros((numSamples))
wTest = np.zeros((numSamples))
aTest = np.zeros((numSamples))
bTest = np.zeros((numSamples))
cTest = np.zeros((numSamples))
dTest = np.zeros((numSamples))
    
for i in range(indexImpEigenvalue+1):

    temp = np.random.normal(0, math.sqrt(w[i]), numSamples).T
        
        
    xTest = xTest + v[0][i]*temp
    yTest = yTest + v[1][i]*temp
    zTest = zTest + v[2][i]*temp
    wTest = wTest + v[3][i]*temp
    aTest = aTest + v[4][i]*temp
    bTest = bTest + v[5][i]*temp
    cTest = cTest + v[6][i]*temp
    dTest = dTest + v[7][i]*temp

            
numNegativeNonCrossing = 0
for i in range(numSamples):
    if ((xTest[i]<=transformedIso1) and (yTest[i]<=transformedIso2) and (zTest[i]<=transformedIso3) and (wTest[i]<=transformedIso4) and (aTest[i]<=transformedIso5) and (bTest[i]<=transformedIso6) and (cTest[i]<=transformedIso7) and (dTest[i]<=transformedIso8)):
        numNegativeNonCrossing = numNegativeNonCrossing + 1
        
numPositiveNonCrossing = 0
for i in range(numSamples):
    if ((xTest[i]>=transformedIso1) and (yTest[i]>=transformedIso2) and (zTest[i]>=transformedIso3) and (wTest[i]>=transformedIso4) and (aTest[i]>=transformedIso5) and (bTest[i]>=transformedIso6) and (cTest[i]>=transformedIso7) and (dTest[i]>=transformedIso8)):
        numPositiveNonCrossing = numPositiveNonCrossing + 1
    
crossingProb = 1.0 - (numNegativeNonCrossing/numSamples) - (numPositiveNonCrossing/numSamples)

print(crossingProb)