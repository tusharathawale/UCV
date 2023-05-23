import numpy as np
from numpy import linalg as LA

A = np.array([
    [8.000000,7.000000,6.000000,5.000000],
    [0.125000,2.125000,6.250000,7.375000],
    [0.375000,-0.294118,5.588235,7.294118],
    [0.875000,-0.058824,0.021053,1.905263]])

gsl_inv_A = np.array([
[-0.314917,0.325967,0.016575,0.027624],
[0.325967,-0.530387,0.193370,-0.011050],
[0.016575,0.193370,-0.685083,0.524862],
[0.027624,-0.011050,0.524862,-0.458564]
])
resultInverse= np.linalg.inv(A)

print(resultInverse)

print(np.matmul(A,resultInverse))

print("test 2")
print(np.matmul(A,gsl_inv_A))