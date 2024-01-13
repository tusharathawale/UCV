import numpy as np
import math

array=[1,2,1,2,3,100,78,34.5,67.2,20,9,3,101]
num_bins = 8
hist1, bin_edges1 = np.histogram(array, bins=num_bins)
density=[]
for i in range(num_bins):
    density.append(0)

print(hist1)
print(bin_edges1)

print("density input",density)

# compute manually
minv = min(array)
maxv = max(array)
slot_size = (maxv-minv)/(1.0*num_bins)
print("slot_size",slot_size)
for v in array:
    #remember to minus minv 
    bin_index = math.floor((v-minv)/slot_size)
    # for the last one
    if(abs(v-maxv)<0.0000001 and bin_index==num_bins):
        bin_index=bin_index-1
    #print(v,bin_index)
    density[bin_index]+=1

print(density)


