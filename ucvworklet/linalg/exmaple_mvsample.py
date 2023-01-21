import numpy as np 
import matplotlib.pyplot as plt
from numpy import linalg as LA

# Define dimension. 
d = 4

# Set mean vector. 
m = np.array([1, 2]).reshape(2, 1)

# Set covariance function. 
K_0 = np.array([[0.2, 0.6, 0.4, 0.8],
                [0.6, 1.8, 1.2, 2.4],
                [0.4, 1.2, 0.8, 1.6],
                [0.8, 2.4, 1.6, 3.2]])

print(np.linalg.eigvals(K_0))


L = np.array([[0, 0, -2.92872e-09, 0.447214],
                [-0, 0, -2.76672e-18, 1.34164],
                [0, -0, 2.92872e-10, 0.894427],
                [0, 0, 5.85744e-10, 1.78885]])

L_transpose = L.transpose()

print("checking results:")
print(np.matmul(L,L_transpose))

K_1 = np.array([[3,1],
                [1,3]])

K_2 = np.array([[1.0,3.0,7.0],
                [3.0,2.0,6.0],
                [7.0,6.0,5.0]])

K_3 = np.array([[1.0,3.0,7.0,8.0],
                [3.0,2.0,6.0,7.0],
                [7.0,6.0,5.0,6.0],
                [8.0,7.0,6.0,5.0]])

w, v = LA.eig(K_2)

# compute eigen values
print("eigen values")
print(w)
#print(w[0])
#print(w[1])

print("eigen vectors shape")
print(np.shape(v))
print(v)
#print(v[0])
#print(v[1])
#print(v[2])
#print(v[3])
# compute eigen vactors

# pay attention here, do not do the transpose
# the column v[:,i] is the eigenvector corresponding to the eigenvalue w[i]
U = v
print("U")
print(U)
sqrt_lbd = np.diag(np.sqrt(w))
print("sqrt_lbd")
print(sqrt_lbd)

A = np.matmul(U,sqrt_lbd)

print("A")
print(A)

print("A*A^t")
print(np.matmul(A, A.transpose()))