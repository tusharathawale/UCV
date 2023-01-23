import numpy as np

# this code is from https://www.andreinc.net/2021/01/25/computing-eigenvalues-and-eigenvectors-using-qr-decomposition
# we use this for testing the correctness of c code

# A is a square random matrix of size n
n = 4
#A = np.random.rand(n, n)
'''
A = np.array([[0.00001,  0.00001,   0.00001,   0.00001],
                [0.00001,   20.000,   60.000,   40.000],
                [0.00001,   60.000 , 180.000 , 120.000],
                [0.00001,   40.000,  120.000,   80.000]])
'''
A = np.array([[0,  0,   0,   0],
                [0,   20.000,   60.000,   40.000],
                [0,   60.000 , 180.000 , 120.000],
                [0,   40.000,  120.000,   80.000]])

print("A=")
print((A))

np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

def eigen_qr_simple(A, iterations=5):
    Ak = np.copy(A)
    n = A.shape[0]
    QQ = np.eye(n)
    for k in range(iterations):
        print(Ak)
        Q, R = np.linalg.qr(Ak)
        print("eigen q")
        print(Q)
        print("eigen r")
        print(R)

        Ak = R @ Q
        QQ = QQ @ Q
        # we "peek" into the structure of matrix A from time to time
        # to see how it looks
        if k%10000 == 0:
            print("A",k,"=")
            print((Ak))
            print("\n")
    return Ak, QQ

# We call the function    
eigen_qr_simple(A)

# We compare our results with the official numpy algorithm
print(np.linalg.eigvals(A))

# compute the eigen vectors based on eigen value