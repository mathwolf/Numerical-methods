# Compute the coefficients
# for m data points (x[i], y[i]) with 0 <= i <= m-1
# where m = 2 ** p and
# x[i] = -PI + i * PI / m

import math             # basic mathematical functions
import cmath            # functions for complex numbers
import numpy as np      # matrix and vector functions

PI = math.pi

def FastFourierTransform(m, p, y_vector):
    # Input parameters
    # m: number of data points
    # p: m = 2 ** p; the data points will be some power of 2
    # y_vector: a vector of length 2m containing the y-values of
    #       the data points

    M = m
    q = p
    zeta = cmath.exp(complex(0,PI) / m)

    # Create an array that holds 2*m complex numbers where
    # we will put our coefficients.  This array will initally
    # be populated with our y_values.
    c_vector = np.zeros((2*m), dtype=complex)
    for j in range(2*m):
        c_vector[j] = y_vector[j]

    print(str(c_vector))


    # Comment to describe xi_vector
    xi_vector = np.zeros((2*M + 1), dtype=complex)
    for j in range(1, M+1):
        xi_vector[j] = zeta ** j
        xi_vector[j + M] = -1 * xi_vector[j]

    K = 0
    xi_vector[0] = 1

    k_vector = np.zeros((p+1))  # Holds the decomposition of K

    # Step 5
    for L in range(1, p+2):
        # Step 6
        if K < 2*m - 1:
            # Step 7
            for j in range(1, M+1):
                # Step 8
                K_temp = K
                # Decompose K
                for n in range(p, -1, -1):
                    if K_temp >= 2**n:
                        k_vector[n] = 1
                        K_temp = K_temp - 2**n
                    else:
                        k_vector[n] = 0
                K1 = 0
                for n in range(p - q + 1):
                    K1 = K1 + k_vector[n] * 2**n
                K2 = 0
                for n in range(q, p+1):
                    K2 = K2 + k_vector[n] * 2**n
                # Step 9
                eta = c_vector[K+M] * xi_vector[K2]
                c_vector[K+M] = c_vector[K] - eta
                c_vector[K] = c_vector[K] + eta
                # Step 10
                K = K + 1
            # Step 11
            K = K + M
        # Step 12
        K = 0
        M = M / 2
        q = q - 1

    # Step 13
    if K < 2*m - 1:
        # Step 14
        K_temp = K
        # Decompose K
        for n in range(p, -1, -1):
            if K_temp >= 2**n:
                k_vector[n] = 1
                K_temp = K_temp - 2**n
            else:
                k_vector[n] = 0
        j = 0
        for n in range(p+1):
            j = j + k_vector[n] * 2**(p-n)
        # Step 15
        if j > K:
            c_temp = c_vector[j]
            c_vector[j] = c_vector[K]
            c_vector[K] = c_temp
        # Step 16
        K = K + 1

    # Step 17
    a_vector = np.zeros((m+1))
    b_vector = np.zeros((m-1))
    a_vector[0] = c_vector[0].real/m
    a_vector[m] = (cmath.exp(complex(0, -1*PI*m)) * c_vector[m] / m).real
    # Step 18
    for j in range(1,m):
        a_vector[j] = (cmath.exp(complex(0, -1*PI*j)) * c_vector[j] / m).real
        b_vector[j-1] = (cmath.exp(complex(0, -1*PI*j)) * c_vector[j] / m).imag

    print(str(c_vector))
    print(str(a_vector))
    print(str(b_vector))
   
    
    

        

        
                 

        


def main():
    m = 4
    p = 2
    y_vector = np.array([1,2,3,4,5,6,7,8])

    FastFourierTransform(m, p, y_vector)


main()
