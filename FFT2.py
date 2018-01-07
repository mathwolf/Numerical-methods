# python 2.7.2

# Steven Wray
# Complex Variables 2
# Implementation of the Fast Fourier Transform

# Compute the coefficients for the trigonometric
# polynomial on m data points (x[i], y[i]) where
# 0 <= i <= m-1, m is a power of 2 and
# x[i] = -PI + i * PI / m

import math             # basic mathematical functions
import cmath            # functions for complex numbers
import numpy as np      # matrix and vector functions

PI = math.pi


def FastFourierTransform(y_vector, N):
    # y_vector holds the y values to be interpolated.
    # The x values must be evenly spaced between -PI and PI
    # N is the length of y_vector.
    # N must be a power of 2
    
    if N == 1:
        # If there is only one y value in the vector,
        # do nothing.
        return(y_vector)
    else:
        # Create two subvectors for terms with even and odd
        # indices in the original vector
        subvector1 = np.zeros((N/2), dtype=complex)
        subvector2 = np.zeros((N/2), dtype=complex)
        for i in range(N/2):
            subvector1[i] = y_vector[2*i]       # even indices
            subvector2[i] = y_vector[2*i + 1]   # odd indices
            
        # Run the FFT on each subvector
        subvector1 = FastFourierTransform(subvector1, N/2)
        subvector2 = FastFourierTransform(subvector2, N/2)

        # Reassemble the original vector from the transformed
        # subvectors.  The FFT of entries with even indices go
        # in positions 0 to N/2-1.  The FFT of the odd indices
        # go in positions N/2 to N-1.
        for i in range(N/2):
            y_vector[i] = subvector1[i]
            y_vector[i + N/2] = subvector2[i]
            
        # Recover the complex coefficients from the reassembled vector
        for i in range(N/2):
            temp = y_vector[i]
            y_vector[i] = temp + cmath.exp(complex(0, -1 * 2 * PI * i / N)) *\
                          y_vector[i + N/2]
            y_vector[i + N/2] = temp - cmath.exp(complex(0, -1 * 2 * PI * i / N)) *\
                          y_vector[i + N/2]

        # Return the complex coefficients
        return y_vector

def PrintRealCoefficients(c_vector, N):
    # c_vector of length N contains the complex output of the FFT
    # subroutine.  Extract the real-valued coefficients a_k and b_k
    # for the trig polynomial from the first N/2 + 1 entries of the
    # vector.  Print the results.

    a_vector = np.zeros((N/2 + 1))
    b_vector = np.zeros((N/2 + 1))
    for i in range(N/2 + 1):
        a_vector[i] = ((-1)**i / (N/2.0)) * c_vector[i].real
        b_vector[i] = -1 * ((-1)**i / (N/2.0)) * c_vector[i].imag

    print("Cosine coefficients")
    print("Power\tCoefficient")
    for i in range(N/2 + 1):
        print(str(i) + "\t" + str(a_vector[i]))
    print("Sine coefficients")
    print("Power\tCoefficient")
    for i in range(1, N/2):
        print(str(i) + "\t" + str(b_vector[i]))


def main():
    while(True):
        print("Choose a sample problem for FFT.")
        print("\t1. Example 2")
        print("\t2. Problem 3b")
        print("\t3. Problem 3d")
        print("\t4. Problem 6")
        print("\t5. Exit")
        s = input("Enter a number: ")
        example = int(s)

        if example == 1:        # Burden & Faires 8.6 ex 2
            y_vector = np.array([0.0, 0.549761275, 1.11909646, \
                                 1.537853226, 1.557407725, 1.069103226, \
                                 0.36909646, -0.106488725], dtype=complex)
            y_vector = FastFourierTransform(y_vector, 8)
            PrintRealCoefficients(y_vector, 8)

        if example == 2:        # Burden & Faires 8.6 prob 3b
            y_vector = np.array([3.141592654, 2.35619449, 1.570796327, \
                                 0.785398163, 0.0, 0.785398163, \
                                 1.570796327, 2.35619449], dtype=complex)
            y_vector = FastFourierTransform(y_vector, 8)
            PrintRealCoefficients(y_vector, 8)

        if example == 3:        # Burden & Faires 8.6 prob 3d
            y_vector = np.array([2.879043276, -1.659010436, 1.430528853, \
                                 -0.231289529, 0.540302306, -0.638149666, \
                                 -0.756029015, -2.786850562], dtype=complex)
            y_vector = FastFourierTransform(y_vector, 8)
            PrintRealCoefficients(y_vector, 8)

        if example == 4:        # Burden & Faires 8.6 prob 6
            y_vector = np.array([-9.869604401, -8.507779734, -6.981217961, \
                                 -5.417424486, -3.925611112, -2.591696361, \
                                 -1.475364878, -0.609228939, 0, \
                                 0.368545901, 0.531131356, 0.535474455, \
                                 0.436179012, 0.288501896, 0.142473836, \
                                 0.037812534, 0, 0.037812354, \
                                 0.142473836, 0.288501896, 0.436179012, \
                                 0.535474455, 0.531131356, 0.368545901, \
                                 0, -0.609228939, -1.475364878, \
                                 -2.591696361, -3.925611112, -5.417424486, \
                                 -6.981217961, -8.507779734], dtype=complex)
            y_vector = FastFourierTransform(y_vector, 32)
            PrintRealCoefficients(y_vector, 32)


        else:
            break               # Leave the program
            



main()
