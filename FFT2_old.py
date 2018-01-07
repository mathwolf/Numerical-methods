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

def main():
    N = 8
    y_vector = np.array([1,2,3,4,5,6,7,8], dtype=complex)

    y_vector = FastFourierTransform(y_vector, N)
    print(str(y_vector))


main()
