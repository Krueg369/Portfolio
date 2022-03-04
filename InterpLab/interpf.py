# Functions Program
import math
import numpy as np
import matplotlib.pyplot as plt

def test1(x):
    y = np.sin(x)
    return y


def test2(p, exact):
    n = len(p)
    sum1 = 0
    sum2 = 0
    for i in p:
        sum1 = sum1 + i
    for i in exact:
        sum2 = sum2 + i

    err = sum2 / sum1 
    return err
def interpfx(x, xdata, fdata):
    p = 0
    n = len(xdata)
    for i in range(0, n):
        L = 1
        for j in range(0, n):
            if j != i:
                numerator = x - xdata[j]
                denomenator = xdata[i] - xdata[j]
                L = (L * (numerator)/(denomenator))
        p = p + L*fdata[i]
    return p
