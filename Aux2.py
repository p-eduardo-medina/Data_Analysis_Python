import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import astropy as ap
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

def fact(n):  
    fa=1
    for i in range(1, n+1):
        fa *= i
    return fa
def combinatoria(n,k):
    C = 1
    C = fact(n)/((fact(k))*(fact(n-k)))
    return C

def BinomialCDF(n, p):
    sum = 0
    for i in range(0,n):
        sum += combinatoria(n,i)*(p**i)*((1-p)**(n-i))
    return round(sum,7)
print(BinomialCDF(5,0.4))