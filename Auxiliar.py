import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import astropy as ap
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

D = np.loadtxt('perlmutter.txt', comments='#')
z = D[:,0]
meff = D[:,1]
m_err = 1

from scipy.optimize import minimize

def model(params, z):
    Olambda , m , Ho = params       
    #print(f'{Olambda}, {m}, {Ho}')
    cosmo = FlatLambdaCDM(H0=Ho* u.km / u.s / u.Mpc, Tcmb0=2.725* u.K, Om0=1-Olambda)
    m_effModel = []     
    for i in range(len(z)):        
        dl = cosmo.luminosity_distance(z[i])
        m_effModel.append(m + 5*math.log10(dl.value*1e5))        
    return m_effModel

def chi_square(params, z_obs, m_obs, m_err):                        
    m_model = model(params, z_obs)
    chi_square = np.sum(((m_obs - m_model) / m_err) ** 2)
    return chi_square

guess = [0.99, -15.0, 65.0]
result = minimize(chi_square, guess, args = (z, meff, m_err),bounds= ((0, 1),(None,None), (0, None)))
Olambda = result.x
print(f"EL valor de Omega_Lambda es de {Olambda}")