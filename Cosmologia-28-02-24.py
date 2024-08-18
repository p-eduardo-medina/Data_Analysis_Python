

#Importaciones
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
#Variables matemáticas
dx = 0.01
#Parametros cosmológicos
Om0 = 0.315
Ol0 = 0.685
Og0 = 9.25e-5
#Constantes
G = 6.67e-11
c = 3e8
eV = 1.6e-19
ec = 4870e6*eV
Gyr = 365*24*60*60*1e9
pc = 3.0857e13
#Funciones
def integral(x,y):
    a = 0
    for i in range(1,len(x)-1):
        a += (x[i+1]-x[i])*y[i]
    return a

def Re(a,w,Om0):
    Ret = ec*Om0*a**(-1-3*w)
    return Ret

def epsilon(a):
    return Re(a,0,Om0) + Re(a,1/3,Og0) + Re(a,-1/3,Ol0)

def ap2(a): 
    return 8*math.pi*G/(3*c**2)*epsilon(a)   

def H_a(a):
    return math.sqrt( ap2(a)/(a**2))

def f_a(a):
    return 1/(math.sqrt(ap2(a)))

t0 = 0.1
tf = 1.5
N = int((tf-t0)/dx)
a = np.linspace(t0, tf, num=N)
 
F = []
for i in range(len(a)):
    F.append(f_a(a[i]))

dt = []
for i in range(len(a)):
    Auxa = a[0:i]
    Auxf = F[0:i]
    c = integral(Auxa,Auxf)/Gyr
    dt.append(c)

plt.plot(dt,a)
plt.grid()
plt.xlabel('t-$t_0$(Gyr)')
plt.ylabel('a(1)')
plt.title('Factor de Escala')
plt.show()

Ha=[]
for i in range(len(a)):
    Ha.append(H_a(a[i]))

plt.figure
plt.plot(dt,Ha)
plt.grid()
plt.xlabel('t-$t_0$(Gyr)')
plt.ylabel('H(km/s*pc)')
plt.title('Factor de Hubble')
plt.show()
    
           





