import numpy as np
import scipy.integrate as scint
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.stats import multivariate_normal


L1 = 1.
L2 =1.
Sigma = np.eye(2)

data=genfromtxt('/home/kt-fitz/human-ergodicObj/DIdkltest.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
tlist = data[0:-1,0]
x1 = data[0:-1,1]#+(L1/2.)
x2 = data[0:-1,3]#+(L2/2.)
X = np.stack([x1,x2],axis=1)



def qfunc(s):
    qtemp = 0.
    Nt = np.shape(X)[0]
    for i in range(0,Nt):
        qtemp = qtemp + np.exp(-0.5*np.transpose(s-X[i]).dot(np.linalg.inv(Sigma)).dot(s-X[i]))
    qtemp = qtemp/Nt
    return qtemp

def eta_calc(f):
    N1 = 50
    N2 = 50
    x0=0.
    y0 = 0.
    xf = L1
    yf = L2
    d1 = L1/N1
    d2 = L2/N2
    total = (d1*d2/4)*(f([x0,y0])+f([x0,yf])+f([xf,y0])+f([xf,yf]))
    for i in range(1,N1):
        total +=(d1*d2/2)*(f([x0+i*d1,y0])+f([x0+i*d1,yf]))
    for j in range(1,N2):
        total+=(d1*d2/2)*(f([x0,y0+j*d2])+f([xf,y0+j*d2]))
    for m in range(1,N2):
        for n in range(1,N1):
            total+=d1*d2*f([x0+m*d1,y0+n*d2]) 
    return 1/total

plt.figure()
plt.plot(tlist,x1)
plt.plot(tlist,x2)

NBINS = 300
xj, yj = np.mgrid[-L1/2:L1/2:NBINS*1j, -L2/2:L2/2:NBINS*1j]
domain = np.transpose(np.vstack([np.linspace(-L1/2,L1/2,NBINS),np.linspace(-L2/2,L2/2,NBINS)]))
ref = 0.3*multivariate_normal.pdf(np.dstack((xj,yj)), mean = [-0.2,0.1], cov=[[0.001,0.0],[0.0,0.001]])\
+0.3*multivariate_normal.pdf(np.dstack((xj,yj)), mean = [0.3,-0.1], cov=[[0.01,0.0],[0.0,0.01]])\
+0.4*multivariate_normal.pdf(np.dstack((xj,yj)), mean = [0.35,0.3], cov=[[0.001,0.0],[0.0,0.001]])
plt.figure()
plt.pcolormesh(xj, yj, ref.reshape(xj.shape))#,norm=colors.Normalize(vmin=0,vmax=1.0))
plt.plot(x1,x2,'k')
#plt.plot(tlist,data[0:-1,5])
#plt.plot(tlist,data[0:-1,6])

plt.figure()
plt.plot(tlist,data[0:-1,7])

plt.show()