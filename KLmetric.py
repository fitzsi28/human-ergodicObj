import numpy as np
import scipy.integrate as scint
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.stats import multivariate_normal


L1 = 1.
L2 =1.
Sigma = 0.01
NBINS = 300
fontsz=14
title_font = {'fontname':'Liberation Sans', 'size':'20', 'color':'black', 'weight':'normal'}#,'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font = {'fontname':'Liberation Sans', 'size':'18'}#'21'}

xj, yj = np.mgrid[-L1/2:L1/2:NBINS*1j, -L2/2:L2/2:NBINS*1j]
xi, yi = np.mgrid[-0.6:0.6:NBINS*1j, -0.6:0.6:NBINS*1j]
domain = np.transpose(np.vstack([np.linspace(-L1/2,L1/2,NBINS),np.linspace(-L2/2,L2/2,NBINS)]))
ref = 0.3*multivariate_normal.pdf(np.dstack((xj,yj)), mean = [-0.2,0.1], cov=[[0.001,0.0],[0.0,0.001]]) \
+0.3*multivariate_normal.pdf(np.dstack((xj,yj)), mean = [0.3,-0.1], cov=[[0.01,0.0],[0.0,0.01]]) \
+0.4*multivariate_normal.pdf(np.dstack((xj,yj)), mean = [0.35,0.3], cov=[[0.001,0.0],[0.0,0.001]])


data=genfromtxt('/home/kt-fitz/human-ergodicObj/DIdkltest.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
tlist = data[0:-1,0]
x1 = data[0:-1,1]#+(L1/2.)
x2 = data[0:-1,3]#+(L2/2.)
X = np.stack([x1,x2],axis=1)

samps=genfromtxt('/home/kt-fitz/human-ergodicObj/Domain_samples.csv',delimiter=",",dtype=float)

Nt = np.shape(X)[0]
x_approx=1./Nt*multivariate_normal.pdf(np.dstack((xi,yi)), mean = [x1[0],x2[0]], cov=[[Sigma,0.0],[0.0,Sigma]])
for i in range(1,Nt):
    x_approx = x_approx +1./Nt*multivariate_normal.pdf(np.dstack((xi,yi)), mean = [x1[i],x2[i]], cov=[[Sigma,0.0],[0.0,Sigma]])

plt.figure()
#plt.plot(tlist,data[0:-1,5])
#plt.plot(tlist,data[0:-1,6])
plt.plot(samps[0],samps[1],'k.')
plt.ylim(-0.5,0.5)
plt.xlim(-0.5,0.5)


plt.figure()
plt.pcolormesh(xj, yj, ref.reshape(xj.shape))#,norm=colors.Normalize(vmin=0,vmax=10.0))
plt.plot(x1,x2,'k')
plt.title("Reference Distribution", **title_font)
plt.xlabel ( r"$x$",**axis_font)
plt.ylabel ( r"$y$",**axis_font)
plt.margins(0)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=fontsz)
cbar.ax.set_ylabel('Density',fontsize=16)

plt.figure()
plt.pcolormesh(xi, yi, x_approx.reshape(xi.shape))#,norm=colors.Normalize(vmin=0,vmax=10.0))
plt.plot(x1,x2,'k')
plt.title("Gaussian Approximation of X(t)", **title_font)
plt.xlabel ( r"$x$",**axis_font)
plt.ylabel ( r"$y$",**axis_font)
plt.margins(0)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=fontsz)
cbar.ax.set_ylabel('Density',fontsize=16)

plt.figure()
plt.plot(tlist,data[0:-1,7])

plt.show()