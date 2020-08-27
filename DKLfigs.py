import numpy as np
import scipy.integrate as scint
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.stats import multivariate_normal
import cv2

L1 = 1.
L2 =1.
Sigma = 0.01
NBINS = 300
fontsz=14
title_font = {'fontname':'Liberation Sans', 'size':'20', 'color':'black', 'weight':'normal'}#,'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font = {'fontname':'Liberation Sans', 'size':'18'}#'21'}

xj, yj = np.mgrid[-L1/2:L1/2:NBINS*1j, -L2/2:L2/2:NBINS*1j]
xi, yi = np.mgrid[-L1/2:L1/2:NBINS*1j, -L2/2:L2/2:NBINS*1j]#np.mgrid[-0.6:0.6:NBINS*1j, -0.6:0.6:NBINS*1j]
domain = np.transpose(np.vstack([np.linspace(-L1/2,L1/2,NBINS),np.linspace(-L2/2,L2/2,NBINS)]))

img = cv2.imread("banana.png",cv2.IMREAD_GRAYSCALE)

            
data=genfromtxt('/home/kt-fitz/data/cpp2/s10_c_set01-clean.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
tlist = data[0:-1,0]
x2 = (data[0:-1,-2])
x1 = -1*(data[0:-1,-3])

xx1 = x1[8050:8500]
xx2 = x2[8050:8500]
X = np.stack([xx1,xx2],axis=1)

Nt = np.shape(X)[0]
x_approx=1./Nt*multivariate_normal.pdf(np.dstack((xi,yi)), mean = [xx1[0],xx2[0]], cov=[[Sigma,0.0],[0.0,Sigma]])
for i in range(1,Nt):
    x_approx = x_approx +1./Nt*multivariate_normal.pdf(np.dstack((xi,yi)), mean = [xx1[i],xx2[i]], cov=[[Sigma,0.0],[0.0,Sigma]])

samps=genfromtxt('/home/kt-fitz/human-ergodicObj/apple_samples.csv',delimiter=",",dtype=float)
Nsamp = np.shape(samps)[1]
phi_approx=samps[2,0]*multivariate_normal.pdf(np.dstack((xi,yi)), mean = [samps[0,0],samps[1,0]], cov=[[Sigma/10,0.0],[0.0,Sigma/10]])
for i in range(1,Nsamp):
    phi_approx = phi_approx +samps[2,i]*multivariate_normal.pdf(np.dstack((xi,yi)), mean = [samps[0,i],samps[1,i]], cov=[[Sigma/10,0.0],[0.0,Sigma/10]])

print(np.count_nonzero(samps[2]))
   

fig = plt.figure()
ax = fig.add_subplot(111)
plt.ylim((0.5,-0.5))
plt.xlim((-0.5,0.5))
plt.pcolormesh(xi, yi, x_approx.reshape(xi.shape),norm=colors.Normalize(vmin=0,vmax=6.0))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=fontsz)
cbar.ax.set_ylabel('Density',fontsize=16)
#ax.set_aspect('equal')
plt.margins(0)
plt.xlabel ( r"$x\ (meters)$",**axis_font)
plt.ylabel ( r"$y\ (meters)$",**axis_font)
plt.title("Baseline Drawing", **title_font)
"""
#s10_banana_c
plt.plot(x1[23050:23500],x2[23050:23500],'ko',markersize=2)
plt.plot(x1[38050:38500],x2[38050:38500],'ko',markersize=2)
plt.plot(x1[18050:18500],x2[18050:18500],'ko',markersize=2)
"""
"""
#s10_apple_v
#plt.plot(x1[50:550],x2[50:550],'ko',markersize=2)
#s10_apple_p
plt.plot(x1[1000:1500],x2[1000:1500],'ko',markersize=2)
plt.plot(x1[3050:3500],x2[3050:3500],'ko',markersize=2)
plt.plot(x1[6050:6500],x2[6050:6500],'ko',markersize=2)
"""
#s10_apple_c
#plt.plot(x1[7050:7500],x2[7050:7500],'ko',markersize=2)
plt.plot(x1[8050:8500],x2[8050:8500],'ko',markersize=2)
#plt.plot(x1[22050:22500],x2[22050:22500],'ko',markersize=2)

plt.show()