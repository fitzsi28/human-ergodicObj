import numpy as np
import scipy.integrate as scint
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.interpolate import interp1d
from scipy.stats import multivariate_normal
from scipy.stats import gaussian_kde

L1 = 20.0
L2 = 20.0
Knum = 5-1
NBINS=300
fontsz=14
title_font = {'fontname':'Arial', 'size':'20', 'color':'black', 'weight':'normal'}#,'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font = {'fontname':'Arial', 'size':'18'}#'21'}

def fourierRecon(xt,yt,cks,hk):
    x=xt+L1/2
    y=yt+L2/2
    z=0.0
    for j in xrange(0,Knum):
        for k in xrange(0,Knum):
            z = z+cks[j,k]*np.cos(j*np.pi*x/L1)*np.cos(k*np.pi*y/L2)/hk[j,k]
    return z

#hk=hkfunc()

data = genfromtxt('/home/kt-fitz/human-ergodicObj/DI_coefficients.csv',delimiter=",",dtype=float)
hk = data[0:5,0:5]
phik = data[5:10]
ck = data[11:16]

data=genfromtxt('/home/kt-fitz/human-ergodicObj/DIergtest.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
x1 = data[0:-1,1]
x2 = data[0:-1,3]

#xi, yi = np.mgrid[-L1/2:L1/2:NBINS*1j, -L2/2:L2/2:NBINS*1j]
xi, yi = np.mgrid[-20:15:NBINS*1j, -10:20:NBINS*1j]
zphik = fourierRecon(xi.flatten(), yi.flatten(),phik,hk)
zck = fourierRecon(xi.flatten(), yi.flatten(),ck,hk)

#domain = np.transpose(np.vstack([np.linspace(-L1/2,L1/2,NBINS),np.linspace(-L2/2,L2/2,NBINS)]))
ref = multivariate_normal.pdf(np.dstack((xi,yi)), mean = [-2,-2], cov=[[1.0,0.0],[0.0,1.]])
#delta = multivariate_normal.pdf(np.dstack((xi,yi)), mean = [0,0], cov=[[0.0025,0.],[0.,0.01]])


pltarea=[-L1/2,L1/2,-L2/2,L2/2]
ax = plt.subplot()
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(16) 
plt.pcolormesh(xi, yi, ref.reshape(xi.shape),norm=colors.Normalize(vmin=0,vmax=0.15))
plt.title("Reference Distribution",**title_font)
plt.xlabel ( r"$x$",**axis_font)
plt.ylabel ( r"$y$",rotation='horizontal',**axis_font)
plt.margins(0)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=fontsz)
cbar.ax.set_ylabel('Density',fontsize=16)
plt.show()

ax = plt.subplot()
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(16) 
    #plt.hexbin(xsob, ysob,gridsize = 40, extent =pltarea,norm=colors.Normalize(vmin=0,vmax=16))
plt.pcolormesh(xi, yi, zphik.reshape(xi.shape))#,norm=colors.Normalize(vmin=0,vmax=0.9))
plt.plot(x1,x2,'k',linewidth=2)
#plt.yticks(np.array([-3,0,3]))
plt.title("Fourier Recon of Phik's", **title_font)
plt.xlabel ( r"$x$",**axis_font)
plt.ylabel ( r"$y$",**axis_font)
plt.margins(0)
plt.axes().set_aspect('equal')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=fontsz)
cbar.ax.set_ylabel('Density',fontsize=16)
plt.show()

ax = plt.subplot()
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(16) 
    #plt.hexbin(xsob, ysob,gridsize = 40, extent =pltarea,norm=colors.Normalize(vmin=0,vmax=16))
plt.pcolormesh(xi, yi, zck.reshape(xi.shape))#,norm=colors.Normalize(vmin=0,vmax=0.9))
plt.plot(x1,x2,'k',linewidth=2)
#plt.yticks(np.array([-3,0,3]))
plt.title("Fourier Recon of Ck's", **title_font)
plt.xlabel ( r"$x$",**axis_font)
plt.ylabel ( r"$y$",**axis_font)
plt.margins(0)
plt.axes().set_aspect('equal')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=fontsz)
cbar.ax.set_ylabel('Density',fontsize=16)
plt.show()
