import numpy as np
import scipy.integrate as scint
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.interpolate import interp1d
from scipy.stats import multivariate_normal
from scipy.stats import gaussian_kde

L1 = 10.0#2200#2*np.pi#
L2 =10.0#2200#
Knum = 30-1
NBINS=300
fontsz=14
title_font = {'fontname':'Liberation Sans', 'size':'20', 'color':'black', 'weight':'normal'}#,'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font = {'fontname':'Liberation Sans', 'size':'18'}#'21'}

def ckfunc(tlist,q,hk):
    x1 = q[:,0]#+L1/2
    x2 = q[:,1]#+L2/2
    fx1 = interp1d(tlist,x1)
    fx2 = interp1d(tlist,x2)
    cklist = np.zeros((Knum,Knum))
    for j in range(0,Knum):
        for k in range(0,Knum):
            fk=lambda t: np.cos(j*np.pi*fx1(t)/L1)*np.cos(k*np.pi*fx2(t)/L2)/hk[j,k]
            ck=scint.quad(fk,0,tlist[-1])[0]
            cklist[j,k] = ck/tlist[-1]
    return cklist

def fourierRecon(xt,yt,cks,hk):
    x=xt#-L1/2.
    y=yt#-L2/2.
    z=0.0
    for j in xrange(0,Knum):
        for k in xrange(0,Knum):
            z = z+cks[j,k]*np.cos(j*np.pi*x/L1)*np.cos(k*np.pi*y/L2)/hk[j,k]
    return z

data = genfromtxt('/home/kt-fitz/human-ergodicObj/DI_coefficients.csv',delimiter=",",dtype=float)
#data = genfromtxt('/home/kt-fitz/human-ergodicObj/CP_coefficients.csv',delimiter=",",dtype=float)
Knum=np.shape(data)[1]-1
hk = data[0:Knum+1,0:Knum+1]
phik = data[Knum+1:Knum+Knum+2]
ck = data[Knum+Knum+2:Knum+Knum+Knum+4]

data=genfromtxt('/home/kt-fitz/human-ergodicObj/DIergtest.csv',delimiter=",",dtype=float)
#data=genfromtxt('/home/kt-fitz/human-ergodicObj/CPergtest.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
tlist = data[0:-1,0]
x1 = data[0:-1,1]+(L1/2.)
x2 = data[0:-1,3]+(L2/2.)#data[0:-1,2]+(L2/2.)#
cost = data[0:-1,-1]

xj, yj = np.mgrid[-L1/2:L1/2:NBINS*1j, -L2/2:L2/2:NBINS*1j]
xi, yi = np.mgrid[0:L1:NBINS*1j, 0:L2:NBINS*1j]
zphik = fourierRecon(xi.flatten(), yi.flatten(),phik,hk)
zck = fourierRecon(xi.flatten(), yi.flatten(),ck,hk)
#zck2 = fourierRecon(xi.flatten(), yi.flatten(),ck2,hk)

#domain = np.transpose(np.vstack([np.linspace(-L1/2,L1/2,NBINS),np.linspace(-L2/2,L2/2,NBINS)]))
ref = multivariate_normal.pdf(np.dstack((xj,yj)), mean = [-2,1], cov=[[0.5,0.0],[0.0,0.5]])+multivariate_normal.pdf(np.dstack((xj,yj)), mean = [3.,-1], cov=[[0.5,0.0],[0.0,0.5]])
#delta = multivariate_normal.pdf(np.dstack((xi,yi)), mean = [0,0], cov=[[0.0025,0.],[0.,0.01]])

plt.figure()
plt.plot(tlist,cost)

plt.figure()
ax = plt.gca()
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(16) 
    #plt.hexbin(xsob, ysob,gridsize = 40, extent =pltarea,norm=colors.Normalize(vmin=0,vmax=16))
plt.pcolormesh(xi, yi, zphik.reshape(xi.shape))#,norm=colors.Normalize(vmin=0,vmax=0.3))
plt.plot(x1,x2,'ko',linewidth=2, markersize=1)
#plt.yticks(np.array([-3,0,3]))
plt.title("Fourier Recon of Phik's", **title_font)
plt.xlabel ( r"$x$",**axis_font)
plt.ylabel ( r"$y$",**axis_font)
plt.margins(0)
plt.axes().set_aspect('equal')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=fontsz)
cbar.ax.set_ylabel('Density',fontsize=16)


plt.figure()
ax = plt.gca()
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(16) 
    #plt.hexbin(xsob, ysob,gridsize = 40, extent =pltarea,norm=colors.Normalize(vmin=0,vmax=16))
plt.pcolormesh(xi, yi, zck.reshape(xi.shape))#,norm=colors.Normalize(vmin=0,vmax=0.3))
plt.plot(x1,x2,'ko',markersize=1)
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
"""
plt.figure()
ax = plt.gca()
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(16) 
    #plt.hexbin(xsob, ysob,gridsize = 40, extent =pltarea,norm=colors.Normalize(vmin=0,vmax=16))
plt.pcolormesh(xi, yi, zck2.reshape(xi.shape))#,norm=colors.Normalize(vmin=0,vmax=0.3))
plt.plot(x1,x2,'ko',markersize=1)
#plt.yticks(np.array([-3,0,3]))
plt.title("Fourier Recon of Ck's", **title_font)
plt.xlabel ( r"$x$",**axis_font)
plt.ylabel ( r"$y$",**axis_font)
plt.margins(0)
plt.axes().set_aspect('equal')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=fontsz)
cbar.ax.set_ylabel('Density',fontsize=16)
"""