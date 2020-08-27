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

img = cv2.imread("apple.png",cv2.IMREAD_GRAYSCALE)

            
data=genfromtxt('/home/kt-fitz/data/cpp2/s10_c_set01-clean.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
tlist = data[0:-1,0]
x2 = (data[0:-1,-2]+0.5)*2200
x1 = 2200-(data[0:-1,-3]+0.5)*2200
X = np.stack([x1,x2],axis=1)

plt.figure()
plt.ylim((2200,0))
plt.xlim((0,2200))
plt.imshow(img, cmap = 'gray', interpolation = 'bicubic')
plt.xlabel ( r"$x\ (pixels)$",**axis_font)
plt.ylabel ( r"$y\ (pixels)$",**axis_font)
plt.title("Baseline Drawing", **title_font)
"""
#s10_banana_c
plt.plot(x1[23050:23500],x2[23050:23500],'-',color='green')
plt.plot(x1[38050:38500],x2[38050:38500],'-',color='green')
plt.plot(x1[18050:18500],x2[18050:18500],'-',color='green')
"""
"""
#s10_apple_v
plt.plot(x1[50:550],x2[50:550],'-',color='green')
#s10_apple_p
plt.plot(x1[1000:1500],x2[1000:1500],'-',color='green')
plt.plot(x1[3050:3500],x2[3050:3500],'-',color='green')
plt.plot(x1[6050:6500],x2[6050:6500],'-',color='green')
"""
#s10_apple_c
plt.plot(x1[7050:7500],x2[7050:7500],'-',color='green')
plt.plot(x1[8050:8500],x2[8050:8500],'-',color='green')
plt.plot(x1[22050:22500],x2[22050:22500],'-',color='green')


plt.show()
