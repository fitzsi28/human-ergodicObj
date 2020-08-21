import numpy as np
import scipy.integrate as scint
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.stats import multivariate_normal
import cv2

fontsz=14
title_font = {'fontname':'Liberation Sans', 'size':'20', 'color':'black', 'weight':'normal'}#,'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font = {'fontname':'Liberation Sans', 'size':'18'}#'21'}


data=genfromtxt('/home/kt-fitz/human-ergodicObj/DIdkltest.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
tlist = data[0:-1,0]
x1 = (data[50:-1,1]+0.5)*2200
x2 = 2200-(data[50:-1,3]+0.5)*2200
X = np.stack([x1,x2],axis=1)

img = cv2.imread("apple.png",cv2.IMREAD_GRAYSCALE)
plt.figure()
plt.imshow(img, cmap = 'gray', interpolation = 'bicubic')
plt.scatter(x1,x2,marker='o')
plt.title("X(t)", **title_font)
plt.xlabel ( r"$x$",**axis_font)
plt.ylabel ( r"$y$",**axis_font)
plt.margins(0)
#cbar = plt.colorbar()
#cbar.ax.tick_params(labelsize=fontsz)
#cbar.ax.set_ylabel('Density',fontsize=16)

#plt.figure()
#plt.plot(tlist,data[0:-1,7])

plt.show()
