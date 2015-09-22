#! /usr/bin/python2
# -*- coding: utf-8 -*-
import numpy as np
from pylab import *
import matplotlib
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm
import math as maths #localisation
import sys
import argparse


#from  matplotlib.mlab import griddata


parser = argparse.ArgumentParser(description='Plot data in histograms');

parser.add_argument('-b',help='number of bins (defaults to 1000)',type=int,default=1000)
parser.add_argument('-i',help='input file(s) (defaults to stdin)',type=argparse.FileType('r'),default=[sys.stdin],nargs='*')
parser.add_argument('-lx',help='lower x limit',type=float,default=-sys.float_info.max)
parser.add_argument('-ux',help='upper x limit',type=float,default=sys.float_info.max)

parser.add_argument('-ly',help='lower y limit',type=float,default=-sys.float_info.max)
parser.add_argument('-uy',help='upper y limit',type=float,default=sys.float_info.max)

args = parser.parse_args()

matplotlib.rcParams['text.latex.unicode']=True #for greek letters

def convertData(data,xLowerLimit,xUpperLimit,yLowerLimit,yUpperLimit):
    xl, yl, z1l = [], [], []
    for thing in data:
        try:
            strings = thing.split()
            x = float(strings[0])
            y = float(strings[1])
            z1 = float(strings[2])
            
            if maths.isnan(x) or maths.isnan(y) or maths.isnan(z1): #or maths.isnan(z2):
                raise ValueError
            if (xLowerLimit < x < xUpperLimit) and (yLowerLimit < y < yUpperLimit):
                xl.append(x)
                yl.append(y)
                z1l.append(z1)
        except ValueError:
            continue
    return xl, yl, z1l

for inFile in args.i:
    X, Y, Z1 = convertData(inFile.readlines(),args.lx,args.ux,args.ly,args.uy) 

    xMin, xMax = min(X), max(X)
    yMin, yMax = min(Y), max(Y)

    xi = np.linspace(xMin,xMax, args.b) 
    yi = np.linspace(yMin,yMax, args.b) 


    zi = griddata(X, Y, Z1, xi, yi)

    zi.swapaxes(0,1)
    
    plt.tick_params(labelsize = 15)

    p = plt.imshow(zi,extent= [xMin, xMax, yMin, yMax], origin='lower', interpolation = 'none')
    colourBar = plt.colorbar(p, shrink = 0.5, aspect = 5)
    colourBar.ax.tick_params(labelsize = 15)
    colourBar.set_label('$V$', fontsize = 20)
   # colourBar.set_clim(vmin=0, vmax=5)
    



    plt.xlabel('$x/[m]$', fontsize = 23)
    plt.ylabel('$y/[m]$', fontsize = 23)
    
    axis('image')
    
    show()
    
    plt.show()

