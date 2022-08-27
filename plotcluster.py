
import pandas as pd
import numpy as np
import re
import matplotlib
import matplotlib.pyplot as plt
import sys


import argparse


def setupArgs():
    parser = argparse.ArgumentParser(description='Generate cluster plot from csv')

    parser.add_argument('input', metavar='file', type=str, help='input file name',nargs='+')
    parser.add_argument('--size', dest='size', action='store', default=5,type=int,\
                        help='size of points')
    parser.add_argument('--colours', dest='colours', action='store', default="red,blue",
                        help='colours to use')
    parser.add_argument('--cartesian', dest='cartesian', action='store_true',\
                        default=False,\
                        help='plot in Cartesian coordinates rather than polar')
    parser.add_argument("--output",dest="output",action='store',required=True)
    args = parser.parse_args()
    return args

args = setupArgs()

default_colours="red,blue,black,magenta,orange,cyan,pink,yellow".split(",")

colours = args.colours.split(",")

if len(colours)<len(args.input):
    colours=colours+(default_colours-colours)

if args.output != "show":
    matplotlib.use('Agg')

ax=plt.subplot()
x="R"
y="G"
    
for i, fname in enumerate(args.input):
    data = pd.read_csv(fname,delim_whitespace=True,header=None,names=["R","G"])
    if args.cartesian:
        ax.scatter(data[x],data[y],s=args.size,c=colours[i],label=fname)        
    else:
        phi = np.arctan2(data[y],data[x])/(3.14159/2)*90
        r   = np.sqrt(np.power(data[x],2) + np.power(data[y],2))
        ax.scatter(phi,r,s=args.size,c=colours[i],label=fname)

output=args.output
ax.legend()
if output=="show":
    plt.show()
else:
    plt.savefig(output)

