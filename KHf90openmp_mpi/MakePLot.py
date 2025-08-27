#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re
import glob

# usage
# ./MakePlot.py snapshots 1

dirname = sys.argv[1]
step_s = int(sys.argv[2])

print("dir ",dirname)
print("step ",step_s)

def makedirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)

makedirs("./images")
fig = plt.figure()  
plt.xlim(0, 1)     
plt.ylim(0, 1)
plt.xlabel("x axis") 
plt.ylabel("y axis") 

xmin = -0.5
xmax =  0.5
ymin = -0.5
ymax =  0.5

for istep in range(step_s,step_s+1):
    filelists = glob.glob(dirname + "/snap*-%05d.xss"%(istep))
    print("lists",filelists)
    for foutname in filelists:  
        print("making plot ",foutname)
        with open(foutname, 'r') as data_file:
            line = data_file.readline();
            attributes1 = re.findall("\d+\.\d+", line)

            line = data_file.readline();
            attributes2 = re.findall("\d+", line)

        time = float(attributes1[0]) 
        nx = int(attributes2[0])
        ny = int(attributes2[1])

        data = np.loadtxt(foutname)

        x   = data[:,0].reshape(ny,nx)
        y   = data[:,1].reshape(ny,nx)
        den = data[:,2].reshape(ny,nx)
        vx  = data[:,3].reshape(ny,nx)
        vy  = data[:,4].reshape(ny,nx)
        pre = data[:,5].reshape(ny,nx)
        #    sca = data[:,10].reshape(ny,nx)
        im=plt.contourf(x,y,den,100)

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    pg00 = plt.text(0.5*(xmin+xmax),ymax*1.1,r"$\mathrm{time}=%.2f$"%(time),horizontalalignment="center")


    if istep == step_s: 
        plt.colorbar(im,orientation="vertical")

    plt.savefig("./images/den%05d.png"%(istep),dpi=300)

#ani = animation.ArtistAnimation(fig, graph_list, interval=200) 
#print("making animation file", fname_anime)
#ani.save(fname_anime, writer="imagemagick")
plt.show()
