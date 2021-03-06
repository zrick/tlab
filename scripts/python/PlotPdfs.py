#!/usr/bin/python3

import numpy as np
import struct
import sys
import matplotlib.pyplot as plt
import os

sizeofdata = 4 # in bytes
# sizeofdata = 1 # for gate files

# etype = ">" # big-endian
etype = "<" # little-endian

dtype = "f" # floating number
# dtype = 'B' # unsigned character, for gate files

# do not edit below this line

# getting data from stdin
if ( len(sys.argv) <= 1 ):
    print("Usage: python $0 dimensions level list-of-files.")
    quit()

ndim = int(sys.argv[1])                 # 1 for univariate, 2 for bivariate...
level = int(sys.argv[2])                # -1 to plot whole map instead of individual levels
setoffiles = sorted(sys.argv[3:])

# obtain size from first file
# the first data is the time and we skip it
# then comes the size and vertical coordinates
# the last level contains the global pdfs
# the last two entries per level are min/max values
fin = open(setoffiles[0], 'rb')
fin.seek( 4 )
raw = fin.read( (1+ndim)*4 )
ny  = struct.unpack((etype+'{}i').format(1+ndim), raw)[0]
nb  = struct.unpack((etype+'{}i').format(1+ndim), raw)[1:]
raw = fin.read( ny*sizeofdata )
y   = np.array(struct.unpack((etype+'{}'+dtype).format(ny), raw))
fin.close()

print("Files with {} bins and {} levels.".format(nb,ny))

# reading data
nb_size = np.prod(list(nb)) + 2 +2*(ndim-1)*nb[0]
a = np.zeros((nb_size*(ny+1)),dtype=float)
for file in setoffiles:
    print("Processing file {}. ".format(file),end="")
    fin = open(file, 'rb')
    raw = fin.read( 4 )                         # Read the time
    t = struct.unpack((etype+'f'), raw)[0]
    print("Time {:4.2f}...".format(t))
    fin.seek( (1+ndim)*4 + ny*sizeofdata, 1 )   # Skip the y-coordinates
    raw = fin.read()
    a   = a +np.array(struct.unpack((etype+'{}'+dtype).format(int(nb_size*(ny+1))), raw))
    fin.close()

a = a /len(setoffiles)
a = a.reshape((ny+1,nb_size))

# processing data
if ndim == 1:
    nb = nb[0]

    # # normalinzing histograms to obtain pdf (s.t. it integrates to 1 using midpoint rule)
    # for j in range(ny+1):
    #     samplesize = np.sum(a[j,:nb])
    #     samplestep = (a[j,nb+1]-a[j,nb]) /( nb -1 )
    #     if samplestep *samplesize > 0: # otherwise the pdf is zero, by construction in Tlab
    #         a[j,:nb] = a[j,:nb] /( samplesize *samplestep )

    if level > 0:
        var1 = np.linspace(a[level-1,nb],a[level-1,nb+1],num=nb)
        plt.plot( var1, a[level-1,:nb])
        plt.xlabel("var1")
        plt.ylabel("pdf, cavg")
        if level <= ny:
            plt.title(setoffiles[0]+" - height {:4.2f}".format(y[level-1]))
        else:
            plt.title(setoffiles[0]+" - global")
        plt.tight_layout(pad=0.1)
        plt.savefig("{}.pdf".format(setoffiles[0]))
        plt.show()

    else:
        # axis information
        var1 = np.zeros((ny,nb),dtype=float)
        y_ex = np.zeros((ny,nb),dtype=float)
        for j in range(ny):
            var1[j,:] = np.linspace(a[j,nb],a[j,nb+1],num=nb)
            y_ex[j,:] = y[j]
        #   var1[j,:] = var1[j,:]-0.5*nbstep # colormesh uses coordinates for the corners of the region

        # choose an interval to define the color range
        levels = 20                 # Default
        # nymin = int(ny/10)        # defined by intermediate interval
        # nymax = nymin *3
        # nymin = 0
        # nymax = int(2*ny/3)
        # levels=np.linspace(np.amin(a[nymin:nymax,:nb]),np.amax(a[nymin:nymax,:nb]),num=20)
        # levels=np.linspace(np.amin(a[:ny,:nb]),np.amax(a[:ny,:nb]),num=20) *0.5

        plt.contourf(var1,y_ex,a[:ny,:nb],levels)
        # plt.pcolormesh(var1,y_ex,a[:ny,:nb],levels)

        plt.xlabel("var1")
        plt.ylabel("height")
        plt.colorbar(label='pdf, cavg',format="%.2g")
        plt.tight_layout(pad=0.1)
        plt.savefig("{}.pdf".format(setoffiles[0]))
        plt.show()

if ndim == 2:
    if level > 0:
        # axis information
        var1 = np.empty((nb[1],nb[0]),dtype=float)
        var2 = np.empty((nb[1],nb[0]),dtype=float)
        for j in range(nb[1]):
            var1[j,:] = np.linspace(a[level-1,nb[0]*nb[1]    ],a[level-1,nb[0]*nb[1]+1        ],num=nb[0])
        for i in range(nb[0]):
            var2[:,i] = np.linspace(a[level-1,nb[0]*nb[1]+2+i],a[level-1,nb[0]*nb[1]+2+i+nb[0]],num=nb[1])

        # var1 = np.linspace( a[level-1,nb[0]*nb[1]],   a[level-1,nb[0]*nb[1]+1], num=nb[0] )
        # var2 = np.linspace( a[level-1,nb[0]*nb[1]+2], a[level-1,nb[0]*nb[1]+3], num=nb[1] )

        z = a[level-1,:nb[0]*nb[1]].reshape(nb[1],nb[0])
        surf=plt.contourf( var1, var2, z )
        #
        # ax = plt.axes(projection='3d')
        # ax.plot_surface( var1, var2, z, cmap='viridis', alpha=0.5 )
        # ax.scatter( var1, var2, color='k')
        # surf=ax.scatter( var1, var2, z, c=z, cmap='viridis' )
        #
        plt.xlabel("var1")
        plt.ylabel("var2")
        if level <= ny:
            plt.title(setoffiles[0]+" - height {:4.2f}".format(y[level-1]))
        else:
            plt.title(setoffiles[0]+" - global")
        plt.colorbar(surf,label='pdf, cavg',format="%.2g")
        plt.tight_layout(pad=0.1)
        plt.savefig("{}.pdf".format(setoffiles[0]))
        plt.show()

    else:
        print("Undveloped option.")
        quit()
