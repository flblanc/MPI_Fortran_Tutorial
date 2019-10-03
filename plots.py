#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt 

import sys

nreplicas = int(sys.argv[1])

niter=10
nsteps=100


fig, ax = plt.subplots()

for i in xrange(nreplicas):

    xx = np.genfromtxt("replica%i.log" %i)
    ax.plot(xx[:,0], xx[:,1])

for k in xrange(1,niter+1):

    ax.axvline(x=k*nsteps)

#ax.set_ylim(4.0,6.0)


fig.savefig("mc.png")
