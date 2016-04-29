import numpy as np 
import matplotlib.pyplot as plt 
import corner

mcmc = np.loadtxt("mychain1.dat")

ntheta = mcmc.shape[1]

fig = plt.figure(1, figsize=(15, 6))

ax = fig.add_subplot(231)
ax.hist(mcmc[:, 0]/np.log(10.0), 100, normed=True, range=(-0.9, -0.1))

ax = fig.add_subplot(232)
ax.hist(mcmc[:, 1]/np.log(10.0), 100, normed=True, range=(0.0, 2.0))

ax = fig.add_subplot(234)
ax.hist(mcmc[:, 2], 100, normed=True, range=(1.0, 2.8))

ax = fig.add_subplot(235)
ax.hist(mcmc[:, 3], 100, normed=True, range=(0.0, 1.2))

ax = fig.add_subplot(236)
ax.hist(mcmc[:, 4], 100, normed=True, range=(5, 13))


mcmc = np.loadtxt("../data/mcmc.txt")

ntheta = mcmc.shape[1]

nb = 20000
fig = plt.figure(2, figsize=(15, 6))

ax = fig.add_subplot(231)
ax.hist( (mcmc[nb:, 1]+0.5*mcmc[nb:, 2]-0.5*np.log(2.0))/np.log(10.0), 100, normed=True, range=(-0.9, -0.1))

ax = fig.add_subplot(232)
ax.hist(mcmc[nb:, 2]/np.log(10), 100, normed=True, range=(0.0, 2.0))

ax = fig.add_subplot(234)
ax.hist(mcmc[nb:, 5], 100, normed=True, range=(1.0, 2.8))

ax = fig.add_subplot(235)
ax.hist(mcmc[nb:, 3], 100, normed=True, range=(0.0, 1.2))

ax = fig.add_subplot(236)
ax.hist(mcmc[nb:, 4], 100, normed=True, range=(5, 13))

plt.show()
