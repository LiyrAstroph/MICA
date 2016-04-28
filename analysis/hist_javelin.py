import numpy as np 
import matplotlib.pyplot as plt 
import corner

mcmc = np.loadtxt("mychain1.dat")

ntheta = mcmc.shape[1]

fig = plt.figure(1)

ax = fig.add_subplot(231)
ax.hist(mcmc[:, 0], 100, normed=True)

ax = fig.add_subplot(232)
ax.hist(mcmc[:, 1], 100, normed=True)

ax = fig.add_subplot(234)
ax.hist(mcmc[:, 2], 100, normed=True)

ax = fig.add_subplot(235)
ax.hist(mcmc[:, 3], 100, normed=True)

ax = fig.add_subplot(236)
ax.hist(mcmc[:, 4], 100, normed=True)


mcmc = np.loadtxt("../data/mcmc.txt")

ntheta = mcmc.shape[1]

nb = 20000
fig = plt.figure(2)

ax = fig.add_subplot(231)
ax.hist( (mcmc[nb:, 1]+0.5*mcmc[nb:, 2]-0.5*np.log(2.0))/np.log(10.0), 100, normed=True)

ax = fig.add_subplot(232)
ax.hist(mcmc[nb:, 2]/np.log(10), 100, normed=True)

ax = fig.add_subplot(234)
ax.hist(mcmc[nb:, 5], 100, normed=True)

ax = fig.add_subplot(235)
ax.hist(mcmc[nb:, 3], 100, normed=True)

ax = fig.add_subplot(236)
ax.hist(mcmc[nb:, 4], 100, normed=True)

plt.show()
