import numpy as np 
import matplotlib.pyplot as plt 

nbuilt = 30000
mcmc = np.loadtxt("../data/mcmc.txt")

ntheta = mcmc.shape[1]

fig = plt.figure()

ax = fig.add_subplot(231)
ax.hist((mcmc[nbuilt:, 1]+0.5*mcmc[nbuilt:, 2]-0.5*np.log(2.0))/np.log(10.0), 100)

ax = fig.add_subplot(232)
ax.hist(mcmc[nbuilt:, 2]/np.log(10.0), 100)

ax = fig.add_subplot(233)
ax.hist(mcmc[nbuilt:, 3], 100)

ax = fig.add_subplot(234)
ax.hist(mcmc[nbuilt:, 4]/np.log(10.0), 100)

ax = fig.add_subplot(235)
ax.hist(mcmc[nbuilt:, 5], 100)

plt.show()