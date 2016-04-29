import numpy as np 
import matplotlib.pyplot as plt 

fig = plt.figure(1, figsize=(15, 6))

mcmc = np.loadtxt("mychain0.dat")

ax = fig.add_subplot(121)
ax.hist(mcmc[:, 0]/np.log(10.0), 100, normed=True, range=(-1.0, 0.2))

lim1 = ax.get_xlim()

ax = fig.add_subplot(122)
ax.hist(mcmc[:, 1]/np.log(10.0), 100, normed=True, range=(0.0, 3.0))

lim2 = ax.get_xlim()

nbuilt = 30000
mcmc = np.loadtxt("../data/mcmc_con.txt")


ntheta = mcmc.shape[1]


fig = plt.figure(2, figsize=(15, 6))

ax = fig.add_subplot(121)
ax.hist((mcmc[nbuilt:, 1]+0.5*mcmc[nbuilt:, 2]-0.5*np.log(2.0))/np.log(10.0), 100, normed=True, range=(-1.0, 0.2))

ax.set_xlim(lim1[0], lim1[1])

ax = fig.add_subplot(122)
ax.hist(mcmc[nbuilt:, 2]/np.log(10.0), 100, normed=True, range=(0.0, 3.0))

ax.set_xlim(lim2[0], lim2[1])

plt.show()
