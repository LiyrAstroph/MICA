import numpy as np 
import matplotlib.pyplot as plt 

nbuilt = 30000
mcmc = np.loadtxt("../data/mcmc_con.txt")

ntheta = mcmc.shape[1]

fig = plt.figure()

ax = fig.add_subplot(121)
ax.hist((mcmc[nbuilt:, 1]+0.5*mcmc[nbuilt:, 2]-0.5*np.log(2.0))/np.log(10.0), 100)

ax = fig.add_subplot(122)
ax.hist(mcmc[nbuilt:, 2]/np.log(10.0), 100)

plt.show()