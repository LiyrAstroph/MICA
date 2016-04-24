import numpy as np 
import matplotlib.pyplot as plt 
import corner

nbuilt = 30000
mcmc = np.loadtxt("../data/mcmc.txt")

ntheta = mcmc.shape[1]

#fig = plt.figure()

corner.corner(mcmc[nbuilt:, 1:])


plt.show()