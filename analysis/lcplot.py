#
# MICA (Multiple and Inhomogeneous Component Analysis)
# A Non-parametric Approach to Constrain the Transfer Function in Reverberation Mapping
#
# Yan-Rong Li, liyanrong@mail.ihep.ac.cn
#
# Reference: Li, Y.-R. et al. 2016, arXiv:1608.03741
#
#
# lcplot.py
# plot the best recovered light curves and transfer function.


import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

plt.rc('text', usetex=True)

con = np.loadtxt("../data/con_test.txt");
hb = np.loadtxt("../data/line_test.txt")

con_rec = np.loadtxt("../data/sall_con.txt")
hb_rec = np.loadtxt("../data/sall_line.txt")

tf = np.loadtxt("../data/transfer.txt")

fig = plt.figure(figsize=(10, 4))

ax1=fig.add_axes((0.1, 0.5, 0.5, 0.4))

ax1.fill_between(con_rec[:, 0], con_rec[:, 1]+con_rec[:, 2], con_rec[:, 1]-con_rec[:, 2], where=None, color='grey')
ax1.errorbar(con[:, 0], con[:, 1], yerr=con[:, 2], ls='none', marker='None', markersize=4, color='k', capsize=1)
ax1.plot(con_rec[:, 0], con_rec[:, 1], linewidth=1, color='b')

[i.set_visible(False) for i in ax1.get_xticklabels()]

ax1.set_ylabel(r'$F_{\rm 5100}$')


ax2=fig.add_axes((0.1, 0.1, 0.5, 0.4))

ax2.fill_between(hb_rec[:, 0], hb_rec[:, 1]+hb_rec[:, 2], hb_rec[:, 1]-hb_rec[:, 2], where=None, color='grey')
ax2.errorbar(hb[:, 0], hb[:, 1], yerr=hb[:, 2], ls='none', marker='None', markersize=4, color='k', capsize=1)
ax2.plot(hb_rec[:, 0], hb_rec[:, 1], color='b')

ax2.set_xlabel(r'$\rm HJD$')
ax2.set_ylabel(r'$F_{\rm H\beta}$')

xlim1 = ax1.get_xlim()
xlim2 = ax2.get_xlim()
xlim =[0,0]
xlim[0] = min(xlim1[0], xlim2[0])
xlim[1] = max(xlim1[1], xlim2[1])

ax2.set_xlim(xlim)
ax1.set_xlim(xlim)



ax3 = fig.add_axes((0.7, 0.1, 0.25, 0.8))

#ax3.plot(tf_input[:, 0], tf_input[:, 1], color='k')
ax3.plot(tf[:, 0], tf[:, 1], linewidth=2)
ax3.fill_between(tf[:, 0], tf[:, 2] , tf[:, 3], where=None, color='grey')


ax3.set_xlabel(r'$\tau$')
ax3.set_ylabel(r'$\rm Transfer\ Function$')

plt.savefig("fig.pdf", bbox_inches='tight')
