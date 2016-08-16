import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

plt.rc('text', usetex=True)

con = np.loadtxt("../data/ngc6814_con_bentz2009.txt");
hb = np.loadtxt("../data/ngc6814_hb_bentz2009.txt")

con_rec = np.loadtxt("../data/sall_con.txt")
hb_rec = np.loadtxt("../data/sall_line.txt")

tf = np.loadtxt("../data/transfer.txt")
#tf_input = np.loadtxt("../data/sim_tf.txt")

fig = plt.figure(figsize=(10, 4))

ax1=fig.add_axes((0.1, 0.5, 0.5, 0.4))

ax1.fill_between(con_rec[:, 0], con_rec[:, 1]+con_rec[:, 2], con_rec[:, 1]-con_rec[:, 2], where=None, color='grey')
ax1.errorbar(con[:, 0], con[:, 1], yerr=con[:, 2], ls='none', marker='None', markersize=4, color='k', capsize=1)
ax1.plot(con_rec[:, 0], con_rec[:, 1], linewidth=1, color='b')

[i.set_visible(False) for i in ax1.get_xticklabels()]

ax1.set_ylabel(r'$F_{\rm 5100}$')
#ax1.yaxis.set_major_locator(MultipleLocator(0.5))
#ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(MultipleLocator(10))

#ax1.text(5420.0, 2.8,r'$\rm PG\ 2130+099$')

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

#ax2.yaxis.set_major_locator(MultipleLocator(0.1))
#ax2.yaxis.set_minor_locator(MultipleLocator(0.05))
ax2.xaxis.set_minor_locator(MultipleLocator(10))

ax3 = fig.add_axes((0.7, 0.1, 0.25, 0.8))

#ax3.plot(tf_input[:, 0], tf_input[:, 1], color='k')
ax3.plot(tf[:, 0], tf[:, 1], linewidth=2)
ax3.fill_between(tf[:, 0], tf[:, 2] , tf[:, 3], where=None, color='grey')
#ax3.plot(tf[:, 0], tf[:, 1:tf.shape[1]-4], ls=':', color='red', linewidth=0.1)
#ax3.plot(tf[:, 0], tf[:, 2])
#ax3.plot(tf[:, 0], tf[:, 3])
#ax3.plot(tf[:, 0], tf[:, 4])
#ax3.plot(tf[:, 0], tf[:, 5])

ax3.set_xlabel(r'$\tau$')
ax3.set_ylabel(r'$\rm Transfer\ Function$')
#ax3.yaxis.set_minor_locator(MultipleLocator(0.005))
#ax3.xaxis.set_minor_locator(MultipleLocator(2.5))

plt.savefig("fig.pdf", bbox_inches='tight')
