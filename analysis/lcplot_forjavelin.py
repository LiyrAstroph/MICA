import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import numpy as np 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

con = np.loadtxt("../data/con_test.txt") 
line = np.loadtxt("../data/line_test.txt")

con_sim = np.loadtxt("../data/sall_con.txt") 
line_sim = np.loadtxt("../data/sall_line.txt")

fig  = plt.figure(figsize=(8, 4))
height = 0.85/2.0

marker="o"
ms=4
ls='None'
lw=1

i = 1
nlc = 2
ax = fig.add_axes([0.10, 0.1+i*height, 0.85, height])
mfc = cm.jet(i/(nlc-1.) if nlc > 1 else 0)
ax.plot(line_sim[:, 0], line_sim[:, 1],
                    color=mfc, ls="-", lw=2)
                    
ax.fill_between(line_sim[:, 0],
                    y1=line_sim[:, 1]+line_sim[:, 2], 
                    y2=line_sim[:, 1]-line_sim[:, 2], 
                    color=mfc, alpha=0.5)
ax.errorbar(line[:, 0], line[:, 1], 
                            yerr=line[:, 2], 
                            ecolor='k', marker=marker, ms=ms, mfc=mfc, mec='k', ls=ls, lw=lw)
ax.set_xlim(495, 630)  
ax.set_ylim(9.0, 19.0)

i = 0
nlc = 2
ax = fig.add_axes([0.10, 0.1+i*height, 0.85, height])
mfc = cm.jet(i/(nlc-1.) if nlc > 1 else 0)
ax.plot(con_sim[:, 0], con_sim[:, 1],
                    color=mfc, ls="-", lw=2)
                    
ax.fill_between(con_sim[:, 0],
                    y1=con_sim[:, 1]+con_sim[:, 2], 
                    y2=con_sim[:, 1]-con_sim[:, 2], 
                    color=mfc, alpha=0.5)
ax.errorbar(con[:, 0], con[:, 1], 
                            yerr=con[:, 2], 
                            ecolor='k', marker=marker, ms=ms, mfc=mfc, mec='k', ls=ls, lw=lw)

ax.set_xlim(495, 630)
ax.set_ylim(1.0, 2.1)                        
plt.show()
