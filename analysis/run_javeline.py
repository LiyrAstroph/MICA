# JAVELINE 
#

from javelin.zylc import get_data
from javelin.lcmodel import Cont_Model, Rmap_Model

javdata_con = get_data(["../data/con_test.txt",], names=["con",])
javdata_rm = get_data(["../data/con_test.txt","../data/line_test.txt"], names=["con","line"])
cont = Cont_Model(javdata_con)
cont.do_mcmc(fchain="mychain0.dat",nwalkers=100, nburn=100, nchain=100)
cont.show_hist(bins=100)
cont.get_hpd()
conthpd = cont.hpd
print(conthpd)

rmap1 = Rmap_Model(javdata_rm)
rmap1.do_mcmc(conthpd=conthpd, fchain="mychain1.dat", laglimit=[[0, 10],], nwalkers=100, nburn=100, nchain=200)
rmap1.show_hist()
rmap1.get_hpd()
rmap1hpd = rmap1.hpd
par_best = rmap1hpd[1,:]
print(par_best)

javdata_best =  rmap1.do_pred(par_best)
javdata_best.plot(set_pred=True, obs=javdata_rm)
