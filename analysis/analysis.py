import numpy as np
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


pdf = PdfPages('hist.pdf')

nbin = 50
ntau = 100
tau1 = 0.0
tau2 = 50.0
nc = 4
nb = 20000
data = np.loadtxt("../data/mcmc.txt", skiprows=nb)
nt=data.shape[0]

#data[:, 1] = (data[:, 1] + 0.5*data[:, 2] - 0.5*np.log(2.0))/np.log(10.0)
#data[:, 2] = data[:, 2]/np.log(10.0)


print(nt)

par=np.zeros((nc+4, 3))

def gaussian(x, amp, cen, wid):
    return (amp) * np.exp(-(x-cen)**2 /(2*wid**2))

def residual(params, x, data, eps_data):
    amp = params['amp'].value
    cen = params['cen'].value
    wid = params['wid'].value

    model = amp * np.exp(-(x-cen)**2 /(2*wid**2))

    return (data-model)/eps_data

for i in range(1,nc+4+1):
  print(i)
  sort=np.sort(data[:, i])
# print(data[40000-1, i], data[40100, i], data[40200, i])
  mean = np.mean(data[:,i])
  std = np.std(data[:,i])
  low = np.percentile(data[:,i], 15.85, interpolation='linear')
  up = np.percentile(data[:,i], 100.0-15.85, interpolation='linear')
  #if((i<=3) or (i==nc+4)):
  #  hrange=(data[:,i].min(),data[:,i].max())
  #else:
  hrange=[data[:,i].min(), data[:,i].max()]
#hrange = (max(mean - 5.0*std, -15.0), mean + 5.0*std)
  if(i==1000):
    print("i-6")
    idxsel = np.where((data[:,i] >= -5.0) & (data[:,i] <= 10.0) )
  else:
    idxsel = np.where((data[:,i] >= hrange[0]) & (data[:,i] <= hrange[1]) )
  

  hist_data = data[idxsel[0],i]
  hist_range = (hist_data.min(), hist_data.max())
  print(i, hist_range)

  hist, bin_edges = np.histogram(hist_data, density=False, bins=nbin, range=hist_range)

  bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
  eps_data = np.zeros(len(hist)) + 1.0

  params = Parameters()
  idx = np.argmax(hist);
  print(i, idx, hist[idx], )
  params.add('amp', value=np.max(hist), min=0.0)
  params.add('cen', value=bin_centres[idx], min=bin_centres[0], max=bin_centres[-1])
  params.add('wid', value=std, min=0.0)
  idx = np.where( hist > np.max(hist)*0.0)
  #if(i==6):
  #  idx = np.where( hist > np.max(hist)*0.6)

  out = minimize(residual, params, args=(bin_centres[idx[0]], hist[idx[0]], eps_data[idx[0]]), method='leastsq')

# Get the fitted curve
  bin_fit = np.linspace(bin_centres[0], bin_centres[-1], nbin)
  hist_fit = gaussian(bin_fit, out.params['amp'].value, out.params['cen'].value,out.params['wid'].value)

  fig=plt.figure()

  #plt.plot(bin_centres, hist, label='Test data')
  plt.hist(hist_data, normed=False, bins=nbin, range=hist_range)
  plt.plot(bin_fit, hist_fit, label='Fitted data')
  
  print(np.min(data[:, i]), np.max(data[:, i]))
  plt.xlim(hist_range)
# Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
  mean = out.params['cen'].value
  par[i-1, 0] = mean
  par[i-1, 1] = mean - low
  par[i-1, 2] = up - mean
  pdf.savefig(fig)
  plt.close()

pdf.close()

print(par)

np.savetxt("../data/par.txt", par, fmt="%f\t%f\t%f")

fig = plt.figure()
ax = fig.add_subplot(111)

tau = np.linspace(tau1, tau2, ntau)
tf = np.zeros(ntau)
tau_grid = np.linspace(tau1, tau2, nc)

width = par[2, 0]
for i in range(nc):
  fk = np.exp(par[3+i, 0])/np.sqrt(2.0*np.pi)/width * np.exp( -0.5*(tau - tau_grid[i])**2/width/width)
  tf += fk
  ax.plot(tau, fk, "r:")


ax.plot(tau, tf, linewidth=2)

plt.show()
