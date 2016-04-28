import numpy as np
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pdf = PdfPages('hist.pdf')

data = np.loadtxt("../data/mcmc.txt")
nb = 50000
nt=data.shape[0]

data[:, 1] = (data[:, 1] + 0.5*data[:, 2] - 0.5*np.log(2.0))/np.log(10.0)
data[:, 2] = data[:, 2]/np.log(10.0)
data[:, 4] = np.exp(data[:, 4])


print(nt)

nc = 2

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
  sort=np.sort(data[nb:nt, i])
# print(data[40000-1, i], data[40100, i], data[40200, i])
  mean = np.mean(data[nb:nt,i])
  std = np.std(data[nb:nt,i])
  low = np.percentile(data[nb:nt,i], 15.85, interpolation='linear')
  up = np.percentile(data[nb:nt,i], 100.0-15.85, interpolation='linear')
  #if((i<=3) or (i==nc+4)):
  #  hrange=(data[nb:nt,i].min(),data[nb:nt,i].max())
  #else:
  hrange=(data[nb:nt,i].min(), data[nb:nt,i].max())
#hrange = (max(mean - 5.0*std, -15.0), mean + 5.0*std)
  hist, bin_edges = np.histogram(data[nb:nt,i], density=False, bins=100, range=hrange)
  bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
  eps_data = np.zeros(len(hist)) + 1.0

  params = Parameters()
  idx = np.argmax(hist);
  print(idx, hist[idx], )
  params.add('amp', value=np.max(hist), min=0.0)
  params.add('cen', value=bin_centres[idx], min=bin_centres[0], max=bin_centres[-1])
  params.add('wid', value=std, min=0.0)
  idx = np.where( hist > np.max(hist)*0.0)
  #if(i==5):
  #  idx = np.where( hist > np.max(hist)*0.5)

  out = minimize(residual, params, args=(bin_centres[idx[0]], hist[idx[0]], eps_data[idx[0]]), method='leastsq')

# Get the fitted curve
  bin_fit = np.linspace(bin_centres[0], bin_centres[-1], 100)
  hist_fit = gaussian(bin_fit, out.params['amp'].value, out.params['cen'].value,out.params['wid'].value)

  fig=plt.figure()

  #plt.plot(bin_centres, hist, label='Test data')
  plt.hist(data[nb:nt,i], normed=False, bins=100, range=hrange)
  plt.plot(bin_fit, hist_fit, label='Fitted data')
  
  print(np.min(data[nb:nt, i]), np.max(data[nb:nt, i]))
  plt.xlim(np.min(data[nb:nt, i]), np.max(data[nb:nt, i]))
# Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
  mean = out.params['cen'].value
  par[i-1, 0] = mean
  par[i-1, 1] = mean - low
  par[i-1, 2] = up - mean
  pdf.savefig(fig)
  plt.close()

pdf.close()

print(par)

np.savetxt("par.txt", par, fmt="%f\t%f\t%f")
