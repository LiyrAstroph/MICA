import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pdf = PdfPages('hist.pdf')

data = np.loadtxt("../data/mcmc.txt")
nb = 20000
nt=data.shape[0]

print(nt)

nc = 3

par=np.zeros((nc+4, 5))

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

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
  hist, bin_edges = np.histogram(data[nb:nt,i], density=True, bins=100, range=hrange)
  bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
  
  p0 = [hist.max(), mean, std]

#  coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)

# Get the fitted curve
#bin_fit = np.linspace(bin_centres[0], bin_centres[-1], 100)
#  hist_fit = gauss(bin_fit, *coeff)

  fig=plt.figure()

  plt.plot(bin_centres, hist, label='Test data')
#  plt.plot(bin_fit, hist_fit, label='Fitted data')
  
  print(np.min(data[nb:nt, i]), np.max(data[nb:nt, i]))
  plt.xlim(np.min(data[nb:nt, i]), np.max(data[nb:nt, i]))
# Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
#  print coeff[1], coeff[2]
#  par[i-1, 0] = coeff[1]
#  par[i-1, 1] = coeff[2]
#  mean = coeff[1]
#  par[i-1, 2] = mean
#  par[i-1, 3] = mean - low
#  par[i-1, 4] = up - mean
  pdf.savefig(fig)
  plt.close()

pdf.close()

print(par)

np.savetxt("par.txt", par[:, 2:5], fmt="%f\t%f\t%f")
