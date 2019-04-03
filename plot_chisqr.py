import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


shear_catalog = np.load('shear_ctlg.npy')
chi_sqr = np.zeros(101)

sample_fg = 10000
sample_bg = 2000
sigma = 0.25
for ii in range(sample_fg):
    noise = np.random.normal(0, sigma, sample_bg)
    g_meas = shear_catalog[50,ii] + noise
    chi2 = np.zeros(101)
    for i in range(101):
        chi2[i] = np.sum((g_meas - shear_catalog[i,ii])**2)/sigma**2
    chi_sqr = chi_sqr + chi2


    


fig, ax = plt.subplots(figsize=[10,10])
x = np.linspace(8,12,101)
ax.plot(x,chi_sqr)
#ax.plot(x,y_2,color='black')
#ax.set_xscale('log')
#ax.set_xticks(np.linspace(0,10,6))
#ax.set_yticks([1e-6,1e-5,1e-4,1e-3])
ax.set_xlabel(r'$\mathrm{Log}(M_{\mathrm{L}}/M_{\odot})$',fontsize=40)
ax.set_ylabel(r'$\chi^2$',fontsize=40)
#ax.set_ylim([2000608,2000610])
#ax.set_yticks([2000608,2000609,2000610])
ax.tick_params(labelsize=25)
ax.legend(fontsize=25,frameon=False)

plt.show()
