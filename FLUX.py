import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
plt.style.use('seaborn-ticks')
mpl.rcParams['mathtext.default'] = 'regular'
font = {'family': 'Times New Roman',
         'style': 'normal',
        'weight': 'normal',
         'color': 'black', 
          'size': 18,
        }
my_path = os.path.abspath('')
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
a = np.loadtxt('flux')
t = a[:,0]
pfe = a[:,1]
pfi = a[:,2]
pfc = a[:,3]
hfe = a[:,4]
hfi = a[:,5]
hfc = a[:,6]
dum = 191616766./191700
gamma = dum*(np.log(hfi[-1]) - np.log(hfi[len(t)//2])) / (t[-1] - t[len(t)//2]) * .5
print('growth rate = %f'%gamma)
ax.semilogy(t,abs(hfi),label=r'$ion heat flux$')
ax.set_xlabel(r'$time$',fontdict = font)
ax.set_ylabel(r'',fontdict = font)
ax.legend(fontsize=13)
ax.tick_params(direction='in',labelsize=13)
plt.show()
