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
a = np.loadtxt('yyre')
t = a[:,0]
v1 = a[:,1]
v2 = a[:,2]
v3 = a[:,3]
v4 = a[:,4]
v5 = a[:,5]
ax.plot(t,v3,label=r'$f_b=0.0$')
ax.set_xlabel(r'$t/\tau_a$',fontdict = font)
ax.set_ylabel(r'island width',fontdict = font)
ax.legend(fontsize=13)
ax.tick_params(direction='in',labelsize=13)
plt.show()
