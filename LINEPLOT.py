import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
plt.style.use('seaborn-ticks')
mpl.rcParams['mathtext.default'] = 'regular'
font = {'family': 'Times New Roman',
         'style': 'normal',
        'weight': 'normal',
         'color': 'black', 
          'size': 24,
        }
my_path = os.path.abspath('')
data = np.loadtxt(my_path+'/out/pol1d.out')
R = data[:,0]
phi = data[:,1]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(R,phi,label=r'$\phi$')
ax.set_xlabel(r'$r$',fontdict = font)
ax.set_ylabel(r'$\phi$',fontdict = font)
ax.legend(fontsize=17)
ax.tick_params(direction='in',labelsize=17)
plt.show()
