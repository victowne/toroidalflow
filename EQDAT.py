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

f = open('eqdat')
lines = f.readlines()
eqdat = np.loadtxt(lines[15:-1])
i = eqdat[:,0]
q = eqdat[:,1]
ni = eqdat[:,2]
ne = eqdat[:,3]
nc = eqdat[:,4]
ti = eqdat[:,5]
tc = eqdat[:,6]
capni = eqdat[:,7]
capnc = eqdat[:,8]
captc = eqdat[:,9]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(i,q)
ax.set_xlabel(r'$t/\tau_a$',fontdict = font)
ax.set_ylabel(r'island width',fontdict = font)
ax.legend(fontsize=13)
ax.tick_params(direction='in',labelsize=13)
plt.show()
