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
ne = eqdat[:,2]
ni = eqdat[:,3]
nc = eqdat[:,4]
te = eqdat[:,5]
ti = eqdat[:,6]
tc = eqdat[:,7]
capne = eqdat[:,8]
capni = eqdat[:,9]
capnc = eqdat[:,10]
capte = eqdat[:,11]
capti = eqdat[:,12]
captc = eqdat[:,13]
omg = eqdat[:,14]
domg = eqdat[:,15]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(i,q)
ax.set_xlabel(r'$t/\tau_a$',fontdict = font)
ax.set_ylabel(r'island width',fontdict = font)
ax.legend(fontsize=13)
ax.tick_params(direction='in',labelsize=13)
plt.show()
