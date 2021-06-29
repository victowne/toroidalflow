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
rin = 0.2
rout = 0.8
lref = 706.71
r = np.linspace(rin,rout,301)

fig = plt.figure()
ax = fig.add_subplot(1,2,1)
ax.plot(r,ti)
ax.set_xlabel(r'$r$',fontdict = font)
ax.set_ylabel(r'$T_i$',fontdict = font)
ax.legend(fontsize=13)
ax.tick_params(direction='in',labelsize=13)
ax.set_xlim((rin,rout))
ax = fig.add_subplot(1,2,2)
ax.plot(r,capti*lref)
ax.set_xlabel(r'$r$',fontdict = font)
ax.set_ylabel(r'$L_{ref}/L_{Ti}$',fontdict = font)
ax.legend(fontsize=13)
ax.tick_params(direction='in',labelsize=13)
ax.set_xlim((rin,rout))
plt.show()
