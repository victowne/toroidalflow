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
data = np.loadtxt(my_path+'/out/polxz.out')
bdy = np.loadtxt(my_path+'/out/bd.out')
dum=bdy.shape[0]/2
bdir = bdy[:dum,0]
bdiz = bdy[:dum,1]
bdor = bdy[dum:,0]
bdoz = bdy[dum:,1]
R = data[:,1].reshape(401,601,order='C')
Z = data[:,2].reshape(401,601,order='C')
polxz = data[:,3].reshape(401,601,order='C')

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.contourf(R,Z,polxz,100,cmap='PiYG')
ax.plot(bdor,bdoz,'k-')
ax.set_xlabel(r'$R/\rho_i$',fontdict = font)
ax.set_ylabel(r'$Z/\rho_i$',fontdict = font)
ax.legend(fontsize=17)
ax.tick_params(direction='in',labelsize=17)
plt.show()
