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
data = np.loadtxt(my_path+'/out/mpol.out')
nm = 25
nr = 5001
r = data[:,1].reshape(nr,nm,order='F')
val = data[:,2].reshape(nr,nm,order='F')
valr = data[:,3].reshape(nr,nm,order='F')
vali = data[:,4].reshape(nr,nm,order='F')

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for i in range(0,nm):
    ax.plot(r[:,i],val[:,i],label='m='+str(i))
    k = np.where(val[:,i] == val[:,i].max())[0][0]
    ax.text(r[k,i],val[k,i],'m='+str(i))
#ax.plot(r[:,22],valr[:,22])
#ax.plot(r[:,22],vali[:,22],'--')
ax.set_xlabel(r'$r$',fontdict = font)
ax.set_ylabel(r'$a. u.$',fontdict = font)
#ax.legend(fontsize=17)
ax.tick_params(direction='in',labelsize=17)
plt.show()
