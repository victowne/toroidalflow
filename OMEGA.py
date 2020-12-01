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
          'size': 18,
        }
my_path = os.path.abspath('')
nr = 7
nz = 64
nf = 640
vf = np.zeros((nf,1))
amp = np.zeros((nf,nr*nz))
tmp = np.zeros((3,nf))
f = open('freq')
lines = f.readlines()
pt = 0
for i in range(nr*nz):
    pt = pt + 6
    if i == 0:
        for j in range(nf):
            vf[j,0] = lines[pt].split()[1]
            amp[j,i] = lines[pt].split()[2]
            pt = pt + 1
    else:
        for j in range(nf):
            amp[j,i] = lines[pt].split()[2]
            pt = pt +1
    pt = pt + 1

i = np.argmax(amp[:,nr*nz-1])
print('frequency = %f'%vf[i,0])

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(vf,amp[:,:])
ax.set_xlabel(r'$frequency$',fontdict = font)
ax.set_ylabel(r'amplitude',fontdict = font)
ax.legend(fontsize=13)
ax.tick_params(direction='in',labelsize=13)
plt.show()
