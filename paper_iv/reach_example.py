import teneto
import numpy as np
import teneto
import matplotlib.pyplot as plt
plt.rcParams['image.cmap'] = 'gist_gray'

A=np.load('./examples/data/reach_example.npy')
A=A[:,:,:8]

#Create subplots
fig,ax=plt.subplots(3,1,figsize=(5,8))

#Plot 0: overview

ax[0]=teneto.plot.slice_plot(A,ax[0],dlabs=range(1,9),vlabs=range(1,6))

ax[0].set_xlabel('time')
ax[0].set_ylabel('nodes')
ax[0].set_title('A',loc='left')


#Plot 1: Longest shortest path from 0

ax[1]=teneto.plot.slice_plot(A,ax[1],dlabs=range(1,9),vlabs=range(1,6))

path=[(0,0),(0,1),(2,1),(2,2),(4,2),(4,3),(5,3),(5,4)]
path=list(zip(*path))
for n in range(0,len(path)-1):
    ax[1].plot(path[n],path[n+1],'#8c000f',linewidth=3)

ax[1].set_xlabel('time')
ax[1].set_ylabel('nodes')
ax[1].set_title('B',loc='left')


#Plot 2: vertical movement illustration

ax[2]=teneto.plot.slice_plot(A,ax[2],dlabs=range(1,9),vlabs=range(1,6))

path=[(1,4),(1,2),(2,2),(2,1)]
path=list(zip(*path))
for n in range(0,len(path)-1):
    ax[2].plot(path[n],path[n+1],'#8c000f',linewidth=3)

path=[(1,4),(6,4),(6,1)]
path=list(zip(*path))
for n in range(0,len(path)-1):
    ax[2].plot(path[n],path[n+1],'#0a437a',linewidth=3)

ax[2].set_xlabel('time')
ax[2].set_ylabel('nodes')
ax[2].set_title('C',loc='left')

fig.tight_layout()

fig.savefig('./examples/figures/shortestpath_example.pdf')
fig.show()
