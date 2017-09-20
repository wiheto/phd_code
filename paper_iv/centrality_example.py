
import numpy as np
import teneto
import matplotlib.pyplot as plt
plt.rcParams['image.cmap'] = 'gist_gray'

A=np.zeros((3,3,20))
A[0,2,0:4]=1
A[0,1,0]=1
A[0,1,5]=1
A[0,1,10]=1
A[0,1,15]=1


fig,ax = plt.subplots(1)


ax = teneto.plot.slice_plot(A,ax,vlabs=range(1,4),dlabs=range(1,21))
ax.set_ylabel('nodes')
ax.set_xlabel('time')
ax.set_ylim(-0.25,2.25)

fig.tight_layout()
fig.show()

fig.savefig('./examples/figures/centrality_examples.pdf')
