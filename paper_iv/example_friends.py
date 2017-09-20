
# <markdowncell>
# ## Shows the plotting tools.

# <markdowncell> Import teneto, numpy and matplotlib
# <codecell>

import teneto
import numpy as np
import matplotlib.pyplot as plt


# <markdowncell> Set color sceme
# <codecell>

plt.rcParams['image.cmap'] = 'gist_gray'

# <markdowncell> Create a 3D network
# <codecell>

A=np.zeros((3,3,3))
A[0,1,:]=1
A[1,0,:]=1
A[0,2,1:]=1
A[2,0,1:]=1
A[1,2,2]=1
A[2,1,2]=1

# <markdowncell> Create a 3D network
# <codecell>

fig, ax = plt.subplots(1,2)


# <markdowncell> Plot circle graph
# <codecell>

ax[0] = teneto.plot.circle_plot(A[:,:,1],ax[0])
ax[0].set_title('A',loc='left')

# <markdowncell> Plot slice graphlet
# <codecell>

ax[1] = teneto.plot.slice_plot(A,ax[1],['Ashley','Blake','Casey','Dylan'],['2014','2015','2016'])
ax[1].set_xlabel('time (years)')
ax[1].set_title('B',loc='left')

# <markdowncell> save figures
# <codecell>
fig.tight_layout()

fig.savefig('./examples/figures/friendexampletst.pdf')
fig.savefig('./examples/figures/friendexample.eps')
