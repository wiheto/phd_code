
# <markdowncell> Import scipy.io for matlab loading, numpy and teneto
# <codecell>
import scipy.io as sio
import scipy.stats as sps
import numpy as np
import teneto
import matplotlib.pyplot as plt
import pandas as pd
import os






# <markdowncell> Set matplotlib color style
# <codecell>
plt.rcParams['image.cmap'] = 'gist_gray'

# <markdowncell> load data
# <codecell>
closeness_eo=np.load('./data/closeness_eo.npy')
closeness_ec=np.load('./data/closeness_ec.npy')
degree_eo=np.load('./data/degree_eo.npy')
degree_ec=np.load('./data/degree_ec.npy')

# <markdowncell> make data <subject x node>
# <codecell>
closeness_eo = np.stack(closeness_eo)
closeness_ec = np.stack(closeness_ec)
degree_eo = np.stack(degree_eo)
degree_ec = np.stack(degree_ec)

# <markdowncell> define figure
# <codecell>
fig, ax = plt.subplots(2,2)


# <markdowncell> Plot temporal degree eo vs ec
# <codecell>
ax[0,0].scatter(np.mean(degree_eo,axis=0),np.mean(degree_ec,axis=0),color='k',s=30,alpha=0.5)
ax[0,0].set_xlabel('Temporal degree: EO')
ax[0,0].set_ylabel('Temporal degree: EC')
ax[0,0].set_title('A',loc='left')
ax[0,0].set_yticks(np.arange(1800,1940.1,40))
ax[0,0].set_ylim(1800,1940)
ax[0,0].set_xticks(np.arange(1800,1940.1,40))
ax[0,0].set_xlim(1800,1940)
x0,x1 = ax[0,0].get_xlim()
y0,y1 = ax[0,0].get_ylim()
ax[0,0].set_aspect((x1-x0)/(y1-y0))

# <markdowncell> Plot temporal closeness eo vs ec
# <codecell>
ax[0,1].scatter(np.mean(closeness_eo,axis=0),np.mean(closeness_ec,axis=0),color='k',s=30,alpha=0.5)
ax[0,1].set_xlabel('Temporal closeness: EO')
ax[0,1].set_ylabel('Temporal closeness: EC')
ax[0,1].set_title('B',loc='left')
ax[0,1].set_yticks(np.arange(0.14,0.181,0.01))
ax[0,1].set_ylim(0.14,0.18)
ax[0,1].set_xticks(np.arange(0.14,0.181,0.01))
ax[0,1].set_xlim(0.14,0.18)
x0,x1 = ax[0,1].get_xlim()
y0,y1 = ax[0,1].get_ylim()
ax[0,1].set_aspect((x1-x0)/(y1-y0))

# <markdowncell> Plot temporal degree vs temporal  closeness for eo condition
# <codecell>
ax[1,0].scatter(np.mean(degree_eo,axis=0),np.mean(closeness_eo,axis=0),color='k',s=30,alpha=0.5)
ax[1,0].set_xlabel('Temporal degree: EO')
ax[1,0].set_ylabel('Temporal closeness: EO')
ax[1,0].set_title('C',loc='left')
ax[1,0].set_xticks(np.arange(1800,1940.1,40))
ax[1,0].set_xlim(1800,1940)
ax[1,0].set_yticks(np.arange(0.14,0.181,0.01))
ax[1,0].set_ylim(0.14,0.18)
x0,x1 = ax[1,0].get_xlim()
y0,y1 = ax[1,0].get_ylim()
ax[1,0].set_aspect((x1-x0)/(y1-y0))


# <markdowncell> Plot temporal degree vs temporal  closeness for eo condition
# <codecell>
ax[1,1].scatter(np.mean(degree_ec,axis=0),np.mean(closeness_ec,axis=0),color='k',s=30,alpha=0.5)
ax[1,1].set_xlabel('Temporal degree: EC')
ax[1,1].set_ylabel('Temporal closeness: EC')
ax[1,1].set_title('D',loc='left')
ax[1,1].set_xticks(np.arange(1800,1940.1,40))
ax[1,1].set_xlim(1800,1940)
ax[1,1].set_yticks(np.arange(0.14,0.181,0.01))
ax[1,1].set_ylim(0.14,0.18)
x0,x1 = ax[1,1].get_xlim()
y0,y1 = ax[1,1].get_ylim()
ax[1,1].set_aspect((x1-x0)/(y1-y0))

# <markdowncell> make the figure tight and save
# <codecell>
fig.tight_layout()
fig.savefig('./examples/figures/TDvsTC.pdf')


# <markdowncell> Get spearman values for each of the subplots
# <codecell>
spear_td=sps.spearmanr(np.mean(degree_eo,axis=0),np.mean(degree_ec,axis=0))
spear_tc=sps.spearmanr(np.mean(closeness_eo,axis=0),np.mean(closeness_ec,axis=0))
spear_tdtceo=sps.spearmanr(np.mean(closeness_eo,axis=0),np.mean(degree_eo,axis=0))
spear_tdtcec=sps.spearmanr(np.mean(closeness_ec,axis=0),np.mean(degree_ec,axis=0))
