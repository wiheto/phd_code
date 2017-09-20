import scipy.stats as sps
import numpy as np
import teneto
import teneto.stats.shufflegroups as shuffle
import matplotlib.pyplot as plt
import seaborn as sns
import math
plt.rcParams['image.cmap'] = 'gist_gray'


#Plot

vol_ec=np.load('./examples/data/vol_ec.npy')
vol_eo=np.load('./examples/data/vol_eo.npy')
fluct_ec=np.load('./examples/data/fluct_ec.npy')
fluct_eo=np.load('./examples/data/fluct_eo.npy')


np.random.seed(2016)
p_vol,pdist_vol = shuffle.shufflegroups(vol_eo,vol_ec,100000,'subjects',2)
p_fluct,pdist_fluct = shuffle.shufflegroups(fluct_eo,fluct_ec,100000,'subjects',2)

bsp_fluct=sps.spearmanr(fluct_eo,fluct_ec)
bsp_vol=sps.spearmanr(vol_eo,vol_ec)


fig, ax = plt.subplots(2, 2)

ax[0,0].scatter(fluct_eo,fluct_ec,color='k',s=30,alpha=0.5)
ax[0,0].set_xlim(0.125,0.155)
ax[0,0].set_ylim(0.125,0.155)
ax[0,0].set_xlabel('Fluctuability (EO)')
ax[0,0].set_ylabel('Fluctuability (EC)')
ax[0,0].annotate('NS',(0.14,0.155),size=12)
ax[0,0].set_title('A',loc='left')

ax[1,0].scatter(vol_eo,vol_ec,color='k',s=30,alpha=0.5)
ax[1,0].set_xlabel('Volatility (EO)')
ax[1,0].set_ylabel('Volatility (EC)')
ax[1,0].annotate('**',(885,1075),size=20)
ax[1,0].set_title('C',loc='left')
ax[1,0].set_xlim(700,1100)
ax[1,0].set_ylim(700,1100)
ax[1,0].set_xticks(range(700,1101,100))
ax[1,0].set_yticks(range(700,1101,100))


sns.violinplot(data=[fluct_eo,fluct_ec], bw=.2, cut=1, linewidth=2,ax=ax[0,1],color=(0.3,0.3,0.3))
sns.swarmplot(data=[fluct_eo,fluct_ec], ax=ax[0,1],color="w",alpha=0.5)
[ ax[0,1].plot([0,1],[fluct_eo[n],fluct_ec[n]],'.-',color=(0.5,0.5,0.5),markerfacecolor=(0.2,0.2,0.2),linewidth=0.5,markersize=0) for n in range(0,len(fluct_eo))]


ax[0,1].set_xlim(-0.5,1.5)
ax[0,1].annotate('**',(0.4,0.155),size=20)
ax[0,1].set_xticks([0,1])
ax[0,1].set_yticks(np.arange(0.125,0.1551,0.005))
ax[0,1].set_ylim(0.125,0.155)
ax[0,1].set_xticklabels(['EO','EC'])
ax[0,1].set_ylabel('Fluctuability')
ax[0,1].set_title('B',loc='left')


sns.violinplot(data=[vol_eo,vol_ec], bw=.2, cut=1, linewidth=2,ax=ax[1,1],color=(0.3,0.3,0.3))
sns.swarmplot(data=[vol_eo,vol_ec], ax=ax[1,1],color="w",alpha=0.5)
[ ax[1,1].plot([0,1],[vol_eo[n],vol_ec[n]],'.-',color=(0.5,0.5,0.5),markerfacecolor=(0.2,0.2,0.2),linewidth=0.5,markersize=0) for n in range(0,len(vol_ec))]


ax[1,1].annotate('NS',(0.45,1100),size=12)
ax[1,1].set_xlim(-0.5,1.5)
ax[1,1].set_xticks([0,1])
ax[1,1].set_xticklabels(['EO','EC'])
ax[1,1].set_ylabel('Volatility')
ax[1,1].set_title('D',loc='left')

for n in range(0,2):
    for m in range(0,2):
        x0,x1 = ax[n,m].get_xlim()
        y0,y1 = ax[n,m].get_ylim()
        ax[n,m].set_aspect((x1-x0)/(y1-y0))

fig.tight_layout()

fig.savefig('./examples/figures/fmri_fluctvol.pdf')

#for n in range(0,9):
#    p_net[n],pdist2 = tegrato.stat_shuffle(vol_eo_n[:,n],vol_ec_n[:,n],1000,'subjects')
