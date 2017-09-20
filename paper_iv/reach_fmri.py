
import scipy.io as sio
import scipy.stats as sps
import numpy as np
import teneto
import matplotlib.pyplot as plt
import pandas as pd
import teneto.stats.shufflegroups as shuffle
plt.rcParams['image.cmap'] = 'gist_gray'



Reo=np.zeros([264,46])
Rec=np.zeros([264,46])
for s in range(0,46):
    print('Calculating for subject: ' + str(s))
    dat=sio.loadmat('./examples/data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c1.mat')['binGraph']
    dat[dat>0]=1
    Reo[:,s]=teneto.reachabilityLatency(dat,1,'nodes')
    dat=sio.loadmat('./examples/data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c2.mat')['binGraph']
    dat[dat>0]=1
    Rec[:,s]=teneto.reachabilityLatency(dat,1,'nodes')


np.save('./examples/data/reachability_eo.npy',Reo)
np.save('./examples/data/reachability_ec.npy',Rec)

Reo=np.load('./examples/data/reachability_eo.npy')
Rec=np.load('./examples/data/reachability_ec.npy')

GReo = np.nanmean(Reo,axis=0)
GRec = np.nanmean(Rec,axis=0)

netid=list(map(int,sio.loadmat('./examples/data/networkassignment')['PowerNetClass']))
network = np.array([1,3,4,5,7,8,9,10,11,12])
netlab = np.array(['SM','CO','AU','DM','V','FP','SA','Sub','VA','DA'])

np.random.seed(2016)
preach,pdist = shuffle.shufflegroups(GReo,GRec,100000)

NetReo = np.zeros(13)
NetRec = np.zeros(13)
sNetRec = np.zeros([13,46])
sNetReo = np.zeros([13,46])
for n in set(netid):
    if n!=-1:
        fid = np.where(np.array(netid)==n)
        NetReo[n] = np.nanmean(np.nanmean(np.squeeze(Reo[fid,:]),axis=0))
        NetRec[n] = np.nanmean(np.nanmean(np.squeeze(Rec[fid,:]),axis=0))

Rnetdif=NetRec-NetReo
Rnetdif=Rnetdif[network]

odr = np.argsort(Rnetdif)[::-1]
netodr = netlab[odr]


fig, ax = plt.subplots(1,3)


sns.violinplot(data=[GReo,GRec], bw=.2, cut=1, linewidth=2,ax=ax[0],color=(0.3,0.3,0.3))
sns.swarmplot(data=[GReo,GRec], ax=ax[0],color="w",alpha=0.5)

[ ax[0].plot([0,1],[GReo[n],GRec[n]],'.-',color=(0.5,0.5,0.5),markerfacecolor=(0.2,0.2,0.2),linewidth=0.5,markersize=0) for n in range(0,len(GRec))]

ax[0].set_xlim(-0.5,1.5)
ax[0].annotate('***',(0.47,34),size=20)
ax[0].set_xticks([0,1])
ax[0].set_xticklabels(['EO','EC'])
ax[0].set_ylabel('Global reachability (time)')
ax[0].set_title('A',loc='left')

ax[1].scatter(range(0,10),Rnetdif[odr],color='k',s=100,alpha=0.5)
ax[1].set_xticks(np.arange(0,10))
ax[1].set_xlim(-0.5,9.5)
ax[1].set_xticklabels(netodr)
ax[1].set_ylabel('Reachability difference per network (EC>EO)')
ax[1].set_xlabel('Network')
ax[1].set_title('B',loc='left')


ax[2].scatter(GReo,GRec,color='k',s=30,alpha=0.5)
ax[2].set_ylabel('Global Reachability EC')
ax[2].set_xlabel('Global Reachability EO')
ax[2].set_title('C',loc='left')
ax[2].annotate('*',(26,33),size=20)

for n in range(0,3):
    x0,x1 = ax[n].get_xlim()
    y0,y1 = ax[n].get_ylim()
    ax[n].set_aspect((x1-x0)/(y1-y0))

sp_avg=sps.spearmanr(GReo,GRec)

fig.tight_layout()

fig.savefig('./examples/figures/reach_fmri.pdf')
