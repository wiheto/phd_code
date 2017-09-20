
import scipy.io as sio
import scipy.stats as sps
import numpy as np
import teneto
import matplotlib.pyplot as plt
import pandas as pd
#import teneto.stats.shufflegroups as shuffle
plt.rcParams['image.cmap'] = 'gist_gray'



Eeo=np.zeros([46])
Eec=np.zeros([46])
for s in range(0,46):
    print('Calculating for subject: ' + str(s))
    dat=sio.loadmat('./data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c1.mat')['binGraph']
    dat[dat>0]=1
    Eeo[s]=teneto.networkmeasures.temporal_efficiency(dat)
    dat=sio.loadmat('./data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c2.mat')['binGraph']
    dat[dat>0]=1
    Eec[s]=teneto.networkmeasures.temporal_efficiency(dat)

CEeo=np.zeros([46,264])
CEec=np.zeros([46,264])
for s in range(27,46):
    print('Calculating for subject: ' + str(s))
    dat=sio.loadmat('./data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c1.mat')['binGraph']
    dat[dat>0]=1
    CEeo[s,:]=teneto.networkmeasures.temporal_efficiency(dat,do='node')
    dat=sio.loadmat('./data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c2.mat')['binGraph']
    dat[dat>0]=1
    CEec[s,:]=teneto.networkmeasures.temporal_efficiency(dat,do='node')

np.save('./data/efficiency_centrality_eo.npy',CEeo)
np.save('./data/efficiency_centrality_ec.npy',CEec)



np.save('./data/efficiency_eo.npy',Eeo)
np.save('./data/efficiency_ec.npy',Eec)

Eeo=np.load('./data/efficiency_eo.npy')
Eec=np.load('./data/efficiency_ec.npy')

Reo=np.load('./data/reachability_eo.npy')
Rec=np.load('./data/reachability_ec.npy')

CEeo=np.load('./data/efficiency_centrality_eo.npy')
CEec=np.load('./data/efficiency_centrality_ec.npy')


GReo = np.nanmean(Reo,axis=0)
GRec = np.nanmean(Rec,axis=0)


np.random.seed(2016)
preach,pdist = shuffle.shufflegroups(Eeo,Eec,100000)


fig, ax = plt.subplots(1,3)


sns.violinplot(data=[Eeo,Eec], bw=.2, cut=1, linewidth=2,ax=ax[0],color=(0.3,0.3,0.3))
sns.swarmplot(data=[Eeo,Eec], ax=ax[0],color="w",alpha=0.5)

[ ax[0].plot([0,1],[Eeo[n],Eec[n]],'.-',color=(0.5,0.5,0.5),markerfacecolor=(0.2,0.2,0.2),linewidth=0.5,markersize=0) for n in range(0,len(Eec))]






ax[0].set_xlim(-0.5,1.5)
ax[0].annotate('**',(0.45,.22),size=20)
ax[0].set_xticks([0,1])
ax[0].set_xticklabels(['EO','EC'])
ax[0].set_ylabel('Temporal Efficiency')
ax[0].set_title('A',loc='left')
ax[0].set_yticks(np.arange(0.05,0.251,0.05))
ax[0].set_ylim(0.075,0.225)


ax[1].scatter(Eeo,GReo,color='k',s=30,alpha=0.5)
ax[1].set_ylabel('Temporal Efficiency (EO)')
ax[1].set_xlabel('Reachability Latency (EO)')
ax[1].set_title('B',loc='left')
ax[1].annotate('***',(.15,33),size=20)
ax[1].set_xticks(np.arange(0.05,0.251,0.05))
ax[1].set_xlim(0.075,0.225)

ax[2].scatter(Eec,GRec,color='k',s=30,alpha=0.5)
ax[2].set_ylabel('Temporal Efficiency (EC)')
ax[2].set_xlabel('Reachability Latency (EC)')
ax[2].set_title('C',loc='left')
ax[2].annotate('***',(.15,33),size=20)
ax[2].set_xticks(np.arange(0.05,0.251,0.05))
ax[2].set_xlim(0.075,0.225)

yth
sp_eo=sps.spearmanr(Eeo,GReo)
sp_ec=sps.spearmanr(Eec,GRec)

fig.tight_layout()

fig.savefig('./examples/figures/eff_fmri.pdf')
