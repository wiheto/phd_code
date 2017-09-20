import scipy.io as sio
import scipy.stats as sps
import numpy as np
import tegrato
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import slice_graph as sg
import pickle
plt.rcParams['image.cmap'] = 'gist_gray'


netid=list(map(int,sio.loadmat('./data/networkassignment')['PowerNetClass']))
nodeOrder = np.argsort(netid)
netdat_eo=pd.DataFrame(columns=['v','dmn','fp','va','co','sa','sn','au','da'],index=range(0,46))
netdat_ec=pd.DataFrame(columns=['v','dmn','fp','va','co','sa','sn','au','da'],index=range(0,46))
netlab=pd.Series([7,5,8,11,12,3,9,1,4,-1,10],index=['V','DM','FP','VA','DA','CO','SA','SM','AU','U','Sub'])
plotorder=np.array([10,1,5,6,0,2,3,4,7,9,8])

GEeo=[]
for s in range(0,46):
    print('Calculating for subject: ' + str(s))
    dat=sio.loadmat('./data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c1.mat')['binGraph']
    GEeo.append(tegrato.efficiency(dat))
pickle.dump(GEeo,open('./data/geff_eo.pickle','wb'))
[ ax[0].plot([0,1],[Ravg_eo[n],Ravg_ec[n]],'.-',color=(0.5,0.5,0.5),markerfacecolor=(0.2,0.2,0.2),linewidth=0.5,markersize=0) for n in range(0,len(Rec))]

GEec=[]
for s in range(0,46):
    print('Calculating for subject: ' + str(s))
    dat=sio.loadmat('./data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c2.mat')['binGraph']
    GEec.append(tegrato.efficiency(dat))

pickle.dump(GEec,open('./data/geff_ec.pickle','wb'))

#Load reachability data
Reo=pickle.load(open('./data/reach_eo.pickle','rb'))
Rec=pickle.load(open('./data/reach_ec.pickle','rb'))
Ravg_eo = np.array([Reo[s].groupby('time').mean().deltat.mean() for s in range(0,len(Reo))])
Ravg_ec = np.array([Rec[s].groupby('time').mean().deltat.mean() for s in range(0,len(Rec))])



f, ax = plt.subplots(1,3,figsize=(12,4))
sns.violinplot(data=[GEeo,GEec], bw=.2, cut=1, linewidth=2,ax=ax[0],color=(0.3,0.3,0.3))
sns.swarmplot(data=[GEeo,GEec], ax=ax[0],color="w",alpha=0.5)
[ ax[0].plot([0,1],[GEeo[n],GEec[n]],'.-',color=(0.5,0.5,0.5),markerfacecolor=(0.2,0.2,0.2),linewidth=0.5,markersize=0) for n in range(0,len(GEec))]
ax[0].set_xlim(-0.5,1.5)
ax[0].annotate('**',(0.47,36),size=20)
ax[0].set_xticks([0,1])
ax[0].set_xticklabels(['EO','EC'])
ax[0].set_ylabel('Global reachability (time)')

ax[1].scatter(Ravg_eo,GEeo,color='k',s=30,alpha=.5)
ax[1].set_xlabel('Global Reachability (EO)')
ax[1].set_ylabel('Temporal Efficiency (EO)')
ax[2].scatter(Ravg_ec,GEec,color='k',s=30,alpha=.5)
ax[2].set_xlabel('Global Reachability (EC)')
ax[2].set_ylabel('Temporal Efficiency (EC)')

r1=sps.spearmanr(Ravg_eo,GEeo)
r2=sps.spearmanr(Ravg_ec,GEec)
f.tight_layout()

f.savefig('./figures/fmri_geff.png')
f.savefig('./figures/fmri_geff.eps')


p_geff,pdist_geff=tegrato.stat_shuffle(np.array(GEeo),np.array(GEec),100000)
