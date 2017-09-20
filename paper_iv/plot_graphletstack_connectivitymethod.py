
# <markdowncell> Import scipy.io for matlab loading, numpy and teneto
# <codecell>
import numpy as np
import teneto
import matplotlib.pyplot as plt
import scipy.io as sio
import pandas as pd

dat=sio.matlab.loadmat('./additionalsfrommatlab/graphletplotdat.mat')['save4tenetoplot']

for n in range(0,264):
    dat[n,n,:]=0

fig,ax=plt.subplots(1,figsize=(20,10))

ax=teneto.plot.graphlet_stack_plot(np.array(dat[:,:,9:100:10]),ax,t0=1,Fs=10,q=20,cmap='seismic',gridcolor='none',sharpen='yes',borderwidth=5)
fig.savefig('./figures/graphletexample_connectivitymethod.pdf',r=600)




netid=list(map(int,sio.loadmat('./data/networkassignment')['PowerNetClass']))
nodeOrder = np.argsort(netid)
netdat_eo=pd.DataFrame(columns=['v','dmn','fp','va','co','sa','sn','au','da'],index=range(0,46))
netdat_ec=pd.DataFrame(columns=['v','dmn','fp','va','co','sa','sn','au','da'],index=range(0,46))
netlab=pd.Series([7,5,8,11,12,3,9,1,4],index=['v','dmn','fp','va','da','co','sa','sn','au'])

dat=sio.loadmat('./data/bingraph_weightcorr_s10_c1.mat')['binGraph']
plt.rcParams['image.cmap'] = 'inferno_r'

v=np.where(netid==netlab['v'])[0]
vdatslice=dat[v,:,:][:,v,:][:,:,80:160]
nodeNames = []
for n in range(0,len(v)):
    nodeNames.append(str(n+1))

fig,ax=plt.subplots(1,1,figsize=(20,8))
ax = teneto.plot.slice_plot(vdatslice,ax,nLabs='',tLabs='',timeunit='',linestyle='k-',nodesize=20,nodealpha=.7)
ax.set_xlabel('time')
ax.set_ylabel('Visual network ROIs')
ax.set_xticks(range(9,80,10))
ax.set_xticklabels(range(90,160,10))
fig.savefig('./figures/visualnet_slice_graph_example_s10.pdf')
fig.savefig('./figures/visualnet_slice_graph_example_s10.png')
fig.show()
