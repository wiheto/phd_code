import scipy.io as sio
import numpy as np
import teneto
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams['image.cmap'] = 'gist_gray'

netid=list(map(int,sio.loadmat('./data/networkassignment')['PowerNetClass']))
nodeOrder = np.argsort(netid)
netdat_eo=pd.DataFrame(columns=['v','dmn','fp','va','co','sa','sn','au','da'],index=range(0,46))
netdat_ec=pd.DataFrame(columns=['v','dmn','fp','va','co','sa','sn','au','da'],index=range(0,46))
netlab=pd.Series([7,5,8,11,12,3,9,1,4],index=['v','dmn','fp','va','da','co','sa','sn','au'])

np.where(netid==netlab['dmn'])[0]

dat=sio.loadmat('./examples/data/bingraph_weightcorr_s30_c1.mat')['binGraph']
v=np.where(netid==netlab['v'])[0]
vdatslice=dat[v,:,:][:,v,:][:,:,0:100]
nodeNames = []
for n in range(0,len(v)):
    nodeNames.append(str(n+1))

fig,ax=plt.subplots(1,1,figsize=(20,8))
ax = teneto.plot.slice_plot(vdatslice,ax,nodeNames,range(0,100),nodesize=20)
ax.set_xlabel('time')
ax.set_ylabel('Visual network ROIs')
ax.set_xticks(range(9,100,10))
ax.set_xticklabels(range(10,101,10))
fig.show()
fig.savefig('./examples/figures/visualnet_slice_graph_example.pdf')
fig.savefig('./examples/figures/visualnet_slice_graph_example.png')



dat=sio.loadmat('./data/bingraph_weightcorr_s30_c1.mat')['binGraph']
teneto.networkmeasures.temporal_degree
