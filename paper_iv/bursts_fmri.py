# <markdowncell>
# Here the intercontact times per edge are concenated, bursty coefficient is calculated and plotted

# <markdowncell> Import scipy.io for matlab loading, scipy.stats for spearman, numpy and teneto, matplotlib and pandas.
# <codecell>
import scipy.io as sio
import scipy.stats as sps
import numpy as np
import teneto
import matplotlib.pyplot as plt
import pandas as pd

# <markdowncell>
### Part 1, data needs to be loaded and concatenated
# <codecell>
#Load intercontact times
ict_eo = np.load('./examples/data/ict_eo.npy')
ict_ec = np.load('./examples/data/ict_ec.npy')


# <markdowncell>
#Concatenate over subjects (a little clunky but works)
# <codecell>
for s in range(0,46):
    if s == 0:
        ict_ec_psub = np.triu(ict_ec[s]['intercontacttimes'],k=1)
        ict_eo_psub = np.triu(ict_eo[s]['intercontacttimes'],k=1)
    else:
        ict_ec_psub=np.dstack([ict_ec_psub,np.triu(ict_ec[s]['intercontacttimes'],k=1)])
        ict_eo_psub=np.dstack([ict_eo_psub,np.triu(ict_eo[s]['intercontacttimes'],k=1)])


# <markdowncell>
#For example (panel A) where all edges are concatenated in eyes open concatenated over all subjects
# <codecell>
triind=np.triu_indices(264,k=1)
ict_eo_concat=np.concatenate(np.concatenate(ict_eo_psub[triind[0],triind[1]]))
# Calculate distribution of ict_eo_comb (i.e. all nodes concatenated)
hist,histlabs,x=plt.hist(ict_eo_concat,bins=range(0,241),normed=True)


# <markdowncell>
# Example edge. Calculates bursty coeffieicnt for each subject and adding ICTs of that subject to those already subject's already processed (and recalculating)
# <codecell>
#Example nodes
i=101 #rmfc
j=94 #rpcc
B_sub_concat = np.zeros(46)
B_sub = np.zeros(46)
#Create dummy dictioanry
ict_concat_dict={}
for s in range(0,46):
    ictoi=ict_eo[s]
    if s == 0:
        ict_concat_dict['intercontacttimes']=ictoi['intercontacttimes'][i,j]
    else:
        ict_concat_dict['intercontacttimes']=np.concatenate([ict_concat_dict['intercontacttimes'],ictoi['intercontacttimes'][i,j]])
    B_sub_concat[s]=teneto.burstycoeff(ict_concat_dict)
    B_sub[s] = teneto.burstycoeff(ictoi,nodes=[i,j])[i,j]



# <markdowncell>
# This collects all icts of each subject into a pandas. Each row is an edge. Each column is a subject
# <codecell>

## Slightly convluted way to concatenate intercontacttimes across subjects. Will make more dedicated function for it
ictSub = pd.DataFrame(index=np.arange(0,len(triind[0])),columns=np.arange(0,46))
for s in range(0,46):
    ictSub[s]=ict_eo_psub[triind[0],triind[1],s]
icttmp={}
b=np.zeros(len(triind[0]))
for n in range(0,len(triind[0])):
    icttmp['intercontacttimes']=np.concatenate(ictSub.iloc[n])
    b[n]=teneto.burstycoeff(icttmp)
Beo = np.zeros([264,264])*np.nan
Beo[triind[0],triind[1]] = b

# <markdowncell>
# Do the same for eyes closed data
# <codecell>
ictSub = pd.DataFrame(index=np.arange(0,len(triind[0])),columns=np.arange(0,46))
for s in range(0,46):
    ictSub[s]=ict_ec_psub[triind[0],triind[1],s]
icttmp={}
b=np.zeros(len(triind[0]))

for n in range(0,len(triind[0])):
    icttmp['intercontacttimes']=np.concatenate(ictSub.iloc[n])
    b[n]=teneto.burstycoeff(icttmp)
Bec = np.zeros([264,264])*np.nan
Bec[triind[0],triind[1]] = b

# <markdowncell>
# Add transpose so entire matrix (not jsut upper tri) has been calculated
# <codecell>
Beo=np.nansum([Beo,Beo.transpose()],axis=0)
Bec=np.nansum([Bec,Bec.transpose()],axis=0)

# <markdowncell>
# save data
# <codecell>
np.save('./examples/data/burst_ec',Bec)
np.save('./examples/data/burst_eo',Beo)




# <markdowncell>
# get network information
# <codecell>
netid=list(map(int,sio.loadmat('./examples/data/networkassignment')['PowerNetClass']))
nodeOrder = np.argsort(netid)
netlab=pd.Series([7,5,8,11,12,3,9,1,4,-1,10],index=['V','DM','FP','VA','DA','CO','SA','SM','AU','U','Sub'])
netlab_10=netlab.drop('U')

# <markdowncell>
# average within and between network bursty coefficient
# <codecell>
Bwithin_eo=np.array([])
Bout_eo=np.array([])
Bwithin_ec=np.array([])
Bout_ec=np.array([])
for n in netlab_10:
    for m in netlab_10:
        ind1=np.where(np.array(netid)==n)[0]
        ind2=np.where(np.array(netid)==m)[0]
        if n==m:
            ind=np.triu_indices(len(ind1),k=1)
            Bwithin_ec=np.hstack([Bwithin_ec,Bec[ind1][:,ind1][ind]])
            Bwithin_eo=np.hstack([Bwithin_eo,Beo[ind1][:,ind1][ind]])
        else:
            Bout_ec=np.hstack([Bout_ec,Bec[ind1][:,ind2].reshape(len(ind1)*len(ind2))])
            Bout_eo=np.hstack([Bout_eo,Beo[ind1][:,ind2].reshape(len(ind1)*len(ind2))])


# <markdowncell>cs
# # PART 2 - plotting
# All one code for now
# <codecell>
fig,ax = plt.subplots(2,3)

ax[0,0].set_yscale("log")
ax[0,0].scatter(histlabs[1:],hist,color='k',s=30,alpha=0.5)
ax[0,0].set_ylim(1e-7,1e0)
ax[0,0].set_xlim(1,250)
ax[0,0].set_xlabel('Inter-contact times (volumes)')
ax[0,0].set_ylabel('P(τ)')

ax[0,0].set_title('A',loc='left')

ax[0,1].plot(range(0,46),B_sub_concat,color='k')
ax[0,1].scatter(range(0,46),B_sub,color='k',s=30,alpha=0.5)
ax[0,1].set_ylim(-1,1)
ax[0,1].set_xlim(0,46)
x0,x1 = ax[0,1].get_xlim()
y0,y1 = ax[0,1].get_ylim()
ax[0,1].set_aspect((x1-x0)/(y1-y0))
ax[0,1].set_ylabel('B')
ax[0,1].set_xlabel('Subject')
ax[0,1].set_title('B',loc='left')


ax[0,2].hist(Beo[triind[0],triind[1]],normed=True,bins=20,color='b',alpha=0.5)
ax[0,2].hist(Bec[triind[0],triind[1]],normed=True,bins=20,color='r',alpha=0.5)
ax[0,2].set_xticks(np.arange(0.1,0.41,0.1))

ax[0,2].set_xlim(0.1,0.4)
x0,x1 = ax[0,2].get_xlim()
y0,y1 = ax[0,2].get_ylim()
ax[0,2].set_aspect((x1-x0)/(y1-y0))
ax[0,2].set_xlabel('B')
ax[0,2].set_ylabel('Number of Edges (%)')
ax[0,2].set_title('C',loc='left')

val,bins,c=np.histogram2d(Beo[triind[0],triind[1]],Bec[triind[0],triind[1]],bins=np.arange(0.2,0.36,0.005))
img=ax[1,2].pcolormesh(bins,bins,val)
x0,x1 = ax[1,2].get_xlim()
y0,y1 = ax[1,2].get_ylim()
ax[1,2].set_yticks(np.arange(0.1,0.41,0.05))
ax[1,2].set_ylim(0.2,0.35)
ax[1,2].set_xticks(np.arange(0.1,0.41,0.05))
ax[1,2].set_xlim(0.2,0.35)
ax[1,2].set_aspect((x1-x0)/(y1-y0))
ax[1,2].set_xlabel('B(Eyes Open)')
ax[1,2].set_ylabel('B(Eyes Closed)')
ax[1,2].set_title('F',loc='left')
ax[1,2].set_frame_on(True)
plt.colorbar(img)

ax[1,0].hist(Bout_eo,normed=True,bins=20,color='b',alpha=0.5)
ax[1,0].hist(Bwithin_eo,normed=True,bins=20,color='r',alpha=0.5)

ax[1,0].set_xticks(np.arange(0.1,0.41,0.1))
ax[1,0].set_xlim(0.1,0.4)
x0,x1 = ax[1,0].get_xlim()
y0,y1 = ax[1,0].get_ylim()
ax[1,0].set_aspect((x1-x0)/(y1-y0))
ax[1,0].set_xlabel('B')
ax[1,0].set_ylabel('Number of Edges (%)')
ax[1,0].set_title('D',loc='left')


ax[1,1].hist(Bout_ec,normed=True,bins=20,color='b',alpha=0.5)
ax[1,1].hist(Bwithin_ec,normed=True,bins=20,color='r',alpha=0.5)
ax[1,1].set_xticks(np.arange(0.1,0.41,0.1))
ax[1,1].set_xlim(0.1,0.4)
x0,x1 = ax[1,1].get_xlim()
y0,y1 = ax[1,1].get_ylim()
ax[1,1].set_aspect((x1-x0)/(y1-y0))
ax[1,1].set_xlabel('B')
ax[1,1].set_ylabel('Number of Edges (%)')
ax[1,1].set_title('E',loc='left')

# <markdowncell>
# make layout tight and save
# <codecell>
fig.tight_layout()
fig.savefig('./examples/figures/bursts_fmri.pdf')
#This figure then got manually editted a bit for aestetics

# <markdowncell>
# Do (quite unnecessary) correlation on panel F data.
# <codecell>
spearman=sps.spearmanr(Beo[triind[0],triind[1]],Bec[triind[0],triind[1]])
