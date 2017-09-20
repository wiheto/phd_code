
# <markdowncell> Import scipy.io for matlab loading, numpy and teneto
# <codecell>
import scipy.io as sio
import scipy.stats as sps
import numpy as np
import teneto
import matplotlib.pyplot as plt
import plot_surf
import pandas as pd
import os
import markdown2
import tenetostats

coord = sio.matlab.loadmat('./data/coord_power264.mat')['coord']
netid=np.array(list(map(int,sio.loadmat('./data/networkassignment')['PowerNetClass'])))
netid[netid==-1]=13
network = np.array([1,3,4,5,7,8,9,10,11,12,13])
netlab = np.array(['SM','CO','AU','DM','V','FP','SA','Sub','VA','DA','U'])
plotnet = [5,1,7,8,9,3,4,10,11,12,13]
plotorder=[]
for n in plotnet:
    fid = np.where(netid==n)[0]
    plotorder=plotorder + list(fid)
plotorder = np.array(plotorder)
PowerInfo = pd.read_csv('./data/PowerInfo.csv')

plotcol = np.zeros([14,3])
plotcol[1,:] = [255,102,0]
plotcol[3,:] = [255,255,0]
plotcol[4,:] = [128,0,128]
plotcol[5,:] = [128,0,0]
plotcol[7,:] = [255,0,100]
plotcol[8,:] = [124,124,220]
plotcol[9,:] = [0,255,0]
plotcol[10,:] = [255,0,255]
plotcol[11,:] = [0,255,255]
plotcol[12,:] = [0,128,0]
plotcol[13,:] = [0,0,0]
plotcol=plotcol/255

colvec=np.zeros([264,3])
for n in range(0,264):
    colvec[n,:]=plotcol[netid[n],:]




# <markdowncell> Set matplotlib color style
# <codecell>
plt.rcParams['image.cmap'] = 'gist_gray'

# <markdowncell> load data
# <codecell>
closeness_eo=np.load('./data/closeness_eo.npy')
closeness_ec=np.load('./data/closeness_ec.npy')
degree_eo=np.load('./data/degree_eo.npy')
degree_ec=np.load('./data/degree_ec.npy')
# Bec=np.load('./data/burst_ec.npy')
# Beo=np.load('./data/burst_eo.npy')
# CEeo=np.load('./data/efficiency_centrality_eo.npy')thon
# CEec=np.load('./data/efficiency_centrality_ec.npy')

closeness_eo = np.stack(closeness_eo)
closeness_ec = np.stack(closeness_ec)
degree_eo = np.stack(degree_eo)
degree_ec = np.stack(degree_ec)

# <markdowncell> Plot on brains
# <codecell>
plotvar=np.mean(degree_eo,axis=0)
plotvar=(plotvar-plotvar.min())/(plotvar.max()-plotvar.min())*5
plot_surf.plot_brain_surface(coord,plotvar,colvec,'./figures/tdegree_brain.png')

plotvar=np.mean(closeness_eo,axis=0)
plotvar=(plotvar-plotvar.min())/(plotvar.max()-plotvar.min())*5
plot_surf.plot_brain_surface(coord,plotvar,colvec,'./figures/closeness_brain.png')

plotvar=np.mean(degree_ec,axis=0)
plotvar=(plotvar-plotvar.min())/(plotvar.max()-plotvar.min())*5
plot_surf.plot_brain_surface(coord,plotvar,colvec,'./figures/tdegree_ec_brain.png')

plotvar=np.mean(closeness_ec,axis=0)
plotvar=(plotvar-plotvar.min())/(plotvar.max()-plotvar.min())*5
plot_surf.plot_brain_surface(coord,plotvar,colvec,'./figures/closeness_ec_brain.png')


# plotvar=np.mean(Beo,axis=0)
# plotvar=(plotvar-plotvar.min())/(plotvar.max()-plotvar.min())*5
# plot_surf.plot_brain_surface(coord,plotvar,colvec,'./figures/bursts_brain.png')
#

# <markdowncell> Plot on brains
# <codecell>
pth=np.zeros([1000,264])
np.random.seed(2017)
for p in range(0,1000):
    porder= np.argsort(np.random.rand(264,46),axis=0)
    pdeg = np.zeros(264)
    for s in range(0,46):
        pdeg += degree_eo[s,porder[:,s]]
    pth[p,:]=pdeg/46
pth=np.sort(pth,axis=0)
thD=pth[950,:].max()

pth=np.zeros([1000,264])
np.random.seed(2017)
for p in range(0,1000):
    porder= np.argsort(np.random.rand(264,46),axis=0)
    pdeg = np.zeros(264)
    for s in range(0,46):
        pdeg += degree_ec[s,porder[:,s]]
    pth[p,:]=pdeg/46
pth=np.sort(pth,axis=0)
thDec=pth[950,:].max()


pth=np.zeros([1000,264])
np.random.seed(2017)
for p in range(0,1000):
    porder= np.argsort(np.random.rand(264,46),axis=0)
    pdeg = np.zeros(264)
    for s in range(0,46):
        pdeg += closeness_ec[s,porder[:,s]]
    pth[p,:]=pdeg/46
pth=np.sort(pth,axis=0)
thCec=pth[950,:].max()


pth=np.zeros([1000,264])
np.random.seed(2017)
for p in range(0,1000):
    porder= np.argsort(np.random.rand(264,46),axis=0)
    pdeg = np.zeros(264)
    for s in range(0,46):
        pdeg += closeness_eo[s,porder[:,s]]
    pth[p,:]=pdeg/46
pth=np.sort(pth,axis=0)
thC=pth[950,:].max()




eD = np.mean(degree_eo,axis=0)
eC = np.mean(closeness_eo,axis=0)
eDec = np.mean(degree_ec,axis=0)
eCec = np.mean(closeness_ec,axis=0)




def print_table_of_node_and_all(val,valname,threshold,sname,title,tag='1'):
    PowerInfo = pd.read_csv('./data/PowerInfo.csv')
    sigVals=val[np.where(val>threshold)[0]]
    results=PowerInfo.iloc[PowerInfo.index[val>threshold]]
    results=results[['coord_x','coord_y','coord_z','network','aal']]
    results.rename(columns={'coord_x':'x','coord_y':'y','coord_z':'z','aal':'AAL','network':'Network'},inplace=True)
    results[valname]=np.around(sigVals,3)
    results.sort_values(by=valname,ascending=False,inplace=True)
    # Get column names
    cols = results.columns
    # Create a new DataFrame with just the markdown
    # strings
    hdrresults = pd.DataFrame([['---',]*len(cols)], columns=cols)
    #Create a new concatenated DataFrame
    results = pd.concat([hdrresults, results])
    results.to_csv(sname + '.md', sep="|", index=False)
    with open(sname + ".md", "a") as myfile:
        myfile.write("\n \pagenumbering{gobble} \n Table " + tag + ": " + title)
    os.system('pandoc ' + sname + '.md -o ' + sname + '.pdf')


print_table_of_node_and_all(eD,'D',thD,'degree_eo_top','Temporal Degree Centrality during eyes open condition. Nodes with degree centrality where p\<0.05, their designated network and corresponding AAL. XYZ are MNI cordinates.','1')
print_table_of_node_and_all(eC,'C',thC,'closeness_eo_top','Closeness Centrality during eyes open condition. Nodes with closeness centrality where p\<0.05, their designated network and corresponding AAL. XYZ are MNI cordinates.','3')

print_table_of_node_and_all(eDec,'D',thDec,'degree_ec_top','Temporal Degree Centrality during eyes closed condition. Nodes with degree centrality where p\<0.05, their designated network and corresponding AAL. XYZ are MNI cordinates.','2')
print_table_of_node_and_all(eCec,'C',thCec,'closeness_ec_top','Closeness Centrality during eyes closed condition. Nodes with closeness centrality where p\<0.05, their designated network and corresponding AAL. XYZ are MNI cordinates.','4')


B=np.mean(Beo,axis=0)
thB = np.sort(B)[-27]
print_table_of_node_and_all(B,'B',thB,'burstiness_eo_top','Burstiness during eyes open condition. Top 10 percent of nodes, their designated network and corresponding AAL. XYZ are MNI cordinates.')
B=np.mean(Bec,axis=0)
thB = np.sort(B)[-27]
print_table_of_node_and_all(B,'B',thB,'burstiness_ec_top','Burstiness during eyes closed condition. Top 10 percent of nodes, their designated network and corresponding AAL. XYZ are MNI cordinates.')



def spornlike_barplot(dat,col,ylim,th,ax):
    sig = np.where(dat>th)[0]
    ax=fig.add_subplot(212)
    ax.bar(range(0,264),dat,edgecolor=[0.45,0.45,0.45],facecolor=[0.45,0.45,0.45],width=1.0)
    ax.bar(range(0,len(sig)),dat[sig],edgecolor=[1,1,0],facecolor=[1,1,0],width=1.0)
    ax.set_xlim([-0.5,264.5])
    ax.set_ylim(ylim)
    ax.plot(range(-1,265),np.zeros(266)+th,linestyle='--',color='k')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax=fig.add_subplot(211)
    ax.scatter(range(0,264),np.zeros(264),s=30,edgecolor='none',c=col)
    ax.set_ylim([-1,1])
    ax.set_xlim([-0.5,264.5])
    #ax[0].set_aspect('equal')
    ax.axis('off')
    return ax


sorted_dat= np.array(list(reversed(np.sort(eC))))
sorted_net = netid[np.array(list(reversed(np.argsort(eC))))]
fig = plt.figure(figsize=(20,5))
spornlike_barplot(sorted_dat,plotcol[sorted_net,:],[.155,.18],thC,ax)
plt.tight_layout()
fig.show()
fig.savefig('./figures/closeness_eo_spornsplot.pdf',r=600)

sorted_dat= np.array(list(reversed(np.sort(eCec))))
sorted_net = netid[np.array(list(reversed(np.argsort(eCec))))]
fig = plt.figure(figsize=(20,5))
spornlike_barplot(sorted_dat,plotcol[sorted_net,:],[.145,.17],thCec,ax)
plt.tight_layout()
fig.show()
fig.savefig('./figures/closeness_ec_spornsplot.pdf',r=600)

sorted_dat= np.array(list(reversed(np.sort(eD))))
sorted_net = netid[np.array(list(reversed(np.argsort(eD))))]
fig = plt.figure(figsize=(20,5))
spornlike_barplot(sorted_dat,plotcol[sorted_net,:],[1800,1950],thD,ax)
plt.tight_layout()
fig.show()
fig.savefig('./figures/tdegree_eo_spornsplot.pdf',r=600)

sorted_dat= np.array(list(reversed(np.sort(eDec))))
sorted_net = netid[np.array(list(reversed(np.argsort(eDec))))]
fig = plt.figure(figsize=(20,5))
spornlike_barplot(sorted_dat,plotcol[sorted_net,:],[1800,1950],thDec,ax)
plt.tight_layout()
fig.show()
fig.savefig('./figures/tdegree_ec_spornsplot.pdf',r=600)



eC_net=np.zeros(11)
eCec_net=np.zeros(11)
for i,net in enumerate(plotnet):
    eC_net[i] = np.mean(eC[netid==net])
    eCec_net[i] = np.mean(eCec[netid==net])


pdegnode=np.zeros(264)
for r in range(0,264):
    ptmp,tmp=tenetostats.shufflegroups(degree_eo[:,r],degree_ec[:,r],pnum=10000)
    pdegnode[r]=ptmp

pdegnode_correct=teneto.misc.correct_pvalues_for_multiple_testing(pdegnode)
sig = np.where(pdegnode_correct<0.05)[0]
degree_eo[:,sig]-degree_ec[:,sig]
np.mean(degree_eo[:,sig],axis=0)-np.mean(degree_ec[:,sig],axis=0)



pclosenode=np.zeros(264)
for r in range(0,264):
    ptmp,tmp=tenetostats.shufflegroups(closeness_eo[:,r],closeness_ec[:,r],pnum=10000)
    pclosenode[r]=ptmp
pclosenode_correct=teneto.misc.correct_pvalues_for_multiple_testing(pclosenode)
sig = np.where(pclosenode_correct<0.05)[0]
np.mean(closeness_eo[:,sig],axis=0)-np.mean(closeness_ec[:,sig],axis=0)
