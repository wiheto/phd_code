# <markdowncell>
# Script calculates the shortest temporal path for all number of subjects. Each shortest path dictionary gets put in a numpy array
# (perhaps not the cleanest way to do this)

# <markdowncell> Import scipy.io for matlab loading, numpy and teneto
# <codecell>
import scipy.io as sio
import numpy as np
import teneto


# <markdowncell>
# Define all parameters
# <codecell>
# number of subjects
N = 46

# <markdowncell>
# Predefine numpy arrays for each subject
# <codecell>
closeness_eo = np.empty(N,dtype=list)
closeness_ec = np.empty(N,dtype=list)

# <markdowncell>
# Loop through subjects (for both conditions). Add dictionary of shortest paths each to numpy array.
# <codecell>
for s in range(0,46):
    print('--- calculating subject: ' + str(s) + ' ---')
    dat=sio.loadmat('./examples/data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c1.mat')['binGraph']
    dat[dat>0]=1
    closeness=teneto.temporalCloseness(dat)
    closeness_eo[s]=closeness

    dat=sio.loadmat('./examples/data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c2.mat')['binGraph']
    dat[dat>0]=1
    closeness=teneto.temporalCloseness(dat)
    closeness_ec[s]=closeness

# <markdowncell>
# Save ouput
# <codecell>
np.save('./examples/data/closeness_eo.npy',closeness_eo)
np.save('./examples/data/closeness_ec.npy',closeness_ec)
