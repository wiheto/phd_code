import scipy.io as sio
import numpy as np
import teneto


vol_eo=np.zeros(46)
vol_ec=np.zeros(46)

fluct_eo=np.zeros(46)
fluct_ec=np.zeros(46)

for s in range(0,46):
    print('Calculating for subject: ' + str(s))
    dat=sio.loadmat('./examples/data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c1.mat')['binGraph']
    dat[dat>0]=1
    fluct_eo[s]=teneto.fluctuability(dat)
    vol_eo[s]=teneto.volatility(dat)

    dat=sio.loadmat('./examples/data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c2.mat')['binGraph']
    dat[dat>0]=1
    fluct_ec[s]=teneto.fluctuability(dat)
    vol_ec[s]=teneto.volatility(dat)

np.save('./examples/data/vol_ec.npy',vol_ec)
np.save('./examples/data/vol_eo.npy',vol_eo)
np.save('./examples/data/fluct_ec.npy',fluct_ec)
np.save('./examples/data/fluct_eo.npy',fluct_eo)
