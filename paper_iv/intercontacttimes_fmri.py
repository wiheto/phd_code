import scipy.io as sio
import numpy as np
import teneto

#Script calculates the intercontact times for the number subjects

ict_eo = np.empty(46,dtype=list)
ict_ec = np.empty(46,dtype=list)

for s in range(0,46):
    print('Calculating for subject: ' + str(s))
    dat=sio.loadmat('./examples/data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c1.mat')['binGraph']
    dat[dat>0]=1
    ict=teneto.intercontacttimes(dat)
    ict_eo[s]=ict

    dat=sio.loadmat('./examples/data/bingraph_weightcorr_2stdth_s' + str(s+1) + '_c2.mat')['binGraph']
    dat[dat>0]=1
    ict=teneto.intercontacttimes(dat)
    ict_ec[s]=ict

np.save('./examples/data/ict_eo.npy',ict_eo)
np.save('./examples/data/ict_ec.npy',ict_ec)

ict_eo=np.load('./examples/data/ict_eo.npy')
ict_ec=np.load('./examples/data/ict_ec.npy')
