import pandas as pd
import numpy as np
import scipy.io as sio
import nibabel

pathToAAL = '/usr/share/mricron/templates/aal.nii.gz'
AAL=nibabel.load(pathToAAL)
coord = sio.matlab.loadmat('./data/coord_power264.mat')['coord']
coordras = (coord-AAL.get_affine()[:-1,-1])

dat=pd.DataFrame(index=np.arange(0,117))
R=['U','Precentral_L',
'Precentral_R',
'Frontal_Sup_L',
'Frontal_Sup_R',
'Frontal_Sup_Orb_L',
'Frontal_Sup_Orb_R',
'Frontal_Mid_L',
'Frontal_Mid_R',
'Frontal_Mid_Orb_L',
'Frontal_Mid_Orb_R',
'Frontal_Inf_Oper_L',
'Frontal_Inf_Tri_L',
'Frontal_Inf_Oper_R',
'Frontal_Inf_Tri_R',
'Frontal_Inf_Orb_L',
'Frontal_Inf_Orb_R',
'Rolandic_Oper_L',
'Supp_Motor_Area_L',
'Rolandic_Oper_R',
'Supp_Motor_Area_R',
'Olfactory_L',
'Olfactory_R',
'Frontal_Sup_Medial_R',
'Frontal_Med_Orb_L',
'Frontal_Sup_Medial_L',
'Frontal_Med_Orb_R',
'Rectus_L',
'Rectus_R',
'Insula_L',
'Cingulum_Ant_L',
'Cingulum_Ant_R',
'Cingulum_Mid_L',
'Insula_R',
'Cingulum_Mid_R',
'Cingulum_Post_L',
'Cingulum_Post_R',
'Hippocampus_L',
'Hippocampus_R',
'ParaHippocampal_L',
'ParaHippocampal_R',
'Amygdala_L',
'Amygdala_R',
'Calcarine_L',
'Calcarine_R',
'Cuneus_L',
'Cuneus_R',
'Lingual_R',
'Lingual_L',
'Occipital_Sup_L',
'Occipital_Sup_R',
'Occipital_Mid_L',
'Occipital_Mid_R',
'Occipital_Inf_L',
'Occipital_Inf_R',
'Fusiform_L',
'Fusiform_R',
'Postcentral_L',
'Postcentral_R',
'Parietal_Sup_L',
'Parietal_Sup_R',
'Parietal_Inf_L',
'Parietal_Inf_R',
'SupraMarginal_L',
'SupraMarginal_R',
'Angular_L',
'Angular_R',
'Precuneus_L',
'Precuneus_R',
'Paracentral_Lobule_L',
'Paracentral_Lobule_R',
'Caudate_R',
'Putamen_L',
'Caudate_L',
'Putamen_R',
'Pallidum_L',
'Pallidum_R',
'Thalamus_L',
'Thalamus_R',
'Heschl_L',
'Heschl_R',
'Temporal_Sup_L',
'Temporal_Sup_R',
'Temporal_Pole_Sup_L',
'Temporal_Pole_Sup_R',
'Temporal_Mid_L',
'Temporal_Mid_R',
'Temporal_Pole_Mid_L',
'Temporal_Pole_Mid_R',
'Temporal_Inf_L',
'Temporal_Inf_R',
'Cerebelum_Crus1_L',
'Cerebelum_Crus1_R',
'Cerebelum_Crus2_L',
'Cerebelum_Crus2_R',
'Cerebelum_3_L',
'Cerebelum_4_5_L',
'Cerebelum_3_R',
'Cerebelum_4_5_R',
'Cerebelum_6_L',
'Cerebelum_6_R',
'Cerebelum_7b_L',
'Cerebelum_7b_R',
'Cerebelum_8_L',
'Cerebelum_8_R',
'Cerebelum_9_L',
'Cerebelum_9_R',
'Cerebelum_10_L',
'Cerebelum_10_R',
'Vermis_1_2',
'Vermis_3',
'Vermis_4_5',
'Vermis_6',
'Vermis_7',
'Vermis_9',
'Vermis_10',
'Vermis_8']
dat['AAL']=R
dat.to_csv('./data/aalinfo.csv')

netid=np.array(list(map(int,sio.loadmat('/home/william/work/teneto/examples_article/data/networkassignment')['PowerNetClass'])))
netid[netid==-1]=13
network = np.array([1,3,4,5,7,8,9,10,11,12,13])
netlab = np.array(['SM','CO','AU','DM','V','FP','SA','Sub','VA','DA','U'])
netvec=[]
for n in range(0,264):
    netvec.append(str(netlab[np.where(network==netid[n])[0][0]]))

PowerInfo=pd.DataFrame(index=range(0,264),columns={'coord_x','coord_y','coord_z','network','aal','note'})
PowerInfo['coord_x']=coord[:,0]
PowerInfo['coord_y']=coord[:,1]
PowerInfo['coord_z']=coord[:,2]
PowerInfo['network']=netvec

aal=[]
for n in range(0,264):
    aal.append(R[AAL.get_data()[int(coordras[n,0]),int(coordras[n,1]),int(coordras[n,2])]])

PowerInfo['aal']=aal

#Then comes manual work with PowerInfo[['coord_x','coord_y','coord_z']].where(PowerInfo['aal']=='U').dropna()
#

for n in range(0,264):
    PowerInfo['note'][n]=''

PowerInfo['aal'][13]=R[33]
PowerInfo['note'][13]='closet-2'

PowerInfo['aal'][27]=R[2]
PowerInfo['note'][27]='closet-2 (equidistant to closet Postcentral_R. Prioritizing Z-axis)'

PowerInfo['aal'][51]=R[29]
PowerInfo['note'][51]='closet-2'

PowerInfo['aal'][56]=R[29]
PowerInfo['note'][56]='closet-3 (equidistant to Caudate_L. Priortizing Z-axis)'

PowerInfo['aal'][59]=R[30]
PowerInfo['note'][59]='closet-1'

PowerInfo['aal'][60]=R[80]
PowerInfo['note'][60]='closet-1'

PowerInfo['aal'][72]=R[79]
PowerInfo['note'][72]='closet-2'

PowerInfo['aal'][90]=R[67]
PowerInfo['note'][90]='closet-1'

PowerInfo['aal'][211]=R[31]
PowerInfo['note'][211]='closet-3'

PowerInfo['aal'][220]=R[34]
PowerInfo['note'][220]='closet-1'

PowerInfo['note'][225]='undefined'

PowerInfo['aal'][227]=R[75]
PowerInfo['note'][227]='closet-2'


PowerInfo.to_csv('./data/PowerInfo.csv')
