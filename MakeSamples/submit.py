import os
from JChaves.Tools.CRABSubmitter import *

ee = False
if ee:
    for i in range(0,5):
        os.system('cp ee/step1_'+str(1000+1000*i)+'_cfg.py step1_cfg.py')
        crab2_submit('nstep',100000,1000+1000*i)
        os.system('rm step1_cfg.py')
else:
    for i in range(0,5):
        os.system('cp mumu/step1_'+str(1000+1000*i)+'_cfg.py step1_cfg.py')
        crab2_submit('nstep',100000,'mumu',1000+1000*i)
        os.system('rm step1_cfg.py')
