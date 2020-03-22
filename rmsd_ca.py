#!/usr/bin/env python
# coding: utf-8

# In[58]:


import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
#get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import pandas as pd


# In[59]:


data_dir = '/Volumes/JGI/AAHD/round_1/'
traj_paths = [data_dir+'{}ns/100ns-production-stripped.xtc'.format(i+1) for i in range(100)]
top_path = 'data/2agy_final_min-stripped_1frame.pdb'
xtal_path = 'data/2agy_c36_state0.pdb'


# Load the crystal structure to align to, I need to get rid of the crystal waters though. 

# In[60]:


xtal = md.load(xtal_path)
xtal = xtal.atom_slice(xtal.top.select('not water'))


# Load the and process data

# In[61]:


alpha_ix = xtal.top.select('backbone and name CA')


# In[ ]:


rmsds = []
for path in traj_paths: 
    print(path)
    traj = md.load(path, top=top_path, stride=10)
    traj.superpose(reference=xtal, atom_indices = alpha_ix, ref_atom_indices = alpha_ix)
    rmsds.append(md.rmsd(target = traj, reference = xtal, atom_indices = alpha_ix))
    


# Convert to data frame and save. 

# In[ ]:


rmsd = np.concatenate(rmsds)
idx = np.concatenate([np.repeat(i+1, rmsds[i].shape[0]) for i in range(len(rmsds))])
pd.DataFrame({'rmsd': rmsd, 'traj_idx': idx}).to_csv('outputs/rmsd_ca.csv', index=False)

