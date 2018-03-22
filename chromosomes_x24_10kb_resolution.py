
# coding: utf-8

# # A Model for the Human Chromosome
# This example shows a simple simulation of a single block copolymer that is a potential model of the Human Chromosome in hoomdd.
# 
# The human genome consists of 24 chromosomes with a total size of 3Gb.  Aiden et al. 2014 produced a 1kB resolution (below single gene resolution) Genome Contact Map.  Contact domains (∼185 kb) segregate into six nuclear subcompartments with distinct histone marks.
# This can be reproduced phenomenologically with this simple block copolymer model.
# 
# from: https://lost-contact.mit.edu/afs//umich.edu/user/j/o/joaander/Public/hoomd-web/doc/page_example_scripts.html
# 
# ## Initialize
# First import hoomd and associated libraries and then initialize.  If no mode is given then the fastest processor cpu or gpu will be selected automatically.

# In[1]:


#import hoomd
import hoomd.md
from hoomd import *
from hoomd.deprecated import *
context.initialize("");
#context.initialize("--mode=cpu");
#context.initialize("--mode=gpu");


# ## Import other libraries and define parameters
# 
# The human genome has 24 chromasomes but for now we will simulate a single chromosome.

# In[2]:


import math
# parameters
phi_P = 0.25
n_poly=1


# ## Define the polymer
# The average length of a human chromosome is 133 Mb

# In[3]:


import numpy as np
genome_size=3079843747
NChr=24 #Mb
Nbeads=genome_size/10000
#chromosome=np.linspace(1,Nchromosomes,Nchromosomes,dtype=intp).tolist()
chromosome_size=np.fix(np.r_[8,7.9,6.5,6.2,5.9,5.5,5.2,4.7,4.6,4.4,4.4,4.3,3.7,3.5,3.3,2.9,2.6,2.5,2.1,2.0,1.5,1.6,5.0,1.9]/100.0*Nbeads)
chromosome=list('cdefghijklmnopqurstuvwxy')
separation_radius='=0.35, '.join(chromosome[:NChr])+'=0.35'
print separation_radius
separation_radius=dict(c=0.35, d=0.35, e=0.35, f=0.35, g=0.35, h=0.35, i=0.35, j=0.35, k=0.35, l=0.35, m=0.35, n=0.35, o=0.35, p=0.35, q=0.35, u=0.35, r=0.35, s=0.35, t=0.35)
separation_radius.update(dict(u=0.35, v=0.35, w=0.35, x=0.35, y=0.35))


# So if we assume that Euchromatin makes up 60% of the genome and heterochromatin makes up 40% then a block copolymer of 10000 beads would have a resolution of approximately 13 kb/bead the average size of a human gene.  This would correspond to approximately 18 beads per enhancer promoter loop.

# In[ ]:


#polymer1 = dict(bond_len=1.2, type=['A']*3000 + ['B']*4000 + ['A']*3000,bond="linear", count=n_poly)
polymer = [[] for i in range(NChr)]
for i in range(NChr): polymer[i] = dict(bond_len=1.2, type=chromosome[i]*int(chromosome_size[i]),bond="linear", count=1)
# perform some simple math to find the length of the box
#N = len(polymer1['type']) * polymer1['count']#+len(polymer2['type']) * polymer2['count']+len(polymer3['type']) * polymer3['count']
N=0
for i in range (NChr): N += len(polymer[i]['type']) * polymer[i]['count']#+len(polymer2['type']) * polymer2['count']+len(polymer3['type']) * polymer3['count']

# generate the polymer system
#init.create_random_polymers(box=data.boxdim(volume=10*math.pi * N / (6.0 * phi_P)), polymers=[polymer1,polymer2,polymer3],separation=separation_radius,seed=12)
init.create_random_polymers(box=data.boxdim(volume=5*math.pi * N / (6.0 * phi_P)), polymers=polymer,separation=separation_radius,seed=12)


# ## Setup the bonds and force fields

# In[ ]:


# force field setup
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('polymer', k=330.0, r0=0.84)
nl = hoomd.md.nlist.cell();
lj = hoomd.md.pair.lj(nlist=nl, r_cut=3.0)
for i in range(NChr): 
    lj.pair_coeff.set(chromosome[i],chromosome[i], epsilon=1.0, sigma=1.0, alpha=1.0)
#lj.pair_coeff.set('b','b', epsilon=1.0, sigma=1.0, alpha=1.0)
#lj.pair_coeff.set('c','c', epsilon=1.0, sigma=1.0, alpha=1.0)
#lj.pair_coeff.set('d','d', epsilon=1.0, sigma=1.0, alpha=1.0)
for i in range(NChr):
    for j in range(i+1,NChr): 
        lj.pair_coeff.set(chromosome[i],chromosome[j], epsilon=1.0, sigma=1.0, alpha=0.0, r_cut=2**(1.0/6.0))
#lj.pair_coeff.set('b','c', epsilon=1.0, sigma=1.0, alpha=0.0, r_cut=2**(1.0/6.0))
#lj.pair_coeff.set('b','d', epsilon=1.0, sigma=1.0, alpha=0.0, r_cut=2**(1.0/6.0))
#lj.pair_coeff.set('c','d', epsilon=1.0, sigma=1.0, alpha=0.0, r_cut=2**(1.0/6.0))
lj.set_params(mode='shift')


# ## Special Pairs will create our Contact Domains (not working yet)

# In[ ]:


#CTCF = hoomd.md.special_pair()
#lj_ctcf = hoomd.md.special_pair.lj(name="CTCF")
#lj_ctcf.pair_coeff.set('A','A', epsilon=1.0, sigma=1.0, r_cut=3)
#lj_ctcf.set_params(mode='shift')


# ## Integrate the simulation

# In[ ]:


all = group.all()
# integrate NVT for a bunch of time steps
hoomd.md.integrate.mode_standard(dt=0.005)
hoomd.md.integrate.nvt(group=all, kT=1.2, tau=0.5)
#hoomd.md.integrate.brownian(group=all, kT=1.2, seed=1)


# The results of this example may be visualized with VMD which can be installed from here: http://www.ks.uiuc.edu/Research/vmd/

# In[ ]:


hoomd.dump.gsd(filename="/home/mferguson/scratch/24xchromosomes_10000.gsd", overwrite=True, period=1000, group=group.all())
#hoomd.dump.gsd("gsd/24xchromosomes_10000_0.gsd", period=None, group=all, overwrite=True, time_step=0);
# setup the IMD server
#hoomd.analyze.imd(port=54321, period=1)
# run a very long time so the simulation can be watched in VMD
run(1e5)


# The results of this example may be visualized with VMD which can be installed from here: http://www.ks.uiuc.edu/Research/vmd/
# 
# While running must execute the following command in a terminal:

# In[ ]:


#!cat points.vmd
#vmd -e points.vmd gsd/24xchromosomes_10000.gsd


# ## GPU cluster
# We happen to have a high performance computing cluster on campus with 5 NVIDIA P100 GPU's.  In order to take advantage of them we need simply need to have an account and to have a little extra knowledge about how to use an HPC cluster.
# 
# R2 is a heterogeneous compute cluster provided by the BSU Office of Research. It consists of 26 compute nodes and 5 GPU nodes, each with dual Intel Xeon E5-2680 CPUs. The GPU nodes each have dual Nvidia P100 GPUs.  More information about the cluster, how to gain access to it and how to use it can be found here: https://rcs.boisestate.edu/r2/
# 
# While it is probably possible to use an ipython notebook on a server (see fast.ai for an example).  We will reduce our python code to a script shown below.

# ![](snapshots/24xchromosomes_10k.png)

# In[4]:


get_ipython().system(u'cat chromosomes_x24_10kb_resolution.py')


# ## SLURM and how to run jobs
# In addition to the python script we need also to have shell script that identifies how many and which nodes to use.  More information about the SLURM scheduler and how to use it can be found here: https://slurm.schedmd.com/quickstart.html
# and here: 
# https://www.rc.fas.harvard.edu/resources/documentation/convenient-slurm-commands/
# 
# Some useful commands are:
# sinfo
# squeue
# scancel
# 
# to start our job type:
# >sbatch hum_gen.bash
# 
# progress can be followed by typing 
# >cat logs/gpu-job-*.o

# In[3]:





# Then copy the output from the scratch folder on R2 using secure copy (scp):
# >scp mferguson@r2.boisestate.edu:scratch/24xchromosomes_10000.gsd ./gsd/

# ![](snapshots/24xchromosomes_10000.png)
