
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
chromosome=list('bcdefghijklmnopqurstuvwxy')
s="["+"]*%d+[".join(map(str, chromosome))+"]*%d" 
t=chromosome_size.astype(int).tolist()
#print t
#type_string='['+''.join(str(e) for e in [val for pair in zip(chromosome,N*']*', chromosome_size.astype(int),N*'+[]') for val in pair])
type_string=['']
for i in range(NChr): type_string += "['%s']*%d" % (chromosome[i],chromosome_size[i]), 
type_string=' + '.join(type_string)
type_string=type_string[3:]
print type_string
#separation_radius=dict(b=0.35, c=0.35)
separation_radius='=0.35, '.join(chromosome) 
print ('%s=0.35' % (separation_radius))
#print dict(b=0.35, c=0.35, d=0.35, e=0.35, f=0.35, g=0.35, h=0.35, i=0.35, j=0.35, k=0.35, l=0.35, m=0.35, n=0.35, o=0.35, p=0.35, q=0.35, u=0.35, r=0.35, s=0.35, t=0.35)
#print dict(u=0.35, v=0.35, w=0.35, x=0.35, y=0.35)


# So if we assume that Euchromatin makes up 60% of the genome and heterochromatin makes up 40% then a block copolymer of 10000 beads would have a resolution of approximately 13 kb/bead the average size of a human gene.  This would correspond to approximately 18 beads per enhancer promoter loop.

# In[ ]:


type_string=['b']*24638 + ['c']*24330 + ['d']*20018 + ['e']*19095 + ['f']*18171 + ['g']*16939 + ['h']*16015 + ['i']*14475 + ['j']*14167 + ['k']*13551 + ['l']*13551 + ['m']*13243 + ['n']*11395 + ['o']*10779 + ['p']*10163 + ['q']*8931 + ['u']*8007 + ['r']*7699 + ['s']*6467 + ['t']*6159 + ['u']*4619 + ['v']*4927 + ['w']*15399 + ['x']*5851
#polymer1 = dict(bond_len=1.2, type=['A']*3000 + ['B']*4000 + ['A']*3000,bond="linear", count=n_poly)
polymer1 = dict(bond_len=1.2, type=type_string,bond="linear", count=n_poly)
#polymer1 = dict(bond_len=1.2, type=['B']*300,bond="linear", count=n_poly)
#polymer2 = dict(bond_len=1.2, type=['C']*400,bond="linear", count=n_poly)
#polymer3 = dict(bond_len=1.2, type=['D']*300,bond="linear", count=n_poly)
# perform some simple math to find the length of the box
N = len(polymer1['type']) * polymer1['count']#+len(polymer2['type']) * polymer2['count']+len(polymer3['type']) * polymer3['count']
# generate the polymer system
separation_radius=dict(b=0.35, c=0.35, d=0.35, e=0.35, f=0.35, g=0.35, h=0.35, i=0.35, j=0.35, k=0.35, l=0.35, m=0.35, n=0.35, o=0.35, p=0.35, q=0.35, u=0.35, r=0.35, s=0.35, t=0.35)
separation_radius.update(dict(u=0.35, v=0.35, w=0.35, x=0.35, y=0.35))
#separation_radius=dict(B=0.35, C=0.35, D=0.35)
#init.create_random_polymers(box=data.boxdim(volume=10*math.pi * N / (6.0 * phi_P)), polymers=[polymer1,polymer2,polymer3],separation=separation_radius,seed=12)
init.create_random_polymers(box=data.boxdim(volume=10*math.pi * N / (6.0 * phi_P)), polymers=[polymer1],separation=separation_radius,seed=12)


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


#!cat imd.vmd
#vmd -e imd.vmd gsd/24xchromosomes_10000_0.gsd


# ![](snapshots/24xchromosomes_10k.png)
