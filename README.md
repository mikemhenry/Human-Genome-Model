
# A Model for the Human Genome
This example shows a simple simulation of 24 block copolymers that is a potential model of the Human Genome in hoomd.

The human genome consists of 24 chromosomes with a total size of 3Gb.  Aiden et al. 2014 produced a 1kB resolution (below single gene resolution) Genome Contact Map.  Contact domains (∼185 kb) segregate into six nuclear subcompartments with distinct histone marks.
This can be reproduced phenomenologically with this simple block copolymer model.

from: https://lost-contact.mit.edu/afs//umich.edu/user/j/o/joaander/Public/hoomd-web/doc/page_example_scripts.html

## Initialize
First import hoomd and associated libraries and then initialize.  If no mode is given then the fastest processor cpu or gpu will be selected automatically.


```python

```


```python
#import hoomd
import hoomd.md
from hoomd import *
from hoomd.deprecated import *
#context.initialize("--mode=cpu");
context.initialize("--mode=gpu");
```

    HOOMD-blue 2.2.4-unknown CUDA (8.0) DOUBLE HPMC_MIXED SSE SSE2 SSE3 SSE4_1 SSE4_2 AVX AVX2 
    Compiled: 03/11/2018
    Copyright 2009-2017 The Regents of the University of Michigan.
    -----
    You are using HOOMD-blue. Please cite the following:
    * J A Anderson, C D Lorenz, and A Travesset. "General purpose molecular dynamics
      simulations fully implemented on graphics processing units", Journal of
      Computational Physics 227 (2008) 5342--5359
    * J Glaser, T D Nguyen, J A Anderson, P Liu, F Spiga, J A Millan, D C Morse, and
      S C Glotzer. "Strong scaling of general-purpose molecular dynamics simulations
      on GPUs", Computer Physics Communications 192 (2015) 97--107
    -----
    notice(2): This system is not compute exclusive, using local rank to select GPUs
    HOOMD-blue is running on the following GPU(s):
     [0]   GeForce GTX 1050 Ti   6 SM_6.1 @ 1.39 GHz, 4095 MiB DRAM, DIS


## Import other libraries and define parameters

The human genome has 24 chromosomes.


```python
import math
# parameters
phi_P = 0.25
n_poly = 24
```

## Define the polymer
The average length of a human chromosome is 133 Mb


```python
3200/24 #Mb
```




    133



So if we assume that Euchromatin makes up 60% of the genome and heterochromatin makes up 40% then a block copolymer of 10000 beads would have a resolution of approximately 13 kb/bead the average size of a human gene.  This would correspond to approximately 18 beads per enhancer promoter loop.


```python
#polymer1 = dict(bond_len=1.2, type=['A']*300 + ['B']*400 + ['A']*300,bond="linear", count=n_poly)
polymer1 = dict(bond_len=1.2, type=['A']*3000 + ['B']*4000 + ['A']*3000,bond="linear", count=n_poly)
# perform some simple math to find the length of the box
N = len(polymer1['type']) * polymer1['count']
# generate the polymer system
init.create_random_polymers(box=data.boxdim(volume=10*math.pi * N / (6.0 * phi_P)), polymers=[polymer1],separation=dict(A=0.35, B=0.35),seed=12)

```

    notice(2): Group "all" created containing 240000 particles





    <hoomd.data.system_data at 0x10ba2b050>



## Setup the bonds and force fields


```python
# force field setup
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('polymer', k=330.0, r0=0.84)
nl = hoomd.md.nlist.cell();
lj = hoomd.md.pair.lj(r_cut=3.0,nlist=nl)

lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0, alpha=0.0,r_cut=2**(1.0/6.0))
lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0, alpha=0.0,r_cut=2**(1.0/6.0))
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0, alpha=1.0,r_cut=3.0)
lj.set_params(mode='shift')
                  

```

## Integrate the simulation


```python
all = group.all()
# integrate NVT for a bunch of time steps
hoomd.md.integrate.mode_standard(dt=0.005)
hoomd.md.integrate.nvt(group=all, kT=1.2, tau=0.5)
#hoomd.md.integrate.brownian(group=all, kT=1.2, seed=1)
```




    <hoomd.md.integrate.nvt at 0x11c966f10>



The results of this example may be visualized with VMD which can be installed from here: http://www.ks.uiuc.edu/Research/vmd/


```python
#hoomd.dump.gsd(filename="hum_gen_24x1000.gsd", overwrite=True, period=None, group=group.all(), time_step=0)
hoomd.dump.gsd("hum_gen_24x10000_0.gsd", period=None, group=all, overwrite=True, time_step=0);
# setup the IMD server
hoomd.analyze.imd(port=54321, period=100)
# run a very long time so the simulation can be watched in VMD
run(1e6)
```

    notice(2): -- Neighborlist exclusion statistics -- :
    notice(2): Particles with 1 exclusions             : 48
    notice(2): Particles with 2 exclusions             : 239952
    notice(2): Neighbors included by diameter          : no
    notice(2): Neighbors excluded when in the same body: no
    ** starting run **
    notice(2): analyze.imd: listening on port 54321
    Time 00:01:17 | Step 509 / 1000000 | TPS 50.8886 | ETA 05:27:20
    Time 00:01:27 | Step 872 / 1000000 | TPS 36.2858 | ETA 07:38:54
    Time 00:01:37 | Step 1401 / 1000000 | TPS 52.8447 | ETA 05:14:56
    Time 00:01:47 | Step 1971 / 1000000 | TPS 56.9455 | ETA 04:52:06
    Time 00:01:57 | Step 2508 / 1000000 | TPS 53.5007 | ETA 05:10:44
    notice(2): analyze.imd: accepted connection
    Time 00:02:08 | Step 3026 / 

    **ERROR**: 

    1000000 | TPS 51.5442 | ETA 05:22:22
    Time 00:02:18 | Step 3493 / 1000000 | TPS 46.6726 | ETA 05:55:50
    Time 00:02:28 | Step 3936 / 1000000 | TPS 44.2787 | ETA 06:14:55
    Time 00:02:38 | Step 4386 / 1000000 | TPS 44.8438 | ETA 06:10:01
    Time 00:02:48 | Step 4824 / 1000000 | TPS 43.7708 | ETA 06:18:56
    Time 00:02:58 | Step 5337 / 1000000 | TPS 51.2616 | ETA 05:23:23
    Time 00:03:08 | Step 5870 / 1000000 | TPS 53.1044 | ETA 05:12:00
    Time 00:03:18 | Step 6395 / 1000000 | TPS 52.441 | ETA 05:15:47
    Time 00:03:28 | Step 6910 / 1000000 | TPS 51.4536 | ETA 05:21:40
    Time 00:03:38 | Step 7421 / 1000000 | TPS 50.9075 | ETA 05:24:57
    Time 00:03:48 | Step 7925 / 1000000 | TPS 50.3913 | ETA 05:28:07
    Time 00:03:58 | Step 8421 / 1000000 | TPS 49.4118 | ETA 05:34:27
    Time 00:04:08 | Step 8907 / 1000000 | TPS 48.5599 | ETA 05:40:09
    Time 00:04:18 | Step 9388 / 1000000 | TPS 48.0685 | ETA 05:43:28
    Time 00:04:28 | Step 9868 / 1000000 | TPS 47.9683 | ETA 05:44:01
    Time 00:04:38 | Step 10341 / 1000000 | TPS 47.2142 | ETA 05:49:21
    Time 00:04:48 | Step 10812 / 1000000 | TPS 47.0739 | ETA 05:50:13
    Time 00:04:58 | Step 1128

The results of this example may be visualized with VMD which can be installed from here: http://www.ks.uiuc.edu/Research/vmd/

While running must execute the following command in a terminal:


```python
#vmd -e imd4.vmd
!cat imd4.vmd
```

    mol load gsd hum_gen_24x10000.gsd
    set sel [atomselect top all]
    $sel set radius 0
    #imd connect localhost 54321


## GPU cluster
We happen to have a high performance computing cluster on campus with 5 NVIDIA P100 GPU's.  In order to take advantage of them we need simply need to have an account and to have a little extra knowledge about how to use an HPC cluster.

R2 is a heterogeneous compute cluster provided by the BSU Office of Research. It consists of 26 compute nodes and 5 GPU nodes, each with dual Intel Xeon E5-2680 CPUs. The GPU nodes each have dual Nvidia P100 GPUs.  More information about the cluster, how to gain access to it and how to use it can be found here: https://rcs.boisestate.edu/r2/

While it is probably possible to use an ipython notebook on a server (see fast.ai for an example).  We will reduce our python code to a script shown below.


```python
!cat hum_gen.py
```

    import hoomd
    import hoomd.md
    from hoomd.deprecated import *
    
    hoomd.context.initialize("");
    import math
    # parameters
    phi_P = 0.25
    n_poly = 24
    T = 1.2
    polymer1 = dict(bond_len=1.2, type=['A']*3000 + ['B']*4000 + ['A']*3000,bond="linear", count=n_poly)
    # perform some simple math to find the length of the box
    N = len(polymer1['type']) * polymer1['count']
    # generate the polymer system
    init.create_random_polymers(box=hoomd.data.boxdim(volume=10*math.pi * N / (6.0 * phi_P)), polymers=[polymer1],separation=dict(A=0.35, B=0.35),seed=12)
    # force field setup
    harmonic = hoomd.md.bond.harmonic()
    harmonic.bond_coeff.set('polymer', k=330.0, r0=0.84)
    nl = hoomd.md.nlist.cell();
    lj = hoomd.md.pair.lj(r_cut=3.0,nlist=nl)
    lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0, alpha=0.0,r_cut=2**(1.0/6.0))
    lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0, alpha=0.0,r_cut=2**(1.0/6.0))
    lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0, alpha=1.0)
    all = hoomd.group.all()
    hoomd.dump.gsd("/home/mferguson/scratch/hum_gen_24x10000.gsd", period=1000, group=all, overwrite=True);
    #hoomd.analyze.imd(port=54321, period=10000)
    # integrate NVT for a bunch of time steps
    hoomd.md.integrate.mode_standard(dt=0.005)
    hoomd.md.integrate.nvt(group=all, kT=1.2, tau=0.5)
    #run(2000)
    # uncomment the next run() command if you have a few hours to spare
    # running this on a GPU the resulting dump files should show the
    # polymers self-assembling into the hex phase
    hoomd.run(1e5)


## SLURM and how to run jobs
In addition to the python script we need also to have shell script that identifies how many and which nodes to use.  More information about the SLURM scheduler and how to use it can be found here: https://slurm.schedmd.com/quickstart.html
and here: 
https://www.rc.fas.harvard.edu/resources/documentation/convenient-slurm-commands/

Some useful commands are:
sinfo
squeue
scancel

to start our job type:
>sbatch hum_gen.bash

progress can be followed by typing 
>cat logs/gpu-job-*.o


```python
!cat hum_gen.bash
```

    #!/bin/bash
    #SBATCH --mail-user=mattferguson@boisestate.edu
    #SBATCH --mail-type=ALL	#This can be set to all, end, or other options
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --time=24:00:00         # 24 hours 
    #SBATCH --output=logs/gpu-job-%j.o      # output file
    #SBATCH --error=logs/gpu-job-%j.e       # error file
    #SBATCH --partition=gpuq          # GPU partition
    #SBATCH --gres=gpu:1            # single GPU job
    
    ulimit -u 9999
    ulimit -s unlimited
    ulimit -v unlimited
    
    
    module load hoomd-blue/gcc/mvapich2/2.1.5
    
    mpirun python hum_gen.py



```python

```
