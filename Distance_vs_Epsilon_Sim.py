import hoomd
from hoomd import *
from hoomd import md
from hoomd import deprecated
import math
import numpy

#hoomd.context.initialize();

#hoomd.context.initialize("--mode=gpu");

def Simulate(eps,temp):
    
    
    hoomd.context.initialize("--mode=cpu")
    
    # parameters (polymer physics)
    phi_P = 0.25 # this has to do with how good of a solvent the polymer is in and is only used to calculate the box size below
    n_poly = 10 # n_poly is the number of polymers in the simulation
    n_beads = 30
    
    # This polymer can be made up of two different types A and B.  It's called a block copolymer.
    polymer1 = dict(bond_len=1.2, type=['B']*1 + ['A']*int(n_beads) + ['B']*1, bond="linear", count=n_poly)
    # perform some simple math to find the length of the box
    N = len(polymer1['type']) * polymer1['count']
    # generate the polymer system
    system=deprecated.init.create_random_polymers(\
                                                  box=data.boxdim(volume=1000*math.pi * N / (6.0 * phi_P)), \
                                                  polymers=[polymer1],\
                                                  separation=dict(A=0.35, B=0.35)\
                                                  ,seed=12)
    
    # force field setup
    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set('polymer', k=330.0, r0=0.84)
    nl = md.nlist.cell();
    lj = md.pair.lj(r_cut=3.0,nlist=nl)
    lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0, alpha=0.0)
    lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0, alpha=0.0)
    lj.pair_coeff.set('B', 'B', epsilon= eps, sigma=1.0, alpha=1.0)
    
    
    
   
    
    # integrate NVT for a bunch of time steps
    all = group.all()
    md.integrate.mode_standard(dt=0.0005)
    md.integrate.nvt(group=all, kT= temp, tau=0.5)
    
    #save a gsd polymer configuration every 100 timesteps
    
    hoomd.dump.gsd("Distance_vs_Epsilon_Sim.gsd", period=100, group=all, overwrite=True);
    
    #log position of particle 0
    polymers=['p1','p2','p3','p4','p5','p6','p7','p8','p9','p10']
    log = hoomd.analyze.log(filename="ree.log",
                            quantities=polymers[0:n_poly],
                            period=100,
                            overwrite=True);
    
    for i in range(n_poly):
        print(polymers[i],(n_beads+2)*i,(n_beads+2)*i+n_beads+1)
        log.register_callback(polymers[i], lambda timestep,i=i: math.sqrt( \
                                                                          (system.particles[(n_beads+2)*i].position[0]-system.particles[(n_beads+2)*i+n_beads+1].position[0])**2\
                                                                          +(system.particles[(n_beads+2)*i].position[1]-system.particles[(n_beads+2)*i+n_beads+1].position[1])**2\
                                                                          +(system.particles[(n_beads+2)*i].position[2]-system.particles[(n_beads+2)*i+n_beads+1].position[2])**2));
        
    run(1e6)


#Run simulation num_runs times incrementing epsilon. 
#Take Temperature as an input parameter  
    
num_runs = 4
def L_vs_eps(temp):
    AvgL_vs_eps = numpy.array([])
    for i in range(num_runs):
        eps = 25 + 2*i
        Simulate(eps,temp)
        ree_data = numpy.genfromtxt(fname='ree.log', skip_header=True)
        AvgL_vs_eps = numpy.append(AvgL_vs_eps,numpy.mean(ree_data[:,1:]))
    return AvgL_vs_eps #Average distance between ends versus epsilon data


AvgL_vs_eps = L_vs_eps(3)

print(AvgL_vs_eps)

from matplotlib import pyplot as plt
%matplotlib notebook

x = numpy.array([])
for i in range(num_runs):
    x = numpy.append(x, 25 + 2*i)
print(x)
y = AvgL_vs_eps
plt.figure(figsize=(4,2.2), dpi=140)
plt.plot(x, y)
plt.xlabel('epsilon')
plt.ylabel('r_ee')
plt.show()

!vmd -e vdw.vmd Distance_vs_Epsilon_Sim.gsd

