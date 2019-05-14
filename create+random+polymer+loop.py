
# coding: utf-8

# # Random polymer
# This example shows you how to run a simple polymer simulation in hoomdd.
# 
# Here is a script that generates a system of bead-spring polymers that self-assemble into a hex phase when run for a few million time steps. The polymers are A6B7A6 block copolymers in an implicit solvent. The script also shows a few examples of how writing python code in the script can be handy: here the concentration phi_P is a parameter and math operations are performed to calculate the length of the box.
# 
# For more information on the model in this script, see
# "Micellar crystals in solution from molecular dynamics simulations"
# J. Chem. Phys. 128, 184906 (2008); DOI:10.1063/1.2913522
# http://link.aip.org/link/?JCPSA6/128/184906/1
# 
# Any of the polymer systems in the paper could be easily run just by changing a few parameters in this script.
# from: https://lost-contact.mit.edu/afs//umich.edu/user/j/o/joaander/Public/hoomd-web/doc/page_example_scripts.html
# 
# ## Initialize
import hoomd
from hoomd import *
from hoomd import md
from hoomd import deprecated
#hoomd.context.initialize();
#hoomd.context.initialize("--mode=cpu");
hoomd.context.initialize("--mode=gpu");


# ## Import other libraries and define parameters
import math
# parameters (polymer physics)
phi_P = 0.25 # this has to do with how good of a solvent the polymer is in and is only used to calculate the box size below
n_poly = 1 # n_poly in the number of polymers in the simulation
n_beads = 100


# ## Define the polymer


# This polymer can be made up of two different types A and B.  It's called a block copolymer.
polymer1 = dict(bond_len=1.2, type=['B']*1 + ['A']*n_beads + ['B']*1,bond="linear", count=n_poly)
# perform some simple math to find the length of the box
N = len(polymer1['type']) * polymer1['count']
# generate the polymer system
system=deprecated.init.create_random_polymers(box=data.boxdim(volume=1000*math.pi * N / (6.0 * phi_P)),polymers=[polymer1], separation=dict(A=0.35, B=0.35),seed=1234)


# ## Setup the bonds and force fields
# force field setup
harmonic = md.bond.harmonic()
harmonic.bond_coeff.set('polymer', k=330.0, r0=0.84)
nl = md.nlist.cell();
lj = md.pair.lj(r_cut=3.0,nlist=nl)
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0, alpha=0.0)
lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0, alpha=0.0)
lj.pair_coeff.set('B', 'B', epsilon=100.0, sigma=1.0, alpha=1.0)


# ## Integrate the simulation
import numpy
import math
# integrate NVT for a bunch of time steps
all = group.all()
md.integrate.mode_standard(dt=0.005)
md.integrate.nvt(group=all, kT=1.2, tau=0.5)

#save a gsd polymer configuration every 100 timesteps
dump.gsd("gsd/random_polymer_loop.gsd", period=100, group=all, overwrite=True);

#log position of particle 0
polymers=['p1','p2','p3','p4','p5','p6','p7','p8','p9','p10']
log = hoomd.analyze.log(filename="ree.log",
                         quantities=polymers[0:n_poly],
                         period=100,
                         overwrite=True);
for i in range(n_poly):
    print(polymers[i])
    print((n_beads+2)*i)
    print((n_beads+2)*i+n_beads+1)
    log.register_callback(polymers[i], lambda timestep: math.sqrt((system.particles[(n_beads+2)*i].position[0]-system.particles[(n_beads+2)*i+n_beads+1].position[0])**2+(system.particles[(n_beads+2)*i].position[1]-system.particles[(n_beads+2)*i+n_beads+1].position[1])**2+(system.particles[(n_beads+2)*i].position[2]-system.particles[(n_beads+2)*i+n_beads+1].position[2])**2));

# running this on a GPU the resulting dump files should show the
# polymers self-assembling into the hexagonal phase
run(1e6)


# ## Plot results and correlation function
from matplotlib import pyplot
#get_ipython().magic(u'matplotlib inline')
#get_ipython().system(u'head ree.log')
data = numpy.genfromtxt(fname='ree.log',skip_header=True);
t=numpy.ndarray.flatten(data[:,0])
r_ee=data[:,1:]#/(101*(0.84**2))
print(numpy.shape(r_ee))
#pyplot.figure(figsize=(4,2.2), dpi=140);
#pyplot.plot(t,r_ee);
for i in range(n_poly):
    pyplot.figure(figsize=(4,2.2), dpi=140);
    pyplot.plot(t,r_ee[:,i]);
    pyplot.xlabel('time step');
    pyplot.ylabel('r_ee');
#r_ee=r_ee[len(r_ee)/2:]
def autocorr(x):
    result = numpy.correlate(x, x, mode='full')/(numpy.average(x)**2)-1.0
    return result[result.size/2:]
#pyplot.figure(figsize=(4,2.2), dpi=140);
for i in range(n_poly):
    pyplot.figure(figsize=(4,2.2), dpi=140);
    g=autocorr(numpy.ndarray.flatten(r_ee[:,i]))
    pyplot.semilogx(t,g)
    pyplot.xlabel('tau');
    pyplot.ylabel('g');


# The results of this example may be visualized with VMD which can be installed from here: http://www.ks.uiuc.edu/Research/vmd/
#!vmd -e vdw.vmd gsd/random_polymer_loop.gsd

