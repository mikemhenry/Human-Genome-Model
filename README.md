
# A Model for the Human Genome
This example shows a simple simulation of 24 block copolymers that is a potential model of the Human Genome in hoomd.

The human genome consists of 24 chromosomes with a total size of 3Gb.  Aiden et al. 2014 produced a 1kB resolution (below single gene resolution) Genome Contact Map.  Contact domains (∼185 kb) segregate into six nuclear subcompartments with distinct histone marks.

## We begin with a random polymer model
from: https://lost-contact.mit.edu/afs//umich.edu/user/j/o/joaander/Public/hoomd-web/doc/page_example_scripts.html

<a href="https://github.com/fergusonml/Human-Genome-Model/blob/master/create%20random%20polymer.ipynb">Random Polymer Model</a>
This simulates a random polymer of 100 beads.

<a href="https://github.com/fergusonml/Human-Genome-Model/blob/master/create%20random%20polymer2.ipynb">Random Polymer Model</a>
This simulates a block copolymer of 1000 beads.  It will run on a cpu or a gpu.  The can be changed via hoomd.context.initialize()

<a href="https://github.com/fergusonml/Human-Genome-Model/blob/master/human_genome_180kb_resolution.ipynb">A Model for the Human Genome at 100kb resolution</a>
This simulates 24 block copolymer of 1000 beads.  It will run on a cpu or a gpu.  The can be changed via hoomd.context.initialize()

<a href="https://github.com/fergusonml/Human-Genome-Model/blob/master/human_genome_10kb_resolution.ipynb">A Model for the Human Genome at 10kb resolution</a>
This simulates 24 block copolymer of 10000 beads.  It will run only on an NVIDIA GPU.

