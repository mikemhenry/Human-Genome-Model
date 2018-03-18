
# A Model for the Human Genome
This example shows a simple simulation of 24 block copolymers that is a good model of the Human Genome in hoomd.

Hoomd can be downloaded from: http://glotzerlab.engin.umich.edu/hoomd-blue/
On Linux and MacOS type:
>$ conda config --add channels glotzer
>$ conda install hoomd

Unfortunatly Hoomd does not currently compile or run on Windows.
If GPU computation is required on MacOS Hoomd must be compiled from source.

![](https://upload.wikimedia.org/wikipedia/commons/6/6e/PLoSBiol3.5.Fig1bNucleus46Chromosomes.jpg)

![](https://upload.wikimedia.org/wikipedia/commons/4/4b/Chromatin_Structures.png)
Adapted from:https://en.wikipedia.org/wiki/Chromosome

The human genome consists of 24 chromosomes with a total size of 3Gb.  Aiden et al. 2014 produced a 1kB resolution Genome Contact Map(below single gene resolution).  Contact domains (∼185 kb) segregate into six nuclear subcompartments with distinct histone marks.

## We begin with a random polymer model
These models are adapted from: 
https://lost-contact.mit.edu/afs//umich.edu/user/j/o/joaander/Public/hoomd-web/doc/page_example_scripts.html

For more information on the model in this script, see "Micellar crystals in solution from molecular dynamics simulations" J. Chem. Phys. 128, 184906 (2008); DOI:10.1063/1.2913522 http://link.aip.org/link/?JCPSA6/128/184906/1

## Models

1. <a href="https://github.com/fergusonml/Human-Genome-Model/blob/master/create%20random%20polymer.ipynb">Random Polymer Model</a>
This simulates a random polymer of 100 beads.
![](snapshots/polymer.png)

2. <a href="https://github.com/fergusonml/Human-Genome-Model/blob/master/create%20random%20polymer2.ipynb">Block CoPolymer Model</a>
This simulates a block copolymer of 1000 beads.  It will run on a cpu or a gpu which can be changed using the command >hoomd.context.initialize(--mode=cpu)
or
>hoomd.context.initialize(--mode=gpu)

Once context is initialized the python kernel must be restarted to change it.

![](snapshots/copolymer.png)

3. <a href="https://github.com/fergusonml/Human-Genome-Model/blob/master/human_genome_180kb_resolution.ipynb">A Model for the Human Genome at 100kb resolution</a>
This simulates 24 block copolymers of 1000 beads each.  Here each bead approximates a single contact domain (e.g. an enhancer promoter loop).  It will run on a cpu or a gpu.  Which can be changed via hoomd.context.initialize() as described above.
![](snapshots/hum_gen_100kb.png)

4. <a href="https://github.com/fergusonml/Human-Genome-Model/blob/master/human_genome_10kb_resolution.ipynb">A Model for the Human Genome at 10kb resolution</a>
This simulates 24 block copolymer of 10000 beads.  Here each bead approximates a single gene.  This model will run only on an NVIDIA GPU.
![](snapshots/hum_gen_10kb_before.png)
![](snapshots/hum_gen_10kb_after.png)
![](snapshots/hum_gen_10kb.png)

## Acknowledgements
I'd like to thank Joshua Anderson of the University of Michigan for writing Hoomd and helpful discussions and Eric Jankowski of Boise State University for helpful discussions.
