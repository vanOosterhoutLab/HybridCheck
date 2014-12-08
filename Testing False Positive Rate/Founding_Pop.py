#! /usr/bin/python

# Generate population and evolve them for 100,000 generations, for testing HybRIDS capabilities.
# Sequences of 50,000 base pairs.

# First command line arg is the base filename for saves.

import sys
import simuOpt
simuOpt.setOptions(numThreads=12, quiet=False, alleleType='long')
import simuPOP as sim
from simuPOP.utils import importPopulation, export


# Initialize the homogenous population and get it to the right starting size.
pop = sim.Population(size=[1], ploidy=2, loci=50000, chromNames=['Chr1'],
  alleleNames=['A', 'C', 'T', 'G'])
pop.setGenotype(1)
pop.resize(1500, propagate=True)
pop.splitSubPop(subPop = 0, sizes = [500,500,500])
print(pop.numSubPop())

# Main evolver function.
pop.evolve(
  initOps=[
    sim.InitSex(maleFreq = 0.5),
    sim.PyOutput('#seg_sites,#fixed_sites,#generations\n', output = ">>FoundingSegSites.txt")
    ],
  preOps=[
    sim.AcgtMutator(rate=10e-6, model='JC69')
  ],
  matingScheme=sim.RandomMating(ops = sim.Recombinator(rates = [10e-6]), subPopSize = [500,500,500]), 
  postOps=[
    sim.Stat(numOfSegSites=sim.ALL_AVAIL,
            vars=['numOfSegSites', 'numOfFixedSites'], step = 100),
    sim.utils.Exporter(format='GENEPOP', output="!'Founding/GENEPOP/Founding%d_GENEPOP' % (gen)", gui=False, step=1000),
    sim.utils.Exporter(format='csv', output="!'Founding/CSV/Founding%d.csv' % (gen)", gui=False, step=1000)
  ],
  gen=10001
)

