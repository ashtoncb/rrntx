#!/usr/bin/env python
# @author <ashtoncberger@utexas.edu>
# ------------------------------------------------
from src import models

#create simulation
sim_args = {
            'mechanism':'brownian_ratchet',
            'multibubbles':False, 
            'folding':False, 
            'algorithm':'nrm'
            } 
sim = models.Simulation(**sim_args)

#create operon to run simulation on
seqinfo = models.fasta_to_seqdic(open('operons/rrnB.fa', 'r')).next()
rrnB = models.Operon(environment=sim, **seqinfo)

#run simulation and write results
sim.run(operon=rrnB, duration=1200)
sim.write_data('tmp')
