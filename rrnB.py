#!/usr/bin/env python
# @author <ashtoncberger@utexas.edu>
# ------------------------------------------------
'''
live streams from this simulation are hosted at:
https://plot.ly/~ashtoncb/0/simulated-polymerase-densities-on-the-rrnb-operon-at-time-00/
'''

from src import models

#create simulation
sim_args = {
            # 'quiet': 1,
            'mechanism':'brownian_ratchet',
            'multibubbles':False, 
            'folding':False, 
            'algorithm':'nrm'
            } 
sim = models.Simulation(**sim_args)

#create operon to run simulation on
seqinfo = models.fasta_to_seqdic(open('operon_fastas/rrnB.fa', 'r')).next()
rrnB = models.Operon(environment=sim, **seqinfo)

#run simulation and write results
sim.run(operon=rrnB, duration=7200)
rrnB.write_data('rrnBt1')
