#!/usr/bin/env python
# @author <ashtoncberger@utexas.edu>
# ------------------------------------------------
'''
TODO:
    *polymerase densities on operon method for operon
    *plotly implementation for viewing density evolution
'''
#--------
# imports
#--------
import bisect
import math
# import numpy as np
import pandas as pd
#import plotly
import random
import string
# import sys
from heapq import *
from itertools import groupby

#----------
# reference
#----------
nn_DNADNA = {
# Watson-Crick nearest neighbor parameters for DNA-DNA duplexes via SantaLucia & Hicks 2004
('AA','TT'):-1.00, ('AT','TA'):-0.88, ('AC','TG'):-1.56, ('AG','TC'):-1.28,
('TA','AT'):-0.58, ('TT','AA'):-1.00, ('TC','AG'):-1.28, ('TG','AC'):-1.44,
('CA','GT'):-1.45, ('CT','GA'):-1.28, ('CC','GG'):-1.84, ('CG','GC'):-2.17,
('GA','CT'):-1.30, ('GT','CA'):-1.44, ('GC','CG'):-2.24, ('GG','CC'):-1.84 
}

# Watson-Crick nearest neighbor parameters for DNA-RNA hybrids via Sugimoto et al 1995
nn_DNARNA = {
('AA','UU'):-0.20, ('AT','UA'):-0.60, ('AC','UG'):-1.60, ('AG','UC'):-1.50, 
('TA','AU'):-0.90, ('TT','AA'):-1.00, ('TC','AG'):-1.80, ('TG','AC'):-2.10,
('CA','GU'):-1.10, ('CT','GA'):-1.30, ('CC','GG'):-2.90, ('CG','GC'):-2.70,
('GA','CU'):-0.90, ('GT','CA'):-0.90, ('GC','CG'):-1.70, ('GG','CC'):-2.10 
}

#----------
# functions
#----------
# general utilities
def DNADNA_dG(duplex):
    seq, rseq = duplex
    seql = len(seq)
    dG = 1.96                               #initiation parameter
    if seq[-1] + rseq[-1] in ['AT', 'TA']:
        dG += 0.05                          #terminal AT penalty
    for i in range(seql - 1):
        dG += nn_DNADNA[(seq[i:i+2],rseq[i:i+2])]
    return dG

def DNARNA_dG(duplex):
    seq, rseq = duplex
    seql = len(seq)
    dG = 3.1                                #initiation parameter
    if seq[-1] + rseq[-1] in ['AU', 'TA']:
        dG += 0.05                          #terminal AT penalty
    for i in range(seql - 1):
        dG += nn_DNARNA[(seq[i:i+2],rseq[i:i+2])]
    return dG

def foldingRNA_dG(nascentseq):
    # eventually want to calcualte the dG of any hairpin sequences in the proximate nascent growing strand which may cause polymerase pauses
    # for now, just make this an arbitrary function of length
    return math.log(math.sqrt(len(nascentseq)))

def RNAP_RNA_dG(nascentseq):
    # a term that is supposed to affect a state's dG based on interactions between an RNAP and the nascent RNA sequence strand
    # supposedly non-changing between states so any use in computing differences between states will render a 0 for its contribution
    # will not be including it below in the RNAP's properties but wanted to have it represented in the code nonetheless
    return 0

def fasta_to_seqdic(openedfile):
    seqinfo = {}
    fiter = (x[1] for x in groupby(openedfile, lambda line: line[0] == ">"))
    for header in fiter:
        infolist = header.next()[1:].strip().split(';')
        for x in infolist:
            k,v = x.split(':')
            if ',' in v:
                seqinfo[k] = tuple(int(n) for n in v.split(','))
            elif v in ['3','5']:
                seqinfo[k] = int(v)
            else:
                seqinfo[k] = v
        seq = "".join(s.strip() for s in fiter.next())
        seqinfo['seq'] = seq
        yield seqinfo

def dnacomp(seq):
    ttable = string.maketrans('AGTC', 'TCAG')
    return seq.upper().translate(ttable)

def revcomp(seq):
    # nucs = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
    ttable = string.maketrans('AGTC', 'TCAG')
    return seq.upper()[::-1].translate(ttable)
    # comp = ''.join([nucs[x] for x in seq.upper()])
    # return comp[::-1]

def rnacomp(seq):
    ttable = string.maketrans('AGTC', 'UCAG')
    return seq.upper().translate(ttable)

def sample_Tau(K_j):
    U_j = random.uniform(0.0, 1.0)
    Tau_j = (1.0/K_j) * math.log(1.0/U_j)
    return Tau_j

def sample_polymerase_binding_promoter():
    return 1

def sample_t1_detachment():
    if random.uniform(0.0, 1.0) <= 0.1:
        return 1

def sample_t2_detachment():
    return 1

# simulation events
def attempt_docking(operon):
    print 'A polymerase is attempting to bind the promoter region...'
    if operon.promoter_available and sample_polymerase_binding_promoter():
        promote_EC(operon)
    operon._sim.heap_push((sample_Tau(operon._sim._params['K_dock']), operon, 'attempt_docking'))

def promote_EC(operon):
    print 'A polymerase has successfully bound the promoter region and formed an elongation complex!'
    RNAP = RNAPII(operon)
    operon._RNAPS.append(RNAP)
    operon._sim.heap_push((sample_Tau(RNAP.Kavg), RNAP, 'translocate_polymerase'))

def translocate_polymerase(RNAP):
    # if not RNAP.active:
        # sys.exit('Inactive polymerase was queued for an event on the simulation heap')

    # if RNAP.operon.promoter_available and 'attempt_docking' not in [event[2] for event in RNAP.operon._sim._heap]:
    #     # RNAP.operon._sim.heap_push((sample_Tau(RNAP.operon._sim._params['K_dock']), RNAP.operon, 'attempt_docking'))
    #     RNAP.operon._sim.heap_push((0, RNAP.operon, 'attempt_docking'))

    if RNAP.operon.promoter_available:
        RNAP.operon._sim.heap_push((0, RNAP.operon, 'attempt_docking'))

    # check polymerase's surroundings to see if another is in its way
    # if RNAP.is_blocked:
    #     RNAP.paused = 1
    #     RNAP.operon._sim.heap_pop()
    #     # for i in range(RNAP.queue):
    #     #     RNAP.operon._sim.heap_pop()
    #     RNAP.operon._sim.heap_push((sample_Tau(RNAP.Kavg), RNAP, 'translocate_polymerase'))
    #     return

    RNAP.translocate_brownian()

    # check to see if polymerase is in terminator regions and if so whether it should detach
    if RNAP.position in range(*RNAP.operon.t1) and sample_t1_detachment():
        RNAP.operon._sim.heap_push((0, RNAP, 'terminate_polymerase'))
        return
    if RNAP.position in range(*RNAP.operon.t2) and sample_t2_detachment():
        RNAP.operon._sim.heap_push((0, RNAP, 'terminate_polymerase'))
        return

    # if function makes it to this point, add polymerase back to the heap for its next event
    RNAP.operon._sim.heap_push((sample_Tau(RNAP.Kavg), RNAP, 'translocate_polymerase'))

def terminate_polymerase(RNAP):
    RNAP.active = 0
    RNAP.operon._RNAPS.remove(RNAP)
    print 'The polymerase at nucleotide position {} has detached from the operon at time {}!'.format(RNAP.position, RNAP.operon._sim._t)


#--------
# classes
#--------
class Simulation(object):
    '''
    Class that contains the parameters for running a simulation 
    Can also specify variables
    
    Parameter ideas: 
        temperature - 37oC?
        salt concentration of environment
        concentration of RNAPs
        time scale/increments
        elongation rate cutoff
        NTP concentration
        per operon initiation frequency (how many RNAP begin per s) (fange et al used 1/s, avg dist bt RNAP ~100bp)

    tadigotla had a cutoff to discriminate paused complexes from nonpaused complexes

    possible simulation types:
        monte carlo gillepsie algorithm?
        read this:
            http://people.math.umass.edu/~markos/697SC/ssa.pdf
        and look at stochpy's source code

    sympy text plot for polymerase positions

    self.K_dock = 4400
        -Dennis et al 2004 said 110 initiations/minute = 1.83333 init/sec is max initiation rate... and RNAP concentration at 1/2 maximal activity is 4.35 microM
        -free RNAP concentration is 
    self.NTP_concentration = 20 uM 
        -tadigotla et al 2006 dissociation constant for NTP binding to RNAP during transcription
    self.T = 310.15 K
        -37 degrees celsius to Kelvin - all NN parameters are defined at this temp
    self.Kb = 0.0019872041
        -Boltzmann's constant in terms of Kcal/mol*K, which are the units of the NN dG calcs
    self.KbT = self.Kb * self.T
        -added this to save me time
    '''

    global _default_params 

    _default_params = {
        'K_dock': 1.833,
        'K_maxA': 50,
        'K_dA': 38,
        'K_maxU': 18,
        'K_dU': 24,
        'K_maxG': 36,
        'K_dG': 62,
        'K_dC': 7,
        'K_maxC': 33,
        'NTPc': 20,
        'T': 310.15,
        'Kb': 0.0019872041,
        'KbT': 0.61633135161,
    }

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        # private
        self._params        = _default_params
        self._log           = ''
        self._t             = 0
        self._timer         = 0
        self._silent        = 0
        self._reactions     = {}
        self._heap          = []

    def __str__(self):
        pass

    def get_rate(self, event, subject):
        if event == 'attempt_docking': 
            return self._params['Kpoly']
        return subject.K

    def heap_push(self, reaction_tuple):
        if reaction_tuple[-1]: 
            heappush(self._heap, reaction_tuple)
            # TODO: self._environment._reactants - Vij (update Sis?)
            heapify(self._heap)

    def heap_pop(self):
        reaction_Tau, reaction_Subject, reaction_Event = heappop(self._heap)
        self._t += reaction_Tau
        self._timer += reaction_Tau
        self._log += 'Reaction event "{}" for subject "{}" at time {}\n'.format(reaction_Event, reaction_Subject.name, self._t)
        # print 'Reaction event "{}" for subject "{}" at time {}'.format(reaction_Event, reaction_Subject.name, self._t)
        self._reactions[reaction_Event](reaction_Subject)
        heapify(self._heap)

    def heap_sort(self):
        heapify(self._heap)

    def run(self, operon, duration=1200, increment=100):
        #initialize settings    
        operon._params = self._params
        # self._t = operon._t
        self._timer = 0
        self._log   = ''
        self._reactions = {
                            'attempt_docking': attempt_docking,
                            'translocate_polymerase': translocate_polymerase,
                            'terminate_polymerase': terminate_polymerase
                            }
        tc = increment

        # Gillespie direct method
        if self.algorithm == 'direct':
            while self._timer < duration:
                pass

        # Next reaction method
        elif self.algorithm == 'nrm':
            # species (reactant Lambda(i,j) matrix)  
                            # i = 4 = free RNAPs, promoters, terminators, tx bubble
            # self._reactants = np.array([[1, 1, 0, 0],  # j = 6 = promoter location/docking,
                                        # [0, 1, 0, 0],  #         promoter melting in,
                                        # [0, 1, 0, 1],  #         abortive initiation,
                                        # [0, 0, 0, 1],  #         elongation,
                                        # [0, 0, 1, 1],  #         terminating sequence encounter,
                                        # [0, 0, 0, 1]]) #         ultimate dissolution of tx bubble

            # initial event to begin simulation is attempt to bind polymerase to promoter
            if not self._heap:
                # self.heap_push((sample_Tau(self._params['K_dock']), operon, 'attempt_docking'))
                self.heap_push((0, operon, 'attempt_docking'))

            while self._timer < duration:
                self.heap_pop()
                if self._timer > tc:
                    # print 'Saving snapshot of RNAP positions at simulation time {}...'.format(round(self._timer,2))
                    for key,val in operon._RNAP_bins.items():
                        operon._data.append([key, val, round(self._timer,2)])
                    tc += increment
                    print operon 
                    print operon._RNAP_bins
                    # print 'RNAP data at time {}:'.format(round(self._timer,2))
                    # print '\tName\tKavg\tPforward\tPaused\tUpstream_Neighbor\tUplimit\tDownstream_Neighbor\tDownlimit'
                    # for r in operon._RNAPS:
                    #     print '\t', r.name, r.Kavg, r.Pforward, r.paused, r.upstream_neighbor, r.downstream_neighbor
            operon._t += self._timer
            operon._log += self._log
            print 'Simulation for operon {} has completed!'.format(operon.name)

        else:
            print 'Invalid simulation algorithm'

class Operon(object):
    '''
    Class for an operon instance to represent a template DNA strand during transcription

    attributes:
        name        = name of operon
        seq         = sequence of operon
        strand      = 5 or 3. Specifies whether the sequence is the 5' to 3' or the 3' to 5' representation
        length      = length of sequence
        motifs      = sequence motifs on the operon. should at the bare minimum contain promoters and terminators.
        simulating  = True or False, represents whether operon is currently behaving in a simulation or 
        elapsed     = time elapsed (in seconds) for simulation
        data        = nested list that will contain all information about the polymerase positions on the operon 
    '''
    def __init__(self, environment, **kwargs):
        for key, val in kwargs.items():
            setattr(self, key, val)
        self._sim               = environment
        self._t                 = 0
        self._RNAPS             = []
        self._data              = []
        self._log               = ''
        if self.strand == 5:
            self.Lseq           = self.seq.upper()
            self.Rseq           = revcomp(self.Lseq)
            self.template       = self.Rseq
        else:
            self.Rseq           = self.seq.upper()
            self.Lseq           = revcomp(self.Rseq)
            self.template       = self.Lseq
        self.length             = len(self.Lseq)
        self.binwidth           = self.length / 20
        # self.motifs = self.IdentifyMotifs()
        # self._promoter_status = self._sim._reactants[0][0]

    @property
    def _RNAP_coords(self):
        return [x.position for x in self._RNAPS]

    @property
    def _RNAP_census(self):
        return len(self._RNAPS)

    @property
    def _RNAP_bins(self):
        bind = {x:0 for x in range(0,self.length,self.binwidth)}
        for x in self._RNAP_coords:
            bin = (x//self.binwidth)
            val = (bin + 1) * self.binwidth
            bind[val] += 1
        bind = {k:float(v)/self._RNAP_census for k,v in bind.items()}
        return bind

    @property
    def promoter_available(self):
        binding_region = range(51)
        for x in self._RNAP_coords:
            if x in binding_region:
                return 0
        return 1

    def __str__(self):
        return 'Operon: {}\nTotal ECs:  {}'.format(self.name, str(len(self._RNAPS)))
  
    def display_local(self, x=30, window=30):
        start, stop = x-window, x+window
        strand = self.seq[start:stop]
        coords = [str(x) if x % 10 == 0 else ' ' for x in range(start,stop)]
        display = "5' " + strand + '\n' + '   ' + '|'*len(strand) + '\n' + "3' " + revcomp(strand) + '\n' + '   ' + ''.join(coords)
        print 'Operon: {}\nActive: {}\nTotal ECs:  {}\n\n{}'.format(self.name, str(self.simulating), str(len(self._RNAPS)), display)

    def reset(self):
        # BUG: if i empty the heap using the commented out line of code below, the next call to simulation.run skips straight to line 202 and throws an error instead of 
        #       adding another attempt_docking event
        # self._sim._heap = [] 
        self._RNAPS = []
        self._t = 0
        self._data = []
        self._log = ''

    def write_data(self, prefix):
        df_file = str(prefix) + '.csv'
        df = pd.DataFrame(self._data, columns=['RNAP', 'x', 't'])
        df.to_csv(prefix+'.csv', index=False)
        print 'Wrote polymerase positions to', df_file
        log_file = str(prefix) + '.log'
        with open(log_file, 'w') as w:
            w.write(self._log)
        w.close()
        print 'Wrote event log to', log_file

class RNAPII(object):
    '''
    Class for a rna polymerase/transcription bubble/elongation complex, which consists of a melted (unhybridized) 11-16 nt DNA duplex template enclosed by an RNAP and 
    stabilized by a nascent RNA molecule

    yager/von hippel 1991 says bubble is 17 +- 1

    moves 3' to 5'

    mechanism options: 
        1) powerstroke
            -'makes use of the energy released during phosphodiester bond formation to drive 
            translocation between states 0 and +1'
        2) brownian ratchet
            -'bidirectional translocational steps occur stochastically as a result of reversible
            thermal (Brownian) motion of the polymerase along DNA'
            -'phosphodiester bond formation rectifies the motion and resets the system into a 
            state with a longer transcript length, poised to incorporate the next NTP'

    master equation?
    arrhenius rate?

    duplex should be 2d list [[5'A,3'T], [G, C], [G, gap], etc]
    
    free energy of the state (m,n,b) of the EC is decomposed by:
        Gmnb = G(DNA-DNA)mnb + G(RNA-DNA)mnb + G(RNA-RNA)mnb + G(NS)

    self.template               = operon.template
    # self.polymerase           = None
    self.growing_strand         = begin with string, eventually make it a nascent seq object? is 9 nucleotides long in Bai's model.
    self.position               = nucleotide coordinate of upstream-most element in tx bubble
    self.m                      = length of RNA transcript
    self.n                      = translocational_state # -1, 0, +1 (+- 2,3,4... w/ brownian ratchet backtracking and hypertranslocation)
    self.x                      = number of template unpaired DNA bases upstream (to left)
    self.y                      = number of template unpaired DNA bases downstream (to right)
    self.h                      = length of template DNA-nascent RNA hybrid
    self.b                      = bubble configuration
    self.sb                     = transcription bubble size in nucleotides
    self.bubble_state           = state of EC
    self.bubble_energy          = self.calculate_dG(self.bubble)
    self.tx_energies            = {1: None,       # free energy minima associated with translocational
                                   0: None,       # states along the template for a fixed length of RNA
                                  -1: None}       # transcript m.
                                        # free energy diff between 2 neighboring translocational states n and n+1
                                        # is associated with (i) breaking DNA-DNA base pair in front and annealing the pair
                                        # behind the moving bubble, (ii) the change in the hybrid base pairing,
                                        # (iii) the folding of the RNA transcript protruding out of the exit channel, and
                                        # (iv) the possible change in RNA-DNA and RNAP-RNA interactions
    '''    
    def __init__(self, operon, position=random.randrange(0,50)):
        self.operon     = operon
        self.simulation = operon._sim
        self.active     = 1
        self.position   = position
        self.m              = 12
        self.n              = 0
        self.x              = 2
        self.y              = 9
        self.h              = 1
        self.b              = (self.x, self.h, self.y)
        self.bsize          = sum(self.b)
        self.EC_state       = (self.m, self.n, self.b)
        self.paused         = 0
        for key in self.operon._sim._params:
            setattr(self, key, self.operon._sim._params[key])

    @property
    def name(self):
        return 'RNAP at position ' + str(self.position)

    @property
    def growing_strand(self):
        return rnacomp(self.operon.template[:self.position])

    # @property
    # def is_blocked(self):
    #     for x in self.operon._RNAP_coords:
    #         if x in range(self.position+1,self.position+13):
    #             return 1
    #     return 0

    @property
    def upstream_neighbor(self):
        if self.position == max(self.operon._RNAP_coords):
            return 0 
        v = sorted(self.operon._RNAP_coords).index(self.position)
        return sorted(self.operon._RNAP_coords)[v+1]

    @property
    def downstream_neighbor(self):
        if self.position == min(self.operon._RNAP_coords):
            return 0 
        v = sorted(self.operon._RNAP_coords).index(self.position)
        return sorted(self.operon._RNAP_coords)[v-1]

    @property
    def frontlimit(self):
        if 0 < self.upstream_neighbor - self.position < 13:
            return self.upstream_neighbor - self.position
        else:
            return 13

    @property
    def backlimit(self):
        if not self.downstream_neighbor:
            return 0
        if 0 < abs(self.downstream_neighbor - self.position) < 10:
            return self.downstream_neighbor  - self.position 
        else:
            return -10

    @property
    def queue(self):
        return len([x for x in self.operon._RNAP_coords if x > self.position])

    # @property
    # def K0n1(self):
    #     return self.calc_Ki(0, -1)

    @property
    def K00(self):
        return self.calc_Ki(0, 0)

    # @property
    # def K01(self):
    #     return self.calc_Ki(0, 1)

    @property
    def Kpause(self):
        K_NTPmax = self.operon._sim._params['K_max'+self.nextbase(0)]
        K_NTPd = self.operon._sim._params['K_d'+self.nextbase(0)]
        NTPc = self.operon._sim._params['NTPc']
        return (K_NTPmax*NTPc)/(K_NTPd*(1+self.K00) + NTPc)

    # @property
    # def Kforward(self):
    #     K_NTPmax = self.operon._sim._params['K_max'+self.nextbase(1)]
    #     K_NTPd = self.operon._sim._params['K_d'+self.nextbase(1)]
    #     NTPc = self.operon._sim._params['NTPc']
    #     return (K_NTPmax*NTPc)/(K_NTPd*(1+self.K01) + NTPc)

    # @property
    # def Kbackward(self):
    #     K_NTPmax = self.operon._sim._params['K_max'+self.nextbase(-1)]
    #     K_NTPd = self.operon._sim._params['K_d'+self.nextbase(-1)]
    #     NTPc = self.operon._sim._params['NTPc']
    #     return (K_NTPmax*NTPc)/(K_NTPd*(1+self.K0n1) + NTPc)

    @property
    def Kforwards(self):
        fs = []
        for i in range(0,self.frontlimit):
            K_NTPmax = self.operon._sim._params['K_max'+self.nextbase(i+1)]
            K_NTPd = self.operon._sim._params['K_d'+self.nextbase(i+1)]
            NTPc = self.operon._sim._params['NTPc']
            Ki = self.calc_Ki(i, i+1)
            val = (K_NTPmax*NTPc)/(K_NTPd*(1+Ki) + NTPc)
            # setattr(self, K+str(i)+str(i+1), Ki)
            fs.append((i+1,val))
        return fs

    @property
    def Kbackwards(self):
        if self.paused or not self.backlimit:
            return [(0,0)]
        bs = []
        for i in range(0,self.backlimit,-1):
            K_NTPmax = self.operon._sim._params['K_max'+self.nextbase(i-1)]
            K_NTPd = self.operon._sim._params['K_d'+self.nextbase(i-1)]
            NTPc = self.operon._sim._params['NTPc']
            Ki = self.calc_Ki(i-1, i)
            val = (K_NTPmax*NTPc)/(K_NTPd*(1+Ki) + NTPc)
            # setattr(self, K+str(i)+'n'+str(i-1), self.calc_Ki(i, i-1))
            bs.append((abs(i-1),val))
        return bs

    @property
    def forward_moves(self):
        moves = sorted([[x[1]/float(sum(zip(*self.Kforwards)[-1])),x[0]] for x in self.Kforwards], reverse=True)
        for i in range(len(moves)-1):
            moves[i+1][0] += moves[i][0]
        return moves

    @property
    def backward_moves(self):
        moves = sorted([[x[1]/float(sum(zip(*self.Kbackwards)[-1])),x[0]] for x in self.Kbackwards], reverse=True)
        for i in range(len(moves)-1):
            moves[i+1][0] += moves[i][0]
        return moves

    @property
    def Klist(self):
        return [self.Kpause, sum(zip(*self.Kforwards)[-1]), sum(zip(*self.Kbackwards)[-1])]

    @property
    def Zc(self):
        return sum(self.Klist)

    @property
    def Pforward(self):
        return sum([x[-1] for x in self.Kforwards])/self.Zc

    @property
    def Pbackward(self):
        return sum([x[-1] for x in self.Kbackwards])/self.Zc

    @property
    def Ppause(self):
        return self.Kpause/self.Zc

    @property
    def Kavg(self):
        # return self.Zc/len(self.Klist)
        return self.Zc

    def calc_Ki(self, n1, n2):
        return math.exp((self.calc_dGi(n2) - self.calc_dGi(n1))/self.operon._sim._params['KbT'])

    def calc_dGi(self, n):
        # return DNADNA_dG(self.generate_bubble(n)) + DNARNA_dG(self.generate_hybrid(n)) + foldingRNA_dG(self.generate_nasceq(n))
        return DNADNA_dG(self.generate_bubble(n)) + DNARNA_dG(self.generate_hybrid(n))

    def generate_nasceq(self, n):
        return rnacomp(self.operon.template[self.position+n:self.position+10+n])

    def generate_bubble(self, n):
        return (dnacomp(self.operon.template[self.position+n:self.position+12+n]), self.operon.template[self.position+n:self.position+12+n])
    
    def generate_hybrid(self, n):
        return (self.operon.template[self.position+n:self.position+10+n], self.generate_nasceq(n))

    def nextbase(self,n):
        return rnacomp(self.operon.template[self.position+10+n])

    def translocate_powerstroke(self):
        pass

    def translocate_brownian(self):
        # try to permit no more than 1 consecutive backtrack
        if self.paused:
            if random.uniform(0.0, 1.0) <= self.Pforward:
                newdraw = random.uniform(0.0,1.0) 
                moves = self.forward_moves
                bisect.insort_left(moves, [newdraw, newdraw])
                self.position += moves[zip(*moves)[1].index(newdraw) + 1][1]
                self.paused = 0
            return 
        else:
            draw = random.uniform(0.0, 1.0)
            if draw <= self.Pforward:
                if len(self.forward_moves) > 1:
                    newdraw = random.uniform(0.0,1.0) 
                    moves = self.forward_moves
                    bisect.insort_left(moves, [newdraw, newdraw])
                    self.position += moves[zip(*moves)[1].index(newdraw) + 1][1]
                else:
                    self.position += 1
                self.paused = 0
                return
            elif self.Pforward < draw <= self.Pforward + self.Pbackward:
                if len(self.backward_moves) > 1:
                    newdraw = random.uniform(0.0,1.0) 
                    moves = self.backward_moves
                    bisect.insort_left(moves, [newdraw, newdraw])
                    self.position -= moves[zip(*moves)[1].index(newdraw) + 1][1]
                else:
                    self.position -= 1
                self.paused = 1
                return 
            else:
                self.paused = 1
                return 

# class NascentSeq(object):
#     '''
#     Class for a 5'-3' nascent sequence molecule growing from an RNAP in a transcription bubble during elongation
#     May eventually be implemented to calculate hairpin loop secondary structures that could cause RNAP pausing
#     '''
#     def __init__(self, length = 8):
#         self.hl = length        # hybrid length
#         self.sequence = None
#         self.hybrid = None
#         self.bpseq = {}         # dictionary containing secondary structure
