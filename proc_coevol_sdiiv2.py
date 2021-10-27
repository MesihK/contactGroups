#!/usr/bin/python
import numpy as np
import itertools
import math
import time
import sys

from .sdii import sdii
from .msa import msa
from scipy.special import binom

#alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
alphabet = []
#alphabet = ['X','Y','Z','U','V','W']
#alphabet = ['X1','X2','X3','X4','X5']

def main():
    global alphabet

    if len(sys.argv) < 4:
        print('Simplified (single process) co-evol sdii calculation procedure. (No weighting / data preprocessing)')
        print('Usage: python proc_coevol_sdiiv2.py msafile target_seq msapos order')
        print('Example: python proc_coevol_sdiiv2.py sca_pdz.fa 77 all 3')
        return

    msafile = sys.argv[1]
    targetHeader = sys.argv[2]
    target = sys.argv[3].lower()
    order = int(sys.argv[4])

    print(('msafile: [%s]' % msafile))
    print(('target msa header: [%s]' % targetHeader))
    print(('target var: [%s]' % target))
    print(('order: [%d]' % order))

    outfile = '%s.%s_%d_sdii' % (msafile, target, order)
    print(('write to [%s]' % outfile))

    m = msa(msafile)
    m.setTarget(targetHeader)
    print(('data dimension: (%d, %d)' % (m.seqNum, m.seqlen)))
    #weight_cutoff = 0.3 # for weighting msa sequence # taken care of in matlab
    drop_cutoff = 0.0
    score, varlist = m.msaboard(drop_cutoff) #, weight_cutoff) # return a compact score
    print(('reduced data dimension: %s' % repr(score.shape)))

    '''
    score: A..C..D.EF
    index: 0123456789
    # after reduction
    score: ACDE
    index: 0123 -> input in sdii calculation
    index: 0368 = varlist = alphabet
    '''

    alphabet = [str(i) for i in varlist]
    #print alphabet
    #m.writeScoreboard('1k2p_PF07714_seed.score')
    if (target != 'all') and (int(target) not in varlist):
        print(('The alignment for var %s is not significant. exit.' % target))
        return

    if target == 'all':
        pk = binom(len(varlist), order)
    else:
        pk = binom(len(varlist), order-1) - len(varlist) - 1

    print(('total calculations: %d' % pk))

    #print 'Loading weight ...'
    #pfam_weight = np.loadtxt(weightfile, delimiter=',')
    #print 'Weight vector: %s' % repr(pfam_weight.shape)

    sdii_core = sdii(score)
    #print 'Applying weight to sdii data ...'
    #sdii_core.setWeight(pfam_weight) # set sequence weight
    print('calculating sdii ...')
    fout = open(outfile, 'w')
    t0 = time.time()
    count = 0
    for s in set(itertools.combinations(list(range(len(alphabet))), order)):
        if (target == 'all') or (alphabet.index(target) in s):
            count+=1
            print(('%d/%d: %s          ' % (count, pk, '-'.join([(alphabet[i]) for i in s]))))
            ret_sdii = sdii_core.calc_sdii(list(s))
            t1 = time.time()
            print(('time used: %d seconds\n' % (t1-t0)))
            fout.write('%s %.15f\n' % ('-'.join([(alphabet[i]) for i in s]), ret_sdii))
            t0 = t1

    fout.close()


if __name__=="__main__":
    main()
