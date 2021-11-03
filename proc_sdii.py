#!/usr/bin/python
import numpy.random as np_random
from scipy.special import binom
from sets import Set
import numpy as np
import itertools
import math
import csv
import time
import sys

from .sdii import sdii


# bootstrap data generation
# index_list: each elemenet is a sampling with replacement
def bootstrap_data(data, index_list, nvar):
    data1 = data[index_list[0], 0]
    for i in range(1, nvar):
        col = data[index_list[i], i]
        data1 = np.c_[data1, col]
    return data1



# indicator function on bootstrap data*
#       count how many T_l >= user input threshold
# t: user threshold
def G_Bn(sdii_obj, bootstrap_indexSet, t, varset, order):
    B = len(bootstrap_indexSet) # number of bootstrap
    count = 0
    for i in range(B):
        #print 'G_Bn()::bootstrap # %d' % i
        t0 = time.time()
        #print 'G_Bn()::bootstrap index[1:10]: %s' % str(bootstrap_index[i][0:10])
        data_1 = bootstrap_data(sdii_obj.data, bootstrap_indexSet[i], len(varset))
        sdii_bootstrap = sdii(data_1) # new hashing object for new data
        '''
        print 'G_Bn()::data_1 shape: %s' % repr(data_1.shape)
        print 'G_Bn()::data_1 : %s' % repr(data_1)
        print
        print 'G_Bn()::data : %s' % repr(data)
        exit()
        '''
        for s in set(itertools.combinations(varset, order)): # generate all variable subset with length of 2
        # varset = Set([2,4,6]), order = 2
        # set([(2, 6), (2, 4), (4, 6)])
            if sdii_bootstrap.T_l(list(s)) >= t: # using the hash table in sdii_bootstrap
                count+=1
        t1 = time.time()

    print(('G_Bn():: # of T >= t : %d, t: %f, count*(1/B): %f' % (count, t, (1.0/B)*count)))
    return (1.0/B)*count


# max(1:pk if Tl>=t count++, 1)
def max_Tl_1(sdii_obj, t, varset, order):
    count = 0
    for s in set(itertools.combinations(varset, order)): # generate all variable subset with length of 2
    # varset = Set([2,4,6]), order = 2
    # set([(2, 6), (2, 4), (4, 6)])
        T = sdii_obj.T_l(list(s)) # hashing for real data
        #print 'max_Tl_1()::T: %f' % T
        if T >= t:
            count+=1

    if count < 1:
        print(('max_Tl_1():: # of T >= t: %d, change to 1' % count))
        count = 1.0
        return count

    print(('max_Tl_1():: # of T >= t: %d' % count))
    return count*1.0


# find threshold with boostrap
def threshold_t_B(sdii_obj, alpha, varset, order, B):
    sk = 4
    pk = binom(len(varset), order)
    #top = 2*math.sqrt(sk*math.log(pk))
    top = 1.0
    #print 'threshold_t_B()::'
    #print (len(alphabet), sk, pk, top)

    final_t = 0.0
    min_diff = sys.float_info.max

    n = sdii_obj.data.shape[0]

    # a list of lists
    bootstrap_indexSet = []
    for b in range(0,B):
        single_var_idx_list = []
        # column-wised samepling, for ith variable
        for i in range(0, len(varset)):
            single_var_idx_list.append(np_random.choice(n, n, replace=True))
        bootstrap_indexSet.append(single_var_idx_list)
        #print 'threshold_t_B()::the %dth bootstrap: %s' % (b, repr(single_var_idx_list))

    # get inf(t<=alpha) from all t
    #for t in np.linspace(0.1,top,10):
    for t in np.linspace(0.0, 0.05, 100):
        v_G = G_Bn(sdii_obj, bootstrap_indexSet, t, varset, order)
        #print 'threshold_t_B()::data: %s' % repr(data[1:10,:])
        v_m = max_Tl_1(sdii_obj, t, varset, order)
        ratio = v_G/v_m
        print(('threshold_t_B():: G/Max_T ratio: %f\n' % ratio))
        if ratio < 1e-10:
            print('zero ratio reached. break')
            break
        diff = alpha - ratio
        if diff > 0 and diff < min_diff:
            min_diff = diff
            final_t = t
        #break # test sampling

    if final_t == 0:
        final_t = top

    print(('threshold_t_B()::final t: %f, v_G: %f, v_m: %f' % (final_t, v_G, v_m)))
    return final_t


# forward selection procedure
# return a set of significant variables (index)
def forward_selection(data, alpha, varset, order, B):
    global alphabet
    ret_varset = Set()

    print(('forward_selection()::varset: %s, order: %d' % (repr(varset), order)))

    sdii_core = sdii(data)
    th = threshold_t_B(sdii_core, alpha, varset, order, B)
    print(('forward_selection()::threshold of order [%d]: %f' % (order, th)))

    return th



#alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
alphabet = []
#alphabet = ['X','Y','Z','U','V','W']
#alphabet = ['X1','X2','X3','X4','X5']

def main():
    global alphabet

    aa_alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T']
    na_alphabet = [
            'AB', 'AE', 'CB', 'CE', 'DB', 'DE', 'EB', 'EE', 'FB', 'FE', 'GB', 'GE', 'HB', 'HE', 'IB', 'IE',
            'KB', 'KE', 'LB', 'LE', 'MB', 'ME', 'NB', 'NE', 'PB', 'PE', 'QB', 'QE', 'RB', 'RE', 'SB', 'SE',
            'TB', 'TE', 'VB', 'VE', 'WB', 'WE', 'YB', 'YE'
    ]

    if len(sys.argv) < 3:
        print('Usage: python proc_sdii.py var_type score_file')
        return

    vartype = sys.argv[1]
    if vartype == 'AA':
        alphabet = aa_alphabet
        print(('use AA varset : %s' % repr(alphabet)))
    elif vartype == 'NA':
        alphabet = na_alphabet
        print(('use NA varset : %s' % repr(alphabet)))

    scorefile = sys.argv[2]
    print(('score file: %s' % scorefile))
    outfile = '%s.sdii' % scorefile
    print(('write to %s' % outfile))

    score = np.loadtxt(scorefile, delimiter=',')
    #print score.shape[0]

    '''
    t1 = time.time()
    varset = range(len(alphabet))
    th2 = forward_selection(score, 0.1, varset, 2, 300)
    th3 = forward_selection(score, 0.1, varset, 3, 300)
    t2 = time.time()

    print 'Threshold of order 2: %f' % th2
    print 'Threshold of order 3: %f' % th3
    print 'use %d seconds' % (t2 - t1)

    return
    '''

    sdii_core = sdii(score)
    fout = open(outfile, 'w')
    print('calculating mutual information ...')
    t0 = time.time()
    for s in set(itertools.combinations(list(range(len(alphabet))), 2)): # generate all variable subset with length of 2
        fout.write('%s %.15f\n' % ('-'.join([(alphabet[i]) for i in s]), sdii_core.calc_sdii(list(s))))

    t1 = time.time()
    print(('MI time: %d seconds' % (t1-t0)))

    print('calculating DeltaK(3) ...')
    for s in set(itertools.combinations(list(range(len(alphabet))), 3)): # generate all variable subset with length of 3
        fout.write('%s %.15f\n' % ('-'.join([(alphabet[i]) for i in s]), sdii_core.calc_sdii(list(s))))
    t2 = time.time()
    print(('DeltaK(3) time: %d seconds' % (t2-t1)))



if __name__=="__main__":
    main()
